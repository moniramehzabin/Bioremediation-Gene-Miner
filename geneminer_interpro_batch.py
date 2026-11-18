"""Fault-isolated batch orchestration around the installed v6 sequence adapter."""
from pathlib import Path
import importlib.util
import json
import time
import csv

TERMINAL = {"FAILURE", "ERROR", "NOT_FOUND", "SUBMITTING", "SUBMISSION_UNCERTAIN"}


def _adapter():
    path = Path(__file__).with_name("geneminer_interpro_v6.py")
    spec = importlib.util.spec_from_file_location("gm_v6_sequence_adapter", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def run_interpro_batch(fasta_path, email, outdir, poll_seconds=30,
                       timeout_minutes=180, *, offline=False, cache_dir=None):
    adapter = _adapter()
    fasta_path, outdir = Path(fasta_path), Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    records, tasks = adapter.read_fasta(fasta_path)
    if not records:
        (outdir / "interpro_web_status.txt").write_text("NOT_RUN\tNo Review candidates.\n", encoding="utf-8")
        return None
    cache = Path(cache_dir) if cache_dir else Path(__file__).parent / "interpro_v6_cache"
    deadline = time.monotonic() + max(1, float(timeout_minutes)) * 60
    manifest = {"adapterVersion": "batch-0.1.0", "endpoint": adapter.ENDPOINT,
                "settings": adapter.SETTINGS, "inputSHA256": adapter.sha(fasta_path),
                "complete": False, "proteins": [],
                "format": "adapter-generated 15-column TSV; authoritative raw JSON retained"}
    rows, versions = [], set()
    service_problem = None
    batch_dir = outdir / "interpro_v6_batch"
    batch_dir.mkdir(exist_ok=True)

    def read_job(key):
        path = cache / (key + ".job.json")
        if not path.exists():
            return {"status": "NOT_SUBMITTED", "settings": adapter.SETTINGS}
        return json.loads(path.read_text(encoding="utf-8"))

    with adapter.cache_lock(outdir / "interpro_batch_lock"):
        for index, (key, task) in enumerate(tasks.items(), 1):
            print("[info] InterPro batch", index, "/", len(tasks), ",".join(task["ids"]), flush=True)
            item = {"ids": task["ids"], "sha256": key, "status": "NOT_SUBMITTED"}
            try:
                job = read_job(key)
                item.update(status=job.get("status", "UNKNOWN"), jobId=job.get("jobId"))
                if job.get("settings") != adapter.SETTINGS:
                    item.update(status="SETTINGS_MISMATCH", note="Different cached analysis settings; not resubmitted.")
                elif item["status"] in TERMINAL:
                    item["note"] = "Existing failed/uncertain job retained; not automatically resubmitted. Other proteins continue."
                elif item["status"] != "FINISHED":
                    if offline:
                        item["note"] = "Offline check: no submission or polling performed."
                    elif service_problem:
                        item["note"] = "Further remote work paused after a connection/nonterminal error: " + service_problem
                    elif time.monotonic() >= deadline:
                        item["note"] = "Batch waiting budget reached; rerun to resume."
                    else:
                        per_dir = batch_dir / key
                        per_dir.mkdir(exist_ok=True)
                        per_fasta = per_dir / "query.faa"
                        per_fasta.write_text("".join(
                            ">" + identifier + "\n" + task["sequence"] + "\n" for identifier in task["ids"]
                        ), encoding="ascii")
                        remaining = max(0.01, (deadline - time.monotonic()) / 60)
                        try:
                            adapter.run_interpro_v6(per_fasta, email, per_dir,
                                poll_seconds, remaining, cache_dir=cache)
                        except Exception as exc:
                            job = read_job(key)
                            item.update(status=job.get("status", "UNKNOWN"), jobId=job.get("jobId"),
                                        note=type(exc).__name__ + ": " + str(exc))
                            # Terminal per-sequence failure is isolated. A connection
                            # problem stops new requests, while cached results still load.
                            if item["status"] not in {"FAILURE", "ERROR", "NOT_FOUND"}:
                                service_problem = item["note"]
                        else:
                            job = read_job(key)
                            item.update(status=job.get("status", "UNKNOWN"), jobId=job.get("jobId"))
                if item["status"] == "FINISHED":
                    result_path = cache / (key + ".result.json")
                    result = json.loads(result_path.read_text(encoding="utf-8"))
                    record, version, release = adapter.validate_result(result, task["sequence"])
                    if (job.get("version"), job.get("release")) != (version, release):
                        raise RuntimeError("Cache version metadata differs from result.")
                    if versions and (version, release) not in versions:
                        item.update(status="VERSION_MISMATCH", note="Completed result retained but not mixed with another data/software release.")
                    else:
                        protein_rows = list(adapter.tsv_rows(record, task["ids"], job["completedUTC"]))
                        versions.add((version, release))
                        rows.extend(protein_rows)
                        item.update(version=version, release=release,
                                    signatureMatches=len(record["matches"]),
                                    resultFile=str(result_path), resultSHA256=adapter.sha(result_path))
            except Exception as exc:
                item.update(status="RESULT_INVALID", note=type(exc).__name__ + ": " + str(exc))
            manifest["proteins"].append(item)
            adapter.write_json(outdir / "interpro_v6_batch_progress.json", manifest)
            print("[info] Protein status:", item["status"], flush=True)

        completed_items = [x for x in manifest["proteins"] if x["status"] == "FINISHED"]
        completed_ids = {identifier for x in completed_items for identifier in x["ids"]}
        counts = {}
        for item in manifest["proteins"]:
            counts[item["status"]] = counts.get(item["status"], 0) + len(item["ids"])
        manifest.update(complete=completed_ids == set(records), proteinCount=len(records),
                        completedProteinCount=len(completed_ids), uniqueSequenceCount=len(tasks),
                        locationRows=len(rows), proteinStatusCounts=counts)
        path = outdir / "interpro_web.tsv"
        temporary = path.with_suffix(".tsv.tmp")
        with temporary.open("w", encoding="utf-8", newline="") as handle:
            csv.writer(handle, delimiter="\t", lineterminator="\n").writerows(rows)
        temporary.replace(path)
        manifest["tsvSHA256"] = adapter.sha(path)
        adapter.write_json(outdir / "interpro_v6_completion.json", manifest)
        with (outdir / "interpro_v6_protein_status.tsv").open("w", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(["Protein", "Status", "JobID", "Note"])
            for item in manifest["proteins"]:
                for identifier in item["ids"]:
                    writer.writerow([identifier, item["status"], item.get("jobId") or "", item.get("note", "")])
        status = "FINISHED" if manifest["complete"] else "PARTIAL" if completed_ids else "INCOMPLETE"
        (outdir / "interpro_web_status.txt").write_text(status + "\t" + json.dumps(counts) + "\n", encoding="utf-8")
        (outdir / "interpro_job_id.txt").write_text("\n".join(
            item["jobId"] for item in manifest["proteins"] if item.get("jobId")) + "\n", encoding="utf-8")
        print("[info] InterPro batch:", status, "completed", len(completed_ids), "of", len(records), flush=True)
        return str(path)
