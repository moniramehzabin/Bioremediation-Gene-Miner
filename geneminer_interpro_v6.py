"""InterProScan 6 adapter for Gene Miner; reviewed JSON remains authoritative.

The public EMBL-EBI MCP tools analyzeSequence/getSequenceAnalysis are used.
No legacy v5 service, browser automation, or precomputed-only substitute.
"""
from __future__ import annotations

import csv
import hashlib
import json
import os
import time
import urllib.error
import urllib.request
from contextlib import contextmanager
from datetime import datetime, timezone
from pathlib import Path

ADAPTER_VERSION = "0.1.0"
ENDPOINT = "https://www.ebi.ac.uk/interpro/mcp"
SETTINGS = {"applications": "DEFAULT", "goTerms": True, "pathways": True}


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def utc():
    return datetime.now(timezone.utc).isoformat()


def write_json(path, data):
    path = Path(path)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(data, indent=2), encoding="utf-8")
    os.replace(temporary, path)


def read_fasta(path):
    from Bio import SeqIO
    records, tasks = {}, {}
    for record in SeqIO.parse(str(path), "fasta"):
        if record.id in records:
            raise RuntimeError("Duplicate FASTA identifier: " + record.id)
        sequence = str(record.seq).upper()
        if not sequence or set(sequence) - set("ACDEFGHIKLMNPQRSTVWYXBZUO"):
            raise RuntimeError("Unsupported sequence characters: " + record.id)
        records[record.id] = sequence
        key = hashlib.sha256(sequence.encode("ascii")).hexdigest()
        tasks.setdefault(key, {"sequence": sequence, "ids": []})["ids"].append(record.id)
    return records, tasks


@contextmanager
def cache_lock(directory):
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    path = directory / "ACTIVE_RUN.lock"
    try:
        fd = os.open(path, os.O_CREAT | os.O_EXCL | os.O_WRONLY)
    except FileExistsError:
        raise RuntimeError(
            "Another adapter run may be active: " + str(path)
            + ". If a process was force-killed, inspect the lock before removing it."
        ) from None
    try:
        with os.fdopen(fd, "w", encoding="utf-8") as handle:
            json.dump({"pid": os.getpid(), "createdUTC": utc()}, handle)
        yield
    finally:
        path.unlink(missing_ok=True)


class MCPClient:
    def __init__(self):
        self.session = None
        self.protocol = "2025-03-26"
        self.counter = 0
        initial = self.rpc("initialize", {
            "protocolVersion": self.protocol, "capabilities": {},
            "clientInfo": {"name": "GeneMiner-InterPro-v6", "version": ADAPTER_VERSION},
        })
        self.protocol = initial.get("protocolVersion", self.protocol)
        self.rpc("notifications/initialized", notification=True)
        names = {t["name"] for t in self.rpc("tools/list").get("tools", [])}
        if not {"analyzeSequence", "getSequenceAnalysis"}.issubset(names):
            raise RuntimeError("Official v6 analysis tools are unavailable.")

    def rpc(self, method, params=None, notification=False):
        payload = {"jsonrpc": "2.0", "method": method}
        if params is not None:
            payload["params"] = params
        if not notification:
            self.counter += 1
            payload["id"] = self.counter
        headers = {"Content-Type": "application/json",
                   "Accept": "application/json, text/event-stream"}
        if self.session:
            headers["Mcp-Session-Id"] = self.session
        if method != "initialize":
            headers["MCP-Protocol-Version"] = self.protocol
        request = urllib.request.Request(
            ENDPOINT, data=json.dumps(payload).encode("utf-8"), headers=headers,
            method="POST",
        )
        # A submission is never blindly retried after transport uncertainty.
        with urllib.request.urlopen(request, timeout=180) as response:
            self.session = response.headers.get("Mcp-Session-Id", self.session)
            body = response.read().decode("utf-8-sig")
        if notification and not body.strip():
            return {}
        try:
            messages = [json.loads(body)]
        except json.JSONDecodeError:
            messages = []
            for block in body.replace("\r\n", "\n").split("\n\n"):
                data = "\n".join(line[5:].lstrip() for line in block.splitlines()
                                 if line.startswith("data:"))
                if data:
                    messages.append(json.loads(data))
        for message in messages:
            if message.get("id") == payload.get("id"):
                if "error" in message:
                    raise RuntimeError(json.dumps(message["error"]))
                return message.get("result", {})
        if notification:
            return {}
        raise RuntimeError("No matching MCP JSON-RPC response.")

    def call(self, name, arguments, evidence_path):
        reply = self.rpc("tools/call", {"name": name, "arguments": arguments})
        write_json(evidence_path, reply)
        if reply.get("isError"):
            raise RuntimeError("MCP tool error; inspect " + str(evidence_path))
        structured = reply.get("structuredContent")
        if isinstance(structured, dict) and any(
            k in structured for k in ("jobId", "status", "error")
        ):
            return structured
        for block in reply.get("content", []):
            if block.get("type") == "text":
                try:
                    value = json.loads(block["text"])
                except (ValueError, TypeError):
                    continue
                if isinstance(value, dict) and any(
                    k in value for k in ("jobId", "status", "error")
                ):
                    return value
        raise RuntimeError("Unrecognized MCP tool result; raw response preserved.")



REST_ENDPOINT = "https://www.ebi.ac.uk/Tools/services/rest/iprscan6"


class RestClient:
    def __init__(self, email):
        self.email = email

    def _request(self, path, form=None):
        from urllib.parse import urlencode
        from urllib.request import Request, urlopen

        body = urlencode(form).encode("ascii") if form is not None else None
        request = Request(REST_ENDPOINT + path, data=body)
        with urlopen(request, timeout=60) as response:
            return response.read().decode("utf-8-sig").strip()

    def call(self, operation, args, evidence_path):
        import re

        if operation == "analyzeSequence":
            job_id = self._request("/run", {
                "email": self.email,
                "sequence": args["sequence"],
                "goterms": str(args["goTerms"]).lower(),
                "pathways": str(args["pathways"]).lower(),
            })
            if not re.fullmatch(r"iprscan6-[A-Za-z0-9_-]+", job_id):
                raise RuntimeError("REST submission did not return a valid v6 job ID.")
            response = {"jobId": job_id, "status": "QUEUED"}
        elif operation == "getSequenceAnalysis":
            job_id = args["jobId"]
            if not re.fullmatch(r"iprscan6-[A-Za-z0-9_-]+", job_id):
                raise RuntimeError("Saved job ID is not a valid v6 job ID.")
            status = self._request("/status/" + job_id).upper()
            response = {"jobId": job_id, "status": status}
            if status == "FINISHED":
                response["result"] = json.loads(self._request("/result/" + job_id + "/json"))
        else:
            raise ValueError(f"Unsupported InterPro operation: {operation}")

        write_json(evidence_path, response)
        return response

    def close(self):
        pass


def validate_result(result, sequence):
    if not isinstance(result, dict):
        raise RuntimeError("Missing InterProScan JSON object.")
    version = str(result.get("interproscan-version", ""))
    release = str(result.get("interpro-version", ""))
    if not version.startswith("6.") or not release:
        raise RuntimeError("Cannot verify v6 version and InterPro release.")
    records = result.get("results")
    if not isinstance(records, list) or len(records) != 1:
        raise RuntimeError("Expected one result sequence per job.")
    record = records[0]
    actual = "".join(str(record.get("sequence", "")).split()).upper()
    if actual != sequence:
        raise RuntimeError("Returned sequence differs from submitted sequence.")
    expected_md5 = hashlib.md5(sequence.encode("ascii")).hexdigest()
    if str(record.get("md5", "")).lower() != expected_md5:
        raise RuntimeError("Result MD5 does not match sequence.")
    matches = record.get("matches")
    if not isinstance(matches, list):
        raise RuntimeError("Missing matches list; not treated as no hits.")
    for match in matches:
        signature = match.get("signature", {})
        if not signature.get("accession") or not signature.get(
            "signatureLibraryRelease", {}
        ).get("library"):
            raise RuntimeError("Missing signature identity/library.")
        locations = match.get("locations")
        if not isinstance(locations, list) or not locations:
            raise RuntimeError("Signature lacks residue locations.")
        for location in locations:
            if not 1 <= int(location["start"]) <= int(location["end"]) <= len(sequence):
                raise RuntimeError("Invalid residue coordinates.")
    return record, version, release


def clean(value):
    return "-" if value is None or value == "" else " ".join(str(value).split())


def tsv_rows(record, identifiers, completed_utc):
    date = datetime.fromisoformat(completed_utc).strftime("%d-%m-%Y")
    for match in record["matches"]:
        signature = match["signature"]
        entry = signature.get("entry") or {}
        gos = sorted({x["id"] for x in entry.get("goXRefs", []) if x.get("id")})
        pathways = sorted({
            str(x.get("databaseName", "")) + ":" + str(x["id"])
            for x in entry.get("pathwayXRefs", []) if x.get("id")
        })
        for location in match["locations"]:
            score = next((v for v in (
                location.get("evalue"), match.get("evalue"),
                location.get("score"), match.get("score")
            ) if v is not None), "-")
            for identifier in identifiers:
                yield [identifier, record["md5"].lower(), len(record["sequence"]),
                       signature["signatureLibraryRelease"]["library"],
                       signature["accession"],
                       clean(signature.get("description") or signature.get("name")),
                       int(location["start"]), int(location["end"]), score, "T", date,
                       clean(entry.get("accession")),
                       clean(entry.get("description") or entry.get("name")),
                       "|".join(gos) or "-", "|".join(pathways) or "-"]


def seed_completed_cache(snapshot, cache_dir):
    """Import validated previously completed jobs, without making HTTP calls."""
    snapshot, cache_dir = Path(snapshot), Path(cache_dir)
    state = json.loads((snapshot / "STATE.json").read_text(encoding="utf-8"))
    source = Path(state["Input"])
    if sha(source) != state["InputSHA256"]:
        raise RuntimeError("Seed FASTA changed since the completed scan.")
    records, tasks = read_fasta(source)
    prepared = []
    for key, task in tasks.items():
        job = state["Jobs"].get(key, {})
        if job.get("Status") != "FINISHED" or set(job.get("IDs", [])) != set(task["ids"]):
            raise RuntimeError("Seed checkpoint does not cover input aliases.")
        result = json.loads((snapshot / (key + ".json")).read_text(encoding="utf-8"))
        record, version, release = validate_result(result, task["sequence"])
        if len(record["matches"]) != job.get("MatchCount"):
            raise RuntimeError("Seed checkpoint match count mismatch.")
        prepared.append((key, result, {
            "settings": SETTINGS, "status": "FINISHED", "jobId": job["JobID"],
            "completedUTC": job["UpdatedUTC"], "version": version, "release": release,
        }))
    with cache_lock(cache_dir):
        for key, result, job in prepared:
            result_path = cache_dir / (key + ".result.json")
            job_path = cache_dir / (key + ".job.json")
            if job_path.exists():
                old = json.loads(job_path.read_text(encoding="utf-8"))
                if old.get("status") != "FINISHED" or old.get("settings") != SETTINGS:
                    raise RuntimeError("Refusing to overwrite an existing pending/different cache job.")
            if result_path.exists():
                old_result = json.loads(result_path.read_text(encoding="utf-8"))
                if old_result != result:
                    raise RuntimeError("Existing cache evidence differs; not overwritten.")
            else:
                write_json(result_path, result)
            write_json(job_path, job)
    return source, len(records)


def run_interpro_v6(fasta_path, email, outdir, poll_seconds=30,
                    timeout_minutes=180, *, cache_dir=None, offline=False):
    """Drop-in automation stage. A partial/failed run never returns a TSV.

    email is retained for Gene Miner's CLI compatibility. The hosted MCP tool
    has no client email argument; this adapter does not send that email.
    """
    fasta_path, outdir = Path(fasta_path), Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    status_file = outdir / "interpro_web_status.txt"
    final_tsv = outdir / "interpro_web.tsv"
    # Prevent an old TSV from looking like this attempt's completed evidence.
    status_file.write_text("IN_PROGRESS\tInterProScan v6 adapter\n", encoding="utf-8")
    if not fasta_path.is_file():
        raise FileNotFoundError(str(fasta_path))
    records, tasks = read_fasta(fasta_path)
    if not records:
        status_file.write_text("NOT_RUN\tNo Review candidates.\n", encoding="utf-8")
        return None
    cache_dir = Path(cache_dir) if cache_dir else Path(__file__).resolve().parent / "interpro_v6_cache"
    deadline = time.monotonic() + max(1, float(timeout_minutes)) * 60
    manifest = {"inputSHA256": sha(fasta_path), "settings": SETTINGS,
                "endpoint": REST_ENDPOINT, "adapterVersion": ADAPTER_VERSION,
                "format": "adapter-generated 15-column TSV; raw v6 JSON authoritative",
                "complete": False, "proteins": []}
    rows, versions = [], set()
    client = None
    try:
        with cache_lock(cache_dir):
            for key, task in tasks.items():
                job_path = cache_dir / (key + ".job.json")
                result_path = cache_dir / (key + ".result.json")
                job = json.loads(job_path.read_text(encoding="utf-8")) if job_path.exists() else {
                    "settings": SETTINGS, "status": "NOT_SUBMITTED"}
                if job.get("settings") != SETTINGS:
                    raise RuntimeError("Cached job has different analysis settings.")
                print("[info] InterPro v6:", ", ".join(task["ids"]), flush=True)
                if job["status"] != "FINISHED":
                    if offline:
                        raise RuntimeError("Offline check requires completed cache: " + key)
                    if job["status"] in {"SUBMITTING", "SUBMISSION_UNCERTAIN"}:
                        raise RuntimeError("Submission uncertain; reconcile before retrying: " + str(job_path))
                    if job["status"] in {"ERROR", "FAILURE", "NOT_FOUND"}:
                        raise RuntimeError("Existing job requires inspection, not automatic resubmission: " + str(job_path))
                    if time.monotonic() >= deadline:
                        raise TimeoutError("Waiting limit reached; rerun to resume saved jobs.")
                    if client is None:
                        client = RestClient(email)
                    evidence_path = cache_dir / (key + ".last_rest_response.json")
                    if job.get("jobId"):
                        response = client.call("getSequenceAnalysis", {"jobId": job["jobId"]}, evidence_path)
                    else:
                        if job["status"] != "NOT_SUBMITTED":
                            raise RuntimeError("Unexpected cache state: " + str(job_path))
                        job["status"] = "SUBMITTING"
                        write_json(job_path, job)
                        try:
                            response = client.call("analyzeSequence", {
                                "sequence": task["sequence"], "goTerms": True, "pathways": True,
                            }, evidence_path)
                        except BaseException as exc:
                            job.update(status="SUBMISSION_UNCERTAIN", error=str(exc), updatedUTC=utc())
                            write_json(job_path, job)
                            raise
                    while True:
                        write_json(cache_dir / (key + ".response.json"), response)
                        if response.get("jobId"):
                            job["jobId"] = response["jobId"]
                        job.update(status=str(response.get("status", "UNKNOWN")).upper(), updatedUTC=utc())
                        if job["status"] != "FINISHED":
                            write_json(job_path, job)
                        if response.get("error"):
                            raise RuntimeError(str(response["error"]))
                        print("[info]", job["status"], job.get("jobId", ""), flush=True)
                        if job["status"] == "FINISHED":
                            result = response.get("result")
                            if isinstance(result, str):
                                result = json.loads(result)
                            _, version, release = validate_result(result, task["sequence"])
                            write_json(result_path, result)
                            job.update(completedUTC=utc(), version=version, release=release)
                            write_json(job_path, job)
                            break
                        if job["status"] not in {"QUEUED", "RUNNING", "PENDING"} or not job.get("jobId"):
                            raise RuntimeError("Unexpected v6 response; evidence saved in cache.")
                        remaining = deadline - time.monotonic()
                        if remaining <= 0:
                            raise TimeoutError("Waiting limit reached; saved jobs will resume on rerun.")
                        time.sleep(min(max(30, float(poll_seconds)), remaining))
                        response = client.call("getSequenceAnalysis", {"jobId": job["jobId"]}, evidence_path)
                result = json.loads(result_path.read_text(encoding="utf-8"))
                record, version, release = validate_result(result, task["sequence"])
                if (job.get("version"), job.get("release")) != (version, release):
                    raise RuntimeError("Cache metadata/version differs from raw result.")
                versions.add((version, release))
                if len(versions) > 1:
                    raise RuntimeError("Mixed v6/data releases; do not combine automatically. Inspect cache.")
                rows.extend(tsv_rows(record, task["ids"], job["completedUTC"]))
                manifest["proteins"].append({
                    "ids": task["ids"], "sha256": key, "jobId": job.get("jobId"),
                    "version": version, "release": release, "status": "FINISHED",
                    "signatureMatches": len(record["matches"]),
                    "resultFile": str(result_path), "resultSHA256": sha(result_path),
                })
                write_json(outdir / "interpro_v6_completion.json", manifest)
            temporary = final_tsv.with_suffix(".tsv.tmp")
            with temporary.open("w", encoding="utf-8", newline="") as handle:
                csv.writer(handle, delimiter="\t", lineterminator="\n").writerows(rows)
            os.replace(temporary, final_tsv)
            manifest.update(complete=True, proteinCount=len(records),
                            uniqueSequenceCount=len(tasks), locationRows=len(rows),
                            tsvSHA256=sha(final_tsv))
            write_json(outdir / "interpro_v6_completion.json", manifest)
            (outdir / "interpro_job_id.txt").write_text(
                "\n".join(x["jobId"] for x in manifest["proteins"]) + "\n", encoding="utf-8")
            status_file.write_text("FINISHED\tInterProScan v6\tinterpro_web.tsv\n", encoding="utf-8")
            print("[info] InterPro v6 complete:", len(records), "proteins;", len(rows), "location rows.")
            return str(final_tsv)
    except BaseException as exc:
        status_file.write_text("FAILED\t" + type(exc).__name__ + ": " + str(exc) + "\n", encoding="utf-8")
        manifest["error"] = str(exc)
        write_json(outdir / "interpro_v6_completion.json", manifest)
        raise
