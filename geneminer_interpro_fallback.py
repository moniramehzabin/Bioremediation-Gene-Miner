"""Optional InterPro orchestration; failures never erase normal core results."""
from pathlib import Path
from datetime import datetime
import hashlib
import json


def fasta_ids(path):
    with Path(path).open(encoding="utf-8") as handle:
        return {line[1:].strip().split()[0] for line in handle if line.startswith(">")}


def save_status(outdir, state):
    path = Path(outdir) / "interpro_enrichment_status.json"
    temporary = path.with_suffix(".json.tmp")
    temporary.write_text(json.dumps(state, indent=2), encoding="utf-8")
    temporary.replace(path)


def status_sheet(report, state, pd):
    values = [["Core_results", "AVAILABLE"], ["InterPro_stage", state["Stage"]],
              ["InterPro_evidence_attached", state.get("EvidenceAttached", False)],
              ["InterPro_selected_proteins", state.get("SelectedReviewProteins", "")],
              ["InterPro_completed_proteins", state.get("CompletedReviewProteins", "")],
              ["InterPro_note", state.get("Note", "")], ["InterPro_error", state.get("Error", "")]]
    for key, value in state.get("ProteinStatusCounts", {}).items():
        values.append(["Protein_status_" + key, value])
    values += [["InterPro_resume", "Rerun the same command/output directory. Completed cache is reused; terminal/uncertain failures are not automatically resubmitted."],
               ["Decision_policy", "Candidate/Review/Reject and sequence-match confidence are unchanged."]]
    with pd.ExcelWriter(report, engine="openpyxl", mode="a", if_sheet_exists="replace") as writer:
        pd.DataFrame(values, columns=["Item", "Value"]).to_excel(writer, sheet_name="Run_status", index=False)


def mark_unresolved(frame, selected, stage, note, completed=None, protein_states=None):
    if frame is None or frame.empty:
        return frame
    frame = frame.copy()
    before = frame[["decision", "confidence"]].copy()
    completed = set(completed or [])
    protein_states = protein_states or {}
    unresolved = set(selected) - completed
    if stage in {"PENDING", "INCOMPLETE", "PARTIAL"}:
        for protein in unresolved:
            value = protein_states.get(protein, "PENDING_ENRICHMENT")
            status = {
                "FAILURE": "FAILED", "ERROR": "FAILED", "NOT_FOUND": "FAILED",
                "SUBMITTING": "SUBMISSION_UNCERTAIN",
                "SUBMISSION_UNCERTAIN": "SUBMISSION_UNCERTAIN",
                "NOT_SUBMITTED": "NOT_SUBMITTED", "RUNNING": "RUNNING", "QUEUED": "QUEUED",
                "RESULT_INVALID": "RESULT_INVALID", "VERSION_MISMATCH": "VERSION_MISMATCH",
                "SETTINGS_MISMATCH": "SETTINGS_MISMATCH",
            }.get(value, "PENDING_ENRICHMENT")
            mask = frame["query"].astype(str).eq(protein)
            frame.loc[mask, "interpro_status"] = status
            frame.loc[mask, "resolver_status"] = "PENDING_INTERPRO"
            frame.loc[mask, "resolver_note"] = "InterPro enrichment unresolved for this protein. " + note
    elif stage == "NOT_REQUESTED":
        mask = frame["query"].astype(str).isin(selected)
        frame.loc[mask, "interpro_status"] = "NOT_REQUESTED"
        frame.loc[mask, "resolver_status"] = ""
        frame.loc[mask, "resolver_note"] = note
    if not frame[["decision", "confidence"]].equals(before):
        raise RuntimeError("InterPro status annotation changed core decisions.")
    return frame


def pending_state(args, selected):
    requested = bool(args.interpro or args.interpro_auto)
    if not selected:
        stage, note = "NOT_NEEDED", "No selected Review sequences require InterPro enrichment."
    elif requested:
        stage, note = "PENDING", "Core results are available. Optional InterPro enrichment has not finished."
    else:
        stage, note = "NOT_REQUESTED", "Core-only run. InterPro was not requested; no negative domain conclusion is made."
    return {"Stage": stage, "EvidenceAttached": False,
            "SelectedReviewProteins": len(selected), "Note": note}


def optional_interpro(args, fasta, outdir, runner):
    import csv
    selected = fasta_ids(fasta)
    state = pending_state(args, selected)
    if args.interpro:
        state.update(Stage="IMPORTED_HITS_ONLY", EvidenceAttached=True,
                     Note="Supplied hit-only TSV imported. Missing rows do not prove completed no-hit scans.")
        save_status(outdir, state)
        return args.interpro, state, set()
    if not args.interpro_auto or not selected:
        save_status(outdir, state)
        return None, state, set()
    try:
        path = runner(fasta, args.interpro_email, outdir,
                      args.interpro_poll_seconds, args.interpro_timeout_minutes)
        if not path or not Path(path).is_file():
            raise RuntimeError("No InterPro batch evidence file returned.")
        completion = json.loads((Path(outdir) / "interpro_v6_completion.json").read_text(encoding="utf-8"))
        if completion.get("inputSHA256") != hashlib.sha256(Path(fasta).read_bytes()).hexdigest():
            raise RuntimeError("Batch manifest refers to different input.")
        items = completion.get("proteins", [])
        all_ids = [identifier for item in items for identifier in item.get("ids", [])]
        if len(all_ids) != len(set(all_ids)) or set(all_ids) != selected:
            raise RuntimeError("Batch manifest does not account for input aliases exactly once.")
        completed = {identifier for item in items if item.get("status") == "FINISHED"
                     for identifier in item.get("ids", [])}
        if bool(completion.get("complete")) != (completed == selected):
            raise RuntimeError("Inconsistent completion flag.")
        if completion.get("tsvSHA256") != hashlib.sha256(Path(path).read_bytes()).hexdigest():
            raise RuntimeError("Batch TSV checksum mismatch.")
        with Path(path).open(encoding="utf-8", newline="") as handle:
            rows = list(csv.reader(handle, delimiter="\t"))
        if len(rows) != completion.get("locationRows"):
            raise RuntimeError("Batch row count mismatch.")
        if any(len(row) != 15 or row[0] not in completed for row in rows):
            raise RuntimeError("TSV contains uncompleted/unexpected protein evidence.")
        for item in items:
            if item.get("status") == "FINISHED":
                raw = Path(item["resultFile"])
                if hashlib.sha256(raw.read_bytes()).hexdigest() != item["resultSHA256"]:
                    raise RuntimeError("Completed raw result changed.")
        stage = "COMPLETE" if completed == selected else "PARTIAL" if completed else "INCOMPLETE"
        state.update(Stage=stage, EvidenceAttached=bool(completed),
                     CompletedReviewProteins=len(completed),
                     ProteinStates={identifier: item["status"] for item in items for identifier in item["ids"]},
                     ProteinStatusCounts=completion.get("proteinStatusCounts", {}),
                     Note="Only verified completed proteins have attached InterPro evidence. Failed, uncertain and unsubmitted proteins remain separate; core decisions are unchanged.")
        save_status(outdir, state)
        return str(path), state, completed
    except Exception as exc:
        state.update(Stage="INCOMPLETE", EvidenceAttached=False, Error=type(exc).__name__ + ": " + str(exc),
                     Note="Core report remains available. InterPro evidence failed validation or the service was unavailable; not labelled no hits.")
        save_status(outdir, state)
        print("[warn] Optional InterPro enrichment incomplete:", exc)
        return None, state, set()


def preserve_old_full_report(outdir, stage):
    if stage not in {"COMPLETE", "IMPORTED_HITS_ONLY"}:
        old = Path(outdir) / "InterPro_Full_Report.xlsx"
        if old.exists():
            archived = old.with_name("InterPro_Full_Report.previous_" + datetime.now().strftime("%Y%m%d_%H%M%S_%f") + ".xlsx")
            old.rename(archived)
