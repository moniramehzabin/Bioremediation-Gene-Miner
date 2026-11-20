"""Foreground recovery around the installed, unmodified InterProScan adapter.

Only the original adapter submits, polls, validates results, and writes its cache.
This supervisor retries failed reads of saved jobs and isolates individual errors.
"""

from contextlib import contextmanager
from datetime import datetime, timezone
import argparse
import hashlib
import http.client
import importlib.util
import json
import os
from pathlib import Path
import socket
import sys
import time
import urllib.error


PENDING = {"QUEUED", "RUNNING", "PENDING"}
BLOCKED = {"SUBMITTING", "SUBMISSION_UNCERTAIN", "ERROR", "FAILURE", "NOT_FOUND"}
READ_HTTP_CODES = {400, 408, 429, 500, 502, 503, 504}


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_json(path):
    return json.loads(Path(path).read_text(encoding="utf-8-sig"))


def write_json(path, value):
    path = Path(path)
    temporary = path.with_name(path.name + ".tmp")
    temporary.write_text(json.dumps(value, indent=2), encoding="utf-8")
    os.replace(temporary, path)


def load_original(repo, filename):
    path = Path(repo) / filename
    spec = importlib.util.spec_from_file_location("original_" + path.stem, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Cannot load the installed adapter: {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@contextmanager
def session_lock(cache):
    """OS-owned lock: closing this process releases it, including on force exit."""
    path = Path(cache) / "launcher_recovery.lock"
    path.parent.mkdir(parents=True, exist_ok=True)
    handle = path.open("a+b")
    acquired = False
    try:
        if path.stat().st_size == 0:
            handle.write(b"0")
            handle.flush()
        handle.seek(0)
        try:
            if os.name == "nt":
                import msvcrt
                msvcrt.locking(handle.fileno(), msvcrt.LK_NBLCK, 1)
            else:
                import fcntl
                fcntl.flock(handle.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
        except OSError as error:
            raise RuntimeError("An InterPro recovery session is already open for this cache.") from error
        acquired = True
        yield
    finally:
        if acquired:
            handle.seek(0)
            if os.name == "nt":
                import msvcrt
                msvcrt.locking(handle.fileno(), msvcrt.LK_UNLCK, 1)
            else:
                import fcntl
                fcntl.flock(handle.fileno(), fcntl.LOCK_UN)
        handle.close()


def process_is_gone(pid, created):
    """True only when the lock owner is demonstrably gone on this computer."""
    if not isinstance(pid, int) or pid <= 0:
        return False
    if os.name != "nt":
        try:
            os.kill(pid, 0)
        except ProcessLookupError:
            return True
        except OSError:
            pass
        return False

    import ctypes
    from ctypes import wintypes
    kernel = ctypes.WinDLL("kernel32", use_last_error=True)
    kernel.OpenProcess.argtypes = [wintypes.DWORD, wintypes.BOOL, wintypes.DWORD]
    kernel.OpenProcess.restype = wintypes.HANDLE
    kernel.GetProcessTimes.argtypes = [wintypes.HANDLE] + [ctypes.POINTER(wintypes.FILETIME)] * 4
    kernel.GetProcessTimes.restype = wintypes.BOOL
    kernel.CloseHandle.argtypes = [wintypes.HANDLE]
    kernel.CloseHandle.restype = wintypes.BOOL
    handle = kernel.OpenProcess(0x1000, False, pid)
    if not handle:
        return ctypes.get_last_error() == 87  # ERROR_INVALID_PARAMETER: no such PID.
    try:
        times = [wintypes.FILETIME() for _ in range(4)]
        if not kernel.GetProcessTimes(handle, *(ctypes.byref(item) for item in times)):
            return False
        started = ((times[0].dwHighDateTime << 32) | times[0].dwLowDateTime) / 10000000 - 11644473600
        locked = datetime.fromisoformat(created.replace("Z", "+00:00")).timestamp()
        return started > locked + 2  # PID has been reused by a newer process.
    finally:
        kernel.CloseHandle(handle)


def check_adapter_lock(directory):
    """Archive a proven stale original lock; retain active or ambiguous locks."""
    path = Path(directory) / "ACTIVE_RUN.lock"
    if not path.exists():
        return
    before = path.read_bytes()
    try:
        owner = json.loads(before)
        gone = process_is_gone(owner.get("pid"), owner.get("createdUTC", ""))
    except (ValueError, TypeError, OSError):
        gone = False
    if not gone or path.read_bytes() != before:
        raise RuntimeError(f"Another run may still own this lock. Close that run before resuming: {path}")
    archived = path.with_name(f"ACTIVE_RUN.stale-{time.time_ns()}.json")
    os.replace(path, archived)
    print(f"[recovery] Previous process has exited; archived its lock: {archived}", flush=True)


def bind_input(repo, run, cache, adapter):
    fasta = run / "interpro_candidates.faa"
    if not fasta.is_file():
        raise FileNotFoundError(f"Select the results folder containing interpro_candidates.faa: {run}")
    expected = sha(fasta)
    state_path = run / "InterPro_resume" / "launcher_recovery_state.json"
    evidence = [
        run / "interpro_v6_completion.json", run / "interpro_v6_batch_progress.json",
        run / "InterPro_resume" / "interpro_v6_completion.json",
        run / "InterPro_resume" / "interpro_v6_batch_progress.json", state_path,
    ]
    found = False
    for path in evidence:
        if path.is_file():
            data = read_json(path)
            if data.get("inputSHA256"):
                found = True
                if data["inputSHA256"] != expected:
                    raise RuntimeError(f"The selected protein FASTA differs from its saved run: {path}")
            if data.get("settings") is not None and data["settings"] != adapter.SETTINGS:
                raise RuntimeError(f"Saved analysis settings differ from the installed adapter: {path}")
    if not found:
        raise RuntimeError("No saved InterPro input fingerprint was found. Select the original results folder.")
    if state_path.is_file():
        saved = read_json(state_path)
        if saved.get("repo") != str(repo) or saved.get("cache") != str(cache):
            raise RuntimeError("This recovery belongs to a different tool or cache folder.")
    return fasta, expected, state_path


def job_state(cache, key, settings):
    path = cache / (key + ".job.json")
    return read_json(path) if path.is_file() else {"status": "NOT_SUBMITTED", "settings": settings}


def log_error(path, key, task, error, before, after):
    item = {
        "timeUTC": datetime.now(timezone.utc).isoformat(), "ids": task["ids"], "sha256": key,
        "jobId": after.get("jobId") or before.get("jobId"), "cachedStatus": after.get("status"),
        "errorType": type(error).__name__, "message": str(error),
    }
    if isinstance(error, urllib.error.HTTPError):
        item.update(httpStatus=error.code, url=error.url)
        try:
            item["responseBody"] = error.read(8192).decode("utf-8", errors="replace")
        except OSError:
            item["responseBody"] = "Response body could not be read."
        finally:
            error.close()
    with path.open("a", encoding="utf-8") as handle:
        handle.write(json.dumps(item) + "\n")


def read_can_retry(error, before, after):
    identifier = after.get("jobId")
    if not identifier or after.get("status") not in PENDING:
        return False
    if before.get("jobId") and before["jobId"] != identifier:
        return False
    if isinstance(error, urllib.error.HTTPError):
        return error.code in READ_HTTP_CODES
    return isinstance(error, (urllib.error.URLError, TimeoutError, ConnectionError, socket.timeout,
                             http.client.IncompleteRead, json.JSONDecodeError))


def recover(repo, run, email, *, max_read_retries=5, max_rounds=3, job_minutes=15,
            total_minutes=720, pause_seconds=30, sleep=time.sleep, clock=time.monotonic,
            adapter=None, batch=None):
    """Run in the calling console; leave all internal source files unchanged."""
    if min(max_read_retries, max_rounds, job_minutes, total_minutes, pause_seconds) <= 0:
        raise ValueError("Recovery limits must be greater than zero.")
    repo, run = Path(repo).resolve(), Path(run).resolve()
    cache = repo / "interpro_v6_cache"
    adapter = adapter or load_original(repo, "geneminer_interpro_v6.py")
    batch = batch or load_original(repo, "geneminer_interpro_batch.py")
    if not hasattr(adapter, "REST_ENDPOINT"):
        raise RuntimeError("Select the original repository with the installed InterProScan v6 REST adapter.")
    fasta, input_hash, state_path = bind_input(repo, run, cache, adapter)
    records, tasks = adapter.read_fasta(fasta)
    out = state_path.parent
    out.mkdir(exist_ok=True)
    diagnostic = out / "launcher_interpro_errors.jsonl"
    source_hashes = {name: sha(repo / name) for name in (
        "geneminer_interpro_v6.py", "geneminer_interpro_batch.py")}
    state = {"inputSHA256": input_hash, "repo": str(repo), "cache": str(cache),
             "settings": adapter.SETTINGS, "originalCodeSHA256": source_hashes,
             "status": "RUNNING", "pid": os.getpid()}
    deadline = clock() + total_minutes * 60
    blocked, failures = set(), {}

    def verify():
        if sha(fasta) != input_hash:
            raise RuntimeError("The saved protein input changed during recovery. Stopped before further requests.")
        if any(sha(repo / name) != digest for name, digest in source_hashes.items()):
            raise RuntimeError("An original adapter changed during recovery. Close other updates and resume again.")

    def aggregate():
        verify()
        check_adapter_lock(out / "interpro_batch_lock")
        batch.run_interpro_batch(fasta, email, out, 30, 1, offline=True, cache_dir=cache)
        if not records:
            return {"complete": True, "proteinCount": 0, "completedProteinCount": 0}
        return read_json(out / "interpro_v6_completion.json")

    print(f"\nInterPro recovery: {len(records)} proteins. This console must stay open.", flush=True)
    print("Existing job IDs and completed results are reused. Ctrl+C stops this session.", flush=True)
    with session_lock(cache):
        check_adapter_lock(cache)
        write_json(state_path, state)
        try:
            manifest = aggregate()
            for round_number in range(1, max_rounds + 1):
                if manifest["complete"] or clock() >= deadline:
                    break
                attempted = False
                print(f"\n[recovery] Pass {round_number}/{max_rounds}", flush=True)
                for key, task in tasks.items():
                    verify()
                    job = job_state(cache, key, adapter.SETTINGS)
                    if job.get("settings") != adapter.SETTINGS or job.get("status") in BLOCKED:
                        blocked.add(key)
                    if key in blocked or job.get("status") == "FINISHED" or clock() >= deadline:
                        continue
                    attempted = True
                    per_dir = out / "launcher_proteins" / key
                    per_dir.mkdir(parents=True, exist_ok=True)
                    per_fasta = per_dir / "query.faa"
                    per_fasta.write_text("".join(
                        f">{identifier}\n{task['sequence']}\n" for identifier in task["ids"]
                    ), encoding="ascii")
                    while clock() < deadline:
                        verify()
                        check_adapter_lock(cache)
                        before = job_state(cache, key, adapter.SETTINGS)
                        print("[recovery] Protein:", ", ".join(task["ids"]),
                              "saved job:", before.get("jobId") or "not yet submitted", flush=True)
                        try:
                            adapter.run_interpro_v6(
                                per_fasta, email, per_dir, 30,
                                min(job_minutes, max(0.01, (deadline - clock()) / 60)), cache_dir=cache,
                            )
                        except Exception as error:
                            after = job_state(cache, key, adapter.SETTINGS)
                            log_error(diagnostic, key, task, error, before, after)
                            print(f"[recovery] {type(error).__name__}: {error}", flush=True)
                            if read_can_retry(error, before, after):
                                failures[key] = failures.get(key, 0) + 1
                                if isinstance(error, TimeoutError) and "Waiting limit reached" in str(error):
                                    print("[recovery] Still waiting. Other proteins continue before checking this job again.", flush=True)
                                    break
                                if failures[key] <= max_read_retries:
                                    delay = min(pause_seconds * 2 ** (failures[key] - 1), 120,
                                                max(0, deadline - clock()))
                                    print(f"[recovery] Retrying the same job in {delay:g} seconds "
                                          f"({failures[key]}/{max_read_retries}).", flush=True)
                                    sleep(delay)
                                    continue
                            blocked.add(key)
                            print("[recovery] This protein is deferred. Continuing with the other proteins.", flush=True)
                        break
                manifest = aggregate()
                print(f"[recovery] Completed {manifest['completedProteinCount']} of {manifest['proteinCount']}", flush=True)
                if manifest["complete"] or not attempted:
                    break
                pending = [key for key in tasks if key not in blocked and
                           job_state(cache, key, adapter.SETTINGS).get("status") != "FINISHED"]
                if not pending:
                    break
                if round_number < max_rounds and clock() < deadline:
                    delay = min(pause_seconds, max(0, deadline - clock()))
                    print(f"[recovery] Next pass in {delay:g} seconds.", flush=True)
                    sleep(delay)
            state.update(status="COMPLETE" if manifest["complete"] else "PARTIAL",
                         completedProteinCount=manifest["completedProteinCount"], proteinCount=len(records))
            write_json(state_path, state)
            if not manifest["complete"]:
                print("[recovery] Unfinished jobs remain saved. Open Resume InterPro.cmd to continue later.", flush=True)
                for item in manifest["proteins"]:
                    if item["status"] != "FINISHED":
                        print("[recovery]", ",".join(item["ids"]), item["status"], item.get("note", ""), flush=True)
                print(f"[recovery] Error details: {diagnostic}", flush=True)
            return manifest
        except BaseException as error:
            state.update(status="INTERRUPTED" if isinstance(error, KeyboardInterrupt) else "STOPPED",
                         message=f"{type(error).__name__}: {error}")
            write_json(state_path, state)
            raise


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", required=True, type=Path)
    parser.add_argument("--run", required=True, type=Path)
    parser.add_argument("--email", required=True)
    args = parser.parse_args()
    if not args.email.strip():
        parser.error("An InterProScan email address is required.")
    try:
        result = recover(args.repo, args.run, args.email.strip())
    except KeyboardInterrupt:
        print("\nStopped. Saved jobs and completed results are retained for Resume InterPro.cmd.", flush=True)
        return 130
    except Exception as error:
        print(f"\nInterPro recovery stopped: {error}", flush=True)
        return 1
    if result["complete"]:
        print(f"\nInterPro complete: {result['completedProteinCount']} of {result['proteinCount']} proteins.", flush=True)
        return 0
    return 2


if __name__ == "__main__":
    sys.exit(main())
