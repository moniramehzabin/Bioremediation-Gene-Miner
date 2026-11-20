"""Select inputs and run the existing Gene Miner command-line workflow."""

from datetime import datetime
import importlib.util
import json
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import tkinter as tk
from tkinter import filedialog, messagebox, simpledialog

from interpro_recovery import sha, write_json


TITLE = "Bioremediation Gene Miner"
HERE = Path(__file__).resolve().parent
REFERENCE_FOLDER = Path("reference_updates") / "2026-11-18_413_sequence_reference"


def run_command(command, working_directory, *, check=True):
    print("\n" + subprocess.list2cmdline([str(item) for item in command]), flush=True)
    result = subprocess.run([str(item) for item in command], cwd=working_directory, check=check)
    if result.returncode in {130, -2, -1073741510, 3221225786}:
        raise KeyboardInterrupt
    return result.returncode


def run_recovery(tool, output, email):
    code = run_command([
        sys.executable, HERE / "interpro_recovery.py", "--repo", tool.parent,
        "--run", output, "--email", email,
    ], tool.parent, check=False)
    if code not in {0, 2}:
        raise RuntimeError("InterPro recovery stopped. The reason is shown in this console.")
    return code == 0


def fingerprint_inputs(tool, reference, diamond, genome):
    paths = [tool, diamond, genome]
    paths += [tool.parent / name for name in (
        "geneminer_interpro_v6.py", "geneminer_interpro_batch.py", "geneminer_interpro_fallback.py")]
    paths += [reference / name for name in ("reference.faa", "reference_metadata.tsv", "annotation_rules.tsv")]
    return {str(path.resolve()): sha(path) for path in paths}


def verify_inputs(fingerprints):
    for name, digest in fingerprints.items():
        if sha(name) != digest:
            raise RuntimeError(f"A saved run input changed. Stopped before reusing alignments: {name}")


def run_analysis(tool, reference, diamond, genome, output, email):
    """Call DIAMOND and Gene Miner through their existing command-line interfaces."""
    fingerprints = fingerprint_inputs(tool, reference, diamond, genome)
    write_json(output / "launcher_run.json", {
        "genome": str(genome), "repo": str(tool.parent), "inputSHA256": fingerprints,
        "created": datetime.now().isoformat(),
    })
    database = output / "reference_413"
    run_command([
        diamond, "makedb", "--in", reference / "reference.faa", "--db", database,
    ], tool.parent)
    command = [
        sys.executable, tool, genome,
        "--reference-db", database,
        "--reference-metadata", reference / "reference_metadata.tsv",
        "--annotation-rules", reference / "annotation_rules.tsv",
        "--diamond", diamond,
        "--outdir", output,
    ]
    if email:
        command.extend(["--interpro-auto", "--interpro-email", email])
    run_command(command, tool.parent)
    if not email:
        return True
    state_file = output / "interpro_enrichment_status.json"
    state = json.loads(state_file.read_text(encoding="utf-8"))
    if state.get("Stage") in {"COMPLETE", "NOT_NEEDED"}:
        return True
    print("\n[launcher] InterPro paused. Continuing saved jobs through foreground recovery.", flush=True)
    if not run_recovery(tool, output, email):
        return False
    # Only after all scans finish, ask the original CLI to read its completed cache.
    # The exact genome and other inputs are checked before any alignment reuse.
    verify_inputs(fingerprints)
    run_command(command + ["--skip-diamond"], tool.parent)
    state = json.loads(state_file.read_text(encoding="utf-8"))
    if state.get("Stage") not in {"COMPLETE", "NOT_NEEDED"}:
        raise RuntimeError("Scans completed, but the original tool did not confirm complete enrichment. See the console.")
    return True


def ask_email(window):
    value = simpledialog.askstring(TITLE, "Enter your email for InterProScan.", parent=window)
    if value is None:
        return None
    value = value.strip()
    if not value:
        raise ValueError("An email address is required for InterProScan.")
    return value


def ensure_dependencies(tool):
    missing = [module for module in ("Bio", "pandas", "openpyxl")
               if importlib.util.find_spec(module) is None]
    if missing:
        run_command([sys.executable, "-m", "pip", "install", "-r", tool.parent / "requirements.txt"],
                    tool.parent)


def resume_interpro(window):
    tool = find_tool(window)
    if tool is None:
        return
    selected = filedialog.askdirectory(
        parent=window, title="Select the original Gene Miner results folder",
        initialdir=str(Path.home() / "Gene Miner Results"),
    )
    if not selected:
        return
    email = ask_email(window)
    if email is None:
        return
    ensure_dependencies(tool)
    complete = run_recovery(tool, Path(selected).resolve(), email)
    if complete:
        messagebox.showinfo(TITLE, "All selected InterPro scans are complete.\n\n"
                            f"Results: {Path(selected) / 'InterPro_resume'}", parent=window)
    else:
        messagebox.showwarning(TITLE, "Some InterPro jobs still need attention. Their job IDs and results are saved.\n\n"
                               "See the unfinished protein names and error details in this console.", parent=window)


def find_tool(window):
    for folder in (HERE, HERE.parent, Path(r"F:\GeneMiner_repo_update_20261031_185718")):
        tool = folder / "bioremediation_gene_miner.py"
        if tool.is_file():
            return tool.resolve()
    selected = filedialog.askopenfilename(
        parent=window, title="Select your existing bioremediation_gene_miner.py",
        filetypes=[("Gene Miner Python script", "bioremediation_gene_miner.py")],
    )
    return Path(selected).resolve() if selected else None


def find_diamond(window, tool):
    candidates = [tool.parent / "tools" / "diamond.exe",
                  Path(r"F:\Tools\EBPT_USER_TEST\diamond-windows\diamond.exe")]
    installed = shutil.which("diamond")
    if installed:
        candidates.append(Path(installed))
    for candidate in candidates:
        if candidate.is_file():
            return candidate.resolve()
    selected = filedialog.askopenfilename(
        parent=window, title="Select your DIAMOND executable",
        filetypes=[("DIAMOND executable", "diamond.exe")],
    )
    return Path(selected).resolve() if selected else None


def main():
    window = tk.Tk()
    window.withdraw()
    try:
        if "--resume-interpro" in sys.argv[1:]:
            resume_interpro(window)
            return
        tool = find_tool(window)
        if tool is None:
            return
        reference = tool.parent / REFERENCE_FOLDER
        for name in ("reference.faa", "reference_metadata.tsv", "annotation_rules.tsv"):
            if not (reference / name).is_file():
                raise FileNotFoundError(f"The 413-sequence reference file is missing.\n{reference / name}")
        diamond = find_diamond(window, tool)
        if diamond is None:
            return
        genome_file = filedialog.askopenfilename(
            parent=window, title="Select your annotated genome",
            filetypes=[("GenBank files", ("*.gbk", "*.gbff", "*.gb", "*.genbank")),
                       ("All files", "*")],
        )
        if not genome_file:
            return
        use_interpro = messagebox.askyesnocancel(
            TITLE, "Include InterProScan v6?\n\nSelected protein sequences will be sent to EMBL-EBI.",
            parent=window, default="yes",
        )
        if use_interpro is None:
            return
        email = ""
        if use_interpro:
            email = ask_email(window)
            if email is None:
                return
        ensure_dependencies(tool)

        results = Path.home() / "Gene Miner Results"
        results.mkdir(exist_ok=True)
        stamp = datetime.now().strftime("%Y-%m-%d_%H%M%S")
        output = Path(tempfile.mkdtemp(prefix=f"Gene Miner {stamp} ", dir=results))
        print(f"\nUsing Gene Miner at {tool}\nResults will be saved in {output}", flush=True)
        complete = run_analysis(tool, reference, diamond, Path(genome_file).resolve(), output, email)
        status_file = output / "interpro_enrichment_status.json"
        if status_file.is_file():
            state = json.loads(status_file.read_text(encoding="utf-8"))
            if use_interpro and state.get("Stage") not in {"COMPLETE", "NOT_NEEDED"}:
                detail = state.get("Error") or state.get("Note") or "See the run log in this window."
                messagebox.showwarning(
                    TITLE, f"InterProScan remains incomplete.\n\n{detail}\n\n"
                    f"Saved jobs can be continued with Resume InterPro.cmd.\n\nResults\n{output}", parent=window,
                )
        if complete and messagebox.askyesno(TITLE, f"Gene Miner has finished.\n\nOpen results?\n{output}", parent=window):
            import os

            os.startfile(output)
    except KeyboardInterrupt:
        print("\nStopped. Open Resume InterPro.cmd to continue saved InterPro jobs later.", flush=True)
    except (OSError, subprocess.CalledProcessError, ValueError, RuntimeError) as error:
        print(f"\nGene Miner stopped. {error}", flush=True)
        messagebox.showerror(TITLE, f"{error}\n\nThe command output is in the open console window.", parent=window)
    finally:
        window.destroy()


if __name__ == "__main__":
    main()
