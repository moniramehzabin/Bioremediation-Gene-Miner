#!/usr/bin/env python3
"""Prepare the saved reference snapshot. Does not search online."""
from pathlib import Path
import csv
import re
import os
import tempfile

ROOT = Path(__file__).resolve().parent
SOURCE = ROOT / "reference_updates" / "2026-10-31_audit_in_progress"
OUT_FASTA = ROOT / "bioremediation_reference_v0.2.1.faa"
OUT_META = ROOT / "bioremediation_reference_v0.2.1_metadata.tsv"

def main():
    source_fasta = SOURCE / "reference.faa"
    source_meta = SOURCE / "reference_metadata.tsv"

    fasta_bytes = source_fasta.read_bytes()
    meta_bytes = source_meta.read_bytes()

    records = {}
    current = None
    for line in fasta_bytes.decode("utf-8-sig").splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            header = line[1:].split()
            if not header:
                raise ValueError("Empty FASTA header")
            current = header[0]
            if current in records:
                raise ValueError(f"Duplicate FASTA ID: {current}")
            records[current] = ""
        else:
            if current is None:
                raise ValueError("Sequence before first FASTA header")
            records[current] += line

    with source_meta.open(encoding="utf-8-sig", newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        required = {
            "Accession", "Length_aa", "Headline_category", "Category_flag",
            "All_family_targets", "Evidence_classes", "Default_confidences",
            "Interpretation_rules", "Source_queries"
        }
        if not required.issubset(reader.fieldnames or []):
            raise ValueError("Required metadata columns are missing")
        rows = list(reader)

    ids = [row["Accession"] for row in rows]
    if len(ids) != len(set(ids)):
        raise ValueError("Duplicate metadata accessions")
    if set(ids) != set(records):
        raise ValueError("FASTA and metadata accessions differ")
    if len(records) != 294:
        raise ValueError("Expected the approved 294-reference snapshot")

    for row in rows:
        acc = row["Accession"]
        sequence = records[acc]
        if not re.fullmatch(r"[ACDEFGHIKLMNPQRSTVWY]+", sequence):
            raise ValueError(f"Invalid protein sequence: {acc}")
        if len(sequence) != int(row["Length_aa"]):
            raise ValueError(f"Sequence length mismatch: {acc}")

    outputs = [(OUT_FASTA, fasta_bytes), (OUT_META, meta_bytes)]
    for destination, data in outputs:
        if destination.exists() and destination.read_bytes() != data:
            raise FileExistsError(
                f"Existing output differs: {destination}. "
                "Preserve or move it before rebuilding."
            )

    temporary = []
    try:
        for destination, data in outputs:
            if destination.exists():
                continue
            with tempfile.NamedTemporaryFile(
                dir=ROOT, prefix=".reference_build_", delete=False
            ) as f:
                temporary.append((Path(f.name), destination))
                f.write(data)
        for temp, destination in temporary:
            os.replace(temp, destination)
    finally:
        for temp, _ in temporary:
            if temp.exists():
                temp.unlink()

    print("Prepared 294 saved reference proteins and matching metadata.")
    print("FASTA:", OUT_FASTA)
    print("Metadata:", OUT_META)
    print("Snapshot status: audit in progress.")
    print("No online searches. Evidence and confidence values unchanged.")
    print("DIAMOND database not built; genome analysis not run.")

if __name__ == "__main__":
    main()
