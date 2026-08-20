#!/usr/bin/env python3
"""
Bioremediation Gene Miner v0.3.2

Unified bacterial WGS screening prototype:
1) screens already-annotated CDSs
2) extracts hypothetical proteins
3) runs DIAMOND against the curated reference DB
4) applies identity/query-coverage/reference-coverage filters
5) optionally integrates InterProScan TSV
6) writes one evidence-ranked Excel/TSV report

This software predicts candidates. It does not experimentally prove
bioremediation activity.
"""
from __future__ import annotations

import argparse
from copy import copy
import csv
import re
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

import pandas as pd

try:
    from Bio import SeqIO
except Exception:
    SeqIO = None

VERSION = "0.3.2"

# ----------------------------
# GenBank parsing
# ----------------------------
def qval(q, key):
    vals = q.get(key, [""])
    return vals[0] if vals else ""

def parse_genbank_biopython(path):
    if SeqIO is None:
        raise RuntimeError("Biopython unavailable")
    cds = []
    for rec in SeqIO.parse(str(path), "genbank"):
        for feat in rec.features:
            if feat.type != "CDS":
                continue
            q = feat.qualifiers
            cds.append({
                "contig": rec.id,
                "locus_tag": qval(q, "locus_tag"),
                "gene": qval(q, "gene"),
                "product": qval(q, "product"),
                "translation": qval(q, "translation").replace(" ", ""),
                "EC_number": ";".join(q.get("EC_number", [])),
                "db_xref": ";".join(q.get("db_xref", [])),
            })
    if not cds:
        raise ValueError("No CDS parsed")
    return cds

def _qualifier(block, key):
    m = re.search(rf'/{re.escape(key)}="([^"]*)"', block, flags=re.S)
    return re.sub(r"\s+", " ", m.group(1)).strip() if m else ""

def parse_genbank_fallback(path):
    # Tolerant parser for some Prokka GenBank files with nonstandard LOCUS spacing.
    text = Path(path).read_text(encoding="utf-8", errors="replace")
    contig_starts = list(re.finditer(r"^LOCUS\s+(\S+)", text, flags=re.M))
    contig_positions = [(m.start(), m.group(1)) for m in contig_starts]
    cds_starts = list(re.finditer(r"^     CDS\s+(.+)$", text, flags=re.M))
    cds = []
    for i, m in enumerate(cds_starts):
        start = m.start()
        end = cds_starts[i+1].start() if i+1 < len(cds_starts) else len(text)
        block = text[start:end]
        contig = ""
        for pos, name in contig_positions:
            if pos <= start:
                contig = name
            else:
                break
        cds.append({
            "contig": contig,
            "locus_tag": _qualifier(block, "locus_tag"),
            "gene": _qualifier(block, "gene"),
            "product": _qualifier(block, "product"),
            "translation": _qualifier(block, "translation").replace(" ", ""),
            "EC_number": _qualifier(block, "EC_number"),
            "db_xref": _qualifier(block, "db_xref"),
        })
    cds = [r for r in cds if r["locus_tag"]]
    if not cds:
        raise ValueError("Fallback parser found no CDS")
    return cds

def parse_genbank(path):
    try:
        return parse_genbank_biopython(path)
    except Exception as e:
        print(f"[info] Biopython GenBank parser failed ({e}); using tolerant fallback parser.")
        return parse_genbank_fallback(path)

# ----------------------------
# Annotation screen
# ----------------------------
def load_annotation_rules(path):
    df = pd.read_csv(path, sep="	").fillna("")
    compiled = []
    for _, r in df.iterrows():
        policy = str(r.get("match_policy", "product_or_gene")).strip() or "product_or_gene"
        compiled.append((r, re.compile(str(r["regex"]), re.I), policy))
    return compiled

def screen_annotated(cds, rules):
    """
    Product-aware annotation screening.

    For rules marked `product_required`, a gene symbol alone is not sufficient.
    This prevents ambiguous symbols such as xylE from overriding an explicit,
    contradictory product annotation (e.g. D-xylose-proton symporter).
    """
    hits = []
    for r in cds:
        product = r["product"] or ""
        if re.search(r"hypothetical protein|uncharacterized protein|unknown protein", product, re.I):
            continue

        gene_text = " | ".join([r["locus_tag"], r["gene"]])
        product_text = " | ".join([product, r["EC_number"], r["db_xref"]])

        for rule, rx, policy in rules:
            gene_match = bool(rx.search(gene_text))
            product_match = bool(rx.search(product_text))

            if policy == "product_required":
                if not product_match:
                    continue
            else:
                if not (gene_match or product_match):
                    continue

            hits.append({
                "locus_tag": r["locus_tag"],
                "gene": r["gene"],
                "product": product,
                "contig": r["contig"],
                "protein_length_aa": len(r["translation"]) if r["translation"] else "",
                "family": rule["family"],
                "major_category": rule["major_category"],
                "evidence_class": rule["evidence_class"],
                "confidence": rule["default_confidence"],
                "evidence_source": "Existing annotation",
                "interpretation": (
                    "Annotation-based candidate; product-aware matching applied. "
                    "Direct biological activity still requires experimental validation."
                ),
            })
    return pd.DataFrame(hits)

def extract_hypotheticals(cds, out_faa):
    rows = []
    with open(out_faa, "w", encoding="utf-8") as f:
        for r in cds:
            product = r["product"] or ""
            if re.search(r"hypothetical protein|uncharacterized protein|unknown protein", product, re.I):
                seq = r["translation"]
                if not seq:
                    continue
                rows.append({
                    "locus_tag": r["locus_tag"],
                    "contig": r["contig"],
                    "product": product,
                    "protein_length_aa": len(seq),
                })
                f.write(f">{r['locus_tag']} | original_product={product} | contig={r['contig']}\n")
                for i in range(0, len(seq), 80):
                    f.write(seq[i:i+80] + "\n")
    return pd.DataFrame(rows)

# ----------------------------
# DIAMOND
# ----------------------------
DIAMOND_COLS = [
    "query","subject","identity_pct","alignment_len","query_len","subject_len",
    "evalue","bitscore","query_coverage_pct"
]

def run_diamond(diamond, query_faa, db, out_tsv, threads=None):
    cmd = [
        diamond, "blastp",
        "--query", str(query_faa),
        "--db", str(db),
        "--out", str(out_tsv),
        "--outfmt", "6", "qseqid", "sseqid", "pident", "length", "qlen", "slen", "evalue", "bitscore", "qcovhsp",
        "--max-target-seqs", "10",
        "--evalue", "1e-5",
        "--sensitive",
    ]
    if threads:
        cmd += ["--threads", str(threads)]
    print("[run]", " ".join(cmd))
    subprocess.run(cmd, check=True)

def load_diamond(path):
    if not Path(path).exists() or Path(path).stat().st_size == 0:
        return pd.DataFrame(columns=DIAMOND_COLS)
    df = pd.read_csv(path, sep="\t", header=None, names=DIAMOND_COLS)
    for c in ["identity_pct","alignment_len","query_len","subject_len","evalue","bitscore","query_coverage_pct"]:
        df[c] = pd.to_numeric(df[c], errors="coerce")
    df["reference_coverage_pct"] = (df["alignment_len"] / df["subject_len"] * 100).round(1)
    df["accession"] = df["subject"].astype(str).str.split("|").str[0]
    return df

# ----------------------------
# Reference metadata + scoring
# ----------------------------
def load_reference_metadata(path):
    return pd.read_csv(path, sep="\t", dtype=str).fillna("")

def base_confidence(row):
    pid = float(row.get("identity_pct", 0) or 0)
    qcov = float(row.get("query_coverage_pct", 0) or 0)
    rcov = float(row.get("reference_coverage_pct", 0) or 0)
    e = float(row.get("evalue", 1) or 1)
    qlen = float(row.get("query_len", 0) or 0)
    slen = float(row.get("subject_len", 0) or 0)

    partial = (rcov < 50) or (qlen < 0.5 * slen)
    if e <= 1e-30 and pid >= 30 and qcov >= 65 and rcov >= 65 and not partial:
        return "High"
    if e <= 1e-10 and pid >= 25 and qcov >= 50 and rcov >= 40:
        return "Moderate"
    return "Weak"

def classify_hypothetical_hits(df, meta):
    if df.empty:
        return pd.DataFrame()
    merged = df.merge(meta, left_on="accession", right_on="Accession", how="left")
    merged["confidence"] = merged.apply(base_confidence, axis=1)
    merged["fragment_or_partial"] = (
        (merged["reference_coverage_pct"] < 50) |
        (merged["query_len"] < 0.5 * merged["subject_len"])
    )

    # Keep best hit per query by bitscore, then e-value.
    best = (
        merged.sort_values(["query","bitscore","evalue"], ascending=[True,False,True])
              .groupby("query", as_index=False)
              .first()
    )

    # Category safeguards from v0.2.1 metadata.
    if "Category_flag" not in best.columns:
        best["Category_flag"] = "REVIEW"
    if "Headline_category" not in best.columns:
        best["Headline_category"] = "Unresolved"

    # Strong rejection rules.
    def decision(r):
        if bool(r["fragment_or_partial"]) and float(r["reference_coverage_pct"]) < 20:
            return "Reject"
        if str(r.get("Category_flag","")).upper() == "REVIEW":
            return "Review"
        if r["confidence"] == "Weak":
            return "Review"
        return "Candidate"

    best["decision"] = best.apply(decision, axis=1)
    best["evidence_source"] = "DIAMOND homology"
    best["interpretation"] = best.apply(
        lambda r: (
            "Short/partial similarity; do not assign full-length function."
            if r["decision"] == "Reject"
            else "Sequence-supported candidate; conserved-domain validation is recommended before a strong functional call."
        ),
        axis=1
    )
    return best

# ----------------------------
# InterPro optional integration
# ----------------------------
def load_interpro(path):
    if not path:
        return pd.DataFrame()
    cols = [
        "protein","md5","length","analysis","signature_accession","signature_description",
        "start","stop","score","status","date","interpro_accession","interpro_description",
        "go_terms","pathways"
    ]
    rows = []
    with open(path, encoding="utf-8") as f:
        for vals in csv.reader(f, delimiter="\t"):
            if len(vals) < 13:
                continue
            vals = vals + [""] * (15-len(vals))
            rows.append(vals[:15])
    return pd.DataFrame(rows, columns=cols)

def add_interpro_support(best, ip):
    if best.empty or ip.empty:
        return best
    support = defaultdict(list)
    for _, r in ip.iterrows():
        for x in [r["signature_description"], r["interpro_description"]]:
            x = str(x).strip()
            if x and x != "-":
                support[str(r["protein"])].append(x)
    def uniq_join(xs):
        out=[]; seen=set()
        for x in xs:
            if x not in seen:
                seen.add(x); out.append(x)
        return "; ".join(out)
    best["interpro_support"] = best["query"].map(lambda q: uniq_join(support.get(str(q), [])))
    best["interpro_detected"] = best["interpro_support"].astype(str).str.len() > 0
    best.loc[best["interpro_detected"], "evidence_source"] = "DIAMOND + InterPro support"
    return best

# ----------------------------
# Report creation
# ----------------------------
CONF_ORDER = {"High":0, "Moderate":1, "Weak":2}

def make_report(out_xlsx, annotated, hyp_best, hyp_meta, diamond_all, interpro, total_cds):

    final_rows = []

    if not annotated.empty:
        for _, r in annotated.iterrows():
            final_rows.append({
                "Locus_tag": r["locus_tag"],
                "Gene": r["gene"],
                "Product_or_prediction": r["product"],
                "Origin": "Annotated",
                "Category": r["major_category"],
                "Evidence_class": r["evidence_class"],
                "Confidence": r["confidence"],
                "Decision": "Candidate",
                "Evidence_source": r["evidence_source"],
                "Identity_pct": "",
                "Query_coverage_pct": "",
                "Reference_coverage_pct": "",
                "Evalue": "",
                "Best_reference": "",
                "Interpretation": r["interpretation"],
            })

    if hyp_best is not None and not hyp_best.empty:
        for _, r in hyp_best.iterrows():
            final_rows.append({
                "Locus_tag": r["query"],
                "Gene": "",
                "Product_or_prediction": r.get("Protein","") or r.get("All_family_targets","") or r.get("subject",""),
                "Origin": "Hypothetical → predicted",
                "Category": r.get("Headline_category","Unresolved"),
                "Evidence_class": r.get("Evidence_classes",""),
                "Confidence": r["confidence"],
                "Decision": r["decision"],
                "Evidence_source": r["evidence_source"],
                "Identity_pct": round(float(r["identity_pct"]),1),
                "Query_coverage_pct": round(float(r["query_coverage_pct"]),1),
                "Reference_coverage_pct": round(float(r["reference_coverage_pct"]),1),
                "Evalue": r["evalue"],
                "Best_reference": r["subject"],
                "Interpretation": r["interpretation"] + (
                    (" InterPro support: " + str(r.get("interpro_support",""))[:500])
                    if str(r.get("interpro_support","")).strip() else ""
                ),
            })

    final = pd.DataFrame(final_rows)
    if not final.empty:
        final["_rank"] = final["Confidence"].map(CONF_ORDER).fillna(9)
        final = final.sort_values(["_rank","Category","Locus_tag"]).drop(columns="_rank")

    high = final[(final["Confidence"]=="High") & (final["Decision"]=="Candidate")] if not final.empty else final
    micro = final[final["Category"].astype(str).str.contains("Microplastics", case=False, na=False)] if not final.empty else final
    supporting = final[final["Evidence_class"].astype(str).str.contains("Supporting", case=False, na=False)] if not final.empty else final
    review = final[final["Decision"].isin(["Review","Reject"])] if not final.empty else final

    summary = pd.DataFrame([
        ["Total CDS", total_cds],
        ["Annotated candidate rows", len(annotated)],
        ["Hypothetical proteins screened", len(hyp_meta)],
        ["Hypothetical queries with DIAMOND hits", 0 if hyp_best is None or hyp_best.empty else hyp_best["query"].nunique()],
        ["High-confidence final candidates", len(high)],
        ["Microplastic/polymer candidates", len(micro)],
        ["Review/rejected hypothetical hits", len(review)],
    ], columns=["Metric","Value"])

    with pd.ExcelWriter(out_xlsx, engine="openpyxl") as writer:
        summary.to_excel(writer, sheet_name="Summary", index=False)
        final.to_excel(writer, sheet_name="All_candidates", index=False)
        high.to_excel(writer, sheet_name="High_confidence", index=False)
        annotated.to_excel(writer, sheet_name="Annotated_candidates", index=False)
        if hyp_best is not None:
            hyp_best.to_excel(writer, sheet_name="Predicted_hypothetical", index=False)
        micro.to_excel(writer, sheet_name="Microplastic_candidates", index=False)
        supporting.to_excel(writer, sheet_name="Supporting_functions", index=False)
        review.to_excel(writer, sheet_name="Review_required", index=False)
        diamond_all.to_excel(writer, sheet_name="All_DIAMOND_hits", index=False)
        if not interpro.empty:
            interpro.to_excel(writer, sheet_name="InterPro_evidence", index=False)

        # Basic formatting
        for ws in writer.book.worksheets:
            ws.freeze_panes = "A2"
            ws.auto_filter.ref = ws.dimensions
            for cell in ws[1]:
                new_font = copy(cell.font)
                new_font.bold = True
                cell.font = new_font
            for col_cells in ws.columns:
                max_len = min(max(len(str(c.value or "")) for c in col_cells[:200]) + 2, 45)
                ws.column_dimensions[col_cells[0].column_letter].width = max(10, max_len)

    # TSV companions
    stem = Path(out_xlsx).with_suffix("")
    final.to_csv(str(stem) + "_all_candidates.tsv", sep="\t", index=False)
    review.to_csv(str(stem) + "_review_required.tsv", sep="\t", index=False)
    micro.to_csv(str(stem) + "_microplastic_candidates.tsv", sep="\t", index=False)

    return final

def export_interpro_candidates(cds, hyp_best, out_faa):
    if hyp_best is None or hyp_best.empty:
        Path(out_faa).write_text("", encoding="utf-8")
        return
    wanted = set(
        hyp_best.loc[
            (hyp_best["decision"].isin(["Candidate","Review"])) &
            (hyp_best["confidence"].isin(["High","Moderate"])),
            "query"
        ].astype(str)
    )
    seqmap = {r["locus_tag"]: r["translation"] for r in cds if r["translation"]}
    with open(out_faa, "w", encoding="utf-8") as f:
        for locus in sorted(wanted):
            seq = seqmap.get(locus, "")
            if not seq:
                continue
            f.write(f">{locus}\n")
            for i in range(0, len(seq), 80):
                f.write(seq[i:i+80] + "\n")

# ----------------------------
# CLI
# ----------------------------
def main():
    ap = argparse.ArgumentParser(description="Mine bacterial WGS annotations for evidence-ranked bioremediation candidates.")
    ap.add_argument("genbank", help="Input GenBank annotation (.gbk/.gbff)")
    ap.add_argument("--reference-db", required=True, help="DIAMOND database basename/path (.dmnd basename)")
    ap.add_argument("--reference-metadata", required=True, help="v0.2.1 reference metadata TSV")
    ap.add_argument("--annotation-rules", default="annotation_rules.tsv", help="Annotated-CDS regex rules TSV")
    ap.add_argument("--diamond", default="diamond", help="DIAMOND executable/path")
    ap.add_argument("--interpro", help="Optional InterProScan TSV for selected hypothetical candidates")
    ap.add_argument("--outdir", default="bioremediation_gene_miner_results")
    ap.add_argument("--threads", type=int)
    ap.add_argument("--skip-diamond", action="store_true", help="Use an existing DIAMOND TSV in outdir")
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    print(f"Bioremediation Gene Miner v{VERSION}")
    cds = parse_genbank(args.genbank)
    print("[info] CDS parsed:", len(cds))

    rules = load_annotation_rules(args.annotation_rules)
    annotated = screen_annotated(cds, rules)
    annotated.to_csv(outdir/"annotated_candidates.tsv", sep="\t", index=False)
    print("[info] annotated candidate rows:", len(annotated))

    hyp_faa = outdir/"hypothetical_proteins.faa"
    hyp_meta = extract_hypotheticals(cds, hyp_faa)
    hyp_meta.to_csv(outdir/"hypothetical_proteins.tsv", sep="\t", index=False)
    print("[info] hypothetical proteins:", len(hyp_meta))

    diamond_tsv = outdir/"hypothetical_vs_reference.tsv"
    if not args.skip_diamond:
        run_diamond(args.diamond, hyp_faa, args.reference_db, diamond_tsv, args.threads)
    diamond_all = load_diamond(diamond_tsv)
    print("[info] DIAMOND alignments:", len(diamond_all))

    refmeta = load_reference_metadata(args.reference_metadata)
    hyp_best = classify_hypothetical_hits(diamond_all, refmeta)

    ip = load_interpro(args.interpro)
    hyp_best = add_interpro_support(hyp_best, ip)

    interpro_faa = outdir/"interpro_candidates.faa"
    export_interpro_candidates(cds, hyp_best, interpro_faa)

    report = outdir/"Bioremediation_Gene_Miner_Report.xlsx"
    final = make_report(report, annotated, hyp_best, hyp_meta, diamond_all, ip, len(cds))

    print("\nDONE")
    print("Report:", report.resolve())
    print("InterPro candidate FASTA:", interpro_faa.resolve())
    if not final.empty:
        print("Final candidate rows:", len(final))
        print("High-confidence:", int(((final["Confidence"]=="High") & (final["Decision"]=="Candidate")).sum()))
        print("Microplastic/polymer:", int(final["Category"].astype(str).str.contains("Microplastics", case=False, na=False).sum()))
    print("\nReminder: predictions are candidates, not experimental confirmation.")

if __name__ == "__main__":
    main()
