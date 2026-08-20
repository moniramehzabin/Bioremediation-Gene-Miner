#!/usr/bin/env python3
"""
Bioremediation Gene Miner v0.3.4

Evidence-ranked bacterial WGS screening:
1) screens already-annotated CDSs using annotation_rules.tsv
2) extracts hypothetical/uncharacterized proteins
3) runs DIAMOND against the curated reference database
4) preserves the best hit PER QUERY AND PER CURATED BIOLOGICAL FAMILY
5) applies identity/query-coverage/reference-coverage safeguards
6) optionally integrates InterProScan TSV evidence
7) writes Excel and TSV reports

Important:
This software predicts candidates. It does not experimentally prove
bioremediation activity.
"""
from __future__ import annotations

import argparse
from copy import copy
import csv
import re
import subprocess
from collections import defaultdict
from pathlib import Path

import pandas as pd

try:
    from Bio import SeqIO
except Exception:
    SeqIO = None

VERSION = "0.3.5-test"


# ============================================================
# GenBank parsing
# ============================================================

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
    """
    Tolerant parser for Prokka/GenBank files that Biopython may reject because
    of nonstandard LOCUS formatting.
    """
    text = Path(path).read_text(encoding="utf-8", errors="replace")

    contig_starts = list(re.finditer(r"^LOCUS\s+(\S+)", text, flags=re.M))
    contig_positions = [(m.start(), m.group(1)) for m in contig_starts]
    cds_starts = list(re.finditer(r"^     CDS\s+(.+)$", text, flags=re.M))

    cds = []
    for i, m in enumerate(cds_starts):
        start = m.start()
        end = cds_starts[i + 1].start() if i + 1 < len(cds_starts) else len(text)
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
    except Exception as exc:
        print(
            f"[info] Biopython GenBank parser failed ({exc}); "
            "using tolerant fallback parser."
        )
        return parse_genbank_fallback(path)


# ============================================================
# Existing-annotation screen
# ============================================================

def load_annotation_rules(path):
    df = pd.read_csv(path, sep="\t").fillna("")
    compiled = []

    for _, r in df.iterrows():
        policy = str(r.get("match_policy", "product_or_gene")).strip()
        policy = policy or "product_or_gene"
        compiled.append((r, re.compile(str(r["regex"]), re.I), policy))

    return compiled


def screen_annotated(cds, rules):
    """
    Product-aware annotation screening.

    Rules marked product_required cannot be triggered by an ambiguous gene
    symbol alone. This retains the v0.3.3 safeguard against misleading calls.
    """
    hits = []

    for r in cds:
        product = r["product"] or ""

        if re.search(
            r"hypothetical protein|uncharacterized protein|unknown protein",
            product,
            re.I,
        ):
            continue

        gene_text = " | ".join([r["locus_tag"], r["gene"]])
        product_text = " | ".join(
            [product, r["EC_number"], r["db_xref"]]
        )

        for rule, rx, policy in rules:
            gene_match = bool(rx.search(gene_text))
            product_match = bool(rx.search(product_text))

            if policy == "product_required":
                if not product_match:
                    continue
            elif not (gene_match or product_match):
                continue

            hits.append({
                "locus_tag": r["locus_tag"],
                "gene": r["gene"],
                "product": product,
                "contig": r["contig"],
                "protein_length_aa": (
                    len(r["translation"]) if r["translation"] else ""
                ),
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

            if not re.search(
                r"hypothetical protein|uncharacterized protein|unknown protein",
                product,
                re.I,
            ):
                continue

            seq = r["translation"]
            if not seq:
                continue

            rows.append({
                "locus_tag": r["locus_tag"],
                "contig": r["contig"],
                "product": product,
                "protein_length_aa": len(seq),
            })

            f.write(
                f">{r['locus_tag']} | original_product={product} "
                f"| contig={r['contig']}\n"
            )
            for i in range(0, len(seq), 80):
                f.write(seq[i:i + 80] + "\n")

    return pd.DataFrame(rows)


# ============================================================
# DIAMOND
# ============================================================

DIAMOND_COLS = [
    "query",
    "subject",
    "identity_pct",
    "alignment_len",
    "query_len",
    "subject_len",
    "evalue",
    "bitscore",
    "query_coverage_pct",
]


def run_diamond(diamond, query_faa, db, out_tsv, threads=None):
    """
    v0.3.4 deliberately requests more reference hits than v0.3.3.

    The curated reference DB contains related reductase/resistance families.
    Keeping only the first 10 DIAMOND targets can prevent a biologically
    relevant family from ever reaching the family-aware classifier.
    """
    cmd = [
        diamond,
        "blastp",
        "--query", str(query_faa),
        "--db", str(db),
        "--out", str(out_tsv),
        "--outfmt", "6",
        "qseqid", "sseqid", "pident", "length", "qlen", "slen",
        "evalue", "bitscore", "qcovhsp",
        "--max-target-seqs", "50",
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

    df = pd.read_csv(
        path,
        sep="\t",
        header=None,
        names=DIAMOND_COLS,
    )

    numeric_cols = [
        "identity_pct",
        "alignment_len",
        "query_len",
        "subject_len",
        "evalue",
        "bitscore",
        "query_coverage_pct",
    ]
    for c in numeric_cols:
        df[c] = pd.to_numeric(df[c], errors="coerce")

    df["reference_coverage_pct"] = (
        df["alignment_len"] / df["subject_len"] * 100
    ).round(1)

    df["accession"] = (
        df["subject"].astype(str).str.split("|", regex=False).str[0]
    )

    return df


# ============================================================
# Reference metadata + evidence scoring
# ============================================================

def load_reference_metadata(path):
    meta = pd.read_csv(path, sep="\t", dtype=str).fillna("")

    if "Accession" not in meta.columns:
        raise ValueError(
            "Reference metadata must contain an 'Accession' column."
        )

    return meta


def base_confidence(row):
    pid = float(row.get("identity_pct", 0) or 0)
    qcov = float(row.get("query_coverage_pct", 0) or 0)
    rcov = float(row.get("reference_coverage_pct", 0) or 0)
    evalue = float(row.get("evalue", 1) or 1)
    qlen = float(row.get("query_len", 0) or 0)
    slen = float(row.get("subject_len", 0) or 0)

    partial = (rcov < 50) or (qlen < 0.5 * slen)

    if (
        evalue <= 1e-30
        and pid >= 30
        and qcov >= 65
        and rcov >= 65
        and not partial
    ):
        return "High"

    if (
        evalue <= 1e-10
        and pid >= 25
        and qcov >= 50
        and rcov >= 40
    ):
        return "Moderate"

    return "Weak"


def _split_families(value):
    """
    Metadata can associate one UniProt accession with several curated targets.
    Exploding these targets prevents one combined metadata string from hiding
    an individual family such as chrR, chrA, azoR, or azoreductase.
    """
    text = str(value or "").strip()

    if not text:
        return ["Unresolved"]

    parts = [
        p.strip()
        for p in re.split(r"\s*;\s*", text)
        if p.strip()
    ]
    return parts or ["Unresolved"]


def classify_hypothetical_hits(df, meta):
    """
    Core v0.3.4 fix.

    v0.3.3 kept ONE global best hit per hypothetical query. v0.3.4 instead:
      1. joins every DIAMOND hit to curated metadata;
      2. expands multi-family metadata into individual family rows;
      3. keeps the best reference hit for EACH query + biological family.

    Therefore a valid ChrR/AzoR/ChrA/etc. family hit is not silently removed
    merely because the same query has a higher-scoring hit to another family.
    """
    if df.empty:
        return pd.DataFrame()

    merged = df.merge(
        meta,
        left_on="accession",
        right_on="Accession",
        how="left",
    )

    # Do not silently turn references missing from metadata into candidates.
    merged["metadata_matched"] = (
        merged["Accession"].astype(str).str.strip() != ""
    )

    if "All_family_targets" not in merged.columns:
        merged["All_family_targets"] = "Unresolved"

    merged["family_target"] = merged["All_family_targets"].apply(_split_families)
    merged = merged.explode("family_target", ignore_index=True)
    merged["family_target"] = (
        merged["family_target"].astype(str).str.strip().replace("", "Unresolved")
    )

    merged["confidence"] = merged.apply(base_confidence, axis=1)

    merged["fragment_or_partial"] = (
        (merged["reference_coverage_pct"] < 50)
        | (merged["query_len"] < 0.5 * merged["subject_len"])
    )

    # Best reference within EACH biological family for EACH query.
    best = (
        merged.sort_values(
            ["query", "family_target", "bitscore", "evalue"],
            ascending=[True, True, False, True],
        )
        .groupby(
            ["query", "family_target"],
            as_index=False,
            dropna=False,
        )
        .first()
    )

    if "Category_flag" not in best.columns:
        best["Category_flag"] = "REVIEW"

    if "Headline_category" not in best.columns:
        best["Headline_category"] = "Unresolved"

    if "Evidence_classes" not in best.columns:
        best["Evidence_classes"] = ""

    def decision(r):
        if not bool(r.get("metadata_matched", False)):
            return "Reject"

        rcov = float(r.get("reference_coverage_pct", 0) or 0)

        if bool(r["fragment_or_partial"]) and rcov < 20:
            return "Reject"

        if str(r.get("Category_flag", "")).upper() == "REVIEW":
            return "Review"

        if r["confidence"] == "Weak":
            return "Review"

        return "Candidate"

    best["decision"] = best.apply(decision, axis=1)
    # Keep weak evidence visible; this label is descriptive, not a filter.
    best["match_strength"] = best["confidence"].astype(str) + " match"
    best["evidence_source"] = "DIAMOND homology"

    def interpretation(r):
        if not bool(r.get("metadata_matched", False)):
            return (
                "DIAMOND reference was not found in the supplied metadata; "
                "functional assignment rejected."
            )

        if r["decision"] == "Reject":
            return (
                "Short/partial or unsupported similarity; "
                "do not assign full-length function."
            )

        return (
            "Sequence-supported family-level candidate; conserved-domain "
            "validation is recommended before a strong functional call."
        )

    best["interpretation"] = best.apply(interpretation, axis=1)

    # Stable, useful ordering.
    conf_rank = {"High": 0, "Moderate": 1, "Weak": 2}
    decision_rank = {"Candidate": 0, "Review": 1, "Reject": 2}
    best["_conf_rank"] = best["confidence"].map(conf_rank).fillna(9)
    best["_decision_rank"] = best["decision"].map(decision_rank).fillna(9)

    best = (
        best.sort_values(
            [
                "query",
                "_decision_rank",
                "_conf_rank",
                "family_target",
                "bitscore",
            ],
            ascending=[True, True, True, True, False],
        )
        .drop(columns=["_conf_rank", "_decision_rank"])
        .reset_index(drop=True)
    )

    return best


# ============================================================
# Optional InterPro integration
# ============================================================

def load_interpro(path):
    if not path:
        return pd.DataFrame()

    cols = [
        "protein",
        "md5",
        "length",
        "analysis",
        "signature_accession",
        "signature_description",
        "start",
        "stop",
        "score",
        "status",
        "date",
        "interpro_accession",
        "interpro_description",
        "go_terms",
        "pathways",
    ]

    rows = []
    with open(path, encoding="utf-8") as f:
        for vals in csv.reader(f, delimiter="\t"):
            if len(vals) < 13:
                continue

            vals = vals + [""] * (15 - len(vals))
            rows.append(vals[:15])

    return pd.DataFrame(rows, columns=cols)


def add_interpro_support(best, ip):
    """
    Add InterPro evidence WITHOUT changing, filtering, hiding, or upgrading
    the original DIAMOND confidence/decision.

    A Weak DIAMOND match remains Weak/Review even when InterPro supplies
    strong domain-level support. InterPro is reported as an independent,
    additive evidence layer.
    """
    if best is None or best.empty:
        return best

    best = best.copy()

    # Always create the columns so report structure is stable even when no
    # InterPro file was supplied.
    best["interpro_support"] = ""
    best["interpro_accessions"] = ""
    best["interpro_analyses"] = ""
    best["interpro_detected"] = False

    if ip is None or ip.empty:
        return best

    descriptions = defaultdict(list)
    accessions = defaultdict(list)
    analyses = defaultdict(list)

    def add_unique(mapping, protein, value):
        value = str(value).strip()
        if value and value != "-" and value not in mapping[protein]:
            mapping[protein].append(value)

    for _, r in ip.iterrows():
        protein = str(r["protein"]).strip()
        if not protein:
            continue

        add_unique(descriptions, protein, r.get("signature_description", ""))
        add_unique(descriptions, protein, r.get("interpro_description", ""))
        add_unique(accessions, protein, r.get("signature_accession", ""))
        add_unique(accessions, protein, r.get("interpro_accession", ""))
        add_unique(analyses, protein, r.get("analysis", ""))

    best["interpro_support"] = best["query"].map(
        lambda q: "; ".join(descriptions.get(str(q), []))
    )
    best["interpro_accessions"] = best["query"].map(
        lambda q: "; ".join(accessions.get(str(q), []))
    )
    best["interpro_analyses"] = best["query"].map(
        lambda q: "; ".join(analyses.get(str(q), []))
    )
    best["interpro_detected"] = (
        best["interpro_support"].astype(str).str.len() > 0
    )

    # Evidence source is additive only. DO NOT alter confidence or decision.
    best.loc[
        best["interpro_detected"],
        "evidence_source",
    ] = "DIAMOND + InterPro (independent domain evidence)"

    return best


# ============================================================
# Report creation
# ============================================================

CONF_ORDER = {
    "High": 0,
    "Moderate": 1,
    "Weak": 2,
}


def make_report(
    out_xlsx,
    annotated,
    hyp_best,
    hyp_meta,
    diamond_all,
    interpro,
    total_cds,
):
    final_rows = []

    if annotated is not None and not annotated.empty:
        for _, r in annotated.iterrows():
            final_rows.append({
                "Locus_tag": r["locus_tag"],
                "Gene": r["gene"],
                "Family_target": r.get("family", ""),
                "Product_or_prediction": r["product"],
                "Origin": "Annotated",
                "Category": r["major_category"],
                "Evidence_class": r["evidence_class"],
                "Confidence": r["confidence"],
                "Match_strength": f"{r['confidence']} match",
                "Decision": "Candidate",
                "Evidence_source": r["evidence_source"],
                "Identity_pct": "",
                "Query_coverage_pct": "",
                "Reference_coverage_pct": "",
                "Evalue": "",
                "Bitscore": "",
                "Best_reference": "",
                "InterPro_detected": "",
                "InterPro_analyses": "",
                "InterPro_accessions": "",
                "InterPro_domain_support": "",
                "Interpretation": r["interpretation"],
            })

    if hyp_best is not None and not hyp_best.empty:
        for _, r in hyp_best.iterrows():
            prediction = (
                r.get("Protein", "")
                or r.get("family_target", "")
                or r.get("subject", "")
            )

            final_rows.append({
                "Locus_tag": r["query"],
                "Gene": "",
                "Family_target": r.get("family_target", ""),
                "Product_or_prediction": prediction,
                "Origin": "Hypothetical → predicted",
                "Category": r.get("Headline_category", "Unresolved"),
                "Evidence_class": r.get("Evidence_classes", ""),
                "Confidence": r["confidence"],
                "Match_strength": f"{r['confidence']} match",
                "Decision": r["decision"],
                "Evidence_source": r["evidence_source"],
                "Identity_pct": round(float(r["identity_pct"]), 1),
                "Query_coverage_pct": round(
                    float(r["query_coverage_pct"]), 1
                ),
                "Reference_coverage_pct": round(
                    float(r["reference_coverage_pct"]), 1
                ),
                "Evalue": r["evalue"],
                "Bitscore": round(float(r["bitscore"]), 1),
                "Best_reference": r["subject"],
                "InterPro_detected": bool(r.get("interpro_detected", False)),
                "InterPro_analyses": r.get("interpro_analyses", ""),
                "InterPro_accessions": r.get("interpro_accessions", ""),
                "InterPro_domain_support": r.get("interpro_support", ""),
                "Interpretation": (
                    r["interpretation"]
                    + (
                        " Independent InterPro domain evidence is present; "
                        "the DIAMOND confidence and decision above are intentionally "
                        "preserved and must be interpreted separately."
                        if bool(r.get("interpro_detected", False))
                        else ""
                    )
                ),
            })

    final = pd.DataFrame(final_rows)

    if not final.empty:
        final["_rank"] = (
            final["Confidence"].map(CONF_ORDER).fillna(9)
        )
        final = (
            final.sort_values(
                [
                    "_rank",
                    "Category",
                    "Family_target",
                    "Locus_tag",
                ]
            )
            .drop(columns="_rank")
            .reset_index(drop=True)
        )

    if not final.empty:
        high = final[
            (final["Confidence"] == "High")
            & (final["Decision"] == "Candidate")
        ].copy()

        micro = final[
            final["Category"].astype(str).str.contains(
                "Microplastics",
                case=False,
                na=False,
            )
        ].copy()

        supporting = final[
            final["Evidence_class"].astype(str).str.contains(
                "Supporting",
                case=False,
                na=False,
            )
        ].copy()

        review = final[
            final["Decision"].isin(["Review", "Reject"])
        ].copy()

        family_candidates = final[
            final["Origin"].eq("Hypothetical → predicted")
            & final["Decision"].isin(["Candidate", "Review"])
        ].copy()

        interpro_supported = final[
            final["InterPro_detected"].eq(True)
        ].copy()
    else:
        high = final.copy()
        micro = final.copy()
        supporting = final.copy()
        review = final.copy()
        family_candidates = final.copy()
        interpro_supported = final.copy()

    hyp_queries = (
        0
        if hyp_best is None or hyp_best.empty
        else hyp_best["query"].nunique()
    )

    family_hits = (
        0
        if hyp_best is None or hyp_best.empty
        else len(hyp_best)
    )

    summary = pd.DataFrame(
        [
            ["Total CDS", total_cds],
            ["Annotated candidate rows", len(annotated)],
            ["Hypothetical proteins screened", len(hyp_meta)],
            ["Hypothetical queries with DIAMOND hits", hyp_queries],
            ["Family-level hypothetical hit rows", family_hits],
            ["InterPro-supported final rows", len(interpro_supported)],
            ["High-confidence final candidates", len(high)],
            ["Microplastic/polymer candidates", len(micro)],
            ["Review/rejected final rows", len(review)],
        ],
        columns=["Metric", "Value"],
    )

    with pd.ExcelWriter(out_xlsx, engine="openpyxl") as writer:
        summary.to_excel(
            writer,
            sheet_name="Summary",
            index=False,
        )

        final.to_excel(
            writer,
            sheet_name="All_candidates",
            index=False,
        )

        high.to_excel(
            writer,
            sheet_name="High_confidence",
            index=False,
        )

        annotated.to_excel(
            writer,
            sheet_name="Annotated_candidates",
            index=False,
        )

        if hyp_best is not None:
            hyp_best.to_excel(
                writer,
                sheet_name="Predicted_hypothetical",
                index=False,
            )

        family_candidates.to_excel(
            writer,
            sheet_name="Family_level_candidates",
            index=False,
        )

        interpro_supported.to_excel(
            writer,
            sheet_name="InterPro_supported",
            index=False,
        )

        micro.to_excel(
            writer,
            sheet_name="Microplastic_candidates",
            index=False,
        )

        supporting.to_excel(
            writer,
            sheet_name="Supporting_functions",
            index=False,
        )

        review.to_excel(
            writer,
            sheet_name="Review_required",
            index=False,
        )

        diamond_all.to_excel(
            writer,
            sheet_name="All_DIAMOND_hits",
            index=False,
        )

        if interpro is not None and not interpro.empty:
            interpro.to_excel(
                writer,
                sheet_name="InterPro_evidence",
                index=False,
            )

        # Basic formatting only; no biological information is altered here.
        for ws in writer.book.worksheets:
            ws.freeze_panes = "A2"
            ws.auto_filter.ref = ws.dimensions

            for cell in ws[1]:
                new_font = copy(cell.font)
                new_font.bold = True
                cell.font = new_font

            for col_cells in ws.columns:
                sampled = col_cells[:200]
                max_len = min(
                    max(len(str(c.value or "")) for c in sampled) + 2,
                    45,
                )
                ws.column_dimensions[
                    col_cells[0].column_letter
                ].width = max(10, max_len)

    stem = Path(out_xlsx).with_suffix("")

    final.to_csv(
        str(stem) + "_all_candidates.tsv",
        sep="\t",
        index=False,
    )

    review.to_csv(
        str(stem) + "_review_required.tsv",
        sep="\t",
        index=False,
    )

    micro.to_csv(
        str(stem) + "_microplastic_candidates.tsv",
        sep="\t",
        index=False,
    )

    family_candidates.to_csv(
        str(stem) + "_family_level_candidates.tsv",
        sep="\t",
        index=False,
    )

    return final


def export_interpro_candidates(cds, hyp_best, out_faa):
    if hyp_best is None or hyp_best.empty:
        Path(out_faa).write_text("", encoding="utf-8")
        return

    wanted = set(
        hyp_best.loc[
            hyp_best["decision"].isin(["Candidate", "Review"])
            & hyp_best["confidence"].isin(["High", "Moderate"]),
            "query",
        ].astype(str)
    )

    seqmap = {
        r["locus_tag"]: r["translation"]
        for r in cds
        if r["translation"]
    }

    with open(out_faa, "w", encoding="utf-8") as f:
        for locus in sorted(wanted):
            seq = seqmap.get(locus, "")
            if not seq:
                continue

            f.write(f">{locus}\n")
            for i in range(0, len(seq), 80):
                f.write(seq[i:i + 80] + "\n")


# ============================================================
# CLI
# ============================================================

def main():
    ap = argparse.ArgumentParser(
        description=(
            "Mine bacterial WGS annotations for evidence-ranked "
            "bioremediation candidates."
        )
    )

    ap.add_argument(
        "genbank",
        help="Input GenBank annotation (.gbk/.gbff)",
    )

    ap.add_argument(
        "--reference-db",
        required=True,
        help="DIAMOND database basename/path (.dmnd basename)",
    )

    ap.add_argument(
        "--reference-metadata",
        required=True,
        help="Curated reference metadata TSV",
    )

    ap.add_argument(
        "--annotation-rules",
        default="annotation_rules.tsv",
        help="Annotated-CDS regex rules TSV",
    )

    ap.add_argument(
        "--diamond",
        default="diamond",
        help="DIAMOND executable/path",
    )

    ap.add_argument(
        "--interpro",
        help="Optional InterProScan TSV for selected hypothetical candidates",
    )

    ap.add_argument(
        "--outdir",
        default="bioremediation_gene_miner_results",
    )

    ap.add_argument(
        "--threads",
        type=int,
    )

    ap.add_argument(
        "--skip-diamond",
        action="store_true",
        help="Use an existing DIAMOND TSV in outdir",
    )

    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    print(f"Bioremediation Gene Miner v{VERSION}")

    cds = parse_genbank(args.genbank)
    print("[info] CDS parsed:", len(cds))

    rules = load_annotation_rules(args.annotation_rules)
    annotated = screen_annotated(cds, rules)

    annotated.to_csv(
        outdir / "annotated_candidates.tsv",
        sep="\t",
        index=False,
    )
    print("[info] annotated candidate rows:", len(annotated))

    hyp_faa = outdir / "hypothetical_proteins.faa"
    hyp_meta = extract_hypotheticals(cds, hyp_faa)

    hyp_meta.to_csv(
        outdir / "hypothetical_proteins.tsv",
        sep="\t",
        index=False,
    )
    print("[info] hypothetical proteins:", len(hyp_meta))

    diamond_tsv = outdir / "hypothetical_vs_reference.tsv"

    if not args.skip_diamond:
        run_diamond(
            args.diamond,
            hyp_faa,
            args.reference_db,
            diamond_tsv,
            args.threads,
        )

    diamond_all = load_diamond(diamond_tsv)
    print("[info] DIAMOND alignments:", len(diamond_all))

    refmeta = load_reference_metadata(args.reference_metadata)

    hyp_best = classify_hypothetical_hits(
        diamond_all,
        refmeta,
    )

    if hyp_best is not None and not hyp_best.empty:
        print(
            "[info] family-level hypothetical rows:",
            len(hyp_best),
        )
        print(
            "[info] distinct families recovered:",
            hyp_best["family_target"].nunique(),
        )

    ip = load_interpro(args.interpro)
    hyp_best = add_interpro_support(hyp_best, ip)

    interpro_faa = outdir / "interpro_candidates.faa"

    export_interpro_candidates(
        cds,
        hyp_best,
        interpro_faa,
    )

    report = outdir / "Bioremediation_Gene_Miner_Report.xlsx"

    final = make_report(
        report,
        annotated,
        hyp_best,
        hyp_meta,
        diamond_all,
        ip,
        len(cds),
    )

    print("\nDONE")
    print("Report:", report.resolve())
    print("InterPro candidate FASTA:", interpro_faa.resolve())

    if not final.empty:
        print("Final candidate rows:", len(final))

        print(
            "High-confidence:",
            int(
                (
                    (final["Confidence"] == "High")
                    & (final["Decision"] == "Candidate")
                ).sum()
            ),
        )

        print(
            "Microplastic/polymer:",
            int(
                final["Category"]
                .astype(str)
                .str.contains(
                    "Microplastics",
                    case=False,
                    na=False,
                )
                .sum()
            ),
        )

    print(
        "\nReminder: predictions are candidates, "
        "not experimental confirmation."
    )


if __name__ == "__main__":
    main()
