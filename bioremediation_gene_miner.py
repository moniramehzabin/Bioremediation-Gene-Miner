#!/usr/bin/env python3
"""
Bioremediation Gene Miner v4.0

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
import time
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from collections import defaultdict
from pathlib import Path

import pandas as pd

try:
    from Bio import SeqIO
except Exception:
    SeqIO = None

VERSION = "4.0"


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

        # Build several matching variants of the annotated gene name.
        # This makes Prokka-style copy/version suffixes transparent to the
        # rule engine without altering the original gene name in the report.
        #
        # Examples:
        #   azoR2_1 -> azoR2_1, azoR2, azoR
        #   catA_2  -> catA_2, catA
        #   copA1   -> copA1, copA
        gene_original = str(r["gene"] or "").strip()
        gene_variants = []

        if gene_original:
            gene_variants.append(gene_original)

            # Remove Prokka copy-number suffix: _1, _2, _3, ...
            gene_no_copy = re.sub(r"_\d+$", "", gene_original)
            if gene_no_copy and gene_no_copy not in gene_variants:
                gene_variants.append(gene_no_copy)

            # Also expose the family root when a terminal number is attached
            # directly to the gene name, e.g. azoR2 -> azoR, copA1 -> copA.
            gene_family_root = re.sub(r"\d+$", "", gene_no_copy)
            if gene_family_root and gene_family_root not in gene_variants:
                gene_variants.append(gene_family_root)

        gene_text = " | ".join([r["locus_tag"]] + gene_variants)

        # Preserve the original product annotation, but also create a
        # punctuation-normalized variant for matching. This allows equivalent
        # forms such as "catechol-2,3-dioxygenase",
        # "catechol_2,3_dioxygenase", and "catechol 2,3 dioxygenase"
        # to be recognized without changing the annotation shown in reports.
        product_variants = [product]
        product_normalized = re.sub(
            r"[-_\u2010\u2011\u2012\u2013\u2014]+",
            " ",
            product,
        )
        product_normalized = re.sub(
            r"\s+",
            " ",
            product_normalized,
        ).strip()

        if product_normalized and product_normalized != product:
            product_variants.append(product_normalized)

        product_text = " | ".join(
            product_variants + [r["EC_number"], r["db_xref"]]
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



def extract_unmatched_annotated(cds, annotated_hits, out_faa):
    """
    Export non-hypothetical proteins that were not recovered by the
    annotation-rule screen.

    These proteins can then be searched against the curated reference DB
    using the existing DIAMOND family-classification workflow.
    """
    matched_loci = set()

    if annotated_hits is not None and not annotated_hits.empty:
        matched_loci = set(
            annotated_hits["locus_tag"].astype(str)
        )

    rows = []

    with open(out_faa, "w", encoding="utf-8") as f:
        for r in cds:
            product = r["product"] or ""

            if re.search(
                r"hypothetical protein|uncharacterized protein|unknown protein",
                product,
                re.I,
            ):
                continue

            locus_tag = str(r["locus_tag"])

            if locus_tag in matched_loci:
                continue

            seq = r["translation"]
            if not seq:
                continue

            rows.append({
                "locus_tag": locus_tag,
                "gene": r["gene"],
                "product": product,
                "contig": r["contig"],
                "protein_length_aa": len(seq),
            })

            f.write(
                f">{locus_tag} | original_product={product} "
                f"| gene={r['gene']} | contig={r['contig']}\n"
            )
            for i in range(0, len(seq), 80):
                f.write(seq[i:i + 80] + "\n")

    return pd.DataFrame(rows)


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
    "subject_start",
    "subject_end",
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
        "evalue", "bitscore", "qcovhsp", "sstart", "send",
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

    df = pd.read_csv(path, sep="\t", header=None)
    if df.shape[1] != len(DIAMOND_COLS):
        raise ValueError(
            "DIAMOND output lacks reference coordinates. "
            "Rerun DIAMOND into a new output directory; "
            "do not reuse old nine-column files."
        )
    df.columns = DIAMOND_COLS

    numeric_cols = [
        "identity_pct",
        "alignment_len",
        "query_len",
        "subject_len",
        "evalue",
        "bitscore",
        "query_coverage_pct",
        "subject_start",
        "subject_end",
    ]
    for c in numeric_cols:
        df[c] = pd.to_numeric(df[c], errors="coerce")

    valid = (
        df["subject_len"].gt(0)
        & df["subject_start"].between(1, df["subject_len"])
        & df["subject_end"].between(1, df["subject_len"])
    )
    if not valid.all():
        raise ValueError("Invalid or missing reference alignment coordinates.")
    df["reference_coverage_pct"] = (
        (df["subject_end"] - df["subject_start"]).abs().add(1)
        / df["subject_len"] * 100
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
    # Zero is a valid E-value; do not replace it with a default.
    try:
        evalue = float(row.get("evalue"))
    except (TypeError, ValueError, OverflowError):
        evalue = float("inf")
    if evalue < 0 or evalue != evalue:
        evalue = float("inf")
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

    # Sequence strength cannot exceed the confidence supported by the
    # curated biological evidence class.
    def apply_evidence_ceiling(r):
        confidence = str(r.get("confidence", "Weak"))
        evidence = str(r.get("Evidence_classes", "")).lower()
        defaults = str(r.get("Default_confidences", "")).lower()

        if (
            "low" in defaults
            and "moderate" not in defaults
            and "high" not in defaults
        ):
            return "Weak"

        if "direct" not in evidence and confidence == "High":
            return "Moderate"

        return confidence

    merged["confidence"] = merged.apply(apply_evidence_ceiling, axis=1)

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

        message = (
            "Sequence-supported family-level candidate; conserved-domain "
            "validation is recommended before a strong functional call."
        )
        reference_note = str(r.get("Interpretation_rules", "") or "").strip()
        if reference_note:
            message += (
                " Curated reference scope/caveat (not experimental "
                "validation of the query): " + reference_note
            )
        return message

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
# Required InterPro workflow stage
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
    Attach InterPro evidence without changing Gene Miner's original
    candidate, confidence, or decision logic.

    Raw one-row-per-hit evidence is retained separately. Common member
    databases are also summarized into dedicated columns.
    """
    if best is None or best.empty:
        return best

    best = best.copy()
    defaults = {
        "interpro_support": "",
        "interpro_accessions": "",
        "interpro_analyses": "",
        "interpro_detected": False,
        "interpro_status": "NOT_SUBMITTED",
        "panther_hits": "",
        "pfam_hits": "",
        "cdd_hits": "",
        "ncbifam_hits": "",
        "prints_hits": "",
        "gene3d_hits": "",
        "superfamily_hits": "",
        "smart_hits": "",
        "prosite_hits": "",
        "other_interpro_member_hits": "",
        "integrated_interpro_entries": "",
        "interpro_go_terms": "",
        "interpro_pathways": "",
        "interpro_coordinates": "",
    }
    for col, default in defaults.items():
        best[col] = default

    if ip is None or ip.empty:
        return best

    descriptions = defaultdict(list)
    accessions = defaultdict(list)
    analyses = defaultdict(list)
    db_hits = defaultdict(lambda: defaultdict(list))
    integrated = defaultdict(list)
    go_terms = defaultdict(list)
    pathways = defaultdict(list)
    coordinates = defaultdict(list)

    def add_unique(d, key, value):
        value = str(value or "").strip()
        if value and value != "-" and value not in d[key]:
            d[key].append(value)

    def norm_db(value):
        return re.sub(r"[^a-z0-9]+", "", str(value or "").lower())

    db_map = {
        "panther": "panther_hits",
        "pfam": "pfam_hits",
        "cdd": "cdd_hits",
        "ncbifam": "ncbifam_hits",
        "prints": "prints_hits",
        "gene3d": "gene3d_hits",
        "cathgene3d": "gene3d_hits",
        "superfamily": "superfamily_hits",
        "smart": "smart_hits",
        "prosite": "prosite_hits",
        "prositepatterns": "prosite_hits",
        "prositeprofiles": "prosite_hits",
    }

    for _, r in ip.iterrows():
        protein = str(r.get("protein", "") or "").strip()
        if not protein:
            continue

        analysis = str(r.get("analysis", "") or "").strip()
        sig_acc = str(r.get("signature_accession", "") or "").strip()
        sig_desc = str(r.get("signature_description", "") or "").strip()
        ipr_acc = str(r.get("interpro_accession", "") or "").strip()
        ipr_desc = str(r.get("interpro_description", "") or "").strip()
        hit_start = str(r.get("start", "") or "").strip()
        hit_stop = str(r.get("stop", "") or "").strip()
        go = str(r.get("go_terms", "") or "").strip()
        pathway = str(r.get("pathways", "") or "").strip()

        add_unique(descriptions, protein, sig_desc)
        add_unique(descriptions, protein, ipr_desc)
        add_unique(accessions, protein, sig_acc)
        add_unique(accessions, protein, ipr_acc)
        add_unique(analyses, protein, analysis)

        target = db_map.get(norm_db(analysis), "other_interpro_member_hits")
        member_hit = " | ".join(
            x for x in (sig_acc, sig_desc, ipr_acc, ipr_desc)
            if x and x != "-"
        )
        if member_hit and member_hit not in db_hits[protein][target]:
            db_hits[protein][target].append(member_hit)

        if ipr_acc and ipr_acc != "-":
            add_unique(
                integrated, protein,
                " | ".join(x for x in (ipr_acc, ipr_desc) if x and x != "-")
            )
        add_unique(go_terms, protein, go)
        add_unique(pathways, protein, pathway)
        if hit_start and hit_stop:
            add_unique(
                coordinates, protein,
                f"{analysis}:{sig_acc or 'NA'}:{hit_start}-{hit_stop}"
            )

    best["interpro_support"] = best["query"].map(
        lambda q: "; ".join(descriptions.get(str(q), []))
    )
    best["interpro_accessions"] = best["query"].map(
        lambda q: "; ".join(accessions.get(str(q), []))
    )
    best["interpro_analyses"] = best["query"].map(
        lambda q: "; ".join(analyses.get(str(q), []))
    )
    best["interpro_detected"] = best["query"].map(
        lambda q: bool(
            descriptions.get(str(q)) or accessions.get(str(q)) or analyses.get(str(q))
        )
    )
    best.loc[
        best["interpro_detected"], "interpro_status"
    ] = "HITS_FOUND"

    for col in [
        "panther_hits", "pfam_hits", "cdd_hits", "ncbifam_hits",
        "prints_hits", "gene3d_hits", "superfamily_hits", "smart_hits",
        "prosite_hits", "other_interpro_member_hits",
    ]:
        best[col] = best["query"].map(
            lambda q, c=col: "; ".join(db_hits.get(str(q), {}).get(c, []))
        )

    best["integrated_interpro_entries"] = best["query"].map(
        lambda q: "; ".join(integrated.get(str(q), []))
    )
    best["interpro_go_terms"] = best["query"].map(
        lambda q: "; ".join(go_terms.get(str(q), []))
    )
    best["interpro_pathways"] = best["query"].map(
        lambda q: "; ".join(pathways.get(str(q), []))
    )
    best["interpro_coordinates"] = best["query"].map(
        lambda q: "; ".join(coordinates.get(str(q), []))
    )

    best.loc[
        best["interpro_detected"], "evidence_source"
    ] = "DIAMOND + InterPro (independent evidence)"
    return best


# ============================================================
# InterPro evidence relationship resolver
# ============================================================
# This resolver ONLY summarizes how independent InterPro evidence relates
# to the DIAMOND candidate. It does not assign pathways, reactions, or a
# forced final protein function.
RESOLVER_GROUPS = {
    "p450": ["cytochrome p450", "cyt_p450", "pf00067", "ipr001128", "ipr002397"],
    "intradiol_dioxygenase": ["intradiol", "ring-cleavage dioxygenase", "pf00775", "ipr000627", "ipr015889"],
    "multicopper": ["multi-copper", "multicopper", "laccase", "pf02578", "ipr038371", "ipr011324"],
    "sdr": ["short-chain dehydrogenase", "short chain dehydrogenase", "sdr family", "rossmann", "pf00106", "ipr002347", "ipr036291"],
    "p_loop_atpase": ["p-loop", "partitioning atpase", "parab", "aaa domain", "soj", "ipr027417", "ipr025669"],
    "fmn_reductase": ["fmn reductase", "flavin reductase", "flavoprotein-like", "pf03358", "ipr005025", "ipr029039", "pf01613", "ipr002563"],
}
FAMILY_EXPECTATIONS = {
    "cytochrome p450": {"p450"}, "gcoa": {"p450"},
    "cata": {"intradiol_dioxygenase"},
    "laccase": {"multicopper"}, "laccase plastic-associated": {"multicopper"},
    "linb": {"sdr"}, "lina": {"sdr"},
    "chrr": {"fmn_reductase"}, "chromate reductase": {"fmn_reductase"},
    "flavin reductase": {"fmn_reductase"},
    "azoreductase": {"fmn_reductase"},
    "arsa": {"p_loop_atpase"},
}
BROAD_RELATED = {"linb", "lina"}

def _resolver_norm(x):
    return re.sub(r"[^a-z0-9]+", " ", str(x or "").lower()).strip()

def _resolver_groups(row):
    blob = "; ".join([
        str(row.get("interpro_support", "")),
        str(row.get("interpro_accessions", "")),
        str(row.get("panther_hits", "")),
        str(row.get("pfam_hits", "")),
        str(row.get("cdd_hits", "")),
        str(row.get("ncbifam_hits", "")),
        str(row.get("prints_hits", "")),
        str(row.get("gene3d_hits", "")),
        str(row.get("superfamily_hits", "")),
        str(row.get("smart_hits", "")),
        str(row.get("prosite_hits", "")),
    ]).lower()
    return {g for g, terms in RESOLVER_GROUPS.items() if any(t in blob for t in terms)}

def resolve_interpro_evidence(best):
    """Evidence-relationship summary only; never changes DIAMOND confidence/decision."""
    if best is None or best.empty:
        return best

    best = best.copy()
    best["resolver_status"] = "UNINFORMATIVE"
    best["resolver_note"] = ""

    for i, r in best.iterrows():
        if not bool(r.get("interpro_detected", False)):
            best.at[i, "resolver_note"] = "No InterPro evidence was returned for this candidate."
            continue

        fam = _resolver_norm(r.get("family_target", ""))
        groups = _resolver_groups(r)
        expected = set()
        for key, vals in FAMILY_EXPECTATIONS.items():
            if fam == key or fam.startswith(key + " "):
                expected |= vals

        matched = expected & groups
        support = str(r.get("interpro_support", "")).lower()

        # Known strong alternative for ArsA-like weak DIAMOND hits.
        if fam == "arsa" and "p_loop_atpase" in groups and any(
            x in support for x in ["partitioning atpase", "parab", "sporulation initiation inhibitor soj"]
        ):
            best.at[i, "resolver_status"] = "CONFLICTING"
            best.at[i, "resolver_note"] = (
                "InterPro returns a ParAB/Soj-like P-loop ATPase interpretation rather "
                "than evidence specifically consistent with the DIAMOND candidate."
            )

        # Compatible but only broad architecture.
        elif matched and fam in BROAD_RELATED:
            best.at[i, "resolver_status"] = "BROAD/RELATED"
            best.at[i, "resolver_note"] = (
                "InterPro supports a broader protein-family/domain architecture related "
                "to the DIAMOND candidate."
            )

        # Mixed laccase/multicopper plus YfiH/CNF1-like evidence.
        elif matched and fam.startswith("laccase") and any(
            x in support for x in ["yfih", "cnf1", "cysteine hydrolase", "peptidoglycan editing"]
        ):
            best.at[i, "resolver_status"] = "BROAD/RELATED"
            best.at[i, "resolver_note"] = (
                "InterPro contains multicopper/laccase-related evidence together with "
                "alternative YfiH/CNF1-like annotations."
            )

        elif matched:
            best.at[i, "resolver_status"] = "SUPPORTING"
            best.at[i, "resolver_note"] = (
                "Independent InterPro family/domain evidence is consistent with the "
                "DIAMOND candidate."
            )

        elif expected and groups:
            best.at[i, "resolver_status"] = "CONFLICTING"
            best.at[i, "resolver_note"] = (
                "InterPro returns recognizable family/domain evidence that is not "
                "consistent with the expected architecture of the DIAMOND candidate."
            )

        elif groups:
            best.at[i, "resolver_status"] = "BROAD/RELATED"
            best.at[i, "resolver_note"] = (
                "InterPro provides related or broader family/domain information but "
                "does not directly support the DIAMOND candidate."
            )

        else:
            best.at[i, "resolver_note"] = (
                "InterPro evidence was returned, but it is not informative enough for "
                "a simple relationship summary."
            )

    return best

# ============================================================
# Automated InterProScan via EMBL-EBI Job Dispatcher
# ============================================================
INTERPRO_REST_BASE = "https://www.ebi.ac.uk/Tools/services/rest/iprscan5"

def _http_text(url, data=None, timeout=60):
    req = urllib.request.Request(url, data=data, headers={
        "User-Agent": "Bioremediation-Gene-Miner/4.0"
    })
    with urllib.request.urlopen(req, timeout=timeout) as response:
        return response.read().decode("utf-8", errors="replace")

def _interpro_result_types(job_id):
    root = ET.fromstring(_http_text(
        f"{INTERPRO_REST_BASE}/resulttypes/{job_id}"
    ))
    found = []
    for elem in root.iter():
        if elem.tag.split("}")[-1] == "type":
            item = {c.tag.split("}")[-1]: (c.text or "").strip() for c in elem}
            if item.get("identifier"):
                found.append(item)
    return found

def run_interpro_web(fasta_path, email, outdir, poll_seconds=5, timeout_minutes=45):
    import importlib.util
    adapter_path = Path(__file__).resolve().with_name("geneminer_interpro_batch.py")
    spec = importlib.util.spec_from_file_location("geneminer_v6_batch_runtime", adapter_path)
    adapter = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(adapter)
    return adapter.run_interpro_batch(fasta_path, email, outdir, poll_seconds, timeout_minutes)

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
    unmatched_best,
    unmatched_meta,
    hyp_best,
    hyp_meta,
    diamond_all,
    interpro,
    total_cds,
    pet_controls=None,
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
                "Sequence match strength": r["confidence"],
                "Match_strength": f"{r['confidence']} match",
                "Decision": "Candidate",
                "Evidence_source": r["evidence_source"],
                "Identity_pct": "",
                "Query_coverage_pct": "",
                "Reference_coverage_pct": "",
                "Evalue": "",
                "Bitscore": "",
                "Best_reference": "",
                "InterPro_status": "NOT_SUBMITTED",
                "InterPro_analyses": "",
                "InterPro_accessions": "",
                "InterPro_domain_support": "",
                "PANTHER_hits": "",
                "Pfam_hits": "",
                "CDD_hits": "",
                "NCBIfam_hits": "",
                "PRINTS_hits": "",
                "Gene3D_hits": "",
                "SUPERFAMILY_hits": "",
                "SMART_hits": "",
                "PROSITE_hits": "",
                "Other_InterPro_member_hits": "",
                "Integrated_InterPro_entries": "",
                "InterPro_GO_terms": "",
                "InterPro_pathway_records": 0,
                "InterPro_pathway_evidence": "",
                "InterPro_coordinates": "",
                "Resolver_status": "",
                "Resolver_note": "",
                "Interpretation": r["interpretation"],
            })

    if unmatched_best is not None and not unmatched_best.empty:
        unmatched_lookup = {}

        if (
            unmatched_meta is not None
            and not unmatched_meta.empty
            and "locus_tag" in unmatched_meta.columns
        ):
            unmatched_lookup = {
                str(r["locus_tag"]): r
                for _, r in unmatched_meta.iterrows()
            }

        for _, r in unmatched_best.iterrows():
            locus_tag = str(r["query"])
            meta = unmatched_lookup.get(locus_tag, {})

            final_rows.append({
                "Locus_tag": locus_tag,
                "Gene": meta.get("gene", ""),
                "Family_target": r.get("family_target", ""),
                "Product_or_prediction": meta.get("product", ""),
                "Origin": "Annotated → sequence rescued",
                "Category": r.get("Headline_category", "Unresolved"),
                "Evidence_class": r.get("Evidence_classes", ""),
                "Sequence match strength": r["confidence"],
                "Match_strength": f"{r['confidence']} match",
                "Decision": r["decision"],
                "Evidence_source": (
                    str(r["evidence_source"])
                    + " (annotated-unmatched rescue)"
                ),
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
                "InterPro_status": r.get("interpro_status", "NOT_SUBMITTED"),
                "InterPro_analyses": r.get("interpro_analyses", ""),
                "InterPro_accessions": r.get("interpro_accessions", ""),
                "InterPro_domain_support": r.get("interpro_support", ""),
                "PANTHER_hits": r.get("panther_hits", ""),
                "Pfam_hits": r.get("pfam_hits", ""),
                "CDD_hits": r.get("cdd_hits", ""),
                "NCBIfam_hits": r.get("ncbifam_hits", ""),
                "PRINTS_hits": r.get("prints_hits", ""),
                "Gene3D_hits": r.get("gene3d_hits", ""),
                "SUPERFAMILY_hits": r.get("superfamily_hits", ""),
                "SMART_hits": r.get("smart_hits", ""),
                "PROSITE_hits": r.get("prosite_hits", ""),
                "Other_InterPro_member_hits": r.get("other_interpro_member_hits", ""),
                "Integrated_InterPro_entries": r.get("integrated_interpro_entries", ""),
                "InterPro_GO_terms": r.get("interpro_go_terms", ""),
                "InterPro_pathway_records": 0 if not str(r.get("interpro_pathways", "") or "").strip()
                    else len([x for x in str(r.get("interpro_pathways", "")).split("; ") if x.strip()]),
                "InterPro_pathway_evidence": "See InterPro_Full_Report.xlsx and interpro_web.tsv"
                    if str(r.get("interpro_pathways", "") or "").strip() else "",
                "InterPro_coordinates": r.get("interpro_coordinates", ""),
                "Resolver_status": r.get("resolver_status", "UNINFORMATIVE"),
                "Resolver_note": r.get("resolver_note", ""),
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
                "Sequence match strength": r["confidence"],
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
                "InterPro_status": r.get("interpro_status", "NOT_SUBMITTED"),
                "InterPro_analyses": r.get("interpro_analyses", ""),
                "InterPro_accessions": r.get("interpro_accessions", ""),
                "InterPro_domain_support": r.get("interpro_support", ""),
                "PANTHER_hits": r.get("panther_hits", ""),
                "Pfam_hits": r.get("pfam_hits", ""),
                "CDD_hits": r.get("cdd_hits", ""),
                "NCBIfam_hits": r.get("ncbifam_hits", ""),
                "PRINTS_hits": r.get("prints_hits", ""),
                "Gene3D_hits": r.get("gene3d_hits", ""),
                "SUPERFAMILY_hits": r.get("superfamily_hits", ""),
                "SMART_hits": r.get("smart_hits", ""),
                "PROSITE_hits": r.get("prosite_hits", ""),
                "Other_InterPro_member_hits": r.get("other_interpro_member_hits", ""),
                "Integrated_InterPro_entries": r.get("integrated_interpro_entries", ""),
                "InterPro_GO_terms": r.get("interpro_go_terms", ""),
                "InterPro_pathway_records": (
                    0 if not str(r.get("interpro_pathways", "") or "").strip()
                    else len([x for x in str(r.get("interpro_pathways", "")).split("; ") if x.strip()])
                ),
                "InterPro_pathway_evidence": (
                    "See InterPro_Full_Report.xlsx and interpro_web.tsv"
                    if str(r.get("interpro_pathways", "") or "").strip() else ""
                ),
                "InterPro_coordinates": r.get("interpro_coordinates", ""),
                "Resolver_status": r.get("resolver_status", "UNINFORMATIVE"),
                "Resolver_note": r.get("resolver_note", ""),
                "Interpretation": (
                    r["interpretation"]
                    + (
                        " Independent InterPro evidence is present; the original DIAMOND "
                        "confidence and decision are preserved unchanged."
                        if bool(r.get("interpro_detected", False))
                        else ""
                    )
                ),
            })

    final = pd.DataFrame(final_rows)
    if not final.empty and pet_controls is not None and not pet_controls.empty:
        notes = {
            "REVIEW_EXACT_ASSAY_NEGATIVE_CONTROL_SEQUENCE":
                "Exact assay-negative sequence match. No PET activity was detected under the tested conditions.",
            "REVIEW_ASSAY_NEGATIVE_SEQUENCE_WITH_X":
                "Matches an assay-negative sequence record containing an unresolved residue (X). "
                "No PET activity was detected under the tested conditions.",
            "REVIEW_CONFLICTING_EXACT_CONTROL_LABELS":
                "Matching sequence records have conflicting PET assay labels.",
        }
        review_by_id = pet_controls.set_index("Locus_tag")["PET_control_flag"].map(notes)
        reasons = final["Locus_tag"].astype(str).map(review_by_id)
        families = final["Family_target"].fillna("").str.strip().str.casefold()
        is_pet = families.eq("petase") | families.str.startswith("pet hydrolase")
        review = is_pet & reasons.notna() & final["Decision"].eq("Candidate")
        final.loc[review, "Decision"] = "Review"
        final.loc[review, "Interpretation"] = (
            final.loc[review, "Interpretation"].fillna("").astype(str).str.strip()
            + " " + reasons[review]
        ).str.strip()

    if not final.empty:
        final["_rank"] = (
            final["Sequence match strength"].map(CONF_ORDER).fillna(9)
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
            (final["Sequence match strength"] == "High")
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
            final["InterPro_status"].eq("HITS_FOUND")
        ].copy()
    else:
        high = final.copy()
        micro = final.copy()
        supporting = final.copy()
        review = final.copy()
        family_candidates = final.copy()
        interpro_supported = final.copy()

    unmatched_queries = (
        0
        if unmatched_best is None or unmatched_best.empty
        else unmatched_best["query"].nunique()
    )

    unmatched_family_hits = (
        0
        if unmatched_best is None or unmatched_best.empty
        else len(unmatched_best)
    )

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
            ["Annotated unmatched proteins screened", len(unmatched_meta)],
            ["Annotated-unmatched queries with DIAMOND hits", unmatched_queries],
            ["Family-level annotated rescue rows", unmatched_family_hits],
            ["Hypothetical proteins screened", len(hyp_meta)],
            ["Hypothetical queries with DIAMOND hits", hyp_queries],
            ["Family-level hypothetical hit rows", family_hits],
            ["InterPro rows with hits", len(interpro_supported)],
            ["High-confidence final candidates", len(high)],
            ["Review/rejected final rows", len(review)],
        ],
        columns=["Metric", "Value"],
    )

    # Main workbook: concise InterPro evidence only.
    # Full evidence remains in InterPro_Full_Report.xlsx + interpro_web.tsv.
    verbose_main_columns = [
        "InterPro_analyses", "InterPro_accessions", "InterPro_domain_support",
        "PRINTS_hits", "Gene3D_hits", "SUPERFAMILY_hits", "SMART_hits",
        "PROSITE_hits", "Other_InterPro_member_hits",
        "InterPro_pathway_records", "InterPro_pathway_evidence",
        "InterPro_coordinates",
    ]

    def _compact_main(df):
        if df is None:
            return df
        return df.drop(columns=verbose_main_columns, errors="ignore")

    final_main = _compact_main(final)
    high_main = _compact_main(high)
    family_candidates_main = _compact_main(family_candidates)
    interpro_supported_main = _compact_main(interpro_supported)
    supporting_main = _compact_main(supporting)
    review_main = _compact_main(review)
    with pd.ExcelWriter(out_xlsx, engine="openpyxl") as writer:
        summary.to_excel(
            writer,
            sheet_name="Summary",
            index=False,
        )

        final_main.to_excel(
            writer,
            sheet_name="All_candidates",
            index=False,
        )

        high_main.to_excel(
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
            hyp_best_excel = hyp_best.drop(
                columns=[
                    "interpro_detected",
                    "interpro_pathways", "interpro_analyses",
                    "interpro_accessions", "interpro_support",
                    "prints_hits", "gene3d_hits", "superfamily_hits",
                    "smart_hits", "prosite_hits",
                    "other_interpro_member_hits", "interpro_coordinates",
                ],
                errors="ignore",
            )
            hyp_best_excel.to_excel(
                writer,
                sheet_name="Predicted_hypothetical",
                index=False,
            )

        family_candidates_main.to_excel(
            writer,
            sheet_name="Family_level_candidates",
            index=False,
        )

        interpro_supported_main.to_excel(
            writer,
            sheet_name="InterPro_supported",
            index=False,
        )


        micro_main = _compact_main(micro)
        micro_main.to_excel(
            writer,
            sheet_name="Microplastic_candidates",
            index=False,
        )

        supporting_main.to_excel(
            writer,
            sheet_name="Supporting_functions",
            index=False,
        )

        review_main.to_excel(
            writer,
            sheet_name="Review_required",
            index=False,
        )

        diamond_all.to_excel(
            writer,
            sheet_name="All_DIAMOND_hits",
            index=False,
        )

        # Full raw/detailed InterPro evidence is written separately to
        # InterPro_Full_Report.xlsx and preserved in interpro_web.tsv.
        # Keeping it out of the main workbook prevents Excel cell overflow.

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



def export_interpro_candidates(cds, hyp_best, out_faa, unmatched_best=None):
    """Export every Review protein from both DIAMOND routes."""
    selected = []
    wanted = set()

    for route, frame in [
        ("Hypothetical", hyp_best),
        ("Annotated sequence rescue", unmatched_best),
    ]:
        if frame is None or frame.empty:
            continue
        for _, row in frame.iterrows():
            if str(row.get("decision", "")).strip() != "Review":
                continue
            query = str(row.get("query", "")).strip()
            if not query:
                raise ValueError("Review row has no protein identifier.")
            wanted.add(query)
            selected.append({
                "Protein": query,
                "Route": route,
                "Family": row.get("family_target", ""),
                "Sequence_match_strength": row.get("confidence", ""),
                "Decision": "Review",
            })

    seqmap = {}
    for row in cds:
        locus = str(row.get("locus_tag", ""))
        sequence = str(row.get("translation", "") or "").strip()
        if locus in seqmap and sequence and seqmap[locus] != sequence:
            raise ValueError("Conflicting sequences for identifier: " + locus)
        if sequence:
            seqmap[locus] = sequence

    missing = sorted(wanted - set(seqmap))
    if missing:
        raise ValueError(
            "Review proteins lack sequences: " + ", ".join(missing)
        )

    with open(out_faa, "w", encoding="utf-8") as handle:
        for query in sorted(wanted):
            handle.write(">" + query + "\n")
            sequence = seqmap[query]
            for start in range(0, len(sequence), 80):
                handle.write(sequence[start:start + 80] + "\n")

    with Path(out_faa).with_suffix(".selection.tsv").open(
        "w", encoding="utf-8", newline=""
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "Protein", "Route", "Family",
                "Sequence_match_strength", "Decision"
            ],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(selected)


def attach_interpro_to_route(best, ip, searched):
    """Attach evidence without changing prediction decisions or strength."""
    if best is None or best.empty:
        return best

    before = best[["decision", "confidence"]].copy()
    result = add_interpro_support(best, ip)
    result["interpro_status"] = "NOT_SUBMITTED"
    result.loc[
        result["query"].astype(str).isin(searched), "interpro_status"
    ] = "NO_HITS"
    result.loc[
        result["interpro_detected"].eq(True), "interpro_status"
    ] = "HITS_FOUND"

    result = resolve_interpro_evidence(result)
    not_submitted = result["interpro_status"].eq("NOT_SUBMITTED")
    result.loc[not_submitted, "resolver_status"] = ""
    result.loc[not_submitted, "resolver_note"] = ""

    no_hits = result["interpro_status"].eq("NO_HITS")
    result.loc[no_hits, "resolver_status"] = "UNINFORMATIVE"
    result.loc[no_hits, "resolver_note"] = (
        "Protein was searched by InterPro, but no usable hits were returned."
    )

    if not result[["decision", "confidence"]].equals(before):
        raise RuntimeError("InterPro attachment changed prediction decisions.")
    return result



# ============================================================
# CLI
# ============================================================


def apply_merA_sequence_safeguard(best, fasta_path):
    """Apply the existing MerA sequence check to either DIAMOND route."""
    if best is None or best.empty:
        return best

    best = best.copy()
    sequences = {
        record.id: str(record.seq).upper()
        for record in SeqIO.parse(str(fasta_path), "fasta")
    }

    for idx, row in best.iterrows():
        if str(row.get("family_target", "")).strip().lower() != "mera":
            continue

        # Preserve rejection safeguards; never promote rejected matches.
        if str(row.get("decision", "")).strip().lower() == "reject":
            continue

        sequence = sequences.get(str(row.get("query", "")), "")
        catalytic_motif = bool(re.search(r"CVN[IV]GC", sequence))
        terminal_cys_pair = "CC" in sequence[-10:]

        if not (catalytic_motif and terminal_cys_pair):
            best.at[idx, "decision"] = "Review"
            existing = str(row.get("interpretation", "") or "").strip()
            note = (
                "MerA sequence safeguard: the existing catalytic-motif "
                "and terminal cysteine-pair criteria were not both met. "
                "Retain for review; do not accept a MerA assignment from "
                "similarity alone. Failure of this screen does not prove "
                "absence of mercury-reduction activity."
            )
            best.at[idx, "interpretation"] = (
                existing + " " + note
            ).strip()

    return best


def _pet_control_flag(active_score, negative_score, exact_labels):
    if len(exact_labels) > 1:
        return "REVIEW_CONFLICTING_EXACT_CONTROL_LABELS"
    if exact_labels:
        if "ACTIVE_IN_ASSAY_PANEL" in exact_labels:
            return "EXACT_ACTIVE_CONTROL_SEQUENCE"
        return "REVIEW_EXACT_ASSAY_NEGATIVE_CONTROL_SEQUENCE"
    if negative_score > active_score:
        return "REVIEW_NEGATIVE_CONTROL_SIMILARITY"
    if negative_score > 0:
        return "REVIEW_BOTH_CONTROL_CLASSES_MATCH"
    if active_score > 0:
        return "ACTIVE_CONTROL_SIMILARITY_ONLY"
    return "NO_QUALIFYING_CONTROL_MATCH"


def _pet_control_evidence(cds, annotated, hyp_best, unmatched_best,
                          outdir, diamond, threads):
    import hashlib
    import json
    import shutil

    root = Path(__file__).resolve().parent / "pet_activity_controls"
    if not root.is_dir():
        raise RuntimeError("Connected PET control database is missing.")

    manifest = json.loads(
        (root / "MANIFEST.json").read_text(encoding="utf-8")
    )
    for name, expected in manifest.items():
        actual = hashlib.sha256((root / name).read_bytes()).hexdigest()
        if actual != expected:
            raise RuntimeError("PET control file changed: " + name)

    # Select existing PET-related candidates from all three routes.
    selected = set()
    for frame, id_col, family_col in [
        (annotated, "locus_tag", "family"),
        (hyp_best, "query", "family_target"),
        (unmatched_best, "query", "family_target"),
    ]:
        if frame is None or frame.empty:
            continue
        if id_col not in frame or family_col not in frame:
            raise RuntimeError("Unexpected PET candidate-table structure.")
        for _, row in frame.iterrows():
            family = str(row[family_col])
            if re.search(r"\b(?:PET(?:ase)?|lipase|cutinase)\b", family, re.I):
                if str(row.get("decision", "")).lower() != "reject":
                    selected.add(str(row[id_col]))

    columns = [
        "Locus_tag", "PET_control_flag",
        "Active_control", "Active_source_ID", "Active_similarity",
        "Assay_negative_control", "Assay_negative_source_ID",
        "Assay_negative_similarity", "Exact_control_labels", "Scope",
    ]
    output = Path(outdir) / "PET_control_comparison.tsv"
    if not selected:
        result = pd.DataFrame(columns=columns)
        result.to_csv(output, sep="\t", index=False)
        return result

    sequences = {}
    for row in cds:
        name = str(row.get("locus_tag", ""))
        if name not in selected:
            continue
        sequence = "".join(
            str(row.get("translation", "")).split()
        ).upper().rstrip("*")
        if not sequence:
            raise RuntimeError("PET query lacks a translation: " + name)
        if name in sequences and sequences[name] != sequence:
            raise RuntimeError("Conflicting PET query identifier: " + name)
        sequences[name] = sequence

    if set(sequences) != selected:
        raise RuntimeError("Not every PET query has a CDS translation.")

    aliases = {f"PETQUERY{i:06d}": name
               for i, name in enumerate(sorted(sequences), 1)}
    query = Path(outdir) / "pet_control_queries.faa"
    query.write_text(
        "".join(">" + alias + "\n" + sequences[name] + "\n"
                for alias, name in aliases.items()),
        encoding="ascii"
    )

    executable = (
        str(diamond) if Path(str(diamond)).is_file()
        else shutil.which(str(diamond))
    )
    if not executable:
        fallback = Path(
            r"F:\Tools\EBPT_USER_TEST\diamond-windows\diamond.exe"
        )
        executable = str(fallback) if fallback.is_file() else None
    if not executable:
        raise RuntimeError("DIAMOND unavailable for PET control comparison.")

    raw = Path(outdir) / "pet_control_alignments.tsv"
    command = [
        executable, "blastp", "--query", str(query),
        "--db", str(root / "controls.dmnd"), "--out", str(raw),
        "--outfmt", "6", "qseqid", "sseqid", "pident",
        "qcovhsp", "slen", "sstart", "send", "evalue",
        "--max-target-seqs", "0", "--max-hsps", "1",
        "--evalue", "1e-5", "--sensitive",
    ]
    if threads:
        command += ["--threads", str(threads)]
    with (Path(outdir) / "PET_control_search.log").open(
        "w", encoding="utf-8"
    ) as log:
        subprocess.run(command, check=True, stdout=log,
                       stderr=subprocess.STDOUT)

    metadata = pd.read_csv(
        root / "controls.tsv", sep="\t", keep_default_na=False
    ).set_index("ControlID").to_dict("index")
    best = {name: {} for name in selected}
    exact = {name: set() for name in selected}
    exact_with_x = set()

    # Use literal sequence hashes for exact matches.
    for name, sequence in sequences.items():
        digest = hashlib.sha256(sequence.encode("ascii")).hexdigest()
        for control in metadata.values():
            if str(control["SequenceSHA256"]).lower() != digest:
                continue
            exact[name].add(control["AssayLabel"])
            if str(control["ContainsX"]).lower() == "true":
                exact_with_x.add(name)

    with raw.open(encoding="utf-8", newline="") as handle:
        for values in csv.reader(handle, delimiter="\t"):
            if not values:
                continue
            if len(values) != 8:
                raise RuntimeError("Unexpected PET control alignment format.")
            alias, control_id = values[:2]
            name = aliases[alias]
            control = metadata[control_id]
            if str(control["ContainsX"]).lower() == "true":
                continue
            pid, qcov, slen, start, end, evalue = map(float, values[2:])
            if slen <= 0:
                raise RuntimeError("Invalid PET control alignment length.")
            rcov = min(100.0, (abs(end-start)+1) / slen * 100)
            # Broad comparison safeguards; not an activity classifier.
            if pid < 25 or qcov < 50 or rcov < 50 or evalue > 1e-5:
                continue
            score = (pid / 100.0) * (min(qcov, rcov) / 100.0)
            label = control["AssayLabel"]
            old = best[name].get(label, (0.0, "", ""))
            item = (score, control_id, control["SourceID"])
            if score > old[0]:
                best[name][label] = item

    rows = []
    for name in sorted(selected):
        active = best[name].get("ACTIVE_IN_ASSAY_PANEL", (0.0, "", ""))
        negative = best[name].get(
            "NO_ACTIVITY_DETECTED_IN_ASSAY_PANEL", (0.0, "", "")
        )
        flag = _pet_control_flag(active[0], negative[0], exact[name])
        if name in exact_with_x:
            if flag == "REVIEW_EXACT_ASSAY_NEGATIVE_CONTROL_SEQUENCE":
                flag = "REVIEW_ASSAY_NEGATIVE_SEQUENCE_WITH_X"
            elif flag == "EXACT_ACTIVE_CONTROL_SEQUENCE":
                flag = "ASSAY_ACTIVE_SEQUENCE_WITH_X"
        rows.append([
            name, flag,
            active[1], active[2], active[0],
            negative[1], negative[2], negative[0],
            ";".join(sorted(exact[name])),
            "Exact assay-negative or conflicting sequence records route PET calls to Review. "
            "Assay-negative means no detected activity under tested conditions.",
        ])
    result = pd.DataFrame(rows, columns=columns)
    result.to_csv(output, sep="\t", index=False)
    print("[info] PET control comparison:", output)
    return result


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
        help="Optional InterProScan TSV for selected Review candidates",
    )
    ap.add_argument(
        "--interpro-auto",
        action="store_true",
        help="Automatically submit Review candidates to EMBL-EBI InterProScan REST.",
    )
    ap.add_argument("--interpro-email", help="Email required by EMBL-EBI Job Dispatcher.")
    ap.add_argument("--interpro-poll-seconds", type=int, default=5)
    ap.add_argument("--interpro-timeout-minutes", type=int, default=180)

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

    # InterPro is optional enrichment; core reports do not require the service.

    if args.interpro:
        supplied_interpro = Path(args.interpro)
        if not supplied_interpro.is_file():
            ap.error("The supplied InterPro TSV does not exist.")
        # A hit-only TSV cannot prove completion for absent proteins.
        # Require readable, nonempty parsed evidence for this import route.
        supplied_evidence = load_interpro(args.interpro)
        if supplied_evidence.empty:
            ap.error(
                "The supplied InterPro TSV has no readable hit records. "
                "An empty file alone cannot establish a completed search."
            )
    elif args.interpro_auto and not args.interpro_email:
        ap.error("--interpro-auto requires --interpro-email.")


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

    unmatched_faa = outdir / "annotated_unmatched.faa"
    unmatched_meta = extract_unmatched_annotated(
        cds,
        annotated,
        unmatched_faa,
    )

    unmatched_meta.to_csv(
        outdir / "annotated_unmatched.tsv",
        sep="\t",
        index=False,
    )
    print("[info] annotated unmatched proteins:", len(unmatched_meta))

    hyp_faa = outdir / "hypothetical_proteins.faa"
    hyp_meta = extract_hypotheticals(cds, hyp_faa)

    hyp_meta.to_csv(
        outdir / "hypothetical_proteins.tsv",
        sep="\t",
        index=False,
    )
    print("[info] hypothetical proteins:", len(hyp_meta))

    # Sequence rescue for annotated proteins not recognized by annotation rules.
    # This is separate from the existing hypothetical-protein rescue pathway.
    unmatched_diamond_tsv = outdir / "annotated_unmatched_vs_reference.tsv"

    if (
        not args.skip_diamond
        and unmatched_faa.exists()
        and unmatched_faa.stat().st_size > 0
    ):
        run_diamond(
            args.diamond,
            unmatched_faa,
            args.reference_db,
            unmatched_diamond_tsv,
            args.threads,
        )

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

    unmatched_diamond_all = load_diamond(unmatched_diamond_tsv)
    print(
        "[info] annotated-unmatched DIAMOND alignments:",
        len(unmatched_diamond_all),
    )

    refmeta = load_reference_metadata(args.reference_metadata)

    unmatched_best = classify_hypothetical_hits(
        unmatched_diamond_all,
        refmeta,
    )

    unmatched_best = apply_merA_sequence_safeguard(
        unmatched_best, unmatched_faa
    )

    hyp_best = classify_hypothetical_hits(
        diamond_all,
        refmeta,
    )

    hyp_best = apply_merA_sequence_safeguard(hyp_best, hyp_faa)

    if hyp_best is not None and not hyp_best.empty:
        print(
            "[info] family-level hypothetical rows:",
            len(hyp_best),
        )
        print(
            "[info] distinct families recovered:",
            hyp_best["family_target"].nunique(),
        )

    pet_controls = _pet_control_evidence(
        cds, annotated, hyp_best, unmatched_best,
        outdir, args.diamond, args.threads
    )

    interpro_faa = outdir / "interpro_candidates.faa"
    export_interpro_candidates(
        cds, hyp_best, interpro_faa, unmatched_best=unmatched_best
    )

    from geneminer_interpro_fallback import (
        fasta_ids, pending_state, save_status, status_sheet, mark_unresolved,
        optional_interpro, preserve_old_full_report,
    )
    selected_interpro_ids = fasta_ids(interpro_faa)
    interpro_state = pending_state(args, selected_interpro_ids)
    save_status(outdir, interpro_state)

    # Always write usable core results BEFORE waiting for an optional service.
    core_ip = pd.DataFrame()
    core_hyp = attach_interpro_to_route(hyp_best, core_ip, set())
    core_unmatched = attach_interpro_to_route(unmatched_best, core_ip, set())
    core_hyp = mark_unresolved(core_hyp, selected_interpro_ids,
                               interpro_state["Stage"], interpro_state["Note"])
    core_unmatched = mark_unresolved(core_unmatched, selected_interpro_ids,
                                     interpro_state["Stage"], interpro_state["Note"])
    core_report = outdir / "Bioremediation_Gene_Miner_Core_Report.xlsx"
    make_report(core_report, annotated, core_unmatched, unmatched_meta,
                core_hyp, hyp_meta, diamond_all, core_ip, len(cds), pet_controls=pet_controls)
    status_sheet(core_report, interpro_state, pd)
    print("[info] Core report available now:", core_report.resolve())

    interpro_input, interpro_state, submitted_interpro = optional_interpro(
        args, interpro_faa, outdir, run_interpro_web
    )
    preserve_old_full_report(outdir, interpro_state["Stage"])
    ip = load_interpro(interpro_input)
    hyp_best = attach_interpro_to_route(hyp_best, ip, submitted_interpro)
    unmatched_best = attach_interpro_to_route(
        unmatched_best, ip, submitted_interpro
    )
    hyp_best = mark_unresolved(hyp_best, selected_interpro_ids,
                                interpro_state["Stage"], interpro_state["Note"],
        completed=submitted_interpro, protein_states=interpro_state.get("ProteinStates", {}))
    unmatched_best = mark_unresolved(unmatched_best, selected_interpro_ids,
                                      interpro_state["Stage"], interpro_state["Note"],
        completed=submitted_interpro, protein_states=interpro_state.get("ProteinStates", {}))

    report = outdir / "Bioremediation_Gene_Miner_Report.xlsx"

    final = make_report(
        report,
        annotated,
        unmatched_best,
        unmatched_meta,
        hyp_best,
        hyp_meta,
        diamond_all,
        ip,
        len(cds),
        pet_controls=pet_controls,
    )

    status_sheet(report, interpro_state, pd)
    print("[info] InterPro enrichment:", interpro_state["Stage"])

    with pd.ExcelWriter(
        report, engine="openpyxl", mode="a",
        if_sheet_exists="replace"
    ) as pet_writer:
        pet_controls.to_excel(
            pet_writer, sheet_name="PET_control_comparison", index=False
        )

    # Separate complete InterPro workbook.
    # The raw interpro_web.tsv remains the lossless source file.
    # Excel has a 32,767-character hard limit per cell, so any exceptionally
    # long field is split into numbered continuation rows without truncation.
    interpro_full_xlsx = outdir / "InterPro_Full_Report.xlsx"

    def _excel_chunk_table(df, max_chars=30000):
        if df is None or df.empty:
            return pd.DataFrame()
        rows = []
        for _, src_row in df.iterrows():
            base = src_row.to_dict()
            long_cols = {
                c: str(v)
                for c, v in base.items()
                if pd.notna(v) and len(str(v)) > max_chars
            }
            if not long_cols:
                row = dict(base)
                row["_continuation_index"] = 1
                row["_continuation_total"] = 1
                rows.append(row)
                continue

            total = max(
                (len(v) + max_chars - 1) // max_chars
                for v in long_cols.values()
            )
            for idx in range(total):
                row = dict(base)
                for c, value in long_cols.items():
                    a = idx * max_chars
                    b = (idx + 1) * max_chars
                    row[c] = value[a:b]
                row["_continuation_index"] = idx + 1
                row["_continuation_total"] = total
                rows.append(row)
        return pd.DataFrame(rows)

    if ip is not None and not ip.empty:
        ip_excel = _excel_chunk_table(ip)
        with pd.ExcelWriter(interpro_full_xlsx, engine="openpyxl") as ip_writer:
            ip_excel.to_excel(
                ip_writer, sheet_name="All_InterPro_hits", index=False
            )

            if "pathways" in ip.columns:
                _p = ip["pathways"].fillna("").astype(str).str.strip()
                ip_path = ip[(_p != "") & (_p != "-")].copy()
                if not ip_path.empty:
                    _excel_chunk_table(ip_path).to_excel(
                        ip_writer, sheet_name="Pathway_metadata", index=False
                    )

            summary_cols = [
                "query", "panther_hits", "pfam_hits", "cdd_hits", "ncbifam_hits",
                "prints_hits", "gene3d_hits", "superfamily_hits", "smart_hits",
                "prosite_hits", "integrated_interpro_entries", "interpro_go_terms",
                "interpro_coordinates", "resolver_status", "resolver_note",
            ]
            summary_frames = [frame for frame in [hyp_best, unmatched_best]
                              if frame is not None and not frame.empty]
            summary_best = (pd.concat(summary_frames, ignore_index=True)
                            if summary_frames else pd.DataFrame())
            available = [c for c in summary_cols if c in summary_best.columns]
            if available:
                _excel_chunk_table(
                    summary_best[available].drop_duplicates()
                ).to_excel(
                    ip_writer, sheet_name="Protein_hit_summary", index=False
                )

        print(f"[info] Full InterPro report: {interpro_full_xlsx}")

    if interpro_full_xlsx.is_file():
        status_sheet(interpro_full_xlsx, interpro_state, pd)

    print("\nDONE")
    print("Report:", report.resolve())
    print("InterPro candidate FASTA:", interpro_faa.resolve())

    if not final.empty:
        print("Final candidate rows:", len(final))

        print(
            "High-confidence:",
            int(
                (
                    (final["Sequence match strength"] == "High")
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
