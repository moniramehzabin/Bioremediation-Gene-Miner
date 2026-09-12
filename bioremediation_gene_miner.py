#!/usr/bin/env python3
"""
Bioremediation Gene Miner v0.3.7-release-candidate-frozen

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

VERSION = "0.3.7-release-candidate"


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
        "User-Agent": "Bioremediation-Gene-Miner/0.3.7-release-candidate"
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
    fasta_path, outdir = Path(fasta_path), Path(outdir)
    status_file = outdir / "interpro_web_status.txt"
    tsv_path = outdir / "interpro_web.tsv"

    if not fasta_path.exists() or fasta_path.stat().st_size == 0:
        status_file.write_text("NOT_RUN\tNo Weak+Review candidates.\n", encoding="utf-8")
        print("[info] InterPro web: no candidates; skipped.")
        return None
    if not email:
        status_file.write_text("NOT_RUN\t--interpro-email missing.\n", encoding="utf-8")
        print("[warn] --interpro-auto requested but --interpro-email is missing.")
        return None

    payload = urllib.parse.urlencode({
        "email": email,
        "title": "Bioremediation Gene Miner weak-review candidates",
        "stype": "p",
        "sequence": fasta_path.read_text(encoding="utf-8"),
        "goterms": "true",
        "pathways": "true",
    }).encode("utf-8")

    try:
        print("[run] submitting Weak+Review candidates to EMBL-EBI InterProScan...")
        job_id = _http_text(f"{INTERPRO_REST_BASE}/run/", data=payload, timeout=120).strip()
        if not job_id:
            raise RuntimeError("Empty InterProScan job ID.")
        (outdir / "interpro_job_id.txt").write_text(job_id + "\n", encoding="utf-8")
        print("[info] InterPro job:", job_id)

        deadline = time.time() + timeout_minutes * 60
        last = None
        while time.time() < deadline:
            status = _http_text(f"{INTERPRO_REST_BASE}/status/{job_id}").strip()
            if status != last:
                print("[info] InterPro status:", status)
                last = status
            if status == "FINISHED":
                break
            if status not in {"QUEUED", "RUNNING", "PENDING"}:
                raise RuntimeError(f"InterProScan ended with status: {status}")
            time.sleep(max(3, int(poll_seconds)))
        else:
            raise TimeoutError(f"InterProScan exceeded {timeout_minutes} minute timeout.")

        ids = {x.get("identifier", "") for x in _interpro_result_types(job_id)}
        if "tsv" not in ids:
            raise RuntimeError("TSV unavailable; result types: " + ", ".join(sorted(ids)))

        req = urllib.request.Request(
            f"{INTERPRO_REST_BASE}/result/{job_id}/tsv",
            headers={"User-Agent": "Bioremediation-Gene-Miner/0.3.7-release-candidate"},
        )
        with urllib.request.urlopen(req, timeout=120) as response:
            tsv_path.write_bytes(response.read())

        status_file.write_text(f"FINISHED\t{job_id}\t{tsv_path.name}\n", encoding="utf-8")
        print("[info] InterPro TSV:", tsv_path.resolve())
        return str(tsv_path)
    except Exception as exc:
        status_file.write_text(f"FAILED\t{type(exc).__name__}: {exc}\n", encoding="utf-8")
        print(f"[warn] InterPro web failed: {exc}")
        print("[warn] continuing without InterPro evidence; FASTA preserved.")
        return None

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
            final["InterPro_status"].eq("HITS_FOUND")
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


def export_interpro_candidates(cds, hyp_best, out_faa):
    if hyp_best is None or hyp_best.empty:
        Path(out_faa).write_text("", encoding="utf-8")
        return

    wanted = set(
        hyp_best.loc[
            (hyp_best["decision"] == "Review")
            & (hyp_best["confidence"] == "Weak"),
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
        "--interpro-auto",
        action="store_true",
        help="Automatically submit Weak+Review candidates to EMBL-EBI InterProScan REST.",
    )
    ap.add_argument("--interpro-email", help="Email required by EMBL-EBI Job Dispatcher.")
    ap.add_argument("--interpro-poll-seconds", type=int, default=5)
    ap.add_argument("--interpro-timeout-minutes", type=int, default=45)

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

    interpro_faa = outdir / "interpro_candidates.faa"
    export_interpro_candidates(cds, hyp_best, interpro_faa)

    interpro_input = args.interpro
    if args.interpro_auto and args.interpro:
        print("[warn] --interpro supplied; skipping automatic web submission.")
    elif args.interpro_auto:
        interpro_input = run_interpro_web(
            interpro_faa, args.interpro_email, outdir,
            args.interpro_poll_seconds, args.interpro_timeout_minutes,
        )

    ip = load_interpro(interpro_input)
    hyp_best = add_interpro_support(hyp_best, ip)

    # Distinguish proteins not submitted to InterPro from submitted proteins
    # that returned no usable hits. This avoids misleading FALSE values.
    submitted_interpro = set()
    if interpro_faa.exists() and interpro_faa.stat().st_size > 0:
        with open(interpro_faa, encoding="utf-8") as _f:
            for _line in _f:
                if _line.startswith(">"):
                    submitted_interpro.add(_line[1:].strip().split()[0])

    if hyp_best is not None and not hyp_best.empty:
        hyp_best.loc[
            hyp_best["query"].astype(str).isin(submitted_interpro),
            "interpro_status",
        ] = "NO_HITS"
        hyp_best.loc[
            hyp_best["interpro_detected"].eq(True),
            "interpro_status",
        ] = "HITS_FOUND"

    hyp_best = resolve_interpro_evidence(hyp_best)

    # Resolver has exactly four categories only for submitted/evaluated proteins.
    # Non-submitted candidates keep resolver fields blank.
    if hyp_best is not None and not hyp_best.empty:
        not_submitted = hyp_best["interpro_status"].eq("NOT_SUBMITTED")
        hyp_best.loc[not_submitted, "resolver_status"] = ""
        hyp_best.loc[not_submitted, "resolver_note"] = ""

        no_hits = hyp_best["interpro_status"].eq("NO_HITS")
        hyp_best.loc[no_hits, "resolver_status"] = "UNINFORMATIVE"
        hyp_best.loc[no_hits, "resolver_note"] = (
            "Protein was submitted to InterPro, but no usable InterPro evidence was returned."
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
            available = [c for c in summary_cols if c in hyp_best.columns]
            if available:
                _excel_chunk_table(
                    hyp_best[available].drop_duplicates()
                ).to_excel(
                    ip_writer, sheet_name="Protein_hit_summary", index=False
                )

        print(f"[info] Full InterPro report: {interpro_full_xlsx}")

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
