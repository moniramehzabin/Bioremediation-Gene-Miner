#!/usr/bin/env python3
from urllib.parse import urlencode
from urllib.request import Request, urlopen
import csv, json, time
from pathlib import Path

BASE = "https://rest.uniprot.org"
KB = Path("bioremediation_knowledge_base_v0.2.1.tsv")
OUT_FASTA = Path("bioremediation_reference_v0.2.1.faa")
CURATED_EXTERNAL = Path("curated_external_references_v0.2.1.tsv")
OUT_META = Path("bioremediation_reference_v0.2.1_metadata.tsv")
MAX_HITS_PER_RULE = 8

def get_json(url, params=None):
    if params:
        url += "?" + urlencode(params)
    req = Request(url, headers={"User-Agent":"Bioremediation-Gene-Miner/0.2.1"})
    with urlopen(req, timeout=90) as r:
        return json.loads(r.read().decode("utf-8"))

def seq(entry):
    return entry.get("sequence",{}).get("value","").replace("\n","").strip()

def protein_name(entry):
    p = entry.get("proteinDescription",{})
    r = p.get("recommendedName",{}).get("fullName",{})
    if isinstance(r,dict) and r.get("value"):
        return r["value"]
    s = p.get("submissionNames",[])
    if s:
        fn = s[0].get("fullName",{})
        if isinstance(fn,dict):
            return fn.get("value","")
    return ""

def genes(entry):
    vals=[]
    for g in entry.get("genes",[]):
        n=g.get("geneName",{})
        if isinstance(n,dict) and n.get("value"):
            vals.append(n["value"])
    return ";".join(vals)

with KB.open(encoding="utf-8") as f:
    rules = list(csv.DictReader(f, delimiter="\t"))

records = {}
provenance = {}
qc_review = []

APPROVED_GENE_ALIASES = {
    "chra": {"chra1"},
    "merb": {"merb1"},
    "nahac": {"ndob"},
    "nahad": {"ndoc"},
}

# Optional protein-name consistency QC for gene symbols that are reused
# for biologically different proteins.
FUNCTION_NAME_QC = {
    "linb": ("haloalkane dehalogenase",),
}

for i, rule in enumerate(rules,1):
    q0 = rule["uniprot_query"].strip()
    query = q0 if q0.startswith("accession:") else f"({q0}) AND (reviewed:true) AND (taxonomy_id:2)"
    print(f"[{i}/{len(rules)}] {rule['family_or_target']}")
    try:
        if q0.startswith("accession:"):
            acc = q0.split(":",1)[1]
            entries = [get_json(f"{BASE}/uniprotkb/{acc}")]
        else:
            data = get_json(f"{BASE}/uniprotkb/search", {
                "query": query, "format":"json", "size":str(MAX_HITS_PER_RULE)
            })
            entries = data.get("results", [])

        for e in entries:
            acc=e.get("primaryAccession","")
            s=seq(e)
            if not acc or not s:
                continue

            # QC test: for gene_exact queries, flag results whose primary
            # gene name does not match the requested gene.
            if q0.startswith("gene_exact:"):
                wanted_gene = q0.split(":", 1)[1].strip().lower()
                returned_genes = [
                    g.strip().lower()
                    for g in genes(e).split(";")
                    if g.strip()
                ]

                alias_ok = any(
                    g in APPROVED_GENE_ALIASES.get(wanted_gene, set())
                    for g in returned_genes
                )

                if wanted_gene not in returned_genes and not alias_ok:
                    qc_review.append({
                        "Accession": acc,
                        "Requested_gene": wanted_gene,
                        "Returned_gene": genes(e),
                        "Protein": protein_name(e),
                        "Family": rule["family_or_target"],
                        "Interpretation": rule["interpretation_rule"],
                        "Source_query": query,
                    })
                    print(
                        f"  [QC-REVIEW] {acc}: wanted={wanted_gene} "
                        f"returned={genes(e)} | {protein_name(e)}"
                    )
                    continue

                # Extra functional-name QC for gene symbols that can refer
                # to biologically different proteins.
                required_names = FUNCTION_NAME_QC.get(wanted_gene)
                if required_names:
                    returned_protein = protein_name(e).lower()
                    if not any(name in returned_protein for name in required_names):
                        qc_review.append({
                            "Accession": acc,
                            "Requested_gene": wanted_gene,
                            "Returned_gene": genes(e),
                            "Protein": protein_name(e),
                            "Family": rule["family_or_target"],
                            "Interpretation": rule["interpretation_rule"],
                            "Source_query": query,
                        })
                        print(
                            f"  [QC-REVIEW-FUNCTION] {acc}: gene={wanted_gene} "
                            f"| {protein_name(e)}"
                        )
                        continue
            if acc not in records:
                records[acc] = {
                    "Accession":acc,
                    "UniProt_ID":e.get("uniProtkbId",""),
                    "Gene":genes(e),
                    "Protein":protein_name(e),
                    "Organism":e.get("organism",{}).get("scientificName",""),
                    "Length_aa":len(s),
                    "Sequence":s
                }
                provenance[acc]=[]
            provenance[acc].append({
                "family_or_target":rule["family_or_target"],
                "major_category":rule["major_category"],
                "pollutant_or_role":rule["pollutant_or_role"],
                "evidence_class":rule["evidence_class"],
                "default_confidence":rule["default_confidence"],
                "interpretation_rule":rule["interpretation_rule"],
                "specificity_class":rule["specificity_class"],
                "headline_category_allowed":rule["headline_category_allowed"],
                "source_query":query
            })
    except Exception as ex:
        print("  ! failed:", ex)
    time.sleep(0.2)


# Add manually curated external references (for example, verified NCBI proteins).
if CURATED_EXTERNAL.exists():
    with CURATED_EXTERNAL.open(encoding="utf-8") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            acc = row["Accession"].strip()
            fasta_source = row["FASTA_file"].strip()

            if not acc or not fasta_source:
                continue

            if fasta_source.startswith(("http://", "https://")):
                req = Request(
                    fasta_source,
                    headers={"User-Agent": "Bioremediation-Gene-Miner/0.2.1"},
                )
                with urlopen(req, timeout=90) as response:
                    fasta_text = response.read().decode("utf-8")
            else:
                fasta_path = Path(fasta_source)
                if not fasta_path.exists():
                    continue
                fasta_text = fasta_path.read_text()

            sequence = "".join(
                line.strip()
                for line in fasta_text.splitlines()
                if not line.startswith(">")
            )

            if not sequence:
                continue

            records[acc] = {
                "Accession": acc,
                "UniProt_ID": "",
                "Gene": row["Gene"],
                "Protein": row["Protein"],
                "Organism": row["Organism"],
                "Length_aa": len(sequence),
                "Sequence": sequence,
            }

            provenance[acc] = [{
                "family_or_target": row["family_or_target"],
                "major_category": row["major_category"],
                "pollutant_or_role": row["pollutant_or_role"],
                "evidence_class": row["evidence_class"],
                "default_confidence": row["default_confidence"],
                "interpretation_rule": row["interpretation_rule"],
                "specificity_class": row["specificity_class"],
                "headline_category_allowed": row["headline_category_allowed"],
                "source_query": row["source_query"],
            }]

            print(f"[CURATED] added {acc} ({row['Gene']})")

def choose_headline(ps):
    # only categories from rules allowed to define a headline are considered.
    eligible = [p for p in ps if p["headline_category_allowed"]=="YES"]
    if not eligible:
        return "Broad/supporting family", "REVIEW"
    cats = sorted({p["major_category"] for p in eligible})
    if len(cats) == 1:
        return cats[0], "OK"
    # If multiple allowed categories point to the same accession, do not merge them
    # into a misleading combined label; flag for review.
    return "Ambiguous multi-category", "REVIEW"

with OUT_FASTA.open("w", encoding="utf-8") as f:
    for acc,r in records.items():
        headline, flag = choose_headline(provenance[acc])
        hdr = f">{acc}|{r['Gene'] or 'NA'}|{headline.replace(' ','_')}|{r['Protein']}"
        f.write(hdr.replace("\n"," ")+"\n")
        s=r["Sequence"]
        for j in range(0,len(s),80):
            f.write(s[j:j+80]+"\n")

fields = [
    "Accession","UniProt_ID","Gene","Protein","Organism","Length_aa",
    "Headline_category","Category_flag","All_family_targets","All_major_categories",
    "Pollutants_or_roles","Evidence_classes","Specificity_classes",
    "Default_confidences","Interpretation_rules","Source_queries"
]
with OUT_META.open("w", newline="", encoding="utf-8") as f:
    w=csv.DictWriter(f, fieldnames=fields, delimiter="\t")
    w.writeheader()
    for acc,r in records.items():
        ps=provenance[acc]
        headline, flag = choose_headline(ps)
        w.writerow({
            "Accession":acc,
            "UniProt_ID":r["UniProt_ID"],
            "Gene":r["Gene"],
            "Protein":r["Protein"],
            "Organism":r["Organism"],
            "Length_aa":r["Length_aa"],
            "Headline_category":headline,
            "Category_flag":flag,
            "All_family_targets":"; ".join(sorted({p["family_or_target"] for p in ps})),
            "All_major_categories":"; ".join(sorted({p["major_category"] for p in ps})),
            "Pollutants_or_roles":"; ".join(sorted({p["pollutant_or_role"] for p in ps})),
            "Evidence_classes":"; ".join(sorted({p["evidence_class"] for p in ps})),
            "Specificity_classes":"; ".join(sorted({p["specificity_class"] for p in ps})),
            "Default_confidences":"; ".join(sorted({p["default_confidence"] for p in ps})),
            "Interpretation_rules":"; ".join(sorted({p["interpretation_rule"] for p in ps})),
            "Source_queries":"; ".join(sorted({p["source_query"] for p in ps}))
        })

with Path("bioremediation_reference_v0.2.1_QC_review.tsv").open(
    "w", newline="", encoding="utf-8"
) as f:
    fields = [
        "Accession",
        "Requested_gene",
        "Returned_gene",
        "Protein",
        "Family",
        "Interpretation",
        "Source_query",
    ]
    w = csv.DictWriter(f, fieldnames=fields, delimiter="\t")
    w.writeheader()
    w.writerows(qc_review)

print("\nDONE")
print("Unique reference proteins:", len(records))
print("Created:", OUT_FASTA.resolve())
print("Created:", OUT_META.resolve())
print("Safeguard: broad/supporting families cannot define microplastic/direct headline categories.")
print("Safeguard: conflicting direct categories are flagged Ambiguous multi-category / REVIEW.")
