"""Look up assay evidence for a Gene Miner hit without a new screening run."""

import argparse
import csv
import hashlib
import json
import re
from concurrent.futures import ThreadPoolExecutor, as_completed
from time import sleep
from datetime import datetime, timezone
from html import escape
from pathlib import Path
from urllib.error import HTTPError, URLError
from urllib.parse import quote, urlencode
from urllib.request import urlopen

EUROPE_PMC = "https://www.ebi.ac.uk/europepmc/webservices/rest/search"
CONTROLS = Path(__file__).resolve().parent / "pet_activity_controls" / "controls.tsv"


def exact_controls(sequence, controls=CONTROLS):
    if not sequence:
        return []
    digest = hashlib.sha256(sequence.encode("ascii")).hexdigest()
    with controls.open(encoding="utf-8-sig", newline="") as handle:
        return [
            record for record in csv.DictReader(handle, delimiter="\t")
            if record["SequenceSHA256"] == digest
        ]


def load_online_json(url, timeout):
    for attempt in range(2):
        try:
            with urlopen(url, timeout=timeout) as response:
                return json.load(response)
        except HTTPError as exc:
            if attempt or exc.code not in {502, 503, 504}:
                raise
        except (TimeoutError, URLError):
            if attempt:
                raise
        sleep(1)


def reference_search_terms(accession, timeout):
    """Build literature terms from a UniProt record when the accession resolves."""
    try:
        entry = load_online_json(
            f"https://rest.uniprot.org/uniprotkb/{quote(accession)}.json", timeout
        )
    except (OSError, ValueError):
        return [], "unavailable"

    def clean(value):
        return re.sub(r'[^\w .-]', ' ', value).strip()

    terms = []
    for item in entry.get("references", []):
        for identifier in item.get("citation", {}).get("citationCrossReferences", []):
            if identifier.get("database") == "PubMed" and identifier["id"].isdigit():
                terms.append(f"SRC:MED AND EXT_ID:{identifier['id']}")

    for gene in entry.get("genes", []):
        for identifier in [gene.get("geneName", {}), *gene.get("synonyms", [])]:
            name = clean(identifier.get("value", ""))
            if name:
                terms.append(f'TITLE_ABS:"{name}"')

    organism = clean(entry.get("organism", {}).get("scientificName", ""))
    protein = clean(entry.get("proteinDescription", {}).get("recommendedName", {})
                    .get("fullName", {}).get("value", ""))
    if organism and protein:
        species = " ".join(organism.split()[:2])
        words = [word for word in protein.split() if word.lower() not in {"protein", "putative"}]
        if words:
            keywords = " AND ".join(f'TITLE_ABS:"{word}"' for word in words)
            terms.append(f'TITLE_ABS:"{species}" AND {keywords}')
    return list(dict.fromkeys(terms)), "resolved"


def literature_search(accession, timeout, family=None, product=None):
    if accession is not None and not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]{0,79}", accession):
        raise ValueError(f"Invalid protein accession: {accession}")
    if not family or re.search(r"\bPET(?:ase)?\b|polyethylene terephthalate", family, re.I):
        activity_terms = "PET OR polyethylene terephthalate"
    else:
        names = dict.fromkeys(
            re.sub(r'[^\w .,+()-]', ' ', name).strip()
            for name in (family, product or "")
        )
        activity_terms = " OR ".join(f'"{name}"' for name in names if name)
    linked = []
    if accession is None:
        terms = " OR ".join(
            f'TITLE_ABS:"{name}"' for name in dict.fromkeys(
                re.sub(r'[^\w .,+()-]', ' ', value).strip()
                for value in (family or "", product or "")
            )
            if name and name.casefold() not in {
                "protein", "hypothetical protein", "uncharacterized protein",
                "unknown", "uncharacterized", "hypothetical",
            }
        )
        if not terms:
            return {"status": "no_searchable_description"}
        metadata_status = "not_applicable"
        search_basis = "family_keywords"
    else:
        terms = f'"{accession}"' + (f" AND ({activity_terms})" if activity_terms else "")
        reference_terms, metadata_status = reference_search_terms(accession, timeout)
        linked = [term for term in reference_terms if term.startswith("SRC:MED AND EXT_ID:")]
        if linked:
            terms = " OR ".join(f"({term})" for term in linked)
        elif reference_terms:
            terms = " OR ".join(f"({term})" for term in [terms, *reference_terms])
        search_basis = "uniprot_linked_pubmed" if linked else "accession_and_keywords"
    parameters = urlencode({
        "query": terms, "format": "json", "resultType": "lite", "pageSize": 10,
    })
    try:
        data = load_online_json(f"{EUROPE_PMC}?{parameters}", timeout)
        articles = data["resultList"]["result"]
        if not isinstance(articles, list) or any(
            not isinstance(article, dict) for article in articles
        ):
            raise ValueError("Unexpected Europe PMC search response.")
    except (URLError, OSError, ValueError, KeyError, TypeError) as exc:
        return {"accession": accession, "status": "search_failed", "error": str(exc)}
    return {
        "accession": accession,
        "search_query": terms,
        "uniprot_metadata": metadata_status,
        "search_basis": search_basis,
        "status": "unverified_publications_found" if articles else "no_publications_found",
        "publications": [
            {
                "title": article.get("title"),
                "doi": article.get("doi"),
                "pmid": article.get("pmid"),
                "pmcid": article.get("pmcid"),
            }
            for article in articles
        ],
    }


def report_search_rows(report):
    with report.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"Locus_tag", "Family_target", "Decision", "Best_reference"}
        if not required.issubset(reader.fieldnames or ()):
            raise ValueError("Expected a Gene Miner family-level candidates TSV.")
        return [row for row in reader if row["Decision"] in {"Candidate", "Review"}]


def search_report_rows(rows, proteins, online, timeout=10, progress=None):
    """Use original CDS accessions when available and label homolog leads separately."""
    results = []
    targets = {}
    for row in rows:
        locus, family, reference = (row[key] for key in ("Locus_tag", "Family_target", "Best_reference"))
        accession, sequence = proteins.get(locus, (None, None))
        pet_family = bool(re.search(r"\bPET(?:ase)?\b|polyethylene terephthalate", family, re.I))
        controls = exact_controls(sequence) if pet_family and sequence else []
        entry = {
            "locus_tag": locus, "family": family, "decision": row["Decision"],
            "reference": reference, "query_accession": accession,
            "local_pet_lookup": ("not_applicable" if not pet_family else
                                 "completed" if sequence else "sequence_not_available"),
            "exact_query_sequence_controls": controls,
            "query_publication_leads": None, "homolog_publication_leads": None,
            "family_publication_leads": None,
        }
        if online:
            reference_accession = (reference or "").partition("|")[0].strip()
            for field, identifier in (
                ("query_publication_leads", accession),
                ("homolog_publication_leads", reference_accession if reference_accession != accession else None),
            ):
                if (not identifier or identifier.startswith(("PETCTRL", "CURATED_")) or
                        not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]{0,79}", identifier)):
                    entry[field] = {"status": "no_searchable_accession"}
                    continue
                key = (identifier, family, row.get("Product_or_prediction"))
                targets.setdefault(key, []).append((entry, field))
            if all(entry[field] is not None for field in (
                "query_publication_leads", "homolog_publication_leads",
            )):
                description = row.get("Product_or_prediction") or row.get("Gene")
                key = (None, family, description)
                targets.setdefault(key, []).append((entry, "family_publication_leads"))
        results.append(entry)

    if targets:
        run_lookups(targets, timeout, progress)
    return {
        "searched_at_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "literature_source": "Europe PMC" if online else None,
        "selected_rows": results,
    }


def retry_failed_lookups(report, timeout=20, progress=None):
    """Retry failed requests in a loaded report, retaining all other results."""
    targets = {}
    for row in report["selected_rows"]:
        for field in ("query_publication_leads", "homolog_publication_leads",
                      "family_publication_leads"):
            lead = row.get(field)
            if lead and lead.get("status") == "search_failed":
                key = (lead.get("accession"), row["family"], lead.get("search_description"))
                targets.setdefault(key, []).append((row, field))
    if targets:
        run_lookups(targets, timeout, progress)
        report["retried_at_utc"] = datetime.now(timezone.utc).isoformat(timespec="seconds")
    return sum(map(len, targets.values()))


def run_lookups(targets, timeout, progress):
    total = len(targets)
    if progress:
        progress(0, total)

    with ThreadPoolExecutor(max_workers=4) as pool:
        futures = {
            pool.submit(literature_search, accession, timeout, family, product):
                (accession, family, product)
            for accession, family, product in targets
        }
        for done_count, future in enumerate(as_completed(futures), start=1):
            key = futures[future]
            try:
                result = future.result()
            except (OSError, ValueError, TypeError, KeyError, AttributeError) as exc:
                result = {"accession": key[0], "status": "search_failed", "error": str(exc)}
            if result["status"] == "search_failed" and key[2]:
                result["search_description"] = key[2]
            for entry, field in targets[key]:
                entry[field] = result
            if progress:
                progress(done_count, total)


def render_report_html(report):
    """Render a standalone summary of a launcher's assay evidence search."""
    sections = []
    for row in report["selected_rows"]:
        parts = [
            f"<h2>{escape(row['locus_tag'])} · {escape(row['family'])}</h2>",
            f"<p><strong>{escape(row['decision'])}</strong> · Reference {escape(row['reference'] or 'unavailable')}</p>",
        ]
        controls = row["exact_query_sequence_controls"]
        if row["local_pet_lookup"] == "completed":
            parts.append(f"<p>Exact PET assay records: {len(controls)}</p>")
        elif row["local_pet_lookup"] == "sequence_not_available":
            parts.append("<p>PET control lookup requires a saved protein sequence.</p>")
        for label, field in (
            ("Selected protein", "query_publication_leads"),
            ("Reference protein", "homolog_publication_leads"),
            ("Protein-family literature", "family_publication_leads"),
        ):
            lead = row.get(field)
            if lead is None:
                continue
            parts.append(f"<h3>{label}</h3>")
            if lead.get("accession"):
                parts.append(f"<p>Accession: {escape(lead['accession'])}</p>")
            if lead["status"] == "unverified_publications_found":
                parts.append("<ul>")
                for paper in lead["publications"]:
                    title = escape(re.sub(r"</?(?:i|b|em|strong|sub|sup)>", "", paper.get("title") or "Untitled paper", flags=re.I))
                    doi, pmid = paper.get("doi"), paper.get("pmid")
                    if doi or pmid:
                        url = ("https://doi.org/" + quote(doi)) if doi else ("https://europepmc.org/article/MED/" + quote(pmid))
                        title = f'<a href="{escape(url, quote=True)}">{title}</a>'
                    parts.append(f"<li>{title}</li>")
                parts.append("</ul>")
            else:
                status = {
                    "no_searchable_description": "No usable protein family or description was available.",
                    "no_searchable_accession": "No searchable accession was available.",
                    "no_publications_found": "The search returned no papers.",
                    "search_budget_exceeded": "Lookup deferred when the search time limit was reached.",
                    "search_failed": "The literature search could not complete.",
                }[lead["status"]]
                parts.append(f"<p>{status}</p>")
        sections.append("<article>" + "".join(parts) + "</article>")

    source = escape(report["literature_source"] or "Literature search off")
    return f"""<!doctype html>
<html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width">
<title>Gene Miner assay evidence</title>
<style>
body {{ font: 1rem/1.55 system-ui, sans-serif; color: #192522; background: #f1f5f3;
        max-width: 58rem; margin: 2rem auto; padding: 0 1rem; }}
article {{ background: white; padding: 1rem 1.5rem; border: 1px solid #ccd9d3;
           border-radius: .5rem; margin: 1rem 0; }}
h2 {{ margin-top: .25rem; }} h3 {{ margin-bottom: .25rem; }}
a {{ color: #145c50; }} li {{ margin-bottom: .5rem; }}
</style>
<h1>Gene Miner assay evidence</h1>
<p>{len(sections)} selected rows · {source} · {escape(report["searched_at_utc"])}</p>
<p>Publication matches are leads that need verification. The JSON export contains the full search details and assay records.</p>
{"".join(sections)}
</html>
"""


def selected_report_evidence(report, loci, fasta):
    """Return evidence for explicitly selected Candidate or Review loci."""
    with report.open(encoding="utf-8-sig", newline="") as handle:
        rows = csv.DictReader(handle, delimiter="\t")
        required = {"Locus_tag", "Family_target", "Decision", "Best_reference"}
        if not required.issubset(rows.fieldnames or ()):
            raise ValueError("Expected a Gene Miner family-level candidates TSV.")
        selected = {locus: [] for locus in loci}
        for row in rows:
            locus = row["Locus_tag"]
            if locus in selected and row["Decision"] in {"Candidate", "Review"}:
                selected[locus].append({field: row[field] for field in (
                    "Family_target", "Decision", "Best_reference",
                )})

    sequences = {}
    for record in fasta.read_text(encoding="utf-8-sig").split(">")[1:]:
        header, *lines = record.splitlines()
        fields = header.split()
        if not fields or fields[0] not in selected:
            continue
        locus = fields[0]
        sequence = "".join("".join(lines).split()).upper()
        if locus in sequences or not re.fullmatch(r"[A-Z]+", sequence):
            raise ValueError(f"Duplicate or invalid protein FASTA record for {locus}.")
        sequences[locus] = sequence

    results = []
    for locus, hits in selected.items():
        if not hits or locus not in sequences:
            raise ValueError(f"{locus} requires a Candidate or Review row and a matching FASTA record.")
        controls = exact_controls(sequences[locus])
        results.append({
            "locus_tag": locus,
            "screening_hits": hits,
            "query_assay_status": "exact_control_found" if controls else "no_verified_assay_evidence",
            "exact_sequence_controls": controls,
        })
    return results


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--accession", help="Query protein accession")
    parser.add_argument("--fasta", type=Path, help="Query protein FASTA for local lookup")
    parser.add_argument("--report", type=Path, help="Gene Miner family-level candidates TSV")
    parser.add_argument("--locus", action="append", help="Selected locus tag, repeatable with --report")
    parser.add_argument("--reference", help="Accession of a homolog from the Gene Miner report")
    parser.add_argument("--online", action="store_true", help="Send supplied accessions to Europe PMC")
    parser.add_argument("--timeout", type=float, default=10)
    args = parser.parse_args()
    if args.report or args.locus:
        if not (args.report and args.locus and args.fasta):
            parser.error("Report selection requires --report, --locus, and --fasta.")
        if args.accession or args.reference or args.online:
            parser.error("Report mode uses the local PET controls. Use single-protein mode for online lookup.")
        try:
            selected = selected_report_evidence(args.report, args.locus, args.fasta)
        except (OSError, ValueError) as exc:
            parser.error(str(exc))
        print(json.dumps({
            "searched_at_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
            "selected_proteins": selected,
        }, indent=2))
        return
    if not (args.accession or args.fasta):
        parser.error("Provide --accession or --fasta.")
    if args.online and not (args.accession or args.reference):
        parser.error("Online search requires at least one accession.")
    if args.timeout <= 0:
        parser.error("--timeout must be positive.")
    for accession in (args.accession, args.reference):
        if accession and not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]{0,79}", accession):
            parser.error(f"Invalid accession: {accession}")
    sequence = None
    if args.fasta:
        try:
            header, *lines = args.fasta.read_text(encoding="utf-8-sig").splitlines()
            sequence = "".join("".join(lines).split()).upper()
        except (OSError, ValueError) as exc:
            parser.error(str(exc))
        if not header.startswith(">") or not re.fullmatch(r"[A-Z]+", sequence):
            parser.error("Expected a single protein FASTA record.")

    controls = exact_controls(sequence)
    result = {
        "searched_at_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "query_accession": args.accession,
        "query_assay_status": (
            "exact_control_found" if controls else
            "no_verified_assay_evidence" if sequence else
            "local_sequence_check_not_performed"
        ),
        "local_control_lookup": "completed" if sequence else "sequence_not_supplied",
        "exact_sequence_controls": controls,
        "literature_source": "Europe PMC" if args.online else None,
        "query_publication_leads": None,
        "homolog_publication_leads": None,
    }
    if args.online and args.accession:
        result["query_publication_leads"] = literature_search(args.accession, args.timeout)
    if args.online and args.reference and args.reference != args.accession:
        result["homolog_publication_leads"] = literature_search(args.reference, args.timeout)
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
