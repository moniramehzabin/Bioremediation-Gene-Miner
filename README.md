# Bioremediation Gene Miner

**Version 0.3.2 — research prototype**

A Python/DIAMOND workflow for mining bacterial whole-genome annotations for
**annotated and previously hypothetical proteins** potentially associated with:

- azo-dye degradation and oxidative dye transformation
- Cr(VI) reduction/resistance and other heavy metals
- aromatic compounds and hydrocarbons
- pesticides and dehalogenation
- nitro compounds
- xenobiotic oxidation / cytochrome P450
- microplastic / polymer biodegradation candidates
- oxidative-stress and UV/DNA-repair support systems

## Why this tool exists

A bacterial genome may contain thousands of CDSs. Some useful bioremediation
genes are already annotated, while others remain labelled `hypothetical protein`.

This prototype combines:

`existing annotation → curated rules`

and

`hypothetical proteins → DIAMOND homology → coverage filters → optional InterPro evidence`

into one evidence-ranked report.

## Important scientific limitation

A sequence match is **not experimental proof** of pollutant degradation.

The tool deliberately distinguishes:

- **Direct**
- **Candidate**
- **Supporting**
- **Review / rejected partial matches**

Broad families such as P450s, generic esterases/lipases, laccases and
peroxidases must not be interpreted as proof of microplastic or xenobiotic
degradation without stronger substrate-specific evidence.

## Requirements

- Windows, Linux or macOS
- Python 3.10+
- DIAMOND
- Python packages in `requirements.txt`

Install Python packages:

```bash
py -m pip install -r requirements.txt
```

On Linux/macOS, replace `py` with `python3` where appropriate.

## Build the reference database

First build the curated reference FASTA:

```bash
py build_reference_db_v0.2.1.py
```

Then create the DIAMOND database:

```bash
diamond makedb --in bioremediation_reference_v0.2.1.faa --db bioremediation_reference_v0.2.1
```

## Run the miner

Example on Windows:

```bash
py bioremediation_gene_miner.py genome.gbk ^
  --reference-db bioremediation_reference_v0.2.1 ^
  --reference-metadata bioremediation_reference_v0.2.1_metadata.tsv ^
  --annotation-rules annotation_rules.tsv ^
  --outdir results
```

The program automatically:

1. parses CDS features from the GenBank file
2. screens already annotated genes
3. extracts hypothetical proteins
4. runs DIAMOND
5. calculates query and reference coverage
6. rejects/downgrades short partial matches
7. creates an Excel report
8. writes `interpro_candidates.faa` for optional conserved-domain validation

## Optional InterPro validation

Upload `results/interpro_candidates.faa` to InterPro and download the TSV output.

Then rerun:

```bash
py bioremediation_gene_miner.py genome.gbk ^
  --reference-db bioremediation_reference_v0.2.1 ^
  --reference-metadata bioremediation_reference_v0.2.1_metadata.tsv ^
  --annotation-rules annotation_rules.tsv ^
  --interpro interpro_results.tsv ^
  --outdir results_with_interpro
```

## Main output workbook

`Bioremediation_Gene_Miner_Report.xlsx`

Sheets include:

- Summary
- All_candidates
- High_confidence
- Annotated_candidates
- Predicted_hypothetical
- Microplastic_candidates
- Supporting_functions
- Review_required
- All_DIAMOND_hits
- InterPro_evidence (when supplied)

## Development validation

The prototype has been technically tested on two bacterial WGS annotations.
These tests demonstrate that the workflow runs and that its filtering can
recover domain-supported hidden candidates while rejecting misleading short
partial hits. They **do not establish biological sensitivity/specificity**.

## Repository status

This is a research prototype. Reference curation, benchmarking, automated
domain integration and formal validation are ongoing.


## v0.3.2 patch

- Fixed CLI high-confidence count to match the Excel report (`High` + `Candidate`).
- Removed the openpyxl deprecation warning from header formatting.


## v0.3.2 patch

- Added **product-aware annotation matching**.
- Ambiguous gene symbols no longer override an explicit contradictory product annotation.
- Example fixed during testing: `xylE` annotated as `D-xylose-proton symporter` is no longer misclassified as catechol 2,3-dioxygenase.
- Selected ambiguous families now require supporting product text rather than gene-symbol-only matches.
