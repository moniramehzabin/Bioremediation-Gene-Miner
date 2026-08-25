
# Bioremediation Gene Miner

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.22046925.svg)](https://doi.org/10.5281/zenodo.22046925)

**Version 0.3.7 — Research Prototype**

Bioremediation Gene Miner is a Python/DIAMOND workflow for mining bacterial whole-genome annotations for **annotated and hypothetical proteins** potentially associated with environmental bioremediation.

The workflow combines existing genome annotation, curated protein-reference similarity searches, evidence-aware candidate classification, and optional InterPro-based domain/family evidence.

> **Important:** Bioremediation Gene Miner identifies computational candidates. A predicted gene or protein should not be interpreted as experimental confirmation of a bioremediation phenotype or biochemical activity.

---

## What the tool screens

The reference library and annotation rules cover bioremediation-associated functions including:

- aromatic compound degradation
- azo-dye degradation and oxidative dye transformation
- Cr(VI) reduction/resistance and other heavy-metal-associated functions
- pesticide/xenobiotic transformation
- dehalogenation
- nitro-compound transformation
- oxidative enzymes
- hydrocarbon-associated functions
- polymer/microplastic-associated candidate functions

The library can be expanded as additional experimentally supported reference proteins and bioremediation families are curated.

---

## Core workflow

```text
Bacterial GenBank genome
        |
        +--> Annotated CDS screening
        |       |
        |       +--> annotation_rules.tsv
        |
        +--> Hypothetical/uncharacterized proteins
                |
                +--> DIAMOND search
                        |
                        +--> curated bioremediation reference database
                        |
                        +--> family-level candidate resolution
                        |
                        +--> evidence/confidence classification
                        |
                        +--> selected Weak + Review candidates
                                |
                                +--> InterProScan
                                        |
                                        +--> PANTHER
                                        +--> Pfam
                                        +--> CDD
                                        +--> NCBIfam
                                        +--> other InterPro member databases
                                        |
                                        +--> evidence resolver
```

The workflow intentionally separates **annotation-derived evidence**, **sequence-homology evidence**, and **protein-family/domain evidence** rather than treating all predictions as equivalent.

---

## Evidence classes

Candidates are evaluated using multiple evidence sources.

### Existing annotation

Already annotated CDS features are screened using product-aware rules from `annotation_rules.tsv`.

### DIAMOND homology

Hypothetical or uncharacterized proteins are compared with the curated reference protein database.

The output retains useful alignment information including:

- sequence identity
- query coverage
- reference coverage
- E-value
- bit score
- best reference

Candidate assignment is performed at the curated biological-family level rather than treating every individual reference alignment as an independent biological prediction.

### InterPro evidence

Selected hypothetical candidates requiring additional evaluation can be submitted to EMBL-EBI InterProScan.

Version 0.3.7 can retain compact evidence from:

- PANTHER
- Pfam
- CDD
- NCBIfam
- integrated InterPro entries
- GO terms
- additional InterPro member databases

Verbose InterPro output is preserved separately so that the main candidate report remains interpretable.

---

## InterPro resolver

InterPro evidence is compared with the DIAMOND-derived family hypothesis.

The resolver can classify evaluated candidates as:

- **SUPPORTING** — independent domain/family evidence is consistent with the DIAMOND candidate
- **BROAD/RELATED** — InterPro supports a broader or related protein family but does not specifically establish the proposed function
- **CONFLICTING** — InterPro evidence favors an alternative interpretation
- **NO_HITS** — the submitted protein produced no usable InterPro evidence

Proteins that were not selected for InterPro evaluation remain explicitly identifiable as **NOT_SUBMITTED** through `InterPro_status`.

The resolver is an evidence-integration layer, not experimental validation.

---

## Confidence and review logic

The program distinguishes stronger candidates from uncertain or conflicting predictions.

Output categories include:

- high-confidence candidates
- annotated candidates
- predicted hypothetical candidates
- family-level candidates
- InterPro-supported candidates
- supporting functions
- candidates requiring review
- polymer/microplastic-associated candidates

Weak or conflicting evidence is retained for inspection rather than silently converted into a positive functional assignment.

---

## Input

Primary input:

```text
annotated bacterial genome in GenBank format (.gbk/.gbff)
```

The genome should contain CDS features and translated protein sequences.

Additional resources used by the workflow include:

```text
annotation_rules.tsv
bioremediation_reference_v0.2.1.dmnd
bioremediation_reference_v0.2.1_metadata.tsv
```

The repository also contains the corresponding curated reference FASTA and database-building resources.

---

## Requirements

- Python 3
- DIAMOND
- pandas
- Biopython
- openpyxl
- requests

Install Python dependencies with:

```bash
pip install -r requirements.txt
```

DIAMOND must be installed separately and available either through the system PATH or supplied using the command-line option.

---

## Basic usage

Example:

```bash
python bioremediation_gene_miner.py genome.gbk \
  --reference-db bioremediation_reference_v0.2.1 \
  --reference-metadata bioremediation_reference_v0.2.1_metadata.tsv \
  --annotation-rules annotation_rules.tsv \
  --outdir results
```

On Windows, an explicit DIAMOND executable can be supplied:

```cmd
py bioremediation_gene_miner.py genome.gbk --reference-db bioremediation_reference_v0.2.1 --reference-metadata bioremediation_reference_v0.2.1_metadata.tsv --annotation-rules annotation_rules.tsv --diamond "C:\path\to\diamond.exe" --outdir results
```

---

## Automated InterProScan

Version 0.3.7 supports automated submission of selected candidates to the EMBL-EBI InterProScan web service.

Example:

```cmd
py bioremediation_gene_miner.py genome.gbk --reference-db bioremediation_reference_v0.2.1 --reference-metadata bioremediation_reference_v0.2.1_metadata.tsv --annotation-rules annotation_rules.tsv --diamond "C:\path\to\diamond.exe" --interpro-auto --interpro-email your_email@example.com --outdir results
```

With automated InterPro enabled, the workflow:

1. identifies candidates selected for additional domain/family evaluation;
2. creates the InterPro candidate FASTA;
3. submits the sequences to EMBL-EBI InterProScan;
4. monitors job status;
5. retrieves the InterPro TSV result;
6. extracts compact family/domain evidence;
7. applies the InterPro resolver;
8. writes both compact and full InterPro reports.

Availability of automated InterPro analysis depends on the external EMBL-EBI service.

---

## Main output report

The primary workbook is:

```text
Bioremediation_Gene_Miner_Report.xlsx
```

It contains:

```text
Bioremediation_Gene_Miner_Report.xlsx
│
├── Summary
├── All_candidates
├── High_confidence
├── Annotated_candidates
├── Predicted_hypothetical
├── Family_level_candidates
├── InterPro_supported
├── Supporting_functions
├── Review_required
├── Microplastic_candidates
└── All_DIAMOND_hits
```

### `Summary`

Provides a compact overview of the genome-screening run.

### `All_candidates`

Main candidate table combining annotation-derived and hypothetical-protein candidates.

### `High_confidence`

Candidates meeting the tool's high-confidence criteria.

### `Annotated_candidates`

Candidates recovered directly from existing genome annotation.

### `Predicted_hypothetical`

Hypothetical/uncharacterized proteins recovered through the sequence-evidence workflow.

### `Family_level_candidates`

Family-resolved hypothetical candidate assignments.

### `InterPro_supported`

Candidates for which InterPro evidence provides supporting or otherwise informative family/domain evidence according to the resolver logic.

### `Supporting_functions`

Additional functions relevant to interpretation of the bioremediation potential.

### `Review_required`

Weak, broad, conflicting, or otherwise uncertain candidates retained for manual evaluation.

### `Microplastic_candidates`

Dedicated output for polymer/microplastic-associated candidate functions. An empty sheet means that no qualifying candidates were detected in that genome; it does not disable this screening category.

### `All_DIAMOND_hits`

Underlying DIAMOND alignment evidence retained for transparency and manual inspection.

---

## InterPro full report

When InterPro results are available, the workflow also produces:

```text
InterPro_Full_Report.xlsx
│
├── All_InterPro_hits
├── Pathway_metadata
└── Protein_hit_summary
```

This workbook preserves verbose InterPro-derived information separately from the main candidate report.

The design keeps the primary workbook compact while retaining the detailed evidence needed for auditing or deeper interpretation.

Additional files can include:

```text
interpro_web.tsv
interpro_candidates.faa
```

---

## Interpreting InterPro fields

The main candidate-facing tables retain compact fields such as:

```text
InterPro_status
PANTHER_hits
Pfam_hits
CDD_hits
NCBIfam_hits
Integrated_InterPro_entries
InterPro_GO_terms
Resolver_status
Resolver_note
```

Large or highly verbose InterPro fields are kept in `InterPro_Full_Report.xlsx` rather than duplicated throughout the main report.

---

## Biological interpretation

The program is designed for **candidate discovery and prioritization**.

A result such as a DIAMOND family match, Pfam domain, PANTHER family, CDD hit, NCBIfam assignment, or integrated InterPro entry provides computational evidence. It does not by itself demonstrate that the organism performs the predicted environmental transformation.

Recommended downstream validation may include:

- manual sequence/domain inspection
- comparison with characterized homologues
- genomic-context analysis
- pathway-level assessment
- expression analysis
- biochemical assays
- phenotype-based bioremediation experiments

---

## Reference library

The bundled reference library is curated specifically around bioremediation-associated biological functions.

Reference metadata can include information such as:

- accession
- protein/enzyme name
- organism
- biological family
- bioremediation category
- pathway/function
- evidence level
- curation information

The reference library should be treated as a versioned scientific resource and expanded conservatively using traceable evidence.

---

## Reproducibility

For reproducible analyses, record:

- Bioremediation Gene Miner version
- reference database version
- reference metadata version
- annotation-rules version
- DIAMOND version
- whether automated InterPro analysis was enabled
- date of analysis

Predictions can change when the reference library, annotation rules, external databases, or evidence thresholds change.

---

## Version 0.3.7

Major changes in v0.3.7 include:

- automated EMBL-EBI InterProScan submission for selected Weak/Review candidates
- structured `InterPro_status`
- InterPro evidence resolver
- explicit `SUPPORTING`, `BROAD/RELATED`, `CONFLICTING`, and `NO_HITS` interpretations
- compact PANTHER, Pfam, CDD, and NCBIfam evidence in the main report
- integrated InterPro entries and GO evidence
- separate `InterPro_Full_Report.xlsx` for verbose InterPro results
- revised candidate-facing workbook organization
- dedicated `Microplastic_candidates` output retained for polymer-associated screening
- improved separation of computational evidence from biological interpretation

---

## Citation

If you use Bioremediation Gene Miner, please cite the software using the repository `CITATION.cff` metadata.

DOI:

**10.5281/zenodo.22046925**

---

## License

This project is distributed under the MIT License. See `LICENSE` for details.

---

## Development status

Bioremediation Gene Miner v0.3.7 is a **research prototype**.

The software is intended to support exploratory bacterial WGS analysis, candidate prioritization, and hypothesis generation. Results should be independently reviewed and experimentally validated before biological conclusions are made.
