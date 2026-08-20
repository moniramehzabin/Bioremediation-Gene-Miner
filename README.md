# Bioremediation Gene Miner

**Version 0.3.5 — Research Prototype**

Bioremediation Gene Miner is a Python/DIAMOND workflow for mining bacterial whole-genome annotations for **annotated and previously hypothetical proteins** potentially associated with environmental bioremediation.

The workflow screens genes related to:

- azo-dye degradation and oxidative dye transformation
- Cr(VI) reduction/resistance and other heavy metals
- aromatic compound degradation
- pesticides and dehalogenation
- nitro compounds
- xenobiotic oxidation and cytochrome P450 systems
- microplastic/polymer biodegradation candidates
- oxidative-stress responses
- UV and DNA-repair systems

---

## Why this tool exists

Bacterial genomes may contain thousands of coding sequences (CDSs). Some potentially important bioremediation genes are already functionally annotated, whereas others remain labelled as `hypothetical protein`.

Bioremediation Gene Miner combines:

```text
Existing annotations
        ↓
Curated annotation rules
```

with:

```text
Hypothetical proteins
        ↓
DIAMOND sequence homology
        ↓
Query + biological-family retention
        ↓
Identity and coverage assessment
        ↓
Optional InterPro domain evidence
        ↓
Evidence-ranked candidate report
```

The workflow is intended for **candidate discovery and prioritization**, not automatic experimental functional assignment.

---

## Important scientific limitation

A sequence or domain match is **not experimental proof** that a protein performs a particular pollutant-degradation reaction.

The workflow therefore separates several concepts in the report, including:

- evidence class
- confidence
- match strength
- decision
- evidence source
- sequence identity
- query coverage
- reference coverage
- InterPro domain evidence

Predicted proteins should be interpreted as **candidates requiring further validation**.

Broad protein families such as cytochrome P450s, esterases, lipases, laccases, peroxidases and oxidoreductases must not automatically be interpreted as evidence of a particular pollutant-degradation phenotype.

---

# Evidence preservation in v0.3.5

A central feature of v0.3.5 is preservation of potentially informative sequence signals.

### Weak matches are retained

Weak DIAMOND signals are **not silently discarded**.

They are explicitly reported as:

```text
Weak match
```

where appropriate.

A Weak match may therefore remain available for later domain analysis, genomic-context analysis, literature comparison or experimental validation.

### Match strength and biological decision are separate

For example:

```text
Confidence: Weak
Match strength: Weak match
Decision: Review
```

does not mean that the sequence similarity is absent.

It means that the available sequence evidence is insufficient for a strong specific functional assignment and should be reviewed.

---

# Family-level DIAMOND retention

Earlier approaches that retained only the single best DIAMOND hit for each hypothetical protein could hide biologically interesting alternative family signals.

v0.3.5 instead retains the best supported reference at the:

```text
query protein + curated biological family
```

level.

Therefore, when one hypothetical protein shows similarity to more than one curated functional family, relevant family-level evidence can remain visible rather than being eliminated solely because another reference produced a better global score.

This is particularly useful when screening multifunctional or evolutionarily related oxidoreductase families.

---

# Independent InterPro evidence

InterProScan results can optionally be integrated into the report.

Importantly:

**InterPro evidence does not overwrite the original DIAMOND classification.**

For example, a protein may have:

```text
DIAMOND:
Weak match / Review

InterPro:
Strong conserved-domain evidence for a broader protein architecture
```

Both observations are retained.

The workflow does **not** automatically convert the Weak DIAMOND assignment to High confidence.

This distinction helps separate:

1. evidence for a **specific curated reference-family assignment**, and
2. evidence for the **broader domain architecture of the protein**.

For multidomain proteins, InterPro may identify combinations such as:

- cytochrome P450 domains
- FAD/NAD(P)-binding domains
- pyridine nucleotide-disulfide oxidoreductase-related domains
- flavodoxin/flavoprotein domains
- other conserved catalytic or structural domains

Such domain evidence can substantially improve biological interpretation while still preserving uncertainty about the exact substrate-specific function.

---

# Requirements

- Python 3.10+
- DIAMOND
- Windows, Linux or macOS
- Python packages listed in `requirements.txt`

Install the Python dependencies:

```bash
py -m pip install -r requirements.txt
```

On Linux/macOS, use `python3` instead of `py` where appropriate.

---

# Building the reference database

First build the curated reference FASTA:

```bash
py build_reference_db_v0.2.1.py
```

Then build the DIAMOND database:

```bash
diamond makedb --in bioremediation_reference_v0.2.1.faa --db bioremediation_reference_v0.2.1
```

---

# Running Bioremediation Gene Miner

Example:

```bash
py bioremediation_gene_miner.py genome.gbk ^
  --reference-db bioremediation_reference_v0.2.1 ^
  --reference-metadata bioremediation_reference_v0.2.1_metadata.tsv ^
  --annotation-rules annotation_rules.tsv ^
  --outdir results
```

The program:

1. parses CDS features from the GenBank annotation
2. screens existing gene/product annotations
3. extracts hypothetical proteins
4. performs DIAMOND similarity searches
5. calculates sequence identity
6. calculates query coverage
7. calculates reference coverage
8. retains relevant query + family-level evidence
9. classifies sequence evidence
10. preserves Weak signals for review
11. produces Excel and TSV reports
12. generates `interpro_candidates.faa` for optional domain analysis

---

# Optional InterProScan integration

After the first run, upload:

```text
results/interpro_candidates.faa
```

to InterProScan.

Download the InterPro TSV result and rerun the miner:

```bash
py bioremediation_gene_miner.py genome.gbk ^
  --reference-db bioremediation_reference_v0.2.1 ^
  --reference-metadata bioremediation_reference_v0.2.1_metadata.tsv ^
  --annotation-rules annotation_rules.tsv ^
  --interpro interpro_results.tsv ^
  --outdir results_with_interpro
```

InterPro evidence is added as an **independent evidence layer**.

The underlying DIAMOND identity, coverage, confidence, match strength and decision remain visible.

---

# Main report

The principal output is:

```text
Bioremediation_Gene_Miner_Report.xlsx
```

The workbook can include sheets such as:

- `Summary`
- `All_candidates`
- `High_confidence`
- `Annotated_candidates`
- `Predicted_hypothetical`
- `Family_level_candidates`
- `InterPro_supported`
- `Microplastic_candidates`
- `Supporting_functions`
- `Review_required`
- `All_DIAMOND_hits`
- `InterPro_evidence`

Companion TSV files are also generated for convenient downstream analysis.

---

# Important report fields

Depending on the evidence source, candidate tables may contain:

```text
Locus_tag
Family_target
Gene
Product_or_prediction
Origin
Category
Evidence_class
Confidence
Match_strength
Decision
Evidence_source
Identity_pct
Query_coverage_pct
Reference_coverage_pct
Evalue
Bitscore
Best_reference
InterPro_detected
InterPro_analyses
InterPro_accessions
InterPro_domain_support
Interpretation
```

These fields are intentionally kept separate so that strong domain evidence does not conceal weak sequence-level evidence, and vice versa.

---

# Development validation

The workflow has been technically tested using two bacterial whole-genome annotations.

During development, regression testing was used to confirm that evidence-integration changes did not unexpectedly alter previously recovered candidate sets.

Testing demonstrated that:

- family-level DIAMOND signals could be retained
- Weak signals remained visible
- optional InterPro evidence could be added without removing underlying DIAMOND rows
- candidate counts remained stable during the final evidence-integration tests

These development tests demonstrate **technical workflow consistency only**.

They do **not** establish biological sensitivity, specificity or experimental confirmation of predicted functions.

---

# Interpretation of predictions

Recommended terminology includes:

> sequence-supported candidate

> family-level candidate

> Weak sequence match requiring review

> InterPro-supported domain architecture

Avoid describing a predicted gene as experimentally confirmed solely on the basis of DIAMOND or InterPro evidence.

Where biologically important, candidates should be validated using approaches such as:

- conserved-domain analysis
- genomic-context analysis
- phylogenetic analysis
- structural prediction
- expression analysis
- enzyme assays
- gene knockout/complementation
- pollutant-transformation experiments

---

# Version 0.3.5

Major changes include:

- preservation of Weak DIAMOND signals
- explicit `Weak match` reporting
- separation of confidence, match strength and decision
- retention of the best reference per query + curated biological family
- improved recovery of alternative family-level evidence
- increased DIAMOND target retention for family-level screening
- independent InterPro evidence integration
- preservation of original DIAMOND classification after InterPro integration
- reporting of InterPro analyses, accessions and domain descriptions
- dedicated `InterPro_supported` output
- continued coverage safeguards against misleading partial matches

---

# Previous updates

## v0.3.3

- corrected total CDS reporting
- retained InterPro signatures as supporting evidence without automatic confidence upgrading
- verified CDS totals during WGS testing

## v0.3.2

- corrected high-confidence CLI reporting
- added product-aware annotation matching
- reduced false assignments caused by ambiguous gene symbols
- improved handling of conflicting gene/product annotations

---

# Repository status

Bioremediation Gene Miner is currently a **research prototype**.

Reference-database curation, benchmarking and formal biological validation are ongoing.

Predictions produced by this software should be interpreted as **computational candidates rather than experimental confirmation**.
