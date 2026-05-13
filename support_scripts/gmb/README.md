# Gene Model Builder

A robust, configurable Python pipeline for generating high-quality consensus gene models by integrating transcriptomic assemblies (Scallop, StringTie), *ab initio* predictions (Helixer), and protein alignment evidence (OrthoDB, UniProt).

Originally developed for eukaryotic genomes, the pipeline has been heavily optimised with a "fungal-friendly" configuration preset, handling compact genomes with single-exon genes and short ORFs accurately.

---

## Features

*   **Evidence Integration:** Merges and reconciles models from Scallop, StringTie, and Helixer.
*   **Protein Evidence Support:** Uses OrthoDB and UniProt alignments to confirm models and identify true open reading frames.
*   **Pre-Genebuild Evidence Filtering:** Radically reduces noise before clustering by filtering out competing fragments, redundant protein alignments, short artifactual Helixer models, and transcriptomic chimeras (e.g., UTR-joined genes).
*   **Configurable Scoring System:** Selects the best isoform for a locus based on a weighted sum of protein evidence, *ab initio* support, and transcriptomic support.
*   **Fungal Optimization:** Defaults cater to fungal biology (e.g., lower `min_codons` threshold of 33, `allow_single_exon` enabled, stringent chimera intron length limits).
*   **100% Python:** No Perl or external script dependencies for FASTA extraction or CDS annotation.
*   **Comprehensive Output:** Generates GFF3 annotations, cDNA/CDS/Protein FASTA files, and detailed summary metrics (JSON/TSV).
*   **Annotation Validation:** Built-in scripts to compare generated consensus annotations against community references (e.g., GenBank), with locus-level classification and visualisation tools.

---

## Installation

Requirements:
*   Python 3.8+
*   `pandas`
*   `pyranges`
*   `biopython`
*   `pyyaml`
*   `matplotlib` (for QC plotting)

```bash
# Option 1: install as an editable package (recommended for development)
cd support_scripts/gmb
pip install -e ".[dev]"

# Option 2: install dependencies only
pip install pandas pyranges biopython pyyaml matplotlib
```

---

## Usage

The canonical way to run GMB is via the Python package CLI:

```bash
python -m gmb.cli.build   --help   # Gene Model Builder
python -m gmb.cli.compare --help   # Annotation comparison
python -m gmb.cli.visualize --help # Disagreement visualisation
```

After `pip install -e .`, the same commands are available as console scripts:

```bash
gmb-build     --help
gmb-compare   --help
gmb-visualize --help
```

**Backward-compatible wrapper scripts** are also available under `scripts/`
(`scripts/gene_model_builder.py`, `scripts/compare_annotations.py`,
`scripts/visualize_disagreements.py`). CLI arguments are identical — these
simply delegate to the package entry points.

### Quickstart — bundled *Z. tritici* example data

The repository includes a complete set of example evidence files for *Zymoseptoria tritici*
(wheat yellow leaf blotch fungus) under `support_scripts/gmb/z_tritici/`, along with a
pre-built DIAMOND database at `support_scripts/gmb/swissprot.dmnd`.  You can run the
full pipeline immediately after cloning — no external downloads required.

```bash
cd support_scripts/gmb

# Install Python dependencies
pip install pandas "pyranges>=0.0.120,<=0.1.4" biopython pyyaml matplotlib

# Fastest smoke-test: pre-subsetted 500 kb region fixture (~8 seconds)
python -m gmb.cli.build \
    --scallop   tests/fixtures/z_tritici_region1/scallop_geneset.gtf \
    --stringtie tests/fixtures/z_tritici_region1/stringtie_geneset.gtf \
    --helixer   tests/fixtures/z_tritici_region1/helixer_remapped.gff3 \
    --orthodb   tests/fixtures/z_tritici_region1/orthodb_geneset.gtf \
    --uniprot   tests/fixtures/z_tritici_region1/uniprot_geneset.gtf \
    --genome    tests/fixtures/z_tritici_region1/genome.fa \
    --output-dir output_region1/ \
    --gene-prefix ZTRITICI

# Full chromosome 1 using --seqname 1 (~25 min; large OrthoDB file loaded in full)
python -m gmb.cli.build \
    --scallop   z_tritici/scallop_geneset.gtf \
    --stringtie z_tritici/stringtie_geneset.gtf \
    --helixer   z_tritici/helixer_remapped.gff3 \
    --orthodb   z_tritici/orthodb_geneset.gtf \
    --uniprot   z_tritici/uniprot_geneset.gtf \
    --genome    z_tritici/zymoseptoria_tritici.fa \
    --assembly-report z_tritici/GCF_000219625.1_MYCGR_v2.0_assembly_report.txt \
    --output-dir output_seqname1/ \
    --gene-prefix ZTRITICI \
    --seqname 1

# Full genome run (all 21 chromosomes, ~15–30 min)
python -m gmb.cli.build \
    --scallop   z_tritici/scallop_geneset.gtf \
    --stringtie z_tritici/stringtie_geneset.gtf \
    --helixer   z_tritici/helixer_remapped.gff3 \
    --orthodb   z_tritici/orthodb_geneset.gtf \
    --uniprot   z_tritici/uniprot_geneset.gtf \
    --genome    z_tritici/zymoseptoria_tritici.fa \
    --assembly-report z_tritici/GCF_000219625.1_MYCGR_v2.0_assembly_report.txt \
    --output-dir output_full/ \
    --gene-prefix ZTRITICI
```

**Outputs in `output_seqname1/`:**

| File | Description |
| :--- | :--- |
| `consensus.gff3` | Final structural annotation (gene / mRNA / exon / CDS / UTR) |
| `cdna.fa` | Spliced transcript nucleotide sequences — one record per mRNA in `consensus.gff3` |
| `cds.fa` | Coding sequences (absent when no CDS features were predicted) |
| `prot.fa` | Translated protein sequences — one record per CDS-bearing mRNA in `consensus.gff3` |
| `summary.json` | Pipeline metrics (gene counts, filtering stats) |
| `summary.tsv` | Same data in tabular form |
| `fasta_qc_report.json` | FASTA QC report (written when `--validate-fasta` is used) |
| `subset_regions.tsv` | Records which regions were selected (when `--seqname` is used) |

> **Note on seqname mapping:** The input GTF/GFF3 files use NCBI/GenBank accession
> numbers (e.g. `CM001642.1`) while the genome FASTA uses short integer chromosome names
> (`1`, `2`, …).  The `--assembly-report` flag handles this remapping automatically.
> Always specify `--seqname` using the *mapped* name (e.g. `--seqname 1`), not the raw
> accession.  See [KI-002](KNOWN_ISSUES.md#ki-002--seqname-mapping-is-applied-globally-before-any-subsetting)
> for details.

---

### Protein validation with bundled DIAMOND database

If DIAMOND and Psauron are installed, you can enable the protein validation stage using
the bundled `swissprot.dmnd` database.  Edit `configs/fungi_default.yaml` (or pass a
custom config) to add:

```yaml
protein_validation:
  enabled: true
  diamond_db: swissprot.dmnd   # path relative to the gmb/ working directory
  diamond_threads: 4
  psauron_threshold: 0.5
```

Then run:

```bash
python -m gmb.cli.build \
    ... \
    --config configs/fungi_default.yaml \
    --seqname 1
```

To verify that DIAMOND and Psauron are on `$PATH` before running:

```bash
python -m gmb.cli.build --check-deps
```

> **CI note:** Protein validation is skipped in CI by default.  To run it locally, set
> `RUN_EXTERNAL_TOOLS=1` before running pytest (see [Testing](#testing) below).

---

### Compare and validate the output

After building the consensus annotation, compare it against the bundled Ensembl Fungi
reference:

```bash
python -m gmb.cli.compare \
    --consensus output_seqname1/consensus.gff3 \
    --reference z_tritici/ensembl_fungi_reference.gff3 \
    --assembly-report z_tritici/GCF_000219625.1_MYCGR_v2.0_assembly_report.txt \
    --seqname 1 \
    --output-dir compare_seqname1/ \
    --scallop   z_tritici/scallop_geneset.gtf \
    --stringtie z_tritici/stringtie_geneset.gtf \
    --helixer   z_tritici/helixer_remapped.gff3 \
    --orthodb   z_tritici/orthodb_geneset.gtf \
    --uniprot   z_tritici/uniprot_geneset.gtf
```

This produces per-locus classifications (`Exact_Match`, `Partial_Match`,
`Structural_Mismatch`, `Missed`, `Novel`) plus sensitivity metrics and locus plots.

---

### Troubleshooting

**"No features found on seqname X"**
> Your evidence files likely use different chromosome/contig names from the genome FASTA.
> Pass `--assembly-report` to remap NCBI accessions to chromosome numbers, or provide a
> custom two-column TSV with `--seqname-map from_name,to_name`.

**Pipeline produces 0 genes**
> Check that `--helixer` points to `helixer_remapped.gff3` (not `helixer_geneset.gff3`).
> The raw Helixer output uses NCBI accessions; the remapped file uses short chr names.

**ImportError on pyranges**
> Install a compatible version: `pip install "pyranges>=0.0.120,<=0.1.4"`.  The pipeline
> does not support pyranges ≥ 0.2.

---

### Inputs expected

The pipeline expects GFF3/GTF files for the evidence tracks (Scallop, StringTie, Helixer, OrthoDB, UniProt) and a standard FASTA file for the genome. All evidence intervals must share coordinate systems with the genome FASTA (use `--assembly-report` if seq-names differ).

### Outputs produced

All output files are generated in the specified `--output-dir` (e.g., `output/`). See the **Main Pipeline** section below for a breakdown.

### Main Pipeline

```bash
python -m gmb.cli.build \
    --scallop scallop_annotation.gtf \
    --stringtie stringtie_annotation.gtf \
    --helixer helixer_remapped.gff3 \
    --orthodb orthodb_annotation.gtf \
    --uniprot uniprot_annotation.gtf \
    --genome genome.fa \
    --config config.yaml \
    --output-dir output/
```

**Output Files (`output/`):**
*   `consensus_genes.gff3`: Final structural annotation with CDS and UTR features.
*   `cdna.fa`: Spliced transcript sequences.
*   `cds.fa`: Coding sequences.
*   `prot.fa`: Translated protein sequences.
*   `report.json` / `report.tsv`: Detailed metrics on retained/filtered evidence and final gene counts.

### Configuration (`configs/fungi_default.yaml`)

The pipeline behaviour is deeply customisable via the YAML configuration file. The default configuration uses the `fungi_default.yaml` preset, which explicitly locks in rules optimised for fungal genes, but the merging logic allows overriding specific keys.

Key configurable areas:
*   `orf.min_codons`: Minimum ORF length (default 33 for fungi).
*   `protein_filter`: Thresholds for dropping fragmented or poorly supported protein alignments. **OrthoDB filters**: now configurable via `min_alignment_coverage`, `min_percent_identity`, and `min_bitscore` here.
*   `transcriptomic_filter`: Thresholds for identifying artificially merged transcript models (e.g., max intron length).
*   `helixer_filter`: Rules for filtering unreliable *ab initio* models.
*   `scoring`: Weights determining how isoforms are selected (e.g., prioritizing protein support/configurable thresholds for keeping Helixer models).
*   `protein_validation`: Enable a batch validation stage (`enabled: true`) that writes translated candidate models to a FASTA, runs DIAMOND and Psauron, and injects a derived `protein_coding_score` back into the candidate gating logic.
*   `utr`: Thresholds and criteria for keeping or trimming Untranslated Regions. Features robust "end support" validation (`require_end_support: true`), allowing you to configure exactly which evidence assemblies must agree on a transcript's start/end coordinates (`end_support_sources`), with a distance tolerance (`end_tolerance_bp`). If the transcript edges lack multi-source agreement or protein validation (depending on `end_support_mode` policy), they can fallback to `drop_utr`, `hard_cap` to a hardcoded base definition, or `drop_transcript` fully based on `fallback_policy_when_unsupported`.

**Config Merge Rules:**
Dicts deep-merge, lists replace entirely, unknown keys raise an error to catch typos.

### How prot.fa and cdna.fa are generated

Both FASTA files are derived from the **post-processed** `consensus.gff3` and are guaranteed to be consistent with it.

**`cdna.fa`**
- One record per `mRNA` feature in `consensus.gff3`.
- Sequence is the spliced cDNA: exonic genomic sequence concatenated in 5′→3′ transcript order, reverse-complemented for minus-strand transcripts.
- FASTA header: `>{transcript_id}` (matches the `ID=` attribute of the mRNA row).

**`prot.fa`**
- One record per `mRNA` that has at least one `CDS` child in `consensus.gff3`.
- Sequence is the translated protein from the CDS interval(s). Terminal stop codons are stripped.
- FASTA header: `>{transcript_id}` (same stable ID as the mRNA row).
- UTR-only transcripts (no CDS rows) are excluded from `prot.fa` but still appear in `cdna.fa`.

**ID guarantee**
Every `>{id}` in `prot.fa` and `cdna.fa` maps to exactly one `mRNA` row in `consensus.gff3`. No transcript appears in the FASTA files unless it survived all GFF post-processing steps (structural validation, deduplication).

**FASTA QC**

Run the standalone QC tool to verify output integrity:

```bash
# Coverage checks only
python -m gmb.pipeline.fasta_qc output/

# Coverage + sequence correctness (reconstructs from genome)
python -m gmb.pipeline.fasta_qc output/ --genome genome.fa

# Or via the pipeline with --validate-fasta
python -m gmb.cli.build ... --validate-fasta
```

This writes `fasta_qc_report.json` summarising:
- transcript/protein/cDNA record counts
- missing or extra FASTA records (relative to GFF)
- duplicate headers
- protein and cDNA sequence reconstruction mismatches (when `--genome` is supplied)

### Output Validation & Comparison

Use `gmb.cli.compare` to evaluate your consensus against a reference (like GenBank).

```bash
python -m gmb.cli.compare \
    --consensus output/consensus_genes.gff3 \
    --reference GenBank_annotation.gff3.gz \
    --assembly-report assembly_report.txt \
    --output-dir validation/ \
    --plots-per-category 3 \
    --scallop scallop_annotation.gtf \
    --stringtie stringtie_annotation.gtf \
    --helixer helixer_remapped.gff3 \
    --orthodb orthodb_annotation.gtf \
    --uniprot uniprot_annotation.gtf
```

This generates:
*   Locus-level classifications (`Exact_Match`, `Partial_Match`, `Structural_Mismatch`, `Missed`, `Novel`).
*   Sensitivity metrics and tabular reports.
*   Visual plots of representative loci comparing all evidence tracks side-by-side.

---

## Project Structure

The codebase is organised as a proper Python package (`gmb/`) with clear
separation between pipeline logic, comparison tools, shared utilities, and CLI
entry points.

```
gmb/                         # Installable Python package
  __init__.py                # Package root (__version__ = "2.0.0")
  pipeline/                  # Core pipeline modules
    builder.py               # Main orchestrator (15-step pipeline)
    config.py                # YAML config loading & dataclass hierarchy
    evidence_filter.py       # Noise removal (fragments, chimeras, etc.)
    scoring.py               # Isoform scoring & selection
    annotate_cds_utrs.py     # ORF/CDS/UTR derivation
    fasta_export.py          # Strand-aware FASTA extraction
    fasta_qc.py              # FASTA QC checks
    gff3_validate.py         # GFF3 structural validation
    dedup_genes.py           # Gene deduplication
    protein_validation.py    # DIAMOND + Psauron scoring
    reporting.py             # Summary metrics (JSON/TSV)
    subset_utils.py          # Region/seqname subsetting
  compare/                   # Annotation comparison tools
    compare_annotations.py   # Consensus vs reference classification
    visualize_disagreements.py  # Locus-level diagnostic plots
    validate_annotation.py   # Standalone annotation validator
  utils/                     # Shared helpers (no duplicated code)
    intervals.py             # Interval overlap & merge functions
    fasta.py                 # Genome loading, FASTA I/O
    gff.py                   # GFF3 parsing utilities
    io.py                    # Directory helpers
    logging.py               # Logging setup
  cli/                       # CLI entry points (installed via pip)
    build.py                 # python -m gmb.cli.build / gmb-build
    compare.py               # python -m gmb.cli.compare / gmb-compare
    visualize.py             # python -m gmb.cli.visualize / gmb-visualize
scripts/                     # Backward-compatible wrapper scripts
  gene_model_builder.py      # -> gmb.pipeline.builder
  compare_annotations.py     # -> gmb.compare.compare_annotations
  visualize_disagreements.py # -> gmb.compare.visualize_disagreements
tools/                       # Optional helper utilities
  remap_helixer.py           # Remap Helixer GFF3 seq IDs via assembly report
  make_fixed_prot.py         # Re-translate Helixer CDS with corrected sort
  retranslate_from_gff3.py   # Re-translate all CDS from consensus GFF3
  check_lengths.py           # Quick transcript length sanity checks
exploratory/                 # Debugging & one-off analysis scripts
  debug_scoring.py           # Interactive scoring inspection
  gen_extra_plots.py         # Generate additional comparison plots
  test_clustering.py         # Ad-hoc clustering experiments
```

---

## Optional Tools

Helper scripts live in `tools/` and are not required for normal pipeline use:

| Script | Purpose |
| :----- | :------ |
| `tools/remap_helixer.py` | Remap Helixer GFF3 sequence IDs using an NCBI assembly report |
| `tools/make_fixed_prot.py` | Re-translate Helixer-sourced CDS with corrected strand sort |
| `tools/retranslate_from_gff3.py` | Re-translate all CDS from a consensus GFF3 + genome |
| `tools/check_lengths.py` | Quick transcript length sanity checks on GFF3/GTF files |

Ad-hoc debugging and plotting scripts are in `exploratory/` — these are not
maintained for general use.

---

## Testing

The pipeline ships with a pytest test suite covering unit logic, synthetic integration,
and a real-data subset integration test using the bundled *Z. tritici* example data.

```bash
cd support_scripts/gmb

# Install dev dependencies
pip install pytest

# Run all fast tests (~30 s)
pytest tests/ -q

# Run only integration tests (includes the z_tritici seqname-1 subset, ~2 min)
pytest tests/ -m integration -v

# Run external-tool tests (requires DIAMOND + Psauron on $PATH)
RUN_EXTERNAL_TOOLS=1 pytest tests/ -m external_tools -v
```

Key test modules:

| Module | What it tests |
| :--- | :--- |
| `tests/test_integration.py` | End-to-end run on a tiny synthetic 501 bp dataset |
| `tests/test_z_tritici_subset.py` | Real-data smoke test on bundled *Z. tritici* chr 1 |
| `tests/test_scoring.py` | Isoform scoring and selection logic |
| `tests/test_evidence_filter.py` | Chimera / protein / Helixer filtering logic |
| `tests/test_annotate_cds_utrs.py` | ORF finding, CDS boundary derivation |
| `tests/test_fasta_export.py` | Strand-aware FASTA sequence extraction |
| `tests/test_config.py` | YAML config loading and preset merging |
| `tests/test_compare_annotations.py` | Locus classification logic |
| `tests/test_protein_validation.py` | Protein scoring (skipped without external tools) |

**Golden regression fixtures** (optional, generated once and stored under
`tests/fixtures/expected/`) enable exact locus-level regression testing on the
*Z. tritici* chr-1 subset.  Generate or refresh them with:

```bash
python tests/generate_golden_fixtures.py
```

---

## Pipeline Logic Diagram

1. **Load Transcriptomics** (Scallop, StringTie)
2. **Load Ab Initio** (Helixer)
3. **Load Protein Aligns** (OrthoDB, UniProt)
4. **Pre-Genebuild Filtering**:
   * *Drop fragmented proteins (competing with better models).*
   * *Drop spurious Helixer models (single-exon, no start/stop, low support).*
   * *Drop transcriptomic chimeras (massive introns).*
5. **Transcript Splitting** *(optional, config-driven)*: Split mega-transcripts into local segments.
6. **Cluster Loci**: Group overlapping evidence using PyRanges.
7. **Score & Select**: Pick top isoforms based on configurable weightings.
8. **Annotate**: Deduce CDS/UTRs for selected models via longest-ORF discovery.
9. **Export**: Write GFF3 and FASTA sets.
10. **Report**: Write summary metrics.

---

## Transcript Splitting

### Why it exists

Transcriptome assemblers (especially StringTie `--merge`) can produce pathological
"mega-transcripts" — merged super-loci like `MSTRG.1.*` spanning tens or hundreds
of kb. Without Helixer evidence, the validation step drops most of these, causing
extremely low gene counts. Transcript splitting breaks them into multiple local
candidate segments using exon geometry, so downstream clustering and selection
operate on realistic segments rather than discarding entire evidence tracks.

### How to enable

In your config YAML (or via override):

```yaml
transcript_splitting:
  split_enabled: true       # enable splitting
  split_gap_bp: 3000        # gap threshold (default = max_intron_length)
  split_on_large_exon_bp: 15000  # drop transcripts with any exon > this
  max_segments_per_transcript: 50  # safety: drop transcript if exceeding this
```

Set `split_enabled: false` (the default for `fungi_default.yaml`) to disable.

All thresholds are config-driven — no hard-coded species constants.

---

## Supported Dependency Versions

| Package    | Base (`requirements.txt`)  | Compat target (`requirements-compat.txt`) |
| :--------- | :------------------------- | :---------------------------------------- |
| pandas     | `>=2.0,<3`                 | `==3.0.0`                                 |
| pyranges   | `>=0.0.120,<=0.1.4`        | `==0.1.4`                                 |
| biopython  | latest                     | latest                                    |
| pyyaml     | latest                     | latest                                    |
| matplotlib | latest                     | latest                                    |

CI runs both the default and compat environments.

---

## Fast Test Mode

All three main commands (`gmb.cli.build`, `gmb.cli.compare`, `gmb.cli.visualize`) support consistent CLI flags for quick, reproducible testing on subsets:

### Region-based subsetting

```bash
# Run builder on a single chromosome
python -m gmb.cli.build ... --seqname '1'

# Run builder on a specific region
python -m gmb.cli.build ... --region '1:100000-200000'

# Use a file with multiple regions
python -m gmb.cli.build ... --regions-file my_regions.txt
```

### Locus-based random sampling

```bash
# Compare on 50 random loci (reproducible with --seed)
python -m gmb.cli.compare ... --sample-loci 50 --seed 42

# Sample from reference genes only
python -m gmb.cli.compare ... --sample-loci 50 --sample-from reference

# Visualise 10 random structural mismatches
python -m gmb.cli.visualize ... --sample-loci 10 \
    --sample-from-category Structural_Mismatch
```

### Seqname mapping

Mapping is applied **before** subsetting, so region names reference the mapped seqnames:

```bash
# Map GenBank accessions to chromosome numbers, then subset chr 1
python -m gmb.cli.compare ... \
    --assembly-report assembly_report.txt \
    --seqname '1'

# Custom mapping (TSV: from_seqname → to_seqname), overrides assembly-report
python -m gmb.cli.build ... \
    --seqname-map custom_mapping.tsv \
    --region '1:100000-200000'
```

When subsetting is active, a `subset_regions.tsv` manifest is written to the output directory recording the selected regions and random seed.

