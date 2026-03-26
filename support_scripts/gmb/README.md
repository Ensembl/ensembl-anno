# Gene Model Builder

A robust, configurable Python pipeline for generating high-quality consensus gene models by integrating transcriptomic assemblies (Scallop, StringTie), *ab initio* predictions (Helixer), and protein alignment evidence (OrthoDB, UniProt).

Originally developed for eukaryotic genomes, the pipeline has been heavily optimized with a "fungal-friendly" configuration preset, handling compact genomes with single-exon genes and short ORFs accurately.

---

## Features

*   **Evidence Integration:** Merges and reconciles models from Scallop, StringTie, and Helixer.
*   **Protein Evidence Support:** Uses OrthoDB and UniProt alignments to confirm models and identify true open reading frames.
*   **Pre-Genebuild Evidence Filtering:** Radically reduces noise before clustering by filtering out competing fragments, redundant protein alignments, short artifactual Helixer models, and transcriptomic chimeras (e.g., UTR-joined genes).
*   **Configurable Scoring System:** Selects the best isoform for a locus based on a weighted sum of protein evidence, *ab initio* support, and transcriptomic support.
*   **Fungal Optimization:** Defaults cater to fungal biology (e.g., lower `min_codons` threshold of 33, `allow_single_exon` enabled, stringent chimera intron length limits).
*   **100% Python:** No Perl or external script dependencies for FASTA extraction or CDS annotation.
*   **Comprehensive Output:** Generates GFF3 annotations, cDNA/CDS/Protein FASTA files, and detailed summary metrics (JSON/TSV).
*   **Annotation Validation:** Built-in scripts to compare generated consensus annotations against community references (e.g., GenBank), with locus-level classification and visualization tools.

---

## 🛠️Installation

Requirements:
*   Python 3.8+
*   `pandas`
*   `pyranges`
*   `biopython`
*   `pyyaml`
*   `matplotlib` (for QC plotting)

```bash
pip install pandas pyranges biopython pyyaml matplotlib
```

---

## Usage

### Quickstart

If you just cloned the repository, we recommend running the provided smoke test to ensure all dependencies are met.

```bash
# Install the package and dependencies
pip install -e .

# Run the smoke test
cd examples
./run_smoke_test.sh
```

### Inputs expected

The pipeline expects GFF3/GTF files for the evidence tracks (Scallop, StringTie, Helixer, OrthoDB, UniProt) and a standard FASTA file for the genome. All evidence intervals must share coordinate systems with the genome FASTA (use `--assembly-report` if seq-names differ).

### Outputs produced

All output files are generated in the specified `--output-dir` (e.g., `output/`). See the **Main Pipeline** section below for a breakdown.

### Main Pipeline

The core entry point is `gene_model_builder.py`.

```bash
python gene_model_builder.py \
    --scallop scallop_annotation.gtf \
    --stringtie stringtie_annotation.gtf \
    --helixer helixer_remapped.gff3 \
    --orthodb orthodb_annotation.gtf \
    --uniprot uniprot_annotation.gtf \
    --genome reference_genome.fa \
    --config default_config.yaml \
    --output-dir output/
```

**Output Files (`output/`):**
*   `consensus_genes.gff3`: Final structural annotation with CDS and UTR features.
*   `cdna.fa`: Spliced transcript sequences.
*   `cds.fa`: Coding sequences.
*   `prot.fa`: Translated protein sequences.
*   `report.json` / `report.tsv`: Detailed metrics on retained/filtered evidence and final gene counts.

### Configuration (`configs/fungi_default.yaml`)

The pipeline behaviour is deeply customizable via the YAML configuration file. The default configuration uses the `fungi_default.yaml` preset, which explicitly locks in rules optimized for fungal genes, but the merging logic allows overriding specific keys.

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

### Output Validation & Comparison

Use `compare_annotations.py` to evaluate your consensus against a reference (like GenBank).

```bash
python compare_annotations.py \
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

| File | Description |
| :--- | :--- |
| `gene_model_builder.py` | Main orchestrator script running the 15-step pipeline. |
| `config.py` & `.yaml` | Clade-specific pipeline configuration and dataclass definitions. |
| `evidence_filter.py` | Logic for purging noise (fragments, chimeras, poor ab initio models). |
| `scoring.py` | Evaluating gene models and selecting the top isoform(s) per locus. |
| `annotate_cds_utrs.py` | Deriving ORF, CDS boundaries, UTR regions, and handling partial edges. |
| `fasta_export.py` | Exporting strand-aware sequences (cDNA, CDS, Protein). |
| `reporting.py` | Compiling filtering stats and pipeline metrics. |
| `compare_annotations.py`| Validation tool for scoring consensus vs reference + Locus plot generation. |
| `visualize_disagreements.py`| Diagnostic tool for visualizing locus models during pipeline dev. |

---

## Testing

The pipeline is fully covered by a custom test suite.

```bash
# Run all tests (requires pytest)
pytest test_*.py
```

Test coverage includes:
*   `test_annotate_cds_utrs.py`: Core ORF finding and CDS logic.
*   `test_evidence_filter.py`: Filtering logic for proteins, chimeras, and Helixer.
*   `test_scoring.py`: Isoform selection and weighting.
*   `test_fasta_export.py`: Reverse-strand sequence extraction and FASTA header formatting.
*   `test_config.py`: YAML parsing and preset overrides.
*   `test_integration.py`: End-to-end synthetic pipeline execution.

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

All three main scripts (`gene_model_builder.py`, `compare_annotations.py`, `visualize_disagreements.py`) support consistent CLI flags for quick, reproducible testing on subsets:

### Region-based subsetting

```bash
# Run builder on a single chromosome
python gene_model_builder.py ... --seqname '1'

# Run builder on a specific region
python gene_model_builder.py ... --region '1:100000-200000'

# Use a file with multiple regions
python gene_model_builder.py ... --regions-file my_regions.txt
```

### Locus-based random sampling

```bash
# Compare on 50 random loci (reproducible with --seed)
python compare_annotations.py ... --sample-loci 50 --seed 42

# Sample from reference genes only
python compare_annotations.py ... --sample-loci 50 --sample-from reference

# Visualize 10 random structural mismatches
python visualize_disagreements.py ... --sample-loci 10 \
    --sample-from-category Structural_Mismatch
```

### Seqname mapping

Mapping is applied **before** subsetting, so region names reference the mapped seqnames:

```bash
# Map GenBank accessions to chromosome numbers, then subset chr 1
python compare_annotations.py ... \
    --assembly-report assembly_report.txt \
    --seqname '1'

# Custom mapping (TSV: from_seqname → to_seqname), overrides assembly-report
python gene_model_builder.py ... \
    --seqname-map custom_mapping.tsv \
    --region '1:100000-200000'
```

When subsetting is active, a `subset_regions.tsv` manifest is written to the output directory recording the selected regions and random seed.

