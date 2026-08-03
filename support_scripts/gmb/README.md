# Gene Model Builder

A configurable Python pipeline for generating consensus eukaryotic gene models by
integrating transcriptomic assemblies (Scallop, StringTie), *ab initio* predictions
(Helixer or Tiberius), and protein alignment evidence (OrthoDB, UniProt, GenBlast).
Optional protein validation (DIAMOND + Psauron) and canonical transcript selection
are available as post-build steps.

---

## Requirements

| Dependency | Version |
| :--------- | :------ |
| Python | **≥ 3.10** (production cluster baseline is 3.10.x) |
| pandas | `>=2.0,<3` (compat-tested to `==3.0.0` on Python ≥ 3.11) |
| pyranges | `>=0.0.120,<=0.1.4` |
| biopython | latest |
| pyyaml | latest |
| matplotlib | latest (for QC plotting) |

External tools are **optional** — the core pipeline runs without them:

| Tool | Used for |
| :--- | :------- |
| DIAMOND | Protein validation (`protein_validation.enabled: true`) |
| Psauron | Protein-coding score in protein validation |
| InterProScan | Canonical-choice resolver only — not required for a normal build |

---

## Installation

```bash
cd support_scripts/gmb

# Recommended: editable install (development)
pip install -e ".[dev]"

# Or: install from a built wheel
pip install gene_model_builder-2.0.0-py3-none-any.whl
```

Entry points installed: `gmb-build`, `gmb-compare`, `gmb-visualize`,
`gmb-longread-consensus`, `gmb-canonical-selection`, `gmb-interpro-review`,
`gmb-interpro-resolve`.

---

## Quickstart — bundled *Z. tritici* fixture

The repository includes a pre-subsetted 500 kb fixture under
`tests/fixtures/z_tritici_region1/` for a fast smoke test (~8 seconds):

```bash
cd support_scripts/gmb

gmb-build \
    --scallop   tests/fixtures/z_tritici_region1/scallop_geneset.gtf \
    --stringtie tests/fixtures/z_tritici_region1/stringtie_geneset.gtf \
    --helixer   tests/fixtures/z_tritici_region1/helixer_remapped.gff3 \
    --orthodb   tests/fixtures/z_tritici_region1/orthodb_geneset.gtf \
    --uniprot   tests/fixtures/z_tritici_region1/uniprot_geneset.gtf \
    --genome    tests/fixtures/z_tritici_region1/genome.fa \
    --output-dir output_region1/ \
    --gene-prefix ZTRITICI
```

The fixture has no corresponding full-genome data in this repository; the subset
fixture is self-contained and sufficient to verify installation.

---

## Basic usage

```bash
gmb-build \
    --scallop   scallop.gtf \
    --stringtie stringtie.gtf \
    --helixer   helixer_remapped.gff3 \
    --orthodb   orthodb.gtf \
    --uniprot   uniprot.gtf \
    --genome    genome.fa \
    --config    my_species.yaml \
    --output-dir output/
```

All evidence files must share coordinate systems with the genome FASTA.
If sequence names differ (e.g. NCBI accessions vs. short chromosome names),
pass `--assembly-report` to remap automatically, or `--seqname-map` for a
custom two-column TSV.

**Subsetting a single chromosome:**

```bash
gmb-build ... --seqname 1
```

---

## Configuration

GMB assembles the effective configuration in layers:

```text
standard.yaml (organism-neutral base, always loaded)
    → clade preset  (--preset fungi | apicomplexa | none)
        → user-supplied --config files (in order)
```

The default preset is `fungi`.  Use `--list-presets` to see what is
installed.

```bash
# Apicomplexa preset with local paths layered on top
gmb-build \
    --preset apicomplexa \
    --config local_cluster_paths.yaml \
    ...

# Standard base only, no clade preset
gmb-build --preset none --config my_species.yaml ...

# List installed clade presets
gmb-build --list-presets
```

**A misspelled or absent `--config` path raises an error immediately** — GMB
never silently falls back.

Every run writes `resolved_config.yaml` and `resolved_config_sha256` to the
output directory so runs are reproducible.

Example config files in `configs/`:

| File | Purpose |
| :--- | :------- |
| `configs/apicomplexa_first_pass.yaml` | Apicomplexa delta (use `--preset apicomplexa` instead) |
| `configs/apicomplexa_chr1_protein_validation.example.yaml` | Protein-validation overlay example |
| `configs/ebi_protein_validation.example.yaml` | EBI cluster protein-validation paths |

For the full layering model, preset descriptions, and config schema see
**[docs/build_and_configuration.md](docs/build_and_configuration.md)**.

---

## Protein validation

DIAMOND and Psauron must be installed separately — they are **not bundled**.
The DIAMOND database must also be provided by the user (`swissprot.dmnd` is **not
included** in the repository or the wheel).

To enable protein validation, set in your config YAML:

```yaml
protein_validation:
  enabled: true
  diamond_db: /path/to/swissprot.dmnd   # required; no default
  diamond_weight: 0.7
  psauron_weight: 0.3
  min_score: 0.7
  policy: penalize   # "drop" | "penalize" (also accepts "penalise")
```

Verify tools are on `$PATH` before a long run:

```bash
gmb-build --check-deps
```

Per-transcript DIAMOND/Psauron results are written to `protein_validation.tsv` in
the output directory.

---

## Long-read consensus preprocessing

Raw Minimap2 per-read GTF alignments are too noisy to use directly as `--minimap2`
evidence. Collapse them first with `gmb-longread-consensus`:

```bash
# A species preset or --config is required for Stage 2 consensus
gmb-longread-consensus \
    --input     raw_minimap2.gtf \
    --output-dir longread_consensus_out/ \
    --preset    pfalciparum_pure
```

This writes `minimap2_consensus.gtf` (source label `Minimap2Consensus`) that can
then be passed to `gmb-build --minimap2 longread_consensus_out/minimap2_consensus.gtf`.

Available presets:

```bash
gmb-longread-consensus --help   # lists available presets
```

See **[docs/longread_consensus.md](docs/longread_consensus.md)** for the full
config schema, preset descriptions, short-read rescue modes, and cluster
(Slurm array) usage.

---

## Canonical transcript selection

For multi-isoform genes, `gmb-canonical-selection` is a non-destructive post-build
step that ranks isoforms and picks one representative per gene:

```bash
gmb-canonical-selection \
    --consensus-gff3        output/consensus.gff3 \
    --evidence-attribution  output/evidence_attribution.tsv \
    --protein-validation    output/protein_validation.tsv \
    --output-dir            output/canonical_selection/
```

Writes `canonical_transcripts.tsv`, `transcript_ranking.tsv`, and
`canonical_selection_summary.json`. Selection uses a deterministic priority order
(complete ORF → protein-validation support → biological evidence-class breadth → GMB
score → CDS length → transcript ID). The `canonical_total_score` is continuous and
reflects *how much better* the winner is, not *what picked it*.

See **[docs/canonical_selection.md](docs/canonical_selection.md)**.

---

## InterProScan resolver (optional second stage)

For the small subset of genes where GMB could not choose confidently,
`gmb-interpro-review` + `gmb-interpro-resolve` provides a second-stage canonical
resolver backed by InterProScan domain evidence. InterProScan is **never required**
for a normal build.

See **[docs/interpro_resolver.md](docs/interpro_resolver.md)** for the full
config schema, replacement policy, and cluster execution examples.

---

## Outputs

All output files are written to `--output-dir`.

| File | Description |
| :--- | :---------- |
| `consensus.gff3` | Final structural annotation (gene / mRNA / exon / CDS / UTR) |
| `cdna.fa` | Spliced transcript sequences — one record per mRNA |
| `cds.fa` | Coding sequences (absent when no CDS features were predicted) |
| `prot.fa` | Translated proteins — one record per CDS-bearing mRNA |
| `summary.json` | Pipeline metrics (gene counts, filtering statistics) |
| `summary.tsv` | Same data in tabular form |
| `evidence_attribution.tsv` | Per-transcript evidence source and score provenance |
| `protein_validation.tsv` | Per-transcript DIAMOND/Psauron results (when enabled) |
| `collapsed_duplicate_transcripts.tsv` | Log of exact-duplicate collapses (when any occur) |
| `fasta_qc_report.json` | FASTA QC report (when `--validate-fasta` is used) |
| `subset_regions.tsv` | Regions selected (when `--seqname` / `--region` is used) |

Every `>id` in `prot.fa` and `cdna.fa` maps to exactly one `mRNA` row in
`consensus.gff3`. See **[docs/output_contracts.md](docs/output_contracts.md)** for
the complete output contract.

**FASTA QC:**

```bash
# Coverage only
python -m gmb.pipeline.fasta_qc output/

# Coverage + sequence reconstruction (requires genome FASTA)
python -m gmb.pipeline.fasta_qc output/ --genome genome.fa
```

---

## Testing

```bash
cd support_scripts/gmb

# Install dev dependencies
pip install -e ".[dev]"

# All fast tests (~30 s)
pytest tests/ -q -m "not integration"

# Integration tests on the bundled z_tritici region fixture (~8 s)
pytest tests/test_z_tritici_subset.py -v -m integration

# Protein validation tests (requires DIAMOND + Psauron on $PATH)
RUN_EXTERNAL_TOOLS=1 pytest tests/ -m external_tools -v
```

Key test modules:

| Module | What it covers |
| :----- | :------------- |
| `tests/test_z_tritici_subset.py` | End-to-end on the bundled region fixture; golden regression |
| `tests/test_integration.py` | Synthetic 501 bp dataset |
| `tests/test_scoring.py` | Isoform scoring and selection |
| `tests/test_evidence_filter.py` | Chimera / protein / backbone filtering |
| `tests/test_annotate_cds_utrs.py` | ORF finding and CDS derivation |
| `tests/test_config.py` | YAML config loading, preset merging, validation |
| `tests/test_longread_rescue.py` | Long-read consensus + short-read rescue |
| `tests/test_canonical_selection.py` | Canonical transcript ranking and selection |
| `tests/test_protein_validation.py` | Protein scoring (skipped without external tools) |

Regenerate golden regression fixtures after intentional output changes:

```bash
python tests/generate_golden_fixtures.py
```

---

## Project structure

```
gmb/                           # Installable Python package
  __init__.py                  # Package root (__version__ = "2.0.0")
  configs/                     # Bundled preset YAML files (package data)
    standard.yaml              # Organism-neutral base (always loaded)
    fungi.yaml                 # Fungal overrides over standard
    apicomplexa.yaml           # Apicomplexa overrides over standard
    longread_consensus/
      pfalciparum_pure.yaml
      pfalciparum_assisted.yaml
  pipeline/                    # Core pipeline logic
    builder.py                 # Main orchestrator (15-step pipeline)
    config.py                  # YAML config loading & dataclass hierarchy
    evidence_filter.py         # Noise removal (fragments, chimeras, etc.)
    scoring.py                 # Isoform scoring & selection
    annotate_cds_utrs.py       # ORF/CDS/UTR derivation
    fasta_export.py            # Strand-aware FASTA extraction
    fasta_qc.py                # FASTA QC checks
    gff3_validate.py           # GFF3 structural validation
    dedup_genes.py             # Gene deduplication
    protein_validation.py      # DIAMOND + Psauron scoring
    canonical_evidence.py      # Evidence-class vocabulary
    canonical_selection.py     # Canonical transcript selector
    reporting.py               # Summary metrics (JSON/TSV)
    subset_utils.py            # Region/seqname subsetting
    longread/                  # Long-read consensus collapsing
      config.py                # LongreadConsensusConfig + preset loading
      consensus.py             # Read collapsing and grouping
      rescue.py                # Short-read-assisted rescue policies
      io.py                    # GTF I/O and split-by-seqname
      reporting.py             # Run manifest + summary outputs
  compare/                     # Annotation comparison tools
    compare_annotations.py
    visualize_disagreements.py
    validate_annotation.py
  utils/                       # Shared helpers
    intervals.py
    fasta.py
    gff.py
    io.py
    logging.py
  cli/                         # CLI entry points (installed by pip)
    build.py
    compare.py
    visualize.py
    longread_consensus.py
    canonical_selection.py
    interpro_review.py
    interpro_resolve.py
configs/                       # User-facing example config files (not installed)
  apicomplexa_first_pass.yaml
  apicomplexa_chr1_protein_validation.example.yaml
  ebi_protein_validation.example.yaml
docs/                          # Detailed documentation
  build_and_configuration.md
  longread_consensus.md
  canonical_selection.md
  output_contracts.md
  interpro_resolver.md
tests/                         # pytest suite
  fixtures/
    z_tritici_region1/         # Bundled 500 kb real-data fixture
tools/                         # Optional standalone helper scripts
  remap_helixer.py
  retranslate_from_gff3.py
  audit_duplicate_transcripts.py
```

---

## Supported dependency versions

| Package    | Base (`requirements.txt`)   | Compat (`requirements-compat.txt`) |
| :--------- | :-------------------------- | :--------------------------------- |
| pandas     | `>=2.0,<3`                  | `==3.0.0` (requires Python ≥ 3.11) |
| pyranges   | `>=0.0.120,<=0.1.4`         | `==0.1.4`                          |
| biopython  | latest                      | latest                             |
| pyyaml     | latest                      | latest                             |
| matplotlib | latest                      | latest                             |

CI runs both environments on Python 3.10 and 3.11 (and 3.13 for unit tests).
