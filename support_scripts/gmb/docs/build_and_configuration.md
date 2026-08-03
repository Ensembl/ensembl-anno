# GMB build configuration reference

Configuration is assembled in layers. This document describes the layering
model, every top-level config section, and how to supply environment-specific
paths without duplicating biological settings.

---

## Configuration layers and precedence

GMB builds the final configuration in four layers, applied in order:

```text
1. standard.yaml        — organism-neutral base (always loaded, bundled)
        ↓
2. clade preset         — species-class overrides (--preset fungi | apicomplexa | ...)
        ↓
3. --config file(s)     — project- or run-specific overrides (user-supplied, in order)
```

Each later layer deep-merges its dict keys on top of the previous state;
list-valued keys are replaced entirely (never concatenated).

### Bundled presets

| Preset | Description |
| :----- | :---------- |
| `fungi` | Compact-genome ascomycete / basidiomycete defaults (default) |
| `apicomplexa` | Apicomplexa defaults (derived from *P. falciparum* GCA_000002765.3) |

```bash
# List installed presets
gmb-build --list-presets

# Use the apicomplexa preset
gmb-build --preset apicomplexa ...

# Standard base only, no clade preset
gmb-build --preset none ...
```

### Typical usage patterns

**Single-species project with organism tuning:**

```bash
gmb-build \
  --preset apicomplexa \
  --config local_cluster_paths.yaml \
  ...
```

Effective configuration:

```text
standard.yaml
    → apicomplexa.yaml (backbone: Tiberius, intron/UTR caps, etc.)
        → local_cluster_paths.yaml (diamond_db, diamond_path, etc.)
```

**Multiple config files — later file wins on shared keys:**

```bash
gmb-build \
  --config my_species.yaml \
  --config protein_validation_overlay.yaml
```

```text
standard.yaml
    → fungi.yaml
        → my_species.yaml
            → protein_validation_overlay.yaml
```

### Splitting biological settings from environment paths

Keep organism-specific biology in one file and site-specific paths in another.
The path file can be kept out of version control:

```yaml
# apicomplexa.yaml  (versioned, shared)
orf:
  min_codons: 33

protein_validation:
  enabled: true
```

```yaml
# local_cluster_paths.yaml  (local only, not committed)
protein_validation:
  diamond_path: /hps/software/diamond
  diamond_db: /hps/nobackup/team/apicomplexa.dmnd
```

```bash
gmb-build \
  --preset apicomplexa \
  --config local_cluster_paths.yaml \
  ...
```

### Merge rules

- **Dicts** deep-merge: keys absent in the overlay are preserved from the layer below.
- **Lists** replace entirely: supplying a list key discards the earlier list.
- **Unknown top-level keys** raise `ValueError` immediately (catches typos).
- **Deprecated keys** (`helixer_filter`, `keep_helixer_without_support`,
  `weights.helixer`) emit a `DeprecationWarning` and are aliased to their
  generic equivalents.
- **Deprecated preset names** (`fungi_default`) emit a `DeprecationWarning`
  and are silently remapped to their current names.

A missing or misspelled `--config` path raises `FileNotFoundError`
immediately — GMB never silently falls back.

### Resolved configuration output

Every run writes two files to `--output-dir`:

| File | Description |
| :--- | :---------- |
| `resolved_config.yaml` | Complete effective configuration (all layers merged) |
| `resolved_config_sha256` | SHA-256 hex digest of the above for integrity checks |

Use `resolved_config.yaml` to reproduce a run exactly by passing it as
`--config resolved_config.yaml --preset none`.


---

## Top-level sections

### `orf`

Controls open-reading-frame detection.

```yaml
orf:
  min_codons: 33               # minimum ORF length (codons)
  allow_partial_5: true        # allow ORFs without a start codon
  allow_partial_3: false       # allow ORFs without a stop codon
```

### `protein_filter`

Controls filtering of protein-alignment evidence tracks (OrthoDB, UniProt,
GenBlast).

```yaml
protein_filter:
  min_alignment_coverage: 0.5  # minimum query coverage fraction
  min_percent_identity: 0.3    # minimum percent identity
  min_bitscore: 50             # minimum bitscore
  min_protein_aa: 30           # minimum protein length (amino acids)
  top_n_per_locus: 3           # keep at most N protein models per locus
```

### `transcriptomic_filter`

Controls filtering of RNA-seq assembly evidence (Scallop, StringTie).

```yaml
transcriptomic_filter:
  max_intron_length: 100000    # drop transcripts with any intron > this bp
  max_transcript_length: null  # drop transcripts longer than this bp (null = no limit)
  allow_single_exon: true      # keep single-exon transcriptomic models
  min_exon_length: 10          # drop exons shorter than this bp
```

### `backbone_filter`

Controls filtering of *ab initio* backbone models (Helixer or Tiberius).
The `backbone` terminology is generic; the legacy keys `helixer_filter`,
`keep_helixer_without_support`, and `weights.helixer` are still accepted
with a deprecation warning.

```yaml
backbone_filter:
  min_exon_count: 1
  max_intron_length: 500000
  require_protein_support: false
  require_transcriptomic_support: false
```

### `scoring`

Weights and thresholds for isoform scoring and selection.

```yaml
scoring:
  max_isoforms_per_locus: 3
  fungal_single_exon_mode: true     # fungi-specific single-exon handling
  keep_backbone_without_support: true  # keep backbone models with no other evidence
  weights:
    backbone: 2.0
    scallop: 1.0
    stringtie: 1.0
    minimap2: 1.0
    orthodb: 1.5
    uniprot: 1.5
    genblast: 1.0
```

### `protein_validation`

Optional DIAMOND + Psauron protein-coding validation. Disabled by default.

```yaml
protein_validation:
  enabled: false
  diamond_path: diamond         # path to diamond executable
  psauron_path: psauron         # path to psauron executable
  diamond_db: null              # REQUIRED when enabled=true and diamond_weight > 0
                                # set to the absolute path of your .dmnd database
  diamond_weight: 0.7
  psauron_weight: 0.3
  min_score: 0.7
  policy: penalize              # "drop" | "penalize" (also accepts "penalise")
  psauron_min_length: 5         # psauron -m/--minimum-length (aa)
  psauron_use_cpu: false        # psauron -c/--use-cpu
  diamond_min_query_coverage: 0.0   # additional DIAMOND hit gate (0-100)
  diamond_min_target_coverage: 0.0  # additional DIAMOND hit gate (0-100)
```

`diamond_db` has no default. Setting `enabled: true` with `diamond_weight > 0`
and no `diamond_db` raises a `ValueError` at config load time — it never silently
runs without a database.

### `utr`

Controls UTR retention and end-support validation.

```yaml
utr:
  require_end_support: true
  end_support_mode: multisource_end_agreement  # "multisource_end_agreement" | "protein_validated"
  end_support_sources:                         # assemblers that must agree on UTR boundaries
    - Scallop
    - StringTie
  end_tolerance_bp: 50
  require_multisource_for_utr_5p: true
  require_multisource_for_utr_3p: true
  fallback_policy_when_unsupported: drop_utr   # "drop_utr" | "hard_cap" | "drop_transcript"
  min_protein_coding_score_for_utr: null       # minimum score to keep UTRs (null = no gate)
  max_end_extension_bp: null                   # cap UTR extension (null = no limit)
```

### `qc`

Per-track quality-control thresholds applied before scoring.

```yaml
qc:
  max_transcripts_per_track: 5
  skip_orf_inference_tracks:
    - OrthoDB
    - UniProt
  parallel: false
  workers: 4
```

### `duplicate_transcript_collapse`

Controls exact-duplicate transcript collapsing within a gene.

```yaml
duplicate_transcript_collapse:
  collapse_exact_duplicates: true
```

### `transcript_splitting`

Controls splitting of pathological "mega-transcripts" produced by some
assemblers.

```yaml
transcript_splitting:
  split_enabled: false          # disabled in fungi_default.yaml
  split_gap_bp: 3000
  split_on_large_exon_bp: 15000
  max_segments_per_transcript: 50
```

### `canonical_selection`

Controls post-build canonical transcript selection and the InterPro
resolver second stage.

```yaml
canonical_selection:
  interpro_resolver:
    enabled: false
    # ... (see docs/interpro_resolver.md for the full schema)
```

### `validation`

Controls GFF3 structural validation and repair.

```yaml
validation:
  mode: drop_transcript         # "error" | "fix" | "drop_transcript"
  log_violations: true
  max_feature_drift_bp: 1500
  feature_outside_exons_policy: trim
  max_exon_len_bp: 15000
  max_exon_len_mode: fixed      # "fixed" | "percentile"
  max_exon_len_percentile: 99.5
```

---

## Config merge rules

- **Dicts** deep-merge: keys present in the override are updated; keys absent
  in the override are preserved from the preset.
- **Lists** replace entirely: if an override sets a list key, the preset's
  list is discarded.
- **Unknown top-level keys** raise `ValueError` immediately (catches typos).
- **Deprecated keys** (`helixer_filter`, `keep_helixer_without_support`,
  `weights.helixer`, `scoring.keep_helixer_without_support`,
  `scoring.weights.helixer`) emit a `DeprecationWarning` and are aliased to
  their generic equivalents.

---

## Checking for external tools

Before a run that uses protein validation or long-read consensus:

```bash
gmb-build --check-deps
```

This detects DIAMOND and Psauron availability, reports the detected Psauron
version, and fails clearly if a required flag is missing from the installed binary.
