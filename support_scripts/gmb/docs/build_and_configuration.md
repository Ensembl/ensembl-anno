# GMB build configuration reference

> **Superseded for the production path — read [`configuration.md`](configuration.md) first.**
>
> This file remains the exhaustive per-section schema reference and is still accurate for
> the sections it covers, but it predates evidence-role weights, the
> `backbone_intron_rescue` applicability gate, `longread_disposition` and the preflight
> thresholds. Use it for depth, not as a starting point.

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
  min_codons: 33               # minimum ORF length (codons); standard 50, fungi 33
```

`allow_partial_5`, `allow_partial_3`, `allow_non_atg_start`, `stop_codon_char` and
`partial_prefix` are accepted but **not applied** (see `known_issues.md`).

### `protein_filter`

Controls filtering of protein-alignment evidence tracks (OrthoDB, UniProt,
GenBlast).

```yaml
protein_filter:
  min_alignment_coverage: 0.80 # minimum coverage FRACTION (needs a Coverage attribute)
  min_percent_identity: 60.0   # minimum identity in PERCENT (needs an Identity attribute)
  min_bitscore: 50.0           # minimum score (needs a numeric score on exon rows)
  min_protein_aa: 30           # alignments spanning < 3 x this many bp are dropped
  max_span_bp: 50000           # alignments spanning more are dropped as artefacts
  redundancy_overlap: 0.80     # reciprocal overlap at which alignments are collapsed
  top_n_per_locus: 3           # rank the N best per locus; with keep_secondary the rest stay
  keep_secondary: true
```

The three score thresholds only act when the alignment GTF carries those attributes on its
exon rows. Ensembl anno genBlastG output does not, so for it they remove nothing.

### `transcriptomic_filter`

Controls filtering of assembled transcript evidence (short- and long-read).

```yaml
transcriptomic_filter:
  max_intron_length: 3000      # drop transcripts with any intron > this bp (fungi 3000)
  max_transcript_length: 20000 # drop transcripts spanning > this bp; null = off (fungi 20000)
  allow_single_exon: true      # false drops single-exon transcripts
  min_intergenic_gap: 500
```

`max_transcript_length` and `allow_single_exon` are applied since 2.0.0 (earlier versions
accepted and ignored them). The build log and `summary.json` report how many transcripts
each rule removed (`chimeras_large_intron`, `chimeras_long_span`, `single_exon_removed`).
`strand_consistency_check` is deprecated and has no effect.

### `backbone_filter`

Controls filtering of *ab initio* backbone models (Helixer or Tiberius).
The `backbone` terminology is generic; the legacy keys `helixer_filter`,
`keep_helixer_without_support`, and `weights.helixer` are still accepted
with a deprecation warning.

```yaml
backbone_filter:
  enabled: true
  min_cds_bp: 90               # drop backbone models with less CDS than this
  max_exons: 50                # flag models with more exons than this
```

### `scoring`

Weights and thresholds for isoform scoring and selection.

```yaml
scoring:
  max_isoforms_per_locus: 3
  fungal_single_exon_mode: true        # single-exon handling; the name is historical
  keep_backbone_without_support: true  # keep backbone models with no other evidence
  locus_clustering: exon_overlap       # or transcript_linked -- see known_issues.md
  weights:                             # keyed by evidence ROLE, not tool name
    backbone: 3.1
    short_read: 1.0
    long_read: 1.0
    protein_alignment: 1.0
    unknown: 1.0
```

Weights are keyed by role (see `configuration.md`). The legacy tool keys `helixer`,
`scallop`, `stringtie` and `minimap2` are accepted with a deprecation warning; there are no
per-protein-track weights, because protein alignments are support, not candidates.

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

**Not applied.** This section configured the plots of the removed `gmb-visualize` tool. It is
still accepted so existing configs load, but nothing reads it.

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
