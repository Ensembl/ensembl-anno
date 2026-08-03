# Long-read consensus preprocessing

`gmb-longread-consensus` collapses a raw per-read Minimap2 GTF into a
deduplicated consensus track suitable for use as `--minimap2` evidence in
`gmb-build`.

---

## Why preprocessing is required

Raw Minimap2 output contains one GTF record per aligned read. Feeding thousands
of overlapping single-read alignments into GMB directly would:

- Inflate cluster sizes (each read becomes a separate isoform candidate).
- Allow low-support singleton models to survive into the final geneset.
- Produce redundant identical structures that inflate isoform counts.

The consensus step groups reads by intron-chain identity (with configurable
boundary tolerance), requires a minimum number of independent supporting reads,
and emits one representative model per group.

---

## Config schema

All thresholds in `LongreadConsensusConfig` were derived from
GCA_000002765.3 (*P. falciparum*) and **are not organism-neutral**. Always
supply `--preset` or `--config` for your target organism.

```yaml
# Per-read rejection
max_span_length: 45000         # reads spanning > this bp are dropped as genomic/chimeric
max_intron_length: 3000        # reads with any intron > this bp are dropped

# Locus clustering
min_intergenic_gap: 500        # minimum gap between loci (bp)

# Consensus grouping
collapse_terminal_variation_bp: 100   # terminal boundary wobble tolerated before grouping
splice_site_tolerance_bp: 10          # intron boundary wobble for chain-signature matching

# Support thresholds
min_read_support_multi_exon: 2   # minimum independent reads for a multi-exon consensus
min_read_support_single_exon: 4  # minimum independent reads for a single-exon consensus

# Post-consensus
suppress_contained_models: true  # suppress models fully contained within a larger model

# Short-read rescue (disabled by default — see "Rescue modes" below)
shortread_rescue:
  enabled: false
```

---

## Bundled presets

Presets are installed as package data in `gmb.configs.longread_consensus`
and are available from any working directory after `pip install`.

| Preset | Description |
| :----- | :---------- |
| `pfalciparum_pure` | *P. falciparum* thresholds, no short-read rescue |
| `pfalciparum_assisted` | *P. falciparum* thresholds, short-read rescue enabled |

```bash
# List installed presets
gmb-longread-consensus --help
```

---

## Usage

**Minimum invocation (Stage 1 + Stage 2):**

```bash
gmb-longread-consensus \
    --input     raw_minimap2.gtf \
    --output-dir longread_out/ \
    --preset    pfalciparum_pure
```

**With a custom config file:**

```bash
gmb-longread-consensus \
    --input     raw_minimap2.gtf \
    --output-dir longread_out/ \
    --config    my_organism_lr.yaml
```

**Preset + per-run override (preset loaded first, file merged on top):**

```bash
gmb-longread-consensus \
    --input     raw_minimap2.gtf \
    --output-dir longread_out/ \
    --preset    pfalciparum_pure \
    --config    min_support_override.yaml
```

**Split-only (Stage 1, no consensus; useful on a cluster before array jobs):**

```bash
gmb-longread-consensus \
    --input     raw_minimap2.gtf \
    --output-dir longread_out/ \
    --split-only
```

`--preset` or `--config` is not required for `--split-only` (no thresholds
are applied during splitting). It is required for any Stage 2 consensus run.

**Slurm array pattern (process one chromosome per job):**

```bash
# Stage 1: split the raw file (run once)
gmb-longread-consensus \
    --input     raw_minimap2.gtf \
    --output-dir longread_out/ \
    --split-only

# Stage 2 per chromosome (one Slurm task per seqname)
gmb-longread-consensus \
    --input-split-gtf longread_out/_split_by_seqname/${SEQNAME}.gtf \
    --seqname         ${SEQNAME} \
    --output-dir      longread_out/ \
    --preset          pfalciparum_pure
```

**Reuse a validated split directory (resume after a partial run):**

```bash
gmb-longread-consensus \
    --reuse-split-dir longread_out/_split_by_seqname/ \
    --output-dir      longread_out/ \
    --preset          pfalciparum_pure
```

---

## Outputs

All outputs are written to `--output-dir`.

| File | Description |
| :--- | :---------- |
| `minimap2_consensus.gtf` | Consensus transcript models (source: `Minimap2Consensus`) |
| `longread_consensus_run_manifest.tsv` | Per-run config and provenance |
| `longread_consensus_summary.tsv` | Per-seqname read counts and filtering statistics |
| `longread_consensus_dropped_records.tsv` | Records dropped and the reason (per threshold) |
| `longread_consensus_rescue_attribution.tsv` | Per-model rescue outcome (when rescue enabled) |
| `_split_by_seqname/` | Per-seqname split files (Stage 1 output; reusable) |

Pass `minimap2_consensus.gtf` to `gmb-build` as:

```bash
gmb-build \
    ... \
    --minimap2 longread_out/minimap2_consensus.gtf
```

---

## Short-read rescue modes

When `shortread_rescue.enabled: true`, reads that fail the minimum-support
threshold can be rescued if short-read evidence corroborates their structure.

### Multi-exon rescue

| Policy | Rule |
| :----- | :--- |
| M1 (default) | Every intron junction confirmed by ≥1 short-read transcript (combined support across SR models allowed) |
| M2 | Complete ordered intron chain matches ≥1 single short-read transcript |
| M3 | Junctions independently confirmed by both Scallop AND StringTie |

`require_both_assemblers: true` also requires M3-level dual confirmation regardless of `policy`.

### Single-exon rescue

| Policy | Rule |
| :----- | :--- |
| S1 (default) | A single-exon SR transcript reciprocally overlaps the candidate by ≥`reciprocal_overlap_min` AND both boundaries are within `boundary_tolerance_bp` |
| S2 | Both Scallop AND StringTie independently satisfy S1 |
| S3 | ≥2 independent SR transcripts each satisfy S1 |

### Example rescue config

```yaml
shortread_rescue:
  enabled: true
  multi_exon:
    enabled: true
    policy: M1
    junction_tolerance_bp: 10
  single_exon:
    enabled: false
```

Supply short-read files at runtime:

```bash
gmb-longread-consensus \
    --input     raw_minimap2.gtf \
    --output-dir longread_out/ \
    --preset    pfalciparum_assisted \
    --scallop   scallop.gtf \
    --stringtie stringtie.gtf
```

---

## Config error handling

| Scenario | Error |
| :------- | :---- |
| Missing `--config` file | `Config error: Config file not found: /path` |
| Unknown config key | `Config error: Unknown longread_consensus config key: 'x'. Known keys: [...]` |
| Invalid rescue policy | `Config validation failed: Unknown multi_exon rescue policy 'MX'. Valid: ['M1', 'M2', 'M3']` |
| Unknown preset | `Config error: Preset 'x' not found. Available: ['pfalciparum_assisted', 'pfalciparum_pure']` |
| No `--preset` or `--config` for Stage 2 | `ERROR: --preset or --config is required for Stage 2 consensus.` |
