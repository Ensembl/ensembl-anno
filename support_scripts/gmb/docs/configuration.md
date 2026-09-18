# GMB configuration

**This is the starting point for anything configuration-related.**

| What you want to do | Where |
|---|---|
| Pick between `standard`, `fungi`, `apicomplexa` | [presets.md](presets.md) — the catalogue, with the evidence behind each |
| See what a preset actually resolves to | `gmb-preflight` prints it, or `build/resolved_config.yaml` |
| Override one value for a run | [Layering](#layering), below |
| **Add a new evidence source** (a tool GMB has not seen) | [Evidence roles](#evidence-roles--the-central-abstraction), below |
| Change how much a source is trusted | [Weights](#weights--keyed-by-role), below |
| Turn a biological policy on or off | [Biological policy](#biological-policy), below |
| Start a brand-new clade | `configs/new_clade_template.yaml` — copy it and fill in `CHANGE_ME` |
| **Create a validated clade preset** | [creating_a_preset.md](creating_a_preset.md) — the 8-step procedure |
| Set DIAMOND/Psauron paths | [Machine paths](#machine-paths), below |

## Layering

```
gmb/configs/standard.yaml      neutral base — ALWAYS loaded
        |
        v
gmb/configs/<preset>.yaml      --preset apicomplexa | fungi | standard
        |
        v
--config overlay1.yaml         run-specific, repeatable
--config overlay2.yaml         each later file wins on shared keys
        |
        v
build/resolved_config.yaml     what actually ran (+ resolved_config_sha256)
```

`--preset standard` (or `none`, or omitting it) loads the neutral base alone.

**Never edit a shipped preset.** They are validated baselines, and changing one invalidates
every past comparison. Put run-specific values in an overlay.

A useful convention is **two overlays**: one for machine paths, one for biological policy.
That lets the policy file be byte-identical across machines and species — which is itself
evidence that no per-species tuning crept in.

## Strict validation

- An **unknown key** raises `ValueError`. Typos fail at load, not silently.
- A **duplicate top-level key** raises `ValueError`. Plain YAML would let a second `scoring:`
  block replace the first outright — a silent, total loss of every setting in the earlier
  block. (This is not hypothetical: it happened during productionisation and quietly reverted
  two clade presets to neutral defaults while the build still succeeded.)
- A **source assigned to two evidence roles** raises `ValueError`.
- An invalid `backbone_intron_rescue` mode raises `ValueError`.
- Configurations that are merely *ineffective* warn rather than fail, so a run is never blocked
  by an evidence track simply being absent.

## Deprecated keys

Old keys still load, each with one `DeprecationWarning`:

| legacy key | current key |
|---|---|
| `helixer_filter` | `backbone_filter` |
| `scoring.keep_helixer_without_support` | `scoring.keep_backbone_without_support` |
| `scoring.weights.helixer` | `scoring.weights.backbone` |
| `scoring.weights.minimap2` | `scoring.weights.long_read` |
| `scoring.weights.scallop`, `.stringtie` | `scoring.weights.short_read` |

`scallop` and `stringtie` both collapse onto `short_read` because they are the same evidence
role. If they are set to the same value it is adopted; if they **disagree**, GMB warns loudly
and takes the maximum, so no evidence class is silently down-weighted. Setting `short_read`
explicitly always wins.

---

## Evidence roles — the central abstraction

Selection logic never tests a tool name. Every source label resolves to a **role**, and all
behaviour — ranking, retention, corroboration counting **and numeric weights** — acts on the
role.

| role | meaning |
|---|---|
| `backbone` | ab initio prediction; exactly one per run |
| `short_read_transcriptomic` | assembled short-read transcripts |
| `long_read_transcriptomic` | long-read transcript models |
| `protein_alignment` | protein-to-genome alignments — support/veto only, **never a candidate structure** |
| `protein_validation` | post-model scoring; not a track |

```yaml
scoring:
  backbone_label: Helixer                    # set by --helixer / --tiberius
  shortread_labels: [Scallop, StringTie]
  longread_label: Minimap2
  protein_alignment_labels: [OrthoDB, GenBlast, UniProt]
```

**Set these if your tools are named anything else.** An assembler called `IsoQuant` or
`AssemblerX` then behaves identically to `StringTie`. Preflight reports the resolved role and
weight for every track, and **fails** a source that lands in no role.

### Adding an evidence source GMB has never seen

Nothing in the selection logic tests a tool name. To add `AssemblerAlpha` as a short-read
source, list its label — that is the entire change:

```yaml
scoring:
  shortread_labels: [Scallop, StringTie, AssemblerAlpha]
```

It is now weighted, ranked, counted for structural corroboration and attributed exactly like
any other short-read track. The same applies to `backbone_label`, `longread_label` and
`protein_alignment_labels`.

The label must match what GMB assigns when it loads the file — which is fixed by the CLI flag
you pass it (`--scallop` → `Scallop`, `--minimap2` → `Minimap2`, and so on). The flag is only a
slot; the role is what decides behaviour.

**Check your work with preflight**: it prints the resolved role and weight for every track, and
**FAILs** any source that lands in no role rather than letting it fall through to the generic
`unknown` weight.

## Weights — keyed by role

```yaml
scoring:
  weights:
    backbone: 2.0
    short_read: 1.0
    long_read: 1.0
    protein_alignment: 1.0
    unknown: 1.0        # a source resolving to no role
```

Relative values matter; the absolute scale does not. Shipped: `standard` 2.0,
`apicomplexa` 2.6, `fungi` 3.1 for the backbone.

## Biological policy

All default **off**. Enable one at a time, for a reason visible in preflight output.

| setting | values | default |
|---|---|---|
| `structural_corroboration` | bool | `false` |
| `protein_support_mode` | `positional` \| `cds_span_compatible` | `positional` |
| `longread_structural_guard` | bool | `false` |
| `longread_disposition` | `primary_structural` \| `support_only` \| `reject` | `primary_structural` |
| `backbone_intron_rescue` | `off` \| `on` \| `auto` | `"off"` |

`backbone_intron_rescue` is **applicability-gated**:

- `off` — never fires.
- `on` — expert override; fires regardless of the measured evidence.
- `auto` — fires only where the backbone measurably under-resolves introns relative to
  credible assembled transcripts. Refuses whenever the evidence is missing, too thin or
  ambiguous.

The gate is reference-free and reports its reasoning in the log, the preflight report and the
run manifest. See `gmb/pipeline/applicability.py` for its calibration.

> **Why it is gated rather than on.** Rescued models were **68.6% CDS-exact** against a
> diatom-trained Tiberius backbone and **1.8%** against a Helixer backbone, where the rule
> destroyed models that were already correct. The same flag, opposite outcomes — the
> discriminator is the evidence state, not the clade.

## Correctness behaviour is not configurable

These are always on and have no switches, by design:

strand/end-aware UTR correction · exon reconstruction after trimming · FASTA regeneration from
the final GFF3 · strict sequence QC · exact duplicate transcript collapse · same-strand protein
support · protein evidence attribution · cross-seqid ID namespacing · gene-boundary
recomputation · gene/transcript containment validation · one canonical transcript per gene ·
run provenance.

A correctness property is never expressed as a tuning switch.

## Preflight thresholds

```yaml
preflight:
  backbone_splice_warn: 0.90
  backbone_splice_fail: 0.70
  shortread_splice_warn: 0.95
  shortread_splice_fail: 0.80
  longread_splice_warn: 0.95
  longread_splice_fail: 0.85      # strictest — this track claims to OBSERVE structure
  min_introns_for_splice_check: 100
  max_unknown_seqid_fraction: 0.05
  max_unstranded_fraction: 0.05
```

Roles are judged differently on purpose. A long-read track asserts it observed the splice
structure directly, so a poor canonical fraction there is a hard failure; a backbone is a
prediction and is judged more leniently; protein alignments are support-only and are not
splice-gated at all.

## Machine paths

Keep these in their own overlay, separate from biological policy:

```yaml
protein_validation:
  enabled: true
  diamond_path: "/path/to/diamond"
  psauron_path: "/path/to/psauron"
  diamond_db: "/path/to/proteins.dmnd"
```

GMB does not install DIAMOND or Psauron. Set `enabled: false` if unavailable — protein
validation only *scores* models that already exist.

## Starting a new clade

Copy `configs/new_clade_template.yaml`. An unedited copy resolves to exactly `standard`, so it
changes nothing until you edit it. See `creating_a_preset.md`.
