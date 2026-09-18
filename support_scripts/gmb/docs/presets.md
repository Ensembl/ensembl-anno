# GMB presets

A preset is a bundle of clade-appropriate settings layered on the neutral base. Pick one with
`--preset`; override anything with `--config`.

```
standard.yaml  ->  <preset>.yaml  ->  --config overlays
```

**Never edit a shipped preset.** They are validated baselines; changing one silently
invalidates every past comparison. Put run-specific values in an overlay.

---

## Which preset should I use?

| preset | intended use | backbone expectation | short-read expectation | long-read policy | `backbone_intron_rescue` |
|---|---|---|---|---|---|
| **`standard`** | a new or unknown clade | unknown | unknown | quality-gated by preflight; no guard by default | **`off`** |
| **`apicomplexa`** | Apicomplexa with a general-purpose ab initio backbone | **may be intron-collapsed** | a genuinely better structural source | guard on; quality-gated | **`auto`** |
| **`fungi`** | fungi with a strong Helixer-like backbone | **strong** | supportive, *not* automatically superior | none assumed; guard off (inert) | **`off`** |

---

## `standard` — the neutral base

**Use when you have not validated an evidence state for your clade.** Also reachable as
`--preset none` or by omitting `--preset`; `standard` is the explicit spelling.

| setting | value |
|---|---|
| `weights.backbone` | 2.0 |
| `weights.short_read` / `long_read` / `protein_alignment` / `unknown` | 1.0 |
| `multi_source_bonus` | 1.0 |
| `structural_corroboration` | `false` |
| `protein_support_mode` | `positional` |
| `longread_structural_guard` | `false` |
| `longread_disposition` | `primary_structural` |
| `backbone_intron_rescue` | `"off"` |

It is **behaviour-neutral**: no optional policy is enabled, and no clade assumption is
embedded. It deliberately assumes neither that assemblies beat the backbone nor the reverse —
which of those is true is exactly what differs between the two validated clades, and preflight
measures it for you.

## `apicomplexa` — validated on P. falciparum and T. gondii

| setting | value | why |
|---|---|---|
| `backbone_label` | `Tiberius` | no Helixer model for these genomes |
| `weights.backbone` | **2.6** | lowered from the neutral 2.0-vs-fungi-3.1 scale: this backbone's multi-exon splice placement is its clear weak point |
| `weights.long_read` | 1.3 | long reads span introns directly *when the track is sound* |
| `multi_source_bonus` | 1.2 | raised: agreement with independent transcript evidence is where backbone errors get corrected |
| `structural_corroboration` | `true` | rank on independent agreement before the numeric score |
| `protein_support_mode` | `cds_span_compatible` | only CDS-span-compatible alignments gate retention |
| `longread_structural_guard` | `true` | fired 332× on P. falciparum; 0× with no long-read track |
| `backbone_intron_rescue` | **`auto`** | applicability-gated, not unconditional |

**Evidence.** P. falciparum GCA_000002765.3, baseline → this policy, no long-read track:

```
CDS exact              2,128 -> 2,430
CDS exact, multi-exon    570 ->   957   (+68%)
Exact Match              447 ->   518
canonical intron frac   87.8% -> 100.0%   (reference 99.94%)
improvements : regressions = 731 : 447
```

Validated independently on T. gondii GCA_000006565.2 with the identical policy file.

**What this preset actually assumes** is an *evidence state*, not a clade: a backbone that
under-calls introns (17.3% and 38.7% multi-exon in the two genomes) alongside assemblies that
do not (95–99%). If a future Apicomplexan arrives with a well-resolved backbone, the `auto`
gate declines and the preset degrades safely.

## `fungi` — validated on Z. tritici

| setting | value | why |
|---|---|---|
| `backbone_label` | `Helixer` | inherited from `standard` |
| `weights.backbone` | **3.1** | Helixer is well-trained on fungi and is the dominant signal |
| `multi_source_bonus` | 0.5 | lowered: agreement is a modest rather than strong signal here |
| `max_intron_length` | 3,000 | fungal introns are short (reference median 62 bp) |
| `max_transcript_length` | 20,000 | compact transcripts; caps chimeras |
| `structural_corroboration` | **`false`** | validated baseline |
| `protein_support_mode` | **`positional`** | validated baseline |
| `longread_structural_guard` | **`false`** | no long-read evidence existed; provably inert |
| `backbone_intron_rescue` | **`"off"`** | see below |

### Why fungi does NOT inherit the Apicomplexa policy

The Apicomplexa policy was applied unchanged to the full Z. tritici genome and made the
annotation **measurably worse**:

```
CDS exact              4,895 -> 4,611   (-284)
CDS exact, multi-exon  2,984 -> 2,877   (-107)
Exact Match            1,881 -> 1,727   (-154)
locus detection       97.24% -> 94.91%
improvements : regressions = 249 : 555   (ratio 0.449)
```

The mechanism is visible in the inputs, with no reference needed:

| | Tiberius (P. falciparum) | **Helixer (Z. tritici)** | Z. tritici reference |
|---|---|---|---|
| backbone multi-exon | 17.3% | **71.9%** | **69.8%** |

The fungal backbone already matches the reference's exon-count distribution. There is nothing
to rescue, so every intervention is a chance to damage a model that was already right — and it
did: rescued models were **1.8% CDS-exact** against 31.9% for backbone-only.

CDS-exact by evidence bucket **inverts** between the two clades:

| bucket | P. falciparum | Z. tritici |
|---|---|---|
| backbone + short-read agreement | 76.7% | 55.5% |
| backbone only | 47.1% | 30.1% |
| **short-read only** | **60.2%** | **5.2%** |

In Apicomplexa the assemblies beat the backbone; in fungi they are six times worse than it.
A policy that promotes assemblies over the backbone must help in one and harm in the other.

### Two fungal settings that are candidates for a future experiment

Both are **off** today because only the baseline has been validated end-to-end, but both look
promising in isolation and are worth a per-flag experiment:

- **`structural_corroboration`** — multi-source models are 42.7% CDS-exact against 28.1% for
  single-source. It enriches more strongly in fungi than in either apicomplexan.
- **`protein_support_mode: cds_span_compatible`** — protein support is **strictly monotone** in
  fungi (strong 38.0% / positional-only 6.8% / none 0.8%), which is exactly what the policy
  assumes. It is *not* monotone in P. falciparum, where "none" beats "positional-only".

### Scope of the fungal validation

One genome. `fungi.yaml`'s own header records that it was tuned against Z. tritici, so the
**preset** is not independent of this organism — only the *selection policy* was tested
independently. Do not read "validated on Z. tritici" as "validated for fungi".

---

## Building a preset for a new clade

See `creating_a_preset.md` and start from `configs/new_clade_template.yaml`.
