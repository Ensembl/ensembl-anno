---
name: gmb-new-clade-config
description: Derive and justify a Gene Model Builder (GMB) configuration for a taxonomic group that has no validated preset, from that group's own annotation evidence. Use when someone wants to create a GMB config or clade preset for a new clade (plants, oomycetes, nematodes, protists, ...), to assess whether the existing standard/fungi/apicomplexa settings suit new evidence, to compare the quality of backbone vs transcript vs protein evidence tracks, or to produce a provisional or validated clade preset. Do NOT use to simply run GMB with an existing preset — that is `gmb-preflight` → `gmb-build` → `gmb-finalise` as documented in support_scripts/gmb/docs/quickstart.md.
---

# GMB new-clade configuration

Take the evidence files for a representative genome of a new clade, measure them, and write a
GMB configuration whose every difference from the neutral `standard` preset is justified by a
measurement. Then run it, QC it, and — only if a reference is supplied, and only after the
configuration is frozen — evaluate it.

The question this skill answers is: **"Given these evidence files, what should my GMB
configuration look like, and why?"** The answer is always a delta from `standard`, never a copy
of `fungi.yaml` or `apicomplexa.yaml`.

Paths below are relative to the GMB root, `support_scripts/gmb/`. `$SKILL` is this directory.

---

## Guardrails — non-negotiable

- **Never use a reference annotation in model selection.** It is never passed to
  `gmb-preflight`, `gmb-build` or `gmb-finalise`, and never listed as evidence.
- **Never dynamically tune production configuration from reference metrics.** No code, loop or
  script may read comparison results and write config values.
- **Never choose a policy purely from the organism or clade name.** "This is an apicomplexan"
  or "this is a fungus" is not a reason. A measurement is.
- **Never assume long-read evidence is structurally reliable.** Measure its splice quality first.
- **Never assume transcript evidence is superior to the backbone** — or the reverse. Measure.
- **Never assume protein support proves exon structure.** A protein hit says a locus is coding;
  it does not say the introns are right.
- **Never change a parameter without recording why**, in `config_decisions.tsv`.
- **Never call one-genome tuning "clade validation".**
- **Never edit the shipped presets** (`gmb/configs/standard.yaml`, `fungi.yaml`,
  `apicomplexa.yaml`) or any existing clade config unless the user explicitly asks.
- **Never commit or push** unless the user explicitly asks.

---

## Facts about GMB you must respect

Verify these against the current code if anything below looks stale (`gmb-build --help`,
`gmb/pipeline/config.py`); do not invent keys.

**The production path is three commands:** `gmb-preflight` → `gmb-build` → `gmb-finalise`.
`gmb-compare` is evaluation only.

**Configuration layers:** `gmb/configs/standard.yaml` → optional `--preset` → each `--config`
overlay in order, last wins. Unknown keys and duplicate YAML keys are **errors**, so loading a
config through `gmb-preflight` is also a schema validation. A new clade config is an **overlay**
used as `--preset standard --config configs/<clade>.yaml`. It lives in `configs/` (next to
`new_clade_template.yaml`), **not** in `gmb/configs/` — that package directory is for promoted,
clade-validated presets only.

**Evidence roles are the abstraction. Tool names are not.** GMB has fixed input slots:

| role | slots (in order) | cap |
|---|---|---|
| `backbone` | `--helixer` or `--tiberius` | exactly 1 |
| `short_read_transcriptomic` | `--scallop`, `--stringtie` | 2 |
| `long_read_transcriptomic` | `--minimap2` | 1 |
| `protein_alignment` | `--orthodb`, `--uniprot`, `--genblast` | 3 |

- A slot is just a slot. `AssemblerA.gtf` goes in `--scallop`; `PredictorFoo.gff3` goes in
  `--helixer`. The slot fixes GMB's **internal label** (e.g. "Scallop"); the built-in role
  labels already map every slot to the right role, so **do not edit `backbone_label`,
  `shortread_labels`, `longread_label` or `protein_alignment_labels`** — `gmb-build` even forces
  `backbone_label` from the slot used.
- Because internal labels are slot names, **every report must carry the slot map** (real tool ↔
  slot ↔ internal label) so a reader never mistakes "Scallop" in a GMB output for the real tool.
- **Parsing is by file extension**: a name ending `.gtf` is read as GTF, anything else as GFF3.
  A gzipped GTF (`*.gtf.gz`) would be misread. Supply evidence uncompressed.
- More tracks than slots cannot be passed. Merge upstream or drop one, and record why.

`scripts/gmb_clade.py args` enforces all of this.

**Things that are correctness, not configuration** — always on, never touched by this skill: UTR
end/strand handling, exon reconstruction, FASTA regeneration from the final GFF3, strict
sequence QC, duplicate collapse, same-strand protein support, cross-seqid ID namespacing,
gene-boundary recomputation, one canonical per gene, provenance.

**`gmb-finalise` uses the build's own `resolved_config.yaml`.** Do not pass `--preset`/`--config`
to it.

---

## Inputs and interaction

Collect a manifest (see `example_manifest.yaml`):

```yaml
clade: nematodes                 # the config to create: configs/nematodes.yaml
genome: genome.fa
machine_overlay: local_paths.yaml   # optional: tool/database paths ONLY (allowlist enforced)
evidence:
  - {path: predictor.gff3,  tool: PredictorFoo,  role: backbone}
  - {path: assembler_a.gtf, tool: AssemblerA,    role: short_read_transcriptomic, rnaseq_source: runs1}
  - {path: assembler_b.gtf, tool: AssemblerB,    role: short_read_transcriptomic, rnaseq_source: runs1}
  - {path: long_reads.gtf,  tool: LongReadsX,    role: long_read_transcriptomic}
  - {path: orthodb.gtf,     tool: ProteinMapperY, role: protein_alignment}
evaluation_only:                 # NEVER passed to preflight/build/finalise
  reference: reference.gff3      # evaluation only; must currently be GFF3
  seqname_map: map.tsv
```

`rnaseq_source` is optional: give two short-read tracks the same value when they were assembled
from the same RNA-seq, so their agreement is not mistaken for independent evidence.

Ask only for what is genuinely missing: clade name, genome, which file is the backbone, which
are short-read / long-read / protein, and whether a reference or a second genome exists. If a
filename makes a role obvious, **state the inference and proceed**; if a role is ambiguous,
**ask** — never guess silently. If the clade name collides with an existing config
(`configs/<clade>.yaml` or `gmb/configs/<clade>.yaml`), ask for a different name.

---

## Workflow

Set `OUT` to a fresh directory per attempt. Builds on large or fragmented genomes can take a day
or more (a 2,263-sequence genome took ~48 h); run them as dedicated jobs. A `--seqname` subset is
fine for **debugging** but never sufficient for a validation status.

### Phase 1 — inventory and validate inputs

```bash
python $SKILL/scripts/gmb_clade.py args manifest.yaml --out $OUT     # slot map + checks
EV=$(python $SKILL/scripts/gmb_clade.py args manifest.yaml --shell)  # CLI argument string
```

This refuses: missing files, unknown roles, a reference listed as evidence, a machine overlay
containing anything but tool/database paths or hardware settings, GTF/GFF3 content
that contradicts the extension, `.gtf.gz`, more tracks than slots, zero or two backbones. It
writes `slot_map.tsv`. Everything else — seqids, strand, collisions, spans, splice quality — is
preflight's job; do not re-implement it.

### Phase 2 — run GMB preflight (twice)

```bash
eval gmb-preflight --preset standard $EV --output-dir $OUT/preflight_standard --allow-fail
eval gmb-preflight --preset standard --config $SKILL/probe_rescue_auto.yaml $EV \
     --output-dir $OUT/preflight_probe --allow-fail
```

The first is the neutral report. The second is a **probe**: it asks the *implemented*
reference-free `backbone_intron_rescue` gate what it would decide for this bundle (under
`standard` the gate is `off` and reports nothing). The probe is never used for a build.

Fix any FAIL in `seqid_compatibility`, `role_resolved`, `file_readable` or `backbone_present`
before continuing. A FAIL in `splice_quality` / `longread_policy` for a long-read track is a
finding, not a blocker — it drives the long-read decision.

### Phase 2b — neutral baseline build (reference-free measurements + the comparator)

```bash
eval gmb-build --preset standard $EV --gene-prefix XXGMB \
     --output-dir $OUT/standard/build --validate-fasta
gmb-finalise --build-dir $OUT/standard/build --genome <genome.fa> --output-dir $OUT/standard/finalise
```

This gives the measurements preflight cannot: backbone ↔ assembly agreement, protein-support
breadth and redundancy, UTR end-support behaviour. It is also the `standard` arm every later
comparison is made against.

### Phase 3 — measure and classify the evidence state

```bash
python $SKILL/scripts/gmb_clade.py measure manifest.yaml \
    --preflight $OUT/preflight_standard --probe $OUT/preflight_probe \
    --baseline-build $OUT/standard/build --out $OUT
```

Read `$OUT/evidence_summary.md`. **Describe the evidence state in words before touching any
setting.** Flags are computed from measurements, never from the organism:

| flag | evidence state | raised when |
|---|---|---|
| **A** | strong backbone relative to assemblies | the implemented rescue gate declines because the backbone is well resolved / not clearly under-resolved |
| **B** | intron-collapsed backbone, strong spliced transcripts | the implemented rescue gate **fires** |
| **C** | strong backbone ↔ transcript agreement | judgement from the baseline's *identical-intron-chain agreement*: materially above the ~1% seen with collapsed backbones (11.4% was seen with a well-resolved one). Descriptive, not a threshold |
| **D** | poor transcript evidence | no short-read track, or a short-read track WARN/FAIL on splice quality |
| **E** | structurally unreliable long reads | long-read splice quality FAIL (`E?` for WARN) |
| **F** | sparse evidence | preflight WARNs on backbone/transcript/protein presence, or too few models for the gate |

Flags combine (B + E is common). Also read, per structural track: model count, multi-exon
fraction, canonical splice fraction, intron and span percentiles; and for protein tracks: input
volume, redundancy collapse, and the share of candidates with support.

### Phase 4 — start from standard

Copy `configs/new_clade_template.yaml` only as a **reference for what can be set**. Write
`configs/<clade>.yaml` containing **only** the settings that differ from `standard`, each with a
one-line comment giving the measurement. An empty overlay is a legitimate result.

Never start from `fungi.yaml` or `apicomplexa.yaml`. Those encode two *opposite* evidence states.

### Phase 5 — decide, setting by setting

Default for every setting: **keep `standard`**. Change it only when a rule below fires, and
write a row in `config_decisions.tsv` either way for each area considered:

```
setting	standard_value	proposed_value	evidence	reason	confidence	validation_required
```

Split settings into two kinds, and treat them differently:

- **Genome-architecture settings** (intron/transcript caps, splitting, ORF length, UTR caps)
  describe the genome. They can be derived reference-free, and `standard` is deliberately
  permissive for them (a 1,000,000 bp intron cap is effectively no chimera filter), so a new
  clade usually *should* change them.
- **Selection-policy settings** (rescue, corroboration, protein mode, weights, isoforms) change
  which model wins. Their benefit **cannot be measured without a reference**. In a provisional
  config, change one only where a reference-free rule below explicitly supports it; otherwise
  record it as a hypothesis for Phase 7.

#### Genome architecture

**`transcriptomic_filter.max_intron_length`** (standard 1,000,000)
The backbone is the track least prone to chimeric joins, so use its intron distribution as the
anchor — with generous headroom, because a backbone under-samples the longest real introns.
Starting rule: round **~3 × the backbone intron p99.9** up to a round number, then check it
against the short-read tracks *after discarding their chimeric tail*. A chimeric tail shows up
as short-read p95/p99 values tens of kb long while the backbone's maximum is ~1 kb. If the
short-read introns and backbone agree in scale, the short-read p99.9 may raise the cap.
Worked examples from two real bundles: backbone p99.9 of 904 bp and 725 bp gave caps of
~2,500–3,000 bp; in the first, short-read p95 was already 21–36 kb (chimeric) and was rightly
ignored. Too tight deletes real genes; too loose admits chimeras — err generous. Clades with
long introns get a large cap automatically, because their backbone p99.9 is large.
Confidence: medium.

**`transcript_splitting.split_enabled` / `split_gap_bp`** (standard false / 3000)
Enable when assembled transcripts are much longer than backbone models — e.g. short-read span
p99 of 115–145 kb against backbone span p99 of 6–10 kb, which is overjoining. Set
`split_gap_bp` equal to the chosen `max_intron_length`. Confidence: medium-high.

**`transcriptomic_filter.max_transcript_length`** (standard 2,000,000)
~1.5–2 × the backbone's maximum span, rounded up (worked examples: backbone max 30 kb → ~35–45
kb; 19 kb → ~20–30 kb). Confidence: medium.

**`orf.min_codons`** (standard 50)
Lower to 33 only if a non-trivial share (≳1%) of *backbone* models have CDS of 33–49 codons —
the predictor itself calls short genes, and 50 would discard them. If that share is ~0, keep 50.
(Real bundles: 1.3% → evidence for 33; 0.0% → keep 50.) Confidence: medium-low.

**UTR caps: `utr.max_5p_bp`, `max_3p_bp`, `max_total_bp`, `max_utr_to_cds_ratio`**
Keep `standard`'s end-support rules unchanged (`require_end_support: true`,
`multisource_end_agreement`, `fallback_policy_when_unsupported: drop_utr`) — they are safety
rules, and UTR handling was a major source of historical error. Consider tightening a cap only
when the baseline's **end-supported (kept)** UTR p99 is far below it; never set a cap below that
p99, and never tune a UTR length against a reference. Note:
- with **one** short-read track, multi-source end agreement cannot be satisfied, so most UTRs
  are dropped — the safe outcome; record it;
- two assemblers of the **same** RNA-seq agree on ends partly by construction — say so.

#### Backbone and transcripts

**`scoring.backbone_intron_rescue`** (standard `"off"`)
Use the probe's verdict, which *is* the implemented gate:
- gate **fires** (flag B) → `auto`. Reason: the backbone measurably under-resolves introns
  relative to credible spliced assemblies. Use `auto`, never `on`, so another genome of the
  clade with a better backbone is protected.
- gate **declines** or is not evaluable → keep `"off"`. If it declined narrowly or for
  insufficient evidence, say so under uncertainties.
Never set `on`. Never decide this from the clade name.

**`scoring.weights.backbone` / `short_read` / `long_read`** (standard 2.0 / 1.0 / 1.0) and
**`scoring.multi_source_bonus`** (1.0)
Weight magnitudes have **no reference-free calibration** — keep `standard` in a provisional
config. If flags A + C suggest the backbone is the better structural authority, record
"raise backbone weight" as a *hypothesis*. If every short-read track shares one `rnaseq_source`,
do **not** raise `multi_source_bonus`: their agreement is not independent.

**`scoring.structural_corroboration`**, **`scoring.protein_support_mode`**,
**`scoring.require_protein_support_for_single_source`**, **`scoring.max_isoforms_per_locus`**,
**`scoring.fungal_single_exon_mode`** (despite its name, generic single-exon handling)
Keep `standard` in a provisional config. Their effect was only ever measured with a reference,
and it **reversed between clades**: the same policy bundle helped a weak-backbone genome
(731 : 447 improved : regressed) and harmed a strong-backbone one (249 : 555). In state B they
are reasonable hypotheses; in state A treat them sceptically.

#### Long reads

Absent → change nothing. The guard cannot fire and nothing special happens.

Present → read its preflight `splice_quality` verdict:

| verdict | `scoring.longread_disposition` | `scoring.longread_structural_guard` |
|---|---|---|
| PASS | `primary_structural` (standard) | standard (`false`), or `true` as a cautious option — justify |
| WARN | `primary_structural` | `true` — demote long-read-only structures where a multi-exon alternative exists |
| FAIL | `support_only` if it plausibly adds loci; `reject` if it also has other problems | `true` |

Explain the choice with the measured canonical fraction. A long-read track once measured
**16.4%** canonical against 99–100% for every other track and, locus by locus, worsened as many
models as it improved while adding coverage — `support_only` fits that profile. Do not
special-case any sequencing technology or aligner.

#### Proteins

Protein alignments are support and veto evidence, **never** candidate structures. From the
baseline, report: input alignments → after redundancy collapse (heavy collapse is normal for
multi-species databases), and the share of selected models with strong / weak / no support.
Low breadth weakens every protein-gated rule; say so. Keep `protein_support_mode: positional`
unless deliberately testing the hypothesis above.

**`protein_validation`** (DIAMOND + Psauron) is distinct from protein-to-genome evidence. If the
tools exist, set `protein_validation.enabled: true` **in the clade config** (it changes scoring,
so it needs a recorded reason) and put only the tool and database paths in the machine overlay.
Standard `policy: drop` deletes models that
score below `min_score`; if the DIAMOND database covers the clade poorly, prefer
`policy: penalize` — reason: poor database coverage, not the clade name.

### Phase 6 — provisional config, build, hard QC

```bash
# schema check + preflight under the proposed config
eval gmb-preflight --preset standard --config configs/<clade>.yaml $EV \
     --output-dir $OUT/preflight_candidate --allow-fail
# build + finalise (whole genome for any validation status)
eval gmb-build --preset standard --config configs/<clade>.yaml $EV --gene-prefix XXGMB \
     --output-dir $OUT/candidate/build --validate-fasta
gmb-finalise --build-dir $OUT/candidate/build --genome <genome.fa> --output-dir $OUT/candidate/finalise
cp $OUT/candidate/build/resolved_config.yaml $OUT/resolved_config.yaml
# reference-free hard QC + diagnostics, for both arms
python $SKILL/scripts/gmb_clade.py summarise-run --build $OUT/standard/build \
    --finalise $OUT/standard/finalise --genome <genome.fa> --label standard --out $OUT
python $SKILL/scripts/gmb_clade.py summarise-run --build $OUT/candidate/build \
    --finalise $OUT/candidate/finalise --genome <genome.fa> --label candidate --out $OUT
```

All hard QC must pass: 0 cDNA/CDS/protein mismatches, 0 internal stops, 0 gene-boundary
violations, 0 UTR-invariant violations, one canonical per gene, no cross-seqid or
opposite-strand parent/child. Any failure: stop, report it, do not proceed to evaluation.

Also compare the arms reference-free: gene count, single/multi-CDS-exon share, canonical intron
fraction of the annotation, very long genes/introns. **A QC PASS is not a quality verdict.**

### Freeze — before any reference is touched

This is the guarantee the whole workflow rests on. The order is fixed:

```
build the candidate (Phase 6)
    ↓
obtain $OUT/candidate/build/resolved_config.yaml
    ↓
freeze the source overlay + the resolved config
    ↓
only then open, inspect or use the reference
    ↓
evaluate (Phase 7)
```

```bash
python $SKILL/scripts/gmb_clade.py freeze configs/<clade>.yaml \
    --resolved-config $OUT/candidate/build/resolved_config.yaml --out $OUT
```

**The frozen resolved config is the authoritative record of what GMB actually ran with** —
`standard` + clade overlay + machine overlay, as resolved by `gmb-build`. The source overlay is
frozen too, as the record of what the developer wrote. `freeze` refuses a resolved config whose
build's `run_manifest.json` does not list this overlay, or whose hash differs from the one that
build recorded, or that lacks any value the overlay sets.

`evaluate` then verifies **both** hashes and refuses, with distinct messages, if the source
overlay or the resolved build config has changed or gone missing. `freeze` refuses once an
`evaluation/` directory exists. Nothing is ever re-frozen silently.

To compare several candidate configs (e.g. architecture-only vs "+ corroboration"): **write,
build and freeze every one of them first**, each in its own output directory, then evaluate.

### Phase 7 — evaluate against the reference (optional, after freeze only)

**The reference must be GFF3.** `evaluate` refuses a GTF (or GTF content under a `.gff3` name)
before writing anything, rather than risk silently wrong metrics; convert it first. `gmb-compare`
itself also accepts GTF, but the skill's exon-count split does not.

```bash
for arm in standard candidate; do
  gmb-compare --query $OUT/$arm/finalise/consensus.gff3 --reference <reference.gff3> \
    --reference-fasta <genome.fa> [--seqname-map <map.tsv>] \
    --evaluation-mode protein_coding --output-dir $OUT/$arm/comparison
done
python $SKILL/scripts/gmb_clade.py evaluate --freeze $OUT/config_freeze.json \
    --comparison $OUT/candidate/comparison --baseline-comparison $OUT/standard/comparison \
    --reference <reference.gff3> --out $OUT
```

Report both arms and the delta: CDS exact, multi-exon and single-exon CDS exact, CDS intron
chain, Exact Match, locus detection, missed, novel, structural and strand mismatch, merges,
splits, and per-reference-gene improvements : regressions. Report **CDS exact / multi-exon CDS
exact** as the headline; Exact Match also measures UTR extent and is not coding accuracy.

**Do not optimise one number.** Look for trade-offs — higher detection with worse structures,
fewer novel genes but more missed, a better gene count with worse CDS accuracy — and require a
biological explanation for any preference. If the candidate is not better than `standard`,
**recommend `standard`**; that is a valid, useful result.

If the evaluation suggests a change, it is a **new hypothesis**: new file name
(`<clade>_v2.yaml`), build it, freeze *its* resolved config in a new output directory, then
evaluate — and it must be judged on a genome it was *not* developed on. Never edit the frozen
config in place.

### Phase 8 — validation status

| status | meaning |
|---|---|
| `PROVISIONAL` | derived from input evidence, passes hard QC, **not** evaluated against a trusted reference |
| `SINGLE_GENOME_VALIDATED` | passes hard QC and performs acceptably against a reference on **one** representative genome |
| `CLADE_VALIDATED` | frozen after development on genome A, then run **unchanged** on at least one genuinely independent genome B of the clade, and acceptable there too |

With one genome, the ceiling is `SINGLE_GENOME_VALIDATED` — however good the numbers. With a
second genome: run Phase 1–6 on B with the **frozen** config (same SHA-256), evaluate, and only
then assign `CLADE_VALIDATED`. Never re-tune on B. Always name the next independent genome to
test when the status is below `CLADE_VALIDATED`.

Promotion to a shipped preset (`gmb/configs/<clade>.yaml`, usable as `--preset <clade>`) happens
only at `CLADE_VALIDATED`, only when the user asks, and needs a test pinning its resolved values
(see `docs/creating_a_preset.md`).

### Phase 9 — outputs

| file | content |
|---|---|
| `configs/<clade>.yaml` | the overlay: only settings differing from `standard`, each commented with its measurement |
| `$OUT/gmb_clade_config_report.md` | the report (template below) |
| `$OUT/config_decisions.tsv` | one row per setting considered, changed or not |
| `$OUT/resolved_config.yaml` | from the candidate build |
| `$OUT/slot_map.tsv` | real tool ↔ GMB slot ↔ internal label |
| `$OUT/evidence_summary.md`, `evidence_measurements.json` | Phase 3 |
| `$OUT/preflight_{standard,probe,candidate}/` | full preflight reports |
| `$OUT/run_summary_{standard,candidate}.json` | hard QC + diagnostics |
| `$OUT/config_freeze.json` | the freeze record: source overlay + resolved build config, each with its SHA-256 |
| `$OUT/evaluation/` | only if a reference was supplied |

Report template — begin exactly like this:

```
GMB NEW-CLADE CONFIGURATION

Clade:                  <clade>
Representative genome:  <genome, assembly accession if known>
Config:                 configs/<clade>.yaml   (sha256 <first 16>)
Status:                 PROVISIONAL | SINGLE_GENOME_VALIDATED | CLADE_VALIDATED

Evidence:
    backbone:   <tool> (slot --helixer/--tiberius)  models, multi-exon %, canonical %
    short read: <tool(s)>  ... ; same RNA-seq source? yes/no
    long read:  <tool> ... splice verdict  |  none
    protein:    <tool(s)> ... breadth of support  |  none

Evidence state:         <flags, with one sentence each>

Hard QC:                PASS | FAIL   (standard arm: PASS | FAIL)

Key configuration differences from standard:
    <setting>
        standard:   <value>
        <clade>:    <value>
        reason:     <measurement-based reason>

Remaining uncertainties:
Recommended next validation:
```

Then: the slot map; settings deliberately **left** at `standard` and why; reference-free
comparison of the two arms; evaluation (if any) with trade-offs stated; the validation status
and what would raise it. Every non-default setting needs a reason. "Copied from apicomplexa" is
never a reason.

---

## How decisions differ by evidence state — self-check

Apply this before finishing. Decisions must follow the measurements, not the organism.

| scenario | measured signature | backbone_intron_rescue | long reads | architecture | status ceiling |
|---|---|---|---|---|---|
| 1. strong backbone, weak assemblies | gate declines (backbone well resolved); assemblies overjoined (span p99 ≫ backbone) | keep `"off"` | — | intron cap from backbone ×3; enable splitting | per genomes supplied |
| 2. intron-poor backbone, strong spliced assemblies | gate **fires**; backbone ↔ assembly agreement ~1% | `auto` | — | same rule; ignore short-read chimeric intron tail | per genomes supplied |
| 3. poor long-read splice quality | long-read `splice_quality` FAIL | unaffected | `support_only` or `reject`, guard `true` | unaffected | per genomes supplied |
| 4. one reference genome only | — | — | — | — | **`SINGLE_GENOME_VALIDATED`**, never `CLADE_VALIDATED` |

If two scenarios with different measurements would produce the same config for a non-obvious
reason, re-read the evidence.

---

## Pitfalls seen in practice

- **GMB's internal labels are slot names.** "Scallop" in an output may be AssemblerA. Always ship
  `slot_map.tsv` with the report.
- **`.gtf.gz` is misparsed** by `gmb-build`; supply uncompressed files.
- **The baseline must use `--preset standard`.** A baseline built with a clade preset is not the
  neutral comparator.
- **A subset build is for debugging.** Status requires a whole-genome run.
- **Structural agreement between two assemblers of the same RNA-seq is not independent.**
- **Heavy protein redundancy collapse is normal**; low *support breadth* is the warning sign.
- **Strong backbone does not mean "raise the backbone weight".** No reference-free calibration
  exists for weight magnitudes; that is a hypothesis for Phase 7.

## References

- `docs/configuration.md` — config mechanics, evidence roles, overrides
- `docs/creating_a_preset.md` — the 8-step procedure this skill implements
- `docs/presets.md` — what `standard`, `fungi`, `apicomplexa` encode, and the evidence for each
- `docs/input_contract.md`, `docs/qc.md`, `docs/pipeline_integration.md`
- `configs/new_clade_template.yaml` — every settable field, commented
- `gmb/pipeline/applicability.py` — the rescue gate and its calibration
