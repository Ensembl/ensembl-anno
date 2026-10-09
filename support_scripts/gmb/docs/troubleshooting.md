# GMB troubleshooting

Symptom first. Where to look, and what it usually means.

---

## Where to look

| question | file |
|---|---|
| were my inputs fit to build from? | `preflight/preflight_report.txt` |
| what did the build do? | `build/gmb.log` |
| did the annotation pass QC? | `finalise/fasta_qc_report.json`, `utr_qc_report.json` |
| why was this model chosen? | `build/evidence_attribution.tsv` |
| what settings actually applied? | `build/resolved_config.yaml` |
| what produced this run? | `build/run_manifest.json` |

---

## Install

**`ImportError` / `AttributeError` from pyranges** — the version is too new. GMB pins
`pyranges<=0.1.4`; the API changed after that. `pip install 'pyranges<=0.1.4'`, or install from
conda-forge.

**`gmb-build: command not found`** — the package is installed but its entry points are not on
`PATH`. Activate the environment, or call `$ENV/bin/gmb-build` directly.

**`gmb-preflight: command not found` after upgrading** — new entry point; re-run
`pip install -e .`.

---

## Configuration

**`ValueError: Unknown configuration key: 'scoring.foo'`** — a typo, or a key removed with a
failed experiment. Strict validation is deliberate. Check `configuration.md`.

**`ValueError: Duplicate key 'scoring' at line N`** — a YAML file defines the same top-level
key twice. Plain YAML would silently discard the earlier block; GMB refuses. **Merge them into
one block.**

**`ValueError: Evidence source(s) [...] assigned to two roles`** — the same label appears
under two of `backbone_label` / `shortread_labels` / `longread_label` /
`protein_alignment_labels`. Each source has exactly one role.

**`ValueError: scoring.backbone_intron_rescue must be one of ['off','on','auto']`** — a typo.
Legacy booleans still work (`true` → `on`, `false` → `off`).

**`DeprecationWarning: Config key 'scoring.weights.scallop' is deprecated`** — weights are now
keyed by evidence role. Use `short_read`. The old key still works.

**A preset seems to have no effect** — check `build/resolved_config.yaml`, which is the truth.
If clade values look like the neutral defaults, suspect a duplicate top-level key in the
overlay (now an error) or an overlay overriding them later in the chain.

**`NOTE: --preset not given; using 'fungi'`** — `gmb-build` and `gmb-preflight` fall back to
`fungi` for backward compatibility. Pass `--preset` explicitly (the same value to both).

**`config key 'X' is deprecated|unsupported and has no effect`** — one of 21 keys that load
but change nothing (`configuration.md`, "Keys that have no effect"). Remove it from your
config. `transcriptomic_filter.max_transcript_length` and `allow_single_exon` *are* applied
since 2.0.0.

**`seqid_compatibility [Helixer]: 21/21 sequence name(s) are absent`** — the backbone uses
GenBank accessions. Rename it with `tools/remap_helixer.py` (`input_contract.md`).

**`file_readable` / `coordinates_valid` / `exon_rows_present` FAIL** — the file is empty or not
GTF/GFF3, has impossible coordinates or was made against another assembly, or a backbone /
transcript track has no exon rows (`input_contract.md`, failure behaviour).

---

## Preflight

**`splice_quality [X]: only N% ... canonical`** — the track carries little usable splice
information. For a **long-read** track this is FAIL: set
`scoring.longread_disposition: support_only` (or `reject`), or fix the upstream alignment.
Splice-aware alignment needs the right preset (`-x splice`, `-x splice:hq` for high-accuracy
reads), strand handling matched to the library, and a sensible `-G` maximum intron size.

> A P. falciparum long-read consensus measured **16.4% canonical** against 99–100% for every
> other track. Locus-by-locus it worsened 233 models and improved 228 — no net structural
> benefit — while contributing ~91 extra detected loci.

**`role_resolved [X]: ... resolves to no configured evidence role`** — FAIL. The source would
get the generic unknown-source weight. Add its label to the right role list.

**`seqid_compatibility [X]: N sequence name(s) absent from the genome`** — remap upstream.
Models on unmatched sequences are unusable. GMB does not guess.

**`cross_seqid_id_collisions [X]: N transcript ID(s) reused`** — WARN only; GMB namespaces
them. A large number means the upstream tool numbers models per sequence (Tiberius does).

**`strand_completeness [X]: N feature row(s) have no strand`** — unstranded features are
excluded from selection. Usually harmless below a few percent.

**`longread_collapsed [X]: N long-read models (… per Mb) look like per-read alignments`** —
WARN. The track is raw read alignments, not transcript models; a build would score every read.
Collapse it (`gmb-longread-consensus`, `longread_consensus.md`) or leave it out.

---

## Build

**Exit code 1 with `ERROR: FASTA QC failed`** — **this is a QC verdict, not a crash.** The
annotation is fully written and `build/run_manifest.json` was written before the exit. Read
`build/fasta_qc_report.json` for the failing check. Do **not** hand this build over.

**`Backbone intron rescue [auto]: disabled -- ...`** — the gate declined and said why. Usually
correct: with a well-resolved backbone the rule destroys models that were already right
(1.8% CDS-exact in fungi). Override with `on` only deliberately.

**Build is slower than expected** — GMB's hot path is single-threaded; extra cores do not
help. Runtime scales with **evidence volume**, not genome size: a 732 MB protein track with
1.69 M alignments dominates a 40 Mb genome. Use `--seqname` or `--sample-loci` for a smoke test.
Check preflight for `longread_collapsed` — a raw per-read long-read track alone makes a build
impractical. Runtimes recorded before 2.0.0 include two costs since removed (see
`known_issues.md`).

**Out of memory** — peak RSS follows the largest evidence track. Measured 0.9–1.4 GB
(P. falciparum) and 2.3–2.5 GB (Z. tritici) for `gmb-build`; the former `gmb-compare` peaked at **4.2 GB**
with a large protein track. Size for the compare stage if you run it.

---

## Output quality

**Far more genes than the reference** — over-prediction/over-splitting. Check the split count
and the intron count against the reference; a much higher intron count means real genes are
being fragmented. (Z. tritici: 16,137 genes against 10,931, and 73% more introns.)

**Far fewer multi-exon genes than expected** — the backbone is probably collapsing introns.
Compare backbone vs assembled multi-exon fractions in the preflight report. If the backbone is
low and assemblies are high and canonically spliced, `backbone_intron_rescue: auto` is the
intended remedy.

**Non-canonical introns in the output** — attribute them. If you supplied a long-read track,
that is the first suspect: one P. falciparum annotation was 94.2% canonical with the track and
**100.0%** without.

**Exact Match looks terrible but CDS exact is fine** — expected, and not a defect. Exact Match
includes UTR extent, and GMB emits fewer UTRs than a curated reference. In one build, 1,206
reference genes (22.7%) had a byte-identical CDS yet failed Exact Match on UTR extent alone.
**Quote CDS exact.**

**Locus detection changed after upgrading** — earlier builds had gene records wider than their
own transcripts, which inflated locus detection (the comparator pairs genes by span overlap).
Corrected figures are lower and right. CDS exact and Exact Match are unaffected.

---

## Something still looks wrong

For a specific locus, start with its `evidence_attribution.tsv` row: it records the evidence
sources, protein support, the selection reason and the rescue flag — normally enough to
reproduce and explain the decision. Then check `resolved_config.yaml` for the settings that
applied, and `run_manifest.json` for the inputs and their hashes.
