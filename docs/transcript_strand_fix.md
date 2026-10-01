# Minimap2 transcript-strand handling in the long-read path

## The bug

The long-read path built transcript models like this:

```
FASTQ -> minimap2 -> SAM -> paftools.js splice2bed -> BED12 -> _bed_to_gtf -> annotation.gtf
```

`paftools.js splice2bed` sets the BED strand column from the SAM FLAG:

```js
a1 = [t[2], parseInt(t[3])-1, null, t[0], 1000, (flag&16)? '-' : '+'];
```

The `ts:A:` tag is never read. The SAM FLAG records **alignment orientation** — whether the read
had to be reverse-complemented to match the reference — which for an unstranded cDNA library is
effectively a coin flip and says nothing about which strand the transcript came from.

`_bed_to_gtf` then copies the BED strand straight into the GTF, so the error propagates unchanged
into every transcript model.

The effect is measurable. On the *P. falciparum* run ERR9652456, of the junctions that are
canonical on *some* strand, 50.5% matched the BED strand and 49.5% were canonical only on the
opposite strand — a coin flip. Scoring the same alignments three ways:

| strand source | canonical splice fraction |
|---|---:|
| motif on either strand (upper bound) | 28.56% |
| Minimap2 `ts` (correct) | 21.90% |
| SAM FLAG (what `splice2bed` does) | 14.38% |

## Minimap2 `ts` semantics

Minimap2 infers the transcript strand from the splice-site motifs it selected and reports it as
`ts:A:+` or `ts:A:-`. In SAM output the value is relative to **the read as stored in the record**,
so the genomic transcript strand is:

```
genomic transcript strand = ts XOR alignment_is_reverse
```

This was verified empirically rather than assumed. Taking as ground truth the genome's own
donor/acceptor dinucleotides at each intron, restricted to single-intron alignments with an
unambiguously canonical motif (n = 36,652):

| rule | agrees, forward alignments | agrees, reverse alignments |
|---|---:|---:|
| `ts` XOR `is_reverse` | 75.45% | 76.18% |
| `ts` used directly | 75.45% | 23.82% |
| SAM FLAG | — | ~51% (chance) |

The XOR rule is stable across both orientations while the direct reading collapses on reverse
alignments, which identifies the convention. The ~24% shortfall is symmetric across orientations,
so it is Minimap2's strand-inference accuracy on that (poor-quality) dataset, not a convention
error.

## What changed

Two new modules and one rewired call site.

**`transcriptomic_annotation/transcript_strand.py`** — the rule as a small pure helper:

```python
transcript_strand(ts_tag, is_reverse, on_missing=MissingTsPolicy.UNSTRANDED) -> "+" | "-" | "."
```

**`transcriptomic_annotation/splice_to_bed.py`** — a generic SAM/BAM → BED12 converter that uses
the helper. It replaces the `paftools.js splice2bed` subprocess, because the defect is inside that
tool and cannot be configured away. Format is detected from file content, so plain SAM, gzipped SAM
and BGZF/BAM all work. No new third-party dependency: SAM text and the BAM binary layout are both
read natively, matching the project's existing `psutil`/`numpy`-only footprint.

**`transcriptomic_annotation/minimap.py`** — the `splice2bed` subprocess call is replaced by
`splice_to_bed(...)`; a new `on_missing_ts` argument and `--on_missing_ts` CLI flag expose the
policy. `paftools_bin` is retained in the signature for backwards compatibility with existing
callers but is no longer used, and is no longer required on `PATH`.

### Behaviour when `ts` is absent or invalid

**The default is to emit `.` (unstranded), never to guess.** Minimap2 only reports `ts` when it has
splice-site evidence, so an unspliced alignment legitimately has no transcript strand; inventing one
from the FLAG is the original bug. Three policies are available:

| `on_missing_ts` | behaviour |
|---|---|
| `unstranded` (default) | emit `.` |
| `error` | raise `MissingTranscriptStrandError` |
| `flag` | fall back to the SAM FLAG — reproduces the old, incorrect behaviour; opt-in only |

Whatever the policy, the conversion is never silent: `ConversionStats` counts how many strands came
from `ts`, how many were unstranded and how many came from the FLAG, and a spliced alignment with no
usable `ts` (the suspicious case) raises a `WARNING`.

**Operational impact.** Most long reads are unspliced and therefore carry no `ts`. On a 25,925-record
sample of ERR9652456, 7,496 records (29%) got a `ts`-derived strand and 18,429 (71%) became `.`;
1,142 of the unstranded ones were spliced. Unstranded features are excluded from selection by
downstream consumers such as GMB, so this default trades a large number of *wrongly* stranded
single-exon models for a large number of *unstranded* ones. That is the honest representation of
what the aligner actually determined, but it is a real change in what downstream tools receive and
should be reviewed before the next production run.

## What deliberately did NOT change

- **Minimap2 parameters.** The command line is untouched (`-G`, `-u b`, `--secondary=no`, `-ax splice`).
- **Alignment filtering policy.** Secondary alignments (FLAG 0x100) are still dropped and
  supplementary alignments (FLAG 0x800) are still emitted as separate BED records, exactly as
  `splice2bed` does without `-m`. `keep_secondary` exists as a separate, explicit knob; the strand
  fix does not touch filtering.
- **BED12 record shape.** Name (including the `/1`, `/2` mate suffix), score 1000, thickStart and
  thickEnd equal to the chrom range, and the `itemRgb` colour convention (`0,128,255` single primary,
  `255,0,0` multi-primary, `0,192,0` secondary) are all preserved.
- **Coordinate semantics.** SAM POS is 1-based inclusive; BED is 0-based half-open
  (`chromStart = POS - 1`); GTF is 1-based closed. `_bed_to_gtf` and `_bed_block_to_exons` are
  unchanged — they were already correct.
- **GMB.** Nothing in `support_scripts/gmb` was touched.

One intentional, documented difference: `splice2bed` advances the reference offset only on CIGAR `M`
and `D`, ignoring `=` and `X`. This converter handles `=` and `X` correctly. Minimap2 emits `M` by
default, so this cannot change output for the current pipeline; it only prevents silently wrong
blocks if a future aligner or option emits the extended CIGAR operations.

## Verification

The converter was checked against an independent `pysam`-based reimplementation on 25,925 real
alignment records from ERR9652456: the BED output is byte-identical, and the native SAM and native
BAM readers agree byte-for-byte with each other.

Run the tests with:

```bash
PYTHONPATH=src/python python -m pytest tests -q
```
