# Targeted fungal validation — *Zymoseptoria tritici* (GCA_000219625.1)

Release preparation for GMB 2.0.0, 2026-10-09. **No whole-genome GMB build was run.** The
evidence below comes from input analyses, component-level runs of GMB functions,
`gmb-preflight`, regional builds of nine windows (22% of the genome), and comparisons of
existing outputs with ensembl-genes `annotation-qc pairwise-compare`.

The community reference (JGI models distributed by Ensembl Fungi) is a benchmark, not truth:
9.4% of its CDSs lack an ATG start and 11.8% a terminal stop, and it has one isoform per gene.
Helixer predictions absent from it are reported with their independent evidence, not counted
as errors.

Scripts, small result files and exact commands: `diagnostics/gmb_release_prep_20261009/`
(alongside this repository).

---

## 1. Inputs

| role | file | verified content |
|---|---|---|
| genome | `reference_data/…GCA_000219625.1.fa.bgz` | 21 chromosomes `1`–`21`, 39,686,251 bp |
| reference (evaluation only) | `reference_data/…GCA_000219625.1.gff3.bgz` | 10,931 protein-coding genes (one mRNA each), 147 ncRNA, 11 pseudogenes |
| **backbone** | `input_data/…_helixer.gff3` | Helixer `fungi_v0.3_a_0100`; 14,757 genes / mRNA; **sequences named by GenBank accession** (`CM001196.1`…) |
| short-read | `…_scallop_annotation.gtf`, `…_stringtie_annotation.gtf` | 13,169 / 18,360 transcripts; both StringTie-format; 0.4% / 1.7% of rows unstranded |
| protein | `…_orthodb_annotation.gtf`, `…_genblast_annotation.gtf` | genBlastG; 1,642,866 / 11,563 alignments (the latter = UniProt) |
| long-read | `…_minimap2_annotation.gtf` | 2,052,331 **per-read** alignments — not usable as GMB input |
| assembly report | `…_assembly_report.txt` | maps `CM001196.1`…`CM001216.1` → `1`…`21`, lengths identical |

**Helixer checks.** Remapped with `tools/remap_helixer.py` (21/21 sequences, 0 unmapped); the
remapped file is record-for-record identical to the bundled 500 kb test fixture over
chr1:1–500,000. GFF3 structure is clean: one mRNA per gene, gene span = mRNA span, every CDS
inside an exon, every CDS a multiple of 3 starting ATG, GFF3 phases correct, no internal
stops, no UTR/CDS overlap, nothing beyond a sequence end; 1 model lacks a terminal stop. CDS
includes the stop codon, as in the reference. Unrenamed, preflight fails the track (21/21
names absent); renamed, the full bundle passes preflight (36 pass, 5 known warnings, 0 fail;
Helixer 71.9% multi-exon, 93.3% canonical splice sites).

## 2. Whole-genome Helixer vs reference

`annotation-qc pairwise-compare --evaluation-mode protein_coding
--reference-transcript-biotypes protein_coding`. **R** = 10,931 reference protein-coding
genes; **Q** = 14,757 Helixer genes. Matching criteria (annotation-qc): *CDS
coordinate-exact* = start, every splice site and stop identical; *CDS structural exact* =
identical CDS intron chain and ≥ 0.8 reciprocal CDS overlap (ends may differ; the historical
GMB "CDS exact"); *CDS intron chain* is evaluated on the 7,629 reference genes with ≥ 2 CDS
segments; a locus is *recovered* when a Helixer gene overlaps it on the same strand.

Reference genes, mutually exclusive outcomes:

| outcome | all (R = 10,931) | core 1–13 (10,277) | accessory 14–21 (654) |
|---|---|---|---|
| CDS coordinate-exact | 3,954 (36.2%) | 3,866 (37.6%) | 88 (13.5%) |
| CDS structural exact, ends differ | 946 (8.7%) | 926 (9.0%) | 20 (3.1%) |
| overlapped, structure differs (incomplete / partial) | 5,668 (51.9%) | 5,225 (50.8%) | 443 (67.7%) |
| only an opposite-strand prediction | 103 (0.9%) | 78 (0.8%) | 25 (3.8%) |
| **missed entirely** | **260 (2.4%)** | 182 (1.8%) | **78 (11.9%)** |
| locus recovered (sum of first three) | 10,568 (96.7%) | 10,017 (97.5%) | 551 (84.3%) |
| CDS intron chain exact (of multi-CDS genes) | 3,205 / 7,629 (42.0%) | 3,162 / 7,203 (43.9%) | 43 / 426 (10.1%) |

Helixer side: 3,279 of 14,757 genes (22.2%) have no same-strand reference counterpart (3,053
overlap no reference gene at all, 217 only one on the opposite strand, 9 a non-coding gene);
core 2,522 of 13,231 (19.1%), accessory 757 of 1,526 (49.6%). 511 reference genes are split
over ≥ 2 Helixer genes; 130 Helixer genes merge ≥ 2 reference genes. One-to-one exact-CDS
precision / recall / F1 = 26.8% / 36.2% / 0.308.

Per chromosome: core chromosomes are uniform (locus recovery 94.8–98.7%, CDS-exact
35.6–40.2%; lowest recovery chr7, 32 missed); every accessory chromosome is far lower
(recovery 81.6–87.9%, CDS-exact 5.2–26.6%). Table:
`results/helixer_wg_strata/per_chromosome.tsv`.

**Was the 500 kb fixture representative?** Partly. Locus recovery (97.3%) matches the core
average (97.5%); exact CDS (29.1%) is below the chr1 (37.3%) and core (37.6%) rates, though
within sampling error for 110 genes (95% CI about ± 9 points); novel predictions are higher
(25.6% vs 19.1%). It is core-only and therefore says nothing about the accessory chromosomes,
which behave very differently. It also overstated locus fragmentation (12.6% of multi-exon
candidates vs 5.7% genome-wide).

## 3. Evidence for missed and potentially novel loci

Same-strand support; *strict* = protein alignment covering ≥ 50% of the CDS, or (multi-exon)
every intron observed exactly in a short-read assembly.

**Missed reference genes (260).** Short (median CDS 269 bp vs 1,071 for all reference genes)
and mostly single-segment (70%). 75.0% have short-read support, 64.6% protein, 44.6% both;
only 13 (5.0%) have neither — so most misses are a backbone limitation, not an evidence gap.
Accessory: 9 of 78 lack both.

**Helixer genes:**

| | n | protein overlap | short-read overlap | both | neither | strict protein | all introns in RNA-seq | neither strict |
|---|---|---|---|---|---|---|---|---|
| novel, core | 2,522 | 84.1% | 94.7% | 81.8% | 2.9% | 75.3% | 29.4% | 22.9% |
| novel, accessory | 757 | 51.8% | 80.8% | 45.3% | 12.7% | 31.6% | 4.9% | **66.4%** |
| reference-overlapping, core | 10,709 | 99.2% | 97.8% | 97.0% | 0.1% | 98.4% | 58.8% | 1.5% |
| reference-overlapping, accessory | 769 | 89.9% | 91.3% | 82.6% | 1.4% | 80.5% | 15.2% | 19.0% |

Most core novel predictions carry substantial protein evidence and are plausible unannotated
genes; most accessory ones are weakly supported. Caveats: positional protein support is
near-saturated (it covers 99% of reference-overlapping genes), and protein homology may come
from transposable-element proteins that reference gene sets exclude — a repeat annotation is
needed to separate them.

Reference coverage by any evidence (earlier analysis): protein 98.1%, short-read 97.4%, both
95.7%, neither 0.2%; accessory protein 90.5%, short-read 89.6%. Only 57.3% of multi-exon
reference genes have every intron in a short-read assembly (accessory 14.8%).

## 4. What GMB added over its backbone before 2.0.0

Same comparator, whole genome: Helixer alone vs the September 2026 GMB build (fungi preset,
pre-2.0.0 code, protein validation on). Existing output; no new build.

| | Helixer | GMB (Sept) |
|---|---|---|
| genes | 14,757 | 16,137 |
| CDS coordinate-exact | 3,954 (36.2%) | 3,997 (36.6%) |
| CDS structural exact | 4,900 (44.8%) | 4,955 (45.3%) |
| CDS intron chain | 3,205 | 3,257 |
| missed | 260 | 224 |
| novel / opposite-strand-only novel | 3,279 / 217 | 3,964 / 441 |
| splits / merges | 511 / 130 | 661 / 239 |
| one-to-one F1 | 0.308 | 0.295 |

In this evidence state GMB's gain over Helixer is small (+43 exact CDS, −36 missed) and came
with more genes, splits and merges. §5 explains why.

## 5. Locus clustering

**Genome-wide, component level** (46,149 filtered candidates, 39,285 multi-exon):

| | `exon_overlap` | `transcript_linked` |
|---|---|---|
| candidates split across loci | 2,223 (5.7% of multi-exon) | 0 |
| loci | 4,408 | 2,602 |
| transcripts per locus, median / p99 / max | 2 / 182 / 475 | 2 / 240 / 494 |
| largest locus span | 393 kb | 403 kb |
| selection-cost proxy Σk² | 5.75 M | 6.39 M (+11%) |
| clustering time | 0.5 s | 0.8 s |

Both modes cluster unstranded, and exon-overlap loci already contain many neighbouring genes
(chained through overlapping UTRs); genes are separated inside a locus by `select_isoforms`
(strand and overlap tests), so `transcript_linked` does not newly merge opposite-strand or
neighbouring genes into one *gene*. Its specific risk is that a chimeric transcript is kept
whole (exon-overlap fragments it) — which is why it must be paired with a chimera filter.

**Read-through chimeras.** After the previous fungal filters, 4,945 short-read transcripts
(15.8%) still spanned > 20 kb (up to 238 kb); 4,859 of them overlapped ≥ 2 reference genes on
the same strand. No intron exceeded 3 kb — they cross short intergenic gaps, often as one long
"exon" — so `max_intron_length` removed none. The fungi preset's 20 kb
`max_transcript_length` was meant to remove them but was never applied (now implemented).

**Regional builds** — nine windows: core 1:2–3 Mb, 3:1–2 Mb, 5:0.5–1.5 Mb, 7:1–2 Mb,
9:1–2 Mb, 12:0.3–1.3 Mb (R = 1,865); accessory chr14, 18, 21 whole (R = 236). Whole
transcripts only; fungi preset; protein validation off.

| core windows | Q | CDS coord-exact | CDS intron chain | missed | splits | merges | build s |
|---|---|---|---|---|---|---|---|
| Helixer alone | 2,322 | 719 | 557 | 44 | 69 | 18 | — |
| pre-2.0.0 default (exon_overlap, no span filter) | 2,517 | 771 | 619 | 29 | 484 | 465 | 51.9 |
| exon_overlap + 20 kb filter | 2,543 | **761** | 612 | 31 | 424 | 362 | 44.5 |
| transcript_linked, no filter | 2,510 | 772 | 620 | 29 | 484 | 459 | 50.6 |
| transcript_linked + 20 kb filter | 2,510 | 774 | 623 | 32 | 406 | 334 | 44.0 |
| … + one isoform per gene | 2,510 | 719 | 557 | 41 | 71 | 18 | 42.4 |
| **2.0.0 fungi preset** (linked + filter + detached-isoform removal) | 2,510 | 770 | 617 | 35 | **327** | **245** | 45.9 |
| rejected: detached isoforms split into new genes | 2,638 | 774 | 623 | 34 | 406 | 247 | 46.5 |

| accessory windows | Q | CDS coord-exact | missed | splits | merges |
|---|---|---|---|---|---|
| Helixer alone | 473 | 25 | 23 | 33 | 3 |
| pre-2.0.0 default | 586 | 27 | 15 | 63 | 16 |
| **2.0.0 fungi preset** | 583 | 29 | 16 | 51 | 10 |

Findings:

1. **The span filter and `transcript_linked` must go together.** Removing the chimeras also
   removes the transcripts that happened to bridge introns, so under `exon_overlap` candidate
   fragmentation rises from 242 to 1,138 and CDS accuracy falls (771 → 761).
2. **Detached isoforms.** Read-through alternates admitted because their long span overlapped
   the primary were trimmed by validation into pieces up to tens of kb away, leaving genes
   whose transcripts do not overlap (244 of 481 merged genes before 2.0.0). Removing them
   cuts merges 334 → 245 and splits 406 → 327 at a cost of 4 exact CDS; of the 7 reference
   genes newly "missed", 6 had 0–3% exon overlap before (artefacts of inflated gene spans) and
   1 is a real partial loss. Re-homing the pieces as new genes instead keeps those 4 matches
   but emits 128 genes that never went through selection and raises splits; rejected.
3. **Selection is backbone-led.** With one isoform per gene the output equals Helixer's. The
   remaining 255 merges are almost all multi-isoform genes (226 "span-only": overlapping
   alternates reaching a neighbour) — an isoform-policy question (`../known_issues.md`).
4. All nine release-candidate window builds passed FASTA QC; `gmb-finalise` on two of them
   passed FASTA and UTR QC with one canonical transcript per gene; peak RSS 0.42 GB.

**Recommendation.** `transcript_linked` with the 20 kb span filter and detached-isoform
removal as the fungal default for the production-scale test, on structural grounds: no
fragment scoring, a documented chimera filter that actually applies, and no gene joining
non-overlapping transcripts — with exact-CDS agreement unchanged, not improved. Risks:
genome-wide behaviour unconfirmed; genuine fungal transcripts > 20 kb are dropped from the
short-read candidates (1 of 10,931 reference genes exceeds 20 kb; the backbone can still call
it); per-locus selection cost rises ~11%. `exon_overlap` remains available. Keeping it as the
default would keep fragment scoring and either a misleading inert setting or, with the filter
applied, lower accuracy.

## 6. Earlier component findings still valid

- Fungal hard thresholds remove no reference gene (CDS ≥ 99 / 150 / 90 bp; UTR caps above every
  reference UTR); 34 of 17,616 reference introns (0.19%) exceed `max_intron_length` 3 kb.
- The protein filter is strand-safe (identical result when run per strand); its score
  thresholds remove nothing for genBlastG input.
- UTR end support used to count ends on other chromosomes (5.8% / 6.7% of 5′ / 3′ ends); fixed.
- Long-read track: 99.1% of 3.2 M read-introns canonical on the declared strand, but per-read.

## 7. What still needs the production-scale run

Genome-wide accuracy and gene count of the 2.0.0 preset; runtime and memory with all tracks;
canonical-only vs all-isoform evaluation; behaviour on accessory chromosomes at scale; the
effect of `transcript_linked` on loci larger than the windows contain.
