# Evidence integration — *Zymoseptoria tritici* (GCA_000219625.1)

Investigation and tuning, branch `gmb/evidence-tuning` (from `feature/gene_model_builder`
e005a86), 2026-10-09. Nothing is committed. **No whole-genome GMB build was run.** The
reference is a benchmark, not truth: it has one isoform per gene and its own errors (see
`z_tritici_validation.md`). Scripts, configurations and per-window results:
`diagnostics/gmb_evidence_tuning_20261009/` (alongside this repository).

**Question.** Given substantial protein and transcriptomic evidence, why does GMB add so little
beyond Helixer, and which demonstrated changes let that evidence contribute without
sacrificing gene-set quality?

**Answer in one paragraph.** Under the fungi preset the evidence cannot change a Helixer gene
model. Helixer with positional protein overlap scores 5.1; the best structure supported only
by short-read transcripts scores 4.5. Transcript structures therefore appear only as
alternate isoforms. Those alternates (and the canonical stage, which counts Scallop and
StringTie as two sources) then *damage* the gene set: on held-out windows the shipped
canonical output is less accurate than Helixer alone. Protein alignments carry spliced CDS
structures but are reduced to "overlaps something on the same strand", which 93% of backbone
candidates do. Three contained selection changes fix this. Junction-supported primary choice
uses where RNA-seq observes the introns. Isoforms must share coding sequence. Canonical
selection stops counting the two short-read assemblers twice. Together they give held-out core
canonical CDS-exact **710 → 733** (Helixer → updated), and **754** with the collapsed long-read
track, which is the strongest evidence available and is currently not usable as delivered.
Splits and merges return to Helixer's level. What the evidence cannot do: recover most genes
Helixer misses (they sit inside long read-through assemblies) or tell incidental antisense
ORFs from genes.

---

## 1. Method

**Windows.** Fixed before any result was seen. *Development* (tuning allowed): the nine
release-prep windows: core 1:2–3 Mb, 3:1–2 Mb, 5:0.5–1.5 Mb, 7:1–2 Mb, 9:1–2 Mb,
12:0.3–1.3 Mb; accessory chr14, 18, 21 whole (R = 1,865 core + 236 accessory reference
protein-coding genes). *Held-out* (first evaluated after every setting was fixed): core
2:1–2 Mb, 4:1–2 Mb, 6:1–2 Mb, 8:1–2 Mb, 10:0.3–1.3 Mb, 13:0.1–1.1 Mb; accessory chr16, 19
whole (R = 1,874 + 175). Inputs keep whole models only; genome = whole chromosome.

**Builds.** `gmb-build --preset fungi` (transcript_linked, 20 kb span filter, detached-isoform
removal), protein validation off unless stated, then `gmb-finalise`. Clustering and
isoform-count policy are held fixed within every comparison.

**Evaluation.** ensembl-genes `annotation-qc pairwise-compare --evaluation-mode protein_coding
--reference-transcript-biotypes protein_coding --region <window>`. Three views of each build:
*all isoforms*, which rewards extra isoforms because any isoform may match; *build primary*
(`.t1`); and *canonical* (`gmb-finalise`, `--query-transcript-selection canonical`), which is
the handover. Definitions as in `z_tritici_validation.md`. *CDS coord-exact*: start, every
splice site and stop identical. *CDS intron chain*: over reference genes with ≥ 2 CDS
segments. *Missed*: no overlapping prediction. *Novel*: query gene with no same-strand
reference counterpart. One-to-one exact-CDS *F1*. Intron precision/recall are over CDS introns.

---

## A. Evidence integration diagnosis

### How each evidence type travels through GMB

| | short-read assemblies (Scallop, StringTie) | long reads (Minimap2) | protein alignments (OrthoDB, UniProt genBlastG) | backbone (Helixer) |
|---|---|---|---|---|
| parsed | GTF exons; no CDS (both are StringTie `--merge` output) | GTF exons, **per read** (2.05 M); unusable without `gmb-longread-consensus` | GTF exons = spliced **coding** segments | GFF3 exons + CDS |
| filtered before candidacy | introns > 3 kb; span > 20 kb (15.8% of transcripts, mostly read-through); gap splitting | same | fragment/span/redundancy filters; score thresholds inert (no per-exon scores) | CDS < 90 bp |
| candidate structures? | **yes** | yes | **no — never** | yes |
| CDS | longest ATG ORF ≥ 33 codons, **one per transcript** | same | — | its own CDS, kept |
| support signal | role weight 1.0 per source, multi-source bonus 0.5 | 1.0 | +2.0 if *any* alignment overlaps on the same strand (86% of all candidates, 93% of backbone ones) | role weight 3.1 |
| affects ranking? | only as alternate; can never outrank backbone | with SR it can (3 sources + protein = 6.0 > 5.1) | almost never discriminates | wins every locus it predicts |
| affects retention? | single-source multi-exon kept only with protein overlap | same | gate for single-source structures | always kept |

### What the evidence actually changes (dev windows, unchanged code)

| arm | core primary CDS-exact | core all-isoform CDS-exact | core missed | genes / transcripts |
|---|---|---|---|---|
| A Helixer alone | 719 | 719 | 44 | 2,322 / 2,322 |
| B GMB, Helixer only | 719 | 719 | 44 | identical to A |
| C + protein | 719 | 719 | 44 | identical to A |
| D1 + Scallop / D2 + StringTie | 719 | 719 | 44 | identical to A |
| D + both short-read | 719 | 722 | 42 | 2,375 / 2,484 |
| D3 + collapsed long reads | 719 | 719 | 44 | identical to A |
| E all evidence (shipped) | **719** | 770 | 36 | 2,510 / 3,960 |

**Every primary transcript at a Helixer locus is Helixer's.** Protein evidence alone, one
transcript track alone, and long reads alone change nothing. A single-source transcript
structure is retained only with protein overlap, and even then it scores 3.0 against
Helixer's 5.1. The only primary-level change is new genes where Helixer has none (missed
44 → 36). The 51 extra CDS-exact genes in the all-isoform view are alternates, and the same
alternates take splits/merges from 69/18 to 327/245.

### Where useful evidence is lost (dev, per-structure audit of `select_isoforms`)

* **Candidate availability caps what reweighting could achieve.** 780 of 1,861 core reference
  genes have a reference-exact candidate structure: 719 from Helixer and 61 (3.3%) from
  transcripts only. Of those 61 genes' 92 exact structures, 65 became alternates and 27 were
  dropped. Perfect reselection could therefore gain at most +61.
* **Selection discards the junction evidence that discriminates.** Where a Helixer primary
  competes with a transcript structure and exactly one is reference-exact (300 pairs),
  current score picks Helixer every time: right 192, wrong 108. "Fewer introns that no
  transcript observes" is right 91 / wrong 16 when it decides, and spliced-protein
  compatibility right 92 / wrong 28. Positional protein support and CDS length do not
  discriminate.
* **Scores are coarse and saturated.** 22 distinct values over 7,225 structures, 64% at
  exactly 3.0 or 5.1. 24% of multi-isoform genes have an alternate tied with the primary;
  the tie is broken by input order.
* **Finalisation undoes the build.** Canonical selection, with protein validation off, ties
  on ORF and protein tiers. Its *named-source count* then prefers a Scallop+StringTie isoform
  (2) over the Helixer primary (1). That swapped 199 of 3,093 dev canonical transcripts (184
  to Scallop+StringTie models): 10 became exact, 17 stopped being exact, 165 changed without
  either matching. Held-out core canonical CDS-exact: **Helixer 710, shipped GMB 695.**
* **Isoform grouping chains neighbouring genes.** An alternate joins a gene by sharing an
  intron or > 15% of span. 41% of alternates (674 of 1,662, held-out) share **no coding base**
  with their primary. They are a neighbouring ORF reached through a UTR or read-through
  overlap. This produces the merges.
* **Inert setting.** `scoring.min_cds_bp` (150) never applied; `structurally_valid()` read a
  field nothing set.
* **Protein validation, if enabled as configured for fungi,** penalises every model without
  a DIAMOND hit (weights 0.7/0.3, `min_score` 0.7) by the hard-coded 5.0. The effect is to
  remove most alternates (dev core transcripts 3,960 → 2,650) and change no primary (719);
  missed 35 → 41.

---

## B. Biological findings from representative loci

| locus (window) | category | candidates (source: score) | stage responsible | class | original → updated |
|---|---|---|---|---|---|
| Mycgr3G66023, 1:2.66 Mb − | Helixer structure wrong, RNA-seq right | Helixer 4.6 (acceptor 4 bp off, 5′ CDS extension; 1 intron unobserved); Scallop 3.0 and StringTie 3.0 both reference-exact | ranking (backbone dominance) | algorithm limitation | Helixer primary → Scallop promoted (`transcript_junction_support`), exact |
| Mycgr3G65706, 1:2.13 Mb − | merge / canonical swap | Helixer 5.1 exact; Scallop+StringTie 4.5 with CDS 3 kb away (different ORF) admitted as isoform by span | gene grouping + canonical source count | defect (policy) | canonical = foreign ORF → Helixer kept; foreign ORF in its own coding gene |
| Mycgr3G103686, 3:1.43 Mb + | same mechanism | Helixer 5.1 exact; Scallop+StringTie 4.5 CDS-disjoint "isoform" | same | defect (policy) | lost → exact; exon-skipping alternates that share CDS are kept |
| Mycgr3G95776, 9:1.21 Mb − | first version of the new rule | promoted read-through alternate had its CDS outside Helixer's; validation trimmed it, Helixer removed as detached | promotion rule | bug found during tuning | fixed by the same-coding-locus test (regression test) |
| Mycgr3G93233, 5:1.30 Mb − | regression | StringTie acceptor 3 bp from Helixer's (NAGNAG-type), observed in RNA-seq; reference follows Helixer | junction rule | uncertain (alternative splice site) | exact → partial |
| Mycgr3G64098, 12:0.39 Mb − | regression | Scallop (retained intron) chosen over a near-exact StringTie model by protein net (−2 vs −7) | junction rule tie-break | limitation of the compatibility measure | ic → partial |
| Mycgr3G97543, 14:0.12 Mb − | long reads + different ORF | Minimap2+Scallop+StringTie 6.0 (ORF 430 bp away) outranks Helixer 5.1 (exact) | ranking + grouping | defect (policy) | Helixer demoted → its own gene, canonical exact |
| Mycgr3G97615, 14:0.73 Mb + | isoform / canonical | Minimap2+Scallop truncated (retained-intron) model has 2 evidence classes vs Helixer's 1 | canonical class breadth | accepted trade-off | exact → partial (class breadth kept: net +6 core on dev) |
| Mycgr3G31653, 1:2.48 Mb − | missed by Helixer | 6 short-read transcripts contain the 246 bp CDS but span 29 kb (read-through, span-filtered or longest ORF elsewhere); 120 protein alignments reproduce the CDS exactly | candidate generation (one ORF per transcript; proteins not candidates) | algorithm limitation | partial overlap only, both versions |
| Mycgr3G95798, 9:1.28 Mb + | missed | 8 transcripts (14.9 kb) contain an ATG ORF equal to the reference CDS; GMB takes the longest ORF elsewhere | ORF inference | algorithm limitation | missed, both |
| Mycgr3G45560, 7:1.68 Mb − | missed, no RNA-seq | 148 protein alignments, 6 reproduce the CDS exactly | proteins never candidates | algorithm limitation | missed, both |
| ZT_00027 / ZT_00038 (protein-validation build), 1:2.14 / 2.18 Mb + | antisense "gene" | Scallop-only 2- and 3-exon models, canonical GT-AG on their strand, antisense to reference genes; CDS 384 / 333 bp; Psauron 0.003–0.07, no DIAMOND hit | retention (ORF ≥ 33 codons + positional protein) | insufficient evidence to reject | retained, both |

**What the loci show.** GMB's real strength is that transcripts are scored *whole* and
merged by intron chain, so where RNA-seq and Helixer agree the model is solid. CDS-exact is
55% in the backbone + short-read agreement bucket against 30% for backbone only. Its
weaknesses are all in what happens when they disagree: the backbone always wins the primary;
"isoform" means "nearby"; and the canonical stage counts assemblers rather than evidence.

### Missed reference genes (section 3A)

Of 44 core genes Helixer misses (dev), 18 have a same-strand short-read transcript containing
the whole CDS. In **none** does GMB's inferred ORF reproduce it: these are small (median CDS
237 bp, 70% single-segment) genes inside **long read-through assemblies** (median transcript
span 13–160 kb). In 13 a different ATG ORF on the same transcript is the reference CDS. 30
have protein overlap and 5 a protein alignment with the exact CDS. The shipped build recovers
8 as genes, none exactly. The "75% have RNA-seq support" figure from release prep is overlap,
not transcripts able to define a coding model. **This is mostly an evidence limitation.**
Recovering these genes needs multi-ORF inference with coding-potential support, or
protein-derived candidates. Tested as an analysis only: gap-filling with protein-defined
complete CDSs (≥ 2 identical alignments, no overlapping gene) would add 73 models per nine
windows, of which 4 are reference-exact and 37 have no reference gene at all. Not
implemented.

### Evidence-supported novel genes (section 3C)

Held-out core, updated build: 432 genes without a same-strand reference counterpart. 305 are
Helixer's own, and 46 are transcript-only intergenic models (13 supported by ≥ 2 sources or
long reads; 41 of the 351 intergenic models have a compatible spliced protein alignment).
Separately, **188 genes lie antisense to a reference gene, 123 of them transcript-only.**
Their introns are canonical on the declared strand (131 of 144 multi-exon), so the
transcription is real. The ORFs look incidental: median Psauron 0.10 (Helixer primaries:
82% ≥ 0.9), no DIAMOND hit, and 1 of 177 agrees with a spliced protein. They are present in
both versions. They cannot be separated from the ~15 per 2,100 reference genes where Helixer
predicted the wrong strand and the transcript model is right, by opposite-strand overlap,
protein coverage or backbone transcript support. Psauron < 0.5 would remove ~80% of them and
~65% of transcript-only models that overlap a reference gene: a policy decision, not made
here.

### Splits, merges and isoforms (sections 3D, 3E, 6)

| held-out core | genes | alternates | CDS-disjoint alternates | UTR-only alternates | splits / merges (all isoforms) |
|---|---|---|---|---|---|
| Helixer alone | 2,355 | 0 | — | — | 64 / 23 |
| original | 2,520 | 1,662 | **674 (41%)** | 457 | 342 / 253 |
| updated | 2,496 | 1,711 | 0 | 864 | 65 / 24 |

`transcript_linked` clustering keeps candidates whole as intended, and the 20 kb filter
removes read-through assemblies. The merges come from **isoform admission**, not clustering:
UTR and read-through overlaps admit a different ORF as an "isoform". Requiring a shared coding
base removes them while keeping 1,711 coding-sharing alternates. Canonical choice is now
independent of those alternates: canonical equals the build primary except where protein
validation or evidence-class breadth (e.g. long reads) says otherwise. Removing alternates
altogether would hide transcript evidence: a single-isoform build equals Helixer. Half of
the remaining alternates differ from the primary only in UTR, which is a decision about what
to ship, not about accuracy.

---

## C. Tuning experiments (all configurations tested)

Dev windows, core reference genes (R = 1,865); canonical unless stated.

| configuration | CDS-exact | intron chain | missed | splits / merges | outcome |
|---|---|---|---|---|---|
| Helixer alone | 719 | 557 | 44 | 69 / 18 | baseline |
| E shipped fungi | 706 (all-iso 770) | 553 | 36 | 89 / 18 (all-iso 327 / 245) | baseline |
| F_cb: `prefer_build_primary` only | 719 | 557 | 36 | 71 / 18 | fixes the canonical regression, no gain |
| F_jp v1: `junction_supported`, no coding-locus test | 704 canonical / 725 primary | — | 40 | 97 / 18 | **rejected**: promoted read-through ORFs, lost 4 genes |
| F_jp: with coding-locus test | 706 canonical / 733 primary | 553 / 574 | 36 | 89 / 18 | gain lost at canonical without F_cb |
| F_jp + F_cb | 733 | 574 | 36 | 70 / 18 | +17 net genes (19 + 6 ic gained, 3 lost) |
| + `min_cds_bp: 150` | 733 | 574 | 36 | 68 / 18 | −46 core / −36 accessory genes, nearly all novel; no reference CDS < 150 bp |
| + `isoform_cds_overlap: new_gene` | 750 | 599 | 34 | 68 / 18 (all-iso 72 / 20) | +124 genes, mostly novel: rejected |
| **+ `isoform_cds_overlap: drop` (= overlay)** | **750** | **599** | 36 | 68 / 18 (all-iso 72 / 20) | adopted |
| overlay, build-primary tier above class breadth | 750 (767 with long reads) | 599 | 36 | 68 / 18 | rejected: −6 core with long reads |
| E + protein validation (fungi weights, penalty 5.0) | 722 | 561 | 41 | 72 / 18 | prunes isoforms only; primary unchanged |
| protein-defined gap-filling candidates (analysis) | — | — | — | — | 4 exact of 73 added: not implemented |
| antisense guards: opposite-strand overlap / in-frame protein (analysis) | — | — | — | — | remove as many correct genes as antisense ones: not implemented |
| E + collapsed long reads | 738 | 590 | 24 | 92 / 19 | long reads help even unchanged |
| **overlay + collapsed long reads** | **773** | **631** | 24 | 69 / 18 | best |

Pairwise-rule simulation on dev (before implementation, per backbone primary, choosing one
replacement): "fewer unobserved introns + complete ORF + no worse protein + CDS ≥ 90%"
predicted +35 / −3 exact for 224 switches. The implementation was then held fixed for the
held-out evaluation.

---

## D. Implemented improvements

All new behaviour is opt-in. Shipped presets produce the same output as before: the golden
fixture passes unchanged. The four settings are combined in
`configs/fungi_evidence_integration.experimental.yaml`.

| change | files | rationale | regression coverage |
|---|---|---|---|
| junction-level evidence: introns no transcript observes; spliced-protein compatible / incompatible alignments | `pipeline/junction_support.py` (new), `builder.py` | the only measured discriminating signals; always reported | `TestJunctionSupport` |
| `scoring.primary_selection: junction_supported` (+ `junction_primary_min_cds_fraction`) | `scoring.py`, `config.py`, `standard.yaml` | backbone always wins otherwise | `TestJunctionSupportedPrimary` (7 cases incl. the read-through defect), fixture build |
| `scoring.isoform_cds_overlap: off \| drop \| new_gene` | `scoring.py`, `config.py`, `standard.yaml` | 41% of alternates encode a different ORF; source of merges | `TestIsoformCdsOverlap` (backbone never dropped) |
| `canonical_selection.prefer_build_primary` | `canonical_selection.py`, `config.py`, `standard.yaml` | named-source count swapped canonicals, net harmful | `TestPreferBuildPrimary` |
| `scoring.min_cds_bp` now applied; shipped value 0 | `scoring.py`, `config.py`, `standard.yaml`, `configs/fungi_default.yaml` | was silently inert; 0 keeps validated output | `TestMinCdsBp` |
| `protein_validation.penalty` (default 5.0) | `scoring.py`, `config.py`, `standard.yaml` | was hard-coded | `TestProteinValidationPenalty` |
| attribution columns `introns_without_transcript_support`, `protein_alignments_compatible`, `protein_alignments_incompatible`; `selection_reason=transcript_junction_support` | `builder.py`, `tests/test_release_prep_regressions.py` | auditability | pinned column list (appended) |

Tests: 889 passed, 3 skipped (was 862 / 3). Compatibility: `evidence_attribution.tsv` gains
three trailing columns; `standard.yaml` shows `min_cds_bp: 0` instead of the never-applied
150; new keys appear in `resolved_config.yaml`. Documentation: `configuration.md`,
`output_contract.md`, `canonical_selection.md` (tier list corrected and extended),
`known_issues.md`.

---

## E. Before/after benchmarking

Canonical = handover output. "Updated" = fungi preset + experimental overlay. Long reads =
the per-read Minimap2 track collapsed per window with `gmb-longread-consensus` (span 20 kb,
intron 3 kb, ≥ 2 reads multi-exon / ≥ 4 single-exon).

### Held-out core (R = 1,874; multi-CDS 1,301)

| | eval | genes | transcripts | locus | CDS coord-exact | CDS intron chain | missed | novel | splits | merges | F1 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| Helixer alone | — | 2,355 | 2,355 | 1,830 (97.7%) | 710 (37.9%) | 599 (46.0%) | 33 | 420 | 64 | 23 | 0.336 |
| GMB original | canonical | 2,520 | 3,936 | 1,850 (98.7%) | 695 (37.1%) | 589 (45.3%) | 20 | 436 | 88 | 22 | 0.316 |
| GMB original | all isoforms | | | | 751 (40.1%) | 646 (49.7%) | 19 | 429 | 342 | 253 | 0.338 |
| **GMB updated** | canonical | 2,496 | 3,960 | 1,852 (98.8%) | **733 (39.1%)** | **636 (48.9%)** | 18 | 432 | 63 | 23 | 0.335 |
| GMB updated | all isoforms | | | | 750 (40.0%) | 654 (50.3%) | 18 | 428 | 65 | 24 | 0.343 |
| GMB original + long reads | canonical | 2,550 | 4,844 | 1,852 (98.8%) | 719 (38.4%) | 625 (48.0%) | 16 | 455 | 81 | 21 | 0.325 |
| **GMB updated + long reads** | canonical | 2,530 | 4,743 | 1,852 (98.8%) | **754 (40.2%)** | **668 (51.3%)** | 17 | 443 | 63 | 20 | **0.342** |
| GMB updated + long reads | all isoforms | | | | 806 (43.0%) | 702 (54.0%) | 17 | 439 | 63 | 23 | 0.366 |

CDS-intron precision / recall (canonical): Helixer 0.572 / 0.661; original 0.575 / 0.654;
updated 0.594 / 0.680; updated + long reads **0.618 / 0.693**.

### Held-out accessory (R = 175)

| | genes | CDS-exact | missed | novel | splits / merges | F1 |
|---|---|---|---|---|---|---|
| Helixer alone | 433 | 16 | 25 | 224 | 29 / 2 | 0.053 |
| GMB original (canonical) | 476 | 19 | 15 | 231 | 33 / 2 | 0.058 |
| GMB updated (canonical) | 466 | 19 | 13 | 223 | 26 / 2 | 0.059 |
| GMB updated + long reads (canonical) | 475 | 20 | 14 | 227 | 27 / 2 | 0.062 |

### Development windows (tuned on; for completeness)

| core (R = 1,865), canonical | CDS-exact | intron chain | missed | splits / merges | F1 |
|---|---|---|---|---|---|
| Helixer alone | 719 | 557 | 44 | 69 / 18 | 0.343 |
| GMB original | 706 | 553 | 36 | 89 / 18 | 0.323 |
| GMB updated | 750 | 599 | 36 | 68 / 18 | 0.346 |
| GMB updated + long reads | 773 | 631 | 24 | 69 / 18 | 0.352 |

Accessory dev canonical CDS-exact: Helixer 25, original 26, updated 28, updated + long reads 27.

### What changed, gene by gene (held-out canonical)

| transition | gained exact | lost exact | partial → intron chain | missed → found |
|---|---|---|---|---|
| Helixer → GMB original | 11 | 23 | 1 | 26 |
| Helixer → GMB updated | 31 | 5 | 15 | 27 |
| Helixer → GMB updated + long reads | 62 | 14 | 21 | 27 |
| GMB original → GMB updated | 48 | 10 | 17 | 5 |

| coding models vs Helixer (held-out canonical) | changed | improved | regressed | uncertain | new loci (no reference CDS overlap) |
|---|---|---|---|---|---|
| GMB original | 103 | 10 | 16 | 77 | 248 (236) |
| GMB updated | 166 | 44 | 6 | 116 | 216 (207) |
| GMB updated + long reads | 323 | 80 | 19 | 224 | 257 (243) |

"Uncertain" = changed but no closer to or further from the reference. On aggregate these
raise intron precision and recall (above), so on balance they are structural improvements,
but individually they are unverified.

### Updated implementation: evidence ablation (canonical, core)

| | dev CDS-exact | held-out CDS-exact |
|---|---|---|
| Helixer only / + protein / + long reads only | 719 / 719 / 719 | 710 / 710 / 710 |
| + short reads | 721 | 711 |
| + short + long reads | 745 | 725 |
| + short reads + protein | 750 | 733 |
| all | 773 | 754 |

A transcript structure needs a second independent agreeing source, or protein support,
before it can matter. The gains come from short reads + protein, and from long reads +
short reads.

### Generalisation, trade-offs, cost

* Generalises: every direction measured on dev holds on held-out (canonical CDS-exact,
  intron chain, splits/merges, missed). Held-out gains are smaller than dev without long
  reads (+23 vs +31 over Helixer), similar with them (+44 vs +54).
* Trade-offs: canonical F1 without long reads is level with Helixer (0.335 vs 0.336). The
  exact-CDS gain is offset by ~140 extra genes, mostly transcript-only antisense/intergenic
  models present in both GMB versions. All-isoform CDS-exact is unchanged (750 vs 751): the
  update moves correct models into the canonical slot rather than adding new ones.
  Accessory gains are within noise.
* Cost: held-out build time 70 s → 64 s (8 windows), peak RSS 424 → 403 MB. Long reads add
  ~13%. Collapsing long reads took < 1 s per window.

---

## F. Prioritised recommendations

**Essential before fungal pipeline handover**

1. **Run the production-scale test with both configurations**: plain `--preset fungi` and
   `--preset fungi --config configs/fungi_evidence_integration.experimental.yaml`. Evaluate
   the *canonical* output against both the reference and Helixer, core and accessory
   separately.
2. **Supply collapsed long reads.** Run `gmb-longread-consensus` (fungal span/intron limits)
   upstream of GMB. It is the most precise and the only independent transcript evidence, and
   it is unusable as delivered.
3. **Do not ship the original canonical output as-is.** Without `prefer_build_primary` it
   is less accurate than Helixer on held-out windows. If the overlay is not adopted, at least
   set `canonical_selection.prefer_build_primary: true`.

**Worth addressing during production-scale testing**

4. Decide the treatment of transcript-only antisense/intergenic loci: a Psauron gate,
   requiring protein-compatible or multi-technology support, or a low-confidence class.
5. Fix the fungal protein-validation settings before enabling them: with `diamond_weight`
   0.7 and `min_score` 0.7 a DIAMOND hit is mandatory. Decide whether validation should gate
   new loci, which it currently does not, rather than prune alternates.
6. Decide how many UTR-only alternates to ship (half of remaining alternates).
7. Canonical-selection housekeeping: ORF completeness and protein-alignment class are read
   only from protein validation, so without it those tiers are inert and every gene is
   "low confidence".

**Longer-term algorithm development**

8. Protein-to-genome candidate models (e.g. miniprot-style) with a repeat annotation, for the
   missed small genes and CDS boundaries proteins define exactly.
9. Multi-ORF inference on read-through transcripts, gated by coding potential.
10. Replace the additive per-source score with junction- and CDS-level evidence throughout;
    the junction rule is a targeted patch on top of a backbone-dominant score.

### Recommendation

**Carry both as configurations in the first production-scale test, with the experimental
overlay as the intended default and long reads collapsed upstream.** It is the only
configuration in which the evidence measurably improves Helixer's gene models (held-out core
canonical CDS-exact +23, or +44 with long reads, intron precision and recall both up,
splits/merges at Helixer's level). It is opt-in and leaves the validated preset unchanged.
The evidence for it is regional (17 windows, 3.7 k reference genes) and one dataset, so it
should be promoted into `fungi.yaml` only after the genome-wide test confirms it. The
original fungi output should not be handed over without at least `prefer_build_primary`.
