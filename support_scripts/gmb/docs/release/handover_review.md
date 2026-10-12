# Final pre-handover review — evidence integration (fungi)

> **Naming note (2026-10-12).** The overlay this report calls
> `configs/fungi_evidence_integration.experimental.yaml` ("configuration B", "the overlay") is now
> the packaged preset **`fungi_latest_settings`** (`gmb/configs/fungi_latest_settings.yaml`,
> identical settings, byte-identical builds), run with `gmb-run --preset fungi_latest_settings`
> ([`../../RUN_FUNGI.md`](../../RUN_FUNGI.md)). Long reads no longer need collapsing by hand:
> `gmb-run --long-read-alignments` does it. The text below records the work as it was done.

Branch `gmb/evidence-tuning`, 2026-10-09. Reviews the four opt-in evidence-integration settings
and the configurable protein-validation penalty (commit be54b99) before pipeline handover. Fixes
made during this review are listed in §9; nothing is pushed, merged or tagged.

Companion documents: [`../../HANDOVER.md`](../../HANDOVER.md) (how to run both fungal
configurations), [`production_scale_test_plan.md`](production_scale_test_plan.md) (the
whole-genome comparison), [`evidence_integration.md`](evidence_integration.md) (the
investigation). Every number below is reproduced from the original input data by
[`../../validation/z_tritici_evidence_integration/`](../../validation/z_tritici_evidence_integration/)
(`results/`).

---

## 1. Implementation review

| setting | parsed / validated | where it acts | default (A) | regression tests |
|---|---|---|---|---|
| `scoring.primary_selection: junction_supported` (+ `junction_primary_min_cds_fraction`) | `ScoringConfig`; invalid value → `ValueError` | `select_isoforms` → `_promote_junction_supported`, after gene grouping | `score` (no change) | 8 unit + fixture build |
| `scoring.isoform_cds_overlap: drop` | `ScoringConfig`; bare YAML `off` accepted | `select_isoforms` gene assignment | `off` | 6 unit + fixture build |
| `canonical_selection.prefer_build_primary: true` | `CanonicalSelectionConfig`; reaches `gmb-finalise` through the build's `resolved_config.yaml` | `_rank_key` tier after evidence-class breadth | `false` | 5 unit |
| `scoring.min_cds_bp: 150` | `ScoringConfig` | retention floor in `select_isoforms` (now reads the candidate's CDS; before this work it read a field nothing set and never applied) | `0` (= previous effective behaviour) | 2 unit + fixture build |
| `protein_validation.penalty` | `ProteinValidationConfig`; negative → `ValueError` | `score_model`, only under `policy: penalize` | `5.0` (the former hard-coded value) | 2 unit |

**Defaults unchanged, verified end to end.** On the bundled fixture, configuration A on this
branch and the baseline commit `feature/gene_model_builder` (e005a86) give byte-identical
`finalise/consensus.gff3`, canonical GFF3, `cdna.fa`, `cds.fa` and `prot.fa`.
`evidence_attribution.tsv` and `canonical_transcripts.tsv` differ only by appended columns
(`validation/smoke_fixture.sh`). On all 17 windows, current code on the original window inputs
reproduces the earlier configuration A and B builds byte for byte (34 of 34 `consensus.gff3`).

**Structures are valid.** Every finalised output of the reproduction (15 GMB arms × 17 windows
= 255) and of the smoke test passes FASTA and UTR QC. An independent checker
(`validation/check_handover.py`, written without GMB code) also finds unique IDs, intact
parent/child links, CDS and UTR inside exons, gene span = union of transcripts, one canonical
per gene, and sequences consistent with the GFF3. Both configurations are byte-for-byte
deterministic on a rerun. An installed wheel gives output identical to the editable install.

### Junction-supported primary selection

The rule replaces a backbone primary only when **all** of the following hold:

* the backbone model has ≥ 1 intron that no transcriptomic model observes exactly (junctions
  are pooled from every loaded short- and long-read model, including read-through models
  removed by the span filter, which makes the backbone *harder* to replace);
* the replacement is an admitted alternate of the same gene: same sequence and strand, having
  passed the retention gate and `min_alternate_score`; unstranded models never reach selection;
* it has fewer unobserved introns (for a transcript-derived model its own introns count as
  observed, so in practice: it is a transcript structure);
* it has a complete ORF (ATG and stop);
* its CDS is ≥ 90% of the backbone's and shares ≥ 50% of the coding bases of both models (the
  same coding locus — this excludes read-through ORFs);
* compatible minus incompatible spliced protein alignments is no lower than the backbone's.

The backbone model is kept as an alternate in every case (483 of 483 promotions on the 17
windows; 645 of 645 with long reads). Opposite-strand transcripts cannot replace it (test).
Chimeras are constrained twice: by the 20 kb span filter before candidacy and by the
coding-locus test.

What "transcriptomic support" means here: 19% of promotions (94 / 483) rest on one assembly
whose introns no other track confirms. Their outcomes are in proportion with the rest (13
improved, 2 regressed, 79 uncertain vs 94 / 9 / 380 overall). Most promotions are
reference-neutral.

| promotions (U_all, 17 windows) | n | improved | regressed | uncertain |
|---|---|---|---|---|
| ≥ 2 transcript sources | 83 | 12 | 0 | 71 |
| one source, introns confirmed by another track | 306 | 69 | 7 | 230 |
| one source, self-supported only | 94 | 13 | 2 | 79 |

"Improved / regressed": the promoted model's CDS is closer to / further from the reference than
the backbone's (exact > exact intron chain > overlap).

One sensitivity: the protein-compatibility tie-break depends on which redundant protein
alignments survive the protein filter (§8). Correcting the window subsets of the OrthoDB file
changed the primary at two loci in one development window (one reference-exact canonical lost,
773 → 772 with long reads). Held-out results were unaffected.

## 2. Canonical transcript selection

Remaining canonical overrides of the build primary, all 17 windows
(`results/canonical_overrides.txt`):

| configuration | overrides / genes | deciding tier | reference: gained / lost / neither |
|---|---|---|---|
| A | 372 / 6,089 | named-source count 340, CDS length 32 | 20 / 35 / 317 |
| B | **0** / 5,988 | — | — |
| A + long reads | 602 / 6,126 | evidence-class breadth 383, named-source count 179, CDS length 40 | 60 / 31 / 511 |
| B + long reads | 243 / 6,095 | evidence-class breadth 243 | 21 / 14 / 208 |
| B + protein validation | 14 / 5,988 | protein validation 14 | 2 / 1 / 11 |

`prefer_build_primary` removes every named-source-count override, the problem it was introduced
for. What can still override the build primary is, by design, the tiers above it. With long
reads, a long-read + short-read isoform (two evidence classes) beats a single-class primary;
this is close to neutral on held-out windows (9 gained, 7 lost) and slightly positive on
development windows. With protein validation, the protein tier decides 14 cases. Both are
scientifically defensible (independent technology, independent protein evidence), so the
order was not changed. The documented priority order matches the code
(`../canonical_selection.md`, corrected in this review: it previously omitted the
internal-stop and named-source tiers).

**With and without protein validation**, configuration B selects the identical canonical CDS
for 94.5% of genes. The differences come almost entirely from the build, not from canonical
selection. The fungal validation weights penalise every model without a DIAMOND hit by 5.0,
which pushes transcript alternates below `min_alternate_score`, so the junction rule has fewer
alternates to promote. Held-out core canonical exact CDS: B 733, B + protein validation 718.
Protein validation stays off in both configurations for the first run.

**Reporting fixed.** `selection_reason` is `LOW_CONFIDENCE_NO_PROTEIN_SUPPORT` for every gene
when protein validation is off, so it never said what decided. `canonical_transcripts.tsv` now
also carries `deciding_tier` and `canonical_is_build_primary`, and the summary JSON carries
`deciding_tier_counts` and `canonical_not_build_primary_by_tier`. One residual gap: ORF
completeness for canonical selection is read only from `protein_validation.tsv`, so without
protein validation the complete-ORF tier is inert. Changing that would change configuration A's
output, so it is documented, not altered.

## 3. Isoforms and transcript-only loci

**`isoform_cds_overlap: drop`.** On the held-out windows it removes all 674 CDS-disjoint
alternates of configuration A while keeping 1,711 alternates that share coding sequence (864 of
them differ from the primary only in UTR). All-isoform splits / merges go from 342 / 253 to
65 / 24 (Helixer 64 / 23). A backbone model is never discarded; a disjoint one becomes its own
gene (test). Gene bounds and parent/child links pass the independent checker. Every discarded
structure is now written to `build/excluded_candidates.tsv` with the primary it would have
joined, and counted in `summary.json` (`isoform_cds_overlap_excluded`). Before this review the
drop left no record.

**Transcript-only loci** are now identifiable without a reference.
`evidence_attribution.tsv` gains `gene_has_backbone_model`, `opposite_strand_cds_overlap`,
`orf_label`, `is_partial_5/3` and `internal_stop_count`. On held-out core windows (B,
canonical) there are 158 transcript-only genes. The flag "transcript-only and CDS overlapping
an opposite-strand CDS" selects 107, of which 95 are antisense to a reference gene. It misses
16 reference-antisense cases (no output CDS on the other strand); that is 89% precision and 86%
recall against the reference view. Of the 107, 100 have complete ORFs, **none** has a
compatible spliced protein alignment, and 14 have ≥ 2 transcript sources. Psauron (where run)
scores most of them < 0.5. They are probably incidental ORFs in antisense transcripts. No
filter was added; the production-scale plan reports them separately (§5 of the plan).

## 4. Long-read integration

Verified workflow (commands in `../../HANDOVER.md` §5):

1. **Stranded alignments upstream.** GMB never re-derives strand. The Ensembl anno Minimap2
   path sets transcript strand from the `ts` tag. The *Z. tritici* track's splice sites are
   canonical on the declared strand (99.1% of 3.2 M read introns genome-wide, release
   preparation), so its strand is correct.
2. **Collapse** with `gmb-longread-consensus --preset fungi` (new preset: span 20 kb and intron
   3 kb, matching the fungi build preset; read support ≥ 2 multi-exon / ≥ 4 single-exon). It
   reproduces the evaluated consensus track byte for byte. 50–100 k reads per Mb collapse to
   ~500–650 models per Mb.
3. **Preflight** with `--minimap2 <consensus>`: on the raw per-read track it WARNs
   `longread_collapsed` (12,771 models per Mb on one window); on the consensus it passes
   `splice_quality` (98.6% canonical).
4. **Build** with `--minimap2 <consensus>`. A per-read track now also prints a `WARNING` in
   `gmb-build` and records `longread_models_per_mb` in `summary.json`. Before this review only
   preflight warned. Neither refuses the track.

For the first production-scale run: supply `minimap2_consensus.gtf` from step 2 as `--minimap2`
in configurations A + long reads and B + long reads, never `*_minimap2_annotation.gtf`.

## 5. Verification of the reported results

**Reproduced from scratch.** `prepare_inputs.sh` regenerates the genome, reference and remapped
Helixer byte-identically from the original files, then subsets windows and collapses long
reads. `run_matrix.sh` then reruns all 16 arms in all 17 windows. **Every held-out figure in
`evidence_integration.md` reproduces exactly.** On development windows one figure moved
(B + long reads: 773 → 772 exact CDS, 631 → 630 intron chains), explained below.

**Matching criteria and denominators.** `annotation-qc pairwise-compare --evaluation-mode
protein_coding --reference-transcript-biotypes protein_coding --region <window>`. Exact CDS =
start, every splice site and stop identical (`cds_coordinate_exact`). Intron chain =
`cds_intron_chain_recovered` over reference genes with ≥ 2 CDS segments (1,301 held-out core).
Denominator R = reference protein-coding genes overlapping the window, pooled over windows
(held-out core 1,874, accessory 175). Only whole models inside a window are given to GMB, so
genes crossing window boundaries are partial or missed in every arm alike. F1 is pooled
`2·TP / (R + Q)` over one-to-one exact-CDS pairs.

**Evaluated outputs.** The canonical rows use
`finalise/canonical/consensus.canonical_annotated.gff3` with `--query-transcript-selection
canonical`, the finalised handover annotation. The all-isoform rows of the original evaluation
read `build/consensus.gff3`; it is byte-identical to `finalise/consensus.gff3` in every checked
window, and the reproduction evaluates the finalised file.

**Reconciliation.** The earlier "31 gained / 5 lost" transition counts pooled core and
accessory windows under a core heading. Core alone: 28 gained, 5 lost, net +23 = 733 − 710.
Accessory: 3 / 0. The model-change table had the same problem. Both are corrected in
`evidence_integration.md`.

**Development / held-out separation.** The windows were fixed in `windows.tsv` before any
result. Held-out windows were first built and evaluated after every setting and threshold was
fixed, and no code changed in response to them. The protein-validation, canonical-override and
promotion analyses in this review are diagnostic, not tuning.

**Input-subset artefact (found and corrected).** genBlastG reuses alignment IDs on different
chromosomes, and the window-subsetting script grouped rows by ID alone. So window OrthoDB files
carried rows from other chromosomes (sequences absent from the window genome; 50 and 38 rows in
the two windows counted), and 0–11 alignments per window were missing. Gene structures were identical in every
checked window and all held-out figures reproduce. The changes were protein attribution
(`ProteinEvidence`) of 1–26 transcripts per window, and the one dev locus above. The script now
keys rows by (sequence, ID).

**Reference agreement vs independent evidence.** The reference is a benchmark with known errors
and one isoform per gene. Evidence independent of it points the same way. On held-out core
windows, CDS-intron precision against the reference rises from 0.572 (Helixer) to 0.594 (B).
More independently, the share of canonical CDS introns confirmed by the collapsed long-read
track, which none of these runs used, rises from 67.6% (Helixer) and 68.8% (A) to 70.9% (B).
Neither test says whether the extra transcript-only loci are genes.

## 6. Reproducibility

`validation/z_tritici_evidence_integration/` (in the repository, outside the installed package)
holds everything needed: `windows.tsv` (with the selection rationale), `arms.tsv` (the arm
matrix), `prepare_inputs.sh`, `run_arm.sh`, `run_matrix.sh`, `make_results.sh`, the analysis
scripts, `diagnosis/` (the exploratory scripts behind the investigation), and `results/` (small
summary tables and `provenance.txt` with commit and tool versions). All paths come from
environment variables (`WORK`, `INPUT_DATA`, `REFERENCE_DATA`, `GMB_BIN`, `AQC`); no
machine-specific path remains. Large intermediate files go only to `$WORK`. The original
`diagnostics/gmb_evidence_tuning_20261009/` directory beside the repository is superseded and
not needed.

## 7. Scientific limitations

1. **Proteins support and veto; they never construct models.** A gene GMB misses but that a
   protein alignment defines exactly (5 of 44 Helixer-missed core genes on development windows)
   cannot be recovered. Gap-filling with protein-defined CDSs was measured, not implemented: 4
   reference-exact of 73 added models, 37 with no reference gene at all (likely
   transposon-derived). Needs a protein-to-genome model builder and a repeat annotation.
2. **One ORF per transcript; read-through assemblies.** Most missed small genes lie inside long
   read-through assemblies whose longest ORF is elsewhere. The span filter removes most such
   transcripts and the coding-locus test stops their ORFs replacing backbone models, but nothing
   recovers the second gene. Multi-ORF inference would need coding-potential support.
3. **Transcript-only antisense coding predictions**, about 4–5% of genes in both
   configurations, flagged but not filtered (§3).
4. **Canonical selection and isoform policy.** Evidence-class breadth can still override the
   build primary when long reads are present (243 / 6,095, near-neutral). The complete-ORF tier
   is inert without protein validation. Half of B's alternates are UTR-only variants. The
   reference has one isoform per gene, so isoform value is not measurable against it.
5. **Protein validation** with the fungi weights makes a DIAMOND hit mandatory (no-hit models
   score ≤ 0.3 < 0.7) and is net harmful combined with B (−15 exact CDS held-out). Weights and
   threshold need recalibration before production use.
6. **Junction rule.** Most promotions are reference-neutral (uncertain). It cannot help where no
   transcript structure exists, and its protein tie-break is sensitive to protein-filter
   tie-breaking (§8).
7. **Regional validation.** 17 windows, 38% of one genome, chosen by hand; accessory
   chromosomes contribute only 411 reference genes, so accessory effects are within noise.

## 8. Readiness decision

| question | decision |
|---|---|
| Is the current baseline (A) ready for production-scale testing? | **Yes.** Unchanged and byte-identical to the validated release candidate. One caution: its canonical set is slightly *less* accurate than Helixer alone on held-out windows (695 vs 710 core). |
| Is the experimental configuration (B) ready for production-scale testing? | **Yes**, as a second, opt-in configuration. Correct, optional, tested, deterministic, reproducible, and better than A and Helixer on held-out windows. |
| Is anything serious enough to block pipeline integration? | **No.** The remaining issues are policy questions (antisense loci, isoform counts, protein-validation calibration) or longer-term algorithm work, and none produces invalid output. |
| Before adopting B as the fungal default | the production-scale test must meet the pre-registered criteria (`production_scale_test_plan.md` §6). Decide how to treat transcript-only antisense loci. Decide whether long-read collapse becomes a pipeline requirement. |
| Before a stable tag | the above; the apicomplexan presets re-validated or the tag scoped to fungi (`../known_issues.md`); protein-validation weights recalibrated or documented as unsupported for fungi. |

## 9. Change summary (this review; on top of be54b99)

| file | change | why |
|---|---|---|
| `gmb/pipeline/scoring.py` | `select_isoforms(..., excluded=)` records structures refused by `isoform_cds_overlap: drop` | auditability |
| `gmb/pipeline/builder.py` | writes `excluded_candidates.tsv` and the `isoform_cds_overlap_excluded` counter; appends `orf_label`, `is_partial_5`, `is_partial_3`, `internal_stop_count`, `gene_has_backbone_model`, `opposite_strand_cds_overlap` to the attribution; `longread_density_warning` (build-time per-read warning, `longread_models_per_mb`) | identify transcript-only and antisense loci without a reference; safe long-read input |
| `gmb/pipeline/canonical_selection.py` | `deciding_tier` and `canonical_is_build_primary` columns; `deciding_tier_counts`, `n_canonical_not_build_primary`, `canonical_not_build_primary_by_tier`, `prefer_build_primary` in the summary | the reported reason hid what decided |
| `gmb/cli/finalise.py` | `__main__` guard | `python -m gmb.cli.finalise` exited 0 without running |
| `gmb/configs/longread_consensus/fungi.yaml` | new `gmb-longread-consensus` preset | documented, reproducible long-read collapse for fungi |
| `tests/test_evidence_integration.py` | +14 tests: strand protection, exclusion audit, deciding tier, class-breadth override, fixture locus context and exclusion sidecar, long-read preset and warning, every CLI under `python -m` | regression coverage |
| `tests/test_release_prep_regressions.py` | pinned attribution columns extended (appended only) | interface contract |
| `HANDOVER.md` | rewritten as the single run document: both configurations, inputs, verified commands, outputs | handover |
| `README.md`, `docs/quickstart.md` | pointers to the two configurations and long-read collapse | handover |
| `docs/input_contract.md`, `docs/longread_consensus.md` | long-read workflow, `fungi` preset, build-time warning | long-read integration |
| `docs/output_contract.md` | new attribution columns, `excluded_candidates.tsv`, canonical columns and counters | interface contract |
| `docs/canonical_selection.md` | remaining overrides, reporting, standalone-command configuration | accuracy |
| `docs/configuration.md` | exclusion audit; protein validation × junction rule | accuracy |
| `docs/known_issues.md` | new P3 entry (protein-filter tie order); interaction note; fixed items | accuracy |
| `docs/release/evidence_integration.md` | transition and model-change tables corrected (core / accessory separate); independent long-read check; input-artefact note | result verification |
| `docs/release/handover_review.md`, `docs/release/production_scale_test_plan.md` | new | deliverables |
| `validation/` | reproducible evaluation (`z_tritici_evidence_integration/`), `smoke_fixture.sh`, `check_handover.py` | reproducibility, smoke test |

Behaviour of configuration A is unchanged (byte-identical outputs; appended columns only).
Configuration B's gene models are unchanged by this review: it adds reporting and audit files
only.

**Tests.** GMB 903 passed, 3 skipped (was 889 / 3); transcript-strand 25 passed
(`PYTHONPATH=src/python`); clade skill 28 passed. Fixture smoke test (both configurations,
finalise, independent checks, determinism, equivalence with e005a86): passed, from the
editable install and from a wheel in a fresh environment. Lint: no new findings (4
pre-existing `ruff` findings in `builder.py` / `scoring.py`; black clean on changed files).
