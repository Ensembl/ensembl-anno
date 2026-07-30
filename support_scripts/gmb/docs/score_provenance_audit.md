# Score-provenance / double-counting audit (canonical transcript selection)

Scope: does every component that influences canonical-transcript selection
carry independent biological information, or is some of it counting the
same underlying evidence twice? This audit walks every component in
`gmb.pipeline.scoring.score_model()` (the model-building/inclusion stage --
produces `gmb_score`) and `gmb.pipeline.canonical_selection.score_transcript()`
(the post-build canonical-choice stage), states where each one's
information ultimately comes from, and records the one MANDATORY change
this audit produced plus a real chr1 before/after comparison.

## 1. Where every component's information comes from

| Component | Stage | Biological basis | Also counted elsewhere? |
| :--- | :--- | :--- | :--- |
| `weights.backbone/scallop/stringtie/minimap2` per named source | `scoring.score_model` (bakes into `gmb_score`) | Named software source presence | Re-derived from the same `evidence_sources` string in canonical_selection's class/flag components -- **deliberate**: `gmb_score` blends these at build time to decide which models are even *retained*; canonical_selection re-reads the same raw string afterwards to decide, among already-retained isoforms, which is *canonical*. Different questions, same raw input -- not a bug. |
| `multi_source_bonus * (n_sources - 1)` | `scoring.score_model` | Raw **named-source** count breadth | Yes -- counted again as `gmb_score_component` (30% weight) in canonical_selection's `evidence_subtotal`, but ALSO the audit's mandated fix (see §2) means canonical_selection's *own* breadth measure (`n_evidence_classes`) no longer uses raw source count. **Known residual inconsistency, deliberately NOT changed here** (out of scope): `scoring.py` governs which models are *built and retained at all* -- changing its multi-source-bonus semantics is a materially larger-blast-radius change (affects the actual candidate set, not just which isoform is marked canonical) than this task's remit. Flagged for a future, separately-scoped pass. |
| `protein_overlap_bonus` | `scoring.score_model` | Binary: any protein alignment overlap | Not directly re-scored; canonical_selection's `protein_alignment_support_weight` flag checks the same `protein_alignment` evidence class, at a different (post-build) point. Deliberate, not redundant in effect (bonus already baked into the `gmb_score` a transcript carries into canonical_selection). |
| `protein_coding_score` bonus/penalty (`val_cfg.policy in (penalize/bonus)`) | `scoring.score_model`, using **top-level** `ProteinValidationConfig` weights | DIAMOND+psauron combined score | **Yes, deliberately.** canonical_selection's `protein_validation_subtotal` re-scores the SAME DIAMOND/psauron result with **independent, separately-tunable weights** (`CanonicalProteinValidationWeights`), because canonical selection needs finer-grained, re-tunable scoring of protein plausibility than the model-retention stage's single scalar bonus/penalty. This is intentional layering, not accidental double counting: the top-level score decides *inclusion*, the canonical-selection score decides *ranking among included isoforms*, and the latter is a strictly more detailed view of the same underlying DIAMOND/psauron numbers (identity, coverage, hit presence individually), not merely a re-application of the same scalar. |
| `noncanonical_splice_penalty` | `scoring.score_model` (needs `genome`) | Splice-site canonicality | Not available at all in canonical_selection (no genome loaded there -- see §3). Baked into `gmb_score` only. |
| `psauron_component`, `diamond_identity_component` | canonical_selection | Protein plausibility (strength) | No duplication -- these are canonical_selection-only, weighted independently from the top-level score above. |
| `diamond_balanced_coverage_component` (was: separate `qcov_component`/`scov_component`) | canonical_selection | Protein plausibility (strength) | **Changed by this audit** -- see §2. |
| `n_evidence_classes` / `class_breadth_component` (was: raw `n_sources`) | canonical_selection | Evidence *breadth* (how many independent kinds of evidence exist) | **Changed by this audit (mandatory) -- see §2.** Distinct from `protein_validation_subtotal` above: one asks "how good is the protein evidence", this asks "how many independent *kinds* of evidence exist at all" (including `protein_validation` itself as one of five classes, per this task's Part 3 contract) -- not double counting, since presence and strength are different questions. |
| `transcriptomic_support` / `protein_alignment_support` / `longread_support` / `backbone_support` flags | canonical_selection | Presence of a specific evidence class | Individually redundant *with* `n_evidence_classes` in the sense that a transcript with 3 classes gets both the breadth score AND each matching flag's bonus -- **pre-existing, unchanged design**: breadth rewards "how many kinds", the flags let specific kinds (e.g. protein_alignment) be weighted more heavily than others if desired. Not touched by this audit (task Part 4 scope is the mandatory breadth-measure change, not a rebalancing of these weights). |
| `has_complete_orf`, `has_internal_stop`, `is_partial` -> `structure_subtotal` | canonical_selection | ORF/CDS structural validity | The **rank-key's** tier-1 gate (`not has_complete_orf`, `has_internal_stop`) is a hard lexicographic decision using the same two booleans `structure_subtotal` scores continuously -- deliberate: the continuous score explains "how much better" for `total_score`/confidence-gap purposes, the boolean tuple is what actually decides ties (see the module's own long-standing documented rationale). Not new to this audit. |
| `domain_subtotal` | canonical_selection | Domain/Pfam/InterPro evidence | Disabled by default (`domains.enabled=False`), all weights 0.0 -- inert on chr1 (no domain adapter has run outside the InterPro *resolver*, which is a separate, non-domain_subtotal code path). Out of scope for this audit. |

## 2. The mandatory change this audit produced

Decision (from the task this audit was written for): **canonical evidence
breadth must use biological evidence classes, not counts of named software
sources.** Scallop and StringTie are two assemblers customarily run over
the SAME short-read libraries; counting them as two independent sources
overstated support relative to one genuinely independent
`short_read_transcriptomic` class.

Implemented in `gmb.pipeline.canonical_selection`:

* `_rank_key()` tier 3 is now `-n_evidence_classes` (was tier 3
  `-n_independent_sources`, a raw named-source count). The raw count is
  **not discarded** -- it is demoted to tier 4, where it still breaks ties
  the class-level measure leaves open (e.g. Scallop+StringTie, one class,
  vs Scallop alone, still the same class -- the raw count of 2 vs 1 still
  matters as a *more granular*, not *more independent*, tie-break).
* `evidence_subtotal`'s breadth term (`independent_source_weight`) now
  scales with `n_evidence_classes` instead of a raw named-source count
  (`class_breadth_component`, capped at `independent_source_cap` classes,
  same as before).
* The evidence-class vocabulary/mapping itself was centralised into
  `gmb.pipeline.canonical_evidence` (previously duplicated, and already
  diverged, between `canonical_selection.py` and `interpro_review.py`), and
  now includes `protein_validation` as a fifth class (previously excluded
  entirely from `interpro_review`'s copy) -- required by this task's Part 3
  data contract.

A second, related change (required by this task's Part 5 hierarchy, not
merely this audit's Part 4 mandate) replaced the separately-weighted
`diamond_query_coverage_weight`/`diamond_target_coverage_weight` terms in
`protein_validation_subtotal` with one `diamond_balanced_coverage_weight`
term using `min(qcov, scov)/100` (`gmb.pipeline.canonical_evidence.balanced_coverage`,
already used by `interpro_review.py`'s ambiguity detection). The default
weight (0.30) equals the sum of the two weights it replaces, so the
subtotal's achievable range is unchanged.

**Not adopted** (reported, not silently implemented): rebalancing
`scoring.py`'s `multi_source_bonus` to also use evidence classes. That
governs which models are *built and retained* at all, a materially
larger-blast-radius change than this task's canonical-selection-only
remit -- left as a known residual inconsistency for a future, separately
scoped pass (see the table above).

## 3. Real chr1 comparison (before vs after, existing data only)

Re-ran ONLY `gmb.pipeline.canonical_selection` (pure Python, no DIAMOND/
psauron/InterProScan invocation) against the existing, already-computed
chr1 build at
`apicomplexa_data/gmb_runs/20260729_141256_chr1_interpro_contract/`
(`consensus.gff3` + `evidence_attribution.tsv` + `protein_validation.tsv`,
none regenerated), with `scoring.backbone_label: Tiberius` (matching how
that build was actually produced -- confirmed via
`evidence_attribution.tsv`'s `evidence_sources` column, which is all
`Tiberius`/`Scallop`/`StringTie` combinations, no `Helixer`).

* **135 genes total, 3 multi-isoform** (the only genes where a tier below
  "only isoform" can matter at all).
* **1 gene's initial canonical choice changed: `PFAL_00010`.**
  * Old winner: `PFAL_00010.t2` (`BEST_DIAMOND_COVERAGE`) -- driven by a
    fractional query-coverage difference (83.8 vs 83.5 vs 83.2 across the
    3 isoforms) while **target coverage was identical (82.9) across all
    three**. Weighting query and target coverage independently let that
    query-side noise decide a tie that the target side said was a dead
    heat.
  * New winner: `PFAL_00010.t1` (`LONGEST_VALID_CDS_TIEBREAK`) -- balanced
    coverage (`min(qcov, scov)`) is now **identical (0.829) across all
    three isoforms**, since target coverage was already tied; the tie
    correctly falls through every subsequent tier (evidence classes tied,
    named-source count tied, `gmb_score` tied at 5.2) down to CDS length,
    where `t1` (3213 bp) is longest.
  * The other two multi-isoform genes (`PFAL_00118`, `PFAL_00119`) did
    **not** change winner -- both are decided at tier 2
    (`protein_validation_subtotal`, via `BEST_PSAURON` /
    `BEST_DIAMOND_COVERAGE` respectively) with real, non-tied differences,
    so the tier-3/4 breadth-measure change never applies to them.
* **Tier 3/4 (evidence-class breadth / named-source count) decided zero
  genes on this chr1 sample** -- with only 3 multi-isoform genes and both
  non-`PFAL_00010` cases resolved at tier 2, there is no real-data instance
  here where the evidence-class-vs-named-source distinction itself flips a
  winner. The change is still correct and independently covered by unit
  tests (`tests/test_canonical_selection.py::TestDeterministicTieBreak`)
  that construct the tied-tier-2 scenario the tiny chr1 multi-isoform
  sample does not happen to contain.

### Important follow-up flagged by this comparison (not resolved here)

`PFAL_00010.t1` -- the new winner -- was **not** one of the two candidates
previously submitted for InterPro review (`max_candidates_per_gene: 2`
selected only `t2`/`t3` under the old ranking), so the existing
InterProScan output (6 proteins, 3 genes) has **no domain evidence at all**
for `t1`'s protein sequence (checksum `ae6bd830...`, confirmed absent from
`interpro_review_candidates.faa`, whose 3 `PFAL_00010`-related sequences
are `t2`'s and `t3`'s only... actually only `t2`+`t3` were submitted, `t1`
never was). This is exactly the kind of case Part 15's chr1 re-validation
must examine explicitly with the new canonical_selection output rather
than assuming the existing InterPro result still applies unchanged -- see
that section of the project report for how it was handled (re-running
`interpro_review`'s preparation stage against the new ranking is cheap and
does not require a new InterProScan invocation unless the new top
candidates are proven to need domain evidence the existing scan does not
cover).
