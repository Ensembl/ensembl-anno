# Evidence-score primary selection (`primary_selection: evidence_score`)

An **opt-in** primary-selection policy. It chooses each gene's primary transcript among the
candidates GMB already has, using quantitative splice-junction evidence from STAR, per-intron
protein agreement and ORF integrity. It creates no new structure, splits or merges no gene, and
admits no extra isoform. No shipped preset enables it. Every other `primary_selection` value
behaves exactly as before, and its outputs are byte-identical.

The rule is "P4d" from the Phase 1.5 experiment (*Z. tritici* development windows). Its
parameters are the values frozen there; do not tune them on held-out data.

## Enabling it

```yaml
# overlay.yaml
scoring:
  primary_selection: evidence_score
```

```bash
gmb-build --preset fungi_latest_settings --config overlay.yaml ... \
    --star-junctions junctions_manifest.tsv --star-junction-roles construction
```

`gmb-run` and `run_gene_model_builder(..., star_junctions=..., star_junction_roles="construction")`
pass the same inputs through.

**Evidence roles are enforced, not defaulted.** The policy uses only manifest libraries whose
role is `construction`. Each of the following stops the run with an explicit error; there is
never a fallback to another policy:

| condition | `gmb-build` / `gmb-run` | `gmb-preflight` | Python API |
|---|---|---|---|
| no `--star-junctions` | exits before any work | `evidence_score_inputs` FAIL | `ValueError` before any stage |
| `--star-junction-roles` other than exactly `construction` (`holdout`, `all`, lists) | exits | FAIL | `ValueError` |
| a construction library that does not load (missing, zero-byte, malformed) | exits | FAIL | (build exits) |
| an index containing any non-construction library | `EvidenceSelectionError` | – | – |

Reference annotations and evaluation results are never inputs to GMB, so they cannot reach the
score.

## How it fits into selection

`select_isoforms` runs unchanged up to gene grouping, the existing `junction_supported`
promotion included. Then, for each gene of the locus:

1. **Candidate pool.** These are every CDS structure of the locus (identical CDS from several
   transcripts or sources form one structure) that:
   * shares at least `same_coding_locus_fraction` of its CDS with the current primary, in both
     directions, on the same strand;
   * has at least `min_candidate_cds_bp` of CDS (and at least `scoring.min_cds_bp`), and passes
     protein validation when its `drop` policy is on (every carrying transcript must pass, see
     [protein_qc.md](protein_qc.md); in-build validation is not recommended for fungi);
   * is carried by a transcript that would pass the build's GFF3 validation limits (merged exon
     length and, if enabled, transcript span);
   * is not used by any other gene of the locus (primary or isoform), so no CDS is duplicated
     across genes.

   The pool is **not** limited to the gene's retained isoforms. In the development windows,
   33 of 111 replacements chose a structure GMB had not kept as an isoform.
2. **Features** for the primary and each pool member are computed (below).
3. **The P4d rule** decides.
4. **If replaced:**
   * the new primary is the transcript carrying the chosen CDS, preferring one already retained
     in the gene, then an admitted structure, then any scored structure;
   * the old primary and the remaining isoforms stay only if they pass the normal isoform
     admission rules against the new primary (`max_isoforms_per_locus`, same gene, shared coding
     base under `isoform_cds_overlap`, `min_alternate_score`, distinct intron chain);
   * gene grouping is otherwise unchanged.

The junction index is loaded once and shared with the Phase 1 diagnostics
(`junction_evidence/`). The protein-intron index is built once and shared with
`protein_compatibility_counts`.

## Features

Per CDS intron, using the Phase 1 status of `JunctionAnnotator`, with
f(g) = min(1, g / `breadth_saturation_groups`) over **condition groups**, so replicate
libraries count once:

| status | junction term J |
|---|---|
| SUPPORTED, usage u (lower of donor and acceptor usage against competing junctions) | f(own) · u − f(best competitor) · (1 − u) |
| ALTERNATIVE_SITE_OBSERVED | −f(competitor) |
| OPPOSITE_STRAND | −f(opposite-strand junction) |
| UNOBSERVED_EXPRESSED_LOCUS | −`unobserved_expressed_weight` × median f of the candidate's own supported introns (0 if none) |
| UNOBSERVED_NO_LOCAL_EVIDENCE, UNOBSERVABLE_* (outside the aligner's intron range, filtered motif), MULTIMAP_ONLY | 0 |

* **Protein term P:**
  * P = `protein_intron_weight` · (2a / (a + c) − 1), where a counts spliced alignments that
    contain exactly this intron and c those with a different, overlapping intron;
  * unspliced or merely overlapping alignments are uninformative (P = 0);
  * the fraction bounds the term, so repeated alignments of one protein cannot add more than
    one alignment's sign. They can still tilt the ratio (see Limitations).
* **Motif:** −`noncanonical_penalty` when the genome motif is not GT/AG, GC/AG or AT/AC and J ≥ 0.
* **Integrity:** −`partial_orf_penalty` for a partial ORF; −`internal_stop_penalty` for an
  internal stop or a CDS length that is not a multiple of 3.
* **Priors:**
  * `backbone_prior` if a backbone source has this exact CDS;
  * `transcript_class_prior` if any transcript source has it. StringTie, Scallop and long reads
    are **one** class, counted once, and never add junction support.
* **Scores:** evidence = ΣJ + ΣP + Σmotif + integrity; total = evidence + priors.

## The P4d rule

The primary B is replaced by the best eligible candidate X only if all of the following hold:

1. **B is not a protected backbone.** A backbone structure with no intron of negative J and
   intact integrity is never replaced. Absence of RNA evidence therefore never removes a
   correct backbone.
2. **X has a complete ORF, no internal stop and an intact frame.**
3. **X keeps B's strong introns.** It does not drop any intron of B with J ≥
   `keep_supported_intron_J`.
4. **X improves by the margin on both scores.** total(X) − total(B) ≥ `replacement_margin` **and**
   evidence(X) − evidence(B) ≥ `replacement_margin`. A prior may protect B but can never be the
   reason to switch.

* **Choosing X:** the eligible candidate with the highest total, then the highest evidence.
* **Ties:** the earlier candidate in locus order. The locus order is that of the candidate exon
  table, which is deterministic for given inputs.

**Deliberately not part of the rule** (rejected in Phase 1.5):
* a penalty for junctions inside a CDS exon;
* a promoting prior;
* junction bonuses added to the source-weight score;
* a CDS-length floor relative to the old primary;
* dropping novel genes.

## Outputs (only in this mode)

* **`consensus.gff3`:** promoted primaries carry `SelectionReason=evidence_score`. No other
  attribute changes meaning.
* **`evidence_selection.tsv`:** one row per gene, with these columns:
  * identity: `gene_id`, `primary_transcript_id`, `policy_version`, `chrom`, `strand`,
    `outcome` (`replaced`, `protected_backbone`, `kept_below_margin`, `kept_no_eligible`,
    `no_cds_primary`), `protected_backbone`;
  * candidates: `old_candidate` / `new_candidate` (input transcript IDs), `old_structure` /
    `new_structure` (1-based CDS segments), `old_sources` / `new_sources`;
  * score components: `old_*` and `new_*` for J, P, motif, integrity, prior_backbone,
    prior_transcript, evidence and total;
  * the change: `delta_total`, `delta_evidence`, `margin`, `protected_introns` (B's introns that
    no replacement may drop), `gained_introns`, `lost_introns`, `old_cds_bp` / `new_cds_bp`,
    `old_protein_net` / `new_protein_net`;
  * the pool: `n_pool`, `n_eligible`, `n_excluded_*` (other gene, invalid transcript, stop or
    frame, partial ORF, supported intron), plus `best_alternative` and its deltas when nothing
    was replaced;
  * how the replacement was carried out: `representative`, `isoforms_dropped`;
  * `concerns`: reviewer flags that never change the decision (`cds_shortened_over_10pct`,
    `start_codon_changed`, `stop_codon_changed`, `protein_compatibility_worse`,
    `new_primary_was_not_a_retained_isoform`, `old_primary_not_kept_as_isoform`).
* **`evidence_selection_summary.json`:** the parameters, policy version, junction manifest path
  and SHA-256, every construction library with its group and SHA-256, index statistics,
  outcome counts and timing.
* **`resolved_config.yaml`:** contains `scoring.evidence_score`. The run manifest records the
  software revision.

A changed build primary is not necessarily the finalised canonical transcript: `gmb-finalise`
canonical selection runs afterwards and is unchanged.

## Secondary-ORF candidates (`evidence_score.secondary_orfs`, opt-in)

GMB annotates one ORF per transcript, its longest. A read-through StringTie, Scallop or
long-read model that spans two genes therefore never offers the second gene's correctly
spliced ORF as a candidate. With

```yaml
scoring:
  primary_selection: evidence_score
  evidence_score:
    secondary_orfs: true
```

GMB also generates those ORFs (`gmb.pipeline.secondary_orfs`) and offers them to the P4d
candidate pool of the **existing** gene they share a coding locus with. No preset enables it,
and it is refused with any other `primary_selection`. The rules were frozen in the Phase 3A
experiment (A_ev) and are not parameters:

* **Discovery** (after ORF annotation). Transcript-evidence models only, the same sources as
  the transcript-class prior. Every complete, maximal ATG..stop ORF of at least
  `orf.min_codons` codons and `max(min_candidate_cds_bp, scoring.min_cds_bp)` bp, lying wholly
  upstream or downstream of the transcript's own ORF on the cDNA. ORFs overlapping the own ORF
  in another frame are not generated. The projected CDS must re-read identically from the
  genome, start with ATG, end with a stop, contain no internal stop and keep a continuous
  frame. Structures are deduplicated across parent transcripts (all parents are kept); a
  structure that already is a candidate's CDS is not generated.
* **Admission (A_ev).** Only structures sharing a coding locus with an existing candidate of
  their locus are evaluated, since no other structure can enter a gene's pool. A structure is
  admitted when at least 50% of its CDS is covered by same-strand protein-alignment exons, or
  when it is spliced and every CDS intron is a construction junction.
* **Selection.** Admitted structures never enter clustering, gene grouping or isoform
  admission. In a gene's pool they follow the existing candidates, under the pool's own rules
  (shared coding locus, CDS floor, not another gene's structure, GFF3 limits), and the
  unchanged P4d rule decides. No new genes are created, and a generated structure that is not
  chosen as primary is never kept as an isoform.
* **Representation.** A chosen secondary ORF becomes the gene's primary transcript with exons
  equal to its CDS segments: no UTR of the read-through parent is carried into a neighbouring
  locus. It keeps its parents' transcript evidence class (`Evidence=StringTie,...`) and is
  marked `CandidateOrigin=secondary_orf;SecondaryOrfParents=<input transcript IDs>`. Its
  junction and protein compatibility counts are computed as for any candidate.

**Input requirement.** The protein clause uses the protein alignments as supplied. Prepare the
protein set **without the target species' own proteins**, otherwise the species' annotation
is aligned back to itself and counted as independent support. GMB cannot recognise them. The
requirement is restated in `secondary_orfs_summary.json`.

**Outputs.**
* `secondary_orfs.tsv`: one row per generated structure. Columns: `candidate_id`, CDS,
  `status` (`admitted`, `not_admitted`, `no_shared_coding_locus`), `admitted_by`,
  `protein_cover`, `protein_alignments`, `introns_junction_supported`, `parents`,
  `parent_sources`, `relations`, `loci`, `eligible_genes` / `eligible_primaries`, the P4d
  outcome at each of those genes, `selected` / `selected_gene`, and, for a selected structure,
  `abandoned_old_cds_bp` (coding bases of the replaced primary outside the new CDS) and
  `abandoned_bp_in_other_genes` (of those, the bases still coded by another output gene).
  The abandoned columns are reported for evaluation only; they never affect selection.
* `secondary_orfs_summary.json`: discovery counters, structures by status, admission clauses,
  eligible and selected counts, abandoned totals and timing.
* `evidence_selection.tsv`: a replacement by a secondary ORF has
  `new_candidate_origin=secondary_orf` and `representative=secondary_orf`.

With the switch off (the default) nothing above runs, and every output, including
`resolved_config.yaml`, is identical to a build without the key.

### Units of filter-removed transcripts (`evidence_score.filtered_transcript_orfs`, experimental)

The transcriptomic filter removes transcript models with an intron longer than
`transcriptomic_filter.max_intron_length` or a span longer than `max_transcript_length` (fungi:
3 kb and 20 kb). Most removed models are read-through assemblies joining adjacent genes. Their
coding units never reach secondary-ORF discovery, which scans only the surviving models. This
experimental switch, a GMB Phase 4 proof of concept, recovers them:

```yaml
scoring:
  primary_selection: evidence_score
  evidence_score:
    secondary_orfs: true            # required
    filtered_transcript_orfs: true
```

* **Units.** Each removed model is annotated like any candidate and scanned with the same
  frozen discovery rule. Its own complete GMB ORF is a unit too, since GMB never saw it. A unit
  is only a CDS; no read-through exon or UTR re-enters clustering.
* **Safeguards.** Each one is a `status` in `secondary_orfs.tsv`:
  * The unit must be a complete, genome-validated ORF.
  * It must itself pass the filter its parent failed (`exceeds_filter`).
  * It must not equal an existing candidate or a secondary ORF of a kept model (counted, not
    generated).
  * It must share a coding locus with existing CDS structures of exactly one locus. Two or more
    loci are rejected as a possible fusion (`multi_locus`). With no existing structure it is
    `no_shared_coding_locus`: **no new gene is ever created**.
  * It must pass the unchanged A_ev admission (protein cover or junction-observed introns).
  * At most 10 admitted units per locus; the units carried by the most removed models are kept
    (`locus_cap`).
* **Selection** is unchanged: an admitted unit joins the existing gene's P4d pool after every
  other candidate and becomes primary only by the frozen P4d rule. P4d never replaces a
  protected backbone primary, so most units can only replace transcript-derived or contradicted
  primaries.
* **Outputs.**
  * Units appear in `secondary_orfs.tsv` with `relations` `removed_own`, `removed_upstream` or
    `removed_downstream`, and candidate IDs `<parent>.ftorfN` / `<parent>.ftownN`.
  * `secondary_orfs_summary.json` gains a `filtered_transcript_orfs` block: discovery,
    admission, status counts, selected count and parent filter reasons.
  * A chosen unit is marked `CandidateOrigin=filtered_transcript_orf`. The same value appears as
    `new_candidate_origin` in `evidence_selection.tsv`.
* **Evaluation status.** Rules were frozen on the Z. tritici development windows and U. maydis
  chromosome 1, where the result was neutral: +1/−1 and 0 exact reference CDS. Data held out
  from development, run once:
  * *A. fumigatus* pilot: +28/−0;
  * *U. maydis* outside chromosome 1: +16/−0;
  * Z. tritici held-out windows: +2/−0;
  * *B. cinerea*: +5/−1, where few models are filter-removed.

  Details are in the Phase 4 workspace report (`PHASE4_RECONSTRUCTION_DESIGN.md`). Not for
  production use without further validation.

## Limitations

* **Assemblies are not construction-only.** The StringTie and Scallop candidate tracks are
  assembled from all libraries. Only the junction index is construction-only, so held-out
  libraries can still shape the candidate pool, and held-out replication is corroborating
  rather than independent.
* **Counts are per mate.** Junction reads are counted per mate in the Z. tritici batch. The
  rule uses group breadth and usage, not raw counts.
* **Protein redundancy.** Protein alignments are counted per alignment; redundant alignments of
  one protein can tilt the per-intron ratio.
* **What it cannot fix.** Start and stop codons are judged only through ORF integrity: junctions
  cannot decide between start codons. Genes whose correct structure is not among the
  candidates, or that are wrongly joined or split, need candidate generation or locus work.
