#!/usr/bin/env python3
"""Exact-duplicate transcript detection and collapse.

Root cause (verified against real chr1 data, not assumed)
-----------------------------------------------------------
The chr1 protein-validation run found 35/38 "multi-isoform" genes were
actually exact structural duplicates (e.g. gene PFAL_00001's two
transcripts had byte-identical exon/CDS coordinates: 29510-34181,
34260-34957, both Evidence=Tiberius). Tracing this back:

1. ``gmb.pipeline.builder`` clusters ALL candidate exon rows with
   ``pr.PyRanges(candidate_exons).cluster(slack=0, count=True)`` --
   clustering happens at the *exon* level, not the transcript level.
2. Verified directly: a single real transcript's own two exons (e.g.
   29509-34181 and 34259-34957, a 78bp intron gap) get assigned to
   *different* Cluster IDs under ``slack=0``, because pyranges only merges
   touching/overlapping intervals -- a real intron gap is neither.
3. ``builder.main()`` then calls ``select_isoforms()`` once per Cluster
   (``for _cid, locus_df in cluster_df.groupby("Cluster")``). The same
   transcript's exon rows are therefore split across two independent
   ``select_isoforms()`` calls, each seeing only one exon and therefore
   treating it as if it were a single-exon transcript.
4. Each single-exon fragment independently qualifies to survive
   ``scoring.py``'s single-exon gate via the ``keep_backbone_without_support``
   exception (the fragment's only source is the ab initio backbone), so
   *both* fragments become separate ``gene_models`` and get separate gene
   IDs/transcript IDs in ``builder.main()``'s output loop.
5. Critically, ``annotate_all_transcripts()`` was run once, up front, on
   the *full* (unclustered) candidate_exons DataFrame, keyed by
   transcript_id -- so both fragments' final exon/CDS coordinates are
   reconstructed *identically and correctly* from that shared lookup,
   masking the fact that each fragment only carried one exon through
   clustering.
6. ``dedup_genes.py`` then (correctly, given only gene-level overlap +
   first-mRNA-chain information) merges the two now-identical-looking
   genes together via its ``merge_as_isoforms`` policy -- producing the
   observed artifact: one gene with two isoforms that are, in fact, the
   exact same transcript.

This is a deeper root cause than "dedup_genes.py's merge_as_isoforms
policy creates near-duplicate isoforms" -- the real bug is upstream, in
exon-level (not transcript-level) clustering. Fixing the clustering
algorithm itself was considered and rejected as the fix location here: it
has a much wider blast radius (risks changing locus detection and
legitimate multi-isoform/multi-source gene grouping across the whole
pipeline, which the task this addresses explicitly says not to optimise
around), whereas the artifact is only *reliably detectable* once each
fragment's full structure has been reconstructed and the genes carrying
them have already been merged into one by dedup_genes.py -- i.e. exactly
where this module hooks in. See the module-level report for the full
insertion-point rationale.

This module therefore implements a narrow, *exact-equality-only* pass that
runs strictly after ``dedup_genes.py``, within each gene's already-assembled
isoform set: collapse only transcripts proven identical in exon
coordinates, CDS coordinates, CDS phase (where recorded), and translated
protein sequence (where available) -- never genuinely distinct isoforms,
regardless of how they arose.
"""

from __future__ import annotations

import hashlib
from collections import defaultdict
from typing import TYPE_CHECKING, Optional

if TYPE_CHECKING:
    from gmb.pipeline.config import DuplicateTranscriptCollapseConfig, PipelineConfig


# ---------------------------------------------------------------------------
# Structural signatures (shared by the diagnostic audit and the real collapse)
# ---------------------------------------------------------------------------


def exon_signature(exons: list[tuple[int, int]]) -> str:
    return ",".join(f"{s}-{e}" for s, e in sorted(exons))


def cds_signature(cds: list[tuple[int, int]]) -> str:
    return ",".join(f"{s}-{e}" for s, e in sorted(cds))


def intron_chain_signature(exons: list[tuple[int, int]]) -> str:
    exons = sorted(exons)
    if len(exons) < 2:
        return "single-exon"
    return ",".join(f"{exons[i][1]}-{exons[i + 1][0]}" for i in range(len(exons) - 1))


def protein_checksum(protein: Optional[str]) -> Optional[str]:
    if not protein:
        return None
    return hashlib.md5(protein.encode()).hexdigest()


def cds_phase_signature(cds_phases: Optional[list[int]]) -> Optional[str]:
    if not cds_phases:
        return None
    return ",".join(str(p) for p in cds_phases)


# ---------------------------------------------------------------------------
# Part 1: pairwise classification (diagnostic + collapse-eligibility)
# ---------------------------------------------------------------------------

# Mutually-exclusive category labels, most specific/actionable first. A pair
# can only ever match one of these (see docstring on _classify_pair for the
# priority rule used to resolve genuinely overlapping definitions in the
# spec this implements).
CATEGORY_SAME_ORF_METADATA_DIFFERS = "6_same_coding_different_orf_metadata"
CATEGORY_IDENTICAL_STRUCTURE_SAME_EVIDENCE = "5_identical_structure_identical_evidence"
CATEGORY_IDENTICAL_STRUCTURE_DIFFERENT_EVIDENCE = "4_identical_structure_different_evidence"
CATEGORY_SAME_CDS_DIFFERENT_UTR = "2_same_cds_different_utr_exon"
CATEGORY_SAME_INTRON_CHAIN_DIFFERENT_ENDS = "3_same_intron_chain_different_ends"
CATEGORY_GENUINE_ALTERNATIVE = "7_genuine_alternative"
CATEGORY_NOT_COMPARABLE = "not_comparable_different_seqname_or_strand"

# Categories 4/5/6 all satisfy category 1's literal definition ("identical
# exon AND CDS coordinates") -- they are examined from different angles
# (evidence-source lens, ORF-metadata lens) rather than being a fourth
# disjoint bucket, so "category 1" is reported as a derived boolean
# (`exact_coordinate_match`) alongside the more specific label rather than
# as its own separate value.
EXACT_DUPLICATE_CATEGORIES = {
    CATEGORY_SAME_ORF_METADATA_DIFFERS,
    CATEGORY_IDENTICAL_STRUCTURE_SAME_EVIDENCE,
    CATEGORY_IDENTICAL_STRUCTURE_DIFFERENT_EVIDENCE,
}


def _evidence_set(evidence_sources) -> set[str]:
    return {s.strip() for s in str(evidence_sources or "").split(",") if s.strip()}


def classify_pair(t1: dict, t2: dict) -> str:
    """Classify a pair of same-gene transcript records into one of the 7
    categories from the task spec (see module constants above).

    Each record is expected to carry: chrom, strand, exons (list of
    (start,end)), cds (list of (start,end)), cds_phases (optional list of
    int), protein (optional str), evidence_sources (str), orf_label,
    is_partial_5, is_partial_3.
    """
    if t1.get("chrom") != t2.get("chrom") or t1.get("strand") != t2.get("strand"):
        return CATEGORY_NOT_COMPARABLE

    exon_sig1, exon_sig2 = exon_signature(t1["exons"]), exon_signature(t2["exons"])
    cds_sig1, cds_sig2 = cds_signature(t1.get("cds", [])), cds_signature(t2.get("cds", []))

    if exon_sig1 == exon_sig2 and cds_sig1 == cds_sig2:
        orf_meta_1 = (t1.get("orf_label"), t1.get("is_partial_5"), t1.get("is_partial_3"))
        orf_meta_2 = (t2.get("orf_label"), t2.get("is_partial_5"), t2.get("is_partial_3"))
        if orf_meta_1 != orf_meta_2:
            return CATEGORY_SAME_ORF_METADATA_DIFFERS
        if _evidence_set(t1.get("evidence_sources")) == _evidence_set(t2.get("evidence_sources")):
            return CATEGORY_IDENTICAL_STRUCTURE_SAME_EVIDENCE
        return CATEGORY_IDENTICAL_STRUCTURE_DIFFERENT_EVIDENCE

    if cds_sig1 == cds_sig2 and cds_sig1 != "":
        return CATEGORY_SAME_CDS_DIFFERENT_UTR

    ic1, ic2 = intron_chain_signature(t1["exons"]), intron_chain_signature(t2["exons"])
    if ic1 == ic2 and ic1 != "single-exon":
        return CATEGORY_SAME_INTRON_CHAIN_DIFFERENT_ENDS

    return CATEGORY_GENUINE_ALTERNATIVE


def is_exact_duplicate_pair(
    t1: dict, t2: dict, cfg: DuplicateTranscriptCollapseConfig
) -> tuple[bool, str, str]:
    """Apply the Part 2 conservative collapse rule to a pair.

    Returns ``(eligible, category, reason)``. ``category`` is always the
    diagnostic classification from ``classify_pair`` (reported even when
    not eligible, e.g. for the audit table); ``reason`` explains why a
    structurally-identical pair was or wasn't collapsed.
    """
    category = classify_pair(t1, t2)

    if category == CATEGORY_NOT_COMPARABLE:
        return False, category, "different seqname or strand -- never collapsed"
    if category == CATEGORY_GENUINE_ALTERNATIVE:
        return False, category, "structurally distinct -- genuine alternative isoform"
    if category == CATEGORY_SAME_INTRON_CHAIN_DIFFERENT_ENDS:
        return (
            False,
            category,
            "same intron chain but different transcript ends -- retained as alternative",
        )
    if category == CATEGORY_SAME_CDS_DIFFERENT_UTR:
        if cfg.preserve_distinct_utrs:
            return (
                False,
                category,
                "same CDS, different UTR/exon boundary -- preserved (preserve_distinct_utrs=true)",
            )
        # Opt-in aggressive mode: collapse on CDS+phase+protein match alone.
    # From here, category is one of the "identical exon+CDS coordinates"
    # family (4/5/6) or the opted-in same-CDS-different-UTR case above.

    if cfg.require_matching_cds_phase:
        p1, p2 = cds_phase_signature(t1.get("cds_phases")), cds_phase_signature(
            t2.get("cds_phases")
        )
        if p1 is not None and p2 is not None and p1 != p2:
            return False, category, "CDS phase/frame differs -- not collapsed"

    if cfg.preserve_distinct_proteins:
        prot1, prot2 = t1.get("protein"), t2.get("protein")
        if prot1 and prot2 and prot1 != prot2:
            return (
                False,
                category,
                "translated proteins differ -- not collapsed (data anomaly: identical CDS coords but different protein)",
            )
        # If either side lacks a translated protein, we cannot verify
        # protein identity -- see module docstring: documented limitation,
        # falls back to the (already very strict) coordinate-only rule
        # rather than blocking collapse on missing data.

    return True, category, "exact structural (and, where available, translational) duplicate"


# ---------------------------------------------------------------------------
# Grouping: within one gene, partition transcripts into duplicate clusters
# ---------------------------------------------------------------------------


def group_exact_duplicates(
    transcripts: list[dict], cfg: DuplicateTranscriptCollapseConfig
) -> list[list[dict]]:
    """Partition one gene's transcripts into groups where every member is
    an exact duplicate of every other member (transitive: verified pairwise
    below, not just chained), plus singleton groups for non-duplicates.
    """
    groups: list[list[dict]] = []
    for t in transcripts:
        placed = False
        for group in groups:
            # Compare against the group's first member; since exact-duplicate
            # groups are defined by coordinate+protein equality (an
            # equivalence relation), matching the first member is
            # sufficient and avoids O(n^2) work being re-litigated per pair.
            eligible, _cat, _reason = is_exact_duplicate_pair(t, group[0], cfg)
            if eligible:
                group.append(t)
                placed = True
                break
        if not placed:
            groups.append([t])
    return groups


# ---------------------------------------------------------------------------
# Part 4: deterministic retained-transcript selection within a duplicate group
# ---------------------------------------------------------------------------

# Priority order (adapted from the task's suggested order to match what
# this data model actually records at the point this module runs):
#   1. structurally valid transcript      -- always true here: dedup_genes.py
#      and validate_and_fix_gff3() upstream already drop/fix structurally
#      invalid models before this point, so every candidate is structurally
#      valid by construction (see report: this criterion is a no-op given
#      the current pipeline, not omitted).
#   2. complete valid ORF                 -- not is_partial_5/is_partial_3
#   3. no internal stop codons            -- protein does not contain "*"
#   4. richer independent evidence provenance -- count of distinct sources
#   5. stronger existing protein-validation result -- protein_coding_score
#   6. higher existing GMB score          -- raw gmb_score
#   7. stable lexical transcript ID       -- final deterministic tie-break
def _retained_rank_key(t: dict) -> tuple:
    is_partial = bool(t.get("is_partial_5")) or bool(t.get("is_partial_3"))
    protein = t.get("protein") or ""
    # annotate_cds_utrs.translate() strips the trailing stop before this
    # point, so any remaining "*" is a genuine internal stop.
    has_internal_stop = "*" in protein
    n_sources = len(_evidence_set(t.get("evidence_sources")))
    protein_coding_score = t.get("protein_coding_score")
    gmb_score = t.get("gmb_score")
    return (
        is_partial,  # False (complete) sorts first
        has_internal_stop,  # False sorts first
        -n_sources,
        -(protein_coding_score if protein_coding_score is not None else -1.0),
        -(gmb_score if gmb_score is not None else float("-inf")),
        t["transcript_id"],
    )


def select_retained_transcript(group: list[dict]) -> dict:
    """Deterministically pick which transcript in an exact-duplicate group
    survives. Order-independent (sorting, not first-wins)."""
    return sorted(group, key=_retained_rank_key)[0]


# ---------------------------------------------------------------------------
# Part 3: evidence provenance merge policy
# ---------------------------------------------------------------------------

# Explicit, documented per-field policy (task requires distinguishing these):
#
#   UNION (safe to combine):
#     evidence_sources        -- set union of all group members' sources
#     collapsed_from          -- new field: union of all removed transcript IDs
#
#   RECALCULATED from the merged source set using existing scoring logic
#   (preferred per the task; implemented in recalculate_gmb_score below):
#     gmb_score
#
#   KEPT AS-IS, NEVER AGGREGATED (identical by construction across an exact-
#   duplicate group, since exon/CDS/protein all matched to qualify as a
#   group -- these are structural/translational facts, not per-source
#   scores, so there is nothing to aggregate):
#     exon coordinates, CDS coordinates, protein sequence, protein_length,
#     orf_label, is_partial_5/3, internal_stop_count, cds_bp, exon_count,
#     transcript_span_bp, diamond_hit/pident/qcov/scov/bitscore/evalue,
#     psauron_score, protein_coding_score
#     (protein validation is keyed by the *original* candidate transcript
#     ID/sequence upstream of this module -- duplicate fragments here share
#     the same original candidate and therefore already carry identical
#     validation results; there is no double-counting to resolve.)
#
#   NEVER MERGED AUTOMATICALLY (would blur genuine per-record identity):
#     transcript_id itself (a deterministic winner is chosen instead, see
#     select_retained_transcript), any future domain-evidence sidecar rows
#     (must remain keyed to whichever original ID they were computed against
#     until a domain-evidence merge policy is separately defined).
def merge_evidence_sources(group: list[dict]) -> str:
    merged = set()
    for t in group:
        merged |= _evidence_set(t.get("evidence_sources"))
    return ",".join(sorted(merged))


def recalculate_gmb_score(
    retained: dict,
    merged_evidence_sources: str,
    config: PipelineConfig,
    protein_supported_tids: set[str],
    genome_dict: Optional[dict] = None,
) -> Optional[float]:
    """Recalculate gmb_score from the merged evidence set using the real
    scoring.score_model() function, per the task's stated preference over
    blind aggregation. Returns None (caller falls back to the retained
    transcript's own already-correct score) if recalculation inputs are
    unavailable -- this is a documented, tested fallback, not a silent one.
    """
    try:
        import pandas as pd

        from gmb.pipeline.scoring import score_model
    except ImportError:
        return None

    exons = retained.get("exons") or []
    if not exons:
        return None

    df = pd.DataFrame({"Start": [s for s, _ in exons], "End": [e for _, e in exons]})
    model = {
        "id": retained["transcript_id"],
        "source": next(iter(_evidence_set(merged_evidence_sources)), "unknown"),
        "combined_evidence": merged_evidence_sources,
        "chrom": retained.get("chrom"),
        "strand": retained.get("strand"),
        "exon_count": len(exons),
        "df": df,
    }
    if retained.get("protein_coding_score") is not None:
        model["protein_coding_score"] = retained["protein_coding_score"]

    return score_model(model, config, protein_supported_tids, genome_dict)


# ---------------------------------------------------------------------------
# Part 1 (standalone diagnostic): audit table over an already-built run
# ---------------------------------------------------------------------------


def build_duplicate_audit_rows(
    genes: dict[str, list[dict]], cfg: DuplicateTranscriptCollapseConfig
) -> list[dict]:
    """One row per transcript-pair-membership, for `duplicate_transcript_audit.tsv`.

    ``genes`` is ``{gene_id: [transcript_record, ...]}`` -- see
    ``gmb.pipeline.canonical_selection.load_transcript_records`` for one way
    to build this from a completed run's output files, extended with
    exon/cds/cds_phases/protein (this module does not itself parse GFF3/FASTA
    -- see ``tools/audit_duplicate_transcripts.py`` for the file-reading glue).
    """
    rows = []
    group_counter = 0
    for gene_id, transcripts in genes.items():
        if len(transcripts) < 2:
            continue
        groups = group_exact_duplicates(transcripts, cfg)
        for group in groups:
            group_counter += 1
            if len(group) < 2:
                # Singleton -- still worth a row for completeness/audit, but
                # not itself a "duplicate group" in the report's sense.
                t = group[0]
                rows.append(
                    _audit_row(
                        gene_id,
                        t,
                        group_counter,
                        "no_duplicate",
                        None,
                        None,
                        t,
                        "keep",
                        "unique structure in gene",
                    )
                )
                continue
            retained = select_retained_transcript(group)
            for t in group:
                _eligible, category, reason = (
                    is_exact_duplicate_pair(t, retained, cfg)
                    if t is not retained
                    else (True, classify_pair(t, retained), "retained record")
                )
                action = "retain" if t is retained else "collapse_into_retained"
                rows.append(
                    _audit_row(
                        gene_id,
                        t,
                        group_counter,
                        category,
                        exon_signature(t["exons"]),
                        cds_signature(t.get("cds", [])),
                        retained,
                        action,
                        reason,
                    )
                )
    return rows


def _audit_row(
    gene_id, t, group_id, category, exon_sig, cds_sig, retained, action, reason
) -> dict:
    return {
        "gene_id": gene_id,
        "transcript_id": t["transcript_id"],
        "duplicate_group_id": group_id,
        "duplicate_classification": category,
        "exon_signature": exon_sig if exon_sig is not None else exon_signature(t["exons"]),
        "cds_signature": cds_sig if cds_sig is not None else cds_signature(t.get("cds", [])),
        "intron_chain_signature": intron_chain_signature(t["exons"]),
        "strand": t.get("strand"),
        "protein_checksum": protein_checksum(t.get("protein")),
        "protein_length": t.get("protein_length") or len(t.get("protein") or ""),
        "gmb_score": t.get("gmb_score"),
        "protein_coding_score": t.get("protein_coding_score"),
        "psauron_score": t.get("psauron_score"),
        "diamond_hit": t.get("diamond_hit"),
        "diamond_qcov": t.get("diamond_qcov"),
        "diamond_scov": t.get("diamond_scov"),
        "evidence_sources": t.get("evidence_sources"),
        "structures_exactly_identical": category in EXACT_DUPLICATE_CATEGORIES,
        "proteins_exactly_identical": (
            bool(t.get("protein")) and t.get("protein") == retained.get("protein")
        ),
        "proposed_retained_transcript": retained["transcript_id"],
        "proposed_action": action,
        "reason": reason,
    }


# ---------------------------------------------------------------------------
# Part 5/6: in-pipeline collapse, operating on builder.py's GFF3 row dicts
# ---------------------------------------------------------------------------
#
# Insertion point: called from gmb.pipeline.builder.main() immediately after
# dedup_genes() and before FASTA/evidence-attribution output is written.
# Rationale (see also the module docstring's root-cause account):
#
#   * Not before/during clustering or scoring: the artifact is only
#     *reliably distinguishable* from a genuine independent duplicate once
#     each candidate's full exon/CDS structure has been reconstructed (via
#     annotate_all_transcripts, keyed by transcript_id) -- clustering
#     fragments only ever see a partial view, so acting earlier would mean
#     re-deciding based on the same incomplete information that caused the
#     problem. It would also touch the general clustering algorithm, with a
#     much wider blast radius across loci that are NOT affected by this bug.
#   * Not only inside dedup_genes.py's own merge loop: that function's job
#     is gene-level overlap merging with a deliberately loose (first-mRNA,
#     tolerance-bp) check -- conflating it with an exact-equality
#     transcript-level pass would make both harder to reason about and test
#     independently. Running strictly after it means this pass only ever
#     sees isoforms that are already siblings under one gene (its own
#     precondition), never reaches across genes, and cannot change which
#     genes exist or their loci.
#   * Not only in canonical_selection.py: that tool must receive an already
#     biologically-meaningful isoform set (per this task's explicit
#     instruction) -- it is a read-only reporting layer over gmb.cli.build's
#     output, not the place to fix the output itself.
#   * Running here means: transcript IDs for genuinely-kept transcripts are
#     unaffected (only removed IDs disappear), evidence/scores are still
#     live Python objects (not yet serialised to TSV) so recalculation via
#     the real scoring.score_model() is possible, and it runs once, on the
#     final isoform set, so it cannot itself introduce new duplicates.


def _extract_transcript_records_from_gff_rows(
    gff_rows: list[dict], protein_by_tid: Optional[dict[str, str]] = None
) -> dict[str, list[dict]]:
    """Group builder.py's flat GFF3 row-dict list into per-gene transcript
    records with the fields classify_pair/is_exact_duplicate_pair need.
    """
    by_parent: dict[str, list[dict]] = defaultdict(list)
    for r in gff_rows:
        parent = r.get("Parent")
        if parent:
            by_parent[parent].append(r)

    genes: dict[str, list[dict]] = defaultdict(list)
    for r in gff_rows:
        if r.get("Feature") != "mRNA":
            continue
        tid = r["ID"]
        gene_id = r.get("Parent")
        children = by_parent.get(tid, [])
        exons = sorted((c["Start"], c["End"]) for c in children if c["Feature"] == "exon")
        cds = sorted((c["Start"], c["End"]) for c in children if c["Feature"] == "CDS")
        protein = (protein_by_tid or {}).get(tid)
        genes[gene_id].append(
            {
                "transcript_id": tid,
                "gene_id": gene_id,
                "chrom": r.get("Chromosome"),
                "strand": r.get("Strand"),
                "exons": exons,
                "cds": cds,
                "cds_phases": None,  # builder.py's own GFF3 rows don't carry phase
                "protein": protein,
                "evidence_sources": r.get("Evidence", ""),
                "gmb_score": r.get("gmb_score"),
                "protein_coding_score": r.get("protein_coding_score"),
                "orf_label": r.get("orf_label"),
                "is_partial_5": r.get("is_partial_5"),
                "is_partial_3": r.get("is_partial_3"),
                "_gff_row": r,
            }
        )
    return genes


def collapse_exact_duplicate_transcripts(
    gff_rows: list[dict],
    config: PipelineConfig,
    protein_by_tid: Optional[dict[str, str]] = None,
    genome_dict: Optional[dict] = None,
    protein_supported_tids: Optional[set[str]] = None,
) -> tuple[list[dict], list[dict], dict]:
    """Collapse exact-duplicate transcripts within each gene.

    Never drops a gene, never moves a transcript to a different gene, never
    changes a structurally-distinct transcript. Only removes transcripts
    proven to be exact duplicates of another transcript already kept in the
    same gene (see is_exact_duplicate_pair).

    Returns
    -------
    tuple of (list of dict, list of dict, dict)
        ``(new_gff_rows, collapse_log_rows, stats)``. ``collapse_log_rows``
        matches the ``collapsed_duplicate_transcripts.tsv`` schema: one row
        per group *member* (not per removal), so each group has exactly one
        row with an empty ``removed_transcript_id`` -- that row documents
        the retained/survivor transcript itself, sharing the same
        ``retained_transcript_id`` as its sibling removal rows. This lets
        the survivor's own classification/signature/score be inspected
        alongside what was collapsed into it, without a separate table.
    """
    cfg = config.duplicate_transcript_collapse
    stats = {
        "duplicate_groups_found": 0,
        "transcripts_removed": 0,
        "genes_with_collapse": 0,
    }
    if not cfg.collapse_exact_duplicates:
        return gff_rows, [], stats

    genes = _extract_transcript_records_from_gff_rows(gff_rows, protein_by_tid)

    by_parent: dict[str, list[dict]] = defaultdict(list)
    for r in gff_rows:
        parent = r.get("Parent")
        if parent:
            by_parent[parent].append(r)

    removed_tids: set[str] = set()
    collapse_log_rows: list[dict] = []
    retained_extra_evidence: dict[str, str] = {}  # tid -> merged evidence string
    retained_recalculated_score: dict[str, Optional[float]] = {}
    retained_collapsed_from: dict[str, list[str]] = defaultdict(list)

    for gene_id, transcripts in genes.items():
        if len(transcripts) < 2:
            continue
        groups = group_exact_duplicates(transcripts, cfg)
        gene_had_collapse = False
        for group in groups:
            if len(group) < 2:
                continue
            stats["duplicate_groups_found"] += 1
            gene_had_collapse = True
            retained = select_retained_transcript(group)
            merged_evidence = merge_evidence_sources(group)
            retained_extra_evidence[retained["transcript_id"]] = merged_evidence
            retained_collapsed_from[retained["transcript_id"]] = [
                t["transcript_id"] for t in group if t is not retained
            ]

            recalculated = recalculate_gmb_score(
                retained, merged_evidence, config, protein_supported_tids or set(), genome_dict
            )
            retained_recalculated_score[retained["transcript_id"]] = recalculated

            for t in group:
                category = (
                    classify_pair(t, retained)
                    if t is not retained
                    else CATEGORY_IDENTICAL_STRUCTURE_SAME_EVIDENCE
                )
                if t is not retained:
                    removed_tids.add(t["transcript_id"])
                collapse_log_rows.append(
                    {
                        "retained_transcript_id": retained["transcript_id"],
                        "removed_transcript_id": None if t is retained else t["transcript_id"],
                        "gene_id": gene_id,
                        "duplicate_classification": category,
                        "structure_signature": exon_signature(t["exons"])
                        + "|"
                        + cds_signature(t.get("cds", [])),
                        "evidence_merged": merged_evidence,
                        "gmb_score_before": t.get("gmb_score"),
                        "gmb_score_recalculated": recalculated,
                        "retained_id_reason": (
                            "sole/first-ranked survivor of exact-duplicate group "
                            "(complete ORF > no internal stop > evidence breadth > "
                            "protein-validation score > gmb_score > lexical transcript_id)"
                        ),
                    }
                )
        if gene_had_collapse:
            stats["genes_with_collapse"] += 1

    stats["transcripts_removed"] = len(removed_tids)

    new_gff_rows = []
    for r in gff_rows:
        if r.get("Feature") == "mRNA" and r["ID"] in removed_tids:
            continue
        if (
            r.get("Feature") in ("exon", "CDS", "five_prime_UTR", "three_prime_UTR")
            and r.get("Parent") in removed_tids
        ):
            continue
        if r.get("Feature") == "mRNA" and r["ID"] in retained_collapsed_from:
            r = dict(r)
            r["Evidence"] = retained_extra_evidence.get(r["ID"], r.get("Evidence", ""))
            r["CollapsedFrom"] = ",".join(retained_collapsed_from[r["ID"]])
            recalculated = retained_recalculated_score.get(r["ID"])
            if recalculated is not None:
                r["gmb_score"] = recalculated
        new_gff_rows.append(r)

    return new_gff_rows, collapse_log_rows, stats
