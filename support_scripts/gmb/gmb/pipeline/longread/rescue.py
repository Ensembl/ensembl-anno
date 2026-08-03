"""Structure-aware short-read rescue for single-read consensus candidates.

Junction matching uses explicit coordinate comparison, not coordinate snapping
or frozensets:

    abs(candidate_donor  - sr_donor)    <= tol_bp   (donor agreement)
    abs(candidate_acceptor - sr_acceptor) <= tol_bp  (acceptor agreement)

This avoids bin-boundary failures: two junctions differing by exactly
``tol_bp`` would snap to different bins when ``tol_bp`` divides one
coordinate cleanly, but the explicit comparison correctly classifies them
as matching.

Donor/acceptor conventions (GTF coordinate space, 1-based inclusive):
    donor    = rightmost base of exon i   (exons[i][1])
    acceptor = leftmost  base of exon i+1 (exons[i+1][0])

Named short-read sources (Scallop, StringTie) are tracked independently
throughout so that policies requiring both assemblers (M3, S2) can enforce
independent confirmation.  They share one biological evidence class
``short_read_transcriptomic`` — two assemblers do not become two classes.
"""

from __future__ import annotations

from typing import Optional

from gmb.utils.intervals import reciprocal_overlap

from .config import LongreadConsensusConfig, ShortreadRescueConfig

# ---------------------------------------------------------------------------
# Junction extraction — ordered, no snapping
# ---------------------------------------------------------------------------


def extract_junctions(exons: list[tuple[int, int]]) -> list[tuple[int, int]]:
    """Return ordered (donor, acceptor) intron pairs from a sorted exon list.

    No coordinate snapping.  Caller must pass sorted exons.
    Returns an empty list for single-exon or empty inputs.
    """
    exons = sorted(exons)
    if len(exons) < 2:
        return []
    return [(exons[i][1], exons[i + 1][0]) for i in range(len(exons) - 1)]


def junction_matches(
    cand_d: int,
    cand_a: int,
    sr_d: int,
    sr_a: int,
    tol_bp: int,
) -> bool:
    """True if candidate and SR junctions agree within ``tol_bp`` on both ends."""
    return abs(cand_d - sr_d) <= tol_bp and abs(cand_a - sr_a) <= tol_bp


# ---------------------------------------------------------------------------
# Multi-exon rescue policies
# ---------------------------------------------------------------------------


def m1_junction_complete(
    cand_exons: list[tuple[int, int]],
    strand: str,
    sr_transcripts: list[dict],
    tol_bp: int,
) -> tuple[bool, list[str]]:
    """M1: every candidate intron junction is confirmed by at least one SR transcript.

    Combined support across multiple SR transcripts is allowed and expected —
    a short SR transcript covering just part of the intron chain is still valid
    evidence for the junctions it contains.

    Returns ``(kept, support_tids)`` where ``support_tids`` lists every SR
    transcript that contributed at least one matching junction.
    """
    cand_junctions = extract_junctions(cand_exons)
    if not cand_junctions:
        return False, []

    # Pool all SR junctions (from matching-strand transcripts)
    all_sr_junctions: list[tuple[int, int]] = []
    for tx in sr_transcripts:
        if tx["strand"] != strand:
            continue
        all_sr_junctions.extend(extract_junctions(tx["exons"]))

    # All candidate junctions must have at least one match in the SR pool
    for cand_d, cand_a in cand_junctions:
        if not any(
            junction_matches(cand_d, cand_a, sr_d, sr_a, tol_bp) for sr_d, sr_a in all_sr_junctions
        ):
            return False, []

    # Collect which transcripts contributed at least one matching junction
    support_tids: list[str] = []
    for tx in sr_transcripts:
        if tx["strand"] != strand:
            continue
        tx_junctions = extract_junctions(tx["exons"])
        if any(
            any(
                junction_matches(cand_d, cand_a, sr_d, sr_a, tol_bp) for sr_d, sr_a in tx_junctions
            )
            for cand_d, cand_a in cand_junctions
        ):
            support_tids.append(tx["tid"])

    return True, support_tids


def m2_full_chain(
    cand_exons: list[tuple[int, int]],
    strand: str,
    sr_transcripts: list[dict],
    tol_bp: int,
) -> tuple[bool, list[str]]:
    """M2: the candidate's full ordered intron chain matches at least one SR transcript.

    Requires the same intron count AND every corresponding (donor, acceptor)
    pair within ``tol_bp``.  Uses ordered list comparison — not frozenset —
    so insertion or deletion of an intron is never silently collapsed.
    """
    cand_junctions = extract_junctions(cand_exons)
    if not cand_junctions:
        return False, []

    support_tids: list[str] = []
    for tx in sr_transcripts:
        if tx["strand"] != strand:
            continue
        tx_junctions = extract_junctions(tx["exons"])
        if len(tx_junctions) != len(cand_junctions):
            continue
        if all(
            junction_matches(cj[0], cj[1], tj[0], tj[1], tol_bp)
            for cj, tj in zip(cand_junctions, tx_junctions)
        ):
            support_tids.append(tx["tid"])

    return bool(support_tids), support_tids


def m3_dual_assembler(
    cand_exons: list[tuple[int, int]],
    strand: str,
    scallop_txs: list[dict],
    stringtie_txs: list[dict],
    tol_bp: int,
) -> tuple[bool, list[str]]:
    """M3: junctions independently confirmed by both Scallop AND StringTie (M1-level each)."""
    ok_s, tids_s = m1_junction_complete(cand_exons, strand, scallop_txs, tol_bp)
    ok_st, tids_st = m1_junction_complete(cand_exons, strand, stringtie_txs, tol_bp)
    if ok_s and ok_st:
        return True, tids_s + tids_st
    return False, []


# ---------------------------------------------------------------------------
# Single-exon rescue policies
# ---------------------------------------------------------------------------


def s1_reciprocal_boundary(
    cand_exons: list[tuple[int, int]],
    strand: str,
    sr_transcripts: list[dict],
    min_overlap: float,
    boundary_tol_bp: int,
    require_single_exon_sr: bool,
) -> tuple[bool, list[str], Optional[str]]:
    """S1: a single-exon SR transcript reciprocally overlaps the candidate.

    When ``require_single_exon_sr`` is True (the default), only SR transcripts
    that are themselves single-exon are evaluated.  A single exon belonging to
    a multi-exon SR transcript does NOT rescue the candidate — instead, if
    that is the only available SR support, the rejection reason
    ``LIKELY_FRAGMENT_OF_MULTIEXON_TRANSCRIPT`` is recorded.

    Boundary agreement: both the 5′ and 3′ boundaries must be within
    ``boundary_tol_bp``.  Set ``boundary_tol_bp=0`` to disable this check
    (reciprocal overlap alone is then sufficient).

    Returns ``(kept, support_tids, rejection_reason)`` where
    ``rejection_reason`` is ``None`` when kept, or a reason code when not.
    """
    if len(cand_exons) != 1:
        return False, [], None
    c_start, c_end = cand_exons[0]

    support_tids: list[str] = []
    has_multiexon_sr_support = False

    for tx in sr_transcripts:
        if tx["strand"] != strand:
            continue
        is_single_exon_sr = len(tx["exons"]) == 1

        if require_single_exon_sr and not is_single_exon_sr:
            # Check if this multi-exon SR transcript could explain the LR
            # candidate as a fragment of one of its exons.
            for e_start, e_end in tx["exons"]:
                if reciprocal_overlap(c_start, c_end, e_start, e_end) >= min_overlap:
                    has_multiexon_sr_support = True
                    break
            continue

        # Evaluate single-exon SR transcript (or any SR if require_single_exon_sr=False)
        for e_start, e_end in tx["exons"]:
            ro = reciprocal_overlap(c_start, c_end, e_start, e_end)
            if ro < min_overlap:
                continue
            if boundary_tol_bp > 0:
                if (
                    abs(c_start - e_start) > boundary_tol_bp
                    or abs(c_end - e_end) > boundary_tol_bp
                ):
                    continue
            support_tids.append(tx["tid"])
            break  # one satisfying exon per transcript is enough

    if support_tids:
        return True, support_tids, None

    if require_single_exon_sr and has_multiexon_sr_support:
        return False, [], "LIKELY_FRAGMENT_OF_MULTIEXON_TRANSCRIPT"

    return False, [], None


def s2_dual_assembler_boundary(
    cand_exons: list[tuple[int, int]],
    strand: str,
    scallop_txs: list[dict],
    stringtie_txs: list[dict],
    min_overlap: float,
    boundary_tol_bp: int,
    require_single_exon_sr: bool,
) -> tuple[bool, list[str], Optional[str]]:
    """S2: S1-level support independently from both Scallop AND StringTie."""
    ok_s, tids_s, _ = s1_reciprocal_boundary(
        cand_exons, strand, scallop_txs, min_overlap, boundary_tol_bp, require_single_exon_sr
    )
    ok_st, tids_st, _ = s1_reciprocal_boundary(
        cand_exons, strand, stringtie_txs, min_overlap, boundary_tol_bp, require_single_exon_sr
    )
    if ok_s and ok_st:
        return True, tids_s + tids_st, None
    return False, [], None


def s3_independent_evidence(
    cand_exons: list[tuple[int, int]],
    strand: str,
    sr_transcripts: list[dict],
    min_overlap: float,
    boundary_tol_bp: int,
    require_single_exon_sr: bool,
    n_required: int = 2,
) -> tuple[bool, list[str], Optional[str]]:
    """S3: at least ``n_required`` independent SR transcripts each satisfying S1."""
    if len(cand_exons) != 1:
        return False, [], None
    c_start, c_end = cand_exons[0]

    has_multiexon_sr_support = False
    support_tids: list[str] = []

    for tx in sr_transcripts:
        if tx["strand"] != strand:
            continue
        is_single_exon_sr = len(tx["exons"]) == 1

        if require_single_exon_sr and not is_single_exon_sr:
            for e_start, e_end in tx["exons"]:
                if reciprocal_overlap(c_start, c_end, e_start, e_end) >= min_overlap:
                    has_multiexon_sr_support = True
                    break
            continue

        for e_start, e_end in tx["exons"]:
            ro = reciprocal_overlap(c_start, c_end, e_start, e_end)
            if ro < min_overlap:
                continue
            if boundary_tol_bp > 0:
                if (
                    abs(c_start - e_start) > boundary_tol_bp
                    or abs(c_end - e_end) > boundary_tol_bp
                ):
                    continue
            support_tids.append(tx["tid"])
            break

    if len(support_tids) >= n_required:
        return True, support_tids, None

    if require_single_exon_sr and has_multiexon_sr_support and not support_tids:
        return False, [], "LIKELY_FRAGMENT_OF_MULTIEXON_TRANSCRIPT"

    return False, [], None


# ---------------------------------------------------------------------------
# Dispatch functions
# ---------------------------------------------------------------------------


def try_rescue_multi_exon(
    cand_exons: list[tuple[int, int]],
    strand: str,
    sr_data: dict,
    cfg: LongreadConsensusConfig,
) -> tuple[bool, Optional[str], list[str]]:
    """Dispatch structure-aware rescue for a single-read multi-exon candidate.

    Returns ``(kept, reason_code, supporting_sr_tids)``.
    """
    rc = cfg.shortread_rescue
    if not rc.enabled or not rc.multi_exon.enabled:
        return False, None, []

    me = rc.multi_exon
    tol = me.junction_tolerance_bp

    if me.policy == "M1":
        kept, tids = m1_junction_complete(cand_exons, strand, sr_data["all"], tol)
        return kept, "M1_junction_complete" if kept else None, tids

    if me.policy == "M2":
        kept, tids = m2_full_chain(cand_exons, strand, sr_data["all"], tol)
        return kept, "M2_full_chain" if kept else None, tids

    if me.policy == "M3":
        kept, tids = m3_dual_assembler(
            cand_exons, strand, sr_data["scallop"], sr_data["stringtie"], tol
        )
        return kept, "M3_dual_assembler" if kept else None, tids

    return False, None, []


def try_rescue_single_exon(
    cand_exons: list[tuple[int, int]],
    strand: str,
    sr_data: dict,
    cfg: LongreadConsensusConfig,
) -> tuple[bool, Optional[str], list[str], Optional[str]]:
    """Dispatch structure-aware rescue for a single-read single-exon candidate.

    Returns ``(kept, reason_code, supporting_sr_tids, rejection_reason)``
    where ``rejection_reason`` is a machine-readable code (e.g.
    ``LIKELY_FRAGMENT_OF_MULTIEXON_TRANSCRIPT``) when rescue fails for a
    specific structural reason, or ``None`` when it simply lacks support.
    """
    rc = cfg.shortread_rescue
    if not rc.enabled or not rc.single_exon.enabled:
        return False, None, [], None

    se = rc.single_exon
    min_ro = se.reciprocal_overlap_min
    btol = se.boundary_tolerance_bp
    req_single = se.require_single_exon_sr_model

    if se.policy == "S1":
        kept, tids, rej = s1_reciprocal_boundary(
            cand_exons, strand, sr_data["all"], min_ro, btol, req_single
        )
        return kept, "S1_reciprocal_boundary" if kept else None, tids, rej

    if se.policy == "S2":
        kept, tids, rej = s2_dual_assembler_boundary(
            cand_exons, strand, sr_data["scallop"], sr_data["stringtie"], min_ro, btol, req_single
        )
        return kept, "S2_dual_assembler" if kept else None, tids, rej

    if se.policy == "S3":
        kept, tids, rej = s3_independent_evidence(
            cand_exons, strand, sr_data["all"], min_ro, btol, req_single
        )
        return kept, "S3_independent_evidence" if kept else None, tids, rej

    return False, None, [], None


def decide_keep_or_rescue(
    support: int,
    min_req: int,
    members: list[dict],
    strand: str,
    cfg: LongreadConsensusConfig,
    sr_data: Optional[dict],
) -> tuple[bool, bool, Optional[str], list[str], Optional[str]]:
    """Central keep/rescue decision for one consensus group.

    Returns ``(kept, rescued, reason_code, support_sr_tids, rejection_reason)``.

    - ``kept``: whether the group survives.
    - ``rescued``: True only when kept via short-read rescue (not own read support).
    - ``reason_code``: which rescue policy fired (None if kept by own support).
    - ``support_sr_tids``: SR transcript IDs that contributed rescue evidence.
    - ``rejection_reason``: structural rejection code when not kept (or None).
    """
    if support >= min_req:
        return True, False, None, [], None

    if support != 1:
        return False, False, None, [], None

    if sr_data is None:
        return False, False, None, [], None

    exons = members[0]["exons"]

    if len(exons) > 1:
        kept, reason, tids = try_rescue_multi_exon(exons, strand, sr_data, cfg)
        return kept, kept, reason, tids, None
    else:
        kept, reason, tids, rej_reason = try_rescue_single_exon(exons, strand, sr_data, cfg)
        return kept, kept, reason, tids, rej_reason
