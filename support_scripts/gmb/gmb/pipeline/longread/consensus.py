"""Core consensus-building logic: clustering, grouping, suppression, gene IDs.

The two-stage design is justified by empirical observations on the P. falciparum
raw Minimap2 GTF (documented in the module docstring of the pipeline module):

* The raw file is NOT seqname-sorted — splitting first (Stage 1) is required.
* At typical long-read coverage depth, locus clustering collapses to ~1 cluster
  per strand for an entire chromosome, making the raw cluster index useless as
  a gene identity proxy.  Gene IDs are therefore re-derived post-hoc from the
  final (small) set of consensus models.

``_intron_chain_str`` uses coordinate snapping for READ GROUPING ONLY — it
groups reads whose splice sites agree within the configured tolerance into the
same intron-chain group so that minor alignment jitter does not fragment genuine
same-transcript reads.  This is a deliberate approximation for grouping, not a
structural matching claim; rescue matching in ``rescue.py`` uses explicit
coordinate comparison instead.
"""

from __future__ import annotations

from collections import defaultdict
from typing import Optional

from gmb.utils.intervals import reciprocal_overlap

from .config import LongreadConsensusConfig
from .io import load_shortread_transcripts, load_transcripts_from_split_file
from .rescue import decide_keep_or_rescue

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _median(values: list[int]) -> int:
    values = sorted(values)
    return values[len(values) // 2]


def _intron_chain_str(exons: list[tuple[int, int]], snap_bp: int = 0) -> str:
    """Intron-chain signature for grouping reads (same convention as scoring module).

    When ``snap_bp`` > 0 each boundary is rounded to the nearest multiple so
    that reads with minor splice-site jitter land in the same group.  This is
    used for GROUPING only — rescue matching uses explicit comparison instead.
    """
    if len(exons) < 2:
        return "single-exon"
    if snap_bp > 0:

        def snap(x: int) -> int:
            return round(x / snap_bp) * snap_bp

        return ",".join(
            f"{snap(exons[i][1])}-{snap(exons[i + 1][0])}" for i in range(len(exons) - 1)
        )
    return ",".join(f"{exons[i][1]}-{exons[i + 1][0]}" for i in range(len(exons) - 1))


def _consensus_exons_multi(members: list[dict]) -> list[tuple[int, int]]:
    """Median consensus exon coordinates for a group sharing one intron chain.

    Every position — including internal junctions, not just outer boundaries —
    is collapsed to its per-position median across the group.
    """
    n = len(sorted(members[0]["exons"]))
    exon_lists = [sorted(m["exons"]) for m in members]
    consensus = []
    for i in range(n):
        med_start = _median([ex[i][0] for ex in exon_lists])
        med_end = _median([ex[i][1] for ex in exon_lists])
        consensus.append((med_start, med_end))

    # Safety net: fall back to first member if per-position medians are non-monotonic
    for i in range(n - 1):
        if consensus[i][1] >= consensus[i + 1][0]:
            return exon_lists[0]
    return consensus


def _cluster_single_exon(reads: list[dict], tol_bp: int) -> list[list[dict]]:
    """Single-linkage clustering of single-exon reads by proximity.

    Uses a start-sorted sweep with retired-group pruning to bound the active
    set — an earlier all-pairs version was O(n²) and slow on real data.
    """
    if not reads:
        return []
    reads = sorted(reads, key=lambda r: r["start"])

    finished: list[list[dict]] = []
    active: list[list[dict]] = []
    active_bounds: list[list[int]] = []

    for r in reads:
        cutoff = r["start"] - tol_bp
        still_active = []
        still_bounds = []
        for grp, bounds in zip(active, active_bounds):
            if bounds[0] < cutoff:
                finished.append(grp)
            else:
                still_active.append(grp)
                still_bounds.append(bounds)
        active, active_bounds = still_active, still_bounds

        placed = False
        for grp, bounds in zip(active, active_bounds):
            if abs(r["start"] - bounds[0]) <= tol_bp and abs(r["end"] - bounds[1]) <= tol_bp:
                grp.append(r)
                bounds[0] = min(bounds[0], r["start"])
                bounds[1] = max(bounds[1], r["end"])
                placed = True
                break
        if not placed:
            active.append([r])
            active_bounds.append([r["start"], r["end"]])

    finished.extend(active)
    return finished


# ---------------------------------------------------------------------------
# Post-processing
# ---------------------------------------------------------------------------


def suppress_contained_models(records: list[dict]) -> list[dict]:
    """Drop any model whose span is fully contained within a higher-support model's span.

    Applied genome-wide (not per-cluster) because locus clustering at typical
    long-read coverage collapses to ~1 cluster per strand per chromosome.
    """
    kept: list[dict] = []
    for strand in ("+", "-"):
        strand_recs = [r for r in records if r["strand"] == strand]
        strand_recs.sort(key=lambda r: r["exons"][0][0])

        active: list[dict] = []
        for r in strand_recs:
            s, e = r["exons"][0][0], r["exons"][-1][1]
            active = [a for a in active if a["exons"][-1][1] >= s]
            contained = any(
                a["exons"][0][0] <= s
                and e <= a["exons"][-1][1]
                and a["read_support"] >= r["read_support"]
                for a in active
                if a is not r
            )
            if not contained:
                kept.append(r)
            active.append(r)
            active.sort(key=lambda a: -a["read_support"])
    return kept


def assign_gene_ids(records: list[dict], seqname: str, min_intergenic_gap: int) -> int:
    """Reassign gene_id by genomic overlap across the final consensus records.

    The raw locus-cluster id used during grouping is not a usable gene identity
    proxy at real long-read coverage depth.  This re-derives gene groups from
    the (much smaller) surviving models.

    Mutates records in place.  Returns the number of distinct genes assigned.
    """
    by_strand: dict[str, list[dict]] = defaultdict(list)
    for rec in records:
        by_strand[rec["strand"]].append(rec)

    n_genes = 0
    for strand, recs in by_strand.items():
        recs.sort(key=lambda r: r["exons"][0][0])
        cur_end = None
        group_idx = -1
        for rec in recs:
            start = rec["exons"][0][0]
            end = rec["exons"][-1][1]
            if cur_end is None or start - cur_end - 1 > min_intergenic_gap:
                group_idx += 1
                n_genes += 1
                cur_end = end
            else:
                cur_end = max(cur_end, end)
            rec["gene_id"] = f"MMConsensus_{seqname}_{strand}_{group_idx}"

    per_gene_counter: dict[str, int] = defaultdict(int)
    for rec in records:
        per_gene_counter[rec["gene_id"]] += 1
        rec["transcript_id"] = f"{rec['gene_id']}.t{per_gene_counter[rec['gene_id']]}"

    return n_genes


# ---------------------------------------------------------------------------
# Stage 2: per-seqname consensus building
# ---------------------------------------------------------------------------


def build_consensus_for_seqname(
    seqname: str,
    split_path: str,
    cfg: LongreadConsensusConfig,
    scallop_paths: Optional[list[str]] = None,
    stringtie_paths: Optional[list[str]] = None,
) -> tuple[list[dict], dict, list[dict]]:
    """Build consensus transcript models for one seqname's split file.

    Parameters
    ----------
    seqname : str
    split_path : str
        Path to the per-seqname GTF produced by Stage 1.
    cfg : LongreadConsensusConfig
    scallop_paths : list[str] or None
        One or more Scallop GTF paths.  Required only when rescue is enabled.
    stringtie_paths : list[str] or None
        One or more StringTie GTF paths.  Required only when rescue is enabled.

    Returns
    -------
    records : list[dict]
        Consensus transcript records.
    stats : dict
        Per-seqname summary statistics.
    dropped : list[dict]
        Per-read/group records for dropped-input attribution.
    """
    stats: dict = {"seqname": seqname}
    dropped: list[dict] = []

    transcripts, malformed = load_transcripts_from_split_file(split_path)
    stats["input_reads"] = len(transcripts)
    stats["malformed_input_rows"] = len(malformed)

    for m in malformed:
        dropped.append(
            {
                "source_tid": "UNKNOWN",
                "seqname": seqname,
                "start": None,
                "end": None,
                "strand": None,
                "n_exons": None,
                "rejection_stage": "INPUT_PARSE",
                "reason_code": m["reason"],
                "metric": None,
                "rescue_status": "N/A",
            }
        )

    # --- Per-read validation ---
    candidates = []
    for tid, rec in transcripts.items():
        if rec.get("mixed_strand"):
            dropped.append(
                {
                    "source_tid": tid,
                    "seqname": seqname,
                    "start": None,
                    "end": None,
                    "strand": rec["strand"],
                    "n_exons": len(rec["exons"]),
                    "rejection_stage": "FILTER_STRAND",
                    "reason_code": "MIXED_STRAND",
                    "metric": None,
                    "rescue_status": "N/A",
                }
            )
            continue

        if rec["strand"] not in ("+", "-"):
            dropped.append(
                {
                    "source_tid": tid,
                    "seqname": seqname,
                    "start": None,
                    "end": None,
                    "strand": rec["strand"],
                    "n_exons": len(rec["exons"]),
                    "rejection_stage": "FILTER_STRAND",
                    "reason_code": "UNRESOLVED_STRAND",
                    "metric": rec["strand"],
                    "rescue_status": "N/A",
                }
            )
            continue

        exons = sorted(rec["exons"])
        if not exons:
            dropped.append(
                {
                    "source_tid": tid,
                    "seqname": seqname,
                    "start": None,
                    "end": None,
                    "strand": rec["strand"],
                    "n_exons": 0,
                    "rejection_stage": "FILTER_INPUT",
                    "reason_code": "NO_EXON_ROWS",
                    "metric": None,
                    "rescue_status": "N/A",
                }
            )
            continue

        span = exons[-1][1] - exons[0][0] + 1
        if span > cfg.max_span_length:
            dropped.append(
                {
                    "source_tid": tid,
                    "seqname": seqname,
                    "start": exons[0][0],
                    "end": exons[-1][1],
                    "strand": rec["strand"],
                    "n_exons": len(exons),
                    "rejection_stage": "FILTER_SPAN",
                    "reason_code": "SPAN_EXCEEDS_MAX",
                    "metric": span,
                    "rescue_status": "N/A",
                }
            )
            continue

        bad_intron = False
        max_intron = 0
        if len(exons) > 1:
            for i in range(len(exons) - 1):
                gap = exons[i + 1][0] - exons[i][1] - 1
                if gap > max_intron:
                    max_intron = gap
                if gap < 0 or gap > cfg.max_intron_length:
                    bad_intron = True
                    break
        if bad_intron:
            dropped.append(
                {
                    "source_tid": tid,
                    "seqname": seqname,
                    "start": exons[0][0],
                    "end": exons[-1][1],
                    "strand": rec["strand"],
                    "n_exons": len(exons),
                    "rejection_stage": "FILTER_INTRON",
                    "reason_code": "INTRON_EXCEEDS_MAX",
                    "metric": max_intron,
                    "rescue_status": "N/A",
                }
            )
            continue

        candidates.append(
            {
                "tid": tid,
                "gene_id": rec["gene_id"],
                "strand": rec["strand"],
                "exons": exons,
                "start": exons[0][0],
                "end": exons[-1][1],
            }
        )

    stats["dropped_no_exons"] = sum(1 for d in dropped if d["reason_code"] == "NO_EXON_ROWS")
    stats["dropped_unresolved_strand"] = sum(
        1 for d in dropped if d["reason_code"] in ("UNRESOLVED_STRAND", "MIXED_STRAND")
    )
    stats["dropped_bad_intron_or_overlap"] = sum(
        1 for d in dropped if d["reason_code"] == "INTRON_EXCEEDS_MAX"
    )
    stats["dropped_oversized_span"] = sum(
        1 for d in dropped if d["reason_code"] == "SPAN_EXCEEDS_MAX"
    )
    stats["surviving_after_basic_filters"] = len(candidates)

    # --- Load short-read transcripts if rescue is enabled ---
    rc = cfg.shortread_rescue
    sr_data: Optional[dict] = None
    if rc.enabled and (rc.multi_exon.enabled or rc.single_exon.enabled):
        sr_data = load_shortread_transcripts(scallop_paths or [], stringtie_paths or [], seqname)
    stats["sr_transcripts_loaded"] = len(sr_data["all"]) if sr_data is not None else 0

    # --- Cluster by strand + genomic overlap ---
    consensus_records: list[dict] = []
    n_clusters = 0
    n_chain_groups_total = 0
    n_dropped_low_support = 0
    n_dropped_no_rescue = 0
    n_suppressed_fragment_of_multiexon = 0
    n_rescued = 0

    for strand in ("+", "-"):
        strand_candidates = sorted(
            (c for c in candidates if c["strand"] == strand), key=lambda c: c["start"]
        )
        clusters: list[list[dict]] = []
        cur_cluster: list[dict] = []
        cur_end = None
        for c in strand_candidates:
            if cur_cluster and c["start"] - cur_end - 1 > cfg.min_intergenic_gap:
                clusters.append(cur_cluster)
                cur_cluster = []
                cur_end = None
            cur_cluster.append(c)
            cur_end = c["end"] if cur_end is None else max(cur_end, c["end"])
        if cur_cluster:
            clusters.append(cur_cluster)

        n_clusters += len(clusters)

        for cluster_idx, cluster in enumerate(clusters):
            cluster_size = len(cluster)
            multi = [c for c in cluster if len(c["exons"]) > 1]
            single = [c for c in cluster if len(c["exons"]) == 1]

            chain_groups: dict[str, list[dict]] = defaultdict(list)
            for c in multi:
                chain_groups[
                    _intron_chain_str(c["exons"], snap_bp=cfg.splice_site_tolerance_bp)
                ].append(c)

            single_groups = _cluster_single_exon(single, cfg.collapse_terminal_variation_bp)

            n_chain_groups_total += len(chain_groups) + len(single_groups)

            group_id = 0
            accepted_multiexon_exons: list[tuple[int, int]] = []

            for _chain, members in chain_groups.items():
                group_id += 1
                support = len(members)
                min_req = cfg.min_read_support_multi_exon

                kept, rescued, reason, sr_tids, rej_reason = decide_keep_or_rescue(
                    support, min_req, members, strand, cfg, sr_data
                )

                if not kept:
                    rc_str = (
                        "DROPPED_LOW_SUPPORT" if support >= 2 else "DROPPED_SINGLE_READ_NO_RESCUE"
                    )
                    for m in members:
                        dropped.append(
                            {
                                "source_tid": m["tid"],
                                "seqname": seqname,
                                "start": m["start"],
                                "end": m["end"],
                                "strand": strand,
                                "n_exons": len(m["exons"]),
                                "rejection_stage": "SUPPORT",
                                "reason_code": rc_str,
                                "metric": support,
                                "rescue_status": rej_reason or "NO_RESCUE",
                            }
                        )
                    if support > 1 or support < min_req:
                        if support != 1:
                            n_dropped_low_support += 1
                        else:
                            n_dropped_no_rescue += 1
                    continue

                if rescued:
                    n_rescued += 1

                exons_consensus = _consensus_exons_multi(members)
                accepted_multiexon_exons.extend(exons_consensus)
                consensus_records.append(
                    {
                        "seqname": seqname,
                        "strand": strand,
                        "gene_id": f"MMConsensus_{seqname}_{strand}_{cluster_idx}",
                        "transcript_id": f"MMConsensus_{seqname}_{strand}_{cluster_idx}_g{group_id}",
                        "exons": exons_consensus,
                        "read_support": support,
                        "cluster_size": cluster_size,
                        "rescued_by_shortread": rescued,
                        "rescue_reason": reason,
                        "rescue_sr_tids": sr_tids,
                        "n_exons": len(exons_consensus),
                    }
                )

            for members in single_groups:
                group_id += 1

                # Suppress single-exon groups that overlap an accepted multi-exon
                # consensus at this locus — likely truncated/partial reads.
                med_start_probe = sorted(m["start"] for m in members)[len(members) // 2]
                med_end_probe = sorted(m["end"] for m in members)[len(members) // 2]
                if any(
                    reciprocal_overlap(med_start_probe, med_end_probe, es, ee) >= 0.5
                    or (es <= med_start_probe and med_end_probe <= ee)
                    for es, ee in accepted_multiexon_exons
                ):
                    n_suppressed_fragment_of_multiexon += 1
                    for m in members:
                        dropped.append(
                            {
                                "source_tid": m["tid"],
                                "seqname": seqname,
                                "start": m["start"],
                                "end": m["end"],
                                "strand": strand,
                                "n_exons": 1,
                                "rejection_stage": "SUPPRESSION",
                                "reason_code": "SUPPRESSED_FRAGMENT_OF_MULTIEXON",
                                "metric": None,
                                "rescue_status": "N/A",
                            }
                        )
                    continue

                support = len(members)
                min_req = cfg.min_read_support_single_exon

                kept, rescued, reason, sr_tids, rej_reason = decide_keep_or_rescue(
                    support, min_req, members, strand, cfg, sr_data
                )

                if not kept:
                    rc_str = (
                        "DROPPED_LOW_SUPPORT" if support >= 2 else "DROPPED_SINGLE_READ_NO_RESCUE"
                    )
                    for m in members:
                        dropped.append(
                            {
                                "source_tid": m["tid"],
                                "seqname": seqname,
                                "start": m["start"],
                                "end": m["end"],
                                "strand": strand,
                                "n_exons": 1,
                                "rejection_stage": "SUPPORT",
                                "reason_code": rc_str,
                                "metric": support,
                                "rescue_status": rej_reason or "NO_RESCUE",
                            }
                        )
                    if support != 1:
                        n_dropped_low_support += 1
                    else:
                        n_dropped_no_rescue += 1
                    continue

                if rescued:
                    n_rescued += 1

                starts = [m["start"] for m in members]
                ends = [m["end"] for m in members]
                starts.sort()
                ends.sort()
                med_start = starts[len(starts) // 2]
                med_end = ends[len(ends) // 2]
                consensus_records.append(
                    {
                        "seqname": seqname,
                        "strand": strand,
                        "gene_id": f"MMConsensus_{seqname}_{strand}_{cluster_idx}",
                        "transcript_id": f"MMConsensus_{seqname}_{strand}_{cluster_idx}_g{group_id}",
                        "exons": [(med_start, med_end)],
                        "read_support": support,
                        "cluster_size": cluster_size,
                        "rescued_by_shortread": rescued,
                        "rescue_reason": reason,
                        "rescue_sr_tids": sr_tids,
                        "n_exons": 1,
                    }
                )

    stats["clusters_found"] = n_clusters
    stats["chain_or_overlap_groups_found"] = n_chain_groups_total
    stats["dropped_low_support"] = n_dropped_low_support
    stats["dropped_single_read_no_rescue"] = n_dropped_no_rescue
    stats["suppressed_fragment_of_multiexon"] = n_suppressed_fragment_of_multiexon
    stats["rescued_by_shortread"] = n_rescued
    stats["consensus_models_before_containment_suppression"] = len(consensus_records)

    n_suppressed_contained = 0
    if cfg.suppress_contained_models:
        before = len(consensus_records)
        suppressed_set = set()
        surviving = suppress_contained_models(consensus_records)
        surviving_ids = {r["transcript_id"] for r in surviving}
        for rec in consensus_records:
            if rec["transcript_id"] not in surviving_ids:
                suppressed_set.add(rec["transcript_id"])
                dropped.append(
                    {
                        "source_tid": rec["transcript_id"],
                        "seqname": seqname,
                        "start": rec["exons"][0][0],
                        "end": rec["exons"][-1][1],
                        "strand": rec["strand"],
                        "n_exons": rec["n_exons"],
                        "rejection_stage": "SUPPRESSION",
                        "reason_code": "SUPPRESSED_CONTAINED",
                        "metric": rec["read_support"],
                        "rescue_status": "N/A",
                    }
                )
        consensus_records = surviving
        n_suppressed_contained = before - len(consensus_records)

    stats["suppressed_contained_in_higher_support_model"] = n_suppressed_contained
    stats["consensus_models_output"] = len(consensus_records)

    n_genes = assign_gene_ids(consensus_records, seqname, cfg.min_intergenic_gap)
    stats["genes_output"] = n_genes

    return consensus_records, stats, dropped
