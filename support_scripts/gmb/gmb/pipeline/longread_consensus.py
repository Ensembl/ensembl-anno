#!/usr/bin/env python3
"""Long-read (Minimap2) raw-alignment-to-consensus preprocessing tool.

Raw long-read GTFs (one record per aligned read, not per transcript) can be
too large to load with pandas/pyranges and too noisy to use as GMB evidence
directly -- see the module docstring sections below for the empirical
justification. This tool converts such a file into a small, provenance-
preserving consensus transcript track (``Source=Minimap2Consensus``) that
*is* safe to pass to ``gmb.cli.build --minimap2``.

Why this exists (empirical findings on GCA_000002765.3's raw Minimap2 GTF)
---------------------------------------------------------------------------
* 14.5 GB / 130.4M lines / 36.2M "transcript" records for a 23 Mb genome --
  unsafe to load wholesale with pandas/pyranges on a 16 GB machine.
* The file is NOT sorted by seqname: 6.59M chromosome-block transitions
  across 36.2M transcript records (a chromosome's records appear in
  hundreds of thousands of scattered blocks, not one contiguous run) --
  a naive "flush on seqname change" single pass is unsafe. This tool
  therefore always splits by seqname first (Stage 1) before any per-locus
  clustering (Stage 2).
* ~2.93% of transcript-level summary lines have swapped start/end
  coordinates (a bug in whatever produced this GTF, not something to
  "fix" -- this tool never trusts the transcript line's own coordinates,
  deriving every transcript's span purely from its own exon rows, which
  were not found to have this problem).
* Raw per-read span/intron lengths vastly exceed real P. falciparum
  biology (verified independently from the reference GFF: real max
  transcript span 32,487bp, real max intron 2,425bp) -- raw data:
  span p99=344,020bp, p99.9=857,452bp; intron p50=1,051bp (already above
  the real p90), p90=38,432bp. Most multi-exon raw records include at
  least one non-biological gap, consistent with chimeric/split-alignment
  artifacts rather than genuine spliced transcripts.
* 63.3% of raw records are single-exon, far above the true ~45.6%
  single-exon gene rate -- most single-exon reads are truncated fragments
  of real multi-exon transcripts, not independent single-exon evidence.
  Single-exon consensus models are held to a stricter support bar.

Two-stage design
-----------------
Stage 1 (``split_by_seqname``): one sequential pass over the raw file,
appending each line to one of N per-seqname temp files (O(1) memory beyond
small write buffers -- safe regardless of input size).

Stage 2 (``build_consensus_for_seqname``): per seqname (now a much smaller
file), reconstruct transcripts from exon rows, reject clearly-chimeric
records (oversized span/intron), cluster by strand + genomic overlap
(splitting on large gaps), group by exact intron chain (or overlap-collapse
for single-exon reads), and emit one consensus record per group meeting a
minimum-read-support bar -- with an optional short-read-corroboration
rescue path for otherwise-unsupported single-read models. Plain Python data
structures throughout (no pandas/pyranges) to keep per-seqname memory
predictable at this record scale.
"""

from __future__ import annotations

import argparse
import json
import os
import re
import sys
import time
from collections import defaultdict
from dataclasses import dataclass
from typing import Optional

import yaml

from gmb.utils.intervals import reciprocal_overlap
from gmb.utils.io import ensure_dir
from gmb.utils.logging import resolve_log_file, setup_logging

ATTR_RE = re.compile(r'(\w+) "([^"]*)"')


# ---------------------------------------------------------------------------
# Config
# ---------------------------------------------------------------------------


@dataclass
class LongreadConsensusConfig:
    """Thresholds for raw-read -> consensus collapsing.

    Defaults are grounded in stats computed directly from GCA_000002765.3's
    reference GFF (real max transcript span 32,487bp, real max intron
    2,425bp) and from the raw Minimap2 GTF itself (see module docstring).
    Not a general-purpose default for arbitrary genomes -- re-derive these
    from the target genome's own biology where possible.
    """

    # Per-read rejection (chimera / genomic-DNA / mis-alignment screening)
    max_span_length: int = 45_000
    max_intron_length: int = 3_000

    # Locus clustering
    min_intergenic_gap: int = 500

    # Consensus grouping. 100bp matches established long-read transcript
    # collapsing practice (e.g. cDNA_Cupcake's default 5'/3' tolerances) for
    # alignment-boundary noise -- an initial 20bp default was found (via a
    # real seqname-1 run) to fragment genuine same-transcript reads into
    # many spurious low-support groups.
    collapse_terminal_variation_bp: int = 100

    # Tolerance (bp) for snapping intron boundaries before grouping by
    # intron-chain signature. Real long-read splice-site calls can jitter by
    # a few bp in low-complexity/AT-rich regions without being a genuinely
    # different splice choice; exact-string chain matching without this
    # fragments true isoforms.
    splice_site_tolerance_bp: int = 10

    # Minimum independent-read support to keep a consensus model
    min_read_support_multi_exon: int = 2
    min_read_support_single_exon: int = 4

    # Rescue path: keep an otherwise-unsupported single-read model only if
    # it is corroborated by short-read (Scallop/StringTie) evidence at the
    # same locus/strand with compatible structure. Defaults to OFF: on a
    # real seqname-1 run this rescue path was found to be the dominant
    # source of remaining over-fragmentation (76% of output models were
    # single-read models rescued this way) -- "any overlapping short-read
    # exon" is not discriminative at typical Scallop/StringTie coverage
    # depth, since almost any genuine (non-chimeric) exonic fragment
    # overlaps *some* short-read exon. Left available for genomes/datasets
    # where short-read coverage is sparser and this signal is more
    # informative.
    require_shortread_support_for_single_read_models: bool = False
    shortread_overlap_fraction: float = 0.5

    # Genome-wide (not just within one raw-read locus cluster) containment
    # suppression: after all clusters are processed, drop any consensus
    # model whose exon span is fully contained within a higher-read-support
    # model's span on the same strand. On real data, raw-read locus
    # clustering collapses to ~1 cluster per strand for an entire
    # chromosome at typical long-read coverage depth (see module
    # docstring), so within-cluster suppression alone is insufficient --
    # this catches partial/truncated reads (fewer introns than the full
    # transcript, so they never share an intron-chain group with it) that
    # would otherwise survive as spurious extra "genes".
    suppress_contained_models: bool = True


def load_longread_consensus_config(path: Optional[str]) -> LongreadConsensusConfig:
    cfg = LongreadConsensusConfig()
    if path and os.path.exists(path):
        with open(path) as fh:
            data = yaml.safe_load(fh) or {}
        for key, val in data.items():
            if not hasattr(cfg, key):
                raise ValueError(f"Unknown longread_consensus config key: '{key}'")
            setattr(cfg, key, val)
    return cfg


# ---------------------------------------------------------------------------
# Stage 1: split the raw file by seqname
# ---------------------------------------------------------------------------


def split_by_seqname(
    input_path: str,
    split_dir: str,
    only_seqname: Optional[str] = None,
) -> dict[str, str]:
    """Single streaming pass: write each line to a per-seqname temp file.

    O(1) memory beyond small per-file write buffers -- safe for arbitrarily
    large input, because we never hold more than one line in memory at a
    time and never seek/re-read. Necessary because (empirically, see module
    docstring) the raw file is not seqname-sorted: a per-seqname flush-on-
    change single pass would be wrong.

    Parameters
    ----------
    input_path : str
        Raw Minimap2 GTF path.
    split_dir : str
        Directory to write ``<seqname>.gtf`` files into.
    only_seqname : str or None
        If given, only that seqname's lines are written (all others are
        skipped without opening a file for them) -- used for fast subset
        testing.

    Returns
    -------
    dict
        ``{seqname: filepath}`` for every seqname actually written.
    """
    ensure_dir(split_dir)
    handles: dict[str, object] = {}
    paths: dict[str, str] = {}
    n_lines = 0
    n_skipped_other_seqname = 0
    t0 = time.time()

    with open(input_path) as fh:
        for line in fh:
            n_lines += 1
            if n_lines % 20_000_000 == 0:
                print(
                    f"    ... split: {n_lines:,} lines read ({time.time() - t0:.0f}s)",
                    file=sys.stderr,
                )
            if not line or line.startswith("#"):
                continue
            tab = line.find("\t")
            if tab < 0:
                continue
            chrom = line[:tab]

            if only_seqname is not None and chrom != only_seqname:
                n_skipped_other_seqname += 1
                continue

            handle = handles.get(chrom)
            if handle is None:
                path = os.path.join(split_dir, f"{chrom}.gtf")
                handle = open(path, "w")
                handles[chrom] = handle
                paths[chrom] = path
            handle.write(line)

    for handle in handles.values():
        handle.close()

    print(
        f"    Split complete: {n_lines:,} lines -> {len(paths)} seqname file(s) "
        f"({time.time() - t0:.0f}s)"
        + (
            f", {n_skipped_other_seqname:,} lines skipped (not {only_seqname})"
            if only_seqname
            else ""
        )
    )
    return paths


# ---------------------------------------------------------------------------
# Stage 2: per-seqname consensus building
# ---------------------------------------------------------------------------


def _parse_attrs(attrs_str: str) -> dict[str, str]:
    return dict(ATTR_RE.findall(attrs_str))


def _intron_chain_str(exons: list[tuple[int, int]], snap_bp: int = 0) -> str:
    """Intron-chain signature for a sorted exon list (same convention as
    ``gmb.pipeline.scoring._get_intron_chain``, adapted for plain tuples).

    When ``snap_bp`` > 0, each intron boundary is rounded to the nearest
    multiple of ``snap_bp`` before building the signature, so reads whose
    splice calls differ by a few bp (real long-read alignment jitter, not a
    genuinely different splice choice) still group together.
    """
    if len(exons) < 2:
        return "single-exon"
    if snap_bp > 0:

        def snap(x):
            return round(x / snap_bp) * snap_bp

        return ",".join(
            f"{snap(exons[i][1])}-{snap(exons[i + 1][0])}" for i in range(len(exons) - 1)
        )
    return ",".join(f"{exons[i][1]}-{exons[i + 1][0]}" for i in range(len(exons) - 1))


def _load_transcripts_from_split_file(split_path: str) -> dict[str, dict]:
    """Reconstruct transcripts from a per-seqname split file's exon rows.

    Deliberately ignores the ``transcript`` feature line's own Start/End
    (empirically ~2.93% have swapped coordinates) -- every transcript's
    span is derived purely from its (clean) exon rows.
    """
    transcripts: dict[str, dict] = {}
    with open(split_path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) != 9:
                continue
            _chrom, _source, feat, start_s, end_s, _score, strand, _frame, attrs = f
            if feat != "exon":
                continue
            try:
                start, end = int(start_s), int(end_s)
            except ValueError:
                continue
            if start > end:
                continue

            attrs_d = _parse_attrs(attrs)
            tid = attrs_d.get("transcript_id")
            if not tid:
                continue
            gid = attrs_d.get("gene_id", tid)

            rec = transcripts.get(tid)
            if rec is None:
                rec = {"strand": strand, "gene_id": gid, "exons": []}
                transcripts[tid] = rec
            rec["exons"].append((start, end))
    return transcripts


def _cluster_single_exon(reads: list[dict], tol_bp: int) -> list[list[dict]]:
    """Single-linkage clustering of single-exon reads by proximity.

    Reads are merged into the same group if their start and end both fall
    within ``tol_bp`` of the group's running (min start, max end) -- this
    both collapses near-identical fragments and tolerates the minor
    terminal variation expected from independent reads of the same
    transcript.

    Reads are processed in start-sorted order, and any group whose min-start
    has fallen more than ``tol_bp`` behind the current read is retired (no
    later read, sorted ascending, could ever match it on start proximity
    again). This keeps the active-group set bounded by local read density
    rather than total read count -- an earlier all-groups-every-time version
    was O(n^2) and took hours on a single chromosome's worth of reads.
    """
    if not reads:
        return []
    reads = sorted(reads, key=lambda r: r["start"])

    finished: list[list[dict]] = []
    active: list[list[dict]] = []
    active_bounds: list[list[int]] = []  # [gs, ge], mutated in place

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


def _load_shortread_exons_for_seqname(
    paths: list[str], seqname: str
) -> list[tuple[int, int, str]]:
    """Load (start, end, strand) exon intervals for one seqname from small
    short-read GTFs (Scallop/StringTie). These files are small enough to
    read directly without pandas."""
    out = []
    for path in paths:
        if not path or not os.path.exists(path):
            continue
        with open(path) as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                f = line.rstrip("\n").split("\t")
                if len(f) != 9 or f[0] != seqname or f[2] != "exon":
                    continue
                try:
                    out.append((int(f[3]), int(f[4]), f[6]))
                except ValueError:
                    continue
    return out


def _has_shortread_support(
    exons: list[tuple[int, int]],
    strand: str,
    shortread_exons: list[tuple[int, int, str]],
    min_overlap: float,
) -> bool:
    for s_start, s_end, s_strand in shortread_exons:
        if s_strand != strand:
            continue
        for e_start, e_end in exons:
            if reciprocal_overlap(e_start, e_end, s_start, s_end) >= min_overlap:
                return True
    return False


def build_consensus_for_seqname(
    seqname: str,
    split_path: str,
    cfg: LongreadConsensusConfig,
    shortread_paths: Optional[list[str]] = None,
) -> tuple[list[dict], dict]:
    """Build consensus transcript models for one seqname's split file.

    Returns
    -------
    tuple of (list of consensus record dicts, stats dict)
    """
    stats: dict = {"seqname": seqname}

    transcripts = _load_transcripts_from_split_file(split_path)
    stats["input_reads"] = len(transcripts)

    # --- Per-read validation ---
    candidates = []
    n_no_exons = 0
    n_bad_intron = 0
    n_bad_span = 0
    for tid, rec in transcripts.items():
        exons = sorted(rec["exons"])
        if not exons:
            n_no_exons += 1
            continue
        span = exons[-1][1] - exons[0][0] + 1
        if span > cfg.max_span_length:
            n_bad_span += 1
            continue
        bad = False
        if len(exons) > 1:
            for i in range(len(exons) - 1):
                gap = exons[i + 1][0] - exons[i][1] - 1
                if gap < 0 or gap > cfg.max_intron_length:
                    bad = True
                    break
        if bad:
            n_bad_intron += 1
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

    stats["dropped_no_exons"] = n_no_exons
    stats["dropped_bad_intron_or_overlap"] = n_bad_intron
    stats["dropped_oversized_span"] = n_bad_span
    stats["surviving_after_basic_filters"] = len(candidates)

    shortread_exons = None
    if cfg.require_shortread_support_for_single_read_models and shortread_paths:
        shortread_exons = _load_shortread_exons_for_seqname(shortread_paths, seqname)
    stats["shortread_exons_loaded"] = len(shortread_exons) if shortread_exons else 0

    # --- Cluster by strand + genomic overlap (split on large gaps) ---
    consensus_records = []
    n_clusters = 0
    n_chain_groups_total = 0
    n_dropped_low_support = 0
    n_rescued_by_shortread = 0
    n_dropped_no_rescue = 0
    n_suppressed_fragment_of_multiexon = 0

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
                keep, rescued = _keep_or_rescue(
                    support, min_req, members, strand, shortread_exons, cfg
                )
                if not keep:
                    if support < min_req and support == 1:
                        n_dropped_no_rescue += 1
                    else:
                        n_dropped_low_support += 1
                    continue
                if rescued:
                    n_rescued_by_shortread += 1

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
                        "n_exons": len(exons_consensus),
                    }
                )

            for members in single_groups:
                group_id += 1

                # Suppress single-exon groups that overlap an already-accepted
                # multi-exon consensus at this same locus: these are almost
                # always truncated/partial-alignment fragments of that
                # transcript (see module docstring -- raw single-exon reads
                # are 63.3% of records vs a true ~45.6% single-exon gene
                # rate), not independent evidence, regardless of their own
                # read support or short-read rescue status.
                med_start_probe = sorted(m["start"] for m in members)[len(members) // 2]
                med_end_probe = sorted(m["end"] for m in members)[len(members) // 2]
                if any(
                    reciprocal_overlap(med_start_probe, med_end_probe, es, ee) >= 0.5
                    or (es <= med_start_probe and med_end_probe <= ee)
                    for es, ee in accepted_multiexon_exons
                ):
                    n_suppressed_fragment_of_multiexon += 1
                    continue

                support = len(members)
                min_req = cfg.min_read_support_single_exon
                keep, rescued = _keep_or_rescue(
                    support, min_req, members, strand, shortread_exons, cfg
                )
                if not keep:
                    if support < min_req and support == 1:
                        n_dropped_no_rescue += 1
                    else:
                        n_dropped_low_support += 1
                    continue
                if rescued:
                    n_rescued_by_shortread += 1

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
                        "n_exons": 1,
                    }
                )

    stats["clusters_found"] = n_clusters
    stats["chain_or_overlap_groups_found"] = n_chain_groups_total
    stats["dropped_low_support"] = n_dropped_low_support
    stats["dropped_single_read_no_rescue"] = n_dropped_no_rescue
    stats["suppressed_fragment_of_multiexon"] = n_suppressed_fragment_of_multiexon
    stats["rescued_by_shortread"] = n_rescued_by_shortread
    stats["consensus_models_before_containment_suppression"] = len(consensus_records)

    n_suppressed_contained = 0
    if cfg.suppress_contained_models:
        before = len(consensus_records)
        consensus_records = _suppress_contained_models(consensus_records)
        n_suppressed_contained = before - len(consensus_records)
    stats["suppressed_contained_in_higher_support_model"] = n_suppressed_contained
    stats["consensus_models_output"] = len(consensus_records)

    n_genes = _assign_gene_ids(consensus_records, seqname, cfg.min_intergenic_gap)
    stats["genes_output"] = n_genes

    return consensus_records, stats


def _suppress_contained_models(records: list[dict]) -> list[dict]:
    """Drop any model whose exon span is fully contained within a higher-
    read-support model's span on the same strand.

    Applied across *all* of a seqname's surviving consensus models (not
    restricted to one raw-read locus cluster, which -- at real long-read
    coverage depth -- collapses to ~1 cluster per strand for an entire
    chromosome and so is not a useful boundary here). This targets
    partial/truncated reads that share fewer introns than the full
    transcript and therefore never land in the same intron-chain group as
    it, and single-exon fragments landing inside a real exon that neither
    overlapped an accepted multi-exon model in their own cluster nor got
    caught by the earlier per-cluster suppression pass.

    O(n log n): sorted by start, a record can only be contained in a
    still-open higher-support interval, so the active set is bounded by
    local density rather than total record count.
    """
    kept: list[dict] = []
    for strand in ("+", "-"):
        strand_recs = [r for r in records if r["strand"] == strand]
        strand_recs.sort(key=lambda r: r["exons"][0][0])

        # Active intervals ordered by descending support so containment
        # checks prefer the strongest covering model.
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


def _assign_gene_ids(records: list[dict], seqname: str, min_intergenic_gap: int) -> int:
    """Reassign gene_id by genomic overlap across the *final* consensus
    records (mutates in place), replacing the coarse raw-read locus-cluster
    id used during grouping.

    The locus-clustering step upstream exists only to bound the search space
    for intron-chain/overlap grouping -- at real long-read coverage depth it
    collapses to ~1 cluster per strand for an entire chromosome (see module
    docstring), so it is not a usable proxy for gene identity. Genuine gene
    boundaries are only recoverable *after* raw reads have been collapsed
    into consensus transcripts, by clustering those (much smaller) surviving
    models on genomic overlap the same way the rest of the GMB pipeline
    treats isoforms of one gene as overlapping models on the same strand.

    Returns the number of distinct genes assigned.
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

    # transcript_id must stay unique even after gene_id collapses group
    # boundaries -- append a per-gene running index.
    per_gene_counter: dict[str, int] = defaultdict(int)
    for rec in records:
        per_gene_counter[rec["gene_id"]] += 1
        rec["transcript_id"] = f"{rec['gene_id']}.t{per_gene_counter[rec['gene_id']]}"

    return n_genes


def _keep_or_rescue(support, min_req, members, strand, shortread_exons, cfg):
    """Decide whether a candidate group survives, and whether via rescue."""
    if support >= min_req:
        return True, False
    if support == 1 and cfg.require_shortread_support_for_single_read_models and shortread_exons:
        exons = members[0]["exons"]
        if _has_shortread_support(exons, strand, shortread_exons, cfg.shortread_overlap_fraction):
            return True, True
    return False, False


def _median(values: list[int]) -> int:
    values = sorted(values)
    return values[len(values) // 2]


def _consensus_exons_multi(members: list[dict]) -> list[tuple[int, int]]:
    """Consensus exon coordinates for a group sharing one (possibly
    splice-site-tolerance-snapped) intron chain.

    All members have the same exon/intron *count* (grouping keys of
    different lengths never collide), but with splice-site snapping, exact
    internal junction coordinates can differ by up to the snap tolerance
    across members -- so every position (internal junctions too, not just
    the outer first-start/last-end) is collapsed to its per-position median
    across the group, rather than copying one arbitrary member's values.
    """
    n = len(sorted(members[0]["exons"]))
    exon_lists = [sorted(m["exons"]) for m in members]

    consensus = []
    for i in range(n):
        med_start = _median([ex[i][0] for ex in exon_lists])
        med_end = _median([ex[i][1] for ex in exon_lists])
        consensus.append((med_start, med_end))

    # Safety net: per-position medians are computed independently, so in a
    # rare heterogeneous-group edge case they could come out non-monotonic
    # or overlapping. Fall back to the first member's exact coordinates
    # rather than emit a structurally invalid transcript.
    for i in range(n - 1):
        if consensus[i][1] >= consensus[i + 1][0]:
            return exon_lists[0]
    return consensus


# ---------------------------------------------------------------------------
# Output writing
# ---------------------------------------------------------------------------


def write_consensus_gtf(
    records: list[dict], output_path: str, source_label: str = "Minimap2Consensus"
) -> None:
    """Write consensus records as a GTF, including explicit ``gene`` rows.

    ``gmb.compare.annotation_loader`` builds its gene table from
    ``Feature == "gene"`` rows only (it does not auto-derive genes from
    transcript/gene_id grouping for GTF input) -- so a gene-less transcript
    GTF loads as 0 genes when used as ``--query``. One gene row per distinct
    ``gene_id`` is emitted, spanning the union of its transcripts' exons.
    """
    by_gene: dict[str, list[dict]] = defaultdict(list)
    for rec in records:
        by_gene[rec["gene_id"]].append(rec)

    with open(output_path, "w") as fh:
        for gid, gene_records in by_gene.items():
            gene_start = min(r["exons"][0][0] for r in gene_records)
            gene_end = max(r["exons"][-1][1] for r in gene_records)
            strand = gene_records[0]["strand"]
            seqname = gene_records[0]["seqname"]
            fh.write(
                f"{seqname}\t{source_label}\tgene\t{gene_start}\t{gene_end}\t.\t"
                f'{strand}\t.\tgene_id "{gid}";\n'
            )
            for rec in gene_records:
                tid = rec["transcript_id"]
                exons = rec["exons"]
                start = exons[0][0]
                end = exons[-1][1]
                common_attrs = (
                    f'gene_id "{gid}"; transcript_id "{tid}"; '
                    f'read_support "{rec["read_support"]}"; '
                    f'cluster_size "{rec["cluster_size"]}"; '
                    f'rescued_by_shortread "{rec["rescued_by_shortread"]}";'
                )
                fh.write(
                    f"{seqname}\t{source_label}\ttranscript\t{start}\t{end}\t.\t"
                    f"{strand}\t.\t{common_attrs}\n"
                )
                for i, (s, e) in enumerate(exons, start=1):
                    fh.write(
                        f"{seqname}\t{source_label}\texon\t{s}\t{e}\t.\t"
                        f'{strand}\t.\tgene_id "{gid}"; transcript_id "{tid}"; '
                        f'exon_number "{i}"; read_support "{rec["read_support"]}";\n'
                    )


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Collapse raw Minimap2 long-read alignment GTFs into a small, "
            "provenance-preserving consensus transcript track safe to pass "
            "to `gmb.cli.build --minimap2`. Never loads the raw file "
            "wholesale -- always splits by seqname first."
        )
    )
    parser.add_argument("--input", required=True, help="Raw Minimap2 alignment GTF")
    parser.add_argument("--output-dir", required=True, help="Output directory")
    parser.add_argument(
        "--seqname", default=None, help="Restrict to a single seqname (fast subset test)"
    )
    parser.add_argument("--config", default=None, help="YAML overrides for consensus thresholds")
    parser.add_argument(
        "--scallop", default=None, help="Scallop GTF, for single-read rescue cross-check"
    )
    parser.add_argument(
        "--stringtie", default=None, help="StringTie GTF, for single-read rescue cross-check"
    )
    parser.add_argument(
        "--log-file", default=None, help="Log path. Default: <output-dir>/longread_consensus.log"
    )
    parser.add_argument("--no-log-file", action="store_true")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    ensure_dir(args.output_dir)
    log_file = resolve_log_file(args.output_dir, args.log_file, args.no_log_file)
    setup_logging(log_file=log_file, capture_stdio=log_file is not None)
    if log_file:
        print(f"Logging to {log_file}")

    if not os.path.exists(args.input):
        sys.exit(f"ERROR: input file not found: {args.input}")

    cfg = load_longread_consensus_config(args.config)

    print(f"Input: {args.input} ({os.path.getsize(args.input):,} bytes)")
    print("Config:", cfg)

    print("Stage 1: splitting raw file by seqname (single streaming pass)...")
    split_dir = os.path.join(args.output_dir, "_split_by_seqname")
    seqname_paths = split_by_seqname(args.input, split_dir, only_seqname=args.seqname)

    if not seqname_paths:
        sys.exit("ERROR: no records found (check --seqname matches the input file's naming).")

    shortread_paths = [p for p in (args.scallop, args.stringtie) if p]

    print(f"Stage 2: building consensus for {len(seqname_paths)} seqname(s)...")
    all_records = []
    all_stats = []
    for seqname in sorted(seqname_paths):
        t0 = time.time()
        records, stats = build_consensus_for_seqname(
            seqname, seqname_paths[seqname], cfg, shortread_paths=shortread_paths
        )
        stats["elapsed_s"] = round(time.time() - t0, 1)
        all_records.extend(records)
        all_stats.append(stats)
        print(
            f"  {seqname}: {stats['input_reads']:,} reads -> "
            f"{stats['consensus_models_output']:,} consensus models "
            f"({stats['elapsed_s']}s)"
        )

    output_gtf = os.path.join(args.output_dir, "minimap2_consensus.gtf")
    write_consensus_gtf(all_records, output_gtf)
    print(f"Wrote {len(all_records):,} consensus transcripts to {output_gtf}")

    summary = {
        "input_path": args.input,
        "input_size_bytes": os.path.getsize(args.input),
        "seqname_filter": args.seqname,
        "config": vars(cfg) if hasattr(cfg, "__dict__") else cfg.__dataclass_fields__,
        "per_seqname": all_stats,
        "totals": {
            k: sum(s.get(k, 0) for s in all_stats)
            for k in [
                "input_reads",
                "dropped_no_exons",
                "dropped_bad_intron_or_overlap",
                "dropped_oversized_span",
                "surviving_after_basic_filters",
                "clusters_found",
                "chain_or_overlap_groups_found",
                "dropped_low_support",
                "dropped_single_read_no_rescue",
                "suppressed_fragment_of_multiexon",
                "rescued_by_shortread",
                "consensus_models_before_containment_suppression",
                "suppressed_contained_in_higher_support_model",
                "consensus_models_output",
                "genes_output",
            ]
        },
    }
    summary_path = os.path.join(args.output_dir, "longread_consensus_summary.json")
    with open(summary_path, "w") as fh:
        json.dump(summary, fh, indent=2, default=str)
    print(f"Summary: {summary_path}")

    with open(os.path.join(args.output_dir, "longread_consensus_summary.tsv"), "w") as fh:
        fh.write("Metric\tValue\n")
        for k, v in summary["totals"].items():
            fh.write(f"{k}\t{v}\n")

    print("Done.")


if __name__ == "__main__":
    main()
