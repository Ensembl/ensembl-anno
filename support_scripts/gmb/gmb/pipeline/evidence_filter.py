#!/usr/bin/env python3
"""Pre-genebuild evidence filtering for Gene Model Builder.

Filters and denoises protein evidence (OrthoDB/UniProt), transcriptomic
assemblies, and Helixer predictions before consensus building.

All thresholds are driven by PipelineConfig (YAML).

DataFrame Schema Contract
-------------------------
Input DataFrames are exon-level rows sharing these columns:

* ``Chromosome`` : str -- sequence/contig name
* ``Start`` : int -- 0-based start (pyranges half-open convention)
* ``End`` : int -- 1-based end
* ``Strand`` : str -- ``"+"`` or ``"-"``
* ``transcript_id`` : str -- unique transcript identifier (source-prefixed)
* ``Source`` : str -- evidence origin, e.g. ``"OrthoDB"``, ``"Scallop"``
* ``Feature`` : str -- always ``"exon"`` for exon DataFrames

Protein DataFrames may additionally carry ``Score``, ``Coverage``,
and ``Identity`` attribute columns parsed from GFF3/GTF.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from gmb.pipeline.config import PipelineConfig

import pandas as pd
import pyranges as pr

from gmb.utils.intervals import reciprocal_overlap as _reciprocal_overlap

# ---------------------------------------------------------------------------
# Protein evidence filtering
# ---------------------------------------------------------------------------


def _transcript_span(df: pd.DataFrame) -> pd.DataFrame:
    """Collapse exon-level DataFrame to transcript-level spans."""
    return df.groupby("transcript_id", as_index=False).agg(
        Chromosome=("Chromosome", "first"),
        Start=("Start", "min"),
        End=("End", "max"),
        Strand=("Strand", "first"),
        exon_count=("Start", "count"),
    )


def filter_protein_evidence(
    protein_df: pd.DataFrame,
    config: PipelineConfig,
    stats: dict | None = None,
    transcriptomic_df: pd.DataFrame | None = None,
) -> pd.DataFrame:
    """Filter protein evidence (OrthoDB / UniProt GTF exon rows).

    Three-stage pipeline:
      0. Remove alignments failing basic thresholds (coverage, identity, bitscore)
      1. Remove short fragments (< min_protein_aa * 3 bp total CDS span)
      2. Collapse redundant mappings (same locus, high overlap)
      3. Locus competition (keep top-N per locus)

    Parameters
    ----------
    protein_df : pd.DataFrame
        Exon-level rows with Chromosome, Start, End, Strand, transcript_id,
        Source columns. May also contain Score, Coverage, Identity attributes.
    config : PipelineConfig
    stats : dict or None
        Mutable dict to record filtering stats.
    transcriptomic_df : pd.DataFrame or None
        If provided, used to boost protein models that overlap transcriptomic
        evidence.

    Returns
    -------
    pd.DataFrame
        Filtered exon-level rows.
    """
    if protein_df.empty:
        return protein_df

    pcfg = config.protein_filter
    if stats is None:
        stats = {}

    n_input_tx = protein_df["transcript_id"].nunique()
    stats["protein_input_transcripts"] = n_input_tx

    # --- Stage 0: Threshold filters (OrthoDB padding/noise) ---
    remove_tids = set()
    # Apply filters if columns exist (usually parsed from attributes)
    for tid, grp in protein_df.groupby("transcript_id"):
        # Usually these attributes are on the transcript-level row, but we might only have exon rows.
        # We take the max or mean across exons for that transcript, or simply the first available.
        if "Coverage" in grp.columns:
            # Handle string percentages
            cov_val = grp["Coverage"].iloc[0]
            if isinstance(cov_val, str):
                cov_val = (
                    float(cov_val.replace("%", "")) / 100.0 if "%" in cov_val else float(cov_val)
                )
            if cov_val < pcfg.min_alignment_coverage:
                remove_tids.add(tid)
                continue

        if "Identity" in grp.columns:
            id_val = grp["Identity"].iloc[0]
            if isinstance(id_val, str):
                id_val = float(id_val.replace("%", ""))
            if id_val < pcfg.min_percent_identity:
                remove_tids.add(tid)
                continue

        if "Score" in grp.columns:
            score_val = grp["Score"].iloc[0]
            if isinstance(score_val, str) and score_val != ".":
                if float(score_val) < pcfg.min_bitscore:
                    remove_tids.add(tid)
                    continue

    stats["protein_threshold_filtered"] = len(remove_tids)

    filtered_0 = protein_df[~protein_df["transcript_id"].isin(remove_tids)].copy()
    if filtered_0.empty:
        stats["protein_after_fragment_filter"] = 0
        stats["protein_after_redundancy"] = 0
        stats["protein_after_competition"] = 0
        return filtered_0

    # --- Stage 1: Remove short fragments ---
    min_span_bp = pcfg.min_protein_aa * 3
    spans = _transcript_span(filtered_0)
    spans["span"] = spans["End"] - spans["Start"]

    # Keep transcripts with sufficient span
    short_tids = set(spans[spans["span"] < min_span_bp]["transcript_id"])
    # Also remove absurdly long spans
    long_tids = set(spans[spans["span"] > pcfg.max_span_bp]["transcript_id"])
    remove_tids_1 = short_tids | long_tids

    stats["protein_short_fragments_removed"] = len(short_tids)
    stats["protein_long_artifacts_removed"] = len(long_tids)

    filtered = filtered_0[~filtered_0["transcript_id"].isin(remove_tids_1)].copy()

    if filtered.empty:
        stats["protein_after_fragment_filter"] = 0
        stats["protein_after_redundancy"] = 0
        stats["protein_after_competition"] = 0
        return filtered

    n_after_frag = filtered["transcript_id"].nunique()
    stats["protein_after_fragment_filter"] = n_after_frag

    # --- Stage 2: Collapse redundant mappings ---
    # Group by chromosome + strand, then cluster overlapping models
    spans2 = _transcript_span(filtered)
    spans2["span"] = spans2["End"] - spans2["Start"]
    spans2_pr = pr.PyRanges(spans2)

    try:
        clustered = spans2_pr.cluster(slack=0, count=True)
    except Exception:
        clustered = spans2_pr.cluster(count=True)

    cluster_col = "Cluster"
    cdf = clustered.df
    # Recompute span in case cluster() dropped it
    if "span" not in cdf.columns:
        cdf["span"] = cdf["End"] - cdf["Start"]
    if "exon_count" not in cdf.columns:
        cdf = cdf.merge(spans2[["transcript_id", "exon_count"]], on="transcript_id", how="left")

    # Within each cluster, collapse models with high reciprocal overlap
    keep_tids = set()
    redundant_removed = 0

    for _cid, grp in cdf.groupby(cluster_col):
        if len(grp) == 1:
            keep_tids.add(grp.iloc[0]["transcript_id"])
            continue

        # Sort by exon_count desc, span desc → best representative first
        grp = grp.sort_values(["exon_count", "span"], ascending=[False, False])
        representatives = []

        for _, row in grp.iterrows():
            tid = row["transcript_id"]
            merged = False
            for rep in representatives:
                ovl = _reciprocal_overlap(row["Start"], row["End"], rep["Start"], rep["End"])
                if ovl >= pcfg.redundancy_overlap:
                    merged = True
                    redundant_removed += 1
                    break
            if not merged:
                representatives.append(row)
                keep_tids.add(tid)

    stats["protein_redundant_collapsed"] = redundant_removed

    filtered = filtered[filtered["transcript_id"].isin(keep_tids)].copy()
    n_after_redundancy = filtered["transcript_id"].nunique()
    stats["protein_after_redundancy"] = n_after_redundancy

    # --- Stage 3: Locus competition ---
    if pcfg.top_n_per_locus > 0:
        spans3 = _transcript_span(filtered)
        spans3["span"] = spans3["End"] - spans3["Start"]
        spans3_pr = pr.PyRanges(spans3)

        try:
            clustered3 = spans3_pr.cluster(slack=0, count=True)
        except Exception:
            clustered3 = spans3_pr.cluster(count=True)

        cdf3 = clustered3.df
        # Recompute span/exon_count in case cluster() dropped them
        if "span" not in cdf3.columns:
            cdf3["span"] = cdf3["End"] - cdf3["Start"]
        if "exon_count" not in cdf3.columns:
            cdf3 = cdf3.merge(
                spans3[["transcript_id", "exon_count"]], on="transcript_id", how="left"
            )

        # Score each model: exon_count * 10 + span / 1000
        # Optionally boost if overlapping transcriptomic evidence
        tx_overlap_tids = set()
        if transcriptomic_df is not None and not transcriptomic_df.empty:
            tx_pr = pr.PyRanges(_transcript_span(transcriptomic_df))
            prot_pr = pr.PyRanges(spans3)
            ovl = prot_pr.overlap(tx_pr)
            if not ovl.df.empty:
                tx_overlap_tids = set(ovl.df["transcript_id"].unique())

        cdf3["score"] = (
            cdf3["exon_count"] * 10
            + cdf3["span"] / 1000
            + cdf3["transcript_id"].isin(tx_overlap_tids).astype(int) * 50
        )

        keep_tids_final = set()
        secondary_final = set()
        competition_removed = 0

        for _cid, grp in cdf3.groupby(cluster_col):
            grp = grp.sort_values("score", ascending=False)
            for i, (_, row) in enumerate(grp.iterrows()):
                if i < pcfg.top_n_per_locus:
                    keep_tids_final.add(row["transcript_id"])
                elif pcfg.keep_secondary:
                    secondary_final.add(row["transcript_id"])
                    keep_tids_final.add(row["transcript_id"])
                    competition_removed += 1
                else:
                    competition_removed += 1

        stats["protein_competition_demoted"] = competition_removed
        all_keep = keep_tids_final
        filtered = filtered[filtered["transcript_id"].isin(all_keep)].copy()

        # Mark secondary models
        if pcfg.keep_secondary and secondary_final:
            filtered["rank"] = filtered["transcript_id"].apply(
                lambda t: "secondary" if t in secondary_final else "primary"
            )
        else:
            filtered["rank"] = "primary"

    n_after_comp = filtered["transcript_id"].nunique()
    stats["protein_after_competition"] = n_after_comp

    print(
        f"    Protein evidence: {n_input_tx} → {n_after_comp} transcripts "
        f"({n_input_tx - n_after_comp} removed)"
    )
    print(
        f"      Fragments: -{stats.get('protein_short_fragments_removed', 0)}"
        f"  Redundant: -{redundant_removed}"
        f"  Competition: -{stats.get('protein_competition_demoted', 0)}"
    )

    return filtered


# ---------------------------------------------------------------------------
# Transcriptomic chimera / artefact filtering
# ---------------------------------------------------------------------------


def filter_chimeras(
    tx_df: pd.DataFrame,
    config: PipelineConfig,
    stats: dict | None = None,
) -> pd.DataFrame:
    """Filter likely chimeric or artefactual transcriptomic models.

    - Remove transcripts with unreasonably large introns
    - Keep single-exon transcripts by default (fungal mode)

    Parameters
    ----------
    tx_df : pd.DataFrame
        Exon-level transcriptomic rows.
    config : PipelineConfig
    stats : dict or None

    Returns
    -------
    pd.DataFrame
        Filtered exon rows.
    """
    if tx_df.empty:
        return tx_df

    tcfg = config.transcriptomic_filter
    if stats is None:
        stats = {}

    n_input = tx_df["transcript_id"].nunique()

    remove_tids = set()

    # Check intron lengths per transcript
    for tid, grp in tx_df.groupby("transcript_id"):
        exons = sorted(zip(grp["Start"].values, grp["End"].values))
        if len(exons) < 2:
            continue  # single-exon → always keep in fungal mode
        for i in range(len(exons) - 1):
            intron_len = exons[i + 1][0] - exons[i][1]
            if intron_len > tcfg.max_intron_length:
                remove_tids.add(tid)
                break

    stats["chimeras_large_intron"] = len(remove_tids)

    filtered = tx_df[~tx_df["transcript_id"].isin(remove_tids)].copy()
    n_after = filtered["transcript_id"].nunique()

    if remove_tids:
        print(
            f"    Chimera filter: {n_input} → {n_after} transcripts "
            f"({len(remove_tids)} removed for intron > "
            f"{tcfg.max_intron_length}bp)"
        )

    return filtered


# ---------------------------------------------------------------------------
# Mega-transcript splitting
# ---------------------------------------------------------------------------


def split_mega_transcripts(
    tx_df: pd.DataFrame,
    config: PipelineConfig,
    stats: dict | None = None,
) -> pd.DataFrame:
    """Split pathological mega-transcripts into local candidate segments.

    Transcriptome GTFs (especially StringTie --merge outputs) can contain
    transcript models spanning huge genomic regions.  This function splits
    them using exon geometry so downstream clustering/selection operates on
    local segments rather than dropping entire transcripts.

    Splitting is controlled by ``config.transcript_splitting``.

    Parameters
    ----------
    tx_df : pd.DataFrame
        Exon-level transcriptomic rows.
    config : PipelineConfig
    stats : dict or None

    Returns
    -------
    pd.DataFrame
        Exon rows with updated transcript_id / gene_id for segments,
        plus ``parent_transcript_id`` and ``parent_gene_id`` provenance columns.
    """
    scfg = config.transcript_splitting
    if stats is None:
        stats = {}

    if tx_df.empty or not scfg.split_enabled:
        return tx_df

    n_input = tx_df["transcript_id"].nunique()

    # QC: log top 20 worst offenders by span before splitting
    span_per_tid = tx_df.groupby("transcript_id").agg(
        span=("Start", lambda x: x.max()),
        span_end=("End", "max"),
    )
    span_per_tid["total_span"] = span_per_tid["span_end"] - span_per_tid["span"]
    top_offenders = span_per_tid.nlargest(20, "total_span")
    if not top_offenders.empty:
        print("    Transcript splitting — top 20 widest transcripts before split:")
        for tid, row in top_offenders.iterrows():
            print(f"      {tid}: {int(row['total_span']):,} bp")

    # Track stats
    transcripts_split = 0
    segments_emitted = 0
    max_segs_per_original = 0
    dropped_large_exon = 0
    dropped_max_segments = 0
    drop_reasons = []  # list of (tid, reason)

    result_rows = []

    for tid, grp in tx_df.groupby("transcript_id"):
        orig_gene_id = grp["gene_id"].iloc[0]

        # --- Large-exon check (before any splitting) ---
        if scfg.split_on_large_exon_bp is not None:
            exon_lens = grp["End"] - grp["Start"]
            if exon_lens.max() > scfg.split_on_large_exon_bp:
                dropped_large_exon += 1
                drop_reasons.append(
                    (
                        tid,
                        f"exon_exceeds_large_exon_limit "
                        f"(max_exon={int(exon_lens.max())} > "
                        f"{scfg.split_on_large_exon_bp})",
                    )
                )
                continue

        # --- Group by (Chromosome, Strand) ---
        contig_strand_groups = list(grp.groupby(["Chromosome", "Strand"]))

        segments = []  # list of DataFrames, one per segment

        for (_chrom, _strand), cs_grp in contig_strand_groups:
            # Sort exons by Start within this contig/strand
            cs_sorted = cs_grp.sort_values("Start")
            starts = cs_sorted["Start"].values
            ends = cs_sorted["End"].values

            # Segment by gap threshold
            seg_indices = [0]  # index into cs_sorted where new segments begin
            for i in range(1, len(starts)):
                gap_bp = max(0, int(starts[i]) - int(ends[i - 1]) - 1)
                if gap_bp > scfg.split_gap_bp:
                    seg_indices.append(i)

            # Emit segments
            for si in range(len(seg_indices)):
                start_idx = seg_indices[si]
                end_idx = seg_indices[si + 1] if si + 1 < len(seg_indices) else len(cs_sorted)
                seg_df = cs_sorted.iloc[start_idx:end_idx].copy()
                segments.append(seg_df)

        # --- Max-segments safety (drop, not cap) ---
        if len(segments) > scfg.max_segments_per_transcript:
            dropped_max_segments += 1
            drop_reasons.append(
                (
                    tid,
                    f"exceeded_max_segments "
                    f"({len(segments)} > {scfg.max_segments_per_transcript})",
                )
            )
            continue

        # --- Sort segments by genomic order (first exon Start) ---
        segments.sort(key=lambda s: s["Start"].iloc[0])

        if len(segments) == 1:
            # No actual splitting needed — pass through with provenance cols
            seg = segments[0].copy()
            seg["parent_transcript_id"] = tid
            seg["parent_gene_id"] = orig_gene_id
            result_rows.append(seg)
        else:
            # Splitting occurred
            transcripts_split += 1
            was_multi_contig = len(contig_strand_groups) > 1
            for seg_num, seg in enumerate(segments, start=1):
                seg = seg.copy()
                seg_label = f"__seg{seg_num:04d}"
                seg["parent_transcript_id"] = tid
                seg["parent_gene_id"] = orig_gene_id
                seg["transcript_id"] = f"{tid}{seg_label}"
                seg["gene_id"] = f"{orig_gene_id}{seg_label}"
                result_rows.append(seg)
            segments_emitted += len(segments)
            max_segs_per_original = max(max_segs_per_original, len(segments))

            if was_multi_contig:
                drop_reasons.append(
                    (tid, f"multi_contig_strand_split " f"(split into {len(segments)} segments)")
                )

    # --- Log drop/split reasons ---
    for tid, reason in drop_reasons:
        print(f"    Split QC: {tid} — {reason}")

    if result_rows:
        result_df = pd.concat(result_rows, ignore_index=True)
    else:
        result_df = tx_df.iloc[:0].copy()
        result_df["parent_transcript_id"] = pd.Series(dtype="object")
        result_df["parent_gene_id"] = pd.Series(dtype="object")

    n_output = result_df["transcript_id"].nunique()

    # Populate stats
    stats["transcripts_split"] = transcripts_split
    stats["segments_emitted"] = segments_emitted
    stats["max_segments_per_original"] = max_segs_per_original
    stats["transcripts_dropped_large_exon"] = dropped_large_exon
    stats["transcripts_dropped_max_segments"] = dropped_max_segments

    if transcripts_split > 0 or dropped_large_exon > 0 or dropped_max_segments > 0:
        print(
            f"    Transcript splitting: {n_input} → {n_output} transcripts "
            f"({transcripts_split} split, {segments_emitted} segments, "
            f"{dropped_large_exon} dropped-large-exon, "
            f"{dropped_max_segments} dropped-max-segments)"
        )

    return result_df


# ---------------------------------------------------------------------------
# Helixer model filtering
# ---------------------------------------------------------------------------


def filter_helixer_models(
    helixer_df: pd.DataFrame,
    helixer_cds: pd.DataFrame,
    config: PipelineConfig,
    stats: dict | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Filter implausible Helixer predictions.

    - Remove models with very short CDS (< min_cds_bp)
    - Flag models with extreme exon count

    Parameters
    ----------
    helixer_df : pd.DataFrame
        Exon-level Helixer rows.
    helixer_cds : pd.DataFrame
        CDS-level Helixer rows.
    config : PipelineConfig
    stats : dict or None

    Returns
    -------
    tuple of (pd.DataFrame, pd.DataFrame)
        ``(filtered_exons, filtered_cds)``.
    """
    if helixer_df.empty:
        return helixer_df, helixer_cds

    hcfg = config.helixer_filter
    if stats is None:
        stats = {}

    if not hcfg.enabled:
        return helixer_df, helixer_cds

    n_input = helixer_df["transcript_id"].nunique()

    remove_tids = set()

    # Check CDS length — use agg() for pandas 2.x / 3.0 compat
    if helixer_cds is not None and not helixer_cds.empty:
        cds_lengths = helixer_cds.groupby("transcript_id").agg(total_cds=("End", "sum")).copy()
        # Compute actual CDS bp per transcript
        cds_lengths["total_cds"] = (
            helixer_cds.groupby("transcript_id")["End"].sum()
            - helixer_cds.groupby("transcript_id")["Start"].sum()
        )
        short_cds = cds_lengths[cds_lengths["total_cds"] < hcfg.min_cds_bp].index
        remove_tids.update(short_cds)
        stats["helixer_short_cds_removed"] = len(short_cds)

    # Check exon count
    exon_counts = helixer_df.groupby("transcript_id").size()
    extreme = exon_counts[exon_counts > hcfg.max_exons].index
    stats["helixer_extreme_exon_flagged"] = len(extreme)
    # Flag but don't remove extreme exon models (they may be real)

    filtered_exons = helixer_df[~helixer_df["transcript_id"].isin(remove_tids)].copy()
    filtered_cds = helixer_cds
    if helixer_cds is not None and not helixer_cds.empty:
        filtered_cds = helixer_cds[~helixer_cds["transcript_id"].isin(remove_tids)].copy()

    n_after = filtered_exons["transcript_id"].nunique()
    removed = n_input - n_after
    if removed > 0:
        print(f"    Helixer filter: {n_input} → {n_after} transcripts " f"({removed} removed)")

    return filtered_exons, filtered_cds
