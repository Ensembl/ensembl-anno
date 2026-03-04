#!/usr/bin/env python3
"""Pre-genebuild evidence filtering for Gene Model Builder.

Filters and denoises protein evidence (OrthoDB/UniProt), transcriptomic
assemblies, and Helixer predictions before consensus building.

All thresholds are driven by PipelineConfig (YAML).
"""

from collections import defaultdict

import pandas as pd
import pyranges as pr


# ---------------------------------------------------------------------------
# Protein evidence filtering
# ---------------------------------------------------------------------------

def _transcript_span(df):
    """Collapse exon-level DataFrame to transcript-level spans."""
    return df.groupby('transcript_id', as_index=False).agg(
        Chromosome=('Chromosome', 'first'),
        Start=('Start', 'min'),
        End=('End', 'max'),
        Strand=('Strand', 'first'),
        exon_count=('Start', 'count'),
    )


def _reciprocal_overlap(s1, e1, s2, e2):
    """Compute reciprocal overlap fraction between two intervals."""
    overlap = max(0, min(e1, e2) - max(s1, s2))
    len1 = e1 - s1
    len2 = e2 - s2
    if len1 == 0 or len2 == 0:
        return 0.0
    return min(overlap / len1, overlap / len2)


def filter_protein_evidence(protein_df, config, stats=None,
                            transcriptomic_df=None):
    """Filter protein evidence (OrthoDB / UniProt GTF exon rows).

    Three-stage pipeline:
      1. Remove short fragments (< min_protein_aa * 3 bp total CDS span)
      2. Collapse redundant mappings (same locus, high overlap)
      3. Locus competition (keep top-N per locus)

    Parameters
    ----------
    protein_df : DataFrame
        Exon-level rows with Chromosome, Start, End, Strand, transcript_id,
        Source columns.
    config : PipelineConfig
    stats : dict or None
        Mutable dict to record filtering stats.
    transcriptomic_df : DataFrame or None
        If provided, used to boost protein models that overlap transcriptomic
        evidence.

    Returns
    -------
    DataFrame : filtered exon-level rows.
    """
    if protein_df.empty:
        return protein_df

    pcfg = config.protein_filter
    if stats is None:
        stats = {}

    n_input_tx = protein_df['transcript_id'].nunique()
    stats['protein_input_transcripts'] = n_input_tx

    # --- Stage 1: Remove short fragments ---
    min_span_bp = pcfg.min_protein_aa * 3
    spans = _transcript_span(protein_df)
    spans['span'] = spans['End'] - spans['Start']

    # Keep transcripts with sufficient span
    short_tids = set(spans[spans['span'] < min_span_bp]['transcript_id'])
    # Also remove absurdly long spans
    long_tids = set(
        spans[spans['span'] > pcfg.max_span_bp]['transcript_id'])
    remove_tids = short_tids | long_tids

    stats['protein_short_fragments_removed'] = len(short_tids)
    stats['protein_long_artifacts_removed'] = len(long_tids)

    filtered = protein_df[
        ~protein_df['transcript_id'].isin(remove_tids)].copy()

    if filtered.empty:
        stats['protein_after_fragment_filter'] = 0
        stats['protein_after_redundancy'] = 0
        stats['protein_after_competition'] = 0
        return filtered

    n_after_frag = filtered['transcript_id'].nunique()
    stats['protein_after_fragment_filter'] = n_after_frag

    # --- Stage 2: Collapse redundant mappings ---
    # Group by chromosome + strand, then cluster overlapping models
    spans2 = _transcript_span(filtered)
    spans2['span'] = spans2['End'] - spans2['Start']
    spans2_pr = pr.PyRanges(spans2)

    try:
        clustered = spans2_pr.cluster(slack=0, count=True)
    except Exception:
        clustered = spans2_pr.cluster(count=True)

    cluster_col = 'Cluster'
    cdf = clustered.df
    # Recompute span in case cluster() dropped it
    if 'span' not in cdf.columns:
        cdf['span'] = cdf['End'] - cdf['Start']
    if 'exon_count' not in cdf.columns:
        cdf = cdf.merge(
            spans2[['transcript_id', 'exon_count']],
            on='transcript_id', how='left')

    # Within each cluster, collapse models with high reciprocal overlap
    keep_tids = set()
    secondary_tids = set()
    redundant_removed = 0

    for _cid, grp in cdf.groupby(cluster_col):
        if len(grp) == 1:
            keep_tids.add(grp.iloc[0]['transcript_id'])
            continue

        # Sort by exon_count desc, span desc → best representative first
        grp = grp.sort_values(
            ['exon_count', 'span'], ascending=[False, False])
        representatives = []

        for _, row in grp.iterrows():
            tid = row['transcript_id']
            merged = False
            for rep in representatives:
                ovl = _reciprocal_overlap(
                    row['Start'], row['End'],
                    rep['Start'], rep['End'])
                if ovl >= pcfg.redundancy_overlap:
                    merged = True
                    redundant_removed += 1
                    break
            if not merged:
                representatives.append(row)
                keep_tids.add(tid)

    stats['protein_redundant_collapsed'] = redundant_removed

    filtered = filtered[filtered['transcript_id'].isin(keep_tids)].copy()
    n_after_redundancy = filtered['transcript_id'].nunique()
    stats['protein_after_redundancy'] = n_after_redundancy

    # --- Stage 3: Locus competition ---
    if pcfg.top_n_per_locus > 0:
        spans3 = _transcript_span(filtered)
        spans3['span'] = spans3['End'] - spans3['Start']
        spans3_pr = pr.PyRanges(spans3)

        try:
            clustered3 = spans3_pr.cluster(slack=0, count=True)
        except Exception:
            clustered3 = spans3_pr.cluster(count=True)

        cdf3 = clustered3.df
        # Recompute span/exon_count in case cluster() dropped them
        if 'span' not in cdf3.columns:
            cdf3['span'] = cdf3['End'] - cdf3['Start']
        if 'exon_count' not in cdf3.columns:
            cdf3 = cdf3.merge(
                spans3[['transcript_id', 'exon_count']],
                on='transcript_id', how='left')

        # Score each model: exon_count * 10 + span / 1000
        # Optionally boost if overlapping transcriptomic evidence
        tx_overlap_tids = set()
        if transcriptomic_df is not None and not transcriptomic_df.empty:
            tx_pr = pr.PyRanges(
                _transcript_span(transcriptomic_df))
            prot_pr = pr.PyRanges(spans3)
            ovl = prot_pr.overlap(tx_pr)
            if not ovl.df.empty:
                tx_overlap_tids = set(ovl.df['transcript_id'].unique())

        cdf3['score'] = (
            cdf3['exon_count'] * 10 +
            cdf3['span'] / 1000 +
            cdf3['transcript_id'].isin(tx_overlap_tids).astype(int) * 50
        )

        keep_tids_final = set()
        secondary_final = set()
        competition_removed = 0

        for _cid, grp in cdf3.groupby(cluster_col):
            grp = grp.sort_values('score', ascending=False)
            for i, (_, row) in enumerate(grp.iterrows()):
                if i < pcfg.top_n_per_locus:
                    keep_tids_final.add(row['transcript_id'])
                elif pcfg.keep_secondary:
                    secondary_final.add(row['transcript_id'])
                    keep_tids_final.add(row['transcript_id'])
                    competition_removed += 1
                else:
                    competition_removed += 1

        stats['protein_competition_demoted'] = competition_removed
        all_keep = keep_tids_final
        filtered = filtered[
            filtered['transcript_id'].isin(all_keep)].copy()

        # Mark secondary models
        if pcfg.keep_secondary and secondary_final:
            filtered['rank'] = filtered['transcript_id'].apply(
                lambda t: 'secondary' if t in secondary_final else 'primary')
        else:
            filtered['rank'] = 'primary'

    n_after_comp = filtered['transcript_id'].nunique()
    stats['protein_after_competition'] = n_after_comp

    print(f"    Protein evidence: {n_input_tx} → {n_after_comp} transcripts "
          f"({n_input_tx - n_after_comp} removed)")
    print(f"      Fragments: -{stats.get('protein_short_fragments_removed', 0)}"
          f"  Redundant: -{redundant_removed}"
          f"  Competition: -{stats.get('protein_competition_demoted', 0)}")

    return filtered


# ---------------------------------------------------------------------------
# Transcriptomic chimera / artefact filtering
# ---------------------------------------------------------------------------

def filter_chimeras(tx_df, config, stats=None):
    """Filter likely chimeric or artefactual transcriptomic models.

    - Remove transcripts with unreasonably large introns
    - Keep single-exon transcripts by default (fungal mode)

    Parameters
    ----------
    tx_df : DataFrame
        Exon-level transcriptomic rows.
    config : PipelineConfig
    stats : dict or None

    Returns
    -------
    DataFrame : filtered exon rows.
    """
    if tx_df.empty:
        return tx_df

    tcfg = config.transcriptomic_filter
    if stats is None:
        stats = {}

    n_input = tx_df['transcript_id'].nunique()

    remove_tids = set()

    # Check intron lengths per transcript
    for tid, grp in tx_df.groupby('transcript_id'):
        exons = sorted(zip(grp['Start'].values, grp['End'].values))
        if len(exons) < 2:
            continue  # single-exon → always keep in fungal mode
        for i in range(len(exons) - 1):
            intron_len = exons[i + 1][0] - exons[i][1]
            if intron_len > tcfg.max_intron_length:
                remove_tids.add(tid)
                break

    stats['chimeras_large_intron'] = len(remove_tids)

    filtered = tx_df[~tx_df['transcript_id'].isin(remove_tids)].copy()
    n_after = filtered['transcript_id'].nunique()

    if remove_tids:
        print(f"    Chimera filter: {n_input} → {n_after} transcripts "
              f"({len(remove_tids)} removed for intron > "
              f"{tcfg.max_intron_length}bp)")

    return filtered


# ---------------------------------------------------------------------------
# Helixer model filtering
# ---------------------------------------------------------------------------

def filter_helixer_models(helixer_df, helixer_cds, config, stats=None):
    """Filter implausible Helixer predictions.

    - Remove models with very short CDS (< min_cds_bp)
    - Flag models with extreme exon count

    Parameters
    ----------
    helixer_df : DataFrame
        Exon-level Helixer rows.
    helixer_cds : DataFrame
        CDS-level Helixer rows.
    config : PipelineConfig
    stats : dict or None

    Returns
    -------
    (filtered_exons, filtered_cds) : (DataFrame, DataFrame)
    """
    if helixer_df.empty:
        return helixer_df, helixer_cds

    hcfg = config.helixer_filter
    if stats is None:
        stats = {}

    if not hcfg.enabled:
        return helixer_df, helixer_cds

    n_input = helixer_df['transcript_id'].nunique()

    remove_tids = set()

    # Check CDS length
    if helixer_cds is not None and not helixer_cds.empty:
        cds_lengths = helixer_cds.groupby('transcript_id').apply(
            lambda g: (g['End'] - g['Start']).sum()
        )
        short_cds = cds_lengths[cds_lengths < hcfg.min_cds_bp].index
        remove_tids.update(short_cds)
        stats['helixer_short_cds_removed'] = len(short_cds)

    # Check exon count
    exon_counts = helixer_df.groupby('transcript_id').size()
    extreme = exon_counts[exon_counts > hcfg.max_exons].index
    stats['helixer_extreme_exon_flagged'] = len(extreme)
    # Flag but don't remove extreme exon models (they may be real)

    filtered_exons = helixer_df[
        ~helixer_df['transcript_id'].isin(remove_tids)].copy()
    filtered_cds = helixer_cds
    if helixer_cds is not None and not helixer_cds.empty:
        filtered_cds = helixer_cds[
            ~helixer_cds['transcript_id'].isin(remove_tids)].copy()

    n_after = filtered_exons['transcript_id'].nunique()
    removed = n_input - n_after
    if removed > 0:
        print(f"    Helixer filter: {n_input} → {n_after} transcripts "
              f"({removed} removed)")

    return filtered_exons, filtered_cds
