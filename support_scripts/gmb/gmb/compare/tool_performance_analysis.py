#!/usr/bin/env python3
"""Downstream tool performance analysis module.

Analyses outputs from multiple completed comparison runs (compare_annotations.py)
and answers biological/technical questions about DL annotation tool performance.

This module is SEPARATE from the core comparison tool. It reads existing
comparison outputs and optionally the reference annotation for stratification.

Usage:
    python -m gmb.compare.tool_performance_analysis \
      --comparisons helixer=compare_helixer annevo=compare_annevo tiberius=compare_tiberius \
      --reference-gff human_files/Homo_sapiens.GRCh38.115.gff3 \
      --evaluation-mode protein_coding \
      --reference-transcript-selection canonical \
      --output-dir tool_performance_analysis

Inputs:
    - comparison_details.tsv from each tool's comparison directory
    - comparison_summary.json from each tool's comparison directory
    - (optional) reference GFF3 for computing stratification features

Outputs:
    - tool_performance_summary.tsv: overall metrics per tool
    - tool_performance_summary.json: full structured results
    - tool_performance_strata.tsv: metrics broken down by feature strata
    - tool_best_matches.tsv: per-gene best match info with structural features
    - tool_error_taxonomy.tsv: error classification counts per tool
    - tool_intron_analysis.tsv: intron chain analysis per tool
    - tool_boundary_analysis.tsv: full-length recovery (start/end accuracy)
    - tool_cross_ranking.tsv: per-metric/stratum tool rankings
"""

import argparse
import json
import os
import sys
from typing import Optional

import numpy as np
import pandas as pd

from gmb.compare.annotation_loader import (
    apply_evaluation_mode,
    load_annotation,
    select_transcripts,
)

# ---------------------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------------------


def load_comparison_details(comparison_dir: str, tool_label: str) -> pd.DataFrame:
    """Load comparison_details.tsv and tag with tool label."""
    path = os.path.join(comparison_dir, "comparison_details.tsv")
    if not os.path.exists(path):
        raise FileNotFoundError(f"Missing {path}")
    df = pd.read_csv(path, sep="\t")
    df["tool"] = tool_label
    return df


def load_comparison_summary(comparison_dir: str) -> dict:
    """Load comparison_summary.json."""
    path = os.path.join(comparison_dir, "comparison_summary.json")
    if not os.path.exists(path):
        raise FileNotFoundError(f"Missing {path}")
    with open(path) as fh:
        return json.load(fh)


def compute_reference_features(
    ref_genes: pd.DataFrame,
    ref_mrna: pd.DataFrame,
    ref_exons: pd.DataFrame,
    ref_cds: pd.DataFrame,
) -> pd.DataFrame:
    """Compute structural features for each reference gene.

    Returns a DataFrame indexed by gene_id with columns:
      - gene_length, cds_length, exon_count, intron_count,
        transcript_count, is_single_exon
    """
    # Gene length
    gene_lengths = ref_genes.set_index("ID" if "ID" in ref_genes.columns else "gene_id")[
        ["Start", "End"]
    ].copy()
    gene_lengths["gene_length"] = gene_lengths["End"] - gene_lengths["Start"]

    # Exon count per gene (using best/representative transcript)
    # Group exons by gene_id
    if "gene_id" in ref_exons.columns:
        exon_counts = ref_exons.groupby("gene_id").size().rename("exon_count")
    else:
        exon_counts = pd.Series(dtype="int64", name="exon_count")

    # CDS length per gene (sum of CDS spans)
    if not ref_cds.empty and "gene_id" in ref_cds.columns:
        ref_cds_len = ref_cds.assign(_len=ref_cds["End"] - ref_cds["Start"])
        cds_lengths = ref_cds_len.groupby("gene_id")["_len"].sum().rename("cds_length")
    else:
        cds_lengths = pd.Series(dtype="int64", name="cds_length")

    # Transcript count per gene
    if not ref_mrna.empty and "gene_id" in ref_mrna.columns:
        tx_counts = (
            ref_mrna.groupby("gene_id")["transcript_id"].nunique().rename("transcript_count")
        )
    else:
        tx_counts = pd.Series(dtype="int64", name="transcript_count")

    # Intron count per gene (exon_count - 1 for single-transcript; approximate for multi)
    # Better: count unique introns from exon boundaries per gene
    if "gene_id" in ref_exons.columns and "transcript_id" in ref_exons.columns:
        # Use first transcript per gene for intron count
        first_tx = (
            ref_mrna.groupby("gene_id")["transcript_id"].first()
            if not ref_mrna.empty
            else pd.Series(dtype="str")
        )
        first_tx_set = set(first_tx.values)
        first_tx_exons = ref_exons[ref_exons["transcript_id"].isin(first_tx_set)]
        exon_per_tx = first_tx_exons.groupby("gene_id").size()
        intron_counts = (exon_per_tx - 1).clip(lower=0).rename("intron_count")
    else:
        intron_counts = pd.Series(dtype="int64", name="intron_count")

    # Build gene_id index
    if "ID" in ref_genes.columns:
        gids = ref_genes["ID"]
    elif "gene_id" in ref_genes.columns:
        gids = ref_genes["gene_id"]
    else:
        gids = pd.Series(dtype="str")

    feat_df = pd.DataFrame({"gene_id": gids.values})
    feat_df = feat_df.set_index("gene_id")
    feat_df["gene_length"] = gene_lengths["gene_length"]
    feat_df["exon_count"] = exon_counts
    feat_df["cds_length"] = cds_lengths
    feat_df["transcript_count"] = tx_counts
    feat_df["intron_count"] = intron_counts
    feat_df = feat_df.fillna(0).astype(
        {
            "gene_length": "int64",
            "exon_count": "int64",
            "cds_length": "int64",
            "transcript_count": "int64",
            "intron_count": "int64",
        }
    )
    feat_df["is_single_exon"] = feat_df["exon_count"] <= 1

    return feat_df


# ---------------------------------------------------------------------------
# Stratification helpers
# ---------------------------------------------------------------------------


def assign_strata(feat_df: pd.DataFrame) -> pd.DataFrame:
    """Assign stratification bins to reference features DataFrame.

    Adds columns: exon_stratum, intron_stratum, gene_length_stratum,
                  cds_length_stratum, transcript_stratum
    """
    df = feat_df.copy()

    # Exon count strata
    bins_exon = [0, 1, 3, 10, float("inf")]
    labels_exon = ["1_exon", "2-3_exon", "4-10_exon", ">10_exon"]
    df["exon_stratum"] = pd.cut(df["exon_count"], bins=bins_exon, labels=labels_exon, right=True)

    # Intron complexity
    bins_intron = [-1, 0, 2, 9, float("inf")]
    labels_intron = ["0_introns", "1-2_introns", "3-9_introns", ">=10_introns"]
    df["intron_stratum"] = pd.cut(
        df["intron_count"], bins=bins_intron, labels=labels_intron, right=True
    )

    # Gene length quartiles
    q_gene = df["gene_length"].quantile([0.25, 0.5, 0.75])
    bins_gene = [0, q_gene[0.25], q_gene[0.5], q_gene[0.75], float("inf")]
    labels_gene = ["Q1_shortest", "Q2", "Q3", "Q4_longest"]
    df["gene_length_stratum"] = pd.cut(
        df["gene_length"], bins=bins_gene, labels=labels_gene, right=True, duplicates="drop"
    )

    # CDS length quartiles (only for genes with CDS)
    cds_mask = df["cds_length"] > 0
    if cds_mask.any():
        q_cds = df.loc[cds_mask, "cds_length"].quantile([0.25, 0.5, 0.75])
        bins_cds = [0, q_cds[0.25], q_cds[0.5], q_cds[0.75], float("inf")]
        labels_cds = ["Q1_shortest", "Q2", "Q3", "Q4_longest"]
        df["cds_length_stratum"] = pd.cut(
            df["cds_length"], bins=bins_cds, labels=labels_cds, right=True, duplicates="drop"
        )
    else:
        df["cds_length_stratum"] = pd.Categorical(["no_cds"] * len(df))

    # Transcript complexity
    bins_tx = [0, 1, 5, 20, float("inf")]
    labels_tx = ["mono_tx", "2-5_tx", "6-20_tx", ">20_tx"]
    df["transcript_stratum"] = pd.cut(
        df["transcript_count"], bins=bins_tx, labels=labels_tx, right=True
    )

    return df


# ---------------------------------------------------------------------------
# Core analyses
# ---------------------------------------------------------------------------


def compute_overall_profile(ref_details: pd.DataFrame) -> dict:
    """Compute overall tool metrics from reference-perspective details."""
    total = len(ref_details)
    if total == 0:
        return {}

    cls = ref_details["classification"].value_counts()
    cls_cds = ref_details["classification_cds"].value_counts()

    exact = cls.get("Exact_Match", 0)
    structural = cls.get("Structural_Mismatch", 0)
    partial = cls.get("Partial_Match", 0)
    missed = cls.get("Missed", 0)
    strand = cls.get("Strand_Mismatch", 0)

    cds_exact = cls_cds.get("Exact_Match", 0)
    cds_any = (
        total
        - cls_cds.get("Missed", 0)
        - cls_cds.get("Strand_Mismatch", 0)
        - cls_cds.get("No_CDS", 0)
    )

    detected = total - missed
    intron_chain = ref_details["intron_chain_match"].sum()

    return {
        "sampled_ref_genes": total,
        "locus_detection_rate": detected / total if total > 0 else 0,
        "exact_transcript_match_rate": exact / total if total > 0 else 0,
        "structural_mismatch_rate": structural / total if total > 0 else 0,
        "partial_match_rate": partial / total if total > 0 else 0,
        "missed_rate": missed / total if total > 0 else 0,
        "strand_mismatch_rate": strand / total if total > 0 else 0,
        "cds_exact_match_rate": cds_exact / total if total > 0 else 0,
        "cds_any_match_rate": cds_any / total if total > 0 else 0,
        "intron_chain_recovery_rate": intron_chain / detected if detected > 0 else 0,
    }


def compute_novel_rate(cons_details: pd.DataFrame) -> float:
    """Compute novel query gene rate from consensus-perspective details."""
    total = len(cons_details)
    if total == 0:
        return 0.0
    novel = (cons_details["classification"] == "Novel").sum()
    return novel / total


def compute_stratified_metrics(ref_details: pd.DataFrame, feat_df: pd.DataFrame) -> pd.DataFrame:
    """Compute metrics stratified by reference structural features.

    Returns a DataFrame with columns: tool, stratum_type, stratum, metric, value
    """
    strata_df = assign_strata(feat_df)

    # Merge features with details
    merged = ref_details.merge(strata_df, left_on="gene_id", right_index=True, how="left")

    stratum_cols = [
        "exon_stratum",
        "intron_stratum",
        "gene_length_stratum",
        "cds_length_stratum",
        "transcript_stratum",
    ]

    rows = []
    for scol in stratum_cols:
        stratum_type = scol.replace("_stratum", "")
        for stratum_val, group in merged.groupby(scol, observed=True):
            if len(group) == 0:
                continue
            metrics = _compute_group_metrics(group)
            for metric_name, metric_val in metrics.items():
                rows.append(
                    {
                        "tool": ref_details["tool"].iloc[0] if not ref_details.empty else "",
                        "stratum_type": stratum_type,
                        "stratum": str(stratum_val),
                        "n_genes": len(group),
                        "metric": metric_name,
                        "value": metric_val,
                    }
                )

    return pd.DataFrame(rows)


def _compute_group_metrics(group: pd.DataFrame) -> dict:
    """Compute standard metrics for a group of reference genes."""
    n = len(group)
    if n == 0:
        return {}
    cls = group["classification"].value_counts()
    missed = cls.get("Missed", 0)
    detected = n - missed
    exact = cls.get("Exact_Match", 0)
    ic_match = group["intron_chain_match"].sum()

    cls_cds = group["classification_cds"].value_counts()
    cds_exact = cls_cds.get("Exact_Match", 0)

    return {
        "locus_detection_rate": detected / n,
        "exact_match_rate": exact / n,
        "cds_exact_match_rate": cds_exact / n,
        "intron_chain_rate": ic_match / detected if detected > 0 else 0,
        "missed_rate": missed / n,
    }


# ---------------------------------------------------------------------------
# Full-length recovery (boundary analysis)
# ---------------------------------------------------------------------------


def compute_boundary_analysis(
    ref_details: pd.DataFrame,
    ref_genes: pd.DataFrame,
    query_genes_by_tool: dict[str, pd.DataFrame],
) -> pd.DataFrame:
    """Compute transcript start/end boundary accuracy.

    For each matched ref gene, compute the offset between ref and query boundaries.
    Report exact matches and matches within tolerance thresholds.
    """
    rows = []
    tool = ref_details["tool"].iloc[0] if not ref_details.empty else ""

    # Get query genes for this tool
    q_genes = query_genes_by_tool.get(tool)
    if q_genes is None or q_genes.empty:
        return pd.DataFrame()

    # Build lookup: query gene_id → (start, end)
    if "ID" in q_genes.columns:
        q_lookup = q_genes.set_index("ID")[["Start", "End"]].to_dict("index")
    elif "gene_id" in q_genes.columns:
        q_lookup = q_genes.set_index("gene_id")[["Start", "End"]].to_dict("index")
    else:
        return pd.DataFrame()

    # Build ref lookup
    if "ID" in ref_genes.columns:
        r_lookup = ref_genes.set_index("ID")[["Start", "End"]].to_dict("index")
    elif "gene_id" in ref_genes.columns:
        r_lookup = ref_genes.set_index("gene_id")[["Start", "End"]].to_dict("index")
    else:
        return pd.DataFrame()

    matched = ref_details[
        (ref_details["classification"] != "Missed")
        & (ref_details["classification"] != "Strand_Mismatch")
        & (ref_details["matched_id"] != "")
    ]

    start_offsets = []
    end_offsets = []

    for _, row in matched.iterrows():
        ref_info = r_lookup.get(row["gene_id"])
        q_info = q_lookup.get(row["matched_id"])
        if ref_info is None or q_info is None:
            continue
        start_offsets.append(abs(ref_info["Start"] - q_info["Start"]))
        end_offsets.append(abs(ref_info["End"] - q_info["End"]))

    if not start_offsets:
        return pd.DataFrame()

    n_matched = len(start_offsets)
    start_arr = np.array(start_offsets)
    end_arr = np.array(end_offsets)

    thresholds = [0, 10, 50, 100, 500]
    for thresh in thresholds:
        rows.append(
            {
                "tool": tool,
                "boundary": "start",
                "threshold_bp": thresh,
                "count_within": int((start_arr <= thresh).sum()),
                "rate_within": float((start_arr <= thresh).sum() / n_matched),
                "n_matched": n_matched,
            }
        )
        rows.append(
            {
                "tool": tool,
                "boundary": "end",
                "threshold_bp": thresh,
                "count_within": int((end_arr <= thresh).sum()),
                "rate_within": float((end_arr <= thresh).sum() / n_matched),
                "n_matched": n_matched,
            }
        )

    # Combined (both start and end within threshold)
    for thresh in thresholds:
        both = int(((start_arr <= thresh) & (end_arr <= thresh)).sum())
        rows.append(
            {
                "tool": tool,
                "boundary": "both",
                "threshold_bp": thresh,
                "count_within": both,
                "rate_within": float(both / n_matched),
                "n_matched": n_matched,
            }
        )

    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# Intron chain analysis
# ---------------------------------------------------------------------------


def compute_intron_chain_analysis(ref_details: pd.DataFrame, feat_df: pd.DataFrame) -> dict:
    """Classify intron chain outcomes for multi-exon reference genes."""
    tool = ref_details["tool"].iloc[0] if not ref_details.empty else ""

    # Only multi-exon genes
    multi_exon = feat_df[feat_df["exon_count"] > 1].index
    multi_ref = ref_details[ref_details["gene_id"].isin(multi_exon)]

    if multi_ref.empty:
        return {"tool": tool, "n_multi_exon": 0}

    n = len(multi_ref)
    ic_exact = multi_ref["intron_chain_match"].sum()

    # Genes not missed and not strand mismatch but no IC match
    detected = multi_ref[
        (multi_ref["classification"] != "Missed")
        & (multi_ref["classification"] != "Strand_Mismatch")
    ]
    n_detected = len(detected)
    ic_no_match = n_detected - int(detected["intron_chain_match"].sum())

    missed = (multi_ref["classification"] == "Missed").sum()
    strand_mm = (multi_ref["classification"] == "Strand_Mismatch").sum()

    return {
        "tool": tool,
        "n_multi_exon_ref": n,
        "exact_intron_chain": int(ic_exact),
        "exact_intron_chain_rate": ic_exact / n_detected if n_detected > 0 else 0,
        "partial_or_shifted": ic_no_match,
        "partial_or_shifted_rate": ic_no_match / n_detected if n_detected > 0 else 0,
        "missed": int(missed),
        "missed_rate": missed / n if n > 0 else 0,
        "strand_mismatch": int(strand_mm),
    }


# ---------------------------------------------------------------------------
# Error taxonomy
# ---------------------------------------------------------------------------


def compute_error_taxonomy(
    ref_details: pd.DataFrame,
    cons_details: pd.DataFrame,
) -> pd.DataFrame:
    """Classify errors into taxonomy categories.

    Categories:
      - exact_transcript: Exact_Match classification
      - exact_cds_utr_differs: CDS exact but exon classification differs
      - exact_intron_chain_boundary_differs: intron chain match but not exact
      - partial_cds: CDS partial match
      - partial_exon_overlap: exon overlap but no chain match
      - wrong_strand: Strand_Mismatch
      - fragmented: one ref gene overlaps multiple query genes
      - merged: one query gene overlaps multiple ref genes
      - missed: Missed classification
    """
    tool = ref_details["tool"].iloc[0] if not ref_details.empty else ""
    n = len(ref_details)

    # Basic classifications from ref perspective
    exact_tx = (ref_details["classification"] == "Exact_Match").sum()
    missed = (ref_details["classification"] == "Missed").sum()
    wrong_strand = (ref_details["classification"] == "Strand_Mismatch").sum()

    # CDS exact but exon differs
    cds_exact_exon_diff = (
        (ref_details["classification_cds"] == "Exact_Match")
        & (ref_details["classification"] != "Exact_Match")
    ).sum()

    # Intron chain match but not exact overall
    ic_not_exact = (
        (ref_details["intron_chain_match"] == True)
        & (ref_details["classification"] != "Exact_Match")
    ).sum()

    # Fragmentation: ref genes whose matched_id appears more than once
    matched_ref = ref_details[ref_details["matched_id"] != ""]
    fragmented_query_ids = matched_ref["matched_id"].value_counts()
    ref_fragmented = matched_ref[
        matched_ref["matched_id"].isin(fragmented_query_ids[fragmented_query_ids > 1].index)
    ]["gene_id"].nunique()

    # Merging: query genes whose matched_id appears more than once
    matched_cons = cons_details[cons_details["matched_id"] != ""]
    merged_ref_ids = matched_cons["matched_id"].value_counts()
    query_merged = matched_cons[
        matched_cons["matched_id"].isin(merged_ref_ids[merged_ref_ids > 1].index)
    ]["gene_id"].nunique()

    # Remaining: partial overlaps
    partial_cds = (
        (ref_details["classification_cds"].isin(["Partial_Match", "Structural_Mismatch"]))
        & (ref_details["classification"] != "Missed")
        & (ref_details["classification"] != "Strand_Mismatch")
    ).sum()

    taxonomy = {
        "tool": tool,
        "total_ref_genes": n,
        "exact_transcript": int(exact_tx),
        "exact_cds_utr_differs": int(cds_exact_exon_diff),
        "exact_intron_chain_boundary_differs": int(ic_not_exact),
        "partial_cds": int(partial_cds),
        "wrong_strand": int(wrong_strand),
        "fragmented_predictions": int(ref_fragmented),
        "merged_predictions": int(query_merged),
        "missed_gene": int(missed),
    }

    return taxonomy


# ---------------------------------------------------------------------------
# Cross-tool ranking
# ---------------------------------------------------------------------------


def compute_cross_ranking(
    all_strata: pd.DataFrame,
    all_profiles: dict[str, dict],
) -> pd.DataFrame:
    """Rank tools per metric and per stratum.

    Returns DataFrame with: metric, stratum_type, stratum, rank_1, rank_2, ...
    """
    rows = []

    # Overall rankings
    tools = list(all_profiles.keys())
    overall_metrics = [
        "locus_detection_rate",
        "exact_transcript_match_rate",
        "cds_exact_match_rate",
        "cds_any_match_rate",
        "intron_chain_recovery_rate",
        "missed_rate",
    ]
    for metric in overall_metrics:
        vals = {t: all_profiles[t].get(metric, 0) for t in tools}
        # For missed_rate, lower is better
        reverse = metric != "missed_rate"
        ranked = sorted(vals.items(), key=lambda x: x[1], reverse=reverse)
        rows.append(
            {
                "metric": metric,
                "stratum_type": "overall",
                "stratum": "all",
                "rank_1": ranked[0][0] if len(ranked) > 0 else "",
                "rank_1_value": ranked[0][1] if len(ranked) > 0 else 0,
                "rank_2": ranked[1][0] if len(ranked) > 1 else "",
                "rank_2_value": ranked[1][1] if len(ranked) > 1 else 0,
                "rank_3": ranked[2][0] if len(ranked) > 2 else "",
                "rank_3_value": ranked[2][1] if len(ranked) > 2 else 0,
            }
        )

    # Stratified rankings
    if not all_strata.empty:
        strata_metrics = [
            "locus_detection_rate",
            "exact_match_rate",
            "cds_exact_match_rate",
            "intron_chain_rate",
            "missed_rate",
        ]
        for metric in strata_metrics:
            metric_data = all_strata[all_strata["metric"] == metric]
            for (stype, sval), group in metric_data.groupby(["stratum_type", "stratum"]):
                vals = {row["tool"]: row["value"] for _, row in group.iterrows()}
                reverse = metric != "missed_rate"
                ranked = sorted(vals.items(), key=lambda x: x[1], reverse=reverse)
                rows.append(
                    {
                        "metric": metric,
                        "stratum_type": stype,
                        "stratum": sval,
                        "rank_1": ranked[0][0] if len(ranked) > 0 else "",
                        "rank_1_value": ranked[0][1] if len(ranked) > 0 else 0,
                        "rank_2": ranked[1][0] if len(ranked) > 1 else "",
                        "rank_2_value": ranked[1][1] if len(ranked) > 1 else 0,
                        "rank_3": ranked[2][0] if len(ranked) > 2 else "",
                        "rank_3_value": ranked[2][1] if len(ranked) > 2 else 0,
                    }
                )

    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# Main orchestration
# ---------------------------------------------------------------------------


def run_analysis(
    comparisons: dict[str, str],
    reference_gff: Optional[str] = None,
    evaluation_mode: str = "protein_coding",
    transcript_selection: str = "canonical",
    output_dir: str = "tool_performance_analysis",
):
    """Run the full downstream analysis pipeline.

    Args:
        comparisons: dict of {tool_label: comparison_directory}
        reference_gff: path to reference GFF3 (optional, for stratification)
        evaluation_mode: evaluation mode to apply to reference
        transcript_selection: transcript selection mode for reference
        output_dir: output directory for analysis results
    """
    os.makedirs(output_dir, exist_ok=True)

    # Load comparison details
    all_ref_details = {}
    all_cons_details = {}
    all_summaries = {}

    for tool_label, comp_dir in comparisons.items():
        print(f"Loading comparison results for {tool_label} from {comp_dir}...")
        details = load_comparison_details(comp_dir, tool_label)
        all_ref_details[tool_label] = details[details["source"] == "reference"].copy()
        all_cons_details[tool_label] = details[details["source"] == "consensus"].copy()
        all_summaries[tool_label] = load_comparison_summary(comp_dir)

    # Load and filter reference if provided
    feat_df = None
    ref_genes = None
    query_genes_by_tool = {}

    if reference_gff:
        print("\nLoading reference annotation for stratification...")
        r_genes, r_mrna, r_exons, r_cds = load_annotation(reference_gff, "Reference")

        if evaluation_mode != "all":
            print(f"  Applying evaluation mode: {evaluation_mode}")
            r_genes, r_mrna, r_exons, r_cds = apply_evaluation_mode(
                evaluation_mode, r_genes, r_mrna, r_exons, r_cds
            )

        if transcript_selection != "all":
            print(f"  Selecting transcripts: {transcript_selection}")
            r_genes, r_mrna, r_exons, r_cds = select_transcripts(
                r_genes,
                r_mrna,
                r_exons,
                r_cds,
                mode=transcript_selection,
            )

        print("  Computing reference structural features...")
        feat_df = compute_reference_features(r_genes, r_mrna, r_exons, r_cds)
        ref_genes = r_genes
        print(f"  {len(feat_df)} genes with features computed")

        # Load query annotations for boundary analysis
        for tool_label, _comp_dir in comparisons.items():
            summary = all_summaries[tool_label]
            query_path = summary.get("query_path")
            if query_path and os.path.exists(query_path):
                print(f"  Loading {tool_label} annotation for boundary analysis...")
                q_genes, _, _, _ = load_annotation(query_path, tool_label)
                query_genes_by_tool[tool_label] = q_genes

    # --- Analysis 1: Overall tool profiles ---
    print("\n--- Computing overall tool profiles ---")
    all_profiles = {}
    for tool_label in comparisons:
        profile = compute_overall_profile(all_ref_details[tool_label])
        profile["novel_rate"] = compute_novel_rate(all_cons_details[tool_label])
        all_profiles[tool_label] = profile
        print(
            f"  {tool_label}: exact={profile.get('exact_transcript_match_rate', 0):.1%}, "
            f"CDS={profile.get('cds_exact_match_rate', 0):.1%}, "
            f"IC={profile.get('intron_chain_recovery_rate', 0):.1%}, "
            f"missed={profile.get('missed_rate', 0):.1%}"
        )

    # Write summary
    summary_rows = []
    for tool_label, profile in all_profiles.items():
        row = {"tool": tool_label}
        row.update(profile)
        summary_rows.append(row)
    summary_df = pd.DataFrame(summary_rows)
    summary_df.to_csv(
        os.path.join(output_dir, "tool_performance_summary.tsv"), sep="\t", index=False
    )

    # --- Analysis 2: Stratified performance ---
    all_strata = pd.DataFrame()
    if feat_df is not None:
        print("\n--- Computing stratified performance ---")
        strata_parts = []
        for tool_label in comparisons:
            strata = compute_stratified_metrics(all_ref_details[tool_label], feat_df)
            strata_parts.append(strata)
        all_strata = pd.concat(strata_parts, ignore_index=True)
        all_strata.to_csv(
            os.path.join(output_dir, "tool_performance_strata.tsv"), sep="\t", index=False
        )
        print(f"  {len(all_strata)} stratum-metric rows written")

    # --- Analysis 3: Boundary analysis ---
    if ref_genes is not None and query_genes_by_tool:
        print("\n--- Computing boundary analysis ---")
        boundary_parts = []
        for tool_label in comparisons:
            if tool_label in query_genes_by_tool:
                bnd = compute_boundary_analysis(
                    all_ref_details[tool_label], ref_genes, query_genes_by_tool
                )
                if not bnd.empty:
                    boundary_parts.append(bnd)
        if boundary_parts:
            boundary_df = pd.concat(boundary_parts, ignore_index=True)
            boundary_df.to_csv(
                os.path.join(output_dir, "tool_boundary_analysis.tsv"), sep="\t", index=False
            )
            print(f"  Boundary analysis: {len(boundary_df)} rows")

    # --- Analysis 4: Intron chain analysis ---
    if feat_df is not None:
        print("\n--- Computing intron chain analysis ---")
        ic_results = []
        for tool_label in comparisons:
            ic = compute_intron_chain_analysis(all_ref_details[tool_label], feat_df)
            ic_results.append(ic)
            print(f"  {tool_label}: IC exact rate = {ic.get('exact_intron_chain_rate', 0):.1%}")
        ic_df = pd.DataFrame(ic_results)
        ic_df.to_csv(os.path.join(output_dir, "tool_intron_analysis.tsv"), sep="\t", index=False)

    # --- Analysis 5: Error taxonomy ---
    print("\n--- Computing error taxonomy ---")
    error_results = []
    for tool_label in comparisons:
        taxonomy = compute_error_taxonomy(
            all_ref_details[tool_label], all_cons_details[tool_label]
        )
        error_results.append(taxonomy)
        print(
            f"  {tool_label}: fragmented={taxonomy['fragmented_predictions']}, "
            f"merged={taxonomy['merged_predictions']}, missed={taxonomy['missed_gene']}"
        )
    error_df = pd.DataFrame(error_results)
    error_df.to_csv(os.path.join(output_dir, "tool_error_taxonomy.tsv"), sep="\t", index=False)

    # --- Analysis 6: Cross-tool ranking ---
    print("\n--- Computing cross-tool rankings ---")
    ranking_df = compute_cross_ranking(all_strata, all_profiles)
    ranking_df.to_csv(os.path.join(output_dir, "tool_cross_ranking.tsv"), sep="\t", index=False)
    print(f"  {len(ranking_df)} ranking rows")

    # --- Analysis 7: Best matches per gene ---
    print("\n--- Writing best match details ---")
    best_match_parts = []
    for tool_label in comparisons:
        ref_det = all_ref_details[tool_label].copy()
        if feat_df is not None:
            ref_det = ref_det.merge(
                assign_strata(feat_df)[
                    [
                        "exon_count",
                        "cds_length",
                        "gene_length",
                        "intron_count",
                        "exon_stratum",
                        "gene_length_stratum",
                    ]
                ],
                left_on="gene_id",
                right_index=True,
                how="left",
            )
        best_match_parts.append(ref_det)
    best_match_df = pd.concat(best_match_parts, ignore_index=True)
    best_match_df.to_csv(os.path.join(output_dir, "tool_best_matches.tsv"), sep="\t", index=False)

    # --- Write full JSON summary ---
    full_results = {
        "evaluation_mode": evaluation_mode,
        "transcript_selection": transcript_selection,
        "reference_gff": reference_gff,
        "tools": list(comparisons.keys()),
        "overall_profiles": all_profiles,
        "error_taxonomy": error_results,
        "intron_chain_analysis": ic_results if feat_df is not None else [],
    }
    with open(os.path.join(output_dir, "tool_performance_summary.json"), "w") as fh:
        json.dump(full_results, fh, indent=2, default=str)

    # --- Print final report ---
    _print_final_report(
        all_profiles, error_results, ic_results if feat_df is not None else [], ranking_df
    )

    print(f"\nAll outputs written to: {output_dir}/")
    return full_results


def _print_final_report(all_profiles, error_results, ic_results, ranking_df):
    """Print human-readable final report to stdout."""
    tools = list(all_profiles.keys())

    print(f"\n\n{'='*80}")
    print("TOOL PERFORMANCE ANALYSIS — FINAL REPORT")
    print(f"{'='*80}")

    # Overall profile table
    print(f"\n{'--- Overall Profile ---'}")
    header = f"{'Metric':<35}"
    for t in tools:
        header += f" {t:>12}"
    print(header)
    print("-" * len(header))

    metrics = [
        ("Locus detection rate", "locus_detection_rate"),
        ("Exact transcript match", "exact_transcript_match_rate"),
        ("CDS exact match", "cds_exact_match_rate"),
        ("CDS any match", "cds_any_match_rate"),
        ("Intron chain recovery", "intron_chain_recovery_rate"),
        ("Missed rate", "missed_rate"),
        ("Strand mismatch rate", "strand_mismatch_rate"),
        ("Novel query rate", "novel_rate"),
    ]
    for label, key in metrics:
        row = f"{label:<35}"
        for t in tools:
            v = all_profiles[t].get(key, 0)
            row += f" {v:>11.1%}"
        print(row)

    # Error taxonomy
    if error_results:
        print(f"\n{'--- Error Taxonomy ---'}")
        header = f"{'Category':<40}"
        for t in tools:
            header += f" {t:>12}"
        print(header)
        print("-" * len(header))

        err_keys = [
            ("Exact transcript", "exact_transcript"),
            ("Exact CDS, UTR differs", "exact_cds_utr_differs"),
            ("Exact intron chain, boundary diff", "exact_intron_chain_boundary_differs"),
            ("Partial CDS", "partial_cds"),
            ("Wrong strand", "wrong_strand"),
            ("Fragmented predictions", "fragmented_predictions"),
            ("Merged predictions", "merged_predictions"),
            ("Missed gene", "missed_gene"),
        ]
        err_lookup = {e["tool"]: e for e in error_results}
        for label, key in err_keys:
            row = f"{label:<40}"
            for t in tools:
                v = err_lookup.get(t, {}).get(key, 0)
                row += f" {v:>12}"
            print(row)

    # Key findings from rankings
    if not ranking_df.empty:
        print(f"\n{'--- Key Rankings (overall) ---'}")
        overall = ranking_df[ranking_df["stratum_type"] == "overall"]
        for _, r in overall.iterrows():
            metric = r["metric"].replace("_", " ")
            best = r["rank_1"]
            val = r["rank_1_value"]
            print(f"  Best {metric}: {best} ({val:.1%})")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def parse_args():
    parser = argparse.ArgumentParser(
        description="Downstream tool performance analysis from comparison outputs"
    )
    parser.add_argument(
        "--comparisons",
        nargs="+",
        required=True,
        metavar="TOOL=DIR",
        help="Tool comparisons in format tool_name=comparison_dir (e.g. helixer=compare_helixer)",
    )
    parser.add_argument(
        "--reference-gff",
        default=None,
        help="Path to reference GFF3 for computing stratification features",
    )
    parser.add_argument(
        "--evaluation-mode",
        choices=["all", "protein_coding", "cds_only", "canonical"],
        default="protein_coding",
        help="Evaluation mode applied to reference (default: protein_coding)",
    )
    parser.add_argument(
        "--reference-transcript-selection",
        choices=["all", "longest_cds", "canonical"],
        default="canonical",
        help="Transcript selection mode for reference (default: canonical)",
    )
    parser.add_argument(
        "--output-dir",
        default="tool_performance_analysis",
        help="Output directory (default: tool_performance_analysis)",
    )
    return parser.parse_args()


def main():
    args = parse_args()

    # Parse comparisons
    comparisons = {}
    for item in args.comparisons:
        if "=" not in item:
            print(f"ERROR: comparison must be in format TOOL=DIR, got: {item}", file=sys.stderr)
            sys.exit(1)
        tool, directory = item.split("=", 1)
        comparisons[tool] = directory

    run_analysis(
        comparisons=comparisons,
        reference_gff=args.reference_gff,
        evaluation_mode=args.evaluation_mode,
        transcript_selection=args.reference_transcript_selection,
        output_dir=args.output_dir,
    )


if __name__ == "__main__":
    main()
