#!/usr/bin/env python3
"""Reporting module for Gene Model Builder.

Generates summary reports (JSON and/or TSV) with stage counts,
locus statistics, and protein-evidence retention metrics.
"""

from __future__ import annotations

import csv
import json
import os
from collections import Counter
from typing import TYPE_CHECKING, Any

if TYPE_CHECKING:
    from config import PipelineConfig


def generate_report(
    stats: dict[str, Any],
    loci: list[tuple[Any, list[dict]]],
    output_dir: str,
    config: PipelineConfig,
) -> dict:
    """Generate summary reports.

    Parameters
    ----------
    stats : dict
        Accumulated statistics from all pipeline stages.
    loci : list of (cluster_id, list of model dicts)
        Final selected isoforms.
    output_dir : str
        Directory to write report files into.
    config : PipelineConfig

    Returns
    -------
    dict
        The assembled report dictionary.
    """
    os.makedirs(output_dir, exist_ok=True)

    # --- Compute locus-level statistics ---
    n_loci = len(loci)
    n_transcripts = sum(len(models) for _, models in loci)

    tx_per_gene = [len(models) for _, models in loci]
    tx_per_gene_counter = Counter(tx_per_gene)

    # Evidence breakdown
    evidence_counter = Counter()
    exon_counts = []

    for _, models in loci:
        for m in models:
            ev = m.get("combined_evidence", m.get("source", "unknown"))
            for s in ev.split(","):
                evidence_counter[s.strip()] += 1
            exon_counts.append(m.get("exon_count", 0))

    single_exon_genes = sum(
        1 for _, models in loci if all(m.get("exon_count", 0) == 1 for m in models)
    )

    # --- Build report dict ---
    report = {
        "pipeline": "GeneModelBuilder",
        "preset": config.preset,
        "summary": {
            "total_loci": n_loci,
            "total_transcripts": n_transcripts,
            "single_exon_genes": single_exon_genes,
            "multi_exon_genes": n_loci - single_exon_genes,
            "transcripts_per_gene_distribution": dict(tx_per_gene_counter),
            "evidence_sources": dict(evidence_counter),
        },
        "filtering": {},
        "export": {},
    }

    # Add filtering stats
    filter_keys = [
        "protein_input_transcripts",
        "protein_short_fragments_removed",
        "protein_long_artifacts_removed",
        "protein_after_fragment_filter",
        "protein_redundant_collapsed",
        "protein_after_redundancy",
        "protein_competition_demoted",
        "protein_after_competition",
        "chimeras_large_intron",
        "helixer_short_cds_removed",
        "helixer_extreme_exon_flagged",
    ]
    for key in filter_keys:
        if key in stats:
            report["filtering"][key] = stats[key]

    # Add export stats
    export_keys = [
        "total_transcripts",
        "with_cds",
        "with_protein",
        "partial_5",
        "partial_3",
    ]
    for key in export_keys:
        if key in stats:
            report["export"][key] = stats[key]

    # Add stage counts
    stage_keys = [
        "scallop_loaded",
        "stringtie_loaded",
        "transcriptomics_deduped",
        "bridging_removed",
        "helixer_loaded",
        "protein_loaded",
        "loci_count",
        "selected_transcripts",
    ]
    stage_data = {}
    for key in stage_keys:
        if key in stats:
            stage_data[key] = stats[key]
    if stage_data:
        report["stages"] = stage_data

    # --- Write JSON ---
    if "json" in config.reporting.formats:
        json_path = os.path.join(output_dir, "summary.json")
        with open(json_path, "w") as fh:
            json.dump(report, fh, indent=2)
        print(f"  Report: {json_path}")

    # --- Write TSV ---
    if "tsv" in config.reporting.formats:
        tsv_path = os.path.join(output_dir, "summary.tsv")
        with open(tsv_path, "w", newline="") as fh:
            writer = csv.writer(fh, delimiter="\t")
            writer.writerow(["metric", "value"])
            writer.writerow(["preset", config.preset])
            writer.writerow(["total_loci", n_loci])
            writer.writerow(["total_transcripts", n_transcripts])
            writer.writerow(["single_exon_genes", single_exon_genes])
            writer.writerow(["multi_exon_genes", n_loci - single_exon_genes])
            for key, val in sorted(report.get("filtering", {}).items()):
                writer.writerow([f"filter.{key}", val])
            for key, val in sorted(report.get("export", {}).items()):
                writer.writerow([f"export.{key}", val])
            for key, val in sorted(report.get("stages", {}).items()):
                writer.writerow([f"stage.{key}", val])
        print(f"  Report: {tsv_path}")

    return report
