#!/usr/bin/env python3
"""
Compare Annotations
===================
Locus-level comparison of a query annotation against a reference annotation.
Supports multiple input formats:
  - Ensembl GFF3 (ID=gene:...; Parent=transcript:...)
  - Helixer GFF3 (standard ID=/Parent=)
  - ANNEVO GFF3 (standard ID=/Parent=, confidence scores in CDS)
  - Tiberius hybrid (gene/transcript with bare attrs, exon/CDS with GTF-style)

Produces:
  - comparison_summary.json / .tsv
  - comparison_details.tsv (one row per gene with classification)
  - Representative locus plots per category (in qc/ subdirectory)

Usage:
    python compare_annotations.py \
        --reference Homo_sapiens.GRCh38.115.gff3.gz \
        --query helixer.gff3 \
        [--reference-fasta Homo_sapiens.GRCh38.dna.toplevel.fa.gz] \
        [--assembly-report assembly_report.txt] \
        [--seqname-map custom_mapping.tsv] \
        --output-dir validation/ \
        [--evidence-gtfs scallop.gtf stringtie.gtf orthodb.gtf uniprot.gtf] \
        [--helixer helixer_remapped.gff3] \
        [--plots-per-category 3] \
        [--sample-loci 50 --seed 42]

Options for Seqname Mapping:
    --assembly-report: Standard NCBI assembly_report.txt to map GenBank accessions to chromosome names.
    --seqname-map: Custom TSV or CSV with 2 columns: from_seqname, to_seqname.
                   Can have headers (from_seqname/to_seqname), ignores # comments.
                   Takes precedence over --assembly-report on collisions.
"""

import argparse
import csv
import json
import os
import sys
from collections import Counter, defaultdict

import matplotlib.patches as patches
import matplotlib.pyplot as plt
import pandas as pd
import pyranges as pr

from gmb.compare.annotation_loader import (
    apply_evaluation_mode,
    filter_by_biotype,
    generate_filter_audit,
    load_annotation,
    select_transcripts,
    validate_against_fasta,
)
from gmb.pipeline.annotate_cds_utrs import load_genome, reverse_complement, translate
from gmb.pipeline.subset_utils import (
    add_subset_args,
    build_mapping,
    remap_df_seqnames,
    resolve_subset_regions,
    subset_df_by_regions,
    write_subset_manifest,
)

# ---------------------------------------------------------------------------
# Seqname remapping — imported from subset_utils
# (load_assembly_mapping, load_seqname_map, remap_df_seqnames are now in
#  subset_utils.py; build_mapping merges them with correct precedence)
# ---------------------------------------------------------------------------


# ---------------------------------------------------------------------------
# Data loading — delegates to annotation_loader for format-agnostic parsing
# ---------------------------------------------------------------------------


def load_gff(path, source_label="Reference"):
    """Load annotation file, return (exons, cds, genes) DataFrames.

    Delegates to the unified annotation_loader for format detection and parsing.
    Preserved for backward compatibility with evidence track loading.
    """
    genes, _mrna, exons, cds = load_annotation(path, source_label)
    return exons, cds, genes


def load_consensus_genes(gff3_path):
    """Load annotation and extract gene/exon/CDS/mRNA data.

    Delegates to the unified annotation_loader for format detection and parsing.
    """
    genes, mrna, exons, cds = load_annotation(gff3_path, "Consensus")
    return genes, exons, cds, mrna


# ---------------------------------------------------------------------------
# Locus-level classification
# ---------------------------------------------------------------------------


def _gene_span(gene_row):
    """Extract (chrom, start, end, strand) from a gene row."""
    return (gene_row["Chromosome"], gene_row["Start"], gene_row["End"], gene_row["Strand"])


def _exon_intervals(exon_df, tid):
    """Get sorted exon intervals for a transcript."""
    rows = exon_df[exon_df["transcript_id"] == tid]
    if rows.empty:
        return []
    return sorted(zip(rows["Start"].values, rows["End"].values))


def _group_by_transcript(exons_df, cds_df=None):
    """Group exons and CDS by transcript_id.
    Returns: dict of {transcript_id: {'exons': [(s,e)...], 'cds': [(s,e)...]}}
    """
    tx_dict = defaultdict(lambda: {"exons": [], "cds": []})
    for _, row in exons_df.iterrows():
        tid = row["transcript_id"]
        tx_dict[tid]["exons"].append((row["Start"], row["End"]))
    if cds_df is not None and not cds_df.empty:
        for _, row in cds_df.iterrows():
            tid = row.get("transcript_id", row.get("Parent", ""))
            if not tid:
                continue
            tx_dict[tid]["cds"].append((row["Start"], row["End"]))

    for _tid, data in tx_dict.items():
        data["exons"] = sorted(data["exons"])
        data["cds"] = sorted(data["cds"])
    return tx_dict


def _intron_chain(exons):
    """Return intron chain signature from sorted exon list."""
    if len(exons) < 2:
        return "single-exon"
    return ",".join(f"{exons[i][1]}-{exons[i+1][0]}" for i in range(len(exons) - 1))


def _exon_overlap_fraction(exons_a, exons_b):
    """Compute fraction of exonic bases in A covered by B."""
    total_a = sum(e - s for s, e in exons_a)
    if total_a == 0:
        return 0.0

    # Flatten B into intervals
    covered = 0
    for sa, ea in exons_a:
        for sb, eb in exons_b:
            overlap = max(0, min(ea, eb) - max(sa, sb))
            covered += overlap

    return covered / total_a


def classify_locus_pairs(
    ref_genes, ref_exons, ref_cds, cons_genes, cons_exons, cons_cds, cons_mrna
):
    """Classify each reference gene and each consensus gene using
    transcript-to-transcript comparisons (isoform-aware).

    Returns:
        ref_results: list of dicts (one per ref gene)
        cons_results: list of dicts (one per consensus gene)
    """
    print("  Classifying loci...")

    # Build ref transcript intervals
    ref_tx = _group_by_transcript(ref_exons, ref_cds)
    # Build cons transcript intervals
    cons_tx = _group_by_transcript(cons_exons, cons_cds)

    # Build ref gene entries
    ref_entries = []
    for _, g in ref_genes.iterrows():
        gid = g.get("ID", g.get("gene_id", g.get("transcript_id", "")))
        # get all transcript ids for this gene
        mask = ref_exons["Parent"] == gid
        if "gene_id" in ref_exons.columns:
            mask = mask | (ref_exons["gene_id"] == gid)
        tids = ref_exons[mask]["transcript_id"].unique()
        # Fallback if no matching standard ids
        if len(tids) == 0:
            gene_span = (g["Chromosome"], g["Start"], g["End"])
            tids = ref_exons[
                (ref_exons["Chromosome"] == gene_span[0])
                & (ref_exons["Start"] >= gene_span[1])
                & (ref_exons["End"] <= gene_span[2])
            ]["transcript_id"].unique()

        ref_entries.append(
            {
                "gene_id": gid,
                "chrom": g["Chromosome"],
                "start": g["Start"],
                "end": g["End"],
                "strand": g["Strand"],
                "tids": list(tids),
            }
        )

    # Build consensus gene entries
    cons_entries = []
    for _, g in cons_genes.iterrows():
        gid = g.get("ID", g.get("gene_id", g.get("transcript_id", "")))
        # Find the mRNA rows for this gene to get evidence info
        evidence = ""
        if cons_mrna is not None and not cons_mrna.empty:
            mrna_rows = cons_mrna[
                (
                    cons_mrna["transcript_id"].str.startswith(gid.split(".")[0])
                    if isinstance(gid, str)
                    else cons_mrna["transcript_id"] == gid
                )
            ]
            # Try to get evidence from attributes
            if not mrna_rows.empty:
                for attr_col in ["Evidence", "evidence"]:
                    if attr_col in mrna_rows.columns:
                        evidence = mrna_rows.iloc[0].get(attr_col, "")
                        break

        # get all transcript ids for this consensus gene
        mask = cons_exons["Parent"] == gid
        if "gene_id" in cons_exons.columns:
            mask = mask | (cons_exons["gene_id"] == gid)
        tids = cons_exons[mask]["transcript_id"].unique()
        if len(tids) == 0:
            gene_span = (g["Chromosome"], g["Start"], g["End"])
            tids = cons_exons[
                (cons_exons["Chromosome"] == gene_span[0])
                & (cons_exons["Start"] >= gene_span[1])
                & (cons_exons["End"] <= gene_span[2])
            ]["transcript_id"].unique()

        cons_entries.append(
            {
                "gene_id": gid,
                "chrom": g["Chromosome"],
                "start": g["Start"],
                "end": g["End"],
                "strand": g["Strand"],
                "evidence": evidence,
                "tids": list(tids),
            }
        )

    # Build PyRanges for overlap detection
    ref_pr = pr.PyRanges(
        pd.DataFrame(
            {
                "Chromosome": [e["chrom"] for e in ref_entries],
                "Start": [e["start"] for e in ref_entries],
                "End": [e["end"] for e in ref_entries],
                "Strand": [e["strand"] for e in ref_entries],
                "gene_id": [e["gene_id"] for e in ref_entries],
            }
        )
    )

    cons_pr = pr.PyRanges(
        pd.DataFrame(
            {
                "Chromosome": [e["chrom"] for e in cons_entries],
                "Start": [e["start"] for e in cons_entries],
                "End": [e["end"] for e in cons_entries],
                "Strand": [e["strand"] for e in cons_entries],
                "gene_id": [e["gene_id"] for e in cons_entries],
            }
        )
    )

    # Find overlaps (unstranded first to catch strand mismatches)
    ovl_unstrand = ref_pr.join(cons_pr, strandedness=False, suffix="_cons")

    # Build lookup: ref_gid → list of (cons_gid, same_strand)
    ref_to_cons = defaultdict(list)
    cons_to_ref = defaultdict(list)

    if not ovl_unstrand.df.empty:
        odf = ovl_unstrand.df
        for _, row in odf.iterrows():
            rgid = row["gene_id"]
            cgid = row["gene_id_cons"]
            same_strand = row["Strand"] == row["Strand_cons"]
            ref_to_cons[rgid].append(
                {
                    "cons_id": cgid,
                    "same_strand": same_strand,
                }
            )
            cons_to_ref[cgid].append(
                {
                    "ref_id": rgid,
                    "same_strand": same_strand,
                }
            )

    # Classify reference genes
    ref_results = []
    for entry in ref_entries:
        gid = entry["gene_id"]
        matches = ref_to_cons.get(gid, [])

        if not matches:
            entry.update(
                {
                    "classification": "Missed",
                    "classification_cds": "Missed",
                    "matched_id": "",
                    "best_ref_transcript_id": "",
                    "exon_overlap": 0.0,
                    "intron_chain_match": False,
                    "cds_overlap": 0.0,
                    "cds_intron_chain_match": False,
                }
            )
            ref_results.append(entry)
            continue

        # Check best match
        same_strand_matches = [m for m in matches if m["same_strand"]]
        diff_strand_matches = [m for m in matches if not m["same_strand"]]

        if not same_strand_matches and diff_strand_matches:
            entry.update(
                {
                    "classification": "Strand_Mismatch",
                    "classification_cds": "Strand_Mismatch",
                    "matched_id": diff_strand_matches[0]["cons_id"],
                    "best_ref_transcript_id": "",
                    "exon_overlap": 0.0,
                    "intron_chain_match": False,
                    "cds_overlap": 0.0,
                    "cds_intron_chain_match": False,
                }
            )
            ref_results.append(entry)
            continue

        best_exon_class, best_cds_class = "Partial_Match", "Partial_Match"
        best_match_gid = same_strand_matches[0]["cons_id"]
        best_match_tid = ""
        best_exon_ovl, best_cds_ovl = 0.0, 0.0
        best_exon_ic, best_cds_ic = False, False
        best_priority = -1  # Higher is better: 3=exact, 2=high_ovl_only, 1=partial

        for m in same_strand_matches:
            cgid = m["cons_id"]
            cons_gene_entry = next((c for c in cons_entries if c["gene_id"] == cgid), None)
            if not cons_gene_entry:
                continue

            for r_tid in entry["tids"]:
                rt_data = ref_tx.get(r_tid, {"exons": [], "cds": []})
                rt_exons, rt_cds = rt_data["exons"], rt_data["cds"]
                if not rt_exons:
                    continue
                rt_chain = _intron_chain(rt_exons)
                rt_cds_chain = _intron_chain(rt_cds)

                for c_tid in cons_gene_entry["tids"]:
                    ct_data = cons_tx.get(c_tid, {"exons": [], "cds": []})
                    ct_exons, ct_cds = ct_data["exons"], ct_data["cds"]
                    if not ct_exons:
                        continue

                    # Exon metrics
                    ovl_ex_a = _exon_overlap_fraction(rt_exons, ct_exons)
                    ovl_ex_b = _exon_overlap_fraction(ct_exons, rt_exons)
                    recip_ex = min(ovl_ex_a, ovl_ex_b)
                    ic_ex_match = rt_chain == _intron_chain(ct_exons)

                    ex_class = "Partial_Match"
                    if recip_ex >= 0.8:
                        ex_class = "Exact_Match" if ic_ex_match else "Structural_Mismatch"

                    # Priority logic for choosing the representative transcript pair
                    priority = 0
                    if ic_ex_match and recip_ex >= 0.8:
                        priority = 3
                    elif recip_ex >= 0.8:
                        priority = 2
                    elif recip_ex > 0:
                        priority = 1

                    if priority > best_priority or (
                        priority == best_priority and recip_ex > best_exon_ovl
                    ):
                        best_priority = priority
                        best_exon_ovl = recip_ex
                        best_exon_ic = ic_ex_match
                        best_exon_class = ex_class
                        best_match_gid = cgid
                        best_match_tid = c_tid

                        # Evaluate CDS metrics for this pair
                        if rt_cds and ct_cds:
                            ovl_cds_a = _exon_overlap_fraction(rt_cds, ct_cds)
                            ovl_cds_b = _exon_overlap_fraction(ct_cds, rt_cds)
                            recip_cds = min(ovl_cds_a, ovl_cds_b)
                            ic_cds_match = rt_cds_chain == _intron_chain(ct_cds)

                            cds_class = "Partial_Match"
                            if recip_cds >= 0.8:
                                cds_class = (
                                    "Exact_Match" if ic_cds_match else "Structural_Mismatch"
                                )

                            best_cds_ovl = recip_cds
                            best_cds_ic = ic_cds_match
                            best_cds_class = cds_class
                        else:
                            best_cds_ovl = 0.0
                            best_cds_ic = False
                            best_cds_class = "No_CDS" if not rt_cds else "Missed"

        entry.update(
            {
                "classification": best_exon_class,
                "classification_cds": best_cds_class,
                "matched_id": best_match_gid,
                "best_ref_transcript_id": best_match_tid,  # Note: this is actually the consensus tid that matched the ref, so naming is slightly overloaded but consistent
                "exon_overlap": best_exon_ovl,
                "intron_chain_match": best_exon_ic,
                "cds_overlap": best_cds_ovl,
                "cds_intron_chain_match": best_cds_ic,
            }
        )
        ref_results.append(entry)

    # Classify consensus genes
    cons_results = []
    _ref_matched_ids = {r["matched_id"] for r in ref_results if r["classification"] != "Missed"}

    for entry in cons_entries:
        gid = entry["gene_id"]
        matches = cons_to_ref.get(gid, [])

        if not matches:
            entry.update(
                {
                    "classification": "Novel",
                    "classification_cds": "Novel",
                    "matched_id": "",
                    "best_ref_transcript_id": "",
                    "exon_overlap": 0.0,
                    "intron_chain_match": False,
                    "cds_overlap": 0.0,
                    "cds_intron_chain_match": False,
                }
            )
            cons_results.append(entry)
            continue

        # Find the best matching ref based on transcripts
        same_strand = [m for m in matches if m["same_strand"]]
        if not same_strand:
            entry.update(
                {
                    "classification": "Strand_Mismatch",
                    "classification_cds": "Strand_Mismatch",
                    "matched_id": matches[0]["ref_id"],
                    "best_ref_transcript_id": "",
                    "exon_overlap": 0.0,
                    "intron_chain_match": False,
                    "cds_overlap": 0.0,
                    "cds_intron_chain_match": False,
                }
            )
            cons_results.append(entry)
            continue

        best_exon_class, best_cds_class = "Partial_Match", "Partial_Match"
        best_match_gid = same_strand[0]["ref_id"]
        best_match_tid = ""
        best_exon_ovl, best_cds_ovl = 0.0, 0.0
        best_exon_ic, best_cds_ic = False, False
        best_priority = -1

        for m in same_strand:
            rgid = m["ref_id"]
            ref_gene_entry = next((r for r in ref_entries if r["gene_id"] == rgid), None)
            if not ref_gene_entry:
                continue

            for c_tid in entry["tids"]:
                ct_data = cons_tx.get(c_tid, {"exons": [], "cds": []})
                ct_exons, ct_cds = ct_data["exons"], ct_data["cds"]
                if not ct_exons:
                    continue
                ct_chain = _intron_chain(ct_exons)
                ct_cds_chain = _intron_chain(ct_cds)

                for r_tid in ref_gene_entry["tids"]:
                    rt_data = ref_tx.get(r_tid, {"exons": [], "cds": []})
                    rt_exons, rt_cds = rt_data["exons"], rt_data["cds"]
                    if not rt_exons:
                        continue

                    # Exon metrics
                    ovl_ex_a = _exon_overlap_fraction(ct_exons, rt_exons)
                    ovl_ex_b = _exon_overlap_fraction(rt_exons, ct_exons)
                    recip_ex = min(ovl_ex_a, ovl_ex_b)
                    ic_ex_match = ct_chain == _intron_chain(rt_exons)

                    ex_class = "Partial_Match"
                    if recip_ex >= 0.8:
                        ex_class = "Exact_Match" if ic_ex_match else "Structural_Mismatch"

                    priority = 0
                    if ic_ex_match and recip_ex >= 0.8:
                        priority = 3
                    elif recip_ex >= 0.8:
                        priority = 2
                    elif recip_ex > 0:
                        priority = 1

                    if priority > best_priority or (
                        priority == best_priority and recip_ex > best_exon_ovl
                    ):
                        best_priority = priority
                        best_exon_ovl = recip_ex
                        best_exon_ic = ic_ex_match
                        best_exon_class = ex_class
                        best_match_gid = rgid
                        best_match_tid = r_tid

                        # CDS metrics
                        if ct_cds and rt_cds:
                            ovl_cds_a = _exon_overlap_fraction(ct_cds, rt_cds)
                            ovl_cds_b = _exon_overlap_fraction(rt_cds, ct_cds)
                            recip_cds = min(ovl_cds_a, ovl_cds_b)
                            ic_cds_match = ct_cds_chain == _intron_chain(rt_cds)

                            cds_class = "Partial_Match"
                            if recip_cds >= 0.8:
                                cds_class = (
                                    "Exact_Match" if ic_cds_match else "Structural_Mismatch"
                                )

                            best_cds_ovl, best_cds_ic, best_cds_class = (
                                recip_cds,
                                ic_cds_match,
                                cds_class,
                            )
                        else:
                            best_cds_ovl, best_cds_ic, best_cds_class = (
                                0.0,
                                False,
                                "No_CDS" if not ct_cds else "Missed",
                            )

        entry.update(
            {
                "classification": best_exon_class if best_priority > -1 else "Matched",
                "classification_cds": best_cds_class,
                "matched_id": best_match_gid,
                "best_ref_transcript_id": best_match_tid,
                "exon_overlap": best_exon_ovl,
                "intron_chain_match": best_exon_ic,
                "cds_overlap": best_cds_ovl,
                "cds_intron_chain_match": best_cds_ic,
            }
        )
        cons_results.append(entry)

    return ref_results, cons_results


# ---------------------------------------------------------------------------
# Summary & detailed output
# ---------------------------------------------------------------------------


def _stratified_counts(ref_results, ref_exons_df=None):
    """Compute stratified breakdowns of reference classification.

    Returns dict with keys: single_exon, multi_exon, cds_containing, non_cds.
    Each value is a Counter of classification labels.
    """
    strata = {
        "single_exon": Counter(),
        "multi_exon": Counter(),
        "cds_containing": Counter(),
        "non_cds": Counter(),
    }

    # Real per-transcript exon counts, when available. A gene is multi-exon
    # if any of its transcripts has >1 exon — genuine intron content, not
    # alternative-splicing isoform count (genes with a single isoform are
    # not necessarily single-exon, and vice versa).
    tx_exon_count = None
    if ref_exons_df is not None and not ref_exons_df.empty:
        tx_exon_count = ref_exons_df.groupby("transcript_id").size()

    for r in ref_results:
        cls = r["classification"]
        tids = r.get("tids", [])
        has_cds = r.get("classification_cds", "") not in ("No_CDS", "Missed", "")

        if tx_exon_count is not None and tids:
            max_exons = max((tx_exon_count.get(tid, 1) for tid in tids), default=1)
            is_single_exon = max_exons <= 1
        else:
            # Fallback when exon data isn't available: transcript-isoform
            # count is a poor proxy (conflates AS complexity with intron
            # content) but better than nothing.
            is_single_exon = len(tids) <= 1

        if is_single_exon:
            strata["single_exon"][cls] += 1
        else:
            strata["multi_exon"][cls] += 1

        if has_cds or r.get("cds_overlap", 0) > 0:
            strata["cds_containing"][cls] += 1
        else:
            strata["non_cds"][cls] += 1

    return strata


def _locus_detection_rate(ref_results):
    """Gene locus detection: any overlap on same strand counts as detected."""
    detected = sum(
        1
        for r in ref_results
        if r["classification"] in ("Exact_Match", "Partial_Match", "Structural_Mismatch")
    )
    total = len(ref_results)
    return detected, total, round(detected / total, 4) if total > 0 else 0.0


def _intron_chain_recovery_rate(ref_results):
    """Fraction of matched ref genes where intron chain was fully recovered."""
    matched = [
        r
        for r in ref_results
        if r["classification"] != "Missed" and r["classification"] != "Strand_Mismatch"
    ]
    if not matched:
        return 0, 0, 0.0
    recovered = sum(1 for r in matched if r.get("intron_chain_match", False))
    return recovered, len(matched), round(recovered / len(matched), 4)


def write_summary(ref_results, cons_results, output_dir, filter_info=None, ref_exons_df=None):
    """Write summary and detailed comparison files."""
    os.makedirs(output_dir, exist_ok=True)

    # Counts
    ref_class_counts = Counter(r["classification"] for r in ref_results)
    cons_class_counts = Counter(r["classification"] for r in cons_results)

    # Per-chromosome breakdown
    chr_breakdown = defaultdict(lambda: Counter())
    for r in ref_results:
        chr_breakdown[r["chrom"]][r["classification"]] += 1

    cds_exact_but_exon_differs = sum(
        1
        for r in ref_results
        if r["classification_cds"] == "Exact_Match" and r["classification"] != "Exact_Match"
    )

    # Locus-level and intron chain metrics
    locus_detected, locus_total, locus_rate = _locus_detection_rate(ref_results)
    ic_recovered, ic_matched, ic_rate = _intron_chain_recovery_rate(ref_results)

    # CDS intron chain recovery
    cds_ic_matched = [
        r
        for r in ref_results
        if r["classification_cds"] in ("Exact_Match", "Partial_Match", "Structural_Mismatch")
    ]
    cds_ic_recovered = sum(1 for r in cds_ic_matched if r.get("cds_intron_chain_match", False))

    # Stratified
    strata = _stratified_counts(ref_results, ref_exons_df=ref_exons_df)

    summary = {
        "total_reference_genes": len(ref_results),
        "total_consensus_genes": len(cons_results),
        "reference_classification": dict(ref_class_counts),
        "reference_classification_cds": dict(
            Counter(r["classification_cds"] for r in ref_results)
        ),
        "consensus_classification": dict(cons_class_counts),
        "per_chromosome": {c: dict(counts) for c, counts in sorted(chr_breakdown.items())},
        "sensitivity": {
            "exact_match_count": ref_class_counts.get("Exact_Match", 0),
            "partial_match_count": ref_class_counts.get("Partial_Match", 0),
            "any_match_count": (
                ref_class_counts.get("Exact_Match", 0)
                + ref_class_counts.get("Partial_Match", 0)
                + ref_class_counts.get("Structural_Mismatch", 0)
            ),
            "missed_count": ref_class_counts.get("Missed", 0),
            "strand_mismatch_count": ref_class_counts.get("Strand_Mismatch", 0),
            "locus_detected_count": locus_detected,
            "locus_detection_rate": locus_rate,
        },
        "sensitivity_cds": {
            "cds_exact_match_count": sum(
                1 for r in ref_results if r["classification_cds"] == "Exact_Match"
            ),
            "cds_any_match_count": sum(
                1
                for r in ref_results
                if r["classification_cds"]
                in ("Exact_Match", "Partial_Match", "Structural_Mismatch")
            ),
            "cds_missed_count": sum(1 for r in ref_results if r["classification_cds"] == "Missed"),
            "cds_exact_but_exon_differs": cds_exact_but_exon_differs,
            "cds_intron_chain_recovered": cds_ic_recovered,
            "cds_intron_chain_matched": len(cds_ic_matched),
        },
        "intron_chain": {
            "exon_intron_chain_recovered": ic_recovered,
            "exon_intron_chain_matched": ic_matched,
            "exon_intron_chain_rate": ic_rate,
        },
        "specificity": {
            "novel_consensus_count": cons_class_counts.get("Novel", 0),
            "matched_consensus_count": cons_class_counts.get("Matched", 0)
            + sum(cons_class_counts.values())
            - cons_class_counts.get("Novel", 0)
            - cons_class_counts.get("Strand_Mismatch", 0),
        },
        "stratified": {k: dict(v) for k, v in strata.items()},
    }

    if filter_info:
        summary["filter_info"] = filter_info

    # Compute rates
    total_ref = len(ref_results)
    if total_ref > 0:
        summary["sensitivity"]["exact_match_rate"] = round(
            ref_class_counts.get("Exact_Match", 0) / total_ref, 4
        )
        summary["sensitivity"]["any_match_rate"] = round(
            summary["sensitivity"]["any_match_count"] / total_ref, 4
        )
        summary["sensitivity"]["missed_rate"] = round(
            ref_class_counts.get("Missed", 0) / total_ref, 4
        )

        # CDS rates
        cds_total = (
            summary["sensitivity_cds"]["cds_any_match_count"]
            + summary["sensitivity_cds"]["cds_missed_count"]
        )
        if cds_total > 0:
            summary["sensitivity_cds"]["cds_exact_match_rate"] = round(
                summary["sensitivity_cds"]["cds_exact_match_count"] / total_ref, 4
            )
            summary["sensitivity_cds"]["cds_any_match_rate"] = round(
                summary["sensitivity_cds"]["cds_any_match_count"] / total_ref, 4
            )

    # JSON
    json_path = os.path.join(output_dir, "comparison_summary.json")
    with open(json_path, "w") as fh:
        json.dump(summary, fh, indent=2)
    print(f"  Summary: {json_path}")

    # TSV summary
    tsv_path = os.path.join(output_dir, "comparison_summary.tsv")
    with open(tsv_path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["metric", "value"])
        w.writerow(["total_reference_genes", len(ref_results)])
        w.writerow(["total_consensus_genes", len(cons_results)])
        for cls, cnt in sorted(ref_class_counts.items()):
            w.writerow([f"ref_{cls}", cnt])
        for cls, cnt in sorted(cons_class_counts.items()):
            w.writerow([f"cons_{cls}", cnt])
        for k, v in summary["sensitivity"].items():
            w.writerow([f"sens_{k}", v])
        for k, v in summary["sensitivity_cds"].items():
            w.writerow([f"sens_{k}", v])
        for k, v in summary["specificity"].items():
            w.writerow([f"spec_{k}", v])
    print(f"  Summary: {tsv_path}")

    # Detailed per-gene TSV
    detail_path = os.path.join(output_dir, "comparison_details.tsv")
    with open(detail_path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(
            [
                "source",
                "gene_id",
                "chrom",
                "start",
                "end",
                "strand",
                "classification",
                "classification_cds",
                "matched_id",
                "best_match_transcript_id",
                "exon_overlap",
                "intron_chain_match",
                "cds_overlap",
                "cds_intron_chain_match",
                "match_basis",
            ]
        )

        for r in ref_results:
            w.writerow(
                [
                    "reference",
                    r["gene_id"],
                    r["chrom"],
                    r["start"],
                    r["end"],
                    r["strand"],
                    r.get("classification", ""),
                    r.get("classification_cds", ""),
                    r.get("matched_id", ""),
                    r.get("best_ref_transcript_id", ""),
                    round(r.get("exon_overlap", 0), 4),
                    r.get("intron_chain_match", False),
                    round(r.get("cds_overlap", 0), 4),
                    r.get("cds_intron_chain_match", False),
                    "exon",
                ]
            )

        for c in cons_results:
            w.writerow(
                [
                    "consensus",
                    c["gene_id"],
                    c["chrom"],
                    c["start"],
                    c["end"],
                    c["strand"],
                    c.get("classification", ""),
                    c.get("classification_cds", ""),
                    c.get("matched_id", ""),
                    c.get("best_ref_transcript_id", ""),
                    round(c.get("exon_overlap", 0), 4),
                    c.get("intron_chain_match", False),
                    round(c.get("cds_overlap", 0), 4),
                    c.get("cds_intron_chain_match", False),
                    "exon",
                ]
            )
    print(f"  Details: {detail_path}")

    return summary


# ---------------------------------------------------------------------------
# Locus plotting
# ---------------------------------------------------------------------------

TRACK_COLORS = {
    "GenBank": "#D9534F",
    "Reference": "#D9534F",
    "Consensus": "#337AB7",
    "Query": "#337AB7",
    "Scallop": "#5CB85C",
    "StringTie": "#F0AD4E",
    "Minimap2": "#E67E22",
    "Helixer": "#9B59B6",
    "Tiberius": "#16A085",
    "OrthoDB": "#999999",
    "UniProt": "#888888",
    "GenBlast": "#7F8C8D",
}

CDS_HEIGHT = 0.4
UTR_HEIGHT = 0.18


def _draw_gene_track(ax, exon_df, cds_df, y, color, label, chrom, plot_start, plot_end, max_tx=5):
    """Draw exon/CDS blocks for one track at a locus."""
    if exon_df is None or exon_df.empty:
        return

    # Filter to region
    region = exon_df[
        (exon_df["Chromosome"] == chrom)
        & (exon_df["End"] > plot_start)
        & (exon_df["Start"] < plot_end)
    ]
    if region.empty:
        return

    # CDS for this region
    region_cds = pd.DataFrame()
    if cds_df is not None and not cds_df.empty:
        region_cds = cds_df[
            (cds_df["Chromosome"] == chrom)
            & (cds_df["End"] > plot_start)
            & (cds_df["Start"] < plot_end)
        ]

    # Build CDS lookup by transcript
    cds_by_tx = defaultdict(list)
    if not region_cds.empty:
        tid_col = "transcript_id" if "transcript_id" in region_cds.columns else "Parent"
        if tid_col in region_cds.columns:
            for _, row in region_cds.iterrows():
                cds_by_tx[row[tid_col]].append((row["Start"], row["End"]))

    tids = list(region["transcript_id"].unique())
    if len(tids) > max_tx:
        tids = tids[:max_tx]

    for tid in tids:
        tx_exons = region[region["transcript_id"] == tid].sort_values("Start")
        tx_cds = cds_by_tx.get(tid, [])

        for _, ex in tx_exons.iterrows():
            ex_s, ex_e = ex["Start"], ex["End"]

            if tx_cds:
                # Thin exon background
                rect = patches.Rectangle(
                    (ex_s, y - UTR_HEIGHT / 2),
                    ex_e - ex_s,
                    UTR_HEIGHT,
                    linewidth=0.5,
                    edgecolor=color,
                    facecolor=color,
                    alpha=0.3,
                )
                ax.add_patch(rect)

                # CDS portions
                for cs, ce in tx_cds:
                    os_ = max(cs, ex_s)
                    oe = min(ce, ex_e)
                    if os_ < oe:
                        rect = patches.Rectangle(
                            (os_, y - CDS_HEIGHT / 2),
                            oe - os_,
                            CDS_HEIGHT,
                            linewidth=1,
                            edgecolor=color,
                            facecolor=color,
                            alpha=0.7,
                        )
                        ax.add_patch(rect)
            else:
                # No CDS → thick exon block
                rect = patches.Rectangle(
                    (ex_s, y - CDS_HEIGHT / 2),
                    ex_e - ex_s,
                    CDS_HEIGHT,
                    linewidth=1,
                    edgecolor=color,
                    facecolor=color,
                    alpha=0.7,
                )
                ax.add_patch(rect)

        # Intron line
        if len(tx_exons) > 1:
            tx_s = tx_exons["Start"].min()
            tx_e = tx_exons["End"].max()
            ax.plot([tx_s, tx_e], [y, y], color=color, linewidth=1, zorder=0)


def plot_comparison_locus(
    chrom, start, end, category, ref_gid, cons_gid, evidence_tracks, output_path, genome=None
):
    """Plot a single comparison locus with all evidence tracks.

    evidence_tracks: dict of {label: (exon_df, cds_df)}
    """
    pad = max(500, int((end - start) * 0.3))
    plot_start = max(0, start - pad)
    plot_end = end + pad

    track_names = list(evidence_tracks.keys())
    n_tracks = len(track_names)

    fig, ax = plt.subplots(figsize=(16, max(4, 1.5 + n_tracks * 0.8)))

    y_positions = {name: i for i, name in enumerate(reversed(track_names))}

    ax.set_xlim(plot_start, plot_end)
    ax.set_ylim(-1, n_tracks)
    ax.set_yticks(list(y_positions.values()))
    ax.set_yticklabels(list(y_positions.keys()))

    title = f"{category}: {chrom}:{start}-{end}"
    if ref_gid:
        title += f"\nRef: {ref_gid}"
    if cons_gid:
        title += f"  Cons: {cons_gid}"
    ax.set_title(title, fontsize=10)
    ax.set_xlabel("Genomic Position")

    # Draw each track
    for name, (exon_df, cds_df) in evidence_tracks.items():
        color = TRACK_COLORS.get(name, "#666666")
        y = y_positions[name]
        _draw_gene_track(ax, exon_df, cds_df, y, color, name, chrom, plot_start, plot_end)

    # Highlight the locus region
    ax.axvspan(start, end, alpha=0.05, color="blue")

    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close()


def _write_qc_index_html(qc_dir, rows):
    """Write a small HTML index linking each QC plot to its locus metadata.

    `rows` is a list of dicts with keys: filename, category, gene_id,
    matched_id, chrom, start, end. This mirrors the fields written per-locus
    to comparison_details.tsv so a reader can cross-reference a plot with its
    row there (matched on gene_id).
    """
    by_category = defaultdict(list)
    for row in rows:
        by_category[row["category"]].append(row)

    html = [
        "<!doctype html><html><head><meta charset='utf-8'>",
        "<title>GMB comparison QC index</title>",
        "<style>",
        "body{font-family:sans-serif;margin:2em;}",
        "table{border-collapse:collapse;margin-bottom:2em;}",
        "td,th{border:1px solid #ccc;padding:6px 10px;text-align:left;font-size:0.9em;}",
        "img{max-width:480px;display:block;}",
        "h2{margin-top:2em;}",
        "</style></head><body>",
        "<h1>GMB comparison QC index</h1>",
        "<p>Each row's gene_id/matched_id can be cross-referenced against "
        "comparison_details.tsv for the full per-locus record.</p>",
    ]
    for category in sorted(by_category):
        html.append(f"<h2>{category} ({len(by_category[category])})</h2>")
        html.append(
            "<table><tr><th>Plot</th><th>Reference gene_id</th>"
            "<th>Matched query id</th><th>Locus</th></tr>"
        )
        for row in by_category[category]:
            html.append(
                "<tr><td><a href='{f}'><img src='{f}'></a></td>"
                "<td>{gid}</td><td>{mid}</td><td>{chrom}:{start}-{end}</td></tr>".format(
                    f=row["filename"],
                    gid=row["gene_id"],
                    mid=row.get("matched_id", "") or "-",
                    chrom=row["chrom"],
                    start=row["start"],
                    end=row["end"],
                )
            )
        html.append("</table>")
    html.append("</body></html>")

    with open(os.path.join(qc_dir, "index.html"), "w") as fh:
        fh.write("\n".join(html))


def generate_comparison_plots(
    ref_results, cons_results, evidence_tracks, output_dir, plots_per_category=3
):
    """Generate representative locus plots per classification category."""
    qc_dir = os.path.join(output_dir, "qc")
    os.makedirs(qc_dir, exist_ok=True)

    # Group ref results by classification
    by_class = defaultdict(list)
    for r in ref_results:
        by_class[r["classification"]].append(r)

    # Also add novel consensus
    for r in cons_results:
        if r["classification"] == "Novel":
            by_class["Novel"].append(r)

    plot_count = 0
    index_rows = []
    for category, entries in by_class.items():
        # Sort by chromosome + position for reproducibility
        entries.sort(key=lambda e: (e["chrom"], e["start"]))
        selected = entries[:plots_per_category]

        for i, entry in enumerate(selected):
            filename = f"{category.lower()}_{i+1}_{entry['chrom']}_{entry['start']}.png"
            filepath = os.path.join(qc_dir, filename)

            plot_comparison_locus(
                entry["chrom"],
                entry["start"],
                entry["end"],
                category,
                entry["gene_id"],
                entry.get("matched_id", ""),
                evidence_tracks,
                filepath,
            )
            plot_count += 1
            index_rows.append(
                {
                    "filename": filename,
                    "category": category,
                    "gene_id": entry["gene_id"],
                    "matched_id": entry.get("matched_id", ""),
                    "chrom": entry["chrom"],
                    "start": entry["start"],
                    "end": entry["end"],
                }
            )

    _write_qc_index_html(qc_dir, index_rows)
    print(f"  Generated {plot_count} comparison plots in {qc_dir}/ (see index.html)")


# ---------------------------------------------------------------------------
# BUSCO helper
# ---------------------------------------------------------------------------


def extract_genbank_proteins(ref_exons, ref_cds, genome, output_path, mapping=None):
    """Extract protein sequences from GenBank CDS coords.

    For BUSCO protein-mode comparison.
    """
    if ref_cds.empty:
        print("  Warning: No CDS features in reference — cannot extract proteins.")
        return 0

    if mapping:
        ref_cds = remap_df_seqnames(ref_cds, mapping)

    # Get unique transcript IDs with CDS
    tid_col = "transcript_id" if "transcript_id" in ref_cds.columns else "Parent"

    count = 0
    with open(output_path, "w") as fh:
        for tid, grp in ref_cds.groupby(tid_col):
            chrom = grp.iloc[0]["Chromosome"]
            strand = grp.iloc[0]["Strand"]

            if chrom not in genome:
                continue

            cds_intervals = sorted(zip(grp["Start"].values, grp["End"].values))

            # Build CDS nucleotide: concatenate in ascending genomic order, then
            # RC for minus strand.  RC(X+Y)==RC(Y)+RC(X) so ascending order
            # places the highest-coordinate (5') exon first after RC.
            cds_nuc = "".join(genome[chrom][cs:ce] for cs, ce in cds_intervals)
            if strand == "-":
                cds_nuc = reverse_complement(cds_nuc)

            protein = translate(cds_nuc)
            if protein.endswith("*"):
                protein = protein[:-1]

            if len(protein) > 10:  # Skip very short
                fh.write(f">{tid}\n")
                for j in range(0, len(protein), 80):
                    fh.write(protein[j : j + 80] + "\n")
                count += 1

    print(f"  Extracted {count} proteins to {output_path}")
    return count


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def main():
    parser = argparse.ArgumentParser(
        description="Compare a query annotation against a reference annotation. "
        "Supports Ensembl GFF3, Helixer GFF3, ANNEVO GFF3, and Tiberius hybrid formats."
    )
    parser.add_argument(
        "--reference", required=True, help="Reference annotation (GFF3/GTF, can be .gz)"
    )
    query_group = parser.add_mutually_exclusive_group(required=True)
    query_group.add_argument("--query", help="Query annotation to compare against reference")
    query_group.add_argument("--consensus", help="(deprecated, use --query) Consensus annotation")
    parser.add_argument(
        "--reference-fasta",
        default=None,
        help="Reference FASTA for coordinate validation (also used for protein extraction)",
    )
    parser.add_argument(
        "--assembly-report", default=None, help="NCBI assembly report for seqname remapping"
    )
    parser.add_argument(
        "--seqname-map",
        default=None,
        help="Custom TSV/CSV mapping for seqnames (from_seqname, to_seqname)",
    )
    parser.add_argument(
        "--genome",
        default=None,
        help="(deprecated, use --reference-fasta) Genome FASTA for protein extraction",
    )
    parser.add_argument("--output-dir", default="validation", help="Output directory")
    parser.add_argument(
        "--query-format",
        choices=["auto", "gff3", "gtf", "tiberius"],
        default="auto",
        help="Format hint for query annotation (default: auto-detect)",
    )
    parser.add_argument(
        "--reference-format",
        choices=["auto", "gff3", "gtf", "tiberius"],
        default="auto",
        help="Format hint for reference annotation (default: auto-detect)",
    )
    parser.add_argument(
        "--evaluation-mode",
        choices=["all", "protein_coding", "cds_only", "canonical"],
        default="all",
        help="Preset evaluation mode for reference filtering (default: all)",
    )
    parser.add_argument(
        "--reference-gene-biotypes",
        default=None,
        help="Comma-separated gene biotypes to keep in reference (e.g. protein_coding)",
    )
    parser.add_argument(
        "--reference-transcript-biotypes",
        default=None,
        help="Comma-separated transcript biotypes to keep in reference (e.g. protein_coding)",
    )
    parser.add_argument(
        "--reference-transcript-selection",
        choices=["all", "longest_cds", "canonical"],
        default="all",
        help="Transcript selection mode for reference (default: all)",
    )
    parser.add_argument("--scallop", default=None, help="Scallop GTF")
    parser.add_argument("--stringtie", default=None, help="StringTie GTF")
    parser.add_argument(
        "--minimap2", default=None, help="Minimap2 long-read transcript alignments (GTF)"
    )
    parser.add_argument("--helixer", default=None, help="Helixer GFF3 (ab initio backbone)")
    parser.add_argument(
        "--tiberius",
        default=None,
        help="Tiberius GTF evidence track (ab initio backbone, alternative to --helixer)",
    )
    parser.add_argument("--orthodb", default=None, help="OrthoDB GTF")
    parser.add_argument("--uniprot", default=None, help="UniProt GTF")
    parser.add_argument("--genblast", default=None, help="GenBlast protein alignment GTF")
    parser.add_argument(
        "--plots-per-category",
        type=int,
        default=3,
        help="Number of representative plots per category",
    )
    parser.add_argument(
        "--extract-ref-proteins",
        action="store_true",
        default=False,
        help="Extract protein FASTA from reference for BUSCO",
    )
    parser.add_argument(
        "--sample-from",
        choices=["reference", "consensus", "union"],
        default="union",
        help="Source of loci for --sample-loci " "(default: union)",
    )
    parser.add_argument(
        "--evidence-attribution",
        type=str,
        default=None,
        help="Path to evidence_attribution.tsv to label with comparison results",
    )
    add_subset_args(parser)
    args = parser.parse_args()

    # Resolve --consensus -> --query backward compat
    query_path = args.query or args.consensus
    # Resolve --genome -> --reference-fasta backward compat
    fasta_path = args.reference_fasta or args.genome

    os.makedirs(args.output_dir, exist_ok=True)

    # --- Load mappings ---
    mapping = build_mapping(
        assembly_report=args.assembly_report,
        seqname_map=getattr(args, "seqname_map", None),
    )

    # --- Load reference ---
    print("Loading reference annotation...")
    ref_genes, _ref_mrna, ref_exons, ref_cds = load_annotation(
        args.reference, "Reference", format_hint=args.reference_format
    )
    if mapping:
        ref_exons = remap_df_seqnames(ref_exons, mapping)
        ref_cds = remap_df_seqnames(ref_cds, mapping)
        ref_genes = remap_df_seqnames(ref_genes, mapping, label="Reference")

    # --- Load query ---
    print("Loading query annotation...")
    cons_genes, cons_mrna, cons_exons, cons_cds = load_annotation(
        query_path, "Query", format_hint=args.query_format
    )
    if mapping:
        cons_exons = remap_df_seqnames(cons_exons, mapping)
        cons_cds = remap_df_seqnames(cons_cds, mapping)
        cons_genes = remap_df_seqnames(cons_genes, mapping, label="Query")

    # --- FASTA validation ---
    if fasta_path:
        print("Validating coordinates against FASTA...")
        ref_diag = validate_against_fasta(ref_genes, ref_exons, ref_cds, fasta_path, "Reference")
        query_diag = validate_against_fasta(cons_genes, cons_exons, cons_cds, fasta_path, "Query")
        all_diag = ref_diag + query_diag
        if all_diag:
            print(f"  WARNING: {len(all_diag)} validation issue(s) detected.")
            for d in all_diag:
                print(f"    {d}")
        else:
            print("  All coordinates validated successfully.")

    # --- Reference filtering ---
    pre_filter_ref_genes = len(ref_genes)
    pre_filter_ref_tx = _ref_mrna["transcript_id"].nunique() if not _ref_mrna.empty else 0

    # Evaluation mode preset
    if args.evaluation_mode != "all":
        print(f"Applying evaluation mode: {args.evaluation_mode}")
        ref_genes, _ref_mrna, ref_exons, ref_cds = apply_evaluation_mode(
            args.evaluation_mode, ref_genes, _ref_mrna, ref_exons, ref_cds
        )

    # Explicit biotype filtering (overrides/stacks with evaluation mode)
    gene_biotypes = None
    tx_biotypes = None
    if args.reference_gene_biotypes:
        gene_biotypes = [b.strip() for b in args.reference_gene_biotypes.split(",")]
    if args.reference_transcript_biotypes:
        tx_biotypes = [b.strip() for b in args.reference_transcript_biotypes.split(",")]
    if gene_biotypes or tx_biotypes:
        print(f"Filtering reference by biotype: genes={gene_biotypes}, transcripts={tx_biotypes}")
        ref_genes, _ref_mrna, ref_exons, ref_cds = filter_by_biotype(
            ref_genes,
            _ref_mrna,
            ref_exons,
            ref_cds,
            gene_biotypes=gene_biotypes,
            transcript_biotypes=tx_biotypes,
        )

    # Transcript selection
    if args.reference_transcript_selection != "all":
        print(f"Selecting transcripts: {args.reference_transcript_selection}")
        ref_genes, _ref_mrna, ref_exons, ref_cds = select_transcripts(
            ref_genes,
            _ref_mrna,
            ref_exons,
            ref_cds,
            mode=args.reference_transcript_selection,
        )

    post_filter_ref_genes = len(ref_genes)
    post_filter_ref_tx = _ref_mrna["transcript_id"].nunique() if not _ref_mrna.empty else 0
    if post_filter_ref_genes != pre_filter_ref_genes:
        print(
            f"  Reference after filtering: {post_filter_ref_genes} genes "
            f"({pre_filter_ref_genes - post_filter_ref_genes} removed), "
            f"{post_filter_ref_tx} transcripts "
            f"({pre_filter_ref_tx - post_filter_ref_tx} removed)"
        )

    # --- Generate filter audit ---
    generate_filter_audit(
        ref_genes,
        _ref_mrna,
        ref_exons,
        ref_cds,
        output_dir=args.output_dir,
        evaluation_mode=args.evaluation_mode,
        transcript_selection=args.reference_transcript_selection,
        gene_biotypes_filter=args.reference_gene_biotypes,
        transcript_biotypes_filter=args.reference_transcript_biotypes,
        pre_filter_genes=pre_filter_ref_genes,
        pre_filter_transcripts=pre_filter_ref_tx,
    )

    # --- Classify ---
    # Apply subsetting before classification
    subset_regions = None
    if getattr(args, "sample_loci", None) is not None:
        # Build loci for sampling
        sample_source = getattr(args, "sample_from", "union")
        if sample_source == "reference":
            loci_df = ref_genes[["Chromosome", "Start", "End"]].copy()
        elif sample_source == "consensus":
            loci_df = cons_genes[["Chromosome", "Start", "End"]].copy()
        else:  # union
            loci_df = pd.concat(
                [
                    ref_genes[["Chromosome", "Start", "End"]],
                    cons_genes[["Chromosome", "Start", "End"]],
                ],
                ignore_index=True,
            )
        subset_regions = resolve_subset_regions(args, loci_df=loci_df)
    else:
        subset_regions = resolve_subset_regions(args)

    if subset_regions:
        print(f"  Subsetting to {len(subset_regions)} region(s)...")
        ref_genes = subset_df_by_regions(ref_genes, subset_regions)
        ref_exons = subset_df_by_regions(ref_exons, subset_regions)
        ref_cds = subset_df_by_regions(ref_cds, subset_regions)
        cons_genes = subset_df_by_regions(cons_genes, subset_regions)
        cons_exons = subset_df_by_regions(cons_exons, subset_regions)
        cons_cds = subset_df_by_regions(cons_cds, subset_regions)
        print(
            f"  After subsetting: {ref_genes.shape[0]} ref genes, "
            f"{cons_genes.shape[0]} cons genes"
        )
        # Write manifest
        manifest_path = os.path.join(args.output_dir, "subset_regions.tsv")
        write_subset_manifest(subset_regions, getattr(args, "seed", 1), manifest_path)

    print("Classifying loci...")
    ref_results, cons_results = classify_locus_pairs(
        ref_genes, ref_exons, ref_cds, cons_genes, cons_exons, cons_cds, cons_mrna
    )

    # Output transcript labels file
    print("Writing consensus transcript labels...")
    labels_path = os.path.join(args.output_dir, "consensus_transcript_labels.tsv")
    labeled_tx = []
    with open(labels_path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(
            [
                "transcript_id",
                "gene_id",
                "classification",
                "best_ref_gene_id",
                "best_ref_transcript_id",
                "best_overlap",
                "best_cds_overlap",
            ]
        )
        for c in cons_results:
            # The classification stored for the consensus gene is 'Matched', 'Strand_Mismatch', 'Novel' etc.
            # Convert classification classes to the 3 base labels
            cls = c["classification"]
            if cls in ("Exact_Match", "Partial_Match", "Structural_Mismatch", "Matched"):
                label = "Matched"
            elif cls == "Strand_Mismatch":
                label = "Strand_Mismatch"
            else:
                label = "Novel"

            for tid in c.get("tids", []):
                # Write a row per transcript
                w.writerow(
                    [
                        tid,
                        c["gene_id"],
                        label,
                        c["matched_id"],
                        c["best_ref_transcript_id"],
                        round(c["exon_overlap"], 4),
                        round(c["cds_overlap"], 4),
                    ]
                )
                labeled_tx.append({"transcript_id": tid, "comparison_label": label})

    # Evidence attribution labelling
    if args.evidence_attribution and os.path.exists(args.evidence_attribution):
        print("Labelling evidence attribution...")
        ea_df = pd.read_csv(args.evidence_attribution, sep="\t")
        labels_df = pd.DataFrame(labeled_tx)

        merged_df = ea_df.merge(labels_df, on="transcript_id", how="left")

        missing_labels = merged_df["comparison_label"].isna().sum()
        if missing_labels > 0:
            print(
                f"  Warning: {missing_labels} transcripts in evidence attribution lack comparison labels (e.g. subsetting)."
            )

        out_ea = os.path.join(args.output_dir, "evidence_attribution_labeled.tsv")
        merged_df.to_csv(out_ea, sep="\t", index=False)
        print(f"  Labeled evidence attribution: {out_ea}")

    # --- Summary ---
    print("Writing summary...")
    filter_info = {
        "evaluation_mode": args.evaluation_mode,
        "transcript_selection": args.reference_transcript_selection,
        "pre_filter_genes": pre_filter_ref_genes,
        "pre_filter_transcripts": pre_filter_ref_tx,
        "post_filter_genes": post_filter_ref_genes,
        "post_filter_transcripts": post_filter_ref_tx,
    }
    if args.reference_gene_biotypes:
        filter_info["gene_biotypes"] = args.reference_gene_biotypes
    if args.reference_transcript_biotypes:
        filter_info["transcript_biotypes"] = args.reference_transcript_biotypes

    # Sampling note
    sample_loci = getattr(args, "sample_loci", None)
    sample_seed = getattr(args, "seed", None)
    sample_source = getattr(args, "sample_from", "union")
    if sample_loci:
        filter_info["sampling"] = {
            "method": f"sampled {sample_loci} loci from {sample_source} (post-filter ref + query)",
            "n_loci": sample_loci,
            "seed": sample_seed,
            "source": sample_source,
            "note": "Loci are sampled AFTER biotype/transcript filtering. "
            "Changing evaluation mode or transcript selection resamples different loci. "
            "Percentages across modes are not directly comparable due to different sample composition.",
        }

    summary = write_summary(
        ref_results, cons_results, args.output_dir, filter_info=filter_info, ref_exons_df=ref_exons
    )
    # Include subset info in summary
    if subset_regions:
        summary["subset_regions"] = [str(r) for r in subset_regions]
        summary["subset_seed"] = sample_seed
        # Re-write the JSON with subset info
        json_path = os.path.join(args.output_dir, "comparison_summary.json")
        with open(json_path, "w") as fh:
            json.dump(summary, fh, indent=2)

    # Print key metrics
    print("\n" + "=" * 60)
    print("COMPARISON SUMMARY")
    print("=" * 60)
    print(f"Evaluation mode:          {args.evaluation_mode}")
    print(f"Transcript selection:     {args.reference_transcript_selection}")
    if pre_filter_ref_genes != post_filter_ref_genes:
        print(
            f"Reference (pre-filter):   {pre_filter_ref_genes} genes, {pre_filter_ref_tx} transcripts"
        )
    print(
        f"Reference (post-filter):  {post_filter_ref_genes} genes, {post_filter_ref_tx} transcripts"
    )
    if sample_loci:
        print(
            f"Sampling:                 {sample_loci} loci from {sample_source} (seed={sample_seed})"
        )
        print("  NOTE: loci sampled AFTER filtering; changing mode resamples different loci")
    print(f"Reference genes (sample): {summary['total_reference_genes']}")
    print(f"Query genes (sample):     {summary['total_consensus_genes']}")
    print()

    total_ref = summary["total_reference_genes"]
    sens = summary["sensitivity"]
    sens_cds = summary["sensitivity_cds"]
    ic = summary.get("intron_chain", {})

    print("--- 1. Locus Recovery ---")
    print(
        f"  Locus detection rate:       {sens.get('locus_detection_rate', 0):.1%}  ({sens.get('locus_detected_count', 0)}/{total_ref})"
    )
    print(
        f"  Missed rate:                {sens.get('missed_rate', 0):.1%}  ({sens.get('missed_count', 0)}/{total_ref})"
    )
    print(f"  Strand mismatch:            {sens.get('strand_mismatch_count', 0)}/{total_ref}")
    print()

    print("--- 2. CDS Accuracy ---")
    print(
        f"  CDS exact match:            {sens_cds.get('cds_exact_match_count', 0)}"
        + (
            f"  ({sens_cds.get('cds_exact_match_rate', 0):.1%})"
            if "cds_exact_match_rate" in sens_cds
            else ""
        )
    )
    print(
        f"  CDS any match:              {sens_cds.get('cds_any_match_count', 0)}"
        + (
            f"  ({sens_cds.get('cds_any_match_rate', 0):.1%})"
            if "cds_any_match_rate" in sens_cds
            else ""
        )
    )
    print(f"  CDS exact, UTR differs:     {sens_cds.get('cds_exact_but_exon_differs', 0)}")
    print()

    print("--- 3. Splice/Intron-Chain Recovery ---")
    print(
        f"  Intron chain recovered:     {ic.get('exon_intron_chain_recovered', 0)}/{ic.get('exon_intron_chain_matched', 0)}  ({ic.get('exon_intron_chain_rate', 0):.1%})"
    )
    cds_ic_rec = sens_cds.get("cds_intron_chain_recovered", 0)
    cds_ic_mat = sens_cds.get("cds_intron_chain_matched", 0)
    cds_ic_rate = cds_ic_rec / cds_ic_mat if cds_ic_mat > 0 else 0
    print(f"  CDS intron chain recovered: {cds_ic_rec}/{cds_ic_mat}  ({cds_ic_rate:.1%})")
    print()

    print("--- 4. Exact Transcript Structure ---")
    print(
        f"  Exact match rate:           {sens.get('exact_match_rate', 0):.1%}  ({sens.get('exact_match_count', 0)}/{total_ref})"
    )
    print(f"  Any match rate:             {sens.get('any_match_rate', 0):.1%}")
    print()

    print("Reference gene classification:")
    for cls in [
        "Exact_Match",
        "Structural_Mismatch",
        "Partial_Match",
        "Missed",
        "Strand_Mismatch",
    ]:
        cnt = summary["reference_classification"].get(cls, 0)
        pct = cnt / total_ref * 100 if total_ref > 0 else 0
        print(f"  {cls:25s}  {cnt:5d}  ({pct:.1f}%)")
    print()

    print("Query gene classification:")
    for cls in [
        "Exact_Match",
        "Structural_Mismatch",
        "Partial_Match",
        "Matched",
        "Novel",
        "Strand_Mismatch",
    ]:
        cnt = summary["consensus_classification"].get(cls, 0)
        if cnt == 0:
            continue
        pct = (
            cnt / summary["total_consensus_genes"] * 100
            if summary["total_consensus_genes"] > 0
            else 0
        )
        print(f"  {cls:25s}  {cnt:5d}  ({pct:.1f}%)")
    print(
        f"\nNovel query genes:          {summary['specificity'].get('novel_consensus_count', 0)}"
    )

    # --- Evidence tracks for plotting ---
    print("\nLoading evidence tracks for visualisation...")
    evidence_tracks = {}

    # Fixed track order for QC plots/legends, independent of CLI flag order:
    # Reference, Query (GMB consensus), ab initio backbone, short-read
    # transcriptomic, long-read consensus, protein evidence.
    evidence_tracks["Reference"] = (ref_exons, ref_cds)
    evidence_tracks["Query"] = (cons_exons, cons_cds)

    if args.helixer and args.tiberius:
        sys.exit(
            "ERROR: pass only one ab initio backbone track: --helixer or --tiberius, not both."
        )
    if args.helixer and os.path.exists(args.helixer):
        hx_exons, hx_cds, _ = load_gff(args.helixer, "Helixer")
        if mapping:
            hx_exons = remap_df_seqnames(hx_exons, mapping)
            if not hx_cds.empty:
                hx_cds = remap_df_seqnames(hx_cds, mapping)
        evidence_tracks["Helixer"] = (hx_exons, hx_cds)
    if args.tiberius and os.path.exists(args.tiberius):
        tb_exons, tb_cds, _ = load_gff(args.tiberius, "Tiberius")
        if mapping:
            tb_exons = remap_df_seqnames(tb_exons, mapping)
            if not tb_cds.empty:
                tb_cds = remap_df_seqnames(tb_cds, mapping)
        evidence_tracks["Tiberius"] = (tb_exons, tb_cds)
    if args.scallop and os.path.exists(args.scallop):
        sc_exons, _, _ = load_gff(args.scallop, "Scallop")
        if mapping:
            sc_exons = remap_df_seqnames(sc_exons, mapping)
        evidence_tracks["Scallop"] = (sc_exons, pd.DataFrame())
    if args.stringtie and os.path.exists(args.stringtie):
        st_exons, _, _ = load_gff(args.stringtie, "StringTie")
        if mapping:
            st_exons = remap_df_seqnames(st_exons, mapping)
        evidence_tracks["StringTie"] = (st_exons, pd.DataFrame())
    if args.minimap2 and os.path.exists(args.minimap2):
        mm_exons, _, _ = load_gff(args.minimap2, "Minimap2")
        if mapping:
            mm_exons = remap_df_seqnames(mm_exons, mapping)
        evidence_tracks["Minimap2"] = (mm_exons, pd.DataFrame())
    if args.orthodb and os.path.exists(args.orthodb):
        od_exons, _, _ = load_gff(args.orthodb, "OrthoDB")
        if mapping:
            od_exons = remap_df_seqnames(od_exons, mapping)
        evidence_tracks["OrthoDB"] = (od_exons, pd.DataFrame())
    if args.uniprot and os.path.exists(args.uniprot):
        up_exons, _, _ = load_gff(args.uniprot, "UniProt")
        if mapping:
            up_exons = remap_df_seqnames(up_exons, mapping)
        evidence_tracks["UniProt"] = (up_exons, pd.DataFrame())
    if args.genblast and os.path.exists(args.genblast):
        gb_exons, _, _ = load_gff(args.genblast, "GenBlast")
        if mapping:
            gb_exons = remap_df_seqnames(gb_exons, mapping)
        evidence_tracks["GenBlast"] = (gb_exons, pd.DataFrame())

    # --- Generate plots ---
    print("Generating comparison plots...")
    generate_comparison_plots(
        ref_results,
        cons_results,
        evidence_tracks,
        args.output_dir,
        plots_per_category=args.plots_per_category,
    )

    # --- Optionally extract ref proteins for BUSCO ---
    if args.extract_ref_proteins and fasta_path:
        print("Extracting reference proteins for BUSCO comparison...")
        genome = load_genome(fasta_path)
        ref_prot_path = os.path.join(args.output_dir, "reference_proteins.fa")
        extract_genbank_proteins(ref_exons, ref_cds, genome, ref_prot_path)

    print("\nDone.")


if __name__ == "__main__":
    main()
