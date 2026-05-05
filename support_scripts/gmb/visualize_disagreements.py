import argparse
import os

import matplotlib.patches as patches
import matplotlib.pyplot as plt
import pyranges as pr

from annotate_cds_utrs import (
    _make_orf_label,
    build_spliced_seq,
    check_frame_continuity,
    check_splice_sites,
    derive_utrs,
    find_best_orf,
    get_start_stop_positions,
    load_genome,
    map_cds_to_genomic,
    reverse_complement,
    translate,
)
from subset_utils import (
    add_subset_args,
    build_mapping,
    remap_df_seqnames,
    resolve_subset_regions,
)

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

FILES = {
    "Reference": "Candida_Auris.gff3",
    "OldCode": "old_code_annotation.gtf",
    "Consensus": "new_consensus_output/consensus.gff3",
    "Scallop": "strict_scallop_annotation.gtf",
    "StringTie": "strict_stringtie_annotation.gtf",
    "Helixer": "helixer_remapped.gff3",
    "OrthoDB": "orthodb_annotation.gtf",
    "UniProt": "uniprot_annotation.gtf",
}

COLORS = {
    "Reference": "#D9534F",
    "OldCode": "#1ABC9C",
    "Consensus": "#337AB7",
    "Scallop": "#5CB85C",
    "StringTie": "#F0AD4E",
    "Helixer": "#9B59B6",
    "OrthoDB": "#999999",
    "UniProt": "#999999",
}

CDS_HEIGHT = 0.4
UTR_HEIGHT = 0.18
MARKER_SIZE = 6

# Tracks where ORF inference is skipped (protein alignments, not transcripts).
# This is the main QC speedup: these tracks have tens of thousands of entries
# and translating them all was the bottleneck.
_SKIP_ORF_INFERENCE = {"OrthoDB", "UniProt"}

# Default max transcripts per track per locus for QC plots.
_DEFAULT_MAX_TX_PER_TRACK = 5


# ---------------------------------------------------------------------------
# Annotation cache — avoids redundant ORF/CDS computation
# ---------------------------------------------------------------------------

_annotation_cache = {}


def _cache_key(exon_rows, strand, chrom):
    """Build a hashable key for annotation caching."""
    return (chrom, strand, tuple(sorted(exon_rows)))


def _cached_annotate_for_vis(
    exon_rows, strand, chrom, genome, cds_rows=None, infer_cds=False, min_orf_aa=100
):
    """Wrapper around _annotate_for_vis with memoisation."""
    key = _cache_key(exon_rows, strand, chrom)
    # Cache only covers ORF prediction; if CDS is provided, compute fresh
    if cds_rows:
        return _annotate_for_vis(
            exon_rows,
            strand,
            chrom,
            genome,
            cds_rows=cds_rows,
            infer_cds=infer_cds,
            min_orf_aa=min_orf_aa,
        )

    if key in _annotation_cache:
        return _annotation_cache[key]

    result = _annotate_for_vis(
        exon_rows,
        strand,
        chrom,
        genome,
        cds_rows=cds_rows,
        infer_cds=infer_cds,
        min_orf_aa=min_orf_aa,
    )
    _annotation_cache[key] = result
    return result


# ---------------------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------------------


def load_data():
    """Load exon and CDS data from all annotation files."""
    data = {}  # name → PyRanges of exons
    cds_data = {}  # name → PyRanges of CDS (where available)
    print("Loading datasets...")
    for name, path in FILES.items():
        if not os.path.exists(path):
            print(f"  Warning: File {path} not found.")
            continue

        print(f"  Loading {name}...")
        if path.endswith(".gtf"):
            gr = pr.read_gtf(path)
        else:
            gr = pr.read_gff3(path)

        df = gr.df

        # Load exons
        exons = df[df["Feature"] == "exon"].copy()
        if "transcript_id" not in exons.columns:
            if "Parent" in exons.columns:
                exons["transcript_id"] = exons["Parent"]
            elif "gene_id" in exons.columns:
                exons["transcript_id"] = exons["gene_id"]
            else:
                exons["transcript_id"] = "Unknown"
        data[name] = pr.PyRanges(exons)

        # Load CDS where available
        cds = df[df["Feature"] == "CDS"].copy()
        if not cds.empty:
            if "transcript_id" not in cds.columns:
                if "Parent" in cds.columns:
                    cds["transcript_id"] = cds["Parent"]
                elif "gene_id" in cds.columns:
                    cds["transcript_id"] = cds["gene_id"]
                else:
                    cds["transcript_id"] = "Unknown"
            cds_data[name] = pr.PyRanges(cds)

    return data, cds_data


# ---------------------------------------------------------------------------
# Per-transcript annotation for visualization
# ---------------------------------------------------------------------------


def _annotate_for_vis(
    exon_rows, strand, chrom, genome, cds_rows=None, infer_cds=False, min_orf_aa=100
):
    """Compute CDS/UTR/QC for one transcript in the plot.

    Parameters
    ----------
    exon_rows : list of (start, end)
        0-based half-open exon intervals.
    strand : str
    chrom : str
    genome : dict or None
    cds_rows : list of (start, end) or None
        Pre-existing CDS from the input file.
    infer_cds : bool
        If True and cds_rows is empty, run ORF prediction.
    min_orf_aa : int
        Minimum ORF length in codons.

    Returns
    -------
    dict with keys: cds, five_prime_utr, three_prime_utr,
                    start_pos, stop_pos, is_partial_5, is_partial_3,
                    splice_sites, frame_ok, orf_label
    """
    empty = {
        "cds": [],
        "five_prime_utr": [],
        "three_prime_utr": [],
        "start_pos": None,
        "stop_pos": None,
        "is_partial_5": False,
        "is_partial_3": False,
        "splice_sites": [],
        "frame_ok": True,
        "orf_label": "",
    }

    exons = sorted(exon_rows)

    # Path A: CDS provided in the input file
    if cds_rows:
        cds_intervals = sorted(cds_rows)
        five_utr, three_utr = derive_utrs(exons, cds_intervals, strand)
        start_pos, stop_pos = get_start_stop_positions(cds_intervals, strand)
        splice = []
        frame_ok = check_frame_continuity(cds_intervals, strand)
        if genome and chrom in genome:
            splice = check_splice_sites(exons, strand, genome[chrom])
            # Build protein for label: ascending genomic order then RC for
            # minus strand (RC(X+Y)==RC(Y)+RC(X) gives correct 5'→3' order).
            cds_nuc = "".join(genome[chrom][cs:ce] for cs, ce in cds_intervals)
            if strand == "-":
                cds_nuc = reverse_complement(cds_nuc)
            prot = translate(cds_nuc)
            if prot.endswith("*"):
                prot = prot[:-1]
            orf_label = _make_orf_label(prot, False, False)
        else:
            orf_label = ""
        return {
            "cds": cds_intervals,
            "five_prime_utr": five_utr,
            "three_prime_utr": three_utr,
            "start_pos": start_pos,
            "stop_pos": stop_pos,
            "is_partial_5": False,
            "is_partial_3": False,
            "splice_sites": splice,
            "frame_ok": frame_ok,
            "orf_label": orf_label,
        }

    # Path B: infer CDS from genome
    if not infer_cds or genome is None or chrom not in genome:
        return empty

    chrom_seq = genome[chrom]
    cdna = build_spliced_seq(exons, strand, chrom_seq)
    splice = check_splice_sites(exons, strand, chrom_seq)

    orf = find_best_orf(cdna, min_codons=min_orf_aa)
    if orf is None:
        return {**empty, "splice_sites": splice, "orf_label": "no ORF"}

    cds_start, cds_end, is_partial_5, is_partial_3 = orf
    cds_intervals = map_cds_to_genomic(cds_start, cds_end, exons, strand)
    five_utr, three_utr = derive_utrs(exons, cds_intervals, strand)

    cds_nuc = cdna[cds_start:cds_end]
    prot = translate(cds_nuc)
    if prot.endswith("*"):
        prot = prot[:-1]

    start_pos, stop_pos = get_start_stop_positions(cds_intervals, strand)
    if is_partial_5:
        start_pos = None
    if is_partial_3:
        stop_pos = None
    frame_ok = check_frame_continuity(cds_intervals, strand)
    orf_label = _make_orf_label(prot, is_partial_5, is_partial_3)

    return {
        "cds": cds_intervals,
        "five_prime_utr": five_utr,
        "three_prime_utr": three_utr,
        "start_pos": start_pos,
        "stop_pos": stop_pos,
        "is_partial_5": is_partial_5,
        "is_partial_3": is_partial_3,
        "splice_sites": splice,
        "frame_ok": frame_ok,
        "orf_label": orf_label,
    }


# ---------------------------------------------------------------------------
# Disagreement detection
# ---------------------------------------------------------------------------


def find_disagreements(data):
    """Find loci where Consensus disagrees with Reference."""
    print("Identifying disagreements...")
    ref = data.get("Reference")
    cons = data.get("Consensus")
    if ref is None or cons is None:
        return []

    any_overlap = ref.overlap(cons, strandedness=False)
    covered_ids = set(any_overlap.df["transcript_id"].unique())
    all_ref_ids = set(ref.df["transcript_id"].unique())
    missed_ids = list(all_ref_ids - covered_ids)

    disagreements = []

    # False negatives
    print(f"  Found {len(missed_ids)} False Negatives (Reference genes missed).")
    for tid in missed_ids[:5]:
        tx_exons = ref[ref.transcript_id == tid].df
        chrom = tx_exons.iloc[0]["Chromosome"]
        start = tx_exons["Start"].min()
        end = tx_exons["End"].max()
        disagreements.append((chrom, start, end, "False_Negative", tid))

    # Structural mismatches
    overlaps_df = ref.join(cons, strandedness=False).df
    structural_mismatches = []
    if not overlaps_df.empty:
        pairs = overlaps_df[["transcript_id", "transcript_id_b"]].drop_duplicates()
        for _, (rid, cid) in pairs.iterrows():
            if len(structural_mismatches) >= 5:
                break
            r_rows = ref[ref.transcript_id == rid].df
            c_rows = cons[cons.transcript_id == cid].df
            if len(r_rows) != len(c_rows):
                chrom = r_rows.iloc[0]["Chromosome"]
                start = min(r_rows["Start"].min(), c_rows["Start"].min())
                end = max(r_rows["End"].max(), c_rows["End"].max())
                structural_mismatches.append(
                    (chrom, start, end, "Structural_Mismatch", f"{rid}_vs_{cid}")
                )
                continue
            r_strand = r_rows.iloc[0]["Strand"]
            c_strand = c_rows.iloc[0]["Strand"]
            if r_strand != c_strand and c_strand != ".":
                chrom = r_rows.iloc[0]["Chromosome"]
                start = min(r_rows["Start"].min(), c_rows["Start"].min())
                end = max(r_rows["End"].max(), c_rows["End"].max())
                structural_mismatches.append(
                    (chrom, start, end, "Strand_Mismatch", f"{rid}_vs_{cid}")
                )

    print(f"  Found {len(structural_mismatches)} Structural Mismatches " "to visualize.")
    disagreements.extend(structural_mismatches)

    # -----------------------------------------------------------------------
    # Solved Overjoins
    # Cases where one OldCode transcript overlaps >=2 Reference genes
    # AND >=2 Consensus genes (implying OldCode fused them, but Consensus
    # successfully split them to mirror the reference).
    # -----------------------------------------------------------------------
    old_code = data.get("OldCode")
    if old_code is not None and ref is not None and cons is not None:
        print("  Looking for solved overjoins (OldCode fused genes resolved by Consensus)...")
        # 1. Overlap OldCode vs Reference
        joined_ref = old_code.join(ref, strandedness=False).df
        # 2. Overlap OldCode vs Consensus
        joined_cons = old_code.join(cons, strandedness=False).df

        solved = []
        if not joined_ref.empty and not joined_cons.empty:
            # Count distinct Reference transcript IDs per OldCode transcript ID
            ref_counts = joined_ref.groupby("transcript_id")["transcript_id_b"].nunique()
            # Count distinct Consensus transcript IDs per OldCode transcript ID
            cons_counts = joined_cons.groupby("transcript_id")["transcript_id_b"].nunique()

            # Find overjoined candidates
            overjoined_old_tids = ref_counts[ref_counts > 1].index

            for tid in overjoined_old_tids:
                has_cons = tid in cons_counts and cons_counts[tid] > 1
                if has_cons:
                    # Look up the bounds for this old code transcript to plot
                    oc_rows = old_code[old_code.transcript_id == tid].df
                    chrom = oc_rows.iloc[0]["Chromosome"]
                    start = oc_rows["Start"].min()
                    end = oc_rows["End"].max()

                    # We can pick a few to show
                    if len(solved) < 5:
                        solved.append((chrom, start, end, "Solved_Overjoin", tid))

        print(f"  Found {len(solved)} Solved Overjoins to visualize.")
        disagreements.extend(solved)

    return disagreements


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------


def plot_locus(
    data,
    chrom,
    start,
    end,
    label,
    out_filename,
    cds_data=None,
    genome=None,
    show_utrs=False,
    show_start_stop=False,
    infer_cds=False,
    min_orf_aa=100,
    show_splice_qc=False,
    show_orf_labels=False,
    show_frame_qc=False,
    max_tx_per_track=_DEFAULT_MAX_TX_PER_TRACK,
):
    """Plot annotation tracks for a genomic region.

    Default behaviour (no toggles) is identical to the original renderer:
    exons as thick coloured blocks, CDS/UTR differentiation only for tracks
    that already have CDS features in their input file.

    Feature toggles:
    - show_utrs: render UTR as lighter, thinner blocks for ALL tracks
    - show_start_stop: triangle/square at start/stop codon positions
    - infer_cds: run ORF prediction for tracks without CDS
    - show_splice_qc: red ✕ at non-canonical splice junctions
    - show_orf_labels: compact label at right edge of transcript
    - show_frame_qc: '!' warning at transcript label on frame breaks

    Performance improvements:
    - max_tx_per_track: limit transcripts plotted per track (sampling)
    - ORF inference is skipped for OrthoDB/UniProt tracks
    - Annotation results are cached by (chrom, strand, exon_chain)
    """
    pad = 1000
    plot_start = max(0, start - pad)
    plot_end = end + pad

    fig, ax = plt.subplots(figsize=(15, 8))

    track_names = list(FILES.keys())
    y_positions = {name: i for i, name in enumerate(reversed(track_names))}

    ax.set_xlim(plot_start, plot_end)
    ax.set_ylim(-1, len(track_names))
    ax.set_yticks(list(y_positions.values()))
    ax.set_yticklabels(list(y_positions.keys()))
    ax.set_title(f"Locus: {chrom}:{start}-{end} ({label})")
    ax.set_xlabel("Genomic Position")

    for name, pr_obj in data.items():
        color = COLORS.get(name, "#333333")
        y = y_positions[name]

        try:
            region = pr_obj[chrom, plot_start:plot_end]
        except Exception:
            continue
        if region.empty:
            continue

        # Pre-existing CDS from the input file for this track
        input_cds_intervals = []
        if cds_data and name in cds_data:
            try:
                cds_region = cds_data[name][chrom, plot_start:plot_end]
                if not cds_region.empty:
                    input_cds_intervals = list(
                        zip(
                            cds_region.df["Start"].values,
                            cds_region.df["End"].values,
                            cds_region.df["transcript_id"].values,
                        )
                    )
            except Exception:
                pass

        # Build CDS lookup by transcript
        cds_by_tx = {}
        for cs, ce, ctid in input_cds_intervals:
            cds_by_tx.setdefault(ctid, []).append((cs, ce))

        df = region.df

        # SPEEDUP: Sampling — limit transcripts per track
        all_tids = list(df["transcript_id"].unique())
        if len(all_tids) > max_tx_per_track:
            # Keep representative transcripts: sample by span diversity
            spans = df.groupby("transcript_id").agg(s=("Start", "min"), e=("End", "max"))
            spans["span"] = spans["e"] - spans["s"]
            spans = spans.sort_values("span", ascending=False)
            sampled_tids = set(spans.head(max_tx_per_track).index)
            df = df[df["transcript_id"].isin(sampled_tids)]

        tx_count = 0
        for tid, group in df.groupby("transcript_id"):
            group = group.sort_values("Start")
            strand = group.iloc[0]["Strand"]
            exon_list = list(zip(group["Start"].values, group["End"].values))

            # Determine CDS/UTR for this transcript
            tx_input_cds = cds_by_tx.get(tid, [])
            has_input_cds = len(tx_input_cds) > 0
            needs_annotation = (
                show_utrs
                or show_start_stop
                or infer_cds
                or show_splice_qc
                or show_orf_labels
                or show_frame_qc
            )

            ann = None
            if needs_annotation:
                # SPEEDUP: Skip ORF inference for protein alignment tracks
                track_infer = infer_cds and name not in _SKIP_ORF_INFERENCE
                # SPEEDUP: Use cached annotation
                ann = _cached_annotate_for_vis(
                    exon_list,
                    strand,
                    chrom,
                    genome,
                    cds_rows=tx_input_cds if has_input_cds else None,
                    infer_cds=track_infer,
                    min_orf_aa=min_orf_aa,
                )

            # Decide rendering mode
            use_cds_utr_rendering = False
            cds_for_rendering = tx_input_cds
            utr5_for_rendering = []
            utr3_for_rendering = []

            if ann and ann["cds"]:
                cds_for_rendering = ann["cds"]
                utr5_for_rendering = ann["five_prime_utr"]
                utr3_for_rendering = ann["three_prime_utr"]
                use_cds_utr_rendering = True
            elif has_input_cds:
                use_cds_utr_rendering = True

            # ~~~~~~~~~~~~ Draw exons ~~~~~~~~~~~~
            for _, exon in group.iterrows():
                ex_s, ex_e = exon["Start"], exon["End"]

                if use_cds_utr_rendering:
                    # CDS portions within this exon
                    cds_in_exon = [
                        (max(cs, ex_s), min(ce, ex_e))
                        for cs, ce in cds_for_rendering
                        if cs < ex_e and ce > ex_s
                    ]

                    if show_utrs and ann:
                        # UTR as lighter fill
                        utr_in_exon = []
                        for us, ue in utr5_for_rendering + utr3_for_rendering:
                            os_, oe = max(us, ex_s), min(ue, ex_e)
                            if os_ < oe:
                                utr_in_exon.append((os_, oe))
                        for us, ue in utr_in_exon:
                            rect = patches.Rectangle(
                                (us, y - UTR_HEIGHT / 2),
                                ue - us,
                                UTR_HEIGHT,
                                linewidth=0.5,
                                edgecolor=color,
                                facecolor=color,
                                alpha=0.3,
                            )
                            ax.add_patch(rect)
                    else:
                        # Thin underlying block (original behaviour)
                        rect = patches.Rectangle(
                            (ex_s, y - UTR_HEIGHT / 2),
                            ex_e - ex_s,
                            UTR_HEIGHT,
                            linewidth=0.5,
                            edgecolor=color,
                            facecolor=color,
                            alpha=0.4,
                        )
                        ax.add_patch(rect)

                    # CDS as thick dark blocks
                    for cs, ce in cds_in_exon:
                        rect = patches.Rectangle(
                            (cs, y - CDS_HEIGHT / 2),
                            ce - cs,
                            CDS_HEIGHT,
                            linewidth=1,
                            edgecolor=color,
                            facecolor=color,
                            alpha=0.7,
                        )
                        ax.add_patch(rect)
                else:
                    # No CDS — draw exon as thick block (default)
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

            # ~~~~~~~~~~~~ Intron line ~~~~~~~~~~~~
            if len(group) > 1:
                tx_min = group["Start"].min()
                tx_max = group["End"].max()
                ax.plot([tx_min, tx_max], [y, y], color=color, linewidth=1, zorder=0)

            # ~~~~~~~~~~~~ Optional glyphs ~~~~~~~~~~~~
            if ann is None:
                continue

            # Start/stop markers
            if show_start_stop:
                if ann["start_pos"] is not None:
                    sp = ann["start_pos"]
                    if plot_start <= sp <= plot_end:
                        ax.plot(
                            sp,
                            y + CDS_HEIGHT / 2 + 0.04,
                            marker="v",
                            markersize=MARKER_SIZE,
                            color="#2ECC40",
                            zorder=5,
                            clip_on=True,
                        )
                elif ann["is_partial_5"]:
                    # Hollow marker at transcript edge for partial
                    edge = exon_list[-1][1] if strand == "-" else exon_list[0][0]
                    if plot_start <= edge <= plot_end:
                        ax.plot(
                            edge,
                            y + CDS_HEIGHT / 2 + 0.04,
                            marker="v",
                            markersize=MARKER_SIZE,
                            color="#2ECC40",
                            markerfacecolor="none",
                            markeredgewidth=1.2,
                            zorder=5,
                            clip_on=True,
                        )

                if ann["stop_pos"] is not None:
                    sp = ann["stop_pos"]
                    if plot_start <= sp <= plot_end:
                        ax.plot(
                            sp,
                            y - CDS_HEIGHT / 2 - 0.04,
                            marker="s",
                            markersize=MARKER_SIZE - 1,
                            color="#FF4136",
                            zorder=5,
                            clip_on=True,
                        )
                elif ann["is_partial_3"]:
                    edge = exon_list[0][0] if strand == "-" else exon_list[-1][1]
                    if plot_start <= edge <= plot_end:
                        ax.plot(
                            edge,
                            y - CDS_HEIGHT / 2 - 0.04,
                            marker="s",
                            markersize=MARKER_SIZE - 1,
                            color="#FF4136",
                            markerfacecolor="none",
                            markeredgewidth=1.2,
                            zorder=5,
                            clip_on=True,
                        )

            # Splice QC
            if show_splice_qc and ann["splice_sites"]:
                for sj in ann["splice_sites"]:
                    if sj["class"] != "canonical":
                        jx = (sj["intron_start"] + sj["intron_end"]) / 2
                        if plot_start <= jx <= plot_end:
                            ax.plot(
                                jx,
                                y,
                                marker="x",
                                markersize=7,
                                markeredgewidth=2,
                                color="red",
                                zorder=6,
                            )

            # ORF label
            if show_orf_labels and ann["orf_label"]:
                tx_right = max(e for _, e in exon_list)
                if tx_right <= plot_end:
                    ax.text(
                        tx_right + (plot_end - plot_start) * 0.005,
                        y,
                        ann["orf_label"],
                        fontsize=6,
                        va="center",
                        color=color,
                        clip_on=True,
                    )

            # Frame QC
            if show_frame_qc and not ann["frame_ok"]:
                tx_left = min(s for s, _ in exon_list)
                ax.text(
                    tx_left - (plot_end - plot_start) * 0.01,
                    y,
                    "!",
                    fontsize=10,
                    fontweight="bold",
                    va="center",
                    ha="right",
                    color="red",
                )

            tx_count += 1

    plt.tight_layout()
    plt.savefig(out_filename, dpi=150)
    plt.close()
    print(f"  Saved plot to {out_filename}")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def main():
    parser = argparse.ArgumentParser(
        description="Visualize annotation disagreements with optional "
        "CDS/UTR, start/stop, and QC overlays"
    )

    # Feature toggles
    parser.add_argument(
        "--genome",
        default=None,
        help="Genome FASTA (required for --infer-cds " "and --show-splice-qc)",
    )
    parser.add_argument(
        "--show-utrs",
        action="store_true",
        default=False,
        help="Render UTR as lighter shade (default: off)",
    )
    parser.add_argument(
        "--show-start-stop",
        action="store_true",
        default=False,
        help="Mark start/stop codon positions (default: off)",
    )
    parser.add_argument(
        "--infer-cds",
        action="store_true",
        default=False,
        help="Infer ORF/CDS for exon-only transcripts " "(default: off)",
    )
    parser.add_argument(
        "--orf-policy",
        choices=["longest_atg_preferred", "longest_any"],
        default="longest_atg_preferred",
        help="ORF selection policy (default: " "longest_atg_preferred)",
    )
    parser.add_argument(
        "--min-orf-aa", type=int, default=100, help="Minimum ORF length in codons (default: 100)"
    )
    parser.add_argument(
        "--allow-partial-orf",
        action="store_true",
        default=True,
        help="Allow partial ORFs (default: true)",
    )
    parser.add_argument(
        "--show-splice-qc",
        action="store_true",
        default=False,
        help="Highlight non-canonical splice junctions " "(default: off)",
    )
    parser.add_argument(
        "--show-orf-labels",
        action="store_true",
        default=False,
        help="Show compact ORF status labels (default: off)",
    )
    parser.add_argument(
        "--show-frame-qc",
        action="store_true",
        default=False,
        help="Show frame-break warnings (default: off)",
    )
    parser.add_argument(
        "--max-tx-per-track",
        type=int,
        default=_DEFAULT_MAX_TX_PER_TRACK,
        help="Max transcripts per track per locus " f"(default: {_DEFAULT_MAX_TX_PER_TRACK})",
    )
    parser.add_argument(
        "--debug",
        action="store_true",
        default=False,
        help="Enable splice-qc + orf-labels + frame-qc",
    )
    parser.add_argument(
        "--sample-from-category",
        default=None,
        help="Sample disagreements from a specific category "
        "(e.g. Structural_Mismatch, False_Negative)",
    )
    add_subset_args(parser)

    args = parser.parse_args()

    # --debug shortcut
    if args.debug:
        args.show_splice_qc = True
        args.show_orf_labels = True
        args.show_frame_qc = True

    # Load genome if needed
    genome = None
    if args.genome:
        print(f"Loading genome from {args.genome}...")
        genome = load_genome(args.genome)
        print(f"  Loaded {len(genome)} sequences")

    data, cds_data = load_data()

    # --- Apply seqname mapping ---
    mapping = build_mapping(
        assembly_report=getattr(args, "assembly_report", None),
        seqname_map=getattr(args, "seqname_map", None),
    )
    if mapping:
        print("Applying seqname mapping...")
        for name in list(data.keys()):
            remapped = remap_df_seqnames(data[name].df, mapping, name)
            data[name] = pr.PyRanges(remapped)
        for name in list(cds_data.keys()):
            remapped = remap_df_seqnames(cds_data[name].df, mapping)
            cds_data[name] = pr.PyRanges(remapped)

    disagreements = find_disagreements(data)

    # --- Apply subsetting / sampling ---
    subset_regions = resolve_subset_regions(args)

    if subset_regions:
        # Filter disagreements to those overlapping subset regions
        disagreements = [
            (chrom, start, end, dtype, label)
            for chrom, start, end, dtype, label in disagreements
            if any(
                r.seqname == chrom and (r.is_whole_contig() or (r.start <= end and r.end >= start))
                for r in subset_regions
            )
        ]
        print(f"  After region subsetting: {len(disagreements)} loci")

    sample_n = getattr(args, "sample_loci", None)
    if sample_n is not None and disagreements:
        import random as _rng

        rng = _rng.Random(getattr(args, "seed", 1))

        category_filter = getattr(args, "sample_from_category", None)
        if category_filter:
            filtered = [d for d in disagreements if d[3] == category_filter]
            if filtered:
                disagreements = filtered
            else:
                print(
                    f"  Warning: no disagreements in category " f"'{category_filter}', using all"
                )

        n = min(sample_n, len(disagreements))
        disagreements = rng.sample(disagreements, n)
        print(f"  Sampled {n} disagreements (seed={args.seed})")

    bed_lines = []

    # Clear annotation cache for this run
    _annotation_cache.clear()

    print(f"Generating plots for {len(disagreements)} loci...")
    for i, (chrom, start, end, dtype, label) in enumerate(disagreements):
        filename = f"disagreement_{i+1}_{dtype}.png"
        plot_locus(
            data,
            chrom,
            start,
            end,
            f"{dtype} {label}",
            filename,
            cds_data=cds_data,
            genome=genome,
            show_utrs=args.show_utrs,
            show_start_stop=args.show_start_stop,
            infer_cds=args.infer_cds,
            min_orf_aa=args.min_orf_aa,
            show_splice_qc=args.show_splice_qc,
            show_orf_labels=args.show_orf_labels,
            show_frame_qc=args.show_frame_qc,
            max_tx_per_track=args.max_tx_per_track,
        )

        bed_lines.append(f"{chrom}\t{start}\t{end}\t{dtype}:{label}\n")

    with open("disagreements.bed", "w") as f:
        f.writelines(bed_lines)
    print("Saved disagreements.bed")


if __name__ == "__main__":
    main()
