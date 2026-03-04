#!/usr/bin/env python3
"""
Compare Annotations
===================
Locus-level comparison of consensus annotation against a reference
(GenBank community annotation).  Produces:
  - comparison_summary.json / .tsv
  - comparison_details.tsv (one row per gene with classification)
  - Representative locus plots per category (in qc/ subdirectory)

Usage:
    python compare_annotations.py \\
        --consensus consensus_genes.gff3 \\
        --reference GCA_002759435.3_Cand_auris_B8441_V3_genomic.gff.gz \\
        --assembly-report assembly_report.txt \\
        --genome candida_auris_softmasked_toplevel.fa \\
        --output-dir validation/ \\
        [--evidence-gtfs scallop.gtf stringtie.gtf orthodb.gtf uniprot.gtf] \\
        [--helixer helixer_remapped.gff3] \\
        [--plots-per-category 3]
"""

import argparse
import csv
import gzip
import json
import os
import sys
from collections import Counter, defaultdict

import matplotlib.pyplot as plt
import matplotlib.patches as patches
import pandas as pd
import pyranges as pr

from annotate_cds_utrs import (
    load_genome,
    build_spliced_seq,
    find_best_orf,
    translate,
    map_cds_to_genomic,
    derive_utrs,
    reverse_complement,
)

# ---------------------------------------------------------------------------
# Seqname remapping
# ---------------------------------------------------------------------------

def load_assembly_mapping(report_path):
    """Parse NCBI assembly report → {GenBank_accession: chromosome_name}."""
    mapping = {}
    with open(report_path) as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if len(parts) >= 5 and parts[2] != 'na':
                mapping[parts[4]] = parts[2]
    return mapping


def remap_df_seqnames(df, mapping):
    """Remap Chromosome column using assembly report mapping."""
    df = df.copy()
    df['Chromosome'] = df['Chromosome'].map(lambda x: mapping.get(x, x))
    return df


# ---------------------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------------------

def load_gff(path, source_label='Reference'):
    """Load GFF3 (optionally gzipped), return exon and CDS DataFrames."""
    print(f"  Loading {source_label} from {path}...")

    if path.endswith('.gz'):
        # Decompress to temp, load
        import tempfile
        tmp = tempfile.NamedTemporaryFile(suffix='.gff3', delete=False)
        with gzip.open(path, 'rt') as gz_in:
            tmp.write(gz_in.read().encode())
        tmp.close()
        try:
            gr = pr.read_gff3(tmp.name)
        finally:
            os.unlink(tmp.name)
    elif path.endswith('.gff3') or path.endswith('.gff'):
        gr = pr.read_gff3(path)
    elif path.endswith('.gtf'):
        gr = pr.read_gtf(path)
    else:
        gr = pr.read_gff3(path)

    df = gr.df
    df['Source_label'] = source_label

    # Normalise transcript_id
    if 'transcript_id' not in df.columns:
        if 'Parent' in df.columns:
            df['transcript_id'] = df['Parent']
        elif 'ID' in df.columns:
            df['transcript_id'] = df['ID']
        elif 'gene_id' in df.columns:
            df['transcript_id'] = df['gene_id']
        else:
            df['transcript_id'] = 'unknown'

    # For GenBank GFF: gene_id from ID attribute on gene rows
    if 'gene_id' not in df.columns:
        df['gene_id'] = df.get('ID', df['transcript_id'])

    exons = df[df['Feature'] == 'exon'].copy()
    cds = df[df['Feature'] == 'CDS'].copy()
    genes = df[df['Feature'] == 'gene'].copy()

    return exons, cds, genes


def load_consensus_genes(gff3_path):
    """Load consensus GFF3 and extract gene-level + exon-level data."""
    df = pr.read_gff3(gff3_path).df

    if 'transcript_id' not in df.columns:
        if 'Parent' in df.columns:
            df['transcript_id'] = df['Parent']
        elif 'ID' in df.columns:
            df['transcript_id'] = df['ID']

    genes = df[df['Feature'] == 'gene'].copy()
    exons = df[df['Feature'] == 'exon'].copy()
    cds = df[df['Feature'] == 'CDS'].copy()
    mrna = df[df['Feature'] == 'mRNA'].copy()

    return genes, exons, cds, mrna


# ---------------------------------------------------------------------------
# Locus-level classification
# ---------------------------------------------------------------------------

def _gene_span(gene_row):
    """Extract (chrom, start, end, strand) from a gene row."""
    return (gene_row['Chromosome'], gene_row['Start'],
            gene_row['End'], gene_row['Strand'])


def _exon_intervals(exon_df, tid):
    """Get sorted exon intervals for a transcript."""
    rows = exon_df[exon_df['transcript_id'] == tid]
    if rows.empty:
        return []
    return sorted(zip(rows['Start'].values, rows['End'].values))


def _intron_chain(exons):
    """Return intron chain signature from sorted exon list."""
    if len(exons) < 2:
        return 'single-exon'
    return ','.join(f'{exons[i][1]}-{exons[i+1][0]}'
                    for i in range(len(exons) - 1))


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


def classify_locus_pairs(ref_genes, ref_exons, cons_genes, cons_exons,
                         cons_mrna):
    """Classify each reference gene and each consensus gene.

    Returns:
        ref_results: list of dicts (one per ref gene)
        cons_results: list of dicts (one per consensus gene)
    """
    print("  Classifying loci...")

    # Build ref gene entries
    ref_entries = []
    for _, g in ref_genes.iterrows():
        gid = g.get('ID', g.get('gene_id', g.get('transcript_id', '')))
        ref_entries.append({
            'gene_id': gid,
            'chrom': g['Chromosome'],
            'start': g['Start'],
            'end': g['End'],
            'strand': g['Strand'],
        })

    # Build consensus gene entries
    cons_entries = []
    for _, g in cons_genes.iterrows():
        gid = g.get('ID', g.get('gene_id', g.get('transcript_id', '')))
        # Find the mRNA rows for this gene to get evidence info
        evidence = ''
        if cons_mrna is not None and not cons_mrna.empty:
            mrna_rows = cons_mrna[
                cons_mrna['transcript_id'].str.startswith(gid.split('.')[0])
                if isinstance(gid, str) else
                cons_mrna['transcript_id'] == gid
            ]
            # Try to get evidence from attributes
            if not mrna_rows.empty:
                for attr_col in ['Evidence', 'evidence']:
                    if attr_col in mrna_rows.columns:
                        evidence = mrna_rows.iloc[0].get(attr_col, '')
                        break

        cons_entries.append({
            'gene_id': gid,
            'chrom': g['Chromosome'],
            'start': g['Start'],
            'end': g['End'],
            'strand': g['Strand'],
            'evidence': evidence,
        })

    # Build PyRanges for overlap detection
    ref_pr = pr.PyRanges(pd.DataFrame({
        'Chromosome': [e['chrom'] for e in ref_entries],
        'Start': [e['start'] for e in ref_entries],
        'End': [e['end'] for e in ref_entries],
        'Strand': [e['strand'] for e in ref_entries],
        'gene_id': [e['gene_id'] for e in ref_entries],
    }))

    cons_pr = pr.PyRanges(pd.DataFrame({
        'Chromosome': [e['chrom'] for e in cons_entries],
        'Start': [e['start'] for e in cons_entries],
        'End': [e['end'] for e in cons_entries],
        'Strand': [e['strand'] for e in cons_entries],
        'gene_id': [e['gene_id'] for e in cons_entries],
    }))

    # Find overlaps (unstranded first to catch strand mismatches)
    ovl_unstrand = ref_pr.join(cons_pr, strandedness=False, suffix='_cons')

    # Build lookup: ref_gid → list of (cons_gid, same_strand)
    ref_to_cons = defaultdict(list)
    cons_to_ref = defaultdict(list)

    if not ovl_unstrand.df.empty:
        odf = ovl_unstrand.df
        for _, row in odf.iterrows():
            rgid = row['gene_id']
            cgid = row['gene_id_cons']
            same_strand = (row['Strand'] == row['Strand_cons'])
            ref_to_cons[rgid].append({
                'cons_id': cgid,
                'same_strand': same_strand,
            })
            cons_to_ref[cgid].append({
                'ref_id': rgid,
                'same_strand': same_strand,
            })

    # Classify reference genes
    ref_results = []
    for entry in ref_entries:
        gid = entry['gene_id']
        matches = ref_to_cons.get(gid, [])

        if not matches:
            entry['classification'] = 'Missed'
            entry['matched_id'] = ''
            entry['exon_overlap'] = 0.0
            ref_results.append(entry)
            continue

        # Check best match
        same_strand_matches = [m for m in matches if m['same_strand']]
        diff_strand_matches = [m for m in matches if not m['same_strand']]

        if not same_strand_matches and diff_strand_matches:
            entry['classification'] = 'Strand_Mismatch'
            entry['matched_id'] = diff_strand_matches[0]['cons_id']
            entry['exon_overlap'] = 0.0
            ref_results.append(entry)
            continue

        # Get exon intervals for ref gene (first transcript)
        ref_tids = ref_exons[
            ref_exons['Chromosome'] == entry['chrom']
        ]['transcript_id'].unique()
        ref_gene_exon_tids = [
            t for t in ref_tids
            if any(gid in str(t) for _ in [1])
        ]
        # Simpler: get exons overlapping the gene span
        ref_gene_exons_df = ref_exons[
            (ref_exons['Chromosome'] == entry['chrom']) &
            (ref_exons['Start'] >= entry['start']) &
            (ref_exons['End'] <= entry['end'])
        ]
        ref_gene_exons = sorted(zip(
            ref_gene_exons_df['Start'].values,
            ref_gene_exons_df['End'].values))

        best_class = 'Partial_Match'
        best_match = same_strand_matches[0]['cons_id']
        best_ovl = 0.0

        for m in same_strand_matches:
            cgid = m['cons_id']
            # Get consensus exons for this gene
            cons_gene_entry = next(
                (c for c in cons_entries if c['gene_id'] == cgid), None)
            if cons_gene_entry is None:
                continue

            cons_gene_exons_df = cons_exons[
                (cons_exons['Chromosome'] == cons_gene_entry['chrom']) &
                (cons_exons['Start'] >= cons_gene_entry['start']) &
                (cons_exons['End'] <= cons_gene_entry['end'])
            ]
            cons_gene_exons = sorted(zip(
                cons_gene_exons_df['Start'].values,
                cons_gene_exons_df['End'].values))

            if ref_gene_exons and cons_gene_exons:
                ovl_a = _exon_overlap_fraction(ref_gene_exons, cons_gene_exons)
                ovl_b = _exon_overlap_fraction(cons_gene_exons, ref_gene_exons)
                recip = min(ovl_a, ovl_b)

                if recip > best_ovl:
                    best_ovl = recip
                    best_match = cgid

                    ref_chain = _intron_chain(ref_gene_exons)
                    cons_chain = _intron_chain(cons_gene_exons)

                    if recip >= 0.8:
                        if ref_chain == cons_chain:
                            best_class = 'Exact_Match'
                        else:
                            best_class = 'Structural_Mismatch'
                    else:
                        best_class = 'Partial_Match'

        entry['classification'] = best_class
        entry['matched_id'] = best_match
        entry['exon_overlap'] = best_ovl
        ref_results.append(entry)

    # Classify consensus genes
    cons_results = []
    ref_matched_ids = {r['matched_id'] for r in ref_results
                       if r['classification'] != 'Missed'}

    for entry in cons_entries:
        gid = entry['gene_id']
        matches = cons_to_ref.get(gid, [])

        if not matches:
            entry['classification'] = 'Novel'
            entry['matched_id'] = ''
            entry['exon_overlap'] = 0.0
            cons_results.append(entry)
            continue

        # Find the best matching ref
        same_strand = [m for m in matches if m['same_strand']]
        if same_strand:
            entry['classification'] = 'Matched'
            entry['matched_id'] = same_strand[0]['ref_id']
        else:
            entry['classification'] = 'Strand_Mismatch'
            entry['matched_id'] = matches[0]['ref_id']
        entry['exon_overlap'] = 0.0
        cons_results.append(entry)

    return ref_results, cons_results


# ---------------------------------------------------------------------------
# Summary & detailed output
# ---------------------------------------------------------------------------

def write_summary(ref_results, cons_results, output_dir):
    """Write summary and detailed comparison files."""
    os.makedirs(output_dir, exist_ok=True)

    # Counts
    ref_class_counts = Counter(r['classification'] for r in ref_results)
    cons_class_counts = Counter(r['classification'] for r in cons_results)

    # Per-chromosome breakdown
    chr_breakdown = defaultdict(lambda: Counter())
    for r in ref_results:
        chr_breakdown[r['chrom']][r['classification']] += 1

    summary = {
        'total_reference_genes': len(ref_results),
        'total_consensus_genes': len(cons_results),
        'reference_classification': dict(ref_class_counts),
        'consensus_classification': dict(cons_class_counts),
        'per_chromosome': {
            c: dict(counts) for c, counts in sorted(chr_breakdown.items())
        },
        'sensitivity': {
            'exact_match_count': ref_class_counts.get('Exact_Match', 0),
            'partial_match_count': ref_class_counts.get('Partial_Match', 0),
            'any_match_count': (
                ref_class_counts.get('Exact_Match', 0) +
                ref_class_counts.get('Partial_Match', 0) +
                ref_class_counts.get('Structural_Mismatch', 0)
            ),
            'missed_count': ref_class_counts.get('Missed', 0),
            'strand_mismatch_count': ref_class_counts.get(
                'Strand_Mismatch', 0),
        },
        'specificity': {
            'novel_consensus_count': cons_class_counts.get('Novel', 0),
            'matched_consensus_count': cons_class_counts.get('Matched', 0),
        },
    }

    # Compute rates
    total_ref = len(ref_results)
    if total_ref > 0:
        summary['sensitivity']['exact_match_rate'] = round(
            ref_class_counts.get('Exact_Match', 0) / total_ref, 4)
        summary['sensitivity']['any_match_rate'] = round(
            summary['sensitivity']['any_match_count'] / total_ref, 4)
        summary['sensitivity']['missed_rate'] = round(
            ref_class_counts.get('Missed', 0) / total_ref, 4)

    # JSON
    json_path = os.path.join(output_dir, 'comparison_summary.json')
    with open(json_path, 'w') as fh:
        json.dump(summary, fh, indent=2)
    print(f"  Summary: {json_path}")

    # TSV summary
    tsv_path = os.path.join(output_dir, 'comparison_summary.tsv')
    with open(tsv_path, 'w', newline='') as fh:
        w = csv.writer(fh, delimiter='\t')
        w.writerow(['metric', 'value'])
        w.writerow(['total_reference_genes', len(ref_results)])
        w.writerow(['total_consensus_genes', len(cons_results)])
        for cls, cnt in sorted(ref_class_counts.items()):
            w.writerow([f'ref_{cls}', cnt])
        for cls, cnt in sorted(cons_class_counts.items()):
            w.writerow([f'cons_{cls}', cnt])
        for k, v in summary['sensitivity'].items():
            w.writerow([f'sens_{k}', v])
        for k, v in summary['specificity'].items():
            w.writerow([f'spec_{k}', v])
    print(f"  Summary: {tsv_path}")

    # Detailed per-gene TSV
    detail_path = os.path.join(output_dir, 'comparison_details.tsv')
    with open(detail_path, 'w', newline='') as fh:
        w = csv.writer(fh, delimiter='\t')
        w.writerow(['source', 'gene_id', 'chrom', 'start', 'end', 'strand',
                     'classification', 'matched_id', 'exon_overlap'])
        for r in ref_results:
            w.writerow(['Reference', r['gene_id'], r['chrom'],
                         r['start'], r['end'], r['strand'],
                         r['classification'], r['matched_id'],
                         f"{r['exon_overlap']:.3f}"])
        for r in cons_results:
            w.writerow(['Consensus', r['gene_id'], r['chrom'],
                         r['start'], r['end'], r['strand'],
                         r['classification'], r['matched_id'],
                         f"{r.get('exon_overlap', 0):.3f}"])
    print(f"  Details: {detail_path}")

    return summary


# ---------------------------------------------------------------------------
# Locus plotting
# ---------------------------------------------------------------------------

TRACK_COLORS = {
    'GenBank': '#D9534F',
    'Consensus': '#337AB7',
    'Scallop': '#5CB85C',
    'StringTie': '#F0AD4E',
    'Helixer': '#9B59B6',
    'OrthoDB': '#999999',
    'UniProt': '#888888',
}

CDS_HEIGHT = 0.4
UTR_HEIGHT = 0.18


def _draw_gene_track(ax, exon_df, cds_df, y, color, label,
                     chrom, plot_start, plot_end, max_tx=5):
    """Draw exon/CDS blocks for one track at a locus."""
    if exon_df is None or exon_df.empty:
        return

    # Filter to region
    region = exon_df[
        (exon_df['Chromosome'] == chrom) &
        (exon_df['End'] > plot_start) &
        (exon_df['Start'] < plot_end)
    ]
    if region.empty:
        return

    # CDS for this region
    region_cds = pd.DataFrame()
    if cds_df is not None and not cds_df.empty:
        region_cds = cds_df[
            (cds_df['Chromosome'] == chrom) &
            (cds_df['End'] > plot_start) &
            (cds_df['Start'] < plot_end)
        ]

    # Build CDS lookup by transcript
    cds_by_tx = defaultdict(list)
    if not region_cds.empty:
        tid_col = 'transcript_id' if 'transcript_id' in region_cds.columns else 'Parent'
        if tid_col in region_cds.columns:
            for _, row in region_cds.iterrows():
                cds_by_tx[row[tid_col]].append((row['Start'], row['End']))

    tids = list(region['transcript_id'].unique())
    if len(tids) > max_tx:
        tids = tids[:max_tx]

    for tid in tids:
        tx_exons = region[region['transcript_id'] == tid].sort_values('Start')
        tx_cds = cds_by_tx.get(tid, [])

        for _, ex in tx_exons.iterrows():
            ex_s, ex_e = ex['Start'], ex['End']

            if tx_cds:
                # Thin exon background
                rect = patches.Rectangle(
                    (ex_s, y - UTR_HEIGHT / 2), ex_e - ex_s, UTR_HEIGHT,
                    linewidth=0.5, edgecolor=color, facecolor=color, alpha=0.3)
                ax.add_patch(rect)

                # CDS portions
                for cs, ce in tx_cds:
                    os_ = max(cs, ex_s)
                    oe = min(ce, ex_e)
                    if os_ < oe:
                        rect = patches.Rectangle(
                            (os_, y - CDS_HEIGHT / 2), oe - os_, CDS_HEIGHT,
                            linewidth=1, edgecolor=color, facecolor=color,
                            alpha=0.7)
                        ax.add_patch(rect)
            else:
                # No CDS → thick exon block
                rect = patches.Rectangle(
                    (ex_s, y - CDS_HEIGHT / 2), ex_e - ex_s, CDS_HEIGHT,
                    linewidth=1, edgecolor=color, facecolor=color, alpha=0.7)
                ax.add_patch(rect)

        # Intron line
        if len(tx_exons) > 1:
            tx_s = tx_exons['Start'].min()
            tx_e = tx_exons['End'].max()
            ax.plot([tx_s, tx_e], [y, y], color=color, linewidth=1, zorder=0)


def plot_comparison_locus(chrom, start, end, category, ref_gid,
                          cons_gid, evidence_tracks, output_path,
                          genome=None):
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
        color = TRACK_COLORS.get(name, '#666666')
        y = y_positions[name]
        _draw_gene_track(ax, exon_df, cds_df, y, color, name,
                         chrom, plot_start, plot_end)

    # Highlight the locus region
    ax.axvspan(start, end, alpha=0.05, color='blue')

    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close()


def generate_comparison_plots(ref_results, cons_results,
                              evidence_tracks, output_dir,
                              plots_per_category=3):
    """Generate representative locus plots per classification category."""
    qc_dir = os.path.join(output_dir, 'qc')
    os.makedirs(qc_dir, exist_ok=True)

    # Group ref results by classification
    by_class = defaultdict(list)
    for r in ref_results:
        by_class[r['classification']].append(r)

    # Also add novel consensus
    for r in cons_results:
        if r['classification'] == 'Novel':
            by_class['Novel'].append(r)

    plot_count = 0
    for category, entries in by_class.items():
        # Sort by chromosome + position for reproducibility
        entries.sort(key=lambda e: (e['chrom'], e['start']))
        selected = entries[:plots_per_category]

        for i, entry in enumerate(selected):
            filename = f"{category.lower()}_{i+1}_{entry['chrom']}_{entry['start']}.png"
            filepath = os.path.join(qc_dir, filename)

            plot_comparison_locus(
                entry['chrom'], entry['start'], entry['end'],
                category,
                entry['gene_id'],
                entry.get('matched_id', ''),
                evidence_tracks,
                filepath)
            plot_count += 1

    print(f"  Generated {plot_count} comparison plots in {qc_dir}/")


# ---------------------------------------------------------------------------
# BUSCO helper
# ---------------------------------------------------------------------------

def extract_genbank_proteins(ref_exons, ref_cds, genome, output_path,
                             mapping=None):
    """Extract protein sequences from GenBank CDS coords.

    For BUSCO protein-mode comparison.
    """
    if ref_cds.empty:
        print("  Warning: No CDS features in reference — cannot extract proteins.")
        return 0

    if mapping:
        ref_cds = remap_df_seqnames(ref_cds, mapping)

    # Get unique transcript IDs with CDS
    tid_col = 'transcript_id' if 'transcript_id' in ref_cds.columns else 'Parent'

    count = 0
    with open(output_path, 'w') as fh:
        for tid, grp in ref_cds.groupby(tid_col):
            chrom = grp.iloc[0]['Chromosome']
            strand = grp.iloc[0]['Strand']

            if chrom not in genome:
                continue

            cds_intervals = sorted(
                zip(grp['Start'].values, grp['End'].values))

            # Build CDS nucleotide
            parts = []
            for cs, ce in (cds_intervals if strand == '+'
                           else reversed(cds_intervals)):
                parts.append(genome[chrom][cs:ce])
            cds_nuc = ''.join(parts)
            if strand == '-':
                cds_nuc = reverse_complement(cds_nuc)

            protein = translate(cds_nuc)
            if protein.endswith('*'):
                protein = protein[:-1]

            if len(protein) > 10:  # Skip very short
                fh.write(f">{tid}\n")
                for j in range(0, len(protein), 80):
                    fh.write(protein[j:j+80] + '\n')
                count += 1

    print(f"  Extracted {count} proteins to {output_path}")
    return count


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description='Compare consensus annotation against a reference '
                    '(GenBank community annotation)')
    parser.add_argument('--consensus', required=True,
                        help='Consensus GFF3')
    parser.add_argument('--reference', required=True,
                        help='Reference GFF3 (can be .gz)')
    parser.add_argument('--assembly-report', default=None,
                        help='NCBI assembly report for seqname remapping')
    parser.add_argument('--genome', default=None,
                        help='Genome FASTA (for protein extraction)')
    parser.add_argument('--output-dir', default='validation',
                        help='Output directory')
    parser.add_argument('--scallop', default=None, help='Scallop GTF')
    parser.add_argument('--stringtie', default=None, help='StringTie GTF')
    parser.add_argument('--helixer', default=None, help='Helixer GFF3')
    parser.add_argument('--orthodb', default=None, help='OrthoDB GTF')
    parser.add_argument('--uniprot', default=None, help='UniProt GTF')
    parser.add_argument('--plots-per-category', type=int, default=3,
                        help='Number of representative plots per category')
    parser.add_argument('--extract-ref-proteins', action='store_true',
                        default=False,
                        help='Extract protein FASTA from reference for BUSCO')
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    # --- Load assembly mapping ---
    mapping = None
    if args.assembly_report:
        mapping = load_assembly_mapping(args.assembly_report)
        print(f"Loaded assembly mapping: {len(mapping)} entries")

    # --- Load reference ---
    print("Loading reference annotation...")
    ref_exons, ref_cds, ref_genes = load_gff(
        args.reference, 'GenBank')
    if mapping:
        ref_exons = remap_df_seqnames(ref_exons, mapping)
        ref_cds = remap_df_seqnames(ref_cds, mapping)
        ref_genes = remap_df_seqnames(ref_genes, mapping)
    print(f"  Reference: {ref_genes.shape[0]} genes, "
          f"{ref_exons['transcript_id'].nunique()} transcripts")

    # --- Load consensus ---
    print("Loading consensus annotation...")
    cons_genes, cons_exons, cons_cds, cons_mrna = load_consensus_genes(
        args.consensus)
    print(f"  Consensus: {cons_genes.shape[0]} genes, "
          f"{cons_exons['transcript_id'].nunique()} transcripts")

    # --- Classify ---
    print("Classifying loci...")
    ref_results, cons_results = classify_locus_pairs(
        ref_genes, ref_exons, cons_genes, cons_exons, cons_mrna)

    # --- Summary ---
    print("Writing summary...")
    summary = write_summary(ref_results, cons_results, args.output_dir)

    # Print key metrics
    print("\n" + "=" * 60)
    print("COMPARISON SUMMARY")
    print("=" * 60)
    print(f"Reference genes:  {summary['total_reference_genes']}")
    print(f"Consensus genes:  {summary['total_consensus_genes']}")
    print()
    print("Reference gene classification:")
    for cls, cnt in sorted(summary['reference_classification'].items()):
        pct = cnt / summary['total_reference_genes'] * 100
        print(f"  {cls:25s}  {cnt:5d}  ({pct:.1f}%)")
    print()
    print("Consensus gene classification:")
    for cls, cnt in sorted(summary['consensus_classification'].items()):
        pct = cnt / summary['total_consensus_genes'] * 100
        print(f"  {cls:25s}  {cnt:5d}  ({pct:.1f}%)")
    print()

    sens = summary['sensitivity']
    print(f"Sensitivity (exact match):  {sens.get('exact_match_rate', 0):.1%}")
    print(f"Sensitivity (any match):    {sens.get('any_match_rate', 0):.1%}")
    print(f"Missed reference genes:     {sens.get('missed_count', 0)}")
    print(f"Novel consensus genes:      "
          f"{summary['specificity'].get('novel_consensus_count', 0)}")

    # --- Evidence tracks for plotting ---
    print("\nLoading evidence tracks for visualization...")
    evidence_tracks = {}

    # Always include ref and consensus
    evidence_tracks['GenBank'] = (ref_exons, ref_cds)
    evidence_tracks['Consensus'] = (cons_exons, cons_cds)

    # Optional evidence tracks
    if args.scallop and os.path.exists(args.scallop):
        sc_exons, _, _ = load_gff(args.scallop, 'Scallop')
        evidence_tracks['Scallop'] = (sc_exons, pd.DataFrame())
    if args.stringtie and os.path.exists(args.stringtie):
        st_exons, _, _ = load_gff(args.stringtie, 'StringTie')
        evidence_tracks['StringTie'] = (st_exons, pd.DataFrame())
    if args.helixer and os.path.exists(args.helixer):
        hx_exons, hx_cds, _ = load_gff(args.helixer, 'Helixer')
        if mapping:
            hx_exons = remap_df_seqnames(hx_exons, mapping)
            if not hx_cds.empty:
                hx_cds = remap_df_seqnames(hx_cds, mapping)
        evidence_tracks['Helixer'] = (hx_exons, hx_cds)
    if args.orthodb and os.path.exists(args.orthodb):
        od_exons, _, _ = load_gff(args.orthodb, 'OrthoDB')
        evidence_tracks['OrthoDB'] = (od_exons, pd.DataFrame())
    if args.uniprot and os.path.exists(args.uniprot):
        up_exons, _, _ = load_gff(args.uniprot, 'UniProt')
        evidence_tracks['UniProt'] = (up_exons, pd.DataFrame())

    # --- Generate plots ---
    print("Generating comparison plots...")
    generate_comparison_plots(
        ref_results, cons_results,
        evidence_tracks, args.output_dir,
        plots_per_category=args.plots_per_category)

    # --- Optionally extract ref proteins for BUSCO ---
    if args.extract_ref_proteins and args.genome:
        print("Extracting reference proteins for BUSCO comparison...")
        genome = load_genome(args.genome)
        ref_prot_path = os.path.join(args.output_dir, 'reference_proteins.fa')
        extract_genbank_proteins(ref_exons, ref_cds, genome, ref_prot_path)

    print("\nDone.")


if __name__ == '__main__':
    main()
