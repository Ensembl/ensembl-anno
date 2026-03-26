#!/usr/bin/env python3
"""CDS & UTR annotation module.

Produces CDS + UTR for each transcript:
  - If input already contains CDS: derive UTR as exonic bases outside CDS.
  - If input lacks CDS: infer CDS by ORF prediction from the genome FASTA.

Behavioral reference: gtf_to_seq.pl  (Ensembl compute_translation logic).
"""

import argparse

import pandas as pd
import pyranges as pr

# ---------------------------------------------------------------------------
# Standard genetic code
# ---------------------------------------------------------------------------

CODON_TABLE = {
    "TTT": "F",
    "TTC": "F",
    "TTA": "L",
    "TTG": "L",
    "CTT": "L",
    "CTC": "L",
    "CTA": "L",
    "CTG": "L",
    "ATT": "I",
    "ATC": "I",
    "ATA": "I",
    "ATG": "M",
    "GTT": "V",
    "GTC": "V",
    "GTA": "V",
    "GTG": "V",
    "TCT": "S",
    "TCC": "S",
    "TCA": "S",
    "TCG": "S",
    "CCT": "P",
    "CCC": "P",
    "CCA": "P",
    "CCG": "P",
    "ACT": "T",
    "ACC": "T",
    "ACA": "T",
    "ACG": "T",
    "GCT": "A",
    "GCC": "A",
    "GCA": "A",
    "GCG": "A",
    "TAT": "Y",
    "TAC": "Y",
    "TAA": "*",
    "TAG": "*",
    "CAT": "H",
    "CAC": "H",
    "CAA": "Q",
    "CAG": "Q",
    "AAT": "N",
    "AAC": "N",
    "AAA": "K",
    "AAG": "K",
    "GAT": "D",
    "GAC": "D",
    "GAA": "E",
    "GAG": "E",
    "TGT": "C",
    "TGC": "C",
    "TGA": "*",
    "TGG": "W",
    "CGT": "R",
    "CGC": "R",
    "CGA": "R",
    "CGG": "R",
    "AGT": "S",
    "AGC": "S",
    "AGA": "R",
    "AGG": "R",
    "GGT": "G",
    "GGC": "G",
    "GGA": "G",
    "GGG": "G",
}

STOP_CODONS = {"TAA", "TAG", "TGA"}
START_CODON = "ATG"

COMPLEMENT = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")


# ---------------------------------------------------------------------------
# FASTA handling
# ---------------------------------------------------------------------------


def load_genome(fasta_path):
    """Load a FASTA file into a dict of {chrom: sequence (uppercase)}."""
    genome = {}
    current_name = None
    chunks = []
    with open(fasta_path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if current_name is not None:
                    genome[current_name] = "".join(chunks).upper()
                current_name = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line)
        if current_name is not None:
            genome[current_name] = "".join(chunks).upper()
    return genome


def reverse_complement(seq):
    """Return the reverse complement of a DNA sequence."""
    return seq.translate(COMPLEMENT)[::-1]


# ---------------------------------------------------------------------------
# Spliced transcript sequence
# ---------------------------------------------------------------------------


def build_spliced_seq(exons, strand, genome):
    """Build the spliced cDNA sequence for a transcript.

    Parameters
    ----------
    exons : list of (start, end)
        0-based half-open exon coordinates, already sorted by Start.
    strand : str
        '+' or '-'.
    genome : dict
        {chrom: seq} from load_genome.  Caller extracts the right chrom.

    Returns
    -------
    str : spliced cDNA (5′→3′ in transcript orientation).
    """
    # Concatenate exonic genomic sequence in genomic order
    parts = []
    for start, end in sorted(exons):
        parts.append(genome[start:end])
    cdna = "".join(parts)

    if strand == "-":
        cdna = reverse_complement(cdna)
    return cdna


# ---------------------------------------------------------------------------
# ORF finding
# ---------------------------------------------------------------------------


def find_best_orf(cdna_seq, min_codons=100):
    """Find the best ORF in a cDNA sequence.

    Policy (matching Ensembl compute_translation):
      1. Prefer ORFs starting with ATG.
      2. Among ATG-initiated ORFs, take the longest.
      3. If no ATG ORF ≥ min_codons, allow 5′-partial (starts at position 0
         in any frame) and 3′-partial (runs to end without stop).
      4. Final fallback: longest ORF at any length.

    Parameters
    ----------
    cdna_seq : str
        Full spliced cDNA (uppercase).
    min_codons : int
        Minimum ORF length in codons.

    Returns
    -------
    (cds_start, cds_end, is_partial_5, is_partial_3) or None
        0-based half-open coordinates on the cDNA.
        cds_start/cds_end are nucleotide positions (NOT codon positions).
    """
    seq_len = len(cdna_seq)
    if seq_len < 3:
        return None

    best_atg = None  # (length, start, end, partial5, partial3)
    best_non_atg = None

    for frame in range(3):
        # Find all ATG positions and stop positions in this frame
        atg_positions = []
        stop_positions = []
        i = frame
        while i + 3 <= seq_len:
            codon = cdna_seq[i : i + 3]
            if codon == START_CODON:
                atg_positions.append(i)
            if codon in STOP_CODONS:
                stop_positions.append(i)
            i += 3

        # Add sentinel at end for 3′-partial ORFs
        # (the stop position is at the very end of the sequence)
        _has_terminal_stop = len(stop_positions) > 0 and stop_positions[-1] + 3 >= seq_len - 2

        # For each ATG, find the next stop
        for atg_pos in atg_positions:
            orf_end = None
            partial_3 = False
            for stop_pos in stop_positions:
                if stop_pos > atg_pos:
                    orf_end = stop_pos + 3  # include stop codon
                    break
            if orf_end is None:
                # 3′-partial: runs to end
                orf_end = seq_len - ((seq_len - frame) % 3)
                # Ensure we don't go past end
                if orf_end <= atg_pos:
                    continue
                partial_3 = True
            length = orf_end - atg_pos
            n_codons = length // 3
            candidate = (n_codons, atg_pos, orf_end, False, partial_3)
            if best_atg is None or n_codons > best_atg[0]:
                best_atg = candidate

        # 5′-partial ORF: starts at frame start, first stop (or end)
        orf_start = frame
        partial_5 = True
        # Check if the frame itself starts with ATG
        if frame + 3 <= seq_len and cdna_seq[frame : frame + 3] == START_CODON:
            partial_5 = False  # Actually has ATG, covered above

        if partial_5:
            orf_end = None
            partial_3 = False
            for stop_pos in stop_positions:
                if stop_pos >= orf_start:
                    orf_end = stop_pos + 3
                    break
            if orf_end is None:
                orf_end = seq_len - ((seq_len - frame) % 3)
                partial_3 = True
            if orf_end > orf_start:
                length = orf_end - orf_start
                n_codons = length // 3
                candidate = (n_codons, orf_start, orf_end, True, partial_3)
                if best_non_atg is None or n_codons > best_non_atg[0]:
                    best_non_atg = candidate

    # Decision: prefer ATG ORFs
    if best_atg and best_atg[0] >= min_codons:
        _, start, end, p5, p3 = best_atg
        return (start, end, p5, p3)

    # Fallback to best non-ATG (partial) if long enough
    if best_non_atg and best_non_atg[0] >= min_codons:
        _, start, end, p5, p3 = best_non_atg
        return (start, end, p5, p3)

    # Final fallback: take whatever is longest
    candidates = [c for c in [best_atg, best_non_atg] if c is not None]
    if candidates:
        best = max(candidates, key=lambda c: c[0])
        if best[0] >= 1:  # At least 1 codon
            _, start, end, p5, p3 = best
            return (start, end, p5, p3)

    return None


def translate(cds_seq):
    """Translate a CDS nucleotide sequence to protein.

    Parameters
    ----------
    cds_seq : str
        CDS sequence (uppercase), length should be multiple of 3.

    Returns
    -------
    str : protein sequence (stop codon represented as '*').
    """
    protein = []
    for i in range(0, len(cds_seq) - 2, 3):
        codon = cds_seq[i : i + 3]
        aa = CODON_TABLE.get(codon, "X")
        protein.append(aa)
    return "".join(protein)


# ---------------------------------------------------------------------------
# Genomic coordinate mapping
# ---------------------------------------------------------------------------


def map_cds_to_genomic(cds_start, cds_end, exons, strand):
    """Map transcript-relative CDS coordinates to genomic CDS intervals.

    Parameters
    ----------
    cds_start, cds_end : int
        0-based half-open CDS coordinates on the spliced cDNA.
    exons : list of (genomic_start, genomic_end)
        0-based half-open, sorted by genomic position (ascending).
    strand : str
        '+' or '-'.

    Returns
    -------
    list of (genomic_start, genomic_end)
        CDS intervals in genomic coordinates, sorted ascending.
    """
    # Build the exon order as seen in the transcript (5′→3′)
    sorted_exons = sorted(exons)
    if strand == "-":
        # Transcript order is reverse of genomic order for minus strand
        tx_exons = list(reversed(sorted_exons))
    else:
        tx_exons = list(sorted_exons)

    # Walk through transcript-ordered exons, tracking cumulative position
    cds_intervals = []
    cum_pos = 0
    for g_start, g_end in tx_exons:
        exon_len = g_end - g_start
        exon_tx_start = cum_pos
        exon_tx_end = cum_pos + exon_len

        # Overlap between CDS and this exon in transcript coords
        overlap_start = max(cds_start, exon_tx_start)
        overlap_end = min(cds_end, exon_tx_end)

        if overlap_start < overlap_end:
            # Convert back to genomic coordinates
            if strand == "+":
                g_cds_start = g_start + (overlap_start - exon_tx_start)
                g_cds_end = g_start + (overlap_end - exon_tx_start)
            else:
                # Minus strand: transcript coords go from right to left
                g_cds_end = g_end - (overlap_start - exon_tx_start)
                g_cds_start = g_end - (overlap_end - exon_tx_start)
            cds_intervals.append((g_cds_start, g_cds_end))

        cum_pos += exon_len

    return sorted(cds_intervals)


def derive_utrs(exons, cds_intervals, strand):
    """Derive 5′ and 3′ UTR intervals from exons and CDS.

    Parameters
    ----------
    exons : list of (start, end)
        Genomic exon coordinates, 0-based half-open, sorted ascending.
    cds_intervals : list of (start, end)
        Genomic CDS coordinates, 0-based half-open, sorted ascending.
    strand : str
        '+' or '-'.

    Returns
    -------
    (five_prime_utrs, three_prime_utrs) : (list, list)
        Each is a list of (start, end) genomic intervals.
    """
    if not cds_intervals:
        return [], []

    sorted_exons = sorted(exons)
    sorted_cds = sorted(cds_intervals)

    cds_min = sorted_cds[0][0]
    cds_max = sorted_cds[-1][1]

    # UTR = exonic regions outside CDS
    # For + strand: 5′ UTR is upstream (left) of CDS, 3′ is downstream (right)
    # For - strand: 5′ UTR is downstream (right) of CDS, 3′ is upstream (left)
    upstream_utrs = []
    downstream_utrs = []

    for ex_start, ex_end in sorted_exons:
        # Portion of exon before CDS
        if ex_start < cds_min:
            upstream_utrs.append((ex_start, min(ex_end, cds_min)))
        # Portion of exon after CDS
        if ex_end > cds_max:
            downstream_utrs.append((max(ex_start, cds_max), ex_end))

        # Handle internal CDS gaps within a single exon (rare but possible)
        for cds_s, cds_e in sorted_cds:
            if cds_s > ex_start and cds_e < ex_end:
                # Check for gaps between CDS intervals within this exon
                pass  # Covered by the min/max logic above

    if strand == "+":
        return upstream_utrs, downstream_utrs
    else:
        return downstream_utrs, upstream_utrs


# ---------------------------------------------------------------------------
# QC helpers
# ---------------------------------------------------------------------------


def get_start_stop_positions(cds_intervals, strand):
    """Return genomic positions of start and stop codons.

    Parameters
    ----------
    cds_intervals : list of (start, end)
        Genomic CDS intervals, 0-based half-open, sorted ascending.
    strand : str
        '+' or '-'.

    Returns
    -------
    (start_pos, stop_pos) : (int or None, int or None)
        Genomic position (0-based) of the first base of the start codon
        and the first base of the stop codon.  None if partial.
    """
    if not cds_intervals:
        return None, None
    sorted_cds = sorted(cds_intervals)
    if strand == "+":
        start_pos = sorted_cds[0][0]  # leftmost CDS start
        stop_pos = sorted_cds[-1][1] - 3  # 3bp before CDS end
    else:
        start_pos = sorted_cds[-1][1] - 3  # rightmost 3bp (start on - strand)
        stop_pos = sorted_cds[0][0]  # leftmost CDS start (stop on -)
    return start_pos, stop_pos


def check_splice_sites(exons, strand, chrom_seq):
    """Classify splice junctions by donor/acceptor dinucleotides.

    Parameters
    ----------
    exons : list of (start, end)
        Genomic exon coords, 0-based half-open, sorted ascending.
    strand : str
        '+' or '-'.
    chrom_seq : str
        Chromosome sequence (uppercase).

    Returns
    -------
    list of dict
        One per intron, with keys:
        'intron_start', 'intron_end', 'donor', 'acceptor', 'class'.
        class is 'canonical', 'known_noncanonical', or 'noncanonical'.
    """
    sorted_exons = sorted(exons)
    results = []
    seq_len = len(chrom_seq)
    for i in range(len(sorted_exons) - 1):
        intron_start = sorted_exons[i][1]
        intron_end = sorted_exons[i + 1][0]
        if intron_end - intron_start < 4:
            continue

        # Raw genomic dinucleotides at intron boundaries
        donor_raw = chrom_seq[intron_start : min(intron_start + 2, seq_len)]
        acceptor_raw = chrom_seq[max(intron_end - 2, 0) : intron_end]

        if strand == "+":
            donor, acceptor = donor_raw, acceptor_raw
        else:
            # Reverse complement for minus strand
            donor = reverse_complement(acceptor_raw)
            acceptor = reverse_complement(donor_raw)

        pair = (donor, acceptor)
        if pair == ("GT", "AG"):
            cls = "canonical"
        elif pair in (("GC", "AG"), ("AT", "AC")):
            cls = "known_noncanonical"
        else:
            cls = "noncanonical"

        results.append(
            {
                "intron_start": intron_start,
                "intron_end": intron_end,
                "donor": donor,
                "acceptor": acceptor,
                "class": cls,
            }
        )
    return results


def check_frame_continuity(cds_intervals, strand):
    """Check whether total CDS length is a multiple of 3.

    Parameters
    ----------
    cds_intervals : list of (start, end)
        Genomic CDS intervals.
    strand : str
        '+' or '-'.

    Returns
    -------
    bool : True if frame is consistent (total CDS bp % 3 == 0).
    """
    if not cds_intervals:
        return True
    total = sum(e - s for s, e in cds_intervals)
    return total % 3 == 0


def _make_orf_label(protein, is_partial_5, is_partial_3):
    """Build a compact ORF status string like 'ORF: 312aa ATG+STOP'."""
    if protein is None:
        return "no ORF"
    aa_len = len(protein)
    parts = []
    if not is_partial_5:
        parts.append("ATG")
    else:
        parts.append("partial5")
    if not is_partial_3:
        parts.append("STOP")
    else:
        parts.append("partial3")
    return f'ORF:{aa_len}aa {"|".join(parts)}'


# ---------------------------------------------------------------------------
# Per-transcript annotation
# ---------------------------------------------------------------------------


def annotate_transcript(exons_df, chrom, strand, genome, cds_df=None, min_codons=100):
    """Annotate a single transcript with CDS, UTR, and QC info.

    Parameters
    ----------
    exons_df : DataFrame
        Exon rows with 'Start', 'End' columns (0-based half-open).
    chrom : str
        Chromosome name.
    strand : str
        '+' or '-'.
    genome : dict
        {chrom: sequence} from load_genome.
    cds_df : DataFrame or None
        If not None, existing CDS rows with 'Start', 'End'.
    min_codons : int
        Minimum ORF length for prediction.

    Returns
    -------
    dict with keys:
        'cds': list of (start, end)
        'five_prime_utr': list of (start, end)
        'three_prime_utr': list of (start, end)
        'cdna': str (spliced cDNA sequence)
        'protein': str or None
        'is_partial_5': bool
        'is_partial_3': bool
        'start_pos': int or None  (genomic position of start codon)
        'stop_pos': int or None   (genomic position of stop codon)
        'splice_sites': list of dict (intron classifications)
        'frame_ok': bool
        'orf_label': str (compact status string)
    """
    exons = sorted(zip(exons_df["Start"].values, exons_df["End"].values))

    if chrom not in genome:
        return _empty_result(exons)

    chrom_seq = genome[chrom]

    # Path A: CDS already provided
    if cds_df is not None and not cds_df.empty:
        cds_intervals = sorted(zip(cds_df["Start"].values, cds_df["End"].values))
        five_utr, three_utr = derive_utrs(exons, cds_intervals, strand)
        cdna = build_spliced_seq(exons, strand, chrom_seq)
        # Build CDS sequence for translation
        cds_seq_parts = []
        for cs, ce in sorted(cds_intervals) if strand == "+" else reversed(sorted(cds_intervals)):
            cds_seq_parts.append(chrom_seq[cs:ce])
        cds_nuc = "".join(cds_seq_parts)
        if strand == "-":
            cds_nuc = reverse_complement(cds_nuc)
        protein = translate(cds_nuc)
        # Remove trailing stop
        if protein.endswith("*"):
            protein = protein[:-1]
        start_pos, stop_pos = get_start_stop_positions(cds_intervals, strand)
        splice = check_splice_sites(exons, strand, chrom_seq)
        frame_ok = check_frame_continuity(cds_intervals, strand)
        return {
            "exons": exons,
            "cds": cds_intervals,
            "five_prime_utr": five_utr,
            "three_prime_utr": three_utr,
            "cdna": cdna,
            "protein": protein,
            "is_partial_5": False,
            "is_partial_3": False,
            "start_pos": start_pos,
            "stop_pos": stop_pos,
            "splice_sites": splice,
            "frame_ok": frame_ok,
            "orf_label": _make_orf_label(protein, False, False),
        }

    # Path B: predict ORF
    cdna = build_spliced_seq(exons, strand, chrom_seq)
    orf = find_best_orf(cdna, min_codons=min_codons)
    splice = check_splice_sites(exons, strand, chrom_seq)
    if orf is None:
        return {
            "exons": exons,
            "cds": [],
            "five_prime_utr": [],
            "three_prime_utr": [],
            "cdna": cdna,
            "protein": None,
            "is_partial_5": False,
            "is_partial_3": False,
            "start_pos": None,
            "stop_pos": None,
            "splice_sites": splice,
            "frame_ok": True,
            "orf_label": "no ORF",
        }

    cds_start, cds_end, is_partial_5, is_partial_3 = orf
    cds_intervals = map_cds_to_genomic(cds_start, cds_end, exons, strand)
    five_utr, three_utr = derive_utrs(exons, cds_intervals, strand)

    cds_nuc = cdna[cds_start:cds_end]
    protein = translate(cds_nuc)
    if protein.endswith("*"):
        protein = protein[:-1]

    start_pos, stop_pos = get_start_stop_positions(cds_intervals, strand)
    if is_partial_5:
        start_pos = None
    if is_partial_3:
        stop_pos = None
    frame_ok = check_frame_continuity(cds_intervals, strand)

    return {
        "exons": exons,
        "cds": cds_intervals,
        "five_prime_utr": five_utr,
        "three_prime_utr": three_utr,
        "cdna": cdna,
        "protein": protein,
        "is_partial_5": is_partial_5,
        "is_partial_3": is_partial_3,
        "start_pos": start_pos,
        "stop_pos": stop_pos,
        "splice_sites": splice,
        "frame_ok": frame_ok,
        "orf_label": _make_orf_label(protein, is_partial_5, is_partial_3),
    }


def _empty_result(exons):
    return {
        "exons": exons,
        "cds": [],
        "five_prime_utr": [],
        "three_prime_utr": [],
        "cdna": "",
        "protein": None,
        "is_partial_5": False,
        "is_partial_3": False,
        "start_pos": None,
        "stop_pos": None,
        "splice_sites": [],
        "frame_ok": True,
        "orf_label": "no ORF",
    }


# ---------------------------------------------------------------------------
# Batch annotation
# ---------------------------------------------------------------------------


def annotate_all_transcripts(exon_df, genome, cds_df=None, min_codons=100):
    """Annotate all transcripts in a DataFrame.

    Parameters
    ----------
    exon_df : DataFrame
        Must have: Chromosome, Start, End, Strand, transcript_id.
    genome : dict
        From load_genome.
    cds_df : DataFrame or None
        If provided, must have: Chromosome, Start, End, Strand, transcript_id
        (Parent column mapped to transcript_id).
    min_codons : int
        Minimum ORF length in codons.

    Returns
    -------
    dict : {transcript_id: annotation_dict}
    """
    results = {}

    # Build CDS lookup
    cds_by_tid = {}
    if cds_df is not None and not cds_df.empty:
        for tid, grp in cds_df.groupby("transcript_id"):
            cds_by_tid[tid] = grp

    for tid, grp in exon_df.groupby("transcript_id"):
        chrom = grp["Chromosome"].iloc[0]
        strand = grp["Strand"].iloc[0]
        existing_cds = cds_by_tid.get(tid)
        results[tid] = annotate_transcript(
            grp, chrom, strand, genome, cds_df=existing_cds, min_codons=min_codons
        )

    return results


# ---------------------------------------------------------------------------
# GFF3 output helpers
# ---------------------------------------------------------------------------


def annotations_to_gff_rows(annotations, transcript_id, parent_tid):
    """Convert annotation dict to GFF3-style row dicts.

    Returns list of dicts with keys: Chromosome, Feature, Start, End, Strand,
    transcript_id, Parent.
    """
    rows = []
    # We don't know the chrom/strand here — caller should supply or we skip
    return rows


def write_augmented_gff3(input_gff3, annotations, output_path):
    """Write an augmented GFF3 with CDS and UTR features added."""
    df = pr.read_gff3(input_gff3).df

    new_rows = []
    for tid, ann in annotations.items():
        if not ann["cds"]:
            continue

        # Find the chromosome and strand for this transcript
        tx_rows = df[
            (df.get("Parent", df.get("transcript_id", "")) == tid) | (df.get("ID", "") == tid)
        ]
        if tx_rows.empty:
            continue
        chrom = tx_rows["Chromosome"].iloc[0]
        strand = tx_rows["Strand"].iloc[0]

        for i, (s, e) in enumerate(ann["cds"]):
            new_rows.append(
                {
                    "Chromosome": chrom,
                    "Source": "GeneBuilder",
                    "Feature": "CDS",
                    "Start": s,
                    "End": e,
                    "Score": ".",
                    "Strand": strand,
                    "Frame": ".",
                    "ID": f"{tid}.cds{i+1}",
                    "Parent": tid,
                }
            )
        for i, (s, e) in enumerate(ann["five_prime_utr"]):
            new_rows.append(
                {
                    "Chromosome": chrom,
                    "Source": "GeneBuilder",
                    "Feature": "five_prime_UTR",
                    "Start": s,
                    "End": e,
                    "Score": ".",
                    "Strand": strand,
                    "Frame": ".",
                    "ID": f"{tid}.5utr{i+1}",
                    "Parent": tid,
                }
            )
        for i, (s, e) in enumerate(ann["three_prime_utr"]):
            new_rows.append(
                {
                    "Chromosome": chrom,
                    "Source": "GeneBuilder",
                    "Feature": "three_prime_UTR",
                    "Start": s,
                    "End": e,
                    "Score": ".",
                    "Strand": strand,
                    "Frame": ".",
                    "ID": f"{tid}.3utr{i+1}",
                    "Parent": tid,
                }
            )

    if new_rows:
        new_df = pd.DataFrame(new_rows)
        df = pd.concat([df, new_df], ignore_index=True)

    df = df.sort_values(["Chromosome", "Start"])
    with open(output_path, "w") as fh:
        fh.write("##gff-version 3\n")
        for _, r in df.iterrows():
            attr = f"ID={r.get('ID', '.')};Parent={r.get('Parent', '.')}"
            fh.write(
                f"{r['Chromosome']}\t{r.get('Source', 'GeneBuilder')}\t"
                f"{r['Feature']}\t{r['Start']}\t{r['End']}\t.\t"
                f"{r['Strand']}\t.\t{attr}\n"
            )

    print(f"Wrote augmented GFF3 to {output_path}")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def main():
    parser = argparse.ArgumentParser(description="Annotate transcripts with CDS and UTR features")
    parser.add_argument("--input", required=True, help="Input GFF3 or GTF with exon features")
    parser.add_argument("--genome", required=True, help="Genome FASTA file")
    parser.add_argument("--output", required=True, help="Output augmented GFF3")
    parser.add_argument(
        "--min-codons", type=int, default=100, help="Minimum ORF length in codons (default: 100)"
    )
    args = parser.parse_args()

    print("Loading genome...")
    genome = load_genome(args.genome)
    print(f"  Loaded {len(genome)} sequences")

    print("Loading annotations...")
    if args.input.endswith(".gtf"):
        gr = pr.read_gtf(args.input)
    else:
        gr = pr.read_gff3(args.input)
    df = gr.df

    # Separate exons and CDS
    exon_df = df[df["Feature"] == "exon"].copy()
    cds_df = df[df["Feature"] == "CDS"].copy()

    # Ensure transcript_id column
    for sub_df in [exon_df, cds_df]:
        if "transcript_id" not in sub_df.columns:
            if "Parent" in sub_df.columns:
                sub_df["transcript_id"] = sub_df["Parent"]

    if cds_df.empty:
        cds_df = None

    print(
        f'  {exon_df["transcript_id"].nunique()} transcripts, '
        f"{len(cds_df) if cds_df is not None else 0} existing CDS features"
    )

    print("Annotating transcripts...")
    annotations = annotate_all_transcripts(
        exon_df, genome, cds_df=cds_df, min_codons=args.min_codons
    )

    n_with_cds = sum(1 for a in annotations.values() if a["cds"])
    n_total = len(annotations)
    print(f"  {n_with_cds}/{n_total} transcripts have CDS")

    print("Writing output...")
    write_augmented_gff3(args.input, annotations, args.output)
    print("Done.")


if __name__ == "__main__":
    main()
