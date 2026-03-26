#!/usr/bin/env python3
"""Shared utilities for region subsetting, seqname mapping, and locus sampling.

Used by gene_model_builder.py, compare_annotations.py, and
visualize_disagreements.py to provide a consistent "fast test" capability.
"""

import csv
import os
import random
from dataclasses import dataclass
from typing import List, Optional

import pandas as pd

try:
    import pyranges as pr
except ImportError:
    pr = None


# ---------------------------------------------------------------------------
# Region dataclass
# ---------------------------------------------------------------------------


@dataclass
class Region:
    """A genomic region: whole-contig if start/end are None."""

    seqname: str
    start: Optional[int] = None
    end: Optional[int] = None

    def is_whole_contig(self):
        return self.start is None and self.end is None

    def __str__(self):
        if self.is_whole_contig():
            return self.seqname
        return f"{self.seqname}:{self.start}-{self.end}"


# ---------------------------------------------------------------------------
# Region parsing
# ---------------------------------------------------------------------------


def parse_region(s: str) -> Region:
    """Parse 'seqname:start-end' or 'seqname' into a Region.

    Examples
    --------
    >>> parse_region('chr1:100-200')
    Region(seqname='chr1', start=100, end=200)
    >>> parse_region('chr1')
    Region(seqname='chr1', start=None, end=None)
    """
    s = s.strip()
    if ":" in s:
        seqname, coords = s.split(":", 1)
        if "-" in coords:
            parts = coords.split("-", 1)
            return Region(seqname, int(parts[0]), int(parts[1]))
        else:
            # Single position — treat as 1-bp window
            pos = int(coords)
            return Region(seqname, pos, pos)
    return Region(s)


def load_regions_file(path: str) -> List[Region]:
    """Load regions from a file, one per line.  # comments and blank lines ignored."""
    regions = []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            regions.append(parse_region(line))
    return regions


# ---------------------------------------------------------------------------
# Seqname mapping (extracted from compare_annotations.py)
# ---------------------------------------------------------------------------


def load_assembly_mapping(report_path: str) -> dict:
    """Parse NCBI assembly report → {GenBank_accession: chromosome_name}."""
    mapping = {}
    with open(report_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) >= 5 and parts[2] != "na":
                mapping[parts[4]] = parts[2]
    return mapping


def load_seqname_map(path: str) -> dict:
    """Parse custom seqname mapping TSV/CSV.

    Expects 2 columns: from_seqname, to_seqname.
    Ignores comments (#) and handles optional headers.
    """
    mapping = {}
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 2:
                parts = line.split(",")

            if len(parts) >= 2:
                from_id = parts[0].strip()
                to_id = parts[1].strip()
                # Skip header rows
                if (
                    from_id.lower() in ("from_seqname", "from", "seqname", "old", "old_name")
                    and "to" in to_id.lower()
                ):
                    continue
                mapping[from_id] = to_id
    return mapping


def build_mapping(assembly_report=None, seqname_map=None) -> dict:
    """Load and merge seqname mappings.  seqname_map takes precedence."""
    mapping = {}
    if assembly_report and os.path.exists(assembly_report):
        report_map = load_assembly_mapping(assembly_report)
        mapping.update(report_map)
        print(f"  Loaded assembly mapping: {len(report_map)} entries")
    if seqname_map and os.path.exists(seqname_map):
        seq_map = load_seqname_map(seqname_map)
        mapping.update(seq_map)  # takes precedence
        print(
            f"  Loaded seqname mapping: {len(seq_map)} entries "
            "(overrides assembly-report on collisions)"
        )
    return mapping


def remap_df_seqnames(df, mapping, label=""):
    """Remap Chromosome column using a mapping dict.

    Prints safety information about uniquely mapped and unmapped seqnames.
    """
    if df is None or df.empty or not mapping:
        return df

    df = df.copy()

    unique_before = set(df["Chromosome"].unique())
    num_before = len(unique_before)

    df["Chromosome"] = df["Chromosome"].map(lambda x: mapping.get(x, x))

    unique_after = set(df["Chromosome"].unique())
    num_after = len(unique_after)

    if label:
        unmapped = sorted(x for x in unique_before if x not in mapping)
        print(f"  {label} seqnames remapped: {num_before} unique -> " f"{num_after} unique")
        if unmapped:
            print(f"  Warning: {label} has unmapped seqnames: " f"{', '.join(unmapped)}")
    return df


# ---------------------------------------------------------------------------
# Subsetting
# ---------------------------------------------------------------------------


def subset_df_by_regions(df, regions: List[Region]) -> pd.DataFrame:
    """Keep rows whose (Chromosome, Start..End) overlaps at least one region.

    Works on any DataFrame with Chromosome, Start, End columns.
    """
    if df is None or df.empty or not regions:
        return df

    masks = []
    for r in regions:
        if r.is_whole_contig():
            masks.append(df["Chromosome"] == r.seqname)
        else:
            masks.append(
                (df["Chromosome"] == r.seqname) & (df["End"] > r.start) & (df["Start"] < r.end)
            )

    combined = masks[0]
    for m in masks[1:]:
        combined = combined | m

    return df[combined].copy()


# ---------------------------------------------------------------------------
# Locus sampling
# ---------------------------------------------------------------------------


def _build_loci_from_exons(exon_df):
    """Quick locus building: cluster overlapping exons into loci.

    Returns DataFrame with columns: Chromosome, Start, End, locus_id.
    """
    if exon_df is None or exon_df.empty:
        return pd.DataFrame(columns=["Chromosome", "Start", "End", "locus_id"])

    if pr is None:
        raise ImportError("pyranges is required for locus sampling")

    # Build minimal PR for clustering
    minimal = exon_df[["Chromosome", "Start", "End"]].copy()
    if "Strand" in exon_df.columns:
        minimal["Strand"] = exon_df["Strand"]

    gr = pr.PyRanges(minimal)
    try:
        clustered = gr.cluster(slack=0, count=True)
    except TypeError:
        clustered = gr.cluster(count=True)

    cdf = clustered.df

    loci = (
        cdf.groupby("Cluster")
        .agg(
            Chromosome=("Chromosome", "first"),
            Start=("Start", "min"),
            End=("End", "max"),
        )
        .reset_index()
    )
    loci = loci.rename(columns={"Cluster": "locus_id"})
    return loci


def sample_loci(
    loci_df,
    n: int,
    strategy: str = "uniform_locus",
    seed: int = 1,
    window_bp: int = 0,
) -> List[Region]:
    """Reproducible locus sampling.

    Parameters
    ----------
    loci_df : DataFrame
        Must have Chromosome, Start, End columns (one row per locus).
    n : int
        Number of loci to sample.
    strategy : str
        'uniform_locus': sample from locus list uniformly.
        'uniform_genome': sample random genomic windows then include
            overlapping loci.
    seed : int
        Random seed for reproducibility.
    window_bp : int
        Expand each sampled locus by this many bp on each side.

    Returns
    -------
    list[Region]
    """
    rng = random.Random(seed)

    if loci_df is None or loci_df.empty:
        return []

    n = min(n, len(loci_df))

    if strategy == "uniform_locus":
        sampled = loci_df.sample(n=n, random_state=seed).copy()
    elif strategy == "uniform_genome":
        # Sample random genomic positions then find overlapping loci
        chrom_lengths = loci_df.groupby("Chromosome")["End"].max().to_dict()
        total_bp = sum(chrom_lengths.values())

        sampled_indices = set()
        attempts = 0
        max_attempts = n * 100

        while len(sampled_indices) < n and attempts < max_attempts:
            # Pick random position in genome
            pos = rng.randint(0, total_bp - 1)
            cumulative = 0
            target_chrom = None
            target_pos = 0
            for chrom, length in sorted(chrom_lengths.items()):
                if cumulative + length > pos:
                    target_chrom = chrom
                    target_pos = pos - cumulative
                    break
                cumulative += length

            if target_chrom is None:
                attempts += 1
                continue

            # Find loci overlapping this position
            overlapping = loci_df[
                (loci_df["Chromosome"] == target_chrom)
                & (loci_df["Start"] <= target_pos)
                & (loci_df["End"] >= target_pos)
            ]

            if not overlapping.empty:
                sampled_indices.add(overlapping.index[0])
            attempts += 1

        sampled = loci_df.loc[list(sampled_indices)].copy()
    else:
        raise ValueError(f"Unknown sampling strategy: {strategy}")

    # Build regions with optional window expansion
    regions = []
    for _, row in sampled.iterrows():
        start = max(0, int(row["Start"]) - window_bp)
        end = int(row["End"]) + window_bp
        regions.append(Region(row["Chromosome"], start, end))

    # Sort for deterministic order
    regions.sort(key=lambda r: (r.seqname, r.start or 0))
    return regions


# ---------------------------------------------------------------------------
# Shared CLI integration
# ---------------------------------------------------------------------------


def add_subset_args(parser):
    """Add consistent subsetting CLI flags to an argparse parser."""
    grp = parser.add_argument_group("Subsetting / Fast Test")

    # Region-based
    grp.add_argument(
        "--seqname",
        default=None,
        help="Single contig/chromosome to subset to " "(after seqname remapping)",
    )
    grp.add_argument("--region", default=None, help="Region to subset to (seqname:start-end)")
    grp.add_argument(
        "--regions-file", default=None, help="File with regions (one per line, # comments ok)"
    )

    # Locus sampling
    grp.add_argument(
        "--sample-loci", type=int, default=None, help="Sample N loci for quick testing"
    )
    grp.add_argument(
        "--sample-bp",
        type=int,
        default=0,
        help="Expand sampled loci by this many bp " "(default: 0)",
    )
    grp.add_argument(
        "--seed", type=int, default=1, help="Random seed for reproducible sampling " "(default: 1)"
    )
    grp.add_argument(
        "--sample-strategy",
        choices=["uniform_locus", "uniform_genome"],
        default="uniform_locus",
        help="Sampling strategy (default: uniform_locus)",
    )

    # Mapping (if not already present — don't add duplicates)
    existing = {a.dest for a in parser._actions}
    if "assembly_report" not in existing:
        grp.add_argument(
            "--assembly-report", default=None, help="NCBI assembly report for seqname remapping"
        )
    if "seqname_map" not in existing:
        grp.add_argument(
            "--seqname-map",
            default=None,
            help="Custom TSV/CSV seqname mapping " "(from_seqname, to_seqname)",
        )


def resolve_subset_regions(args, loci_df=None) -> Optional[List[Region]]:
    """Resolve CLI args → list of Regions (or None = full dataset).

    Parameters
    ----------
    args : argparse.Namespace
        Parsed CLI args (must include subset flags added by add_subset_args).
    loci_df : DataFrame or None
        Loci DataFrame with Chromosome, Start, End for sampling.
        Required if --sample-loci is used.

    Returns
    -------
    list[Region] or None
        None means run on full dataset.
    """
    regions = []

    # Explicit regions
    if getattr(args, "seqname", None):
        regions.append(Region(args.seqname))
    if getattr(args, "region", None):
        regions.append(parse_region(args.region))
    if getattr(args, "regions_file", None):
        regions.extend(load_regions_file(args.regions_file))

    # Random sampling
    sample_n = getattr(args, "sample_loci", None)
    if sample_n is not None:
        if loci_df is None or loci_df.empty:
            print("  Warning: --sample-loci specified but no loci available " "for sampling")
            return regions if regions else None

        sampled = sample_loci(
            loci_df,
            n=sample_n,
            strategy=getattr(args, "sample_strategy", "uniform_locus"),
            seed=getattr(args, "seed", 1),
            window_bp=getattr(args, "sample_bp", 0),
        )
        regions.extend(sampled)
        print(
            f"  Sampled {len(sampled)} loci (seed={args.seed}, "
            f"strategy={args.sample_strategy})"
        )

    return regions if regions else None


def write_subset_manifest(regions: List[Region], seed: int, output_path: str):
    """Write subset_regions.tsv documenting the selected regions and seed."""
    os.makedirs(os.path.dirname(output_path) or ".", exist_ok=True)
    with open(output_path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["# seed", str(seed)])
        w.writerow(["seqname", "start", "end"])
        for r in regions:
            w.writerow(
                [
                    r.seqname,
                    r.start if r.start is not None else ".",
                    r.end if r.end is not None else ".",
                ]
            )
    print(f"  Subset manifest: {output_path}")
