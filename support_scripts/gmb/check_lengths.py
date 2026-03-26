import pandas as pd
import pyranges as pr


def check_lengths(path, name):
    print(f"Checking {name} ({path})...")
    df = pr.read_gtf(path).df if path.endswith("gtf") else pr.read_gff3(path).df
    if "Feature" in df.columns:
        df = df[df["Feature"] == "transcript"]  # GTF usually has transcript lines

    if df.empty:
        # Reconstruct from exons
        df = pr.read_gtf(path).df if path.endswith("gtf") else pr.read_gff3(path).df
        df = df[df["Feature"] == "exon"]
        if "transcript_id" not in df.columns:
            df["transcript_id"] = df["Parent"] if "Parent" in df.columns else "Unknown"

        df = df.groupby("transcript_id").agg({"Start": "min", "End": "max"})

    length = df["End"] - df["Start"]
    print(f"  Count: {len(length)}")
    print(f"  Min Length: {length.min()}")
    print(f"  Max Length: {length.max()}")
    print(f"  Avg Length: {length.mean()}")
    print(f"  >10kb count: {sum(length > 10000)}")
    print(f"  >50kb count: {sum(length > 50000)}")


check_lengths("scallop_annotation.gtf", "Scallop")
check_lengths("stringtie_annotation.gtf", "StringTie")
check_lengths("helixer_remapped.gff3", "Helixer")
