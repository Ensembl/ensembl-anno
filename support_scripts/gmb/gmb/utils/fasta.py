"""FASTA read/write/validation helpers shared across the pipeline."""

from __future__ import annotations

import gzip
from typing import IO


def load_genome(fasta_path: str) -> dict[str, str]:
    """Load a FASTA file into a dict of {chrom: sequence (uppercase)}.

    Handles both plain and gzip-compressed files.
    """
    genome: dict[str, str] = {}
    current_name = None
    chunks: list[str] = []
    opener = gzip.open if fasta_path.endswith(".gz") else open
    with opener(fasta_path, "rt") as fh:
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


def parse_fasta_ids(fasta_path: str) -> list[str]:
    """Extract record IDs (first whitespace-delimited token after '>') from FASTA."""
    import os

    ids: list[str] = []
    if not os.path.exists(fasta_path):
        return ids
    with open(fasta_path) as fh:
        for line in fh:
            if line.startswith(">"):
                ids.append(line[1:].split()[0].strip())
    return ids


def parse_fasta_records(fasta_path: str) -> dict[str, str]:
    """Parse FASTA into {id: sequence} dict."""
    import os

    records: dict[str, str] = {}
    if not os.path.exists(fasta_path):
        return records
    current_id = None
    parts: list[str] = []
    with open(fasta_path) as fh:
        for line in fh:
            if line.startswith(">"):
                if current_id is not None:
                    records[current_id] = "".join(parts)
                current_id = line[1:].split()[0].strip()
                parts = []
            else:
                parts.append(line.strip())
    if current_id is not None:
        records[current_id] = "".join(parts)
    return records


def write_seq(fh: IO[str], seq: str, line_width: int = 80) -> None:
    """Write a sequence in wrapped FASTA format."""
    for i in range(0, len(seq), line_width):
        fh.write(seq[i : i + line_width] + "\n")
