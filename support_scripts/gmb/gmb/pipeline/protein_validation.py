#!/usr/bin/env python3
"""Protein-coding validation stage for Gene Model Builder.

Runs DIAMOND and Psauron against candidate translated sequences
in batch. Deduplicates strict sequences within run to optimise throughput,
and returns a score dict mapping sequences or structural hashes to their results.
"""

from __future__ import annotations

import os
import subprocess
import sys
import tempfile
from collections import defaultdict
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from gmb.pipeline.config import PipelineConfig, ProteinValidationConfig


def check_dependencies(val_cfg: ProteinValidationConfig) -> None:
    """Ensure tools and databases are present. Exits if missing.

    Parameters
    ----------
    val_cfg : ProteinValidationConfig
        Protein validation configuration section.
    """
    if not val_cfg.enabled:
        return

    try:
        subprocess.run(
            [val_cfg.diamond_path, "help"],
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
    except FileNotFoundError:
        print(
            f"ERROR: Protein validation enabled but diamond binary '{val_cfg.diamond_path}' not found in PATH.",
            file=sys.stderr,
        )
        sys.exit(1)

    try:
        subprocess.run(
            [val_cfg.psauron_path, "--help"],
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
    except FileNotFoundError:
        print(
            f"ERROR: Protein validation enabled but psauron binary '{val_cfg.psauron_path}' not found in PATH.",
            file=sys.stderr,
        )
        sys.exit(1)

    if not os.path.exists(val_cfg.diamond_db) and not os.path.exists(f"{val_cfg.diamond_db}.dmnd"):
        print(f"ERROR: DIAMOND database '{val_cfg.diamond_db}' not found.", file=sys.stderr)
        sys.exit(1)


def batch_score_proteins(
    protein_dict: dict[str, str],
    config: PipelineConfig,
) -> dict[str, float]:
    """Run DIAMOND and Psauron in batch on unique sequences.

    Parameters
    ----------
    protein_dict : dict
        Mapping of arbitrary keys (e.g. transcript ID or struct hash) to
        protein sequences (strings).
    config : PipelineConfig

    Returns
    -------
    dict
        Mapping of the same keys to float ``protein_coding_score``.
    """
    val_cfg = config.protein_validation
    if not val_cfg.enabled or not protein_dict:
        return {k: 0.0 for k in protein_dict}

    # Deduplicate strictly identical sequences.
    # We map sequence -> list of keys that produced it.
    seq_to_keys = defaultdict(list)
    for k, seq in protein_dict.items():
        if seq and len(seq) > 0:
            seq_to_keys[seq].append(k)

    if not seq_to_keys:
        return {k: 0.0 for k in protein_dict}

    with tempfile.TemporaryDirectory() as td:
        fasta_path = os.path.join(td, "candidates.fa")
        # Ensure we write a clean fasta with safe temporary identifiers
        seq_to_id = {}
        with open(fasta_path, "w") as fh:
            for i, seq in enumerate(seq_to_keys.keys()):
                sid = f"seq_{i}"
                seq_to_id[seq] = sid
                fh.write(f">{sid}\n{seq}\n")

        # -----------------------------
        # DIAMOND
        # -----------------------------
        diamond_out = os.path.join(td, "diamond.tsv")
        # max-target-seqs 1 to just get back the best hit per query
        cmd_d = [
            val_cfg.diamond_path,
            "blastp",
            "-d",
            val_cfg.diamond_db,
            "-q",
            fasta_path,
            "-o",
            diamond_out,
            "--outfmt",
            "6",
            "qseqid",
            "sseqid",
            "pident",
            "length",
            "bitscore",
            "evalue",
            "--max-target-seqs",
            "1",
            "--quiet",
        ]

        try:
            subprocess.run(cmd_d, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        except subprocess.CalledProcessError as e:
            print(f"DIAMOND extraction failed: {e.stderr}", file=sys.stderr)
            sys.exit(1)

        diamond_scores = {}
        if os.path.exists(diamond_out):
            with open(diamond_out) as fh:
                for line in fh:
                    parts = line.strip().split("\t")
                    if len(parts) >= 5:
                        qseqid = parts[0]
                        bitscore = float(parts[4])
                        # Normalise between 0 and 1, assuming max bitscore around 1000
                        norm_bs = min(1.0, bitscore / 1000.0)
                        diamond_scores[qseqid] = norm_bs

        # -----------------------------
        # Psauron
        # -----------------------------
        # Output is commonly a JSON or tabular format depending on Psauron version.
        # Assume it can dump a json score map for ids, or simple TSV
        psauron_out = os.path.join(td, "psauron.txt")
        cmd_p = [
            val_cfg.psauron_path,
            "-i",
            fasta_path,
            "-o",
            psauron_out,
            "-m",
            val_cfg.psauron_model,
        ]

        try:
            subprocess.run(cmd_p, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        except subprocess.CalledProcessError as e:
            # Just fail loudly if Psauron is missing somehow
            print(f"Psauron extraction failed: {e.stderr}", file=sys.stderr)
            sys.exit(1)

        psauron_scores = {}
        if os.path.exists(psauron_out):
            with open(psauron_out) as fh:
                for line in fh:
                    if line.startswith("#") or not line.strip():
                        continue
                    parts = line.strip().split()
                    if len(parts) >= 2:
                        qseqid = parts[0]
                        try:
                            score = float(parts[1])
                            # Assume Psauron output is already a probability/score between 0 and 1
                            psauron_scores[qseqid] = score
                        except ValueError:
                            pass

        # Calculate combined scores
        final_scores = {}

        for seq, keys in seq_to_keys.items():
            sid = seq_to_id[seq]
            d_score = diamond_scores.get(sid, 0.0)
            p_score = psauron_scores.get(sid, 0.0)

            comb_score = (d_score * val_cfg.diamond_weight) + (p_score * val_cfg.psauron_weight)
            for k in keys:
                final_scores[k] = comb_score

    # Re-apply defaults for missing/empty
    for k in protein_dict:
        if k not in final_scores:
            final_scores[k] = 0.0

    return final_scores
