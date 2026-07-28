#!/usr/bin/env python3
"""Protein-coding validation stage for Gene Model Builder.

Runs DIAMOND and Psauron against candidate translated sequences
in batch. Deduplicates strict sequences within run to optimise throughput,
and returns both a combined score dict (for the existing scoring gate) and
a per-key detail dict preserving the individual DIAMOND/Psauron component
values (for reporting and canonical-transcript selection).
"""

from __future__ import annotations

import csv
import os
import re
import subprocess
import sys
import tempfile
from collections import defaultdict
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from gmb.pipeline.config import PipelineConfig, ProteinValidationConfig


# ---------------------------------------------------------------------------
# Psauron capability detection
# ---------------------------------------------------------------------------
#
# GMB constructs a psauron command using exactly five flags: -i, -o, -m, -p,
# -c (see batch_score_proteins below). Rather than hard-coding behaviour to
# a specific psauron release, this checks -- by parsing `psauron --help` --
# that the installed binary actually exposes each of those flags, and fails
# with a clear, actionable message if it does not, instead of either
# crashing on an unrecognised-argument error deep inside a subprocess call
# or (worse) silently misinterpreting a changed CLI.
#
# Investigated (2026-07) against the psauron GitHub repository
# (salzberg-lab/PSAURON) release history, from the first public tag v1.0.0
# (2024-07) through the latest at the time, v1.1.3 (2026-05): the CLI
# accepted by GMB (-i/-o/-m/-p/-c) and the CSV output shape it depends on
# (a `description,psauron_is_protein,in-frame_score` header, in -p/protein
# mode) have been identical in every release across that entire span. There
# is no evidence a `--model`/model-selection flag has ever existed in
# psauron at any point -- the tool has always shipped exactly one embedded
# TCN checkpoint, so a model-selection argument would not correspond to
# anything the tool actually does. `-m` has always meant `--minimum-length`
# (an amino-acid length filter), never a model selector. This capability
# check exists as a safety net against a *future* release or a
# non-standard/patched build changing that, not because such a change is
# currently known to exist -- see the project report for the full
# version-history investigation this is based on.
_PSAURON_VERSION_RE = re.compile(r"PSAURON version (\S+)", re.IGNORECASE)


def detect_psauron_capabilities(psauron_path: str) -> dict:
    """Parse ``psauron --help`` to detect its version and flag availability.

    Returns a dict with ``version`` (str or None -- not parseable from
    every conceivable build, handled as "unknown" by callers rather than
    raising), ``help_text``, and one ``has_<flag>`` bool per flag GMB
    actually uses (``has_input_flag``, ``has_output_flag``,
    ``has_minimum_length``, ``has_protein_flag``, ``has_cpu_flag``).

    Never raises for a missing/unparseable version string or a missing
    flag -- this is pure detection; callers (check_dependencies) decide
    whether an absent flag is fatal.
    """
    try:
        result = subprocess.run(
            [psauron_path, "--help"],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )
    except FileNotFoundError:
        return {
            "version": None,
            "help_text": "",
            "has_input_flag": False,
            "has_output_flag": False,
            "has_minimum_length": False,
            "has_protein_flag": False,
            "has_cpu_flag": False,
        }

    # psauron prints its version banner to stdout even when --help alone
    # triggers a non-zero exit in some builds -- check both streams rather
    # than assuming which one carries it.
    help_text = (result.stdout or "") + "\n" + (result.stderr or "")

    version_match = _PSAURON_VERSION_RE.search(help_text)

    return {
        "version": version_match.group(1) if version_match else None,
        "help_text": help_text,
        "has_input_flag": bool(re.search(r"-i[,\s]|--input-fasta", help_text)),
        "has_output_flag": bool(re.search(r"-o[,\s]|--output-path", help_text)),
        "has_minimum_length": bool(re.search(r"-m[,\s]|--minimum-length", help_text)),
        "has_protein_flag": bool(re.search(r"-p,\s*--protein|--protein\b", help_text)),
        "has_cpu_flag": bool(re.search(r"-c,\s*--use-cpu|--use-cpu\b", help_text)),
    }


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

    # Capability (not version-number) detection: verify the installed
    # psauron actually exposes the flags GMB is about to construct a
    # command line from, rather than assuming any particular release.
    caps = detect_psauron_capabilities(val_cfg.psauron_path)
    required_flags = [
        ("-i/--input-fasta", "has_input_flag"),
        ("-o/--output-path", "has_output_flag"),
        ("-m/--minimum-length", "has_minimum_length"),
        ("-p/--protein", "has_protein_flag"),
        ("-c/--use-cpu", "has_cpu_flag"),
    ]
    missing = [label for label, key in required_flags if not caps[key]]
    if missing:
        print(
            f"ERROR: psauron at '{val_cfg.psauron_path}' "
            f"(detected version: {caps['version'] or 'unknown, could not parse --help banner'}) "
            f"does not expose the flag(s) GMB requires: {', '.join(missing)}. "
            "This installed psauron's CLI differs from every publicly released "
            "version investigated (v1.0.0 through v1.1.3) -- refusing to guess "
            "at a compatible command line. See support_scripts/gmb docs for the "
            "capability-detection rationale.",
            file=sys.stderr,
        )
        sys.exit(1)
    print(f"  Detected psauron version: {caps['version'] or 'unknown'}")

    if not os.path.exists(val_cfg.diamond_db) and not os.path.exists(f"{val_cfg.diamond_db}.dmnd"):
        print(f"ERROR: DIAMOND database '{val_cfg.diamond_db}' not found.", file=sys.stderr)
        sys.exit(1)


# DIAMOND outfmt 6 column order requested below -- kept in sync with the
# parsing loop in batch_score_proteins(). qcovhsp/scovhsp (query/target
# coverage as a percentage) were added so canonical-transcript selection
# (see gmb.pipeline.canonical_selection) can use coverage rather than raw
# bitscore, which scales with protein length.
_DIAMOND_OUTFMT_COLUMNS = [
    "qseqid",
    "sseqid",
    "pident",
    "length",
    "bitscore",
    "evalue",
    "qcovhsp",
    "scovhsp",
]


def batch_score_proteins(
    protein_dict: dict[str, str],
    config: PipelineConfig,
) -> tuple[dict[str, float], dict[str, dict]]:
    """Run DIAMOND and Psauron in batch on unique sequences.

    Parameters
    ----------
    protein_dict : dict
        Mapping of arbitrary keys (e.g. transcript ID or struct hash) to
        protein sequences (strings).
    config : PipelineConfig

    Returns
    -------
    tuple of (dict, dict)
        ``(combined_scores, details)``. ``combined_scores`` maps each key to
        the float ``protein_coding_score`` (unchanged shape/semantics from
        before -- this is what ``gmb.pipeline.scoring`` consumes).
        ``details`` maps each key to a dict of the individual DIAMOND/Psauron
        values (``diamond_hit``, ``diamond_pident``, ``diamond_qcov``,
        ``diamond_scov``, ``diamond_bitscore``, ``diamond_evalue``,
        ``psauron_score``, ``protein_length``), each ``None``/``0`` when
        not available (no hit, sequence too short, empty protein, etc.).
    """
    val_cfg = config.protein_validation

    def _empty_detail(seq: str) -> dict:
        return {
            "diamond_hit": None,
            "diamond_pident": None,
            "diamond_qcov": None,
            "diamond_scov": None,
            "diamond_bitscore": None,
            "diamond_evalue": None,
            "psauron_score": None,
            "protein_length": len(seq) if seq else 0,
        }

    if not val_cfg.enabled or not protein_dict:
        return (
            {k: 0.0 for k in protein_dict},
            {k: _empty_detail(seq) for k, seq in protein_dict.items()},
        )

    # Deduplicate strictly identical sequences.
    # We map sequence -> list of keys that produced it.
    seq_to_keys = defaultdict(list)
    for k, seq in protein_dict.items():
        if seq and len(seq) > 0:
            seq_to_keys[seq].append(k)

    if not seq_to_keys:
        return (
            {k: 0.0 for k in protein_dict},
            {k: _empty_detail(seq) for k, seq in protein_dict.items()},
        )

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
            *_DIAMOND_OUTFMT_COLUMNS,
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
        diamond_detail_by_sid = {}
        if os.path.exists(diamond_out):
            with open(diamond_out) as fh:
                for line in fh:
                    parts = line.rstrip("\n").split("\t")
                    if len(parts) < len(_DIAMOND_OUTFMT_COLUMNS):
                        continue
                    row = dict(zip(_DIAMOND_OUTFMT_COLUMNS, parts))
                    qseqid = row["qseqid"]
                    bitscore = float(row["bitscore"])
                    qcov = float(row["qcovhsp"])
                    scov = float(row["scovhsp"])
                    if qcov < val_cfg.diamond_min_query_coverage:
                        continue
                    if scov < val_cfg.diamond_min_target_coverage:
                        continue
                    # Normalise bitscore between 0 and 1, assuming max
                    # bitscore around 1000 (a bounded, if crude, way to keep
                    # a single very-long-protein hit from dominating the
                    # combined score -- see canonical_selection for a
                    # coverage/identity-based alternative that does not
                    # scale with protein length at all).
                    norm_bs = min(1.0, bitscore / 1000.0)
                    diamond_scores[qseqid] = norm_bs
                    diamond_detail_by_sid[qseqid] = {
                        "diamond_hit": row["sseqid"],
                        "diamond_pident": float(row["pident"]),
                        "diamond_qcov": qcov,
                        "diamond_scov": scov,
                        "diamond_bitscore": bitscore,
                        "diamond_evalue": float(row["evalue"]),
                    }

        # -----------------------------
        # Psauron
        # -----------------------------
        # Real psauron 1.0.8 output (verified against the installed binary,
        # not assumed): a CSV file whose first two non-empty lines are the
        # invoked command and a "psauron score: <mean>" summary line, THEN a
        # header row `description,psauron_is_protein,in-frame_score`,
        # followed by one data row per input sequence. An earlier version of
        # this code assumed a simple whitespace-separated `id score` format
        # and split on whitespace with no header handling -- that silently
        # matched zero real data rows (comma-separated, no whitespace), so
        # every psauron score was always 0.0 regardless of the real result.
        #
        # psauron 1.0.8 also has no --model flag (checked via --help); the
        # relevant flags are -m/--minimum-length (aa) and -p/--protein (the
        # input here is always a protein FASTA, never spliced CDS, so -p is
        # required for a correct score -- omitting it was a second bug).
        psauron_out = os.path.join(td, "psauron.csv")
        cmd_p = [
            val_cfg.psauron_path,
            "-i",
            fasta_path,
            "-o",
            psauron_out,
            "-m",
            str(val_cfg.psauron_min_length),
            "-p",
        ]
        if val_cfg.psauron_use_cpu:
            cmd_p.append("-c")

        try:
            subprocess.run(cmd_p, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        except subprocess.CalledProcessError as e:
            # Just fail loudly if Psauron is missing somehow
            print(f"Psauron extraction failed: {e.stderr}", file=sys.stderr)
            sys.exit(1)

        psauron_scores = {}
        if os.path.exists(psauron_out):
            with open(psauron_out, newline="") as fh:
                lines = fh.readlines()
            # Locate the header row rather than assuming a fixed line
            # count, since the number of preamble lines is not part of
            # psauron's documented/stable output contract.
            header_idx = None
            for i, line in enumerate(lines):
                if line.startswith("description,"):
                    header_idx = i
                    break
            if header_idx is not None:
                reader = csv.DictReader(lines[header_idx:])
                if reader.fieldnames and "in-frame_score" not in reader.fieldnames:
                    # Fail loudly rather than silently returning 0.0 for
                    # every sequence -- this exact failure mode (a changed
                    # column name matching nothing, no error, every score
                    # quietly 0.0) was the original bug this module fixed.
                    # See detect_psauron_capabilities()/check_dependencies()
                    # for the corresponding pre-flight CLI-flag check.
                    sys.exit(
                        "ERROR: psauron output CSV header does not contain the "
                        f"expected 'in-frame_score' column (found: {reader.fieldnames}). "
                        "The installed psauron's CSV output format differs from every "
                        "version investigated (v1.0.0-v1.1.3) -- refusing to silently "
                        "produce all-zero protein-validation scores."
                    )
                for row in reader:
                    qseqid = row.get("description")
                    score_str = row.get("in-frame_score")
                    if qseqid and score_str:
                        try:
                            psauron_scores[qseqid] = float(score_str)
                        except ValueError:
                            pass

        # Calculate combined scores
        final_scores = {}
        details = {}

        for seq, keys in seq_to_keys.items():
            sid = seq_to_id[seq]
            d_score = diamond_scores.get(sid, 0.0)
            p_score = psauron_scores.get(sid, 0.0)

            comb_score = (d_score * val_cfg.diamond_weight) + (p_score * val_cfg.psauron_weight)
            detail = _empty_detail(seq)
            detail["protein_length"] = len(seq)
            if sid in diamond_detail_by_sid:
                detail.update(diamond_detail_by_sid[sid])
            if sid in psauron_scores:
                detail["psauron_score"] = psauron_scores[sid]
            for k in keys:
                final_scores[k] = comb_score
                details[k] = detail

    # Re-apply defaults for missing/empty
    for k, seq in protein_dict.items():
        if k not in final_scores:
            final_scores[k] = 0.0
        if k not in details:
            details[k] = _empty_detail(seq)

    return final_scores, details
