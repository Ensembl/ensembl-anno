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
import hashlib
import os
import re
import subprocess
import sys
import tempfile
from collections import defaultdict
from typing import TYPE_CHECKING, Optional

if TYPE_CHECKING:
    from gmb.pipeline.config import PipelineConfig, ProteinValidationConfig


# ---------------------------------------------------------------------------
# Psauron capability detection
# ---------------------------------------------------------------------------
#
# GMB constructs a psauron command from up to five flags: -i, -o, -m, -p,
# and (only when CPU mode is requested) -c. Rather than hard-coding
# behaviour to a specific psauron release, this checks -- by parsing
# `psauron --help` once -- that the installed binary actually exposes each
# flag GMB is about to use, and fails with a clear, actionable message if it
# does not, instead of either crashing on an unrecognised-argument error
# deep inside a subprocess call or (worse) silently misinterpreting a
# changed CLI.
#
# Psauron uses one bundled model checkpoint and exposes no model-selection
# option: -m specifies minimum protein length, not a model. This has been
# true across the tool's entire public release history (see the project
# report for the version-history investigation), so GMB never needs a
# version-specific code path -- only feature detection, as a safety net
# against a future release or non-standard/patched build changing the CLI.
_PSAURON_VERSION_RE = re.compile(r"PSAURON version (\S+)", re.IGNORECASE)


# ---------------------------------------------------------------------------
# Protein identity / provenance (the "compute once, consume many" contract)
# ---------------------------------------------------------------------------


def protein_sha256(seq: Optional[str]) -> Optional[str]:
    """Stable SHA-256 of a protein sequence, or None for an empty sequence.

    This is the reuse key for the whole protein-validation contract:
    identical sequences validate once and every consumer (model scoring,
    duplicate collapse, canonical selection, InterPro review) refers to the
    same result through this checksum. It deliberately does NOT replace
    transcript IDs -- both are carried in the sidecar, because transcript
    IDs are renamed late in the build (candidate ID -> final PFAL_xxxxx.tN)
    while the checksum stays constant across that rename.

    Normalised to uppercase with surrounding whitespace stripped, so the
    same peptide written with different casing/wrapping hashes identically.
    """
    if not seq:
        return None
    return hashlib.sha256(seq.strip().upper().encode()).hexdigest()


def _tool_version(path: str, args: list) -> Optional[str]:
    """Best-effort single-token version string for a tool, or None."""
    try:
        result = subprocess.run(
            [path, *args], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True
        )
    except (FileNotFoundError, OSError):
        return None
    text = (result.stdout or "") + " " + (result.stderr or "")
    match = re.search(r"(\d+\.\d+(?:\.\d+)*)", text)
    return match.group(1) if match else None


def build_validation_provenance(val_cfg: ProteinValidationConfig) -> dict:
    """Tool/database provenance recorded alongside every validation result.

    Lets a consumer decide whether a cached/sidecar result is still valid
    without rerunning anything: a different DIAMOND database or a different
    tool version invalidates reuse, and that is only detectable if the
    producing run recorded what it used.
    """
    return {
        "diamond_version": _tool_version(val_cfg.diamond_path, ["version"]),
        "psauron_version": detect_psauron_capabilities(val_cfg.psauron_path).get("version"),
        "diamond_db": val_cfg.diamond_db,
    }


# One dict per flag GMB may construct a command from: label for error
# messages, and the detect_psauron_capabilities() key that reports it.
_ALWAYS_REQUIRED_FLAGS = [
    ("-i/--input-fasta", "has_input_flag"),
    ("-o/--output-path", "has_output_flag"),
    ("-m/--minimum-length", "has_minimum_length"),
    ("-p/--protein", "has_protein_flag"),
]
_CPU_FLAG = ("-c/--use-cpu", "has_cpu_flag")


def detect_psauron_capabilities(psauron_path: str) -> dict:
    """Run ``psauron --help`` once and detect its version and flag availability.

    Returns a dict with ``found`` (bool -- False only if the executable
    itself could not be run), ``returncode`` (int or None if not found),
    ``version`` (str or None -- not parseable from every conceivable build,
    handled as "unknown" by callers rather than raising), ``help_text``,
    and one ``has_<flag>`` bool per flag GMB can use (``has_input_flag``,
    ``has_output_flag``, ``has_minimum_length``, ``has_protein_flag``,
    ``has_cpu_flag``).

    Never raises for a missing binary, a non-zero exit, or an unparseable
    version string -- this is pure detection; callers (check_dependencies)
    decide what to do with the result.
    """
    empty = {
        "found": False,
        "returncode": None,
        "version": None,
        "help_text": "",
        "has_input_flag": False,
        "has_output_flag": False,
        "has_minimum_length": False,
        "has_protein_flag": False,
        "has_cpu_flag": False,
    }
    try:
        result = subprocess.run(
            [psauron_path, "--help"],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )
    except FileNotFoundError:
        return empty

    # psauron prints its version banner to stdout even when --help alone
    # triggers a non-zero exit in some builds -- check both streams rather
    # than assuming which one carries it.
    help_text = (result.stdout or "") + "\n" + (result.stderr or "")
    version_match = _PSAURON_VERSION_RE.search(help_text)

    return {
        **empty,
        "found": True,
        "returncode": result.returncode,
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

    # Capability (not version-number) detection: run `psauron --help` once
    # and verify the installed binary actually exposes the flags GMB is
    # about to construct a command line from, rather than assuming any
    # particular release. This replaces a separate existence-only check --
    # a missing executable is just one more capability-detection outcome.
    caps = detect_psauron_capabilities(val_cfg.psauron_path)
    if not caps["found"]:
        print(
            f"ERROR: Protein validation enabled but psauron binary '{val_cfg.psauron_path}' not found in PATH.",
            file=sys.stderr,
        )
        sys.exit(1)
    if caps["returncode"] != 0:
        print(
            f"ERROR: '{val_cfg.psauron_path} --help' exited with status "
            f"{caps['returncode']} instead of 0. Cannot verify this psauron "
            "installation's CLI is usable.",
            file=sys.stderr,
        )
        sys.exit(1)

    # -c/--use-cpu is only ever added to the real invocation (see
    # batch_score_proteins) when psauron_use_cpu is set, so it is only
    # required here in that case.
    required_flags = list(_ALWAYS_REQUIRED_FLAGS)
    if val_cfg.psauron_use_cpu:
        required_flags.append(_CPU_FLAG)

    missing = [label for label, key in required_flags if not caps[key]]
    if missing:
        print(
            f"ERROR: psauron at '{val_cfg.psauron_path}' "
            f"(detected version: {caps['version'] or 'unknown, could not parse --help banner'}) "
            f"does not expose the flag(s) GMB requires: {', '.join(missing)}. "
            "This installed psauron's CLI differs from every publicly released "
            "version investigated -- refusing to guess at a compatible command "
            "line. See support_scripts/gmb docs for the capability-detection "
            "rationale.",
            file=sys.stderr,
        )
        sys.exit(1)
    print(f"  Detected psauron version: {caps['version'] or 'unknown'}")

    if val_cfg.diamond_db is None:
        if val_cfg.diamond_weight > 0:
            print(
                "ERROR: protein_validation.diamond_db must be set when diamond_weight > 0.",
                file=sys.stderr,
            )
            sys.exit(1)
    elif not os.path.exists(val_cfg.diamond_db) and not os.path.exists(
        f"{val_cfg.diamond_db}.dmnd"
    ):
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
        not available (no hit, sequence too short, empty protein, etc.), plus
        provenance: ``protein_sha256``, ``validation_status``,
        ``validation_reason``, tool/database versions, and
        ``protein_validation_source`` (``"calculated"``/``"reused"``/``None``)
        with ``protein_validation_reused_from`` recording which key's
        calculation is being reused when the sequence is a duplicate of one
        already scored in this same batch.
    """
    val_cfg = config.protein_validation
    provenance = build_validation_provenance(val_cfg)

    def _empty_detail(seq: str, status: str, reason: str) -> dict:
        return {
            "protein_sha256": protein_sha256(seq),
            "diamond_hit": None,
            "diamond_pident": None,
            "diamond_qcov": None,
            "diamond_scov": None,
            "diamond_bitscore": None,
            "diamond_evalue": None,
            "psauron_score": None,
            "protein_length": len(seq) if seq else 0,
            "validation_status": status,
            "validation_reason": reason,
            # No DIAMOND/psauron invocation happened for this key at all
            # (validation disabled, or no translatable protein), so
            # "calculated vs reused" does not apply -- see the real
            # computation below for the populated values.
            "protein_validation_source": None,
            "protein_validation_reused_from": None,
            **provenance,
        }

    if not val_cfg.enabled:
        return (
            {k: 0.0 for k in protein_dict},
            {
                k: _empty_detail(seq, "not_run", "protein_validation.enabled is false")
                for k, seq in protein_dict.items()
            },
        )

    if not protein_dict:
        return ({}, {})

    # Deduplicate strictly identical sequences -- this is the "compute once"
    # half of the compute-once/consume-many contract: DIAMOND and Psauron
    # each see one copy of each distinct sequence, and every transcript
    # sharing that sequence receives the same result object.
    # We map sequence -> list of keys that produced it.
    seq_to_keys = defaultdict(list)
    for k, seq in protein_dict.items():
        if seq and len(seq) > 0:
            seq_to_keys[seq].append(k)

    if not seq_to_keys:
        return (
            {k: 0.0 for k in protein_dict},
            {
                k: _empty_detail(seq, "no_protein", "empty or missing translated sequence")
                for k, seq in protein_dict.items()
            },
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
        # Real psauron output: a CSV file whose first two non-empty lines
        # are the invoked command and a "psauron score: <mean>" summary
        # line, THEN a header row `description,psauron_is_protein,
        # in-frame_score`, followed by one data row per input sequence.
        #
        # The input here is always a protein FASTA, never spliced CDS, so
        # -p/--protein is required for a correct score; -m/--minimum-length
        # sets the minimum protein length, not a model (psauron has no
        # model-selection option -- see detect_psauron_capabilities above).
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

        # Fail loudly rather than silently returning 0.0 for every sequence
        # at each of these steps -- an earlier version of this code matched
        # zero real data rows (wrong delimiter, no header handling) and
        # every psauron score was silently always 0.0 regardless of the
        # real result. See detect_psauron_capabilities()/check_dependencies()
        # for the corresponding pre-flight CLI-flag check.
        if not os.path.exists(psauron_out):
            sys.exit(
                f"ERROR: psauron did not produce an output file at '{psauron_out}'. "
                "Refusing to silently score every protein 0.0."
            )
        with open(psauron_out, newline="") as fh:
            lines = fh.readlines()
        # Locate the header row rather than assuming a fixed line count,
        # since the number of preamble lines is not part of psauron's
        # documented/stable output contract.
        header_idx = next(
            (i for i, line in enumerate(lines) if line.startswith("description,")), None
        )
        if header_idx is None:
            sys.exit(
                "ERROR: could not find a 'description,...' header row in psauron's "
                f"output '{psauron_out}'. The installed psauron's CSV output format "
                "differs from every version investigated -- refusing to silently "
                "produce all-zero protein-validation scores."
            )

        reader = csv.DictReader(lines[header_idx:])
        required_columns = ("description", "in-frame_score")
        missing_columns = [c for c in required_columns if c not in (reader.fieldnames or [])]
        if missing_columns:
            sys.exit(
                f"ERROR: psauron output CSV header is missing column(s) {missing_columns} "
                f"(found: {reader.fieldnames}). The installed psauron's CSV output format "
                "differs from every version investigated -- refusing to silently "
                "produce all-zero protein-validation scores."
            )

        psauron_scores = {}
        malformed_rows = 0
        for row in reader:
            qseqid = row.get("description")
            score_str = row.get("in-frame_score")
            if not qseqid or not score_str:
                malformed_rows += 1
                continue
            try:
                psauron_scores[qseqid] = float(score_str)
            except ValueError:
                malformed_rows += 1
        if malformed_rows:
            print(
                f"WARNING: {malformed_rows} psauron output row(s) were malformed "
                "(missing description/score, or non-numeric score) and were skipped.",
                file=sys.stderr,
            )
        submitted_sids = set(seq_to_id.values())
        if not (psauron_scores.keys() & submitted_sids):
            sys.exit(
                f"ERROR: psauron produced no score rows matching any of the "
                f"{len(submitted_sids)} submitted sequence ID(s) (got: "
                f"{sorted(psauron_scores)[:5]}...). Refusing to silently score "
                "every protein 0.0."
            )

        # Calculate combined scores
        final_scores = {}
        details = {}

        for seq, keys in seq_to_keys.items():
            sid = seq_to_id[seq]
            d_score = diamond_scores.get(sid, 0.0)
            p_score = psauron_scores.get(sid, 0.0)

            comb_score = (d_score * val_cfg.diamond_weight) + (p_score * val_cfg.psauron_weight)
            has_diamond = sid in diamond_detail_by_sid
            has_psauron = sid in psauron_scores
            if has_diamond and has_psauron:
                status, reason = "ok", ""
            elif has_psauron:
                status, reason = "no_diamond_hit", "no DIAMOND hit passing coverage gates"
            elif has_diamond:
                status, reason = "no_psauron_score", "psauron returned no score for this sequence"
            else:
                status, reason = "no_evidence", "neither DIAMOND nor psauron returned a result"

            detail = _empty_detail(seq, status, reason)
            detail["protein_length"] = len(seq)
            if has_diamond:
                detail.update(diamond_detail_by_sid[sid])
            if has_psauron:
                detail["psauron_score"] = psauron_scores[sid]
            # Every transcript sharing this sequence gets the *same* result;
            # copied per key so a later per-transcript annotation (e.g. the
            # final transcript ID) cannot leak across transcripts. Exactly
            # one key per unique sequence is deterministically designated
            # "calculated" (the lexicographically-first key -- an arbitrary
            # but stable and reproducible choice); every other key sharing
            # that sequence is "reused", with reused_from recording which
            # key's result it is sharing. This makes the compute-once/
            # consume-many contract auditable per transcript rather than
            # only inferable from a shared protein_sha256.
            calculated_key = min(keys)
            for k in keys:
                final_scores[k] = comb_score
                row = dict(detail)
                if k == calculated_key:
                    row["protein_validation_source"] = "calculated"
                    row["protein_validation_reused_from"] = None
                else:
                    row["protein_validation_source"] = "reused"
                    row["protein_validation_reused_from"] = calculated_key
                details[k] = row

    # Re-apply defaults for missing/empty
    for k, seq in protein_dict.items():
        if k not in final_scores:
            final_scores[k] = 0.0
        if k not in details:
            details[k] = _empty_detail(seq, "no_protein", "empty or missing translated sequence")

    return final_scores, details
