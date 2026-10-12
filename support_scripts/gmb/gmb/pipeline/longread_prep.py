"""Automatic preparation of per-read long-read alignments for a GMB run.

Per-read alignments (one GTF transcript per aligned read, strand already set
upstream) are collapsed into consensus transcript models by the existing
``gmb-longread-consensus`` command before preflight and build. This module only
validates the input, decides whether a previous consensus can be reused, runs
that command and records provenance; it implements no consensus logic itself.

Layout inside the run directory::

    longread_consensus/
        minimap2_consensus.gtf               the evidence passed to preflight and build
        longread_preparation.json            input identity, settings, versions, status
        longread_consensus_summary.{json,tsv}, longread_consensus_dropped.tsv,
        longread_consensus_run_manifest.json  (written by gmb-longread-consensus)
    logs/longread_consensus.log

A previous consensus is reused only when the input SHA-256, the resolved
consensus settings and the GMB version all match and the consensus file is
unchanged; otherwise it is rebuilt from an emptied directory. The preparation
manifest is written last, only on success, so a failed or interrupted run leaves
no manifest and can never be reused.
"""

from __future__ import annotations

import dataclasses
import hashlib
import json
import os
import shutil
import subprocess
import time

#: Build presets whose long-read consensus preset is known. Any other build
#: preset must name a consensus preset or config explicitly.
CONSENSUS_PRESET_FOR_BUILD_PRESET = {
    "fungi": "fungi",
    "fungi_latest_settings": "fungi",
}

#: gmb-longread-consensus output signature (source column / attribute).
CONSENSUS_SOURCE = "Minimap2Consensus"
CONSENSUS_ATTRIBUTE = 'read_support "'

PREPARATION_MANIFEST = "longread_preparation.json"
CONSENSUS_GTF = "minimap2_consensus.gtf"


class LongReadPreparationError(RuntimeError):
    """The long-read input is unusable or the consensus step failed."""


def _genome_seqnames(genome: str) -> set[str]:
    names = set()
    with open(genome) as fh:
        for line in fh:
            if line.startswith(">"):
                names.add(line[1:].split()[0])
    return names


def scan_alignments(path: str, genome_names: set[str]) -> dict:
    """One streaming pass: SHA-256 plus structural, strand and naming statistics.

    Never holds the file in memory; only counters and the set of sequence names.
    """
    h = hashlib.sha256()
    stats = {
        "lines": 0,
        "data_rows": 0,
        "exon_rows": 0,
        "transcript_rows": 0,
        "malformed_rows": 0,
        "first_malformed_line": None,
        "exon_rows_without_transcript_id": 0,
        "invalid_strand_rows": 0,
        "unknown_seqname_rows": 0,
        "consensus_signature_rows": 0,
        "seqnames": set(),
    }
    with open(path, "rb") as fh:
        for raw in fh:
            h.update(raw)
            stats["lines"] += 1
            if raw.startswith(b"#") or not raw.strip():
                continue
            f = raw.decode("utf-8", "replace").rstrip("\n").split("\t")
            if len(f) != 9:
                stats["malformed_rows"] += 1
                stats["first_malformed_line"] = stats["first_malformed_line"] or stats["lines"]
                continue
            stats["data_rows"] += 1
            stats["seqnames"].add(f[0])
            if f[0] not in genome_names:
                stats["unknown_seqname_rows"] += 1
            if f[1] == CONSENSUS_SOURCE or CONSENSUS_ATTRIBUTE in f[8]:
                stats["consensus_signature_rows"] += 1
            if f[6] not in ("+", "-"):
                stats["invalid_strand_rows"] += 1
            if f[2] == "exon":
                stats["exon_rows"] += 1
                if 'transcript_id "' not in f[8]:
                    stats["exon_rows_without_transcript_id"] += 1
            elif f[2] == "transcript":
                stats["transcript_rows"] += 1
    stats["sha256"] = h.hexdigest()
    stats["size_bytes"] = os.path.getsize(path)
    stats["seqnames"] = sorted(stats["seqnames"])
    return stats


def validate_alignment_stats(
    path: str, stats: dict, max_unstranded_fraction: float, max_unknown_seqid_fraction: float
) -> list[str]:
    """Reasons the per-read input is unusable (empty list = usable)."""
    problems = []
    if stats["malformed_rows"]:
        problems.append(
            f"{stats['malformed_rows']:,} row(s) do not have the 9 GTF columns "
            f"(first at line {stats['first_malformed_line']}); is this a GTF file?"
        )
    if not stats["exon_rows"]:
        problems.append("no exon rows: not a per-read alignment GTF")
    if stats["exon_rows_without_transcript_id"]:
        problems.append(
            f"{stats['exon_rows_without_transcript_id']:,} exon row(s) lack "
            "a transcript_id attribute"
        )
    if stats["consensus_signature_rows"]:
        problems.append(
            "this file is gmb-longread-consensus output (source "
            f"'{CONSENSUS_SOURCE}' / read_support attribute), not per-read alignments. "
            "Collapsing it again would discard most models: supply it as long-read "
            "transcript models (--long-read-models) instead."
        )
    n = max(stats["data_rows"], 1)
    unstranded = stats["invalid_strand_rows"] / n
    if unstranded > max_unstranded_fraction:
        problems.append(
            f"{stats['invalid_strand_rows']:,} of {stats['data_rows']:,} rows "
            f"({unstranded:.1%}) have no + / - strand (limit "
            f"{max_unstranded_fraction:.0%}). Long-read strand must be set upstream "
            "(e.g. from the Minimap2 ts tag); GMB does not infer it."
        )
    unknown = stats["unknown_seqname_rows"] / n
    if unknown > max_unknown_seqid_fraction:
        problems.append(
            f"{stats['unknown_seqname_rows']:,} rows ({unknown:.1%}) are on sequences "
            f"absent from the genome FASTA (limit {max_unknown_seqid_fraction:.0%}); "
            "rename sequences to match the genome."
        )
    return problems


def _file_sha256(path: str) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def _gmb_version() -> str:
    from gmb import __version__

    return __version__


def resolve_consensus_settings(preset: str | None, config_path: str | None) -> dict:
    from gmb.pipeline.longread.config import load_longread_consensus_config

    if not preset and not config_path:
        raise LongReadPreparationError(
            "no long-read consensus settings: name a consensus preset or config"
        )
    cfg = load_longread_consensus_config(path=config_path, preset=preset)
    return dataclasses.asdict(cfg)


def _reusable(final_dir: str, identity: dict) -> tuple[bool, str]:
    manifest_path = os.path.join(final_dir, PREPARATION_MANIFEST)
    gtf = os.path.join(final_dir, CONSENSUS_GTF)
    if not os.path.exists(manifest_path):
        return False, "no previous consensus"
    try:
        with open(manifest_path) as fh:
            prev = json.load(fh)
    except (OSError, ValueError):
        return False, "previous preparation manifest unreadable"
    if prev.get("status") != "complete":
        return False, "previous consensus incomplete"
    for key in ("input_sha256", "consensus_settings", "gmb_version"):
        if prev.get(key) != identity[key]:
            return False, f"{key.replace('_', ' ')} changed"
    if not os.path.exists(gtf) or _file_sha256(gtf) != prev.get("consensus_gtf_sha256"):
        return False, "previous consensus file missing or modified"
    return True, "input, settings and GMB version match"


def prepare_long_read_alignments(
    alignments: str,
    genome: str,
    run_dir: str,
    consensus_preset: str | None,
    consensus_config: str | None = None,
    bin_dir: str | None = None,
    max_unstranded_fraction: float = 0.05,
    max_unknown_seqid_fraction: float = 0.05,
    log=print,
) -> dict:
    """Validate per-read alignments and produce (or reuse) consensus transcript models.

    Returns the preparation record (also written to
    ``<run_dir>/longread_consensus/longread_preparation.json``); its
    ``consensus_gtf`` is the file to pass to preflight and build as long-read
    evidence. Raises :class:`LongReadPreparationError` on any problem.
    """
    alignments = os.path.abspath(alignments)
    if not os.path.isfile(alignments) or os.path.getsize(alignments) == 0:
        raise LongReadPreparationError(f"long-read alignments missing or empty: {alignments}")
    final_dir = os.path.join(run_dir, "longread_consensus")
    logs = os.path.join(run_dir, "logs")
    os.makedirs(logs, exist_ok=True)

    settings = resolve_consensus_settings(consensus_preset, consensus_config)

    log(f"Long-read preparation: scanning {alignments}")
    t0 = time.time()
    stats = scan_alignments(alignments, _genome_seqnames(genome))
    problems = validate_alignment_stats(
        alignments, stats, max_unstranded_fraction, max_unknown_seqid_fraction
    )
    if problems:
        raise LongReadPreparationError(
            "long-read alignments rejected:\n  - " + "\n  - ".join(problems)
        )
    n_reads = stats["transcript_rows"] or None
    log(
        f"  {stats['size_bytes']:,} bytes, {stats['exon_rows']:,} exon rows"
        f"{f', {n_reads:,} reads' if n_reads else ''}, "
        f"{stats['invalid_strand_rows']:,} unstranded rows; sha256 {stats['sha256'][:12]}… "
        f"({time.time() - t0:.0f}s)"
    )

    identity = {
        "input_sha256": stats["sha256"],
        "consensus_settings": {
            "preset": consensus_preset,
            "config": os.path.abspath(consensus_config) if consensus_config else None,
            "resolved": settings,
        },
        "gmb_version": _gmb_version(),
    }
    ok, why = _reusable(final_dir, identity)
    if ok:
        with open(os.path.join(final_dir, PREPARATION_MANIFEST)) as fh:
            record = json.load(fh)
        record["action"] = "reused"
        # The verified file in THIS run directory (the directory may have been moved or
        # copied since the manifest was written).
        record["consensus_gtf"] = os.path.join(final_dir, CONSENSUS_GTF)
        log(f"Long-read consensus: REUSED ({why}) -> {record['consensus_gtf']}")
        return record

    log(
        f"Long-read consensus: RUNNING ({why}); preset={consensus_preset} "
        f"config={consensus_config}"
    )
    shutil.rmtree(final_dir, ignore_errors=True)
    os.makedirs(final_dir)
    split_dir = os.path.join(final_dir, "_split_by_seqname")
    exe = os.path.join(bin_dir, "gmb-longread-consensus") if bin_dir else "gmb-longread-consensus"
    cmd = [
        exe,
        "--input",
        alignments,
        "--output-dir",
        final_dir,
        "--split-dir",
        split_dir,
        "--no-log-file",
    ]
    if consensus_preset:
        cmd += ["--preset", consensus_preset]
    if consensus_config:
        cmd += ["--config", os.path.abspath(consensus_config)]
    log_path = os.path.join(logs, "longread_consensus.log")
    t0 = time.time()
    with open(log_path, "w") as fh:
        proc = subprocess.run(cmd, stdout=fh, stderr=subprocess.STDOUT, text=True)
    # The split files duplicate the (possibly multi-GB) input; never keep them.
    shutil.rmtree(split_dir, ignore_errors=True)
    gtf = os.path.join(final_dir, CONSENSUS_GTF)
    if proc.returncode != 0 or not os.path.exists(gtf):
        raise LongReadPreparationError(
            f"gmb-longread-consensus failed (exit {proc.returncode}); see {log_path}. "
            f"Partial output (not reusable) left in {final_dir}."
        )
    n_models = 0
    with open(gtf) as fh:
        for line in fh:
            n_models += "\ttranscript\t" in line
    if n_models == 0:
        raise LongReadPreparationError(
            f"gmb-longread-consensus produced no consensus models from {alignments}; "
            f"see {log_path} and {final_dir}/longread_consensus_dropped.tsv."
        )

    summary = {}
    try:
        with open(os.path.join(final_dir, "longread_consensus_summary.json")) as fh:
            summary = json.load(fh).get("totals", {})
    except (OSError, ValueError):
        pass
    record = {
        "status": "complete",
        "action": "ran",
        "input": alignments,
        "input_size_bytes": stats["size_bytes"],
        "input_sha256": stats["sha256"],
        "input_scan": {k: v for k, v in stats.items() if k not in ("sha256", "size_bytes")},
        "consensus_settings": identity["consensus_settings"],
        "gmb_version": identity["gmb_version"],
        "command": cmd,
        "elapsed_s": round(time.time() - t0, 1),
        "consensus_models": n_models,
        "consensus_summary": summary,
        "consensus_gtf": gtf,
        "consensus_gtf_sha256": _file_sha256(gtf),
        "log": log_path,
        "completed_at": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
    }
    # Written last: its presence marks a complete, reusable consensus.
    with open(os.path.join(final_dir, PREPARATION_MANIFEST), "w") as fh:
        json.dump(record, fh, indent=2)
    log(
        f"Long-read consensus: RAN -> {n_models:,} consensus models in "
        f"{record['elapsed_s']}s -> {record['consensus_gtf']}"
    )
    return record
