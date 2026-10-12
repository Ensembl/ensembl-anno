#!/usr/bin/env python3
"""Stable Python entry point for orchestration pipelines.

``run_gene_model_builder()`` is the **only** supported way to drive GMB from
Python. Everything under ``gmb.pipeline`` is internal and may change without
notice; an orchestration layer that imports ``gmb.pipeline.scoring`` or
``gmb.pipeline.builder`` has taken a dependency it is not entitled to.

The function runs the same stages as the CLI, in the same order, with the same
configuration layering::

    [long-read consensus]  ->  preflight  ->  build  ->  finalise

The long-read consensus stage runs only when per-read long-read alignments are
supplied (``long_read_alignments``); it calls ``gmb-longread-consensus`` and the
resulting consensus models are the long-read evidence that preflight validates
and the build uses. ``gmb-run`` is a thin command-line wrapper around this
function, so the command and the API apply the same policy.

Example
-------
::

    from gmb import run_gene_model_builder

    result = run_gene_model_builder(
        genome="/data/genome.fa",
        backbone="/data/backbone.gff3",
        backbone_kind="helixer",
        short_read=["/data/scallop.gtf", "/data/stringtie.gtf"],
        protein_alignment=["/data/orthodb.gtf"],
        long_read_alignments="/data/minimap2_reads.gtf",   # optional, per-read
        preset="fungi_latest_settings",
        output_dir="/work/gmb",
        gene_prefix="XXGMB",
    )
    if result.ok:
        ship(result.handover_dir)
"""

from __future__ import annotations

import json
import os
import shutil
import subprocess
import sys
from dataclasses import dataclass, field

__all__ = ["GmbResult", "run_gene_model_builder"]

# CLI flag to use for each short-read / protein track, in order. GMB ships one
# flag per historically-used tool; the flag is only a slot -- the evidence ROLE
# is what determines behaviour (see docs/input_contract.md).
_SHORT_READ_FLAGS = ("--scallop", "--stringtie")
_PROTEIN_FLAGS = ("--orthodb", "--uniprot", "--genblast")
_LONG_READ_FLAGS = ("--minimap2",)


@dataclass
class GmbResult:
    """Outcome of one GMB run."""

    output_dir: str
    preflight_verdict: str | None = None
    preflight_report: dict | None = None
    qc_passed: bool | None = None
    qc_report: dict | None = None
    build_dir: str = ""
    handover_dir: str = ""
    annotation: str | None = None
    annotation_all_isoforms: str | None = None
    proteins: str | None = None
    cds: str | None = None
    cdna: str | None = None
    evidence_attribution: str | None = None
    resolved_config: str | None = None
    run_manifest: str | None = None
    handover_manifest: str | None = None
    stages_run: list = field(default_factory=list)
    returncodes: dict = field(default_factory=dict)
    long_read_models: str | None = None
    longread_preparation: dict | None = None
    failed_stage: str | None = None
    failure_reason: str | None = None

    @property
    def ok(self) -> bool:
        """True when every stage that ran succeeded, QC passed and finalisation ran."""
        if self.failed_stage is not None:
            return False
        if self.preflight_verdict == "FAIL":
            return False
        if self.qc_passed is False:
            return False
        if "finalise" not in self.stages_run:
            return False
        return all(rc == 0 for rc in self.returncodes.values())

    def to_dict(self) -> dict:
        d = dict(self.__dict__)
        d.pop("preflight_report", None)
        d.pop("qc_report", None)
        d["ok"] = self.ok
        return d


def _flags_for(tracks, flags, kind: str) -> list:
    tracks = [t for t in (tracks or []) if t]
    if len(tracks) > len(flags):
        raise ValueError(
            f"{len(tracks)} {kind} track(s) supplied but GMB exposes only "
            f"{len(flags)} {kind} input slot(s) ({', '.join(flags)}). Merge the "
            f"tracks upstream, or add a slot."
        )
    out = []
    for flag, path in zip(flags, tracks):
        out += [flag, str(path)]
    return out


def _run(cmd: list, log_path: str) -> int:
    os.makedirs(os.path.dirname(log_path), exist_ok=True)
    with open(log_path, "w") as log:
        proc = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                              text=True)
        log.write(proc.stdout or "")
    return proc.returncode


def _read_json(path: str):
    try:
        with open(path) as fh:
            return json.load(fh)
    except Exception:
        return None


def run_gene_model_builder(
    genome: str,
    backbone: str | None = None,
    backbone_kind: str = "helixer",
    backbone_label: str | None = None,
    short_read=None,
    long_read=None,
    protein_alignment=None,
    preset: str = "standard",
    config=None,
    output_dir: str = "gmb_run",
    gene_prefix: str = "GMB",
    run_preflight: bool = True,
    stop_on_preflight_fail: bool = True,
    validate_fasta: bool = True,
    bin_dir: str | None = None,
    long_read_alignments: str | None = None,
    longread_consensus_preset: str | None = None,
    longread_consensus_config: str | None = None,
    log=print,
) -> GmbResult:
    """Run preflight, build and finalise, and return the handover paths.

    Parameters
    ----------
    genome : str
        Genome FASTA. Defines the sequence names every other input must use.
    backbone : str or None
        Ab initio backbone annotation. One per run.
    backbone_kind : {"helixer", "tiberius", "backbone"}
        Which CLI slot the backbone occupies. This also sets
        ``scoring.backbone_label``, which is how the backbone role resolves.
        ``"backbone"`` is the generic slot for any other predictor; its label
        comes from the annotation's own source column unless *backbone_label*
        is given.
    backbone_label : str or None
        Override the source label of the backbone track.
    short_read, long_read, protein_alignment : sequence of str or None
        Evidence tracks. ``long_read`` takes already-collapsed long-read
        transcript models and uses them as they are. A track that looks like
        per-read alignments (preflight ``longread_collapsed``) stops the run when
        preflight runs. Long reads are optional: with none supplied there is no
        special case and the long-read guard cannot fire.
    long_read_alignments : str or None
        Per-read long-read alignments (GTF, one transcript per read, strand set
        upstream). Collapsed into consensus models by ``gmb-longread-consensus``
        in ``<output_dir>/longread_consensus/`` before preflight; a matching
        earlier consensus (same input SHA-256, settings and GMB version) is
        reused. Mutually exclusive with ``long_read``.
    longread_consensus_preset, longread_consensus_config : str or None
        Consensus settings. When both are None the preset follows the build
        preset (``fungi`` and ``fungi_latest_settings`` -> ``fungi``); any other
        build preset must name one.
    preset : str
        ``"standard"`` (neutral), ``"apicomplexa"``, ``"fungi"`` (historical
        fungal baseline) or ``"fungi_latest_settings"`` (recommended for fungi).
    config : sequence of str or None
        YAML overlays applied after the preset, in order; last wins.
    output_dir : str
        Root. ``longread_consensus/`` (when needed), ``preflight/``, ``build/``,
        ``finalise/`` (the handover), ``logs/`` and ``run_summary.json`` are
        created beneath it.
    stop_on_preflight_fail : bool
        When True (default) a preflight FAIL stops the run before the expensive
        build. Set False to record the verdict and build anyway.
    validate_fasta : bool
        Run sequence QC at the end of the build. A QC failure makes
        ``result.qc_passed`` False and ``result.ok`` False; the outputs are still
        written, and so is the run manifest.
    bin_dir : str or None
        Directory holding the ``gmb-*`` executables. Defaults to the directory of
        the running interpreter, which is correct inside a virtualenv or conda env.
    log : callable
        Receives one progress line at a time (default ``print``).

    Returns
    -------
    GmbResult
    """
    if backbone_kind not in ("helixer", "tiberius", "backbone"):
        raise ValueError(
            "backbone_kind must be 'helixer', 'tiberius' or 'backbone'")
    if long_read_alignments and long_read:
        raise ValueError(
            "supply long-read evidence once: per-read alignments (long_read_alignments, "
            "collapsed automatically) or prepared transcript models (long_read), not both")

    bin_dir = bin_dir or os.path.dirname(sys.executable)
    out = os.path.abspath(output_dir)
    logs = os.path.join(out, "logs")
    os.makedirs(logs, exist_ok=True)

    overlays = [str(c) for c in (config or [])]
    config_flags: list = []
    for c in overlays:
        config_flags += ["--config", c]

    result = GmbResult(output_dir=out,
                       build_dir=os.path.join(out, "build"),
                       handover_dir=os.path.join(out, "finalise"))

    # Stage outputs from an earlier run in the same directory are removed, so no
    # stale file can sit beside new results. longread_consensus/ is kept: it is
    # reused only when its recorded input and settings match (longread_prep).
    stale = [d for d in ("preflight", "build", "finalise") if os.path.isdir(os.path.join(out, d))]
    if stale:
        log(f"Removing previous {', '.join(stale)} output(s) in {out}")
        for d in stale:
            shutil.rmtree(os.path.join(out, d))

    def finish(stage=None, reason=None):
        if stage:
            result.failed_stage, result.failure_reason = stage, reason
            log(f"FAILED at {stage}: {reason}")
        with open(os.path.join(out, "run_summary.json"), "w") as fh:
            json.dump(result.to_dict(), fh, indent=2, default=str)
        return result

    # ---- long-read consensus (only for per-read alignments) ---------------
    long_read_tracks = list(long_read or [])
    if long_read_alignments:
        from gmb.pipeline.config import load_config
        from gmb.pipeline.longread_prep import (
            CONSENSUS_PRESET_FOR_BUILD_PRESET,
            LongReadPreparationError,
            prepare_long_read_alignments,
        )

        lr_preset = longread_consensus_preset
        if not lr_preset and not longread_consensus_config:
            lr_preset = CONSENSUS_PRESET_FOR_BUILD_PRESET.get(preset)
            if lr_preset is None:
                raise ValueError(
                    f"no long-read consensus preset is defined for build preset "
                    f"{preset!r}; pass longread_consensus_preset or longread_consensus_config")
        pf_cfg = load_config(overlays or None, preset).preflight
        result.stages_run.append("longread_consensus")
        try:
            record = prepare_long_read_alignments(
                long_read_alignments, genome, out, lr_preset, longread_consensus_config,
                bin_dir=bin_dir,
                max_unstranded_fraction=pf_cfg.max_unstranded_fraction,
                max_unknown_seqid_fraction=pf_cfg.max_unknown_seqid_fraction,
                log=log)
        except LongReadPreparationError as exc:
            result.returncodes["longread_consensus"] = 1
            return finish("longread_consensus", str(exc))
        result.returncodes["longread_consensus"] = 0
        result.longread_preparation = record
        long_read_tracks = [record["consensus_gtf"]]
    else:
        result.longread_preparation = {"action": "skipped" if not long_read else
                                       "not_needed (prepared models supplied)"}
        log("Long-read consensus: SKIPPED ("
            + ("prepared long-read models supplied" if long_read else
               "no long-read alignments supplied") + ")")
    result.long_read_models = long_read_tracks[0] if long_read_tracks else None

    evidence: list = []
    if backbone:
        evidence += [f"--{backbone_kind}", str(backbone)]
        if backbone_label:
            evidence += ["--backbone-label", str(backbone_label)]
    evidence += _flags_for(short_read, _SHORT_READ_FLAGS, "short-read")
    evidence += _flags_for(long_read_tracks, _LONG_READ_FLAGS, "long-read")
    evidence += _flags_for(protein_alignment, _PROTEIN_FLAGS, "protein-alignment")

    # ---- preflight --------------------------------------------------------
    if run_preflight:
        pf_dir = os.path.join(out, "preflight")
        cmd = ([os.path.join(bin_dir, "gmb-preflight"), "--preset", preset]
               + config_flags + ["--genome", str(genome)] + evidence
               + ["--output-dir", pf_dir, "--allow-fail"])
        log("Preflight: validating the prepared inputs")
        rc = _run(cmd, os.path.join(logs, "preflight.log"))
        result.returncodes["preflight"] = rc
        result.stages_run.append("preflight")
        report = _read_json(os.path.join(pf_dir, "preflight_report.json"))
        if report:
            result.preflight_report = report
            result.preflight_verdict = report.get("verdict")
        log(f"  verdict: {result.preflight_verdict} (report: {pf_dir}/preflight_report.txt)")
        per_read = [c for c in (report or {}).get("checks", [])
                    if c.get("name") == "longread_collapsed" and c.get("verdict") != "PASS"]
        if per_read:
            return finish("preflight", per_read[0].get("message", "") + " Per-read long-read "
                          "alignments must be supplied as long_read_alignments, which are "
                          "collapsed automatically.")
        if result.preflight_verdict == "FAIL" and stop_on_preflight_fail:
            fails = [c for c in report.get("checks", []) if c.get("verdict") == "FAIL"]
            return finish("preflight", "; ".join(
                f"{c.get('name')} [{c.get('target') or ''}]: {c.get('message')}" for c in fails)
                or "preflight FAIL")
        if report is None and rc != 0:
            return finish("preflight", f"gmb-preflight exited {rc} without a report; "
                          f"see {logs}/preflight.log")

    # ---- build ------------------------------------------------------------
    cmd = ([os.path.join(bin_dir, "gmb-build"), "--preset", preset]
           + config_flags + ["--genome", str(genome)] + evidence
           + ["--gene-prefix", gene_prefix, "--output-dir", result.build_dir])
    if validate_fasta:
        cmd.append("--validate-fasta")
    log("Build: integrating evidence and selecting gene models")
    rc = _run(cmd, os.path.join(logs, "build.log"))
    result.returncodes["build"] = rc
    result.stages_run.append("build")

    qc = _read_json(os.path.join(result.build_dir, "fasta_qc_report.json"))
    if qc is not None:
        result.qc_report = qc
        result.qc_passed = bool(qc.get("pass"))

    for attr, rel in (("evidence_attribution", "evidence_attribution.tsv"),
                      ("resolved_config", "resolved_config.yaml"),
                      ("run_manifest", "run_manifest.json")):
        path = os.path.join(result.build_dir, rel)
        if os.path.exists(path):
            setattr(result, attr, path)

    # A build that failed QC has still written a complete annotation, but it must
    # not be finalised into a handover -- see docs/qc.md.
    if result.qc_passed is False:
        return finish("sequence_qc", f"FASTA QC failed: {qc.get('failed_checks')}; "
                      f"see {result.build_dir}/fasta_qc_report.json")
    if rc != 0:
        return finish("build", f"gmb-build exited {rc}; see {logs}/build.log")

    # ---- finalise ---------------------------------------------------------
    # No --preset/--config: gmb-finalise runs under the build's own
    # resolved_config.yaml, so finalisation cannot use settings that never
    # applied to the build.
    cmd = [os.path.join(bin_dir, "gmb-finalise"),
           "--build-dir", result.build_dir, "--genome", str(genome),
           "--output-dir", result.handover_dir]
    log("Finalise: regenerating sequences, QC, canonical transcripts, handover manifest")
    rc = _run(cmd, os.path.join(logs, "finalise.log"))
    result.returncodes["finalise"] = rc
    result.stages_run.append("finalise")

    for attr, rel in (
        ("annotation", os.path.join("canonical", "consensus.canonical_annotated.gff3")),
        ("annotation_all_isoforms", "consensus.gff3"),
        ("proteins", "prot.fa"),
        ("cds", "cds.fa"),
        ("cdna", "cdna.fa"),
        ("handover_manifest", "handover_manifest.json"),
    ):
        path = os.path.join(result.handover_dir, rel)
        if os.path.exists(path):
            setattr(result, attr, path)
    if rc != 0:
        return finish("finalise", f"gmb-finalise exited {rc}; see {logs}/finalise.log")
    return finish()
