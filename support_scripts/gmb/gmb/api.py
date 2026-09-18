#!/usr/bin/env python3
"""Stable Python entry point for orchestration pipelines.

``run_gene_model_builder()`` is the **only** supported way to drive GMB from
Python. Everything under ``gmb.pipeline`` is internal and may change without
notice; an orchestration layer that imports ``gmb.pipeline.scoring`` or
``gmb.pipeline.builder`` has taken a dependency it is not entitled to.

The function runs the same three stages as the CLI, in the same order, with the
same configuration layering::

    preflight  ->  build  ->  finalise

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
        preset="fungi",
        output_dir="/work/gmb",
        gene_prefix="XXGMB",
    )
    if result.ok:
        ship(result.handover_dir)
"""

from __future__ import annotations

import json
import os
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

    @property
    def ok(self) -> bool:
        """True when every stage that ran succeeded and QC passed."""
        if self.preflight_verdict == "FAIL":
            return False
        if self.qc_passed is False:
            return False
        return all(rc == 0 for rc in self.returncodes.values())


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
) -> GmbResult:
    """Run preflight, build and finalise, and return the handover paths.

    Parameters
    ----------
    genome : str
        Genome FASTA. Defines the sequence names every other input must use.
    backbone : str or None
        Ab initio backbone annotation. One per run.
    backbone_kind : {"helixer", "tiberius"}
        Which CLI slot the backbone occupies. This also sets
        ``scoring.backbone_label``, which is how the backbone role resolves.
    short_read, long_read, protein_alignment : sequence of str or None
        Evidence tracks. ``long_read`` is genuinely optional: with none supplied
        there is no special case and the long-read guard cannot fire.
    preset : str
        ``"standard"`` (neutral), ``"apicomplexa"``, ``"fungi"``.
    config : sequence of str or None
        YAML overlays applied after the preset, in order; last wins.
    output_dir : str
        Root. ``preflight/``, ``build/``, ``finalise/`` and ``logs/`` are created
        beneath it.
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

    Returns
    -------
    GmbResult
    """
    if backbone_kind not in ("helixer", "tiberius"):
        raise ValueError("backbone_kind must be 'helixer' or 'tiberius'")

    bin_dir = bin_dir or os.path.dirname(sys.executable)
    out = os.path.abspath(output_dir)
    logs = os.path.join(out, "logs")
    os.makedirs(logs, exist_ok=True)

    overlays = [str(c) for c in (config or [])]
    config_flags: list = []
    for c in overlays:
        config_flags += ["--config", c]

    evidence: list = []
    if backbone:
        evidence += [f"--{backbone_kind}", str(backbone)]
    evidence += _flags_for(short_read, _SHORT_READ_FLAGS, "short-read")
    evidence += _flags_for(long_read, _LONG_READ_FLAGS, "long-read")
    evidence += _flags_for(protein_alignment, _PROTEIN_FLAGS, "protein-alignment")

    result = GmbResult(output_dir=out,
                       build_dir=os.path.join(out, "build"),
                       handover_dir=os.path.join(out, "finalise"))

    # ---- preflight --------------------------------------------------------
    if run_preflight:
        pf_dir = os.path.join(out, "preflight")
        cmd = ([os.path.join(bin_dir, "gmb-preflight"), "--preset", preset]
               + config_flags + ["--genome", str(genome)] + evidence
               + ["--output-dir", pf_dir, "--allow-fail"])
        rc = _run(cmd, os.path.join(logs, "preflight.log"))
        result.returncodes["preflight"] = rc
        result.stages_run.append("preflight")
        report = _read_json(os.path.join(pf_dir, "preflight_report.json"))
        if report:
            result.preflight_report = report
            result.preflight_verdict = report.get("verdict")
        if result.preflight_verdict == "FAIL" and stop_on_preflight_fail:
            return result

    # ---- build ------------------------------------------------------------
    cmd = ([os.path.join(bin_dir, "gmb-build"), "--preset", preset]
           + config_flags + ["--genome", str(genome)] + evidence
           + ["--gene-prefix", gene_prefix, "--output-dir", result.build_dir])
    if validate_fasta:
        cmd.append("--validate-fasta")
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
        return result
    if rc != 0:
        return result

    # ---- finalise ---------------------------------------------------------
    # No --preset/--config: gmb-finalise runs under the build's own
    # resolved_config.yaml, so finalisation cannot use settings that never
    # applied to the build.
    cmd = [os.path.join(bin_dir, "gmb-finalise"),
           "--build-dir", result.build_dir, "--genome", str(genome),
           "--output-dir", result.handover_dir]
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
    return result
