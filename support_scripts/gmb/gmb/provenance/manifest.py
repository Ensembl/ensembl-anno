#!/usr/bin/env python3
"""Build and write the GMB run manifest.

A run should be reproducible from ``run_manifest.json`` + the input files + the
recorded software version. Three fields earn their place from experience:

``git_dirty``
    A commit id alone does not identify what ran if the working tree was dirty.
``input sha256``
    Filenames get reused. The hash is what proves two runs saw the same bytes.
``resolved_policy``
    Some policies are decided at run time from the evidence (the
    ``backbone_intron_rescue`` applicability gate). Recording the *configured*
    mode without the *resolved* outcome would hide what actually happened.
"""

from __future__ import annotations

import hashlib
import json
import os
import platform
import shutil
import socket
import subprocess
import sys
import time
from dataclasses import asdict, dataclass, field


def _sha256(path: str | None) -> str | None:
    if not path or not os.path.exists(path) or os.path.isdir(path):
        return None
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(1 << 22), b""):
            h.update(block)
    return h.hexdigest()


def _gmb_version() -> str:
    try:
        from importlib.metadata import version

        return version("gene-model-builder")
    except Exception:
        return "unknown"


def _git_state(package_dir: str) -> dict:
    """Commit and dirtiness of the checkout GMB is running from, if it is one."""
    out = {"commit": None, "branch": None, "dirty": None, "dirty_files": None}
    if not shutil.which("git"):
        return out

    def git(*args):
        try:
            r = subprocess.run(["git", "-C", package_dir, *args],
                               capture_output=True, text=True, timeout=10)
            return r.stdout.strip() if r.returncode == 0 else None
        except Exception:
            return None

    if git("rev-parse", "--is-inside-work-tree") != "true":
        return out
    out["commit"] = git("rev-parse", "HEAD")
    out["branch"] = git("rev-parse", "--abbrev-ref", "HEAD")
    status = git("status", "--short")
    if status is not None:
        out["dirty"] = bool(status)
        out["dirty_files"] = [ln.strip() for ln in status.splitlines()][:50] or None
    return out


def _tool_version(path: str | None, flag: str = "--version") -> str | None:
    if not path:
        return None
    exe = shutil.which(path) or (path if os.path.exists(path) else None)
    if not exe:
        return None
    try:
        r = subprocess.run([exe, flag], capture_output=True, text=True, timeout=20)
        text = (r.stdout or r.stderr or "").strip().splitlines()
        return text[0][:200] if text else None
    except Exception:
        return None


def _peak_rss_kb() -> int | None:
    """Peak RSS of this process, in KB. None where the platform cannot report it."""
    try:
        import resource

        peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        # Linux reports KB; macOS/BSD report bytes.
        return int(peak) if sys.platform.startswith("linux") else int(peak / 1024)
    except Exception:
        return None


@dataclass
class RunManifest:
    """Everything needed to identify and reproduce one GMB run."""

    stage: str = "build"
    started_at: str = ""
    finished_at: str = ""
    runtime_seconds: float | None = None
    peak_rss_kb: int | None = None

    hostname: str = ""
    platform: str = ""
    cpu_count: int | None = None
    python_version: str = ""
    gmb_version: str = ""
    git: dict = field(default_factory=dict)
    external_tools: dict = field(default_factory=dict)

    command_line: str = ""
    preset: str = ""
    config_overlays: list = field(default_factory=list)
    resolved_config_path: str = ""
    resolved_config_sha256: str | None = None

    inputs: list = field(default_factory=list)
    evidence_roles: dict = field(default_factory=dict)
    evidence_weights: dict = field(default_factory=dict)
    resolved_policy: dict = field(default_factory=dict)
    preflight: dict = field(default_factory=dict)

    outputs: dict = field(default_factory=dict)
    qc: dict = field(default_factory=dict)

    def to_dict(self) -> dict:
        return asdict(self)


def build_run_manifest(
    config,
    args=None,
    inputs: dict | None = None,
    output_dir: str = "",
    started_at: float | None = None,
    rescue_decision=None,
    preflight_report=None,
    qc: dict | None = None,
    stage: str = "build",
) -> RunManifest:
    """Assemble a run manifest.

    Parameters
    ----------
    config : PipelineConfig
    args : argparse.Namespace or None
        Used for the preset and config-overlay list.
    inputs : dict or None
        ``{label: path}`` for every evidence file actually used.
    rescue_decision : RescueDecision or None
        The resolved applicability-gate outcome, so the manifest records what
        happened and not merely what was asked for.
    """
    from gmb.pipeline.applicability import normalise_rescue_mode
    from gmb.pipeline.canonical_evidence import EvidenceRoles
    from gmb.pipeline.scoring import weights_for_role

    scfg = config.scoring
    roles = EvidenceRoles.from_config(scfg)
    now = time.time()
    package_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

    m = RunManifest(stage=stage)
    m.started_at = (time.strftime("%Y-%m-%dT%H:%M:%S%z", time.localtime(started_at))
                    if started_at else "")
    m.finished_at = time.strftime("%Y-%m-%dT%H:%M:%S%z", time.localtime(now))
    m.runtime_seconds = round(now - started_at, 1) if started_at else None
    m.peak_rss_kb = _peak_rss_kb()

    m.hostname = socket.gethostname()
    m.platform = platform.platform()
    m.cpu_count = os.cpu_count()
    m.python_version = platform.python_version()
    m.gmb_version = _gmb_version()
    m.git = _git_state(package_dir)

    pv = getattr(config, "protein_validation", None)
    if pv is not None and getattr(pv, "enabled", False):
        m.external_tools = {
            "diamond": _tool_version(getattr(pv, "diamond_path", None)),
            "psauron": _tool_version(getattr(pv, "psauron_path", None)),
            "diamond_db": getattr(pv, "diamond_db", None),
            "diamond_db_sha256": None,  # large; hashed only on request
        }
    else:
        m.external_tools = {"protein_validation": "disabled"}

    m.command_line = " ".join(sys.argv)
    m.preset = getattr(args, "preset", None) or getattr(config, "preset", "") or "none"
    overlays = getattr(args, "config", None) or []
    m.config_overlays = list(overlays) if isinstance(overlays, (list, tuple)) else [overlays]
    if output_dir:
        resolved = os.path.join(output_dir, "resolved_config.yaml")
        m.resolved_config_path = resolved
        m.resolved_config_sha256 = _sha256(resolved)

    for label, path in (inputs or {}).items():
        if not path:
            continue
        entry = {
            "label": label,
            "path": os.path.abspath(path),
            "sha256": _sha256(path),
            "size_bytes": os.path.getsize(path) if os.path.exists(path) else None,
        }
        if label == "genome":
            # The genome is the coordinate system, not an evidence track: it has
            # no evidence role and no weight. Recording a role for it would
            # misreport "other" as if the genome had failed role resolution.
            entry["resolved_role"] = "genome_reference_sequence"
            entry["resolved_weight"] = None
        else:
            role = roles.role_of(label)
            entry["resolved_role"] = role
            entry["resolved_weight"] = weights_for_role(scfg.weights, role)
            m.evidence_roles[label] = role
            m.evidence_weights[label] = weights_for_role(scfg.weights, role)
        m.inputs.append(entry)

    m.resolved_policy = {
        "structural_corroboration": bool(getattr(scfg, "structural_corroboration", False)),
        "protein_support_mode": getattr(scfg, "protein_support_mode", "positional"),
        "longread_structural_guard": bool(getattr(scfg, "longread_structural_guard", False)),
        "longread_disposition": getattr(scfg, "longread_disposition", "primary_structural"),
        "backbone_intron_rescue_mode": normalise_rescue_mode(
            getattr(scfg, "backbone_intron_rescue", "off")),
    }
    if rescue_decision is not None:
        m.resolved_policy["backbone_intron_rescue"] = rescue_decision.as_dict()
    if preflight_report is not None:
        m.preflight = {
            "verdict": preflight_report.verdict,
            "counts": preflight_report.counts(),
        }
    if qc:
        m.qc = qc

    for fname in ("consensus.gff3", "cdna.fa", "cds.fa", "prot.fa",
                  "evidence_attribution.tsv", "protein_validation.tsv"):
        fpath = os.path.join(output_dir, fname) if output_dir else None
        if fpath and os.path.exists(fpath):
            m.outputs[fname] = {"size_bytes": os.path.getsize(fpath),
                                "sha256": _sha256(fpath)}
    return m


def write_run_manifest(manifest: RunManifest, output_dir: str) -> tuple[str, str]:
    """Write run_manifest.json and a flattened run_manifest.tsv."""
    os.makedirs(output_dir, exist_ok=True)
    data = manifest.to_dict()

    json_path = os.path.join(output_dir, "run_manifest.json")
    with open(json_path, "w") as fh:
        json.dump(data, fh, indent=2, default=str)

    def flatten(obj, prefix=""):
        if isinstance(obj, dict):
            for k, v in obj.items():
                yield from flatten(v, f"{prefix}{k}." if prefix or True else k)
        elif isinstance(obj, list):
            for i, v in enumerate(obj):
                yield from flatten(v, f"{prefix}{i}.")
        else:
            yield prefix.rstrip("."), obj

    tsv_path = os.path.join(output_dir, "run_manifest.tsv")
    with open(tsv_path, "w") as fh:
        fh.write("field\tvalue\n")
        for key, value in flatten(data):
            text = "" if value is None else str(value).replace("\t", " ").replace("\n", " ")
            fh.write(f"{key}\t{text}\n")
    return json_path, tsv_path
