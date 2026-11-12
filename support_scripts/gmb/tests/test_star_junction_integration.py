"""Integration tests: STAR junction diagnostics on the bundled Z. tritici region fixture.

The junction fixture (``fixtures/star_junctions_region1/``) holds real STAR junctions for
chr 1:1-500,000 from three construction and two held-out libraries, plus the real
zero-byte file of a failed run. Builds take ~2 s each.
"""

import csv
import json
import os
import subprocess
import sys
from pathlib import Path

import pytest

_HERE = Path(__file__).parent
_GMB_DIR = _HERE.parent
_R = _HERE / "fixtures" / "z_tritici_region1"
_J = _HERE / "fixtures" / "star_junctions_region1" / "manifest.tsv"

# Outputs that must be byte-identical with and without junction diagnostics.
_ANNOTATION_OUTPUTS = [
    "consensus.gff3",
    "cdna.fa",
    "cds.fa",
    "prot.fa",
    "evidence_attribution.tsv",
    "summary.json",
    "summary.tsv",
    "resolved_config.yaml",
    "resolved_config_sha256",
    "excluded_candidates.tsv",
    "collapsed_duplicate_transcripts.tsv",
]
_JUNCTION_TABLES = [
    "libraries.tsv",
    "candidates_introns.tsv",
    "candidates_models.tsv",
    "selected_introns.tsv",
    "selected_models.tsv",
]

pytestmark = pytest.mark.integration


def _inputs():
    return [
        "--scallop",
        str(_R / "scallop_geneset.gtf"),
        "--stringtie",
        str(_R / "stringtie_geneset.gtf"),
        "--helixer",
        str(_R / "helixer_remapped.gff3"),
        "--orthodb",
        str(_R / "orthodb_geneset.gtf"),
        "--uniprot",
        str(_R / "uniprot_geneset.gtf"),
        "--genome",
        str(_R / "genome.fa"),
    ]


def _build(out: Path, preset: str, *extra) -> Path:
    cmd = [
        sys.executable,
        "-m",
        "gmb.cli.build",
        *_inputs(),
        "--preset",
        preset,
        "--gene-prefix",
        "ZTTEST",
        "--no-log-file",
        "--output-dir",
        str(out),
        *extra,
    ]
    r = subprocess.run(cmd, capture_output=True, text=True, cwd=str(_GMB_DIR))
    assert r.returncode == 0, r.stdout[-3000:] + r.stderr[-3000:]
    return out


@pytest.fixture(scope="module")
def runs(tmp_path_factory):
    if not _R.exists() or not _J.exists():
        pytest.skip("region fixtures missing")
    base = tmp_path_factory.mktemp("junction_runs")
    out = {}
    for preset in ("fungi_latest_settings", "fungi"):
        out[(preset, "none")] = _build(base / f"{preset}_none", preset)
        out[(preset, "construction")] = _build(
            base / f"{preset}_c", preset, "--star-junctions", str(_J)
        )
    out["repeat"] = _build(base / "repeat", "fungi_latest_settings", "--star-junctions", str(_J))
    out["holdout"] = _build(
        base / "holdout",
        "fungi_latest_settings",
        "--star-junctions",
        str(_J),
        "--star-junction-roles",
        "holdout",
    )
    return out


def _rows(path):
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


@pytest.mark.parametrize("preset", ["fungi_latest_settings", "fungi"])
def test_annotation_outputs_unchanged_by_diagnostics(runs, preset):
    a, b = runs[(preset, "none")], runs[(preset, "construction")]
    for name in _ANNOTATION_OUTPUTS:
        if (a / name).exists() or (b / name).exists():
            assert (a / name).read_bytes() == (b / name).read_bytes(), name
    assert not (a / "junction_evidence").exists()
    assert (b / "junction_evidence" / "selected_models.tsv").exists()


def test_diagnostics_are_deterministic(runs):
    a = runs[("fungi_latest_settings", "construction")] / "junction_evidence"
    b = runs["repeat"] / "junction_evidence"
    for name in _JUNCTION_TABLES:
        assert (a / name).read_bytes() == (b / name).read_bytes(), name


def test_holdout_excluded_from_construction_diagnostics(runs):
    libs = {
        r["library_id"]: r
        for r in _rows(
            runs[("fungi_latest_settings", "construction")] / "junction_evidence" / "libraries.tsv"
        )
    }
    assert {k for k, r in libs.items() if r["status"] == "LOADED"} == {
        "SRR17893325",
        "SRR22570518",
        "SRR18918656",
    }
    assert libs["SRR22428702"]["status"] == libs["SRR31679459"]["status"] == "NOT_SELECTED"
    assert libs["SRR6941156"]["status"] == "NOT_SELECTED"  # failed run, role 'excluded'
    held = {
        r["library_id"]: r["status"]
        for r in _rows(runs["holdout"] / "junction_evidence" / "libraries.tsv")
    }
    assert {k for k, s in held.items() if s == "LOADED"} == {"SRR22428702", "SRR31679459"}
    # A separate held-out evaluation still leaves the annotation unchanged.
    assert (runs["holdout"] / "consensus.gff3").read_bytes() == (
        runs[("fungi_latest_settings", "none")] / "consensus.gff3"
    ).read_bytes()


def test_selected_diagnostics_cover_every_cds_transcript(runs):
    d = runs[("fungi_latest_settings", "construction")]
    with open(d / "consensus.gff3") as fh:
        tx = {
            l.split("\t")[8].split("Parent=")[1].split(";")[0].strip()
            for l in fh
            if "\tCDS\t" in l
        }
    models = _rows(d / "junction_evidence" / "selected_models.tsv")
    assert {m["model_id"] for m in models} == tx
    assert all(m["extra"].startswith("candidate_transcript_id=") for m in models)
    s = json.load(open(d / "junction_evidence" / "junction_evidence_summary.json"))
    assert s["index"]["libraries_loaded"] == 3 and s["index"]["roles"] == ["construction"]
    assert any("not independent validation" in n for n in s["notes"])


def test_long_introns_never_contradicted(runs):
    rows = _rows(
        runs[("fungi_latest_settings", "construction")]
        / "junction_evidence"
        / "candidates_introns.tsv"
    )
    long_rows = [r for r in rows if int(r["length"]) > 1000]
    assert all(r["evidence_class"] in ("unobservable", "contradicted") for r in long_rows)
    assert all(r["status"] in ("UNOBSERVABLE_LONG", "OPPOSITE_STRAND") for r in long_rows)


def test_cli_matches_build_representation(runs, tmp_path):
    d = runs[("fungi_latest_settings", "construction")]
    out = tmp_path / "cli"
    cmd = [
        sys.executable,
        "-m",
        "gmb.cli.junction_evidence",
        "--star-junctions",
        str(_J),
        "--genome",
        str(_R / "genome.fa"),
        "--annotation",
        f"GMB={d / 'consensus.gff3'}",
        "--output-dir",
        str(out),
    ]
    r = subprocess.run(cmd, capture_output=True, text=True, cwd=str(_GMB_DIR))
    assert r.returncode == 0, r.stderr
    assert (out / "GMB_introns.tsv").read_bytes() == (
        d / "junction_evidence" / "selected_introns.tsv"
    ).read_bytes()


def test_preflight_reports_junction_checks(tmp_path):
    out = tmp_path / "pf"
    cmd = [
        sys.executable,
        "-m",
        "gmb.cli.preflight",
        "--preset",
        "fungi_latest_settings",
        *_inputs(),
        "--star-junctions",
        str(_J),
        "--output-dir",
        str(out),
        "--no-log-file",
        "--allow-fail",
    ]
    r = subprocess.run(cmd, capture_output=True, text=True, cwd=str(_GMB_DIR))
    assert r.returncode == 0, r.stderr
    checks = {
        c["name"]: c["verdict"] for c in json.load(open(out / "preflight_report.json"))["checks"]
    }
    assert checks["star_junction_libraries"] == "PASS"
    assert checks["star_junction_count_unit"] == "WARN"


def test_gmb_run_accepts_junction_arguments():
    from gmb.cli.run import build_parser

    a = build_parser().parse_args(
        [
            "--preset",
            "fungi_latest_settings",
            "--genome",
            "g.fa",
            "--helixer",
            "h.gff3",
            "--output-dir",
            "o",
            "--star-junctions",
            "m.tsv",
            "--star-junction-roles",
            "construction",
        ]
    )
    assert a.star_junctions == "m.tsv" and a.star_junction_roles == "construction"
    import inspect

    from gmb.api import run_gene_model_builder

    assert "star_junctions" in inspect.signature(run_gene_model_builder).parameters
