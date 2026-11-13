"""Integration tests: ``primary_selection: evidence_score`` on the bundled region fixture."""

import csv
import json
import subprocess
import sys
from pathlib import Path

import pytest

_HERE = Path(__file__).parent
_GMB_DIR = _HERE.parent
_R = _HERE / "fixtures" / "z_tritici_region1"
_J = _HERE / "fixtures" / "star_junctions_region1" / "manifest.tsv"

pytestmark = pytest.mark.integration

_INPUTS = ["--scallop", "scallop_geneset.gtf", "--stringtie", "stringtie_geneset.gtf",
           "--helixer", "helixer_remapped.gff3", "--orthodb", "orthodb_geneset.gtf",
           "--uniprot", "uniprot_geneset.gtf", "--genome", "genome.fa"]


def _inputs():
    return [a if a.startswith("--") else str(_R / a) for a in _INPUTS]


def _overlay(tmp):
    p = Path(tmp) / "evidence_score.yaml"
    p.write_text("scoring:\n  primary_selection: evidence_score\n")
    return str(p)


def _run(module, *args):
    return subprocess.run([sys.executable, "-m", module, *args], capture_output=True, text=True,
                          cwd=str(_GMB_DIR))


def _build(out, *extra, ok=True):
    r = _run("gmb.cli.build", *_inputs(), "--preset", "fungi_latest_settings", "--gene-prefix",
             "ZTTEST", "--no-log-file", "--output-dir", str(out), *extra)
    if ok:
        assert r.returncode == 0, r.stdout[-3000:] + r.stderr[-3000:]
    return r


@pytest.fixture(scope="module")
def runs(tmp_path_factory):
    if not _R.exists() or not _J.exists():
        pytest.skip("region fixtures missing")
    base = tmp_path_factory.mktemp("evidence_runs")
    ov = _overlay(base)
    out = {}
    for name in ("a", "b"):
        _build(base / name, "--config", ov, "--star-junctions", str(_J))
        out[name] = base / name
    _build(base / "legacy", "--star-junctions", str(_J))
    out["legacy"] = base / "legacy"
    out["overlay"] = ov
    out["base"] = base
    return out


def _rows(path):
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def test_decision_table_and_provenance(runs):
    rows = _rows(runs["a"] / "evidence_selection.tsv")
    assert rows
    assert {r["outcome"] for r in rows} <= {"replaced", "protected_backbone", "kept_below_margin",
                                            "kept_no_eligible", "no_cds_primary"}
    for r in rows:
        if r["outcome"] == "replaced":
            assert float(r["delta_total"]) >= 0.5 and float(r["delta_evidence"]) >= 0.5
            assert r["new_candidate"] and r["new_structure"] != r["old_structure"]
    s = json.loads((runs["a"] / "evidence_selection_summary.json").read_text())
    assert s["junction_roles"] == ["construction"] and s["policy_version"] == "p4d-1"
    assert {lib["library_id"] for lib in s["junction_libraries"]} == {
        "SRR17893325", "SRR22570518", "SRR18918656"}
    assert all(len(lib["sha256"]) == 64 for lib in s["junction_libraries"])
    assert "evidence_score" in (runs["a"] / "resolved_config.yaml").read_text()
    assert "evidence_score" not in (runs["legacy"] / "resolved_config.yaml").read_text()
    assert not (runs["legacy"] / "evidence_selection.tsv").exists()


def test_deterministic(runs):
    for name in ("consensus.gff3", "prot.fa", "evidence_selection.tsv", "evidence_attribution.tsv"):
        assert (runs["a"] / name).read_bytes() == (runs["b"] / name).read_bytes(), name


def test_every_gene_keeps_a_valid_primary(runs):
    genes, t1 = set(), set()
    for line in (runs["a"] / "consensus.gff3").read_text().splitlines():
        p = line.split("\t")
        if len(p) == 9 and p[2] == "mRNA":
            a = dict(kv.split("=", 1) for kv in p[8].split(";") if "=" in kv)
            genes.add(a["Parent"])
            if a["ID"] == a["Parent"] + ".t1":
                t1.add(a["Parent"])
    assert genes == t1
    proteins = (runs["a"] / "prot.fa").read_text().split(">")[1:]
    assert all("*" not in "".join(p.splitlines()[1:]).rstrip("*") for p in proteins)


def test_phase1_diagnostics_still_written(runs):
    d = runs["a"] / "junction_evidence"
    assert (d / "selected_models.tsv").exists()
    notes = json.loads((d / "junction_evidence_summary.json").read_text())["notes"]
    assert any("evidence_score" in n for n in notes)


def test_refused_without_junctions(runs):
    r = _build(runs["base"] / "nojunctions", "--config", runs["overlay"], ok=False)
    assert r.returncode != 0 and "requires STAR junction evidence" in r.stdout + r.stderr


@pytest.mark.parametrize("roles", ["holdout", "construction,holdout", "all"])
def test_refused_with_non_construction_roles(runs, roles):
    r = _build(runs["base"] / f"roles_{roles.replace(',', '_')}", "--config", runs["overlay"],
               "--star-junctions", str(_J), "--star-junction-roles", roles, ok=False)
    assert r.returncode != 0 and "construction" in r.stdout + r.stderr


def test_preflight_fails_without_junctions_and_passes_with(runs):
    common = [*_inputs(), "--preset", "fungi_latest_settings", "--config", runs["overlay"],
              "--no-log-file", "--json"]
    bad = _run("gmb.cli.preflight", *common)
    assert bad.returncode == 1
    good = _run("gmb.cli.preflight", *common, "--star-junctions", str(_J))
    checks = {c["name"]: c["verdict"] for c in json.loads(good.stdout)["checks"]}
    assert checks["evidence_score_inputs"] == "PASS"
