"""Integration tests: ``scoring.evidence_score.secondary_orfs`` on the bundled region fixture."""

import csv
import filecmp
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

_INPUTS = [
    "--scallop",
    "scallop_geneset.gtf",
    "--stringtie",
    "stringtie_geneset.gtf",
    "--helixer",
    "helixer_remapped.gff3",
    "--orthodb",
    "orthodb_geneset.gtf",
    "--uniprot",
    "uniprot_geneset.gtf",
    "--genome",
    "genome.fa",
]
_OUT_FILES = [
    "consensus.gff3",
    "cdna.fa",
    "cds.fa",
    "prot.fa",
    "evidence_attribution.tsv",
    "evidence_selection.tsv",
    "resolved_config.yaml",
]


def _inputs():
    return [a if a.startswith("--") else str(_R / a) for a in _INPUTS]


def _yaml(tmp, name, text):
    p = Path(tmp) / name
    p.write_text(text)
    return str(p)


def _build(out, *extra, ok=True):
    r = subprocess.run(
        [
            sys.executable,
            "-m",
            "gmb.cli.build",
            *_inputs(),
            "--preset",
            "fungi_latest_settings",
            "--gene-prefix",
            "ZTTEST",
            "--no-log-file",
            "--output-dir",
            str(out),
            *extra,
        ],
        capture_output=True,
        text=True,
        cwd=str(_GMB_DIR),
    )
    if ok:
        assert r.returncode == 0, r.stdout[-3000:] + r.stderr[-3000:]
    return r


@pytest.fixture(scope="module")
def runs(tmp_path_factory):
    if not _R.exists() or not _J.exists():
        pytest.skip("region fixtures missing")
    base = tmp_path_factory.mktemp("sorf_runs")
    on = _yaml(
        base,
        "on.yaml",
        "scoring:\n  primary_selection: evidence_score\n"
        "  evidence_score:\n    secondary_orfs: true\n",
    )
    off = _yaml(
        base,
        "off.yaml",
        "scoring:\n  primary_selection: evidence_score\n"
        "  evidence_score:\n    secondary_orfs: false\n",
    )
    p4d = _yaml(base, "p4d.yaml", "scoring:\n  primary_selection: evidence_score\n")
    out = {"base": base}
    for name, cfg in (("on_a", on), ("on_b", on), ("off", off), ("p4d", p4d)):
        _build(base / name, "--config", cfg, "--star-junctions", str(_J))
        out[name] = base / name
    out["on_cfg"] = on
    return out


def _rows(path):
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def test_attribution_written(runs):
    d = runs["on_a"]
    rows = _rows(d / "secondary_orfs.tsv")
    summary = json.loads((d / "secondary_orfs_summary.json").read_text())
    assert summary["discovery"]["unique_structures"] == len(rows)
    assert set(summary["structures_by_status"]) <= {
        "admitted",
        "not_admitted",
        "no_shared_coding_locus",
    }
    for r in rows:
        assert r["parents"] and int(r["n_parents"]) == len(r["parents"].split(","))
        if r["status"] == "admitted":
            assert r["admitted_by"] in ("protein", "junctions", "protein+junctions")
        if r["selected"] == "True":
            assert r["selected_gene"].startswith("ZTTEST_") and int(r["eligible_genes"]) >= 1
    assert "target species" in " ".join(summary["input_requirements"])


def test_generated_transcripts_are_valid_cds_only_models(runs):
    d = runs["on_a"]
    gff = (d / "consensus.gff3").read_text().splitlines()
    gen = {
        line.split("ID=")[1].split(";")[0]
        for line in gff
        if "\tmRNA\t" in line and "CandidateOrigin=secondary_orf" in line
    }
    selected = [r for r in _rows(d / "secondary_orfs.tsv") if r["selected"] == "True"]
    assert len(gen) == len(selected)
    exons, cds = {}, {}
    for line in gff:
        p = line.split("\t")
        if len(p) == 9 and p[2] in ("exon", "CDS"):
            parent = p[8].split("Parent=")[1].split(";")[0]
            if parent in gen:
                (exons if p[2] == "exon" else cds).setdefault(parent, []).append((p[3], p[4]))
        if len(p) == 9 and p[2] in ("five_prime_UTR", "three_prime_UTR"):
            assert p[8].split("Parent=")[1].split(";")[0] not in gen
    for t in gen:
        assert sorted(exons[t]) == sorted(cds[t])
    prot = (d / "prot.fa").read_text().split(">")[1:]
    seqs = {r.split("\n", 1)[0]: r.split("\n", 1)[1].replace("\n", "") for r in prot}
    for t in gen:
        assert seqs[t].startswith("M") and "*" not in seqs[t].rstrip("*")


def test_deterministic(runs):
    for f in _OUT_FILES + ["secondary_orfs.tsv"]:
        assert filecmp.cmp(runs["on_a"] / f, runs["on_b"] / f, shallow=False), f


def test_switched_off_is_identical_to_plain_evidence_score(runs):
    for f in _OUT_FILES + ["evidence_selection_summary.json"]:
        a, b = runs["off"] / f, runs["p4d"] / f
        if f.endswith(".json"):
            ja, jb = json.loads(a.read_text()), json.loads(b.read_text())
            ja.pop("timing"), jb.pop("timing")
            assert ja == jb
        else:
            assert filecmp.cmp(a, b, shallow=False), f
    assert not (runs["off"] / "secondary_orfs.tsv").exists()
    assert "secondary_orfs" not in (runs["off"] / "resolved_config.yaml").read_text()
    assert "secondary_orfs: true" in (runs["on_a"] / "resolved_config.yaml").read_text()


def test_gene_grouping_unchanged(runs):
    def genes(d):
        return sum(
            1 for line in (d / "consensus.gff3").read_text().splitlines() if "\tgene\t" in line
        )

    assert genes(runs["on_a"]) == genes(runs["p4d"])


def test_refused_with_another_policy(runs):
    cfg = _yaml(
        runs["base"], "bad.yaml", "scoring:\n  evidence_score:\n    secondary_orfs: true\n"
    )
    r = _build(runs["base"] / "bad", "--config", cfg, ok=False)
    assert r.returncode != 0
    assert "requires scoring.primary_selection: evidence_score" in r.stdout + r.stderr


def test_refused_without_junctions(runs):
    r = _build(runs["base"] / "nojx", "--config", runs["on_cfg"], ok=False)
    assert r.returncode != 0 and "--star-junctions" in r.stdout + r.stderr
