"""gmb-protein-qc: standalone post-build DIAMOND/PSAURON QC (tools mocked)."""

import csv
import json
import os
import shutil
import subprocess
from pathlib import Path
from unittest.mock import patch

import pytest

from gmb.cli import protein_qc
from gmb.pipeline.config import load_config

GMB_PKG = Path(protein_qc.__file__).resolve().parent.parent

_PSAURON_HELP = """PSAURON version 1.0.8
  -i INPUT_FASTA, --input-fasta INPUT_FASTA
  -o OUTPUT_PATH, --output-path OUTPUT_PATH
  -m MINIMUM_LENGTH, --minimum-length MINIMUM_LENGTH
  -c, --use-cpu         OPTIONAL
  -p, --protein         OPTIONAL
"""

GFF3 = """##gff-version 3
1\tGMB\tgene\t100\t900\t.\t+\t.\tID=G1
1\tGMB\tmRNA\t100\t900\t.\t+\t.\tID=G1.t1;Parent=G1;Evidence=Scallop,StringTie;tag=Ensembl_canonical
1\tGMB\tmRNA\t100\t900\t.\t+\t.\tID=G1.t2;Parent=G1;Evidence=Scallop
1\tGMB\tgene\t2000\t2600\t.\t+\t.\tID=G2
1\tGMB\tmRNA\t2000\t2600\t.\t+\t.\tID=G2.t1;Parent=G2;Evidence=Scallop;CandidateOrigin=secondary_orf;tag=Ensembl_canonical
"""
PROTEINS = {"G1.t1": "MSTAVLKQAAW", "G1.t2": "MSTAVLKQAAW", "G2.t1": "MKKLLPQ"}


def _fake_tools(calls):
    """DIAMOND hits the first submitted sequence only; PSAURON scores all."""

    def run(cmd, **kwargs):
        tool, rest = os.path.basename(cmd[0]), cmd[1:]
        out, text = "", ""
        if tool == "diamond" and rest[:1] == ["blastp"]:
            calls.append("diamond")
            with open(cmd[cmd.index("-o") + 1], "w") as fh:
                fh.write("seq_0\tSP_HIT\t55.0\t100\t420.0\t1e-50\t95.0\t90.0\n")
        elif tool == "diamond" and rest[:1] == ["version"]:
            out = "diamond version 2.2.4\n"
        elif tool == "diamond" and rest[:1] == ["dbinfo"]:
            out = "         Sequences  42\n"
        elif tool == "psauron" and "--help" in rest:
            out = _PSAURON_HELP
        elif tool == "psauron" and "-i" in rest:
            calls.append("psauron")
            ids = [ln[1:].strip() for ln in open(cmd[cmd.index("-i") + 1]) if ln.startswith(">")]
            with open(cmd[cmd.index("-o") + 1], "w") as fh:
                fh.write("psauron -i x\npsauron score: 0.5\n")
                fh.write("description,psauron_is_protein,in-frame_score\n")
                for i, sid in enumerate(ids):
                    fh.write(f"{sid},{'True' if i == 0 else 'False'},{0.95 if i == 0 else 0.2}\n")
        elif tool != "diamond":
            return subprocess.CompletedProcess(cmd, 1, "", "")
        return subprocess.CompletedProcess(cmd, 0, out, text)

    return run


@pytest.fixture
def inputs(tmp_path):
    (tmp_path / "consensus.gff3").write_text(GFF3)
    (tmp_path / "prot.fa").write_text("".join(f">{k}\n{v}\n" for k, v in PROTEINS.items()))
    (tmp_path / "db.dmnd").write_text("")
    return tmp_path


def _run(inputs, cfg=None, **kw):
    cfg = cfg or load_config(None, "fungi")
    args = type(
        "A",
        (),
        dict(
            build_dir=None,
            preset="fungi",
            config=None,
            diamond_db=str(inputs / "db.dmnd"),
            diamond_path="diamond",
            psauron_path="psauron",
            psauron_min_length=None,
            psauron_cpu=True,
        ),
    )()
    with patch("gmb.cli.protein_qc.load_config", return_value=cfg):
        config, source = protein_qc.resolve_config(args)
    calls = []
    with patch("subprocess.run", side_effect=_fake_tools(calls)):
        summary = protein_qc.run_protein_qc(
            str(inputs / "prot.fa"),
            str(inputs / "consensus.gff3"),
            config,
            str(inputs / "qc"),
            source,
            **kw,
        )
    with open(inputs / "qc" / protein_qc.QC_TSV) as fh:
        rows = {r["transcript_id"]: r for r in csv.DictReader(fh, delimiter="\t")}
    return summary, rows, calls


def test_runs_with_validation_disabled_in_config_without_mutating_it(inputs):
    cfg = load_config(None, "fungi")
    assert cfg.protein_validation.enabled is False
    _summary, rows, calls = _run(inputs, cfg)
    assert cfg.protein_validation.enabled is False  # the caller's config is untouched
    assert calls == ["diamond", "psauron"]  # one run each, on unique sequences
    assert rows["G1.t1"]["validation_status"] == "ok"


def test_schema_scores_and_secondary_orfs(inputs):
    summary, rows, _ = _run(inputs)
    with open(inputs / "qc" / protein_qc.QC_TSV) as fh:
        assert fh.readline().rstrip("\n").split("\t") == protein_qc.COLUMNS
    assert set(rows) == set(PROTEINS)
    sorf = rows["G2.t1"]
    assert sorf["candidate_origin"] == "secondary_orf"
    assert float(sorf["psauron_score"]) == pytest.approx(0.2)
    assert sorf["psauron_is_protein"] == "False"
    assert sorf["validation_status"] == "no_diamond_hit"
    assert rows["G1.t1"]["diamond_hit"] == "SP_HIT"
    assert rows["G1.t1"]["psauron_is_protein"] == "True"
    # identical proteins share one result, recorded as reuse
    assert rows["G1.t2"]["protein_sha256"] == rows["G1.t1"]["protein_sha256"]
    assert rows["G1.t2"]["protein_validation_source"] == "reused"
    assert rows["G1.t1"]["canonical"] == "True" and rows["G1.t2"]["canonical"] == "False"
    assert summary["unique_proteins"] == 2 and summary["genes"] == 2
    assert summary["groups"]["origin:secondary_orf"]["transcripts"] == 1
    assert summary["groups"]["canonical"]["transcripts"] == 2
    assert summary["groups"]["all"]["psauron_is_protein"] == {"true": 2, "false": 1, "missing": 0}
    assert summary["diamond"]["version"] == "2.2.4"
    assert summary["diamond"]["database"]["sequences"] == "42"
    assert summary["inputs"]["proteins"]["sha256"]
    on_disk = json.loads((inputs / "qc" / protein_qc.QC_SUMMARY).read_text())
    assert on_disk["validation_status"] == summary["validation_status"]


def test_refuses_to_overwrite_and_unpaired_inputs(inputs):
    _run(inputs)
    with pytest.raises(SystemExit):
        _run(inputs)
    _run(inputs, overwrite=True)
    (inputs / "prot.fa").write_text(">OTHER.t1\nMAAA\n")
    with pytest.raises(SystemExit):
        _run(inputs, overwrite=True)


def test_missing_database_fails_before_running(inputs):
    os.remove(inputs / "db.dmnd")
    with pytest.raises(SystemExit):
        _run(inputs)


def test_no_pipeline_module_reads_protein_qc_output():
    """gmb-build / gmb-finalise / canonical selection must never consume QC files."""
    readers = [
        p.relative_to(GMB_PKG).as_posix()
        for p in GMB_PKG.rglob("*.py")
        if "protein_qc" in p.read_text() and p.name != "protein_qc.py"
    ]
    assert readers == []
    assert protein_qc.QC_TSV != "protein_validation.tsv"


@pytest.mark.skipif(
    not (os.environ.get("GMB_TEST_DIAMOND") and os.environ.get("GMB_TEST_PSAURON")),
    reason="set GMB_TEST_DIAMOND and GMB_TEST_PSAURON to run with real tools",
)
def test_real_tools_small_dataset(inputs):
    diamond, psauron = os.environ["GMB_TEST_DIAMOND"], os.environ["GMB_TEST_PSAURON"]
    ref = inputs / "ref.fa"
    ref.write_text(">REF1\nMSTAVLKQAAWLLPRRGGEEFFKKLLMMNNPPQQ\n")
    subprocess.run(
        [diamond, "makedb", "--in", str(ref), "-d", str(inputs / "ref")],
        check=True,
        capture_output=True,
    )
    cfg, _ = protein_qc.resolve_config(
        type(
            "A",
            (),
            dict(
                build_dir=None,
                preset="fungi",
                config=None,
                diamond_db=str(inputs / "ref.dmnd"),
                diamond_path=diamond,
                psauron_path=psauron,
                psauron_min_length=None,
                psauron_cpu=True,
            ),
        )()
    )
    s = protein_qc.run_protein_qc(
        str(inputs / "prot.fa"), str(inputs / "consensus.gff3"), cfg, str(inputs / "qc")
    )
    assert s["groups"]["all"]["psauron_is_protein"]["missing"] == 0
    assert shutil.which(diamond) or os.path.exists(diamond)
