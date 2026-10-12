"""gmb-run and automatic long-read consensus, on the bundled Z. tritici fixture.

The fixture's minimap2_reads.gtf holds the 10,596 real per-read long-read
alignments lying wholly inside chr1:1-500,000 (Ensembl anno Minimap2 output,
strand from the ts tag).
"""

import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

from gmb.pipeline.config import dump_config, list_build_presets, load_config

GMB_DIR = Path(__file__).resolve().parent.parent
FIXTURE = GMB_DIR / "tests" / "fixtures" / "z_tritici_region1"
READS = FIXTURE / "minimap2_reads.gtf"
BIN = Path(sys.executable).parent

EVIDENCE = [
    "--genome",
    str(FIXTURE / "genome.fa"),
    "--helixer",
    str(FIXTURE / "helixer_remapped.gff3"),
    "--scallop",
    str(FIXTURE / "scallop_geneset.gtf"),
    "--stringtie",
    str(FIXTURE / "stringtie_geneset.gtf"),
    "--orthodb",
    str(FIXTURE / "orthodb_geneset.gtf"),
    "--uniprot",
    str(FIXTURE / "uniprot_geneset.gtf"),
]


def gmb_run(out, *extra, preset="fungi_latest_settings", cwd=None):
    """Run `python -m gmb.cli.run` from a directory outside the repository."""
    cmd = [
        sys.executable,
        "-m",
        "gmb.cli.run",
        "--preset",
        preset,
        *EVIDENCE,
        "--gene-prefix",
        "ZT",
        "--output-dir",
        str(out),
        *map(str, extra),
    ]
    return subprocess.run(cmd, capture_output=True, text=True, cwd=str(cwd or out.parent))


def summary(out):
    return json.loads((out / "run_summary.json").read_text())


@pytest.fixture(scope="module")
def auto_run(tmp_path_factory):
    out = tmp_path_factory.mktemp("auto") / "run"
    r = gmb_run(out, "--long-read-alignments", READS)
    assert r.returncode == 0, r.stdout[-3000:] + r.stderr[-3000:]
    return out, r


@pytest.fixture(scope="module")
def no_lr_run(tmp_path_factory):
    out = tmp_path_factory.mktemp("nolr") / "run"
    r = gmb_run(out)
    assert r.returncode == 0, r.stdout[-3000:] + r.stderr[-3000:]
    return out, r


# ---------------------------------------------------------------------------
# fungi_latest_settings
# ---------------------------------------------------------------------------


def test_fungi_latest_settings_is_a_packaged_preset_layered_on_fungi():
    from importlib.resources import files

    assert "fungi_latest_settings" in list_build_presets()
    packaged = files("gmb.configs").joinpath("fungi_latest_settings.yaml")
    assert packaged.is_file()
    latest = load_config(None, "fungi_latest_settings")
    via_overlay = load_config([str(packaged)], "fungi")
    a, b = dump_config(latest), dump_config(via_overlay)
    assert a.pop("preset") == "fungi_latest_settings"
    b.pop("preset")
    assert a == b
    assert latest.scoring.primary_selection == "junction_supported"
    assert latest.scoring.isoform_cds_overlap == "drop"
    assert latest.scoring.min_cds_bp == 150
    assert latest.canonical_selection.prefer_build_primary is True
    # inherited from fungi, unchanged
    fungi = load_config(None, "fungi")
    assert latest.transcriptomic_filter == fungi.transcriptomic_filter
    assert latest.scoring.weights == fungi.scoring.weights


def test_no_conflicting_copy_of_latest_settings():
    assert not (GMB_DIR / "configs" / "fungi_evidence_integration.experimental.yaml").exists()
    copies = [
        p
        for p in GMB_DIR.rglob("fungi_latest_settings.yaml")
        if "build" not in p.parts and ".egg-info" not in str(p)
    ]
    assert copies == [GMB_DIR / "gmb" / "configs" / "fungi_latest_settings.yaml"]


@pytest.mark.integration
def test_preset_name_and_old_overlay_route_give_identical_builds(tmp_path):
    from importlib.resources import files

    overlay = str(files("gmb.configs").joinpath("fungi_latest_settings.yaml"))
    outs = {}
    for name, args in (
        ("preset", ["--preset", "fungi_latest_settings"]),
        ("overlay", ["--preset", "fungi", "--config", overlay]),
    ):
        outs[name] = tmp_path / name
        r = subprocess.run(
            [
                sys.executable,
                "-m",
                "gmb.cli.build",
                *args,
                *EVIDENCE,
                "--gene-prefix",
                "ZT",
                "--output-dir",
                str(outs[name]),
                "--no-log-file",
            ],
            capture_output=True,
            text=True,
            cwd=str(tmp_path),
        )
        assert r.returncode == 0, r.stdout[-2000:]
    for f in ("consensus.gff3", "evidence_attribution.tsv", "prot.fa"):
        assert (outs["preset"] / f).read_bytes() == (outs["overlay"] / f).read_bytes(), f


# ---------------------------------------------------------------------------
# gmb-run: no long reads / automatic consensus / prepared models
# ---------------------------------------------------------------------------


@pytest.mark.integration
def test_no_long_read_input_skips_preparation(no_lr_run):
    out, r = no_lr_run
    assert "Long-read consensus: SKIPPED (no long-read alignments supplied)" in r.stdout
    assert not (out / "longread_consensus").exists()
    s = summary(out)
    assert s["ok"] and s["stages_run"] == ["preflight", "build", "finalise"]
    assert s["long_read_models"] is None
    pf = json.loads((out / "preflight" / "preflight_report.json").read_text())
    assert not [t for t in pf["tracks"] if t.get("role") == "long_read_transcriptomic"]


@pytest.mark.integration
def test_per_read_alignments_are_collapsed_and_used(auto_run):
    out, r = auto_run
    assert "Long-read consensus: RUNNING" in r.stdout and "RAN ->" in r.stdout
    prep = json.loads((out / "longread_consensus" / "longread_preparation.json").read_text())
    gtf = out / "longread_consensus" / "minimap2_consensus.gtf"
    assert prep["status"] == "complete" and prep["action"] == "ran"
    assert prep["input"] == str(READS)
    assert prep["input_sha256"] == _sha256(READS)
    assert prep["consensus_gtf"] == str(gtf) and prep["consensus_gtf_sha256"] == _sha256(gtf)
    assert prep["consensus_models"] > 0
    assert prep["gmb_version"]
    assert not (out / "longread_consensus" / "_split_by_seqname").exists()
    s = summary(out)
    assert s["ok"]
    assert s["stages_run"] == ["longread_consensus", "preflight", "build", "finalise"]
    assert s["long_read_models"] == str(gtf)
    # the evidence preflight validated is the evidence the build used
    pf = json.loads((out / "preflight" / "preflight_report.json").read_text())
    assert any(t.get("path") == str(gtf) for t in pf["tracks"])
    manifest = json.loads((out / "build" / "run_manifest.json").read_text())
    assert str(gtf) in json.dumps(manifest)
    assert "Minimap2" in (out / "build" / "evidence_attribution.tsv").read_text()


@pytest.mark.integration
def test_consensus_uses_the_fungal_settings(auto_run):
    import dataclasses

    from gmb.pipeline.longread.config import load_longread_consensus_config

    out, _ = auto_run
    prep = json.loads((out / "longread_consensus" / "longread_preparation.json").read_text())
    resolved = prep["consensus_settings"]["resolved"]
    assert prep["consensus_settings"]["preset"] == "fungi"
    assert resolved == json.loads(
        json.dumps(dataclasses.asdict(load_longread_consensus_config(preset="fungi")))
    )
    build = load_config(None, "fungi_latest_settings").transcriptomic_filter
    assert resolved["max_span_length"] == build.max_transcript_length == 20000
    assert resolved["max_intron_length"] == build.max_intron_length == 3000


@pytest.mark.integration
def test_automatic_path_matches_manual_consensus_path(auto_run, tmp_path):
    out, _ = auto_run
    manual = tmp_path / "manual_lr"
    r = subprocess.run(
        [
            str(BIN / "gmb-longread-consensus"),
            "--input",
            str(READS),
            "--preset",
            "fungi",
            "--output-dir",
            str(manual),
            "--no-log-file",
        ],
        capture_output=True,
        text=True,
        cwd=str(tmp_path),
    )
    assert r.returncode == 0, r.stdout[-2000:]
    auto_gtf = out / "longread_consensus" / "minimap2_consensus.gtf"
    assert (manual / "minimap2_consensus.gtf").read_bytes() == auto_gtf.read_bytes()
    build = tmp_path / "build"
    r = subprocess.run(
        [
            sys.executable,
            "-m",
            "gmb.cli.build",
            "--preset",
            "fungi_latest_settings",
            *EVIDENCE,
            "--minimap2",
            str(manual / "minimap2_consensus.gtf"),
            "--gene-prefix",
            "ZT",
            "--output-dir",
            str(build),
            "--validate-fasta",
            "--no-log-file",
        ],
        capture_output=True,
        text=True,
        cwd=str(tmp_path),
    )
    assert r.returncode == 0, r.stdout[-2000:]
    for f in ("consensus.gff3", "prot.fa", "cds.fa"):
        assert (build / f).read_bytes() == (out / "build" / f).read_bytes(), f


@pytest.mark.integration
def test_rerun_reuses_matching_consensus_and_rebuilds_on_new_input(auto_run, tmp_path):
    out, _ = auto_run
    run = tmp_path / "run"
    shutil.copytree(out / "longread_consensus", run / "longread_consensus")
    # the copied preparation is valid for this input and settings: reused
    r = gmb_run(run, "--long-read-alignments", READS)
    assert r.returncode == 0 and "Long-read consensus: REUSED" in r.stdout
    # reuse points at the consensus inside THIS run directory, not the original one
    assert summary(run)["long_read_models"] == str(
        run / "longread_consensus" / "minimap2_consensus.gtf"
    )
    assert (run / "finalise" / "consensus.gff3").read_bytes() == (
        out / "finalise" / "consensus.gff3"
    ).read_bytes()
    # different reads (one removed): never reuse
    changed = tmp_path / "reads_changed.gtf"
    lines = READS.read_text().splitlines(keepends=True)
    first_tid = lines[0].split('transcript_id "')[1].split('"')[0]
    changed.write_text("".join(l for l in lines if f'transcript_id "{first_tid}"' not in l))
    r = gmb_run(run, "--long-read-alignments", changed)
    assert r.returncode == 0
    assert "RUNNING (input sha256 changed)" in r.stdout
    prep = json.loads((run / "longread_consensus" / "longread_preparation.json").read_text())
    assert prep["input_sha256"] == _sha256(changed)
    # a modified consensus file is not trusted
    gtf = run / "longread_consensus" / "minimap2_consensus.gtf"
    gtf.write_text(gtf.read_text() + "\n")
    r = gmb_run(run, "--long-read-alignments", changed)
    assert "RUNNING (previous consensus file missing or modified)" in r.stdout


@pytest.mark.integration
def test_prepared_models_are_used_directly(auto_run, tmp_path):
    out, _ = auto_run
    models = out / "longread_consensus" / "minimap2_consensus.gtf"
    run = tmp_path / "run"
    r = gmb_run(run, "--long-read-models", models)
    assert r.returncode == 0, r.stdout[-2000:]
    assert "SKIPPED (prepared long-read models supplied)" in r.stdout
    assert not (run / "longread_consensus").exists()
    assert summary(run)["long_read_models"] == str(models)
    assert (run / "build" / "consensus.gff3").read_bytes() == (
        out / "build" / "consensus.gff3"
    ).read_bytes()


# ---------------------------------------------------------------------------
# rejected inputs: clear failure, no gene build
# ---------------------------------------------------------------------------


def _assert_failed(r, run, stage, code, text):
    assert r.returncode == code, r.stdout[-3000:] + r.stderr[-2000:]
    assert text in r.stdout + r.stderr
    s = summary(run)
    assert s["failed_stage"] == stage and not s["ok"]
    assert not (run / "build").exists()


@pytest.mark.integration
def test_consensus_output_supplied_as_alignments_is_not_collapsed_again(auto_run, tmp_path):
    out, _ = auto_run
    run = tmp_path / "run"
    r = gmb_run(
        run, "--long-read-alignments", out / "longread_consensus" / "minimap2_consensus.gtf"
    )
    _assert_failed(r, run, "longread_consensus", 3, "is gmb-longread-consensus output")


@pytest.mark.integration
def test_per_read_alignments_supplied_as_models_are_refused(tmp_path):
    run = tmp_path / "run"
    r = gmb_run(run, "--long-read-models", READS)
    _assert_failed(r, run, "preflight", 4, "look like per-read alignments")


def _rewrite_strand(src, dst, fn):
    with open(src) as fi, open(dst, "w") as fo:
        for line in fi:
            f = line.split("\t")
            f[6] = fn(f[6], f)
            fo.write("\t".join(f))


@pytest.mark.integration
def test_unstranded_long_reads_are_rejected(tmp_path):
    bad = tmp_path / "unstranded.gtf"
    _rewrite_strand(READS, bad, lambda s, f: ".")
    run = tmp_path / "run"
    r = gmb_run(run, "--long-read-alignments", bad)
    _assert_failed(r, run, "longread_consensus", 3, "have no + / - strand")
    assert not (run / "longread_consensus" / "longread_preparation.json").exists()


@pytest.mark.integration
def test_wrongly_stranded_long_reads_fail_preflight_splice_check(tmp_path):
    flipped = tmp_path / "flipped.gtf"
    _rewrite_strand(READS, flipped, lambda s, f: {"+": "-", "-": "+"}.get(s, s))
    run = tmp_path / "run"
    r = gmb_run(run, "--long-read-alignments", flipped)
    _assert_failed(r, run, "preflight", 4, "splice_quality")


@pytest.mark.integration
def test_long_reads_on_unknown_sequences_are_rejected(tmp_path):
    bad = tmp_path / "renamed.gtf"
    with open(READS) as fi, open(bad, "w") as fo:
        for line in fi:
            fo.write("CM001196.1" + line[1:])
    run = tmp_path / "run"
    r = gmb_run(run, "--long-read-alignments", bad)
    _assert_failed(r, run, "longread_consensus", 3, "absent from the genome FASTA")


@pytest.mark.integration
def test_non_gtf_long_read_input_is_rejected(tmp_path):
    bad = tmp_path / "reads.bam.txt"
    bad.write_text("this is not a GTF\nline two\n")
    run = tmp_path / "run"
    r = gmb_run(run, "--long-read-alignments", bad)
    _assert_failed(r, run, "longread_consensus", 3, "do not have the 9 GTF columns")


def test_failed_consensus_stops_before_build(tmp_path):
    from gmb.api import run_gene_model_builder

    fake_bin = tmp_path / "bin"
    fake_bin.mkdir()
    script = fake_bin / "gmb-longread-consensus"
    script.write_text("#!/bin/sh\necho 'simulated failure' >&2\nexit 1\n")
    script.chmod(0o755)
    called = []
    result = run_gene_model_builder(
        genome=str(FIXTURE / "genome.fa"),
        backbone=str(FIXTURE / "helixer_remapped.gff3"),
        long_read_alignments=str(READS),
        preset="fungi_latest_settings",
        output_dir=str(tmp_path / "run"),
        bin_dir=str(fake_bin),
        log=called.append,
    )
    assert not result.ok and result.failed_stage == "longread_consensus"
    assert "gmb-longread-consensus failed (exit 1)" in result.failure_reason
    assert result.stages_run == ["longread_consensus"]
    assert not (tmp_path / "run" / "longread_consensus" / "longread_preparation.json").exists()
    assert (
        "simulated failure" in (tmp_path / "run" / "logs" / "longread_consensus.log").read_text()
    )


def test_api_rejects_both_long_read_inputs_and_unknown_consensus_preset(tmp_path):
    from gmb.api import run_gene_model_builder

    with pytest.raises(ValueError, match="supply long-read evidence once"):
        run_gene_model_builder(
            genome="g.fa",
            long_read=["m.gtf"],
            long_read_alignments="r.gtf",
            output_dir=str(tmp_path),
        )
    with pytest.raises(ValueError, match="no long-read consensus preset"):
        run_gene_model_builder(
            genome=str(FIXTURE / "genome.fa"),
            long_read_alignments=str(READS),
            preset="apicomplexa",
            output_dir=str(tmp_path / "x"),
            log=lambda m: None,
        )


def test_gmb_run_argument_errors_exit_2(tmp_path):
    r = subprocess.run(
        [sys.executable, "-m", "gmb.cli.run", "--genome", "x.fa"],
        capture_output=True,
        text=True,
        cwd=str(tmp_path),
    )
    assert r.returncode == 2
    r = subprocess.run(
        [
            sys.executable,
            "-m",
            "gmb.cli.run",
            "--preset",
            "nonexistent",
            *EVIDENCE,
            "--output-dir",
            str(tmp_path / "o"),
        ],
        capture_output=True,
        text=True,
        cwd=str(tmp_path),
    )
    assert r.returncode == 2 and "unknown preset" in r.stderr
    r = subprocess.run(
        [
            sys.executable,
            "-m",
            "gmb.cli.run",
            "--preset",
            "fungi_latest_settings",
            *EVIDENCE,
            "--long-read-alignments",
            str(READS),
            "--long-read-models",
            str(READS),
            "--output-dir",
            str(tmp_path / "o"),
        ],
        capture_output=True,
        text=True,
        cwd=str(tmp_path),
    )
    assert r.returncode == 2 and "not allowed with argument" in r.stderr


# ---------------------------------------------------------------------------
# finalisation, QC and provenance
# ---------------------------------------------------------------------------


@pytest.mark.integration
def test_handover_outputs_qc_and_provenance(auto_run):
    import yaml

    out, r = auto_run
    fin = out / "finalise"
    for f in (
        "canonical/consensus.canonical_annotated.gff3",
        "consensus.gff3",
        "cdna.fa",
        "cds.fa",
        "prot.fa",
        "fasta_qc_report.json",
        "utr_qc_report.json",
        "handover_manifest.json",
    ):
        assert (fin / f).exists(), f
    assert json.loads((fin / "fasta_qc_report.json").read_text())["pass"] is True
    assert json.loads((fin / "utr_qc_report.json").read_text())["invariant_ok"] is True
    resolved = yaml.safe_load((out / "build" / "resolved_config.yaml").read_text())
    assert resolved["preset"] == "fungi_latest_settings"
    assert resolved["scoring"]["primary_selection"] == "junction_supported"
    assert resolved["canonical_selection"]["prefer_build_primary"] is True
    s = summary(out)
    assert s["annotation"] == str(fin / "canonical" / "consensus.canonical_annotated.gff3")
    assert s["qc_passed"] is True and s["returncodes"] == {
        "longread_consensus": 0,
        "preflight": 0,
        "build": 0,
        "finalise": 0,
    }
    assert "Share these finalised files" in r.stdout
    result = subprocess.run(
        [
            sys.executable,
            str(GMB_DIR / "validation" / "check_handover.py"),
            str(fin),
            str(out / "build"),
        ],
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, result.stdout
    assert (out / "logs" / "gmb-run.log").exists()


@pytest.mark.integration
def test_historical_fungi_preset_still_runs_and_is_flagged(tmp_path):
    run = tmp_path / "run"
    r = gmb_run(run, preset="fungi")
    assert r.returncode == 0, r.stdout[-2000:]
    assert "historical fungal baseline" in r.stdout


def _sha256(path):
    import hashlib

    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


@pytest.mark.integration
@pytest.mark.slow
def test_wheel_install_runs_outside_the_repository(tmp_path, auto_run):
    """Build a wheel, install it into a fresh venv and run gmb-run from elsewhere."""
    wheel_dir, venv = tmp_path / "wheel", tmp_path / "venv"
    r = subprocess.run(
        [
            sys.executable,
            "-m",
            "pip",
            "wheel",
            str(GMB_DIR),
            "--no-deps",
            "--no-build-isolation",
            "-q",
            "-w",
            str(wheel_dir),
        ],
        capture_output=True,
        text=True,
        cwd=str(tmp_path),
    )
    if r.returncode != 0:
        pytest.skip(f"cannot build a wheel here: {r.stderr[-500:]}")
    subprocess.run([sys.executable, "-m", "venv", "--system-site-packages", str(venv)], check=True)
    vpy = venv / "bin" / "python"
    subprocess.run(
        [str(vpy), "-m", "pip", "install", "-q", "--no-deps", str(next(wheel_dir.glob("*.whl")))],
        check=True,
        capture_output=True,
    )
    where = subprocess.run(
        [str(vpy), "-c", "import gmb; print(gmb.__file__)"],
        capture_output=True,
        text=True,
        cwd=str(tmp_path),
    ).stdout
    assert str(venv) in where
    run = tmp_path / "wheel_run"
    cmd = [
        str(venv / "bin" / "gmb-run"),
        "--preset",
        "fungi_latest_settings",
        *EVIDENCE,
        "--long-read-alignments",
        str(READS),
        "--gene-prefix",
        "ZT",
        "--output-dir",
        str(run),
    ]
    r = subprocess.run(cmd, capture_output=True, text=True, cwd=str(tmp_path))
    assert r.returncode == 0, r.stdout[-3000:]
    out, _ = auto_run
    assert (run / "finalise" / "consensus.gff3").read_bytes() == (
        out / "finalise" / "consensus.gff3"
    ).read_bytes()
