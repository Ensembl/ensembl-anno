"""Hardening tests for the gmb-new-clade-config helper.

Cover the two safety guarantees the skill depends on:

* evaluation is tied to the exact resolved configuration the candidate build ran
  with (and to the overlay the developer wrote), frozen before any reference use;
* a machine overlay cannot carry biology, and a reference must be GFF3.

Everything runs on tiny temp files through the real CLI. No GMB build is run.
"""

from __future__ import annotations

import hashlib
import json
import os
import subprocess
import sys
from dataclasses import asdict

import pytest
import yaml

SCRIPT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                      "scripts", "gmb_clade.py")


def run(*args):
    r = subprocess.run([sys.executable, SCRIPT, *map(str, args)],
                       capture_output=True, text=True)
    return r.returncode, r.stdout + r.stderr


def sha(path):
    return hashlib.sha256(open(path, "rb").read()).hexdigest()


# ---------------------------------------------------------------------------
# fixtures
# ---------------------------------------------------------------------------

REF_GFF3 = """##gff-version 3
1\tref\tgene\t100\t400\t.\t+\t.\tID=gene:G1
1\tref\tmRNA\t100\t400\t.\t+\t.\tID=transcript:T1;Parent=gene:G1
1\tref\tCDS\t100\t400\t.\t+\t0\tParent=transcript:T1
1\tref\tgene\t1000\t2000\t.\t+\t.\tID=gene:G2
1\tref\tmRNA\t1000\t2000\t.\t+\t.\tID=transcript:T2;Parent=gene:G2
1\tref\tCDS\t1000\t1200\t.\t+\t0\tParent=transcript:T2
1\tref\tCDS\t1800\t2000\t.\t+\t0\tParent=transcript:T2
"""

REF_GTF = ('1\tref\tCDS\t100\t400\t.\t+\t0\tgene_id "G1"; transcript_id "T1";\n'
           '1\tref\tCDS\t1000\t1200\t.\t+\t0\tgene_id "G2"; transcript_id "T2";\n')


@pytest.fixture
def workspace(tmp_path):
    """A candidate 'build' as gmb-build would leave it, plus a comparison dir."""
    from gmb.pipeline.config import load_config

    overlay = tmp_path / "configs" / "testclade.yaml"
    overlay.parent.mkdir()
    overlay.write_text("transcriptomic_filter:\n  max_intron_length: 3000\n")

    build = tmp_path / "out" / "candidate" / "build"
    build.mkdir(parents=True)
    resolved = build / "resolved_config.yaml"
    resolved.write_text(yaml.safe_dump(
        asdict(load_config(str(overlay), preset="standard")), sort_keys=True))
    (build / "run_manifest.json").write_text(json.dumps({
        "config_overlays": [str(overlay)],
        "resolved_config_sha256": sha(resolved)}))

    cmp_dir = tmp_path / "out" / "candidate" / "comparison"
    cmp_dir.mkdir()
    (cmp_dir / "comparison_summary.tsv").write_text(
        "metric\tvalue\ntotal_reference_genes\t2\ntotal_consensus_genes\t2\n"
        "sens_cds_exact_match_count\t1\nsens_locus_detection_rate\t1.0\n")
    (cmp_dir / "comparison_details.tsv").write_text(
        "source\tgene_id\tclassification_cds\tmatched_id\n"
        "reference\tG1\tExact_Match\tQ1\nreference\tG2\tPartial_Match\tQ2\n"
        "consensus\tQ1\tExact_Match\tG1\nconsensus\tQ2\tPartial_Match\tG2\n")

    ref = tmp_path / "reference.gff3"
    ref.write_text(REF_GFF3)
    return {"root": tmp_path, "overlay": overlay, "resolved": resolved, "build": build,
            "comparison": cmp_dir, "reference": ref, "out": tmp_path / "out"}


def freeze(w):
    return run("freeze", w["overlay"], "--resolved-config", w["resolved"], "--out", w["out"])


def evaluate(w, reference=None):
    return run("evaluate", "--freeze", w["out"] / "config_freeze.json",
               "--comparison", w["comparison"], "--reference", reference or w["reference"],
               "--out", w["out"])


def no_evaluation_output(w):
    return not (w["out"] / "evaluation").exists()


# ---------------------------------------------------------------------------
# A-D: the resolved-config freeze
# ---------------------------------------------------------------------------

class TestResolvedConfigFreeze:
    def test_A_freeze_records_both_hashes_and_evaluation_is_allowed(self, workspace):
        w = workspace
        rc, out = freeze(w)
        assert rc == 0, out
        rec = json.load(open(w["out"] / "config_freeze.json"))
        assert rec["source_overlay_sha256"] == sha(w["overlay"])
        assert rec["resolved_config_sha256"] == sha(w["resolved"])
        assert rec["resolved_config"] == str(w["resolved"])
        assert rec["build_dir"] == str(w["build"])
        rc, out = evaluate(w)
        assert rc == 0, out
        result = json.load(open(w["out"] / "evaluation" / "evaluation.json"))
        assert result["resolved_config_sha256"] == sha(w["resolved"])

    def test_B_source_overlay_edited_after_freeze_is_refused(self, workspace):
        w = workspace
        assert freeze(w)[0] == 0
        w["overlay"].write_text(w["overlay"].read_text() + "orf:\n  min_codons: 33\n")
        rc, out = evaluate(w)
        assert rc != 0 and "SOURCE OVERLAY CHANGED" in out
        assert no_evaluation_output(w)

    def test_C_resolved_config_edited_after_freeze_is_refused(self, workspace):
        w = workspace
        assert freeze(w)[0] == 0
        w["resolved"].write_text(w["resolved"].read_text() + "# edited\n")
        rc, out = evaluate(w)
        assert rc != 0 and "RESOLVED BUILD CONFIG CHANGED" in out
        assert "SOURCE OVERLAY" not in out      # the two failures are distinguishable
        assert no_evaluation_output(w)

    def test_D_freeze_refuses_a_missing_resolved_config(self, workspace):
        w = workspace
        w["resolved"].unlink()
        rc, out = freeze(w)
        assert rc != 0 and "resolved config not found" in out
        assert not (w["out"] / "config_freeze.json").exists()

    def test_D_evaluation_refuses_when_resolved_config_disappears(self, workspace):
        w = workspace
        assert freeze(w)[0] == 0
        w["resolved"].unlink()
        rc, out = evaluate(w)
        assert rc != 0 and "RESOLVED BUILD CONFIG missing" in out
        assert no_evaluation_output(w)

    def test_freeze_is_required_before_evaluation(self, workspace):
        rc, out = evaluate(workspace)
        assert rc != 0 and "no freeze record" in out
        assert no_evaluation_output(workspace)

    def test_old_overlay_only_freeze_record_is_refused(self, workspace):
        w = workspace
        (w["out"] / "config_freeze.json").write_text(json.dumps(
            {"config": str(w["overlay"]), "config_sha256": sha(w["overlay"])}))
        rc, out = evaluate(w)
        assert rc != 0 and "predates resolved-config freezing" in out

    def test_resolved_config_from_a_build_without_this_overlay_is_refused(self, workspace):
        w = workspace
        (w["build"] / "run_manifest.json").write_text(json.dumps({
            "config_overlays": ["/somewhere/else.yaml"],
            "resolved_config_sha256": sha(w["resolved"])}))
        rc, out = freeze(w)
        assert rc != 0 and "was not run with" in out

    def test_resolved_config_edited_after_its_build_is_refused_at_freeze(self, workspace):
        w = workspace
        w["resolved"].write_text(w["resolved"].read_text() + "# edited\n")
        rc, out = freeze(w)
        assert rc != 0 and "run_manifest.json" in out

    def test_resolved_config_lacking_the_overlay_values_is_refused(self, workspace):
        w = workspace
        (w["build"] / "run_manifest.json").unlink()        # force the value check
        w["overlay"].write_text("transcriptomic_filter:\n  max_intron_length: 4242\n")
        rc, out = freeze(w)
        assert rc != 0 and "transcriptomic_filter.max_intron_length" in out

    def test_freeze_is_refused_once_the_reference_has_been_consulted(self, workspace):
        w = workspace
        (w["out"] / "evaluation").mkdir()
        rc, out = freeze(w)
        assert rc != 0 and "already been consulted" in out


# ---------------------------------------------------------------------------
# E-F: machine overlay guard
# ---------------------------------------------------------------------------

@pytest.fixture
def manifest_dir(tmp_path):
    (tmp_path / "genome.fa").write_text(">1\n" + "ACGT" * 50 + "\n")
    (tmp_path / "backbone.gff3").write_text(
        "1\tp\tgene\t1\t90\t.\t+\t.\tID=g1\n"
        "1\tp\tmRNA\t1\t90\t.\t+\t.\tID=t1;Parent=g1\n"
        "1\tp\texon\t1\t90\t.\t+\t.\tParent=t1\n")
    return tmp_path


def write_manifest(d, machine_overlay_text):
    (d / "machine.yaml").write_text(machine_overlay_text)
    (d / "manifest.yaml").write_text(
        "clade: testclade\ngenome: genome.fa\nmachine_overlay: machine.yaml\n"
        "evidence:\n  - {path: backbone.gff3, tool: PredictorFoo, role: backbone}\n")
    return d / "manifest.yaml"


class TestMachineOverlayGuard:
    def test_E_execution_only_keys_are_accepted(self, manifest_dir):
        m = write_manifest(manifest_dir,
                           "protein_validation:\n"
                           "  diamond_path: /opt/diamond\n  psauron_path: /opt/psauron\n"
                           "  diamond_db: /data/proteins.dmnd\n  psauron_use_cpu: true\n"
                           "qc:\n  workers: 8\n")
        rc, out = run("args", m)
        assert rc == 0, out

    @pytest.mark.parametrize("text,key", [
        ("scoring:\n  weights:\n    backbone: 9.0\n", "scoring.weights.backbone"),
        ("scoring:\n  backbone_intron_rescue: 'on'\n", "scoring.backbone_intron_rescue"),
        ("scoring:\n  structural_corroboration: true\n", "scoring.structural_corroboration"),
        ("scoring:\n  protein_support_mode: cds_span_compatible\n", "scoring.protein_support_mode"),
        ("scoring:\n  longread_disposition: reject\n", "scoring.longread_disposition"),
        ("transcriptomic_filter:\n  max_intron_length: 500\n", "transcriptomic_filter.max_intron_length"),
        ("orf:\n  min_codons: 20\n", "orf.min_codons"),
        ("utr:\n  require_end_support: false\n", "utr.require_end_support"),
        # enabling a scoring stage is biology, even though it sits beside the paths
        ("protein_validation:\n  diamond_path: /opt/d\n  enabled: true\n", "protein_validation.enabled"),
    ])
    def test_F_selection_affecting_keys_are_rejected(self, manifest_dir, text, key):
        rc, out = run("args", write_manifest(manifest_dir, text))
        assert rc != 0
        assert key in out and "not machine/execution" in out

    def test_missing_machine_overlay_is_rejected(self, manifest_dir):
        m = write_manifest(manifest_dir, "qc:\n  workers: 2\n")
        (manifest_dir / "machine.yaml").unlink()
        rc, out = run("args", m)
        assert rc != 0 and "machine_overlay not found" in out

    def test_allowlist_names_only_real_schema_keys(self):
        """The allowlist must not invent keys: every entry exists in the current schema."""
        import importlib.util

        from gmb.pipeline.config import load_config

        spec = importlib.util.spec_from_file_location("gmb_clade", SCRIPT)
        mod = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(mod)
        real = set(mod.flatten(asdict(load_config(preset="standard"))))
        assert mod.MACHINE_OVERLAY_ALLOWED <= real, mod.MACHINE_OVERLAY_ALLOWED - real


# ---------------------------------------------------------------------------
# reference format
# ---------------------------------------------------------------------------

class TestReferenceFormat:
    def test_gtf_reference_is_refused_before_any_output(self, workspace):
        w = workspace
        assert freeze(w)[0] == 0
        gtf = w["root"] / "reference.gtf"
        gtf.write_text(REF_GTF)
        rc, out = evaluate(w, reference=gtf)
        assert rc != 0 and "must currently be GFF3" in out
        assert no_evaluation_output(w)

    def test_gtf_content_under_a_gff3_name_is_still_refused(self, workspace):
        w = workspace
        assert freeze(w)[0] == 0
        disguised = w["root"] / "disguised.gff3"
        disguised.write_text(REF_GTF)
        rc, out = evaluate(w, reference=disguised)
        assert rc != 0 and "must currently be GFF3" in out
        assert no_evaluation_output(w)

    def test_missing_reference_is_refused(self, workspace):
        w = workspace
        assert freeze(w)[0] == 0
        rc, out = evaluate(w, reference=w["root"] / "nope.gff3")
        assert rc != 0 and "reference not found" in out
        assert no_evaluation_output(w)

    def test_gff3_without_coding_structure_is_refused(self, workspace):
        w = workspace
        assert freeze(w)[0] == 0
        empty = w["root"] / "genes_only.gff3"
        empty.write_text("##gff-version 3\n1\tref\tgene\t1\t100\t.\t+\t.\tID=g1\n")
        rc, out = evaluate(w, reference=empty)
        assert rc != 0 and "not a usable protein-coding GFF3" in out
        assert no_evaluation_output(w)

    def test_isoforms_do_not_inflate_the_exon_count(self, tmp_path):
        """Two single-exon isoforms must leave the gene single-exon."""
        import importlib.util

        spec = importlib.util.spec_from_file_location("gmb_clade", SCRIPT)
        mod = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(mod)
        ref = tmp_path / "iso.gff3"
        ref.write_text(
            "1\tr\tgene\t1\t300\t.\t+\t.\tID=gene:G\n"
            "1\tr\tmRNA\t1\t300\t.\t+\t.\tID=a;Parent=gene:G\n"
            "1\tr\tCDS\t1\t300\t.\t+\t0\tParent=a\n"
            "1\tr\tmRNA\t1\t240\t.\t+\t.\tID=b;Parent=gene:G\n"
            "1\tr\tCDS\t1\t240\t.\t+\t0\tParent=b\n")
        assert mod._ref_cds_exon_counts(str(ref))["G"] == 1
