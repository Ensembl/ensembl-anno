"""Regression tests for the release-preparation clean-up.

Covers: import integrity after removing gmb.compare, the declared runtime
dependencies, the consistent --preset default, UTR end support (index, same
sequence only), locus clustering modes, inert configuration keys, and an
in-process build through the generic --backbone slot.
"""

from __future__ import annotations

import importlib
import inspect
import os
import pkgutil
import re
import subprocess
import sys
from pathlib import Path

import pandas as pd
import pytest
import yaml

import gmb
from gmb.pipeline import builder
from gmb.pipeline.builder import (
    build_transcript_end_index,
    cluster_candidate_loci,
    compute_utr_end_support,
)
from gmb.pipeline.config import DEFAULT_PRESET, PipelineConfig, load_config

GMB_DIR = Path(__file__).resolve().parent.parent
FIXTURE = GMB_DIR / "tests" / "fixtures" / "z_tritici_region1"


# ---------------------------------------------------------------------------
# Import integrity and packaging
# ---------------------------------------------------------------------------


def _all_gmb_modules():
    return sorted(m.name for m in pkgutil.walk_packages(gmb.__path__, "gmb."))


@pytest.mark.parametrize("module", _all_gmb_modules())
def test_every_module_imports(module):
    importlib.import_module(module)


def test_comparison_package_is_gone():
    with pytest.raises(ModuleNotFoundError):
        importlib.import_module("gmb.compare")


def test_production_modules_do_not_import_matplotlib():
    """Plotting belonged to the removed comparison tools; the build must not need it."""
    code = (
        "import importlib, pkgutil, sys, gmb\n"
        "for m in pkgutil.walk_packages(gmb.__path__, 'gmb.'):\n"
        "    importlib.import_module(m.name)\n"
        "print('matplotlib' in sys.modules)\n"
    )
    out = subprocess.run(
        [sys.executable, "-c", code], capture_output=True, text=True, cwd=str(GMB_DIR), check=True
    )
    assert out.stdout.strip() == "False"


def _pyproject_section(name: str) -> str:
    text = (GMB_DIR / "pyproject.toml").read_text()
    m = re.search(rf"^\[{re.escape(name)}\]\n(.*?)(?=^\[)", text, re.S | re.M)
    assert m, f"[{name}] missing from pyproject.toml"
    return m.group(1)


def test_console_scripts_resolve():
    scripts = dict(
        re.findall(r'^([\w-]+)\s*=\s*"([^"]+)"', _pyproject_section("project.scripts"), re.M)
    )
    assert "gmb-visualize" not in scripts
    for name, target in scripts.items():
        mod, func = target.split(":")
        assert callable(getattr(importlib.import_module(mod), func)), name


def test_declared_dependencies_match_imports():
    """Every third-party import is declared, and nothing declared is unused."""
    import ast

    declared = set(
        re.findall(
            r'"([A-Za-z0-9_-]+)',
            re.search(
                r"^dependencies = \[(.*?)\]", (GMB_DIR / "pyproject.toml").read_text(), re.S | re.M
            ).group(1),
        )
    )
    import_name = {"pyyaml": "yaml"}
    declared_imports = {import_name.get(d.lower(), d.lower()) for d in declared}

    used = set()
    for path in (GMB_DIR / "gmb").rglob("*.py"):
        for node in ast.walk(ast.parse(path.read_text())):
            if isinstance(node, ast.Import):
                used.update(a.name.split(".")[0] for a in node.names)
            elif isinstance(node, ast.ImportFrom) and node.level == 0 and node.module:
                used.add(node.module.split(".")[0])
    third_party = {u for u in used if u not in sys.stdlib_module_names and u != "gmb"}
    assert third_party <= declared_imports, third_party - declared_imports
    assert declared_imports <= third_party, declared_imports - third_party


def test_gmb_compare_stub_points_to_annotation_qc(capsys):
    from gmb.cli.compare import main

    assert main([]) == 2
    err = capsys.readouterr().err
    assert "annotation-qc pairwise-compare" in err
    assert "ensembl-genes" in err


# ---------------------------------------------------------------------------
# One default preset for every entry point
# ---------------------------------------------------------------------------


def test_load_config_default_is_shared_constant():
    assert inspect.signature(load_config).parameters["preset"].default == DEFAULT_PRESET


def test_build_and_preflight_share_the_default_preset(monkeypatch):
    from gmb.cli.preflight import build_parser

    monkeypatch.setattr(sys, "argv", ["gmb-build", "--output-dir", "x"])
    assert builder.parse_args().preset is None  # resolved to DEFAULT_PRESET in main()
    assert build_parser().parse_args(["--genome", "g.fa"]).preset is None


def test_preflight_resolves_omitted_preset_to_build_default(tmp_path, capsys):
    from gmb.cli import preflight

    genome = tmp_path / "g.fa"
    genome.write_text(">1\nACGT\n")
    preflight.main(["--genome", str(genome), "--json", "--allow-fail", "--no-log-file"])
    assert f"using '{DEFAULT_PRESET}'" in capsys.readouterr().err


def test_api_passes_one_preset_to_every_stage(monkeypatch, tmp_path):
    from gmb import api

    seen = []

    def fake_run(cmd, log_path):
        seen.append(cmd)
        return 0

    monkeypatch.setattr(api, "_run", fake_run)
    api.run_gene_model_builder(
        genome="g.fa",
        backbone="b.gff3",
        preset="fungi",
        output_dir=str(tmp_path),
        validate_fasta=False,
    )
    presets = [c[c.index("--preset") + 1] for c in seen if "--preset" in c]
    assert presets == ["fungi", "fungi"]  # preflight and build; finalise uses the build's config


# ---------------------------------------------------------------------------
# UTR end support
# ---------------------------------------------------------------------------


def _utr_cfg():
    cfg = PipelineConfig()
    cfg.utr.require_end_support = True
    cfg.utr.end_support_mode = "multisource_end_agreement"
    cfg.utr.end_support_sources = ["Scallop", "StringTie"]
    cfg.utr.end_tolerance_bp = 10
    cfg.utr.require_multisource_for_utr_5p = True
    cfg.utr.require_multisource_for_utr_3p = True
    cfg.utr.fallback_policy_when_unsupported = "drop_utr"
    return cfg


def _exons(rows):
    return pd.DataFrame(
        rows, columns=["transcript_id", "Source", "Chromosome", "Strand", "Start", "End"]
    )


def test_end_on_another_sequence_is_not_support():
    """Regression: the builder used to pass the genome-wide table with no chromosome
    check, so an end at the same coordinate on another contig counted as agreement
    (5.8% of 5' and 6.7% of 3' short-read ends on Z. tritici)."""
    df = _exons([("sc_1", "Scallop", "2", "+", 1000, 2000)])
    model = {"id": "st_1", "chrom": "1", "strand": "+", "start": 1000, "end": 2000}
    res = compute_utr_end_support(model, df, _utr_cfg())
    assert res["supported_5p"] is False and res["supported_3p"] is False

    df = _exons([("sc_1", "Scallop", "1", "+", 1003, 1995)])
    res = compute_utr_end_support(model, df, _utr_cfg())
    assert res["supported_5p"] is True and res["supported_3p"] is True


def test_model_never_supports_itself():
    df = _exons([("st_1", "StringTie", "1", "+", 1000, 2000)])
    model = {"id": "st_1", "chrom": "1", "strand": "+", "start": 1000, "end": 2000}
    res = compute_utr_end_support(model, df, _utr_cfg())
    assert res["supported_5p"] is False and res["supported_3p"] is False


def test_minus_strand_ends_and_source_filter():
    df = _exons(
        [
            ("sc_1", "Scallop", "1", "-", 500, 900),
            ("sc_1", "Scallop", "1", "-", 1000, 2005),
            ("hx_1", "Helixer", "1", "-", 1000, 2000),  # not an end-support source
        ]
    )
    model = {"id": "st_1", "chrom": "1", "strand": "-", "start": 495, "end": 2000}
    res = compute_utr_end_support(model, df, _utr_cfg())
    # 5' of a minus-strand model is its max coordinate; 3' its min.
    assert res["supported_5p"] is True and res["supported_3p"] is True
    res = compute_utr_end_support(dict(model, start=400), df, _utr_cfg())
    assert res["supported_3p"] is False and res["action_3p"] == "dropped"


def test_prebuilt_index_matches_per_call_scan():
    rng = pd.Series(range(60))
    df = _exons(
        [
            (
                f"t{i}",
                "Scallop" if i % 2 else "StringTie",
                str(i % 3 + 1),
                "+-"[i % 2],
                1000 + 7 * i,
                3000 + 11 * i,
            )
            for i in rng
        ]
    )
    cfg = _utr_cfg()
    index = build_transcript_end_index(df, cfg.utr.end_support_sources)
    for i in rng:
        model = {
            "id": f"m{i}",
            "chrom": str(i % 3 + 1),
            "strand": "+-"[i % 2],
            "start": 1000 + 7 * i + (i % 5) * 4,
            "end": 3000 + 11 * i - (i % 4) * 5,
        }
        assert compute_utr_end_support(model, df, cfg) == compute_utr_end_support(
            model, df, cfg, end_index=index
        )


# ---------------------------------------------------------------------------
# Locus clustering
# ---------------------------------------------------------------------------


def _fixture_locus():
    """Three 2-exon models that agree on one intron (Z. tritici chr1:246,121-247,117).

    No candidate spans the shared intron, so exon-overlap clustering puts the
    two exons of every model into different loci. The Scallop/StringTie CDS
    matches reference Mycgr3T88624 exactly; the Helixer CDS does not.
    """
    import contextlib
    import io

    from gmb.pipeline.builder import load_evidence

    with contextlib.redirect_stdout(io.StringIO()):
        frames = [
            load_evidence(str(FIXTURE / f), label)
            for f, label in (
                ("scallop_geneset.gtf", "Scallop"),
                ("stringtie_geneset.gtf", "StringTie"),
                ("helixer_remapped.gff3", "Helixer"),
            )
        ]
    keep = {"Scallop_MSTRG.22.1", "StringTie_MSTRG.32.1", "Helixer__CM001196.1_000028.1"}
    exons = pd.concat([f[0] for f in frames], ignore_index=True)
    cds = pd.concat(
        [f[1] for f in frames if f[1] is not None and not f[1].empty], ignore_index=True
    )
    return exons[exons.transcript_id.isin(keep)], cds[cds.transcript_id.isin(keep)]


def test_exon_overlap_clustering_splits_agreeing_multi_exon_models():
    """Documents the validated-baseline behaviour (scoring.locus_clustering default)."""
    exons, _ = _fixture_locus()
    clusters = cluster_candidate_loci(exons, "exon_overlap")
    assert (clusters.groupby("transcript_id")["Cluster"].nunique() == 2).all()


def test_transcript_linked_clustering_keeps_models_whole():
    exons, _ = _fixture_locus()
    clusters = cluster_candidate_loci(exons, "transcript_linked")
    assert clusters["Cluster"].nunique() == 1
    assert (clusters["Count"] == len(clusters)).all()


def test_transcript_linked_keeps_unlinked_nested_gene_separate():
    exons = _exons(
        [
            ("outer", "Scallop", "1", "+", 100, 200),
            ("outer", "Scallop", "1", "+", 900, 1000),
            ("nested", "StringTie", "1", "-", 400, 600),
        ]
    )
    clusters = cluster_candidate_loci(exons, "transcript_linked")
    assert clusters.groupby("transcript_id")["Cluster"].nunique().max() == 1
    assert clusters["Cluster"].nunique() == 2


def test_transcript_linked_lets_the_reference_exact_model_be_selected():
    """On real data the correct spliced short-read model is only selectable whole."""
    from gmb.pipeline.annotate_cds_utrs import annotate_all_transcripts, load_genome
    from gmb.pipeline.scoring import select_isoforms

    exons, cds = _fixture_locus()
    genome = load_genome(str(FIXTURE / "genome.fa"))
    cfg = load_config(None, "fungi")
    ann = annotate_all_transcripts(exons, genome, cds, min_codons=cfg.orf.min_codons)
    cand_cds = {t: a.get("cds") or [] for t, a in ann.items()}
    reference_cds = [(246302, 246499), (246563, 246675)]  # Mycgr3T88624, 0-based
    assert sorted(cand_cds["Scallop_MSTRG.22.1"]) == reference_cds
    protein = set(exons.transcript_id)

    def selected(mode):
        out = set()
        for _cid, locus in cluster_candidate_loci(exons, mode).groupby("Cluster"):
            for gene in select_isoforms(locus, cfg, protein, genome, candidate_cds=cand_cds):
                out.update(m["id"] for m in gene)
        return out

    assert not {"Scallop_MSTRG.22.1", "StringTie_MSTRG.32.1"} & selected("exon_overlap")
    assert {"Scallop_MSTRG.22.1", "StringTie_MSTRG.32.1"} & selected("transcript_linked")


def test_locus_clustering_config():
    # fungi: whole-candidate scoring, paired with its transcript length filter.
    # Neutral base and apicomplexa keep the validated exon_overlap behaviour.
    assert load_config(None, "standard").scoring.locus_clustering == "exon_overlap"
    assert load_config(None, "apicomplexa").scoring.locus_clustering == "exon_overlap"
    assert load_config(None, "fungi").scoring.locus_clustering == "transcript_linked"


def test_invalid_locus_clustering_is_fatal(tmp_path):
    overlay = tmp_path / "bad.yaml"
    overlay.write_text("scoring:\n  locus_clustering: span\n")
    with pytest.raises(ValueError, match="locus_clustering"):
        load_config(str(overlay), "fungi")


# ---------------------------------------------------------------------------
# Inert configuration keys
# ---------------------------------------------------------------------------

from gmb.pipeline.config import INERT_CONFIG_KEYS  # noqa: E402


def _unread_config_keys():
    import dataclasses

    src = "\n".join(
        p.read_text() for p in (GMB_DIR / "gmb").rglob("*.py") if p.name != "config.py"
    )
    unread = set()

    def walk(dc, prefix):
        for f in dataclasses.fields(dc):
            value = getattr(dc, f.name)
            if dataclasses.is_dataclass(value):
                walk(value, f"{prefix}{f.name}.")
            elif not re.search(rf"\.{f.name}\b|[\"']{f.name}[\"']", src):
                unread.add(prefix + f.name)

    walk(PipelineConfig(), "")
    return unread


def test_inert_config_keys_are_exactly_the_documented_set():
    """If you wire a key up, remove it from config.INERT_CONFIG_KEYS and its
    "NOT APPLIED" note in standard.yaml; if you add an unread key, list it."""
    assert _unread_config_keys() == set(INERT_CONFIG_KEYS)
    assert set(INERT_CONFIG_KEYS.values()) == {"deprecated", "unsupported"}


def test_inert_config_keys_are_marked_in_standard_yaml():
    text = (GMB_DIR / "gmb" / "configs" / "standard.yaml").read_text()
    data = yaml.safe_load(text)
    for key in INERT_CONFIG_KEYS:
        node = data
        for part in key.split("."):
            node = node[part]  # still present, so existing configs keep loading
    assert text.count("NOT APPLIED") >= 8


@pytest.mark.parametrize(
    "key,category",
    [
        ("orf.allow_partial_5", UserWarning),
        ("qc.workers", FutureWarning),
        ("export.write_cds", FutureWarning),
    ],
)
def test_setting_an_inert_key_warns(tmp_path, key, category):
    section, name = key.split(".")
    overlay = tmp_path / "o.yaml"
    overlay.write_text(f"{section}:\n  {name}: false\n")
    with pytest.warns(category, match=key.replace(".", r"\.")):
        load_config(str(overlay), "standard")


def test_reloading_a_resolved_config_does_not_warn(tmp_path):
    import warnings

    from gmb.pipeline.config import dump_config

    data = dump_config(load_config(None, "fungi"))
    data.pop("preset", None)
    path = tmp_path / "resolved_config.yaml"
    path.write_text(yaml.safe_dump(data))
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        load_config(str(path), preset=None, warn_inert=False)
    with pytest.warns(FutureWarning):
        load_config(str(path), preset=None)


def test_shipped_presets_and_examples_set_no_inert_keys():
    import warnings

    with warnings.catch_warnings():
        warnings.simplefilter("error")
        for preset in ("standard", "fungi", "apicomplexa"):
            load_config(None, preset)
        for name in ("fungi_default.yaml", "apicomplexa_first_pass.yaml"):
            load_config(str(GMB_DIR / "configs" / name), "standard")


# ---------------------------------------------------------------------------
# Transcriptomic filter: span and single-exon rules
# ---------------------------------------------------------------------------


def _filter(rows, **tfilter):
    from gmb.pipeline.evidence_filter import filter_chimeras

    cfg = load_config(None, "standard")
    for k, v in tfilter.items():
        setattr(cfg.transcriptomic_filter, k, v)
    stats = {}
    kept = filter_chimeras(_exons(rows), cfg, stats)
    return set(kept["transcript_id"]), stats


READ_THROUGH = [  # 25 kb, no intron > 3 kb: one long "exon" crosses the intergenic gap
    ("rt", "Scallop", "1", "+", 1000, 1200),
    ("rt", "Scallop", "1", "+", 1300, 16000),
    ("rt", "Scallop", "1", "+", 16100, 26000),
    ("ok", "Scallop", "1", "+", 30000, 30500),
    ("ok", "Scallop", "1", "+", 31000, 32000),
    ("se", "StringTie", "1", "-", 40000, 41000),
]


def test_max_transcript_length_removes_read_through_chimeras():
    kept, stats = _filter(READ_THROUGH, max_transcript_length=20000, max_intron_length=3000)
    assert kept == {"ok", "se"}
    assert stats["chimeras_long_span"] == 1 and stats["chimeras_large_intron"] == 0


def test_max_transcript_length_null_disables():
    kept, stats = _filter(READ_THROUGH, max_transcript_length=None, max_intron_length=3000)
    assert kept == {"rt", "ok", "se"} and stats["chimeras_long_span"] == 0


def test_allow_single_exon_false_drops_single_exon_transcripts():
    kept, stats = _filter(READ_THROUGH, max_transcript_length=None, allow_single_exon=False)
    assert "se" not in kept and stats["single_exon_removed"] == 1


def test_fungi_preset_applies_20kb_and_apicomplexa_disables_it():
    assert load_config(None, "fungi").transcriptomic_filter.max_transcript_length == 20000
    assert load_config(None, "apicomplexa").transcriptomic_filter.max_transcript_length is None


# ---------------------------------------------------------------------------
# Detached isoforms
# ---------------------------------------------------------------------------


def _gene(gid, *spans):
    rows = [
        {
            "Feature": "gene",
            "ID": gid,
            "Parent": "",
            "Start": min(s for s, _ in spans),
            "End": max(e for _, e in spans),
        }
    ]
    for i, (s, e) in enumerate(spans, start=1):
        tid = f"{gid}.t{i}"
        rows.append({"Feature": "mRNA", "ID": tid, "Parent": gid, "Start": s, "End": e})
        rows.append({"Feature": "exon", "ID": f"{tid}.exon1", "Parent": tid, "Start": s, "End": e})
    return rows


def test_detached_isoform_is_removed_and_primary_kept():
    """Regression for a read-through alternate trimmed by validation to a piece 27 kb
    from its gene's primary (Z. tritici ZT_00017): the gene must not join both loci."""
    from gmb.pipeline.gff3_validate import drop_detached_isoforms, recompute_gene_bounds

    rows = _gene("G1", (1000, 2000), (28000, 30000), (1500, 2500)) + _gene("G2", (5000, 6000))
    rows, stats = drop_detached_isoforms(rows)
    ids = {r["ID"] for r in rows}
    assert "G1.t2" not in ids and "G1.t2.exon1" not in ids
    assert {"G1.t1", "G1.t3", "G2.t1"} <= ids
    assert stats == {"genes_with_detached_isoforms": 1, "detached_isoforms_removed": 1}
    rows, _ = recompute_gene_bounds(rows)
    g1 = next(r for r in rows if r["ID"] == "G1")
    assert (g1["Start"], g1["End"]) == (1000, 2500)


def test_detached_check_keeps_chained_and_touching_cases_correct():
    from gmb.pipeline.gff3_validate import drop_detached_isoforms

    chained = _gene("G", (0, 100), (90, 200), (190, 300))  # all connected via t2
    assert drop_detached_isoforms(chained)[1]["detached_isoforms_removed"] == 0
    touching = _gene("H", (0, 100), (100, 200))  # half-open: touching does not overlap
    assert drop_detached_isoforms(touching)[1]["detached_isoforms_removed"] == 1


# ---------------------------------------------------------------------------
# In-process build through the generic --backbone slot
# ---------------------------------------------------------------------------


def _write_inputs(tmp_path):
    utr5, coding, utr3 = "AAA" * 10, "ATG" + "GCT" * 49 + "TAA", "CCC" * 20
    (tmp_path / "genome.fa").write_text(f">1\n{utr5 + coding + utr3 + 'N' * 258}\n")
    (tmp_path / "backbone.gff3").write_text(
        "##gff-version 3\n"
        "1\tAUGUSTUS\tgene\t1\t243\t.\t+\t.\tID=g1\n"
        "1\tAUGUSTUS\tmRNA\t1\t243\t.\t+\t.\tID=t1;Parent=g1\n"
        "1\tAUGUSTUS\texon\t1\t243\t.\t+\t.\tID=t1.e1;Parent=t1\n"
        "1\tAUGUSTUS\tCDS\t31\t183\t.\t+\t0\tID=t1.c1;Parent=t1\n"
    )
    for name, src in (("scallop", "Scallop"), ("stringtie", "StringTie")):
        (tmp_path / f"{name}.gtf").write_text(
            f'1\t{src}\ttranscript\t1\t243\t.\t+\t.\tgene_id "{name}_g"; transcript_id "{name}_t";\n'
            f'1\t{src}\texon\t1\t243\t.\t+\t.\tgene_id "{name}_g"; transcript_id "{name}_t";\n'
        )
    (tmp_path / "orthodb.gtf").write_text(
        '1\tOrthoDB\texon\t31\t183\t.\t+\t.\tgene_id "p"; transcript_id "p1";\n'
    )
    return tmp_path


def test_generic_backbone_build_in_process(tmp_path, monkeypatch):
    inputs = _write_inputs(tmp_path)
    out = tmp_path / "out"
    calls = []
    real = builder.select_isoforms

    def counting(*a, **k):
        calls.append(1)
        return real(*a, **k)

    monkeypatch.setattr(builder, "select_isoforms", counting)
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "gmb-build",
            "--preset",
            "fungi",
            "--config",
            str(GMB_DIR / "tests" / "fixtures" / "test_no_pv.yaml"),
            "--genome",
            str(inputs / "genome.fa"),
            "--backbone",
            str(inputs / "backbone.gff3"),
            "--scallop",
            str(inputs / "scallop.gtf"),
            "--stringtie",
            str(inputs / "stringtie.gtf"),
            "--orthodb",
            str(inputs / "orthodb.gtf"),
            "--output-dir",
            str(out),
            "--no-log-file",
            "--gene-prefix",
            "T",
        ],
    )
    builder.main()

    resolved = yaml.safe_load((out / "resolved_config.yaml").read_text())
    assert resolved["scoring"]["backbone_label"] == "AUGUSTUS"
    attribution = pd.read_csv(out / "evidence_attribution.tsv", sep="\t")
    assert attribution["evidence_sources"].str.contains("AUGUSTUS").all()
    assert attribution["protein_alignment_sources"].eq("OrthoDB").all()
    # One locus, one selection pass (a second, discarded pass used to run).
    assert len(calls) == 1

    # evidence_attribution gene_id must be the GFF3 Parent of the same mRNA.
    parents = {}
    for line in (out / "consensus.gff3").read_text().splitlines():
        f = line.split("\t")
        if len(f) == 9 and f[2] == "mRNA":
            at = dict(kv.split("=", 1) for kv in f[8].split(";"))
            parents[at["ID"]] = at["Parent"]
    assert dict(zip(attribution["transcript_id"], attribution["gene_id"])) == parents
    assert list(attribution.columns) == EVIDENCE_ATTRIBUTION_COLUMNS
    assert os.path.exists(out / "run_manifest.json")


# ---------------------------------------------------------------------------
# Fixture build under each locus-clustering mode (same 500 kb fixture as
# test_z_tritici_subset; ~8 s each)
# ---------------------------------------------------------------------------


def _region_build(out_dir, overlay=None):
    import json

    cmd = [
        sys.executable,
        "-m",
        "gmb.cli.build",
        "--preset",
        "fungi",
        "--scallop",
        str(FIXTURE / "scallop_geneset.gtf"),
        "--stringtie",
        str(FIXTURE / "stringtie_geneset.gtf"),
        "--helixer",
        str(FIXTURE / "helixer_remapped.gff3"),
        "--orthodb",
        str(FIXTURE / "orthodb_geneset.gtf"),
        "--uniprot",
        str(FIXTURE / "uniprot_geneset.gtf"),
        "--genome",
        str(FIXTURE / "genome.fa"),
        "--output-dir",
        str(out_dir),
        "--gene-prefix",
        "ZT",
        "--no-log-file",
    ]
    if overlay:
        cmd += ["--config", str(overlay)]
    result = subprocess.run(cmd, capture_output=True, text=True, cwd=str(GMB_DIR))
    assert result.returncode == 0, result.stdout[-2000:] + result.stderr[-2000:]
    rows = [
        line.split("\t")
        for line in (out_dir / "consensus.gff3").read_text().splitlines()
        if line and not line.startswith("#")
    ]
    summary = json.loads((out_dir / "summary.json").read_text())["summary"]
    return rows, summary


def _assert_structurally_valid(rows):
    ids = [dict(kv.split("=", 1) for kv in r[8].split(";"))["ID"] for r in rows]
    assert len(ids) == len(set(ids)), "duplicate feature IDs"
    mrnas = {i for r, i in zip(rows, ids) if r[2] == "mRNA"}
    with_exon = {
        dict(kv.split("=", 1) for kv in r[8].split(";"))["Parent"] for r in rows if r[2] == "exon"
    }
    assert mrnas and mrnas <= with_exon, "mRNA without exons"


def _assert_gene_transcripts_overlap(rows):
    """No gene may hold transcripts that are not connected by overlaps."""
    by_gene = {}
    for r in rows:
        if r[2] == "mRNA":
            at = dict(kv.split("=", 1) for kv in r[8].split(";"))
            by_gene.setdefault(at["Parent"], []).append((int(r[3]), int(r[4])))
    for gene, spans in by_gene.items():
        spans.sort()
        reach = spans[0][1]
        for s, e in spans[1:]:
            assert s <= reach, f"{gene} has a transcript detached from the rest"
            reach = max(reach, e)


@pytest.mark.integration
def test_region_build_exon_overlap_reports_split_candidates(tmp_path):
    overlay = tmp_path / "exon.yaml"
    overlay.write_text("scoring:\n  locus_clustering: exon_overlap\n")
    rows, summary = _region_build(tmp_path / "exon", overlay)
    _assert_structurally_valid(rows)
    _assert_gene_transcripts_overlap(rows)
    assert summary["candidates_split_across_loci"] > 0


@pytest.mark.integration
def test_region_build_fungi_default(tmp_path):
    rows, summary = _region_build(tmp_path / "fungi")
    _assert_structurally_valid(rows)
    _assert_gene_transcripts_overlap(rows)
    assert summary["candidates_split_across_loci"] == 0
    assert summary["chimeras_long_span"] > 0
    assert sum(r[2] == "gene" for r in rows) > 100


# ---------------------------------------------------------------------------
# Preflight: uncollapsed long-read input
# ---------------------------------------------------------------------------


def _longread_preflight(tmp_path, n_models):
    from gmb.preflight import run_preflight

    (tmp_path / "g.fa").write_text(">1\n" + "ACGT" * 250 + "\n")  # 1 kb
    rows = "".join(
        f'1\tminimap\texon\t{10 + i}\t{60 + i}\t.\t+\t.\tgene_id "r{i}"; transcript_id "r{i}";\n'
        for i in range(n_models)
    )
    (tmp_path / "lr.gtf").write_text(rows)
    report = run_preflight(
        {
            "genome": str(tmp_path / "g.fa"),
            "tracks": [{"label": "Minimap2", "path": str(tmp_path / "lr.gtf")}],
        },
        load_config(None, "standard"),
    )
    return [c for c in report.checks if c.name == "longread_collapsed"], report


def test_preflight_warns_on_per_read_long_read_track(tmp_path):
    hits, report = _longread_preflight(tmp_path, 3)  # 3,000 models per Mb
    assert len(hits) == 1 and hits[0].verdict == "WARN"
    assert any("Collapse the long-read track" in r for r in report.recommendations)


def test_preflight_accepts_collapsed_long_read_track(tmp_path):
    hits, _ = _longread_preflight(tmp_path, 1)  # 1,000 models per Mb
    assert hits == []


def test_regenerated_fasta_minus_strand_multi_exon(tmp_path):
    """Ported from the removed fasta_export tests: the live exporter must splice and
    reverse-complement a minus-strand multi-exon CDS, keep the stop in cds.fa and drop it
    from prot.fa."""
    from gmb.pipeline.annotate_cds_utrs import reverse_complement, translate
    from gmb.pipeline.builder import regenerate_final_fasta

    coding = "ATG" + "GCT" * 10 + "AAA" * 9 + "TAA"  # 63 bp, M A.. K.. *
    plus = reverse_complement(coding)  # as it sits on the + strand
    seq = "C" * 20 + plus[:30] + "GT" + "A" * 26 + "AG" + plus[30:] + "C" * 20
    genome = {"1": seq}
    a, b = (20, 50), (50 + 30, 50 + 30 + 33)  # two CDS segments, genomic order
    rows = [
        {
            "Feature": "mRNA",
            "ID": "t1",
            "Chromosome": "1",
            "Strand": "-",
            "Start": a[0],
            "End": b[1],
        }
    ]
    for s_, e_ in (a, b):
        rows.append({"Feature": "exon", "Parent": "t1", "Start": s_, "End": e_})
        rows.append({"Feature": "CDS", "Parent": "t1", "Start": s_, "End": e_})
    regenerate_final_fasta(rows, genome, str(tmp_path))
    read = lambda n: "".join((tmp_path / n).read_text().splitlines()[1:])  # noqa: E731
    assert read("cds.fa") == coding
    assert read("cdna.fa") == coding
    assert read("prot.fa") == translate(coding)[:-1] and not read("prot.fa").endswith("*")


def test_list_presets_needs_no_output_dir():
    """Regression: --output-dir was argparse-required, so --list-presets always errored."""
    ok = subprocess.run(
        [sys.executable, "-m", "gmb.cli.build", "--list-presets"],
        capture_output=True,
        text=True,
        cwd=str(GMB_DIR),
    )
    assert ok.returncode == 0 and {"fungi", "apicomplexa"} <= set(ok.stdout.split())
    missing = subprocess.run(
        [sys.executable, "-m", "gmb.cli.build", "--genome", "g.fa"],
        capture_output=True,
        text=True,
        cwd=str(GMB_DIR),
    )
    assert missing.returncode == 2 and "--output-dir" in missing.stderr


# ---------------------------------------------------------------------------
# Preflight: malformed inputs
# ---------------------------------------------------------------------------


def _preflight_one(tmp_path, label, name, text):
    from gmb.preflight import run_preflight

    (tmp_path / "g.fa").write_text(">1\n" + "ACGT" * 250 + "\n")  # 1 kb
    (tmp_path / name).write_text(text)
    report = run_preflight(
        {
            "genome": str(tmp_path / "g.fa"),
            "tracks": [{"label": label, "path": str(tmp_path / name)}],
        },
        load_config(None, "fungi"),
    )
    return {c.name: c.verdict for c in report.checks if c.target == label}


GENE = "1\tH\tgene\t1\t{e}\t.\t+\t.\tID=g1\n1\tH\tmRNA\t1\t{e}\t.\t+\t.\tID=t1;Parent=g1\n"


@pytest.mark.parametrize(
    "label,name,text,check,verdict",
    [
        ("Helixer", "empty.gff3", "", "file_readable", "FAIL"),
        ("Helixer", "junk.gff3", "not\ta gff\n", "file_readable", "FAIL"),
        (
            "Scallop",
            "inv.gtf",
            '1\tS\texon\t200\t100\t.\t+\t.\ttranscript_id "a";\n',
            "file_readable",
            "FAIL",
        ),
        (
            "Helixer",
            "beyond.gff3",
            GENE.format(e=5000) + "1\tH\texon\t1\t5000\t.\t+\t.\tParent=t1\n",
            "coordinates_valid",
            "FAIL",
        ),
        (
            "Helixer",
            "cds_only.gff3",
            GENE.format(e=100) + "1\tH\tCDS\t10\t99\t.\t+\t0\tParent=t1\n",
            "exon_rows_present",
            "FAIL",
        ),
        (
            "Helixer",
            "no_cds.gff3",
            GENE.format(e=100) + "1\tH\texon\t1\t100\t.\t+\t.\tParent=t1\n",
            "backbone_cds_present",
            "WARN",
        ),
        (
            "OrthoDB",
            "prot_cds_only.gtf",
            '1\tO\tCDS\t10\t99\t.\t+\t0\ttranscript_id "p";\n',
            "exon_rows_present",
            None,
        ),
    ],
)
def test_preflight_rejects_malformed_inputs(tmp_path, label, name, text, check, verdict):
    assert _preflight_one(tmp_path, label, name, text).get(check) == verdict


# ---------------------------------------------------------------------------
# Public interface stability (2.0.0). Removing or renaming any of these is a
# breaking change: update the release notes, then this test.
# ---------------------------------------------------------------------------

EVIDENCE_ATTRIBUTION_COLUMNS = [
    "gene_id",
    "transcript_id",
    "evidence_sources",
    "protein_alignment_sources",
    "protein_support_strength",
    "structural_support_sources",
    "n_structural_support_sources",
    "backbone_shortread_agreement",
    "longread_structural_role",
    "backbone_intron_rescue",
    "selection_reason",
    "exon_count",
    "cds_bp",
    "utr_5p_bp",
    "utr_3p_bp",
    "max_exon_len_bp",
    "max_intron_len_bp",
    "transcript_span_bp",
    "gmb_score",
    "utr_5p_supported",
    "utr_3p_supported",
    "utr_5p_action",
    "utr_3p_action",
    "utr_5p_reason",
    "utr_3p_reason",
    # junction-level evidence, appended after 2.0.0
    "introns_without_transcript_support",
    "protein_alignments_compatible",
    "protein_alignments_incompatible",
]


def test_python_api_signature_is_stable():
    from gmb.api import GmbResult, run_gene_model_builder

    params = inspect.signature(run_gene_model_builder).parameters
    assert {k: v.default for k, v in params.items() if k != "genome"} == {
        "backbone": None,
        "backbone_kind": "helixer",
        "backbone_label": None,
        "short_read": None,
        "long_read": None,
        "protein_alignment": None,
        "preset": "standard",
        "config": None,
        "output_dir": "gmb_run",
        "gene_prefix": "GMB",
        "run_preflight": True,
        "stop_on_preflight_fail": True,
        "validate_fasta": True,
        "bin_dir": None,
    }
    import dataclasses

    assert {f.name for f in dataclasses.fields(GmbResult)} >= {
        "output_dir",
        "preflight_verdict",
        "qc_passed",
        "build_dir",
        "handover_dir",
        "annotation",
        "annotation_all_isoforms",
        "proteins",
        "cds",
        "cdna",
        "evidence_attribution",
        "resolved_config",
        "run_manifest",
        "handover_manifest",
        "stages_run",
        "returncodes",
    }


@pytest.mark.parametrize(
    "command,flags",
    [
        (
            "gmb.cli.build",
            {
                "--config",
                "--preset",
                "--list-presets",
                "--check-deps",
                "--genome",
                "--scallop",
                "--stringtie",
                "--minimap2",
                "--helixer",
                "--tiberius",
                "--backbone",
                "--backbone-label",
                "--orthodb",
                "--uniprot",
                "--genblast",
                "--output-dir",
                "--gene-prefix",
                "--log-file",
                "--no-log-file",
                "--validate-fasta",
                "--seqname",
                "--region",
                "--regions-file",
                "--sample-loci",
                "--assembly-report",
                "--seqname-map",
            },
        ),
        (
            "gmb.cli.preflight",
            {
                "--config",
                "--preset",
                "--genome",
                "--scallop",
                "--stringtie",
                "--minimap2",
                "--helixer",
                "--tiberius",
                "--backbone",
                "--backbone-label",
                "--orthodb",
                "--uniprot",
                "--genblast",
                "--output-dir",
                "--json",
                "--allow-fail",
                "--log-file",
                "--no-log-file",
            },
        ),
    ],
)
def test_cli_flags_are_stable(command, flags):
    out = subprocess.run(
        [sys.executable, "-m", command, "--help"], capture_output=True, text=True, cwd=str(GMB_DIR)
    ).stdout
    assert flags <= set(re.findall(r"--[a-z0-9-]+", out))


def _remap(tmp_path, gff_text, report_rows, *extra):
    (tmp_path / "in.gff3").write_text(gff_text)
    (tmp_path / "report.txt").write_text(
        "# Sequence-Name\tSequence-Role\tAssigned-Molecule\tAML\tGenBank-Accn\n"
        + "".join(
            f"chr{a}\tassembled-molecule\t{a}\tChromosome\t{acc}\n" for a, acc in report_rows
        )
    )
    return subprocess.run(
        [
            sys.executable,
            str(GMB_DIR / "tools" / "remap_helixer.py"),
            "--input",
            str(tmp_path / "in.gff3"),
            "--assembly-report",
            str(tmp_path / "report.txt"),
            "--output",
            str(tmp_path / "out.gff3"),
            *extra,
        ],
        capture_output=True,
        text=True,
    )


def test_remap_tool_renames_records_and_sequence_region_headers(tmp_path):
    res = _remap(
        tmp_path,
        "##gff-version 3\n##sequence-region CM1.1 1 100\n"
        "CM1.1\tHelixer\tgene\t1\t50\t.\t+\t.\tID=g\n",
        [("1", "CM1.1")],
    )
    assert res.returncode == 0, res.stderr
    out = (tmp_path / "out.gff3").read_text()
    assert "##sequence-region 1 1 100" in out and "\n1\tHelixer\tgene" in out


def test_remap_tool_fails_on_unmapped_sequences(tmp_path):
    gff = "CM2.1\tHelixer\tgene\t1\t50\t.\t+\t.\tID=g\n"
    assert _remap(tmp_path, gff, [("1", "CM1.1")]).returncode != 0
    assert _remap(tmp_path, gff, [("1", "CM1.1")], "--allow-unmapped").returncode == 0
