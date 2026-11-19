"""End-to-end build: a coding unit of a filter-removed read-through transcript is recovered.

The synthetic locus of ``test_secondary_orfs_e2e`` (a StringTie read-through carrying a long ORF
and, downstream, a second spliced ORF whose Helixer model uses an unobserved alternative acceptor)
is built with ``transcriptomic_filter.max_transcript_length`` below the transcript's span, so GMB's
filter removes the transcript. With secondary ORFs alone nothing can use it and the Helixer model
stays primary. With ``filtered_transcript_orfs`` the downstream unit is admitted by the frozen
A_ev rule (its intron is observed in five construction groups) and replaces the Helixer primary
of the existing gene. Writing the switch as false is byte-identical to leaving it out.
"""

import csv
import json
import subprocess
import sys
from pathlib import Path

import pytest

pytestmark = pytest.mark.integration
_GMB_DIR = Path(__file__).parent.parent

F = "GCC"
ORF1 = "ATG" + F * 80 + "TAA"
ORF2_A, ORF2_B = "ATG" + F * 30, F * 30 + "TAA"
INTRON = "GT" + "C" * 53 + "AG" + "C" + "AG"
LEFT = F * 20
E1 = F * 10 + ORF1 + F * 5 + ORF2_A
E2 = ORF2_B + F * 10
SEQ = LEFT + E1 + INTRON + E2 + F * 200
I0, I1 = len(LEFT) + len(E1), len(LEFT) + len(E1) + len(INTRON)
TX = [(len(LEFT), I0), (I1, I1 + len(E2))]  # span 567 bp
ORF2_CDS = [(I0 - len(ORF2_A), I0), (I1, I1 + len(ORF2_B))]  # span 246 bp
HX_CDS = [ORF2_CDS[0], (I1 - 3, I1 + len(ORF2_B))]
MAX_SPAN = 400  # removes the transcript, keeps the unit


def _write(d: Path):
    (d / "genome.fa").write_text(
        ">1\n" + "\n".join(SEQ[i : i + 60] for i in range(0, len(SEQ), 60)) + "\n"
    )
    hx = [
        "##gff-version 3",
        f"1\tHelixer\tgene\t{HX_CDS[0][0] + 1}\t{HX_CDS[1][1]}\t.\t+\t.\tID=hxg",
        f"1\tHelixer\tmRNA\t{HX_CDS[0][0] + 1}\t{HX_CDS[1][1]}\t.\t+\t.\tID=hxg.1;Parent=hxg",
    ]
    phase = 0
    for i, (a, b) in enumerate(HX_CDS):
        hx.append(f"1\tHelixer\texon\t{a + 1}\t{b}\t.\t+\t.\tID=hxg.1.exon.{i + 1};Parent=hxg.1")
        hx.append(
            f"1\tHelixer\tCDS\t{a + 1}\t{b}\t.\t+\t{phase}\tID=hxg.1.CDS.{i + 1};Parent=hxg.1"
        )
        phase = (3 - ((b - a) - phase) % 3) % 3
    (d / "helixer.gff3").write_text("\n".join(hx) + "\n")
    st = [
        f'1\tStringTie\ttranscript\t{TX[0][0] + 1}\t{TX[1][1]}\t1000\t+\t.\tgene_id "S1"; transcript_id "S1.1";'
    ]
    st += [
        f'1\tStringTie\texon\t{a + 1}\t{b}\t1000\t+\t.\tgene_id "S1"; transcript_id "S1.1"; exon_number "{i + 1}";'
        for i, (a, b) in enumerate(TX)
    ]
    (d / "stringtie.gtf").write_text("\n".join(st) + "\n")
    lines = [
        "library_id\tpath\trole\tstudy\tstrain\tcondition\tcount_unit\tmax_intron\tmin_intron"
        "\tnoncanonical_filtered"
    ]
    for g in range(5):
        (d / f"L{g}.sj.tab").write_text(f"1\t{I0 + 1}\t{I1}\t1\t1\t1\t60\t0\t30\n")
        lines.append(f"L{g}\tL{g}.sj.tab\tconstruction\tS{g}\tX\tc{g}\tmate\t1000\t21\ttrue")
    (d / "manifest.tsv").write_text("\n".join(lines) + "\n")
    head = (
        f"transcriptomic_filter:\n  max_transcript_length: {MAX_SPAN}\n"
        "scoring:\n  primary_selection: evidence_score\n  evidence_score:\n"
        "    secondary_orfs: true\n"
    )
    (d / "sorf.yaml").write_text(head)
    (d / "on.yaml").write_text(head + "    filtered_transcript_orfs: true\n")
    (d / "off.yaml").write_text(head + "    filtered_transcript_orfs: false\n")


def _build(d: Path, cfg: str, out: str):
    r = subprocess.run(
        [
            sys.executable,
            "-m",
            "gmb.cli.build",
            "--genome",
            str(d / "genome.fa"),
            "--helixer",
            str(d / "helixer.gff3"),
            "--stringtie",
            str(d / "stringtie.gtf"),
            "--preset",
            "fungi_latest_settings",
            "--config",
            str(d / cfg),
            "--star-junctions",
            str(d / "manifest.tsv"),
            "--gene-prefix",
            "E2E",
            "--no-log-file",
            "--output-dir",
            str(d / out),
        ],
        capture_output=True,
        text=True,
        cwd=str(_GMB_DIR),
    )
    assert r.returncode == 0, r.stdout[-3000:] + r.stderr[-3000:]
    return d / out


@pytest.fixture(scope="module")
def world(tmp_path_factory):
    d = tmp_path_factory.mktemp("ft_e2e")
    _write(d)
    return {k: _build(d, f"{k}.yaml", k) for k in ("sorf", "on", "off")}


def _rows(out):
    rows = {}
    for line in (out / "consensus.gff3").read_text().splitlines():
        p = line.split("\t")
        if len(p) == 9:
            rows.setdefault(p[2], []).append(p)
    return rows


def _cds(rows, tid):
    return sorted((int(p[3]) - 1, int(p[4])) for p in rows["CDS"] if f"Parent={tid}" in p[8])


def test_without_the_option_the_removed_transcript_is_lost(world):
    rows = _rows(world["sorf"])
    assert not any("CandidateOrigin" in p[8] for p in rows["mRNA"])
    assert _cds(rows, "E2E_00001.t1") == HX_CDS


def test_unit_of_the_removed_transcript_replaces_the_contradicted_model(world):
    on, rows = world["on"], _rows(world["on"])
    assert len(rows["gene"]) == len(_rows(world["sorf"])["gene"])  # no new gene
    gen = [p for p in rows["mRNA"] if "CandidateOrigin=filtered_transcript_orf" in p[8]]
    assert len(gen) == 1 and "SecondaryOrfParents=StringTie_S1.1" in gen[0][8]
    tid = gen[0][8].split("ID=")[1].split(";")[0]
    assert tid.endswith(".t1") and _cds(rows, tid) == ORF2_CDS
    sel = [
        r
        for r in csv.DictReader(open(on / "evidence_selection.tsv"), delimiter="\t")
        if r["outcome"] == "replaced"
    ]
    assert len(sel) == 1 and sel[0]["new_candidate_origin"] == "filtered_transcript_orf"
    so = {
        r["candidate_id"]: r
        for r in csv.DictReader(open(on / "secondary_orfs.tsv"), delimiter="\t")
    }
    unit = so["StringTie_S1.1.ftorf1"]
    assert unit["status"] == "admitted" and unit["admitted_by"] == "junctions"
    assert unit["relations"] == "removed_downstream" and unit["selected"] == "True"
    assert so["StringTie_S1.1.ftown1"]["status"] == "no_shared_coding_locus"  # never a new gene
    summ = json.load(open(on / "secondary_orfs_summary.json"))["filtered_transcript_orfs"]
    assert summ["selected_as_primary"] == 1 and summ["parent_filter_reasons"] == {"span": 2}


def test_switch_written_false_is_byte_identical(world):
    for f in (
        "consensus.gff3",
        "cds.fa",
        "prot.fa",
        "evidence_selection.tsv",
        "secondary_orfs.tsv",
        "resolved_config.yaml",
        "evidence_attribution.tsv",
    ):
        assert (world["off"] / f).read_bytes() == (world["sorf"] / f).read_bytes(), f
    a = json.load(open(world["off"] / "secondary_orfs_summary.json"))
    assert "filtered_transcript_orfs" not in a
