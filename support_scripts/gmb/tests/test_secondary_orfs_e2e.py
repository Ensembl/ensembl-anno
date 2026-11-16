"""End-to-end build on a synthetic read-through locus where a secondary ORF must be selected.

One StringTie read-through transcript carries a long ORF (its GMB ORF) and, downstream, a
second spliced ORF. The Helixer model of the second gene uses an alternative acceptor 3 bp
upstream that no junction observes, while the transcript's intron is observed in five
construction groups. With ``secondary_orfs`` on, the second ORF must replace the Helixer
primary of that existing gene as a CDS-only transcript; with it off, nothing changes.
"""

import csv
import json
import subprocess
import sys
from pathlib import Path

import pytest

pytestmark = pytest.mark.integration
_GMB_DIR = Path(__file__).parent.parent

F = "GCC"  # no ATG or stop in any frame
ORF1 = "ATG" + F * 80 + "TAA"
ORF2_A, ORF2_B = "ATG" + F * 30, F * 30 + "TAA"
# Intron with a second canonical acceptor (AG) 3 bp upstream of the real one.
INTRON = "GT" + "C" * 53 + "AG" + "C" + "AG"
LEFT = F * 20
E1 = F * 10 + ORF1 + F * 5 + ORF2_A
E2 = ORF2_B + F * 10
SEQ = LEFT + E1 + INTRON + E2 + F * 200
I0, I1 = len(LEFT) + len(E1), len(LEFT) + len(E1) + len(INTRON)  # 0-based half-open intron
TX = [(len(LEFT), I0), (I1, I1 + len(E2))]
ORF2_CDS = [(I0 - len(ORF2_A), I0), (I1, I1 + len(ORF2_B))]
HX_CDS = [ORF2_CDS[0], (I1 - 3, I1 + len(ORF2_B))]  # alternative acceptor, in frame


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
    for name, on in (("on", "true"), ("off", "false")):
        (d / f"{name}.yaml").write_text(
            "scoring:\n  primary_selection: evidence_score\n"
            f"  evidence_score:\n    secondary_orfs: {on}\n"
        )


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
    d = tmp_path_factory.mktemp("sorf_e2e")
    _write(d)
    return {"on": _build(d, "on.yaml", "on"), "off": _build(d, "off.yaml", "off")}


def _mrnas(out):
    rows = {}
    for line in (out / "consensus.gff3").read_text().splitlines():
        p = line.split("\t")
        if len(p) == 9:
            rows.setdefault(p[2], []).append(p)
    return rows


def test_off_keeps_the_helixer_primary(world):
    rows = _mrnas(world["off"])
    assert not any("CandidateOrigin" in p[8] for p in rows["mRNA"])
    assert not (world["off"] / "secondary_orfs.tsv").exists()


def test_secondary_orf_replaces_contradicted_helixer_model(world):
    out = world["on"]
    rows = _mrnas(out)
    assert len(rows["gene"]) == len(_mrnas(world["off"])["gene"])  # no new gene
    gen = [p for p in rows["mRNA"] if "CandidateOrigin=secondary_orf" in p[8]]
    assert (
        len(gen) == 1
        and "Evidence=StringTie" in gen[0][8]
        and "SecondaryOrfParents=StringTie_S1.1" in gen[0][8]
    )
    tid = gen[0][8].split("ID=")[1].split(";")[0]
    assert tid.endswith(".t1")  # it is the primary
    cds = sorted((int(p[3]) - 1, int(p[4])) for p in rows["CDS"] if f"Parent={tid}" in p[8])
    exons = sorted((int(p[3]) - 1, int(p[4])) for p in rows["exon"] if f"Parent={tid}" in p[8])
    assert cds == ORF2_CDS and exons == ORF2_CDS  # CDS-only transcript
    phases = [
        p[7] for p in sorted(rows["CDS"], key=lambda p: int(p[3])) if f"Parent={tid}" in p[8]
    ]
    assert phases == ["0", str((3 - len(ORF2_A) % 3) % 3)]
    prot = (out / "prot.fa").read_text()
    assert f">{tid}\n" in prot and "M" + "A" * 60 in prot.replace("\n", "")
    sel = list(csv.DictReader(open(out / "evidence_selection.tsv"), delimiter="\t"))
    rec = [r for r in sel if r["outcome"] == "replaced"]
    assert len(rec) == 1 and rec[0]["new_candidate_origin"] == "secondary_orf"
    assert (
        rec[0]["representative"] == "secondary_orf"
        and rec[0]["new_candidate"] == "StringTie_S1.1.sorf1"
    )
    so = list(csv.DictReader(open(out / "secondary_orfs.tsv"), delimiter="\t"))
    chosen = [r for r in so if r["selected"] == "True"]
    assert len(chosen) == 1 and chosen[0]["admitted_by"] == "junctions"
    assert (
        chosen[0]["relations"] == "downstream" and chosen[0]["selected_gene"] == rec[0]["gene_id"]
    )
    # the Helixer model's 3 bp in-frame acceptor extension is the abandoned sequence; no other
    # gene covers it (reported for evaluation only)
    assert float(chosen[0]["abandoned_old_cds_bp"]) == 3.0
    assert float(chosen[0]["abandoned_bp_in_other_genes"]) == 0.0
    summary = json.loads((out / "secondary_orfs_summary.json").read_text())
    assert summary["selected_as_primary"] == 1
