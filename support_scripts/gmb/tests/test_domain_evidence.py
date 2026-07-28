"""Tests for gmb.pipeline.domain_evidence adapters and summarisation."""

from gmb.pipeline.domain_evidence import (
    parse_hmmer_domtblout,
    parse_interproscan_tsv,
    summarize_transcript_domain_metrics,
)


def test_parse_hmmer_domtblout(tmp_path):
    path = tmp_path / "domtblout.txt"
    # Realistic column layout (whitespace-separated, 22 columns); only the
    # columns this adapter reads are given meaningful values.
    header = "# target name accession tlen query name accession qlen E-value score bias # of c-Evalue i-Evalue score bias from to from to from to acc description\n"
    row = (
        "PF00069.28 PF00069.28 260 tx1 - 400 1.2e-40 130.5 0.3 1 1 "
        "3.4e-42 1.5e-38 129.8 0.3 5 255 10 260 8 262 0.95 Protein kinase domain\n"
    )
    path.write_text(header + row)

    rows = parse_hmmer_domtblout(str(path))
    assert len(rows) == 1
    r = rows[0]
    assert r["transcript_id"] == "tx1"
    assert r["feature_accession"] == "PF00069.28"
    assert r["start"] == 10
    assert r["end"] == 260
    assert r["coverage"] is not None


def test_parse_interproscan_tsv(tmp_path):
    path = tmp_path / "interpro.tsv"
    path.write_text(
        "tx1\tmd5\t300\tPfam\tPF00069\tProtein kinase domain\t10\t260\t1.2E-40\tT\t01-01-2024\tIPR000719\tKinase\n"
        "tx1\tmd5\t300\tPANTHER\tPTHR24347\tSerine/threonine kinase\t5\t295\t0.0\tT\t01-01-2024\tIPR000719\tKinase\n"
    )
    rows = parse_interproscan_tsv(str(path))
    assert len(rows) == 2
    assert {r["provider"] for r in rows} == {"Pfam", "PANTHER"}
    assert rows[0]["transcript_id"] == "tx1"
    assert rows[0]["start"] == 10
    assert rows[0]["end"] == 260


def test_summarize_transcript_domain_metrics_basic():
    rows = [
        {"transcript_id": "tx1", "provider": "Pfam", "feature_accession": "PF00069",
         "feature_name": "Kinase", "start": 10, "end": 260, "score": 130.0,
         "evalue": 1e-40, "coverage": 0.9, "status": None},
        {"transcript_id": "tx1", "provider": "PANTHER", "feature_accession": "PTHR24347",
         "feature_name": "Kinase-like", "start": 5, "end": 295, "score": None,
         "evalue": None, "coverage": None, "status": None},
    ]
    summary = summarize_transcript_domain_metrics(rows, protein_lengths={"tx1": 300})
    s = summary["tx1"]
    assert s["n_significant_domains"] == 2  # second has no evalue -> always counted
    assert s["cross_provider_agreement"] is True
    assert s["n_nonoverlapping_domains"] >= 1
    assert s["protein_coverage_fraction"] is not None


def test_summarize_no_domains_is_not_a_penalty_signal():
    # No rows at all for this transcript -- summary dict simply absent,
    # canonical_selection.py must treat that as "no support", not "penalised".
    summary = summarize_transcript_domain_metrics([], protein_lengths={"tx1": 300})
    assert summary == {}


def test_duplicate_domain_detected():
    rows = [
        {"transcript_id": "tx1", "provider": "Pfam", "feature_accession": "PF00069",
         "feature_name": "Kinase", "start": 10, "end": 100, "score": 50.0,
         "evalue": 1e-10, "coverage": None, "status": None},
        {"transcript_id": "tx1", "provider": "Pfam", "feature_accession": "PF00069",
         "feature_name": "Kinase", "start": 150, "end": 240, "score": 48.0,
         "evalue": 1e-9, "coverage": None, "status": None},
    ]
    summary = summarize_transcript_domain_metrics(rows)
    assert summary["tx1"]["has_duplicate_domain"] is True
