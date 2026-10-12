#!/usr/bin/env python3
"""Genes whose canonical transcript is not the build primary (<gene>.t1), with the first ranking
tier separating the two and a reference verdict. Usage: canonical_overrides.py <work> <windows.tsv> <arm> ...
"""
import sys
from collections import defaultdict
import pandas as pd
work, wfile, *arms = sys.argv[1:]
W = pd.read_csv(wfile, sep="\t", dtype=str)
TIERS = [("has_complete_orf", True, "complete ORF"), ("has_internal_stop", False, "no internal stop"),
         ("protein_validation_subtotal", True, "protein validation"), ("n_evidence_classes", True, "evidence-class breadth"),
         ("n_independent_sources", True, "named-source count"), ("gmb_score", True, "GMB score"), ("cds_bp", True, "CDS length")]
ref = defaultdict(list); rs = {}
for l in open(f"{work}/w/ref.gff3"):
    f = l.rstrip("\n").split("\t")
    if len(f) == 9 and f[2] == "CDS":
        p = f[8].split("Parent=")[1].split(";")[0]; ref[p].append((int(f[3]), int(f[4]))); rs[p] = (f[0], f[6])
refsig = {(rs[t][0], rs[t][1], tuple(sorted(v))) for t, v in ref.items()}
for arm in arms:
    rows = []
    for _, w in W.iterrows():
        d = f"{work}/runs/{arm}/{w['name']}/fin"
        try:
            can = pd.read_csv(f"{d}/canonical/canonical_transcripts.tsv", sep="\t")
            rk = pd.read_csv(f"{d}/canonical/transcript_ranking.tsv", sep="\t").set_index("transcript_id")
        except FileNotFoundError:
            continue
        cds = defaultdict(list); st = {}
        for l in open(f"{d}/consensus.gff3"):
            f = l.rstrip("\n").split("\t")
            if len(f) == 9 and f[2] == "CDS":
                p = f[8].split("Parent=")[1].split(";")[0]; cds[p].append((int(f[3]), int(f[4]))); st[p] = (f[0], f[6])
        ex = lambda t: (st[t][0], st[t][1], tuple(sorted(cds[t]))) in refsig if t in st else False
        for _, g in can.iterrows():
            p = f"{g.gene_id}.t1"; c = g.canonical_transcript_id
            if c == p or p not in rk.index: continue
            a, b = rk.loc[c], rk.loc[p]
            tier = "lexical"
            for col, higher, name in TIERS:
                va, vb = a[col], b[col]
                if pd.isna(va) and pd.isna(vb): continue
                if va != vb: tier = name; break
            rows.append(dict(set=w["set"], window=w["name"], gene=g.gene_id, canonical=c, tier=tier,
                             canonical_sources=a.evidence_sources, primary_sources=b.evidence_sources,
                             verdict=("gain" if ex(c) and not ex(p) else "loss" if ex(p) and not ex(c) else "neither/both")))
    df = pd.DataFrame(rows)
    n_genes = sum(len(pd.read_csv(f"{work}/runs/{arm}/{w}/fin/canonical/canonical_transcripts.tsv", sep="\t"))
                  for w in W["name"] if __import__("os").path.exists(f"{work}/runs/{arm}/{w}/fin/canonical/canonical_transcripts.tsv"))
    print(f"== {arm}: {len(df)} canonical != build primary of {n_genes} genes")
    if len(df):
        print(pd.crosstab([df.set, df.tier], df.verdict, margins=True).to_string())
        df.to_csv(f"{work}/canonical_overrides_{arm}.tsv", sep="\t", index=False)
