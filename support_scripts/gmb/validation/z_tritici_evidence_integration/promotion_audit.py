#!/usr/bin/env python3
"""Audit junction-supported promotions: independence of the promoted model's junction support
and reference verdict versus the backbone model it displaced.

Usage: promotion_audit.py <work> <arm> <windows.tsv>
"""
import contextlib, io, sys
from collections import defaultdict
import pandas as pd
from gmb.pipeline.builder import load_evidence

work, arm, wfile = sys.argv[1:4]
W = pd.read_csv(wfile, sep="\t", dtype=str)
TRACKS = {"Scallop": "scallop.gtf", "StringTie": "stringtie.gtf", "Minimap2": "minimap2_consensus.gtf"}

def introns(ivs):
    ivs = sorted(ivs); return [(ivs[i][1], ivs[i + 1][0]) for i in range(len(ivs) - 1)]

ref_cds = defaultdict(list); ref_str = {}
for l in open(f"{work}/w/ref.gff3"):
    f = l.rstrip("\n").split("\t")
    if len(f) == 9 and f[2] == "CDS":
        p = f[8].split("Parent=")[1].split(";")[0]; ref_cds[p].append((int(f[3]) - 1, int(f[4]))); ref_str[p] = (f[0], f[6])
by_key = defaultdict(list)
for t, v in ref_cds.items(): by_key[ref_str[t]].append(tuple(sorted(v)))

def rank(chrom, strand, c):
    best = 0
    for r in by_key[(chrom, strand)]:
        if r == c: return 3
        if len(c) > 1 and introns(r) == introns(c): best = max(best, 2)
        elif r[0][0] < c[-1][1] and c[0][0] < r[-1][1]: best = max(best, 1)
    return best

rows = []
for _, w in W.iterrows():
    d = f"{work}/runs/{arm}/{w['name']}/build"
    try:
        att = pd.read_csv(f"{d}/evidence_attribution.tsv", sep="\t")
    except FileNotFoundError:
        continue
    prom = att[att.selection_reason == "transcript_junction_support"]
    if prom.empty: continue
    obs = {}
    with contextlib.redirect_stdout(io.StringIO()):
        for src, fn in TRACKS.items():
            try:
                df, _ = load_evidence(f"{work}/win/{w['name']}/{fn}", src)
            except Exception:
                df = None
            s = set()
            if df is not None:
                for _t, g in df.groupby("transcript_id"):
                    s.update((str(g.Strand.iloc[0]), *i) for i in introns(zip(g.Start, g.End)))
            obs[src] = s
    exons, cds, strand = defaultdict(list), defaultdict(list), {}
    for l in open(f"{d}/consensus.gff3"):
        f = l.rstrip("\n").split("\t")
        if len(f) != 9: continue
        a = dict(kv.split("=", 1) for kv in f[8].split(";") if "=" in kv)
        if f[2] in ("exon", "CDS"):
            (exons if f[2] == "exon" else cds)[a["Parent"]].append((int(f[3]) - 1, int(f[4]))); strand[a["Parent"]] = f[6]
    for _, p in prom.iterrows():
        t = p.transcript_id; srcs = set(p.evidence_sources.split(","))
        ints = introns(exons[t])
        indep = sum(1 for i in ints if any((strand[t], *i) in obs[s] for s in obs if s not in srcs))
        bb = att[(att.gene_id == p.gene_id) & att.evidence_sources.str.contains("Helixer")]
        bb_rank = max((rank(w["chrom"], strand[b], tuple(sorted(cds[b]))) for b in bb.transcript_id), default=-1)
        rows.append(dict(set=w["set"], window=w["name"], cls=w["class"], tid=t, sources=p.evidence_sources,
                         n_sources=len(srcs), n_introns=len(ints), introns_confirmed_by_other_track=indep,
                         all_introns_independent=(len(ints) > 0 and indep == len(ints)),
                         backbone_kept=len(bb) > 0, new_rank=rank(w["chrom"], strand[t], tuple(sorted(cds[t]))),
                         backbone_rank=bb_rank))
df = pd.DataFrame(rows)
df["verdict"] = ["improved" if n > b else "regressed" if n < b else "uncertain" for n, b in zip(df.new_rank, df.backbone_rank)]
df["support"] = ["multi-source" if n > 1 else ("single-source, introns confirmed by another track" if a else "single-source, self-supported only")
                 for n, a in zip(df.n_sources, df.all_introns_independent)]
out = f"{work}/promotion_audit_{arm}.tsv"; df.to_csv(out, sep="\t", index=False)
print(arm, "promotions:", len(df), " backbone model retained in gene:", int(df.backbone_kept.sum()))
print(pd.crosstab([df.set, df.support], df.verdict, margins=True).to_string())
