#!/usr/bin/env python3
"""Genes whose isoforms have CDSs that do not overlap the primary's CDS. Usage: <runs> <arm> <wins csv> <ref.gff3>"""
import sys
from collections import defaultdict, Counter
runs, arm, wins, ref = sys.argv[1:5]
rc = defaultdict(list); rs = {}
for l in open(ref):
    f = l.rstrip("\n").split("\t")
    if len(f) == 9 and f[2] == "CDS":
        p = f[8].split("Parent=")[1].split(";")[0]; rc[p].append((int(f[3]), int(f[4]))); rs[p] = (f[0], f[6])
refsig = {(rs[t][0], rs[t][1], tuple(sorted(v))) for t, v in rc.items()}
tot = Counter()
for w in wins.split(","):
    tx = {}; cds = defaultdict(list); genes = defaultdict(list)
    for l in open(f"{runs}/{arm}/{w}/build/consensus.gff3"):
        f = l.rstrip("\n").split("\t")
        if len(f) != 9: continue
        a = dict(kv.split("=", 1) for kv in f[8].split(";") if "=" in kv)
        if f[2] == "mRNA": tx[a["ID"]] = (f[0], f[6], a["Parent"], a.get("Evidence", "")); genes[a["Parent"]].append(a["ID"])
        elif f[2] == "CDS": cds[a["Parent"]].append((int(f[3]), int(f[4])))
    for g, ts in genes.items():
        tot["genes"] += 1
        p = f"{g}.t1"
        if p not in tx or not cds.get(p): continue
        pc = cds[p]
        dis = [t for t in ts if t != p and cds.get(t) and not any(min(e1, e2) >= max(s1, s2) for s1, e1 in pc for s2, e2 in cds[t])]
        if not dis: continue
        tot["genes_with_cds_disjoint_isoform"] += 1
        p_bb = "Helixer" in tx[p][3]
        tot[f"primary_backbone={p_bb}"] += 1
        p_ex = (tx[p][0], tx[p][1], tuple(sorted(pc))) in refsig
        d_ex = any((tx[t][0], tx[t][1], tuple(sorted(cds[t]))) in refsig for t in dis)
        tot[f"primary_bb={p_bb} primary_exact={p_ex} disjoint_alt_exact={d_ex}"] += 1
for k, v in sorted(tot.items()): print(f"{k:60s} {v}")
