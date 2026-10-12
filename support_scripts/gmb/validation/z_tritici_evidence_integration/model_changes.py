#!/usr/bin/env python3
"""Classify each canonical transcript against Helixer's models, with a reference verdict.

unchanged   : CDS identical to a Helixer model's CDS
changed     : differs from every Helixer CDS but overlaps one on the same strand
              -> verdict vs the overlapped Helixer model: improved / regressed / uncertain
new_locus   : no same-strand Helixer CDS overlap
Verdict ranks: ref CDS exact (3) > ref CDS intron chain (2) > overlaps ref CDS (1) > none (0).
Usage: model_changes.py <runs> <win root> <ref.gff3> <wins csv> <arm> ...
"""
import sys
from collections import defaultdict, Counter
runs, wroot, ref, wins, *arms = sys.argv[1:]

def read(path, feat="CDS"):
    cds, strand = defaultdict(list), {}
    par = {}
    for l in open(path):
        f = l.rstrip("\n").split("\t")
        if len(f) != 9: continue
        a = dict(kv.split("=", 1) for kv in f[8].split(";") if "=" in kv)
        if f[2] == "CDS":
            cds[a["Parent"]].append((int(f[3]), int(f[4]))); strand[a["Parent"]] = (f[0], f[6])
        elif f[2] == "mRNA":
            par[a["ID"]] = a.get("Parent"); 
            if a.get("Ensembl_canonical") or "Ensembl_canonical" in f[8]: par[a["ID"] + "#canon"] = True
    return {t: tuple(sorted(v)) for t, v in cds.items()}, strand, par

rc, rs, _ = read(ref)
ref_by_strand = defaultdict(list)
for t, c in rc.items(): ref_by_strand[rs[t]].append(c)
def chain(c): return tuple((c[i][1], c[i + 1][0]) for i in range(len(c) - 1))
def rank(c, key):
    best = 0
    for r in ref_by_strand[key]:
        if r == c: return 3
        if len(c) > 1 and chain(r) == chain(c): best = max(best, 2)
        elif r[0][0] <= c[-1][1] and c[0][0] <= r[-1][1]: best = max(best, 1)
    return best
def ov(a, b): return any(min(e1, e2) >= max(s1, s2) for s1, e1 in a for s2, e2 in b)

for arm in arms:
    tot = Counter()
    for w in wins.split(","):
        hc, hs, _ = read(f"{wroot}/{w}/helixer.gff3")
        hx = defaultdict(list)
        for t, c in hc.items(): hx[hs[t]].append(c)
        gff = f"{runs}/{arm}/{w}/fin/canonical/consensus.canonical_annotated.gff3"
        cc, cs, par = read(gff)
        canon = [t for t in cc if par.get(t + "#canon")]
        for t in canon:
            c, key = cc[t], cs[t]
            if c in set(hx[key]): tot["unchanged"] += 1; continue
            hits = [h for h in hx[key] if ov(h, c)]
            if not hits: tot["new_locus"] += 1; tot[f"new_locus_rank{rank(c, key)}"] += 1; continue
            tot["changed"] += 1
            rn = rank(c, key); ro = max(rank(h, key) for h in hits)
            tot["changed_improved" if rn > ro else "changed_regressed" if rn < ro else "changed_uncertain"] += 1
    print(arm, dict(sorted(tot.items())))
