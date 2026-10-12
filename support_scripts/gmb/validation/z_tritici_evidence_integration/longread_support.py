#!/usr/bin/env python3
"""Independent check: fraction of canonical CDS introns confirmed by the collapsed long-read track,
for arms that did NOT use long reads (reference-free). Usage: longread_support.py <work> <wins csv> <arm> ..."""
import sys
from collections import defaultdict
work, wins, *arms = sys.argv[1:]

def introns_by_tx(path, feature, canonical_only=False):
    ivs, strand, canon = defaultdict(list), {}, set()
    for l in open(path):
        f = l.rstrip("\n").split("\t")
        if len(f) != 9:
            continue
        a = dict(kv.split("=", 1) for kv in f[8].split(";") if "=" in kv)
        if f[2] == "mRNA" and "Ensembl_canonical" in a.get("tag", ""):
            canon.add(a["ID"])
        key = a.get("Parent") or a.get("transcript_id")
        if f[2] == feature:
            if key is None:
                key = f[8].split('transcript_id "')[1].split('"')[0]
            ivs[key].append((int(f[3]), int(f[4]))); strand[key] = (f[0], f[6])
    out = {}
    for t, v in ivs.items():
        if canonical_only and t not in canon:
            continue
        v.sort(); c, s = strand[t]
        out[t] = [(c, s, v[i][1], v[i + 1][0]) for i in range(len(v) - 1)]
    return out

for arm in arms:
    tot = conf = 0
    for w in wins.split(","):
        lr = set()
        lrx = introns_by_tx(f"{work}/win/{w}/minimap2_consensus.gtf", "exon")
        for v in lrx.values():
            lr.update(v)
        if arm == "A_helixer":
            q = introns_by_tx(f"{work}/win/{w}/helixer.gff3", "CDS")
        else:
            q = introns_by_tx(f"{work}/runs/{arm}/{w}/fin/canonical/consensus.canonical_annotated.gff3", "CDS", True)
        for v in q.values():
            tot += len(v); conf += sum(i in lr for i in v)
    print(f"{arm:12s} canonical CDS introns {tot:6d}  confirmed by long reads {conf:6d} ({100 * conf / tot:.1f}%)")
