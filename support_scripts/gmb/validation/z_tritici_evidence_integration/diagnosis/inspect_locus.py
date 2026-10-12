#!/usr/bin/env python3
"""Show a reference gene's CDS and the overlapping GMB transcripts in several arms.

Usage: inspect_locus.py <runs> <ref.gff3> <window> <ref_gene_id> <arm> [<arm> ...]
"""
import sys
from collections import defaultdict

import pandas as pd

runs, ref, win, gene, *arms = sys.argv[1:]


def parse(path):
    tx, cds = {}, defaultdict(list)
    for line in open(path):
        f = line.rstrip("\n").split("\t")
        if len(f) != 9:
            continue
        a = dict(kv.split("=", 1) for kv in f[8].split(";") if "=" in kv)
        if f[2] == "mRNA":
            tx[a["ID"]] = (f[0], int(f[3]), int(f[4]), f[6], a.get("Parent"), a.get("Evidence", ""))
        elif f[2] == "CDS":
            cds[a["Parent"]].append((int(f[3]), int(f[4])))
    return tx, {k: sorted(v) for k, v in cds.items()}


rtx, rcds = parse(ref)
rid = next(t for t, v in rtx.items() if v[4] in (gene, f"gene:{gene}"))
chrom, s, e, strand, *_ = rtx[rid]
print(f"REF {gene} {chrom}:{s}-{e} {strand}  CDS {rcds.get(rid)}")
for arm in arms:
    d = f"{runs}/{arm}/{win}/build"
    tx, cds = parse(f"{d}/consensus.gff3")
    att = pd.read_csv(f"{d}/evidence_attribution.tsv", sep="\t").set_index("transcript_id")
    print(f"-- {arm}")
    for t, (c, ts, te, st, par, ev) in sorted(tx.items(), key=lambda kv: kv[1][1]):
        if c != chrom or te < s or ts > e or st != strand:
            continue
        a = att.loc[t] if t in att.index else {}
        mark = "EXACT" if cds.get(t) == rcds.get(rid) else ""
        print(f"  {t:14s} {ev:28s} score={a.get('gmb_score')} reason={a.get('selection_reason')} "
              f"unsup_introns={a.get('introns_without_transcript_support')} "
              f"prot={a.get('protein_alignments_compatible')}/{a.get('protein_alignments_incompatible')} "
              f"{mark}\n      CDS {cds.get(t)}")
