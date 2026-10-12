#!/usr/bin/env python3
"""In-frame protein coverage of new-locus (no backbone transcript) genes in a GMB build (dev)."""
import contextlib, io, sys
from collections import defaultdict
import pandas as pd
from gmb.pipeline.builder import load_evidence
runs, wroot, arm, wins = sys.argv[1:5]

def frames(ivs, strand):
    """{genomic_pos: codon phase} for CDS-like intervals (phase 0 at the 5' end)."""
    out = {}; cum = 0
    order = sorted(ivs) if strand == "+" else sorted(ivs, reverse=True)
    for s, e in order:
        rng = range(s, e) if strand == "+" else range(e - 1, s - 1, -1)
        for p in rng:
            out[p] = cum % 3; cum += 1
    return out

rows = []
for w in wins.split(","):
    with contextlib.redirect_stdout(io.StringIO()):
        od, _ = load_evidence(f"{wroot}/{w}/orthodb.gtf", "OrthoDB"); up, _ = load_evidence(f"{wroot}/{w}/genblast.gtf", "UniProt")
    pr = pd.concat([od, up])
    prot = defaultdict(list)
    for t, g in pr.groupby("transcript_id"):
        s = str(g.Strand.iloc[0]); ivs = sorted(zip(g.Start.astype(int), g.End.astype(int)))
        prot[s].append((ivs[0][0], ivs[-1][1], ivs))
    lab = pd.read_csv(f"{runs}/{arm}/{w}/cmp_primary/consensus_transcript_labels.tsv", sep="\t").set_index("transcript_id")
    tx, cds, gene_src = {}, defaultdict(list), defaultdict(set)
    for line in open(f"{runs}/{arm}/{w}/build/consensus.gff3"):
        f = line.rstrip("\n").split("\t")
        if len(f) != 9: continue
        a = dict(kv.split("=", 1) for kv in f[8].split(";") if "=" in kv)
        if f[2] == "mRNA":
            tx[a["ID"]] = (f[6], a["Parent"]); gene_src[a["Parent"]].update(a.get("Evidence", "").split(","))
        elif f[2] == "CDS":
            cds[a["Parent"]].append((int(f[3]) - 1, int(f[4])))
    for t, (strand, gid) in tx.items():
        if not t.endswith(".t1") or "Helixer" in gene_src[gid] or t not in lab.index or not cds.get(t): continue
        c = sorted(cds[t]); cf = frames(c, strand); lo, hi = c[0][0], c[-1][1]
        cint = {(c[i][1], c[i + 1][0]) for i in range(len(c) - 1)}
        cov = set()
        for ps, pe, ivs in prot[strand]:
            if pe <= lo or ps >= hi: continue
            a, b = max(lo, ps), min(hi, pe)
            pint = {(ivs[i][1], ivs[i + 1][0]) for i in range(len(ivs) - 1)}
            if {i for i in pint if i[0] >= a and i[1] <= b} != {i for i in cint if i[0] >= a and i[1] <= b}: continue
            pf = frames(ivs, strand)
            cov.update(p for p, ph in pf.items() if cf.get(p) == ph)
        r = lab.loc[t]
        grp = r.novel_category if isinstance(r.novel_category, str) else r.classification
        grp = {"Overlaps_evaluated_reference_opposite_strand": "antisense", "Strand_Mismatch": "antisense",
               "Novel_no_reference_overlap": "novel_intergenic", "Matched": "matched_sense"}.get(grp, grp)
        rows.append(dict(w=w, tid=t, grp=grp, cds_bp=sum(e - s for s, e in c), inframe_frac=len(cov) / sum(e - s for s, e in c),
                         cds_exact=r.best_cds_overlap == 1.0))
d = pd.DataFrame(rows); d.to_csv(f"{runs}/../newlocus_inframe_{arm}.tsv", sep="\t", index=False)
d["bin"] = pd.cut(d.inframe_frac, [-0.01, 0, 0.25, 0.5, 0.8, 1.0], labels=["0", "(0,.25]", "(.25,.5]", "(.5,.8]", "(.8,1]"])
print(pd.crosstab(d.grp, d.bin, margins=True))
