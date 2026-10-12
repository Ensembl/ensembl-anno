#!/usr/bin/env python3
"""Potential of protein-alignment-defined gene candidates in gaps of the GMB output (dev only)."""
import contextlib, io, os, sys
from collections import defaultdict, Counter
import pandas as pd
from gmb.pipeline.annotate_cds_utrs import load_genome, reverse_complement, translate
from gmb.pipeline.builder import load_evidence
runs, wroot, ref, wins = sys.argv[1:5]
W = {l.split("\t")[0]: l.split("\t") for l in open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "windows.tsv")) if not l.startswith("name")}
rcds = defaultdict(list); rstr = {}
for line in open(ref):
    f = line.rstrip("\n").split("\t")
    if len(f) == 9 and f[2] == "CDS":
        p = f[8].split("Parent=")[1].split(";")[0]; rcds[p].append((int(f[3]) - 1, int(f[4]))); rstr[p] = (f[0], f[6])
refsig = {(rstr[t][0], rstr[t][1], tuple(sorted(v))) for t, v in rcds.items()}
refspan = defaultdict(list)
for t, v in rcds.items(): refspan[rstr[t]].append((min(a for a, _ in v), max(b for _, b in v)))
tot = Counter()
for w in wins.split(","):
    chrom = W[w][2]
    with contextlib.redirect_stdout(io.StringIO()):
        od, _ = load_evidence(f"{wroot}/{w}/orthodb.gtf", "OrthoDB"); up, _ = load_evidence(f"{wroot}/{w}/genblast.gtf", "UniProt")
        seq = load_genome(f"{wroot}/{w}/genome.fa")[chrom]
    pr = pd.concat([od, up])
    sig = Counter(); src = defaultdict(set)
    for t, g in pr.groupby("transcript_id"):
        k = (str(g.Strand.iloc[0]), tuple(sorted(zip(g.Start.astype(int), g.End.astype(int)))))
        sig[k] += 1; src[k].add(g.Source.iloc[0])
    # GMB E output genes (all transcripts) spans per strand
    genes = defaultdict(list)
    for line in open(f"{runs}/E_all/{w}/build/consensus.gff3"):
        f = line.split("\t")
        if len(f) == 9 and f[2] == "gene": genes[f[6]].append((int(f[3]) - 1, int(f[4])))
    for (strand, cds), n in sig.items():
        if n < 2: continue
        nuc = "".join(seq[s:e] for s, e in cds); nuc = reverse_complement(nuc) if strand == "-" else nuc
        if len(nuc) % 3 or len(nuc) < 99: continue
        prot = translate(nuc)
        if not (nuc.startswith("ATG") and prot.endswith("*") and "*" not in prot[:-1]): continue
        lo, hi = cds[0][0], cds[-1][1]
        if any(s < hi and e > lo for s, e in genes[strand]): continue
        tot["candidates"] += 1
        if (chrom, strand, cds) in refsig: tot["ref_exact"] += 1
        elif any(s < hi and e > lo for s, e in refspan[(chrom, strand)]): tot["ref_overlap"] += 1
        else: tot["novel"] += 1
        if n >= 5: tot["n>=5"] += 1
print(dict(tot))
