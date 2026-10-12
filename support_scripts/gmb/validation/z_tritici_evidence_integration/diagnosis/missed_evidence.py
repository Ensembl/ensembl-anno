#!/usr/bin/env python3
"""For reference genes Helixer misses: what raw evidence overlaps the reference CDS?

Per gene: same-strand raw short-read transcripts whose exons overlap the reference CDS, how
many contain the whole reference CDS within their exons, whether GMB's inferred ORF on that
transcript overlaps the reference CDS, whether a transcript ORF reproduces the reference CDS
exactly, and same-strand protein alignments overlapping the CDS (and any exact CDS match).
"""
import contextlib, io, os, sys
from collections import defaultdict
import pandas as pd
from gmb.pipeline.annotate_cds_utrs import annotate_all_transcripts, load_genome, build_spliced_seq
from gmb.pipeline.builder import load_evidence

runs, wroot, ref, wins, out = sys.argv[1:6]
W = {l.split("\t")[0]: l.split("\t") for l in open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "windows.tsv")) if not l.startswith("name")}
rcds = defaultdict(list); rstrand = {}; tx_gene = {}
for line in open(ref):
    f = line.rstrip("\n").split("\t")
    if len(f) != 9: continue
    a = dict(kv.split("=", 1) for kv in f[8].split(";") if "=" in kv)
    if f[2] == "mRNA": tx_gene[a["ID"]] = a["Parent"].replace("gene:", "")
    elif f[2] == "CDS": rcds[a["Parent"]].append((int(f[3]) - 1, int(f[4]))); rstrand[a["Parent"]] = f[6]
gene_cds = {tx_gene[t]: (sorted(v), rstrand[t]) for t, v in rcds.items() if t in tx_gene}

def ov(a, b): return sum(max(0, min(e1, e2) - max(s1, s2)) for s1, e1 in a for s2, e2 in b)
rows = []
for w in wins.split(","):
    hx = pd.read_csv(f"{runs}/A_helixer/{w}/cmp_consensus/comparison_details.tsv", sep="\t", low_memory=False)
    missed = hx[(hx.source == "reference") & (hx.classification == "Missed")].gene_id.tolist()
    with contextlib.redirect_stdout(io.StringIO()):
        sc, _ = load_evidence(f"{wroot}/{w}/scallop.gtf", "Scallop")
        st, _ = load_evidence(f"{wroot}/{w}/stringtie.gtf", "StringTie")
        od, _ = load_evidence(f"{wroot}/{w}/orthodb.gtf", "OrthoDB")
        up, _ = load_evidence(f"{wroot}/{w}/genblast.gtf", "UniProt")
        genome = load_genome(f"{wroot}/{w}/genome.fa")
    tx = pd.concat([sc, st]); pr_ = pd.concat([od, up])
    txg = {t: (str(g.Strand.iloc[0]), sorted(zip(g.Start, g.End))) for t, g in tx.groupby("transcript_id")}
    prg = {t: (str(g.Strand.iloc[0]), sorted(zip(g.Start, g.End))) for t, g in pr_.groupby("transcript_id")}
    chrom = W[w][2]; seq = genome[chrom]
    for gid in missed:
        cds, strand = gene_cds[gid]
        lo, hi = cds[0][0], cds[-1][1]; cbp = sum(e - s for s, e in cds)
        hits = [(t, ex) for t, (s, ex) in txg.items() if s == strand and ex[0][0] < hi and ex[-1][1] > lo and ov(ex, cds) > 0]
        contain = orf_on = orf_exact = alt_orf_exact = 0
        spans = []
        for t, ex in hits:
            if ov(ex, cds) == cbp: contain += 1
            a = annotate_all_transcripts(pd.DataFrame({"Chromosome": chrom, "Start": [s for s, _ in ex], "End": [e for _, e in ex], "Strand": strand, "transcript_id": t}), genome, None, min_codons=33)[t]
            c = sorted(a.get("cds") or [])
            if c and ov(c, cds) > 0: orf_on += 1
            if c == cds: orf_exact += 1
            spans.append(ex[-1][1] - ex[0][0])
            # would ANY ATG..stop ORF of >= 33 codons on this transcript reproduce the reference CDS?
            cdna = build_spliced_seq(ex, strand, seq)
            from gmb.pipeline.annotate_cds_utrs import map_cds_to_genomic
            for fr in range(3):
                i = fr
                while i + 3 <= len(cdna):
                    if cdna[i:i + 3] == "ATG":
                        j = i
                        while j + 3 <= len(cdna) and cdna[j:j + 3] not in ("TAA", "TAG", "TGA"): j += 3
                        if j + 3 <= len(cdna) and sorted(map_cds_to_genomic(i, j + 3, ex, strand)) == cds:
                            alt_orf_exact += 1; break
                    i += 3
        ph = [t for t, (s, ex) in prg.items() if s == strand and ex[0][0] < hi and ex[-1][1] > lo and ov(ex, cds) > 0]
        p_exact = sum(1 for t in ph if prg[t][1] == cds)
        p_chain = sum(1 for t in ph if len(cds) > 1 and [(prg[t][1][i][1], prg[t][1][i + 1][0]) for i in range(len(prg[t][1]) - 1)] == [(cds[i][1], cds[i + 1][0]) for i in range(len(cds) - 1)])
        rows.append(dict(window=w, ref_gene=gid, cds_bp=cbp, n_cds=len(cds), sr_overlapping=len(hits), sr_containing_cds=contain,
                         sr_orf_overlaps_ref=orf_on, sr_orf_exact=orf_exact, sr_alt_atg_orf_exact=min(alt_orf_exact, len(hits)),
                         median_sr_span=int(pd.Series(spans).median()) if spans else None,
                         prot_overlapping=len(ph), prot_cds_exact=p_exact, prot_chain_exact=p_chain))
df = pd.DataFrame(rows); df.to_csv(out, sep="\t", index=False)
df["acc"] = df.window.isin(["d14", "d18", "d21"])
for acc, d in df.groupby("acc"):
    print("accessory" if acc else "core", len(d), "missed by Helixer")
    print("  any same-strand short-read overlap:", (d.sr_overlapping > 0).sum(),
          "| transcript contains whole ref CDS:", (d.sr_containing_cds > 0).sum(),
          "| GMB ORF overlaps ref CDS:", (d.sr_orf_overlaps_ref > 0).sum(),
          "| GMB ORF exact:", (d.sr_orf_exact > 0).sum(),
          "| some ATG ORF on a transcript is exact:", (d.sr_alt_atg_orf_exact > 0).sum())
    print("  protein overlap:", (d.prot_overlapping > 0).sum(), "| protein CDS exact:", (d.prot_cds_exact > 0).sum(),
          "| protein intron chain exact (multi-CDS):", (d.prot_chain_exact > 0).sum(), "of", (d.n_cds > 1).sum())
    print("  no SR and no protein:", ((d.sr_overlapping == 0) & (d.prot_overlapping == 0)).sum(), "| single-CDS-segment:", (d.n_cds == 1).sum(), "| median CDS bp:", int(d.cds_bp.median()))
