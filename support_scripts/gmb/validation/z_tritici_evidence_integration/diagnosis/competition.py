#!/usr/bin/env python3
"""Pairwise competition between each backbone primary and overlapping alternative CDS structures.

Usage: competition.py <audit_dir> <window names comma-sep> <out.tsv>
"""
import sys
import pandas as pd
AD, WINS, OUT = sys.argv[1:4]
rows = []
for w in WINS.split(","):
    s = pd.read_csv(f"{AD}/{w}_structures.tsv", sep="\t")
    s = s[s.cds_bp > 0]
    s["w"] = w
    for cid, g in s.groupby("cluster"):
        prim = g[(g.role == "primary") & g.has_backbone]
        for _, p in prim.iterrows():
            comp = g[(g.key_strand == p.key_strand) & (g.cds_sig != p.cds_sig) & ~g.has_backbone
                     & (g.cds_start < p.cds_end) & (g.cds_end > p.cds_start)]
            for _, c in comp.iterrows():
                rows.append(dict(w=w, cluster=cid, p_rep=p.rep, c_rep=c.rep, c_sources=c.sources, c_role=c.role,
                    p_exact=p.ref_label == "ref_cds_exact", c_exact=c.ref_label == "ref_cds_exact",
                    p_label=p.ref_label, c_label=c.ref_label, p_score=p.score, c_score=c.score,
                    p_comp=p.prot_compatible, p_incomp=p.prot_incompatible,
                    c_comp=c.prot_compatible, c_incomp=c.prot_incompatible,
                    p_nint=p.n_introns, c_nint=c.n_introns,
                    p_sr=p.sr_introns_observed, c_sr=c.sr_introns_observed,
                    p_sr_miss=p.n_introns - p.sr_introns_observed, c_sr_miss=c.n_introns - c.sr_introns_observed,
                    p_cds=p.cds_bp, c_cds=c.cds_bp, p_ncds=p.n_cds, c_ncds=c.n_cds,
                    c_partial=bool(c.partial5) or bool(c.partial3), c_orf=c.orf, p_orf=p.orf,
                    ref_gene=p.ref_gene if p.ref_label != "no_ref" else c.ref_gene))
df = pd.DataFrame(rows); df.to_csv(OUT, sep="\t", index=False)
dec = df[df.p_exact != df.c_exact]
print("pairs", len(df), "decisive (exactly one exact)", len(dec),
      " backbone exact:", int(dec.p_exact.sum()), " alternative exact:", int(dec.c_exact.sum()))
print("decisive pairs, distinct backbone primaries:", dec.p_rep.nunique())
def test(name, pref):  # pref>0 favours competitor, <0 favours backbone, 0 tie
    good = ((pref > 0) & dec.c_exact) | ((pref < 0) & dec.p_exact)
    bad = ((pref > 0) & dec.p_exact) | ((pref < 0) & dec.c_exact)
    sw_good = ((pref > 0) & dec.c_exact).sum(); sw_bad = ((pref > 0) & dec.p_exact).sum()
    print(f"{name:42s} correct {good.sum():4d}  wrong {bad.sum():4d}  tie {(pref==0).sum():4d} | "
          f"would switch: to-exact {sw_good}, away-from-exact {sw_bad}")
netp = (dec.c_comp - dec.c_incomp) - (dec.p_comp - dec.p_incomp)
test("score (current)", dec.c_score - dec.p_score)
test("protein net compatible", netp)
fr = lambda c, i: c / (c + i).where((c + i) > 0, 1)
test("protein compatible fraction", fr(dec.c_comp, dec.c_incomp) - fr(dec.p_comp, dec.p_incomp))
test("fewer SR-unsupported introns", dec.p_sr_miss - dec.c_sr_miss)
test("longer CDS", dec.c_cds - dec.p_cds)
# conservative composite: switch only if protein net strictly better AND SR no worse AND complete ORF
comp = ((netp > 0) & (dec.c_sr_miss <= dec.p_sr_miss) & ~dec.c_partial).astype(int)
test("composite (prot>, SR>=, complete ORF)", comp)
comp2 = ((netp >= 2) & (dec.c_sr_miss <= dec.p_sr_miss) & ~dec.c_partial).astype(int)
test("composite net>=2", comp2)
