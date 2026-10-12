#!/usr/bin/env python3
"""Trace reference genes that Helixer misses through GMB's candidate stages (dev windows).

Uses annotation-qc's Missed classification for raw Helixer, then for each such gene looks at
(1) raw short-read transcripts overlapping its CDS on the same strand, incl. those filtered,
(2) merged candidate structures from locus_audit (ORF, protein support, score, role),
(3) whether the all-evidence GMB build recovers it.
"""
import sys, pandas as pd
runs, audit, wins, out = sys.argv[1:5]
rows = []
for w in wins.split(","):
    hx = pd.read_csv(f"{runs}/A_helixer/{w}/cmp_consensus/comparison_details.tsv", sep="\t", low_memory=False)
    hx = hx[(hx.source == "reference") & (hx.classification == "Missed")]
    e = pd.read_csv(f"{runs}/E_all/{w}/cmp_primary/comparison_details.tsv", sep="\t", low_memory=False)
    e = e[e.source == "reference"].set_index("gene_id")
    st = pd.read_csv(f"{audit}/{w}_structures.tsv", sep="\t")
    fi = pd.read_csv(f"{audit}/{w}_filtered_inputs.tsv", sep="\t") if True else None
    for _, g in hx.iterrows():
        gid = g.gene_id
        sub = st[st.ref_gene.isin([gid, f"gene:{gid}"])]
        lost = fi[fi.ref_gene.isin([gid, f"gene:{gid}"])] if len(fi) else fi
        best = None
        if len(sub):
            order = {"ref_cds_exact": 0, "ref_cds_ic": 1, "ref_overlap": 2}
            sub = sub.assign(o=sub.ref_label.map(order).fillna(9)).sort_values(["o", "score"], ascending=[True, False])
            best = sub.iloc[0]
        ecls = e.loc[gid, "classification"] if gid in e.index else "?"
        e_exact = e.loc[gid, "cds_coordinate_exact"] if gid in e.index else "?"
        rows.append(dict(window=w, ref_gene=gid, span=g.end - g.start + 1,
            n_struct=len(sub), filtered_inputs=len(lost),
            filtered_exact=int((lost.ref_label == "ref_cds_exact").sum()) if len(lost) else 0,
            best_label=None if best is None else best.ref_label,
            best_sources=None if best is None else best.sources,
            best_exons=None if best is None else best.exon_count,
            best_orf=None if best is None else best.orf,
            best_prot=None if best is None else best.protein_positional,
            best_score=None if best is None else best.score,
            best_role=None if best is None else best.role,
            gmb_classification=ecls, gmb_cds_exact=e_exact))
df = pd.DataFrame(rows); df.to_csv(out, sep="\t", index=False)
pd.set_option("display.width", 250); pd.set_option("display.max_rows", 200)
print(df.drop(columns=["window"]).to_string(index=False))
