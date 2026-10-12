#!/usr/bin/env python3
"""Novel-gene stratification and isoform behaviour of canonical/all output. Usage: <runs> <wins csv> <arm> ..."""
import sys
from collections import defaultdict
import pandas as pd
runs, wins, *arms = sys.argv[1:]
pd.set_option("display.width", 250)
out = []
for arm in arms:
    rows = []
    for w in wins.split(","):
        lab = pd.read_csv(f"{runs}/{arm}/{w}/cmp_canonical/consensus_transcript_labels.tsv", sep="\t")
        if arm == "A_helixer":
            lab["evidence_sources"] = "Helixer"; lab["protein_alignments_compatible"] = None
        else:
            att = pd.read_csv(f"{runs}/{arm}/{w}/build/evidence_attribution.tsv", sep="\t")
            lab = lab.merge(att[[c for c in att.columns if c in ("transcript_id", "evidence_sources", "protein_alignments_compatible", "exon_count")]], on="transcript_id", how="left")
        lab["w"] = w; rows.append(lab)
    d = pd.concat(rows)
    d["acc"] = d.w.isin(["h16", "h19", "d14", "d18", "d21"])
    cat = d.novel_category.where(d.novel_category.notna(), d.classification)
    d["cat"] = cat.replace({"Overlaps_evaluated_reference_opposite_strand": "antisense", "Strand_Mismatch": "antisense",
                            "Novel_no_reference_overlap": "novel_intergenic", "Matched": "matched"})
    d["has_backbone"] = d.evidence_sources.fillna("").str.contains("Helixer")
    d["n_src"] = d.evidence_sources.fillna("").str.count(",") + 1
    d["has_lr"] = d.evidence_sources.fillna("").str.contains("Minimap2")
    d["prot_spliced"] = d.protein_alignments_compatible.fillna(0) > 0
    for acc, g in d.groupby("acc"):
        r = dict(arm=arm, cls="accessory" if acc else "core", genes=len(g))
        for c in ("matched", "antisense", "novel_intergenic"):
            sub = g[g.cat == c]
            r[c] = len(sub)
            if c != "matched":
                r[f"{c}_backbone"] = int(sub.has_backbone.sum())
                r[f"{c}_tx_only"] = int((~sub.has_backbone).sum())
                r[f"{c}_tx_only_multisrc_or_lr"] = int(((~sub.has_backbone) & ((sub.n_src > 1) | sub.has_lr)).sum())
                r[f"{c}_spliced_protein"] = int(sub.prot_spliced.sum())
        out.append(r)
    # isoform behaviour (all transcripts)
    if arm != "A_helixer":
        iso = []
        for w in wins.split(","):
            tx, cds = {}, defaultdict(list)
            for l in open(f"{runs}/{arm}/{w}/build/consensus.gff3"):
                f = l.rstrip("\n").split("\t")
                if len(f) != 9: continue
                a = dict(kv.split("=", 1) for kv in f[8].split(";") if "=" in kv)
                if f[2] == "mRNA": tx[a["ID"]] = a["Parent"]
                elif f[2] == "CDS": cds[a["Parent"]].append((int(f[3]), int(f[4])))
            by = defaultdict(list)
            for t, g in tx.items(): by[g].append(t)
            for g, ts in by.items():
                p = f"{g}.t1"; pc = sorted(cds.get(p, []))
                for t in ts:
                    if t == p: continue
                    c = sorted(cds.get(t, []))
                    iso.append("same_cds_utr_only" if c == pc else ("shares_coding" if any(min(e1, e2) >= max(s1, s2) for s1, e1 in pc for s2, e2 in c) else "cds_disjoint"))
            iso.append(("genes", len(by), sum(len(v) > 1 for v in by.values())))
        genes = sum(x[1] for x in iso if isinstance(x, tuple)); multi = sum(x[2] for x in iso if isinstance(x, tuple))
        alts = [x for x in iso if not isinstance(x, tuple)]
        print(f"{arm:10s} isoforms: genes {genes} multi-isoform {multi} ({100*multi/genes:.1f}%) alternates {len(alts)} "
              + " ".join(f"{k}={v}" for k, v in pd.Series(alts).value_counts().items()))
print(pd.DataFrame(out).to_string(index=False))
