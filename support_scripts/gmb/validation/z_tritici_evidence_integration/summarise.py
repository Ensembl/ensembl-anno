#!/usr/bin/env python3
"""Aggregate annotation-qc results per arm over a window set (core / accessory pooled).

Usage: summarise.py <runs_dir> <windows.tsv> <set: dev|heldout> <arm> [<arm> ...] [--out file.tsv]
Metrics are summed over windows (R and Q are disjoint across windows).
"""
import json, os, re, sys
import pandas as pd

args = sys.argv[1:]
out = None
if "--out" in args:
    i = args.index("--out"); out = args[i + 1]; args = args[:i] + args[i + 2:]
runs, wfile, wset, *arms = args
W = pd.read_csv(wfile, sep="\t", dtype=str)
W = W[W["set"] == wset]


def n_tx(gff):
    g = t = 0
    if not os.path.exists(gff):
        return None, None
    for l in open(gff):
        f = l.split("\t")
        if len(f) == 9:
            g += f[2] == "gene"; t += f[2] == "mRNA"
    return g, t


rows = []
for arm in arms:
    for _, w in W.iterrows():
        for ev in ("consensus", "primary", "canonical"):
            d = f"{runs}/{arm}/{w['name']}"
            f = f"{d}/cmp_{ev}/comparison_summary.json"
            if not os.path.exists(f):
                continue
            s = json.load(open(f))
            sc, sp, sm, o2o = s["sensitivity_cds"], s["specificity"], s["split_merge"], s["cds_exact_one_to_one"]
            r = dict(arm=arm, window=w["name"], cls=w["class"], eval=ev,
                     R=s["total_reference_genes"], Q=s["total_consensus_genes"],
                     locus=s["sensitivity"]["locus_detected_count"],
                     cds_coord_exact=sc["cds_coordinate_exact_count"],
                     cds_struct_exact=sc["cds_exact_match_count"],
                     cds_ic=sc["cds_intron_chain_recovered"], cds_ic_eval=sc["cds_intron_chain_evaluated"],
                     partial=s["sensitivity"]["partial_match_count"],
                     missed=s["sensitivity"]["missed_count"], novel=sp["novel_consensus_count"],
                     splits=sm["gene_split_count"], merges=sm["gene_merge_count"],
                     o2o_tp=o2o["true_positives"])
            gff = {"canonical": f"{d}/fin/canonical/consensus.canonical_annotated.gff3",
                   "consensus": f"{d}/fin/consensus.gff3",
                   "primary": f"{d}/fin/primary.gff3"}[ev] if arm != "A_helixer" else f"{runs}/../win/{w['name']}/helixer.gff3"
            r["genes"], r["transcripts"] = n_tx(gff)
            log = f"{d}/build.log"
            if os.path.exists(log) and ev == "consensus":
                t = open(log).read()
                m = re.search(r"([\d.]+) real", t); r["build_s"] = float(m.group(1)) if m else None
                m = re.search(r"(\d+)\s+maximum resident", t); r["rss_mb"] = round(int(m.group(1)) / 1e6) if m else None
            rows.append(r)
df = pd.DataFrame(rows)
num = ["R", "Q", "locus", "cds_coord_exact", "cds_struct_exact", "cds_ic", "cds_ic_eval", "partial",
       "missed", "novel", "splits", "merges", "o2o_tp", "genes", "transcripts"]
agg = df.groupby(["cls", "eval", "arm"], sort=False)[num].sum()
if "build_s" in df:
    agg["build_s"] = df.groupby(["cls", "eval", "arm"], sort=False)["build_s"].sum()
agg["cds_exact_%"] = (100 * agg["cds_coord_exact"] / agg["R"]).round(1)
agg["ic_%"] = (100 * agg["cds_ic"] / agg["cds_ic_eval"]).round(1)
agg["f1"] = (2 * agg["o2o_tp"] / (agg["R"] + agg["Q"])).round(3)
pd.set_option("display.width", 250); pd.set_option("display.max_columns", 40)
print(agg.drop(columns=["cds_ic_eval"]).to_string())
if out:
    agg.reset_index().to_csv(out, sep="\t", index=False)
    df.to_csv(out.replace(".tsv", "_per_window.tsv"), sep="\t", index=False)
