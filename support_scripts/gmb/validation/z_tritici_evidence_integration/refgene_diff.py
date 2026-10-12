#!/usr/bin/env python3
"""Per-reference-gene outcome changes between two arms (same eval).

Usage: refgene_diff.py <runs> <armA> <armB> <eval> <windows csv> [out.tsv]
"""
import sys, pandas as pd
runs, A, B, ev, wins = sys.argv[1:6]
out = sys.argv[6] if len(sys.argv) > 6 else None
def load(arm, w):
    d = pd.read_csv(f"{runs}/{arm}/{w}/cmp_{ev}/comparison_details.tsv", sep="\t", low_memory=False)
    d = d[d.source == "reference"].set_index("gene_id")
    return d[["classification", "cds_coordinate_exact", "cds_intron_chain_match", "best_cds_match_transcript_id", "start", "end", "strand"]]
rows = []
for w in wins.split(","):
    a, b = load(A, w), load(B, w)
    j = a.join(b, lsuffix="_A", rsuffix="_B")
    j["w"] = w
    rows.append(j)
j = pd.concat(rows)
def state(r, s):
    if r[f"classification_{s}"] == "Missed": return "missed"
    if r[f"cds_coordinate_exact_{s}"] in (True, "True"): return "exact"
    if r[f"cds_intron_chain_match_{s}"] in (True, "True"): return "ic"
    return "partial"
j["sA"] = j.apply(lambda r: state(r, "A"), axis=1); j["sB"] = j.apply(lambda r: state(r, "B"), axis=1)
ch = j[j.sA != j.sB]
print(pd.crosstab(ch.sA, ch.sB))
if out:
    ch.to_csv(out, sep="\t")
print(ch[ch.sB == "missed"][["w", "start_A", "end_A", "strand_A", "best_cds_match_transcript_id_A", "sA"]])
