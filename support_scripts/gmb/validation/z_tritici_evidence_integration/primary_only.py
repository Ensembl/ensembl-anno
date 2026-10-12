#!/usr/bin/env python3
"""Keep only each gene's primary transcript (<gene>.t1) and its children from a GMB GFF3."""
import sys
src, dst = sys.argv[1:3]
keep_tx = set()
lines = [l for l in open(src)]
for l in lines:
    f = l.split("\t")
    if len(f) == 9 and f[2] == "mRNA":
        a = dict(kv.split("=", 1) for kv in f[8].strip().split(";") if "=" in kv)
        if a["ID"] == a["Parent"] + ".t1":
            keep_tx.add(a["ID"])
with open(dst, "w") as out:
    for l in lines:
        f = l.split("\t")
        if len(f) != 9:
            out.write(l); continue
        a = dict(kv.split("=", 1) for kv in f[8].strip().split(";") if "=" in kv)
        if f[2] == "gene" or a.get("ID") in keep_tx or a.get("Parent") in keep_tx:
            out.write(l)
