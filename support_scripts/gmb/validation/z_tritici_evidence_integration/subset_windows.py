#!/usr/bin/env python3
"""Write per-window subsets of the Z. tritici GMB inputs in one pass per file.

A model is kept only when its whole transcript (GTF) / gene (GFF3) lies inside the
window; all its rows are then written. The genome subset is the whole chromosome so
coordinates are unchanged. Input files are only read.

Usage: subset_windows.py <input_data_dir> <helixer_remapped.gff3> <genome.fa> <windows.tsv> <outroot>
"""
import re
import sys
from collections import defaultdict
from pathlib import Path

IN, HELIXER, GENOME, WIN, OUT = sys.argv[1:6]
wins = []
for i, line in enumerate(open(WIN)):
    if i == 0:
        continue
    name, _set, chrom, s, e, _cls = line.rstrip("\n").split("\t")
    wins.append((name, chrom, int(s), int(e)))
by_chrom = defaultdict(list)
for w in wins:
    by_chrom[w[1]].append(w)
for name, *_ in wins:
    (Path(OUT) / name).mkdir(parents=True, exist_ok=True)


def window_of(chrom, s, e):
    for name, _c, ws, we in by_chrom.get(chrom, ()):
        if s >= ws and e <= we:
            return name
    return None


TID = re.compile(r'transcript_id "([^"]+)"')
tracks = {"scallop": "scallop", "stringtie": "stringtie", "orthodb": "orthodb",
          "genblast": "genblast", "minimap2": "minimap2"}
for track in tracks:
    path = f"{IN}/GCA_000219625.1_{track}_annotation.gtf"
    keep = {}           # tid -> window (transcript row inside a window)
    rows = defaultdict(list)
    for line in open(path):
        if line.startswith("#"):
            continue
        f = line.split("\t", 8)
        if len(f) < 9 or f[0] not in by_chrom:
            continue
        m = TID.search(f[8])
        if not m:
            continue
        # genBlastG reuses alignment IDs on different sequences: key by (seq, ID)
        # so another sequence's rows are never copied into this window.
        tid = (f[0], m.group(1))
        w = window_of(f[0], int(f[3]), int(f[4]))
        if w is None:
            continue
        if f[2] == "transcript":
            keep[tid] = w
        rows[tid].append(line)
    handles = {name: open(Path(OUT) / name / f"{track}.gtf", "w") for name, *_ in wins}
    n = defaultdict(int)
    for tid, w in keep.items():
        for line in rows[tid]:
            handles[w].write(line)
        n[w] += 1
    for h in handles.values():
        h.close()
    print(track, dict(n), flush=True)

gene_of = re.compile(r"^(_[A-Za-z0-9.]+?_\d{6})")
genes, rows = {}, defaultdict(list)
for line in open(HELIXER):
    if line.startswith("#"):
        continue
    f = line.split("\t")
    if len(f) < 9 or f[0] not in by_chrom:
        continue
    gid = gene_of.match(f[8].split(";")[0].split("=", 1)[1]).group(1)
    if f[2] == "gene":
        w = window_of(f[0], int(f[3]), int(f[4]))
        if w:
            genes[gid] = w
    rows[gid].append(line)
handles = {name: open(Path(OUT) / name / "helixer.gff3", "w") for name, *_ in wins}
for h in handles.values():
    h.write("##gff-version 3\n")
n = defaultdict(int)
for gid, w in genes.items():
    for line in rows[gid]:
        handles[w].write(line)
    n[w] += 1
for h in handles.values():
    h.close()
print("helixer", dict(n))

seqs, cur = {}, None
for line in open(GENOME):
    if line.startswith(">"):
        cur = line[1:].split()[0]
        seqs[cur] = [line]
    else:
        seqs[cur].append(line)
for name, chrom, *_ in wins:
    with open(Path(OUT) / name / "genome.fa", "w") as fh:
        fh.writelines(seqs[chrom])
