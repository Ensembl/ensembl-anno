#!/usr/bin/env python3
"""Independent structural check of a GMB handover (finalise/ + build/).

Checks, without using GMB code:
  * GFF3: unique IDs; every mRNA has an existing gene parent on the same sequence/strand;
    every exon/CDS/UTR has an existing mRNA parent; CDS and UTR segments lie inside exons;
    gene span == union of its mRNAs; mRNA span == union of its exons
  * exactly one Ensembl_canonical transcript per gene in the canonical GFF3, and the canonical
    GFF3 holds the same transcripts as finalise/consensus.gff3
  * prot.fa / cds.fa: one record per coding mRNA, CDS length multiple of 3, no internal stop
  * evidence_attribution.tsv: one row per mRNA, same IDs
  * QC reports pass
Usage: check_handover.py <finalise_dir> <build_dir>     (exit 1 on any failure)
"""
import json
import sys
from collections import defaultdict
from pathlib import Path

fin, build = Path(sys.argv[1]), Path(sys.argv[2])
errors = []


def err(msg):
    errors.append(msg)


def parse(path):
    feats = []
    for line in open(path):
        if line.startswith("#") or not line.strip():
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) != 9:
            err(f"{path.name}: line with {len(f)} columns")
            continue
        attrs = dict(kv.split("=", 1) for kv in f[8].split(";") if "=" in kv)
        feats.append((f[0], f[2], int(f[3]), int(f[4]), f[6], attrs))
    return feats


def check_gff(path):
    feats = parse(path)
    ids = [a["ID"] for *_x, a in feats if "ID" in a]
    if len(ids) != len(set(ids)):
        err(f"{path.name}: duplicate IDs")
    genes = {a["ID"]: (c, s, e, st) for c, t, s, e, st, a in feats if t == "gene"}
    mrnas = {a["ID"]: (c, s, e, st, a.get("Parent")) for c, t, s, e, st, a in feats if t == "mRNA"}
    kids = defaultdict(lambda: defaultdict(list))
    for c, t, s, e, st, a in feats:
        if t in ("exon", "CDS", "five_prime_UTR", "three_prime_UTR"):
            p = a.get("Parent")
            if p not in mrnas:
                err(f"{path.name}: {t} with missing parent {p}")
                continue
            kids[p][t].append((s, e))
    by_gene = defaultdict(list)
    for m, (c, s, e, st, g) in mrnas.items():
        if g not in genes:
            err(f"{path.name}: mRNA {m} parent {g} missing")
            continue
        gc, gs, ge, gst = genes[g]
        if (c, st) != (gc, gst):
            err(f"{path.name}: mRNA {m} on other sequence/strand than gene {g}")
        by_gene[g].append((s, e))
        ex = sorted(kids[m]["exon"])
        if not ex:
            err(f"{path.name}: mRNA {m} without exons")
            continue
        if (ex[0][0], max(x[1] for x in ex)) != (s, e):
            err(f"{path.name}: mRNA {m} span != exon union")
        for t in ("CDS", "five_prime_UTR", "three_prime_UTR"):
            for a, b in kids[m][t]:
                if not any(xs <= a and b <= xe for xs, xe in ex):
                    err(f"{path.name}: {t} {a}-{b} of {m} outside exons")
    for g, spans in by_gene.items():
        if (min(s for s, _ in spans), max(e for _, e in spans)) != genes[g][1:3]:
            err(f"{path.name}: gene {g} span != union of its mRNAs")
    for g in genes:
        if g not in by_gene:
            err(f"{path.name}: gene {g} without mRNA")
    return feats, mrnas, kids


feats, mrnas, kids = check_gff(fin / "consensus.gff3")
cfeats, cmrnas, _ = check_gff(fin / "canonical" / "consensus.canonical_annotated.gff3")
if set(cmrnas) != set(mrnas):
    err("canonical GFF3 and consensus.gff3 hold different transcripts")
canon = defaultdict(int)
for c, t, s, e, st, a in cfeats:
    if t == "mRNA" and "Ensembl_canonical" in a.get("tag", "").split(","):
        canon[a["Parent"]] += 1
genes = {a["ID"] for c, t, s, e, st, a in cfeats if t == "gene"}
bad = [g for g in genes if canon.get(g, 0) != 1]
if bad:
    err(f"{len(bad)} gene(s) without exactly one canonical transcript (e.g. {bad[:3]})")


def fasta(path):
    seqs, cur = {}, None
    for line in open(path):
        if line.startswith(">"):
            cur = line[1:].split()[0]
            seqs[cur] = []
        else:
            seqs[cur].append(line.strip())
    return {k: "".join(v) for k, v in seqs.items()}


coding = {m for m in mrnas if kids[m]["CDS"]}
cds, prot = fasta(fin / "cds.fa"), fasta(fin / "prot.fa")
if set(cds) != coding or set(prot) != coding:
    err("cds.fa / prot.fa records differ from coding mRNAs in consensus.gff3")
for m, seq in cds.items():
    if len(seq) % 3:
        err(f"cds.fa {m}: length not a multiple of 3")
    if sum(e - s + 1 for s, e in kids[m]["CDS"]) != len(seq):
        err(f"cds.fa {m}: length differs from GFF3 CDS")
for m, seq in prot.items():
    if "*" in seq:
        err(f"prot.fa {m}: internal stop")

import csv
with open(build / "evidence_attribution.tsv") as fh:
    att = [r["transcript_id"] for r in csv.DictReader(fh, delimiter="\t")]
if sorted(att) != sorted(mrnas):
    err("evidence_attribution.tsv rows differ from mRNAs")
fq = json.load(open(fin / "fasta_qc_report.json"))
if fq.get("pass") is not True:
    err(f"fasta_qc_report.json: pass={fq.get('pass')} {fq.get('failed_checks')}")
uq = json.load(open(fin / "utr_qc_report.json"))
if uq.get("invariant_ok") is not True or uq.get("violations"):
    err(f"utr_qc_report.json: invariant_ok={uq.get('invariant_ok')} violations={uq.get('violations')}")

n_genes = len({a["ID"] for c, t, s, e, st, a in feats if t == "gene"})
print(f"{fin}: {n_genes} genes, {len(mrnas)} transcripts, {len(coding)} coding -- "
      f"{'OK' if not errors else f'{len(errors)} problem(s)'}")
for e in errors[:20]:
    print("  FAIL", e)
sys.exit(1 if errors else 0)
