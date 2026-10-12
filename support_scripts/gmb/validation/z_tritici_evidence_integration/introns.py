#!/usr/bin/env python3
"""Pooled intron-level precision/recall from annotation-qc. Usage: introns.py <runs> <wins csv> <arm:eval> ..."""
import json, sys
runs, wins, *specs = sys.argv[1:]
for spec in specs:
    arm, ev = spec.split(":")
    t = dict(cr=0, cq=0, cs=0, er=0, eq=0, es=0)
    for w in wins.split(","):
        s = json.load(open(f"{runs}/{arm}/{w}/cmp_{ev}/comparison_summary.json"))["intron_support"]
        t["cr"] += s["cds"]["reference_introns"]; t["cq"] += s["cds"]["query_introns"]; t["cs"] += s["cds"]["shared_introns"]
        t["er"] += s["exon"]["reference_introns"]; t["eq"] += s["exon"]["query_introns"]; t["es"] += s["exon"]["shared_introns"]
    print(f"{spec:22s} CDS introns query {t['cq']:5d} shared {t['cs']:5d} P {t['cs']/t['cq']:.3f} R {t['cs']/t['cr']:.3f}"
          f" | exon introns query {t['eq']:5d} P {t['es']/t['eq']:.3f} R {t['es']/t['er']:.3f}")
