#!/usr/bin/env python3
"""Per-structure audit of GMB candidate generation, scoring and selection in one window.

Runs GMB's own loaders, filters, ORF annotation, locus clustering and select_isoforms()
(the build's selection stage, without the GFF3 post-processing) with the fungi preset
plus optional --config overlays, and writes one row per *merged candidate structure*
(the unit select_isoforms scores) with:

  sources, n members, exon/CDS counts, ORF integrity (ATG/stop/partial/internal stops),
  numeric score and its components, protein support (positional / CDS-span / CDS-intron
  agreement from the spliced genBlastG alignments), short-read junction support,
  retention-gate outcome, ranking position, final role (primary / alternate / dropped),
  and the reference label of its CDS (exact / intron-chain / overlap / none) plus the
  reference gene it overlaps.

It also writes one row per reference gene summarising what candidates were available
and what was selected. The reference is used ONLY for labelling after selection.

Usage: locus_audit.py <window_dir> <ref.gff3> <chrom> <start> <end> <outprefix> [config.yaml ...]
"""
import contextlib
import io
import sys
from collections import defaultdict

import pandas as pd

from gmb.pipeline import scoring
from gmb.pipeline.annotate_cds_utrs import annotate_all_transcripts, load_genome
from gmb.pipeline.builder import cluster_candidate_loci, load_evidence
from gmb.pipeline.canonical_evidence import EvidenceRoles
from gmb.pipeline.config import load_config
from gmb.pipeline.evidence_filter import (filter_backbone_models, filter_chimeras,
                                          filter_protein_evidence, split_mega_transcripts)
from gmb.utils.intervals import cds_span_compatible_ids, same_strand_overlap_ids

WD, REF, CHROM, WS, WE, OUT, *CONFIGS = sys.argv[1:]
WS, WE = int(WS), int(WE)
cfg = load_config(CONFIGS or None, "fungi")
cfg.scoring.backbone_label = "Helixer"
roles = EvidenceRoles.from_config(cfg.scoring)
st = {}
with contextlib.redirect_stdout(io.StringIO()):
    sc, _ = load_evidence(f"{WD}/scallop.gtf", "Scallop")
    sr, _ = load_evidence(f"{WD}/stringtie.gtf", "StringTie")
    hx, hx_cds = load_evidence(f"{WD}/helixer.gff3", "Helixer")
    od, _ = load_evidence(f"{WD}/orthodb.gtf", "OrthoDB")
    up, _ = load_evidence(f"{WD}/genblast.gtf", "UniProt")
    tx_raw = pd.concat([sc, sr], ignore_index=True)
    tx_f1 = filter_chimeras(tx_raw, cfg, st)
    tx = split_mega_transcripts(tx_f1, cfg, st) if cfg.transcript_splitting.split_enabled else tx_f1
    he, hc = filter_backbone_models(hx, hx_cds, cfg, st)
    prot_raw = pd.concat([od, up], ignore_index=True)
    prot = filter_protein_evidence(prot_raw, cfg, st, tx)
    genome = load_genome(f"{WD}/genome.fa")
    cand = pd.concat([tx, he], ignore_index=True)
    ann = annotate_all_transcripts(cand, genome, hc, min_codons=cfg.orf.min_codons)

# ---------------------------------------------------------------- input fates
fate_in = {}
for tid in tx_raw["transcript_id"].unique():
    fate_in[tid] = "kept"
removed_f1 = set(tx_raw["transcript_id"]) - set(tx_f1["transcript_id"])
for tid in removed_f1:
    g = tx_raw[tx_raw.transcript_id == tid]
    span = g.End.max() - g.Start.min()
    fate_in[tid] = "filtered_long_span" if (cfg.transcriptomic_filter.max_transcript_length
                                            and span > cfg.transcriptomic_filter.max_transcript_length) \
        else "filtered_long_intron"

# ---------------------------------------------------------------- evidence indexes
def introns_of(ivs):
    ivs = sorted(ivs)
    return [(ivs[i][1], ivs[i + 1][0]) for i in range(len(ivs) - 1)]


def build_intron_index(df):
    """(strand, start, end) -> set of sources, from exon rows grouped by transcript."""
    idx = defaultdict(set)
    for (src, tid), g in df.groupby(["Source", "transcript_id"], observed=True):
        strand = str(g["Strand"].iloc[0])
        for i in introns_of(zip(g.Start.astype(int), g.End.astype(int))):
            idx[(strand, *i)].add(src)
    return idx


sr_introns_raw = build_intron_index(tx_raw)     # every observed short-read junction
prot_introns = build_intron_index(prot)          # spliced protein alignments (post-filter)
hx_introns = build_intron_index(he)

prot_spans = [(str(g.Strand.iloc[0]), int(g.Start.min()), int(g.End.max()), tid,
               introns_of(zip(g.Start.astype(int), g.End.astype(int))))
              for tid, g in prot.groupby("transcript_id")]

psup = same_strand_overlap_ids(cand, prot)
cand_spans = {t: (str(g.Chromosome.iloc[0]), str(g.Strand.iloc[0]), int(g.Start.min()), int(g.End.max()))
              for t, g in cand.groupby("transcript_id")}
pspans = {t: (str(g.Chromosome.iloc[0]), str(g.Strand.iloc[0]), int(g.Start.min()), int(g.End.max()))
          for t, g in prot.groupby("transcript_id")}
with_cds = {t for t, a in ann.items() if a and a.get("cds")}
pspan_ok = cds_span_compatible_ids(cand_spans, pspans, with_cds)
cand_cds = {t: (a.get("cds") or []) for t, a in ann.items()}


def protein_cds_agreement(strand, cds):
    """CDS-intron agreement with spliced protein alignments.

    confirmed  : candidate CDS introns observed exactly in >= 1 protein alignment
    conflicting: protein-alignment introns lying inside the candidate CDS span that the
                 candidate does not have (the protein says 'splice here', the model does not)
    best_chain : True if some protein alignment's whole intron chain is a contiguous
                 sub-chain of the candidate's CDS intron chain (and >= 1 intron)
    """
    if not cds:
        return 0, 0, False, 0
    ci = introns_of(cds)
    cset = set(ci)
    lo, hi = cds[0][0], cds[-1][1]
    confirmed = sum(1 for i in ci if (strand, *i) in prot_introns)
    conflicting = set()
    best_chain = False
    n_contained = 0
    for s, ps, pe, _tid, pin in prot_spans:
        if s != strand or pe <= lo or ps >= hi:
            continue
        for i in pin:
            if lo < i[0] and i[1] < hi and i not in cset:
                conflicting.add(i)
        if ps >= lo and pe <= hi:
            n_contained += 1
            if pin and all(i in cset for i in pin):
                best_chain = True
    return confirmed, len(conflicting), best_chain, n_contained


def protein_compatibility(strand, cds):
    """Count overlapping spliced protein alignments compatible / incompatible with the CDS.

    Over the overlap of an alignment with the candidate CDS span, the alignment is
    compatible when its introns there equal the candidate's CDS introns there and it has
    >= 1 intron in the overlap (spliced, informative); incompatible when they differ.
    Unspliced overlaps (no intron on either side) count as neither.
    """
    if not cds:
        return 0, 0
    ci = introns_of(cds)
    lo, hi = cds[0][0], cds[-1][1]
    comp = incomp = 0
    for s, ps, pe, _tid, pin in prot_spans:
        if s != strand or pe <= lo or ps >= hi:
            continue
        a, b = max(lo, ps), min(hi, pe)
        if b - a < 30:
            continue
        p_in = {i for i in pin if i[0] >= a and i[1] <= b}
        c_in = {i for i in ci if i[0] >= a and i[1] <= b}
        if not p_in and not c_in:
            continue
        if p_in == c_in:
            comp += 1
        else:
            incomp += 1
    return comp, incomp


# ---------------------------------------------------------------- reference
ref_cds = defaultdict(list); ref_strand = {}; ref_gene = {}; tx_gene = {}; pc_genes = set()
for line in open(REF):
    if line.startswith("#"):
        continue
    f = line.rstrip("\n").split("\t")
    if len(f) < 9 or f[0] != CHROM:
        continue
    a = dict(kv.split("=", 1) for kv in f[8].split(";") if "=" in kv)
    if f[2] == "gene" and a.get("biotype") == "protein_coding" and int(f[3]) >= WS and int(f[4]) <= WE:
        pc_genes.add(a["ID"])
    elif f[2] == "mRNA":
        tx_gene[a["ID"]] = a["Parent"]
    elif f[2] == "CDS":
        ref_cds[a["Parent"]].append((int(f[3]) - 1, int(f[4])))
        ref_strand[a["Parent"]] = f[6]
ref_cds = {t: sorted(v) for t, v in ref_cds.items() if tx_gene.get(t) in pc_genes}
ref_by_sig = {(ref_strand[t], tuple(v)): tx_gene[t] for t, v in ref_cds.items()}
ref_by_chain = {(ref_strand[t], tuple(introns_of(v))): tx_gene[t] for t, v in ref_cds.items() if len(v) > 1}
ref_spans = [(ref_strand[t], v[0][0], v[-1][1], tx_gene[t]) for t, v in ref_cds.items()]


def ref_label(strand, cds):
    if not cds:
        return "no_cds", ""
    cds = tuple(sorted(cds))
    if (strand, cds) in ref_by_sig:
        return "ref_cds_exact", ref_by_sig[(strand, cds)]
    if len(cds) > 1 and (strand, tuple(introns_of(cds))) in ref_by_chain:
        return "ref_cds_ic", ref_by_chain[(strand, tuple(introns_of(cds)))]
    lo, hi = cds[0][0], cds[-1][1]
    best, bestg = 0, ""
    for s, b, e, g in ref_spans:
        if s == strand:
            ov = min(hi, e) - max(lo, b)
            if ov > best:
                best, bestg = ov, g
    if best > 0:
        return "ref_overlap", bestg
    if any(min(hi, e) > max(lo, b) for s, b, e, g in ref_spans):
        return "ref_antisense", ""
    return "no_ref", ""


# ---------------------------------------------------------------- run selection, instrumented
scored = {}
orig_score = scoring.score_model


def spy(model, *a, **k):
    s = orig_score(model, *a, **k)
    scored[model["id"]] = s
    return s


scoring.score_model = spy
clusters = cluster_candidate_loci(cand, cfg.scoring.locus_clustering)
role = {}
gene_of = {}
gidx = 0
for cid, locus in clusters.groupby("Cluster"):
    genes = scoring.select_isoforms(locus, cfg, psup, genome,
                                    protein_cds_span_tids=pspan_ok, candidate_cds=cand_cds,
                                    canonical_intron_tids=None)
    for gene in genes:
        gidx += 1
        for i, m in enumerate(gene):
            role[m["id"]] = "primary" if i == 0 else "alternate"
            gene_of[m["id"]] = gidx
scoring.score_model = orig_score

# ---------------------------------------------------------------- per-structure table
rows = []
for cid, locus in clusters.groupby("Cluster"):
    merged = {}
    for (src, tid), g in locus.groupby(["Source", "transcript_id"], observed=True):
        ex = sorted(zip(g.Start.astype(int), g.End.astype(int)))
        strand = str(g.Strand.iloc[0])
        chain = introns_of(ex)
        key = (strand, tuple(chain)) if chain else (strand, ex[0][0], ex[-1][1])
        merged.setdefault(key, []).append((src, tid, ex, strand))
    for key, members in merged.items():
        srcs = sorted({m[0] for m in members})
        ids = [m[1] for m in members]
        rep = next((i for i in ids if i in scored), ids[0])
        r_role = next((role[i] for i in ids if i in role), "dropped")
        a = ann.get(rep) or {}
        strand = members[0][3]
        cds = sorted(a.get("cds") or [])
        prot_seq = a.get("protein") or ""
        lab, rg = ref_label(strand, cds)
        conf, confl, pchain, ncont = protein_cds_agreement(strand, cds)
        pcomp, pincomp = protein_compatibility(strand, cds)
        ex = members[0][2]
        ints = introns_of(ex)
        cints = introns_of(cds)
        sr_ok = sum(1 for i in ints if (strand, *i) in sr_introns_raw)
        sr2 = sum(1 for i in ints if len(sr_introns_raw.get((strand, *i), ())) >= 2)
        hx_ok = sum(1 for i in ints if (strand, *i) in hx_introns)
        rows.append(dict(
            cluster=cid, key_strand=strand, start=ex[0][0], end=ex[-1][1], sources=",".join(srcs),
            n_members=len(members), rep=rep, members=";".join(ids),
            has_backbone="Helixer" in srcs, exon_count=len(ex), n_introns=len(ints),
            cds_bp=sum(e - s for s, e in cds), n_cds=len(cds), n_cds_introns=len(cints),
            orf=a.get("orf_label", ""), partial5=a.get("is_partial_5"), partial3=a.get("is_partial_3"),
            internal_stops=prot_seq.count("*"),
            score=scored.get(rep), protein_positional=any(i in psup for i in ids),
            protein_cds_span=any(i in pspan_ok for i in ids),
            prot_cds_introns_confirmed=conf, prot_conflicting_introns=confl,
            prot_chain_contained=pchain, prot_alignments_inside=ncont,
            prot_compatible=pcomp, prot_incompatible=pincomp,
            cds_start=cds[0][0] if cds else None, cds_end=cds[-1][1] if cds else None,
            cds_sig=";".join(f"{a}-{b}" for a, b in cds),
            sr_introns_observed=sr_ok, sr_introns_both_tracks=sr2, hx_introns_shared=hx_ok,
            role=r_role, gene=gene_of.get(next((i for i in ids if i in gene_of), None)),
            ref_label=lab, ref_gene=rg,
            input_fate=";".join(sorted({fate_in.get(i.split("__seg")[0], "backbone") for i in ids})),
        ))
df = pd.DataFrame(rows)
df.to_csv(f"{OUT}_structures.tsv", sep="\t", index=False)

# ---------------------------------------------------------------- dropped-before-candidacy
# transcripts removed by the chimera filter that would have had a reference-exact CDS
lost = []
if removed_f1:
    ann_lost = annotate_all_transcripts(tx_raw[tx_raw.transcript_id.isin(removed_f1)], genome, None,
                                        min_codons=cfg.orf.min_codons)
    for tid, a in ann_lost.items():
        g = tx_raw[tx_raw.transcript_id == tid]
        lab, rg = ref_label(str(g.Strand.iloc[0]), sorted(a.get("cds") or []))
        lost.append(dict(tid=tid, fate=fate_in[tid], span=int(g.End.max() - g.Start.min()),
                         ref_label=lab, ref_gene=rg))
pd.DataFrame(lost).to_csv(f"{OUT}_filtered_inputs.tsv", sep="\t", index=False)

# ---------------------------------------------------------------- per-reference-gene table
order = {"ref_cds_exact": 0, "ref_cds_ic": 1, "ref_overlap": 2}
gr = []
for g in sorted(pc_genes):
    sub = df[df.ref_gene == g]
    hx_sub = sub[sub.has_backbone]
    best_any = min((order[l] for l in sub.ref_label if l in order), default=9)
    best_hx = min((order[l] for l in hx_sub.ref_label if l in order), default=9)
    best_sel = min((order[l] for l, r in zip(sub.ref_label, sub.role) if l in order and r != "dropped"),
                   default=9)
    best_pri = min((order[l] for l, r in zip(sub.ref_label, sub.role) if l in order and r == "primary"),
                   default=9)
    ex = sub[sub.ref_label == "ref_cds_exact"]
    gr.append(dict(ref_gene=g, n_structures=len(sub), best_any=best_any, best_backbone=best_hx,
                   best_selected=best_sel, best_primary=best_pri,
                   exact_sources=";".join(ex.sources), exact_roles=";".join(ex.role),
                   exact_scores=";".join(f"{s}" for s in ex.score),
                   primary_score=";".join(f"{s}" for s in sub[sub.role == "primary"].score)))
pd.DataFrame(gr).to_csv(f"{OUT}_refgenes.tsv", sep="\t", index=False)
print(f"{WD}: {len(df)} structures, {len(pc_genes)} ref genes, "
      f"{(df.role == 'primary').sum()} primary, {(df.role == 'alternate').sum()} alternate")
