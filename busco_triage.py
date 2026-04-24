#!/usr/bin/env python3
"""
BUSCO triage: why does GMB miss 1483 BUSCOs present in Helixer?
Tasks A-D per specification.
"""

import sys
import re
import gzip
from collections import defaultdict
from pathlib import Path

# ── paths ────────────────────────────────────────────────────────────────────
# Script lives in the worktree; data lives in the parent gmb directory
WORKTREE = Path(__file__).parent
BASE     = WORKTREE.parents[2]   # .claude/worktrees/NAME -> gmb
ANNOT    = BASE / "z_tritici/annotations"
OUTPUT   = BASE / "output_full"

GMB_BUSCO    = ANNOT / "busco_core_protein_mode_gmb/run_mycosphaerellaceae_odb12/full_table.tsv"
HEL_BUSCO    = ANNOT / "busco_core_protein_mode_helixer/run_mycosphaerellaceae_odb12/full_table.tsv"
HEL_GFF      = ANNOT / "helixer.gff3"
GMB_GFF      = OUTPUT / "consensus.gff3"
EVIDENCE_TSV = OUTPUT / "evidence_attribution.tsv"
GMB_PROT     = ANNOT / "gmb_prot.fa"        # protein FASTA for Task D
ASSEMBLY_REP = BASE  / "z_tritici/GCF_000219625.1_MYCGR_v2.0_assembly_report.txt"

OUT_DIR = WORKTREE  # write outputs into the worktree

# ── chromosome mapping (CM... -> integer) ────────────────────────────────────
def build_chr_map(assembly_report: Path) -> dict:
    """Return {CM_id: int_str} e.g. CM001196.1 -> '1'"""
    mapping = {}
    with open(assembly_report) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) < 7:
                continue
            # col 4 = GenBank accession, col 2 = chromosome number
            genbank = parts[4]
            chrom   = parts[2]
            mapping[genbank] = chrom
    return mapping


# ── BUSCO table parser ────────────────────────────────────────────────────────
def parse_busco_table(path: Path) -> dict:
    """
    Returns {busco_id: {'status': str, 'hits': [{'seq': str, 'score': float, 'length': int}]}}
    Duplicated BUSCOs appear multiple times; we collect all hits.
    """
    result = {}
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                continue
            bid    = parts[0]
            status = parts[1]
            hit = {}
            if len(parts) >= 5 and parts[2]:
                hit = {"seq": parts[2], "score": float(parts[3]), "length": int(parts[4])}
            if bid not in result:
                result[bid] = {"status": status, "hits": []}
            # Duplicated entries: status is Duplicated for every row
            if status != "Missing":
                result[bid]["status"] = status  # last write wins (all same for Complete/Frag)
            if hit:
                result[bid]["hits"].append(hit)
    return result


# ── GFF3 parser helpers ───────────────────────────────────────────────────────
def parse_attrs(attr_str: str) -> dict:
    d = {}
    for tok in attr_str.split(";"):
        tok = tok.strip()
        if "=" in tok:
            k, v = tok.split("=", 1)
            d[k] = v
    return d


def parse_mrna_coords(gff_path: Path, feature_type="mRNA") -> dict:
    """
    Returns {transcript_id: {'chrom': str, 'start': int, 'end': int, 'strand': str,
                              'gene_id': str}}
    """
    records = {}
    open_fn = gzip.open if str(gff_path).endswith(".gz") else open
    with open_fn(gff_path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            if parts[2] != feature_type:
                continue
            chrom  = parts[0]
            start  = int(parts[3])
            end    = int(parts[4])
            strand = parts[6]
            attrs  = parse_attrs(parts[8])
            tid    = attrs.get("ID", "")
            parent = attrs.get("Parent", "")
            records[tid] = {"chrom": chrom, "start": start, "end": end,
                             "strand": strand, "gene_id": parent,
                             "evidence": attrs.get("Evidence", "")}
    return records


# ── interval index for overlap queries ───────────────────────────────────────
class IntervalIndex:
    """Bin-based overlap index. Bin size 10 kbp."""
    BIN = 10_000

    def __init__(self):
        self._bins = defaultdict(list)   # (chrom, bin) -> [(start,end,data)]

    def add(self, chrom, start, end, data):
        b0 = start // self.BIN
        b1 = end   // self.BIN
        for b in range(b0, b1 + 1):
            self._bins[(chrom, b)].append((start, end, data))

    def query(self, chrom, start, end):
        b0 = start // self.BIN
        b1 = end   // self.BIN
        seen = set()
        out  = []
        for b in range(b0, b1 + 1):
            for rec in self._bins.get((chrom, b), []):
                rid = id(rec)
                if rid in seen:
                    continue
                seen.add(rid)
                rs, re, data = rec
                if rs <= end and re >= start:
                    out.append(data)
        return out


# ── Task A ────────────────────────────────────────────────────────────────────
def task_a(gmb_busco: dict, hel_busco: dict, out_dir: Path):
    print("\n=== TASK A: Missing BUSCO triage ===")

    present_statuses = {"Complete", "Duplicated", "Fragmented"}
    M = {
        bid for bid, rec in gmb_busco.items()
        if rec["status"] == "Missing"
        and hel_busco.get(bid, {}).get("status", "Missing") in present_statuses
    }
    print(f"|M| (GMB Missing, Helixer present) = {len(M)}")

    # Write TSV
    out_path = out_dir / "missing_buscos_present_in_helixer.tsv"
    with open(out_path, "w") as fh:
        fh.write("busco_id\thelixer_status\thelixer_sequence_ids\thelixer_scores\thelixer_lengths\n")
        for bid in sorted(M):
            hrec = hel_busco[bid]
            seqs   = ",".join(h["seq"]             for h in hrec["hits"]) or "."
            scores = ",".join(str(h["score"])      for h in hrec["hits"]) or "."
            lens   = ",".join(str(h["length"])     for h in hrec["hits"]) or "."
            fh.write(f"{bid}\t{hrec['status']}\t{seqs}\t{scores}\t{lens}\n")
    print(f"  -> {out_path.name}")

    # GMB duplication analysis
    gmb_dup_ids = [bid for bid, r in gmb_busco.items() if r["status"] == "Duplicated"]
    print(f"\nGMB Duplicated BUSCOs: {len(gmb_dup_ids)}")
    # check within-gene vs across-gene duplication
    within_gene  = 0
    across_gene  = 0
    for bid in gmb_dup_ids:
        hits  = gmb_busco[bid]["hits"]
        genes = set()
        for h in hits:
            # ZTRITICI_00009.t2 -> gene ZTRITICI_00009
            genes.add(h["seq"].rsplit(".", 1)[0])
        if len(genes) == 1:
            within_gene += 1
        else:
            across_gene += 1
    print(f"  Within-gene isoforms only: {within_gene}")
    print(f"  Across different genes:    {across_gene}")

    # Status distribution summary
    status_counts = defaultdict(int)
    for rec in gmb_busco.values():
        status_counts[rec["status"]] += 1
    print(f"\nGMB BUSCO status distribution:")
    for st, cnt in sorted(status_counts.items(), key=lambda x: -x[1]):
        print(f"  {st:15s}: {cnt}")

    hel_status_counts = defaultdict(int)
    for rec in hel_busco.values():
        hel_status_counts[rec["status"]] += 1
    print(f"\nHelixer BUSCO status distribution:")
    for st, cnt in sorted(hel_status_counts.items(), key=lambda x: -x[1]):
        print(f"  {st:15s}: {cnt}")

    return M


# ── Task B ────────────────────────────────────────────────────────────────────
def task_b(M: set, hel_busco: dict, hel_mrnas: dict, gmb_index: IntervalIndex,
           chr_map: dict, gmb_mrnas: dict, evidence_attr: dict, out_dir: Path):
    print("\n=== TASK B: Map Helixer loci → GMB overlap ===")

    # For each BUSCO in M, get Helixer mRNA coords and query GMB index
    rows = []
    for bid in sorted(M):
        hrec = hel_busco[bid]
        for hit in hrec["hits"]:
            seq_id = hit["seq"]
            hmrna  = hel_mrnas.get(seq_id)
            if hmrna is None:
                rows.append({
                    "busco_id": bid,
                    "helixer_seq": seq_id,
                    "helixer_chrom_cm": "?",
                    "helixer_start": "?",
                    "helixer_end": "?",
                    "gmb_overlapping": ".",
                    "classification": "Helixer_seq_not_in_gff",
                    "gmb_evidence": ".",
                })
                continue
            cm_chrom = hmrna["chrom"]
            int_chrom = chr_map.get(cm_chrom, cm_chrom)
            h_start  = hmrna["start"]
            h_end    = hmrna["end"]

            overlaps = gmb_index.query(int_chrom, h_start, h_end)
            if not overlaps:
                classification = "Dropped_locus"
                gmb_ids   = "."
                gmb_evids = "."
            else:
                gmb_ids   = ",".join(o["tid"] for o in overlaps)
                gmb_evids = ",".join(o.get("evidence", ".") or "." for o in overlaps)

                # Check if any overlap is from Helixer source
                all_evid = [o.get("evidence", "") or "" for o in overlaps]
                has_helixer = any("Helixer" in e for e in all_evid)
                has_transcriptome = any(
                    any(s in e for s in ["Scallop","StringTie","scallop","stringtie"])
                    for e in all_evid
                )

                # Check protein lengths of overlapping GMB models
                short_proteins = []
                for o in overlaps:
                    plen = o.get("prot_len", 0)
                    if plen and plen < hit["length"] * 0.5:
                        short_proteins.append(o["tid"])

                if has_helixer:
                    classification = "Present_Helixer_model"
                elif has_transcriptome:
                    classification = "Replaced_transcriptome_model"
                else:
                    classification = "Replaced_other_model"

            rows.append({
                "busco_id": bid,
                "helixer_seq": seq_id,
                "helixer_chrom_cm": cm_chrom,
                "helixer_chrom_int": int_chrom,
                "helixer_start": h_start,
                "helixer_end": h_end,
                "helixer_busco_score": hit["score"],
                "helixer_busco_length": hit["length"],
                "gmb_overlapping": gmb_ids,
                "gmb_evidence": gmb_evids,
                "classification": classification,
            })

    # Summary counts
    class_counts = defaultdict(int)
    for r in rows:
        class_counts[r["classification"]] += 1

    print(f"  Total BUSCO-hit rows analysed: {len(rows)}")
    print(f"  Classification breakdown:")
    for cls, cnt in sorted(class_counts.items(), key=lambda x: -x[1]):
        pct = 100 * cnt / max(len(rows), 1)
        print(f"    {cls:40s}: {cnt:5d}  ({pct:.1f}%)")

    out_path = out_dir / "busco_missing_root_cause_breakdown.tsv"
    fieldnames = [
        "busco_id","helixer_seq","helixer_chrom_cm","helixer_chrom_int",
        "helixer_start","helixer_end","helixer_busco_score","helixer_busco_length",
        "gmb_overlapping","gmb_evidence","classification"
    ]
    with open(out_path, "w") as fh:
        fh.write("\t".join(fieldnames) + "\n")
        for r in rows:
            fh.write("\t".join(str(r.get(f, ".")) for f in fieldnames) + "\n")
    print(f"  -> {out_path.name}")

    return rows, class_counts


# ── Task C ────────────────────────────────────────────────────────────────────
def task_c(rows: list, evidence_attr: dict, class_counts: dict):
    print("\n=== TASK C: Why is Helixer losing at replaced loci? ===")

    replaced_rows = [r for r in rows if "Replaced" in r["classification"]]
    print(f"  Replaced loci rows: {len(replaced_rows)}")
    if not replaced_rows:
        print("  (no replaced loci - dropped locus dominates)")
        return

    # Analyse evidence sources at replaced loci
    evid_counter = defaultdict(int)
    for r in replaced_rows:
        evids = r["gmb_evidence"].split(",")
        for e in evids:
            evid_counter[e.strip()] += 1

    print("  Evidence at replaced loci (GMB winner sources):")
    for e, cnt in sorted(evid_counter.items(), key=lambda x: -x[1]):
        print(f"    {e:40s}: {cnt}")

    # UTR support analysis from evidence_attribution
    utr5_no_agree  = 0
    utr3_no_agree  = 0
    total_replaced = 0
    for r in replaced_rows:
        if r["gmb_overlapping"] == ".":
            continue
        for tid in r["gmb_overlapping"].split(","):
            if tid in evidence_attr:
                ea = evidence_attr[tid]
                total_replaced += 1
                if ea.get("utr_5p_action") == "dropped":
                    utr5_no_agree += 1
                if ea.get("utr_3p_action") == "dropped":
                    utr3_no_agree += 1

    if total_replaced:
        print(f"\n  UTR support at replaced loci ({total_replaced} transcripts checked):")
        print(f"    5' UTR dropped (no end agreement): {utr5_no_agree} ({100*utr5_no_agree/total_replaced:.1f}%)")
        print(f"    3' UTR dropped (no end agreement): {utr3_no_agree} ({100*utr3_no_agree/total_replaced:.1f}%)")

    total_loci = sum(class_counts.values())
    dropped    = class_counts.get("Dropped_locus", 0)
    replaced   = sum(v for k, v in class_counts.items() if "Replaced" in k)
    print(f"\n  --- Root cause summary ---")
    print(f"  Dropped locus   (no GMB model at Helixer coords): {dropped:5d}  ({100*dropped/max(total_loci,1):.1f}%)")
    print(f"  Replaced locus  (GMB chose non-Helixer model):    {replaced:5d}  ({100*replaced/max(total_loci,1):.1f}%)")
    print(f"  Present (Helixer model kept, BUSCO still missing):{class_counts.get('Present_Helixer_model',0):5d}")

    print("""
  ── Recommendations ──
  PRIMARY ISSUE (Dropped locus):
    If >50% are dropped loci, the Helixer gene was not emitted in consensus.gff3.
    Possible causes:
      1. Locus fell below protein_coding_score threshold after competition.
      2. Helixer model was outcompeted by a transcriptome model that covers only
         part of the locus (partial overlap).
      3. Helixer model was removed by dedup_genes (merged into neighbour).
    Action: Lower Helixer minimum protein_coding_score threshold OR add a
    "BUSCO-protection" pass that forces inclusion of high-scoring Helixer models
    with no transcriptome competition.

  SECONDARY ISSUE (Replaced locus):
    GMB selected a transcriptome-derived model at the same locus, displacing
    the Helixer model. Often the transcriptome model is truncated or multi-exon
    assembly artefact lacking strong coding signal.
    Action: When a locus has Helixer AND transcriptome candidates, gate the
    decision on protein_coding_score delta: retain Helixer unless transcriptome
    model has score >= helixer_score - MARGIN (suggest MARGIN=0.05).
    Also check that end-support logic (utr_5p/3p_action) is not penalising
    Helixer models for missing UTR evidence they structurally lack.
""")


# ── Task D ────────────────────────────────────────────────────────────────────
def task_d(gmb_prot_fa: Path, evidence_attr: dict, out_dir: Path):
    print("\n=== TASK D: Canonical proteome (one protein per gene) ===")

    # Parse protein FASTA - header like >ZTRITICI_00001.t1
    seqs = {}
    current_id = None
    current_seq = []
    if not gmb_prot_fa.exists():
        print(f"  WARNING: {gmb_prot_fa} not found - skipping Task D")
        return

    with open(gmb_prot_fa) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                if current_id:
                    seqs[current_id] = "".join(current_seq)
                current_id = line[1:].split()[0]
                current_seq = []
            else:
                current_seq.append(line)
    if current_id:
        seqs[current_id] = "".join(current_seq)

    print(f"  Total protein sequences: {len(seqs)}")

    # Group by gene
    gene_transcripts = defaultdict(list)
    for tid, seq in seqs.items():
        gene = tid.rsplit(".", 1)[0]
        gene_transcripts[gene].append((tid, seq))

    # For each gene, pick canonical = longest protein
    # Tie-break: prefer Helixer evidence source (from evidence_attr)
    canonical = {}
    for gene, transcripts in gene_transcripts.items():
        # Sort by length desc, then prefer Helixer
        def sort_key(item):
            tid, seq = item
            ea = evidence_attr.get(tid, {})
            evid = ea.get("evidence_sources", "")
            is_helixer = 1 if "Helixer" in evid else 0
            return (len(seq), is_helixer)
        best_tid, best_seq = max(transcripts, key=sort_key)
        canonical[gene] = (best_tid, best_seq)

    out_path = out_dir / "gmb_canonical_proteins.fa"
    with open(out_path, "w") as fh:
        for gene in sorted(canonical):
            tid, seq = canonical[gene]
            fh.write(f">{gene} canonical={tid}\n")
            for i in range(0, len(seq), 60):
                fh.write(seq[i:i+60] + "\n")

    print(f"  Canonical genes written: {len(canonical)}")
    print(f"  -> {out_path.name}")

    # Stats
    lens = [len(seq) for _, seq in canonical.values()]
    print(f"  Protein length stats: min={min(lens)} median={sorted(lens)[len(lens)//2]} max={max(lens)}")


# ── evidence_attribution loader ───────────────────────────────────────────────
def load_evidence_attr(path: Path) -> dict:
    """Returns {transcript_id: {col: val}}"""
    result = {}
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            rec = dict(zip(header, parts))
            tid = rec.get("transcript_id", "")
            if tid:
                result[tid] = rec
    return result


# ── main ──────────────────────────────────────────────────────────────────────
def main():
    print("Loading chromosome mapping...")
    chr_map = build_chr_map(ASSEMBLY_REP)
    print(f"  {len(chr_map)} chromosomes mapped (e.g. CM001196.1 -> {chr_map.get('CM001196.1','?')})")

    print("Loading BUSCO tables...")
    gmb_busco = parse_busco_table(GMB_BUSCO)
    hel_busco = parse_busco_table(HEL_BUSCO)
    print(f"  GMB:    {len(gmb_busco)} BUSCO entries")
    print(f"  Helixer:{len(hel_busco)} BUSCO entries")

    print("Loading evidence attribution...")
    evidence_attr = load_evidence_attr(EVIDENCE_TSV)
    print(f"  {len(evidence_attr)} transcript records")

    print("Parsing Helixer GFF3 mRNA coordinates...")
    hel_mrnas = parse_mrna_coords(HEL_GFF, feature_type="mRNA")
    print(f"  {len(hel_mrnas)} Helixer mRNAs")

    print("Parsing GMB consensus GFF3 and building interval index...")
    gmb_mrnas_raw = {}
    open_fn = gzip.open if str(GMB_GFF).endswith(".gz") else open
    with open_fn(GMB_GFF, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9 or parts[2] != "mRNA":
                continue
            chrom  = parts[0]
            start  = int(parts[3])
            end    = int(parts[4])
            strand = parts[6]
            attrs  = parse_attrs(parts[8])
            tid    = attrs.get("ID", "")
            parent = attrs.get("Parent", "")
            evid   = attrs.get("Evidence", "")
            gmb_mrnas_raw[tid] = {"chrom": chrom, "start": start, "end": end,
                                   "strand": strand, "gene_id": parent,
                                   "evidence": evid, "tid": tid}

    gmb_index = IntervalIndex()
    for tid, rec in gmb_mrnas_raw.items():
        gmb_index.add(rec["chrom"], rec["start"], rec["end"], rec)
    print(f"  {len(gmb_mrnas_raw)} GMB mRNAs indexed")

    # ── run tasks ─────────────────────────────────────────────────────────────
    M = task_a(gmb_busco, hel_busco, OUT_DIR)

    rows, class_counts = task_b(M, hel_busco, hel_mrnas, gmb_index,
                                chr_map, gmb_mrnas_raw, evidence_attr, OUT_DIR)

    task_c(rows, evidence_attr, class_counts)

    task_d(GMB_PROT, evidence_attr, OUT_DIR)

    print("\nDone. Output files written to:", OUT_DIR)


if __name__ == "__main__":
    main()
