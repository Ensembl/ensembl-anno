# GMB BUSCO Triage Report
## Why GMB Misses 1483 BUSCOs Present in Helixer

**Species:** *Zymoseptoria tritici* IPO323  
**BUSCO lineage:** mycosphaerellaceae_odb12 (n = 5686)  
**GMB output:** `output_full/consensus.gff3`, `prot.fa`  
**Date:** 2026-04-23

---

## Summary numbers

| Tool     | Complete | Duplicated | Fragmented | Missing |
|----------|----------|------------|------------|---------|
| Helixer  | 5651     | 3          | 13         | 19      |
| GMB      | 3222     | 895        | 69         | 1500    |
| Anno v1  | (see files) | – | – | – |
| JGI      | (see files) | – | – | – |

**Set M** (GMB Missing ∩ Helixer present) = **1483 BUSCOs**

---

## Task A — Missing BUSCO triage

### Helixer-status breakdown of the 1483

| Helixer status | Count |
|----------------|-------|
| Complete       | 1478  |
| Fragmented     | 5     |

→ `missing_buscos_present_in_helixer.tsv`

### GMB duplication analysis (895 Duplicated BUSCO groups)

| Category | Count | % |
|----------|-------|---|
| Within-gene isoforms (same gene, multiple transcripts) | 456 | 50.9% |
| Across different genes (true paralogues or redundant loci) | 439 | 49.1% |

The ~49% across-gene duplication is elevated and warrants follow-up; some of these are likely legitimate gene-family expansions in *Z. tritici*, but others may reflect gene-model fragmentation or mis-assembly.

---

## Task B — Root cause breakdown

→ `busco_missing_root_cause_breakdown.tsv`  
(includes `refined_classification` and `note` columns)

| Class | Count | % | Meaning |
|-------|-------|---|---------|
| **Spatial_overlap_wrong_gene** | 591 | 39.9% | A GMB gene overlaps the Helixer locus spatially but is a *different* gene; the actual Helixer BUSCO gene was dropped from the GMB proteome |
| **Replaced_transcriptome_model** | 434 | 29.3% | GMB chose a Scallop/StringTie model at this locus instead of Helixer |
| **Helixer_kept_protein_bad_sequence** | 401 | 27.0% | GMB kept the Helixer exon structure but the translated protein has the wrong reading frame → premature stop codons, BUSCO cannot find the gene |
| **Dropped_locus_no_model** | 40 | 2.7% | No GMB model at the Helixer locus at all |
| **Helixer_kept_protein_truncated** | 17 | 1.1% | GMB kept the Helixer model but the protein is < 50% expected length |

---

## Task C — Root causes in detail

### C1 — Wrong protein (bad_sequence + truncated, ~28%) ← highest priority

**23.6% of all 15399 GMB proteins contain internal stop codons**  
(3629/15399 total; 82.5% of bad_sequence Helixer models; 0% of transcriptome-only models).

**Root cause: bug in `_build_cds_seq_from_intervals` (fasta_export.py) and the
inline CDS-building block in `annotate_cds_utrs.py`.**

For **minus-strand multi-exon CDS**, both sites reverse-sort the exon intervals
and then call `reverse_complement()` on the *concatenated* sequence. Because
`RC(X + Y) = RC(Y) + RC(X)`, this inverts the exon order in the final CDS
string — placing the 3'-most exon first — producing a completely wrong protein
with multiple internal stop codons.

Confirmation: for `ZTRITICI_00392.t1` (same exon coordinates as Helixer gene
`_CM001196.1_002208.1`, BUSCO `10004at93133`):
- Helixer protein: `MSASSLKYDFFNVSIPVEYVAVVEVDRVKK…` (correct)
- GMB (buggy): `DWTMSLFGMQQCYKRKT*RTHYSRRKRRGRR…` (9 internal stops)
- Fixed: `MSASSLKYDFFNVSIPVEYVAVVEVDRVKK…` (matches Helixer exactly)

**Fix applied in this branch:**
- `fasta_export.py`: `_build_cds_seq_from_intervals` — sort ascending, not reversed
- `annotate_cds_utrs.py`: inline CDS loop at line ~585 — same fix
- **Regression test added:** `tests/test_fasta_export.py::TestHelixerCdsExport::test_helixer_cds_minus_strand_multiexon`

Expected impact: fixing this should convert ≈400–600 Helixer BUSCO hits from
Missing to Complete in GMB.

### C2 — Replaced by transcriptome model (434, 29.3%)

GMB selects a Scallop or StringTie assembly model at the locus in preference to
the Helixer model. UTR support analysis:

| UTR status | Count | % of 434 |
|------------|-------|----------|
| 5′ UTR dropped (no end agreement) | 205 | 47.2% |
| 3′ UTR dropped (no end agreement) | 223 | 51.4% |

Nearly half of the winning transcriptome models lack 5′/3′ end support from
other transcriptome data, which means they are likely truncated or fragmented
assemblies. GMB is preferring these over Helixer despite their lower coding
confidence.

**Recommendation:**  
At any locus where Helixer provides a candidate, only displace the Helixer model
with a transcriptome model if the transcriptome model satisfies **both**:
1. Has supported 5′ **and** 3′ ends (or `protein_coding_score` ≥ Helixer's score)
2. `protein_coding_score` delta ≤ 0.05 below the Helixer model

Avoid using UTR end-support failure as a disqualification signal against Helixer
models specifically — Helixer does not depend on polyA or TSS signals and its
UTRs will almost never have transcriptome endpoint agreement.

### C3 — Dropped / wrong-gene loci (591 + 40, 42.6%)

The Helixer gene is absent from the GMB proteome. The `Spatial_overlap_wrong_gene`
cases (591) have a different GMB gene at the same genomic coordinates — the Helixer
gene was lost during:
- protein competition/deduplication (a neighbouring gene outcompeted it), or
- `dedup_genes` merging (984 merges reported in `summary.json`)

The 40 `Dropped_locus_no_model` cases had no model at all.

**Recommendation:**  
Investigate whether these 591 loci are genuinely adjacent genes (no fix needed)
or whether the Helixer model was incorrectly merged into a neighbouring locus.
A post-fix BUSCO run (after C1 is resolved) will clarify how many of these
remain as true dropped loci vs. loci that now have correct proteins.

---

## Task D — Canonical proteome

→ `gmb_canonical_proteins.fa`  
**13587 genes**, one protein each (longest isoform; Helixer-source preferred on ties).

| Stat | Value |
|------|-------|
| Sequences | 13587 |
| Min length | 4 aa |
| Median length | 344 aa |
| Max length | 6124 aa |

**Caution:** 23.6% of proteins in the current GMB FASTA have internal stop
codons (see C1). After the pipeline bug is fixed and the run is repeated, this
canonical proteome should be regenerated from the corrected `prot.fa`.

---

## Priority action list

| # | Priority | Action | Expected BUSCO gain |
|---|----------|--------|---------------------|
| 1 | **Critical** | Apply minus-strand CDS fix to `fasta_export.py` + `annotate_cds_utrs.py` and re-run GMB | +400–600 Complete |
| 2 | High | Tighten Helixer-vs-transcriptome selection: require full end support *or* protein_coding_score ≥ Helixer before displacing Helixer model | +300–400 Complete |
| 3 | Medium | Investigate spatial_overlap_wrong_gene cases post-fix to identify true merging errors | +0–100 Complete |
| 4 | Low | Investigate KI-006 single-exon asymmetry (may affect a small number of Helixer single-exon genes) | +minimal |

---

## Output files

| File | Description |
|------|-------------|
| `missing_buscos_present_in_helixer.tsv` | All 1483 BUSCO IDs with Helixer status and sequence IDs |
| `busco_missing_root_cause_breakdown.tsv` | Per-BUSCO mapping + refined classification |
| `gmb_canonical_proteins.fa` | One protein per GMB gene (longest isoform) |
| `busco_triage.py` | Full analysis script (re-runnable) |
