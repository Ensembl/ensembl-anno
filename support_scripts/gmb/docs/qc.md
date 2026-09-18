# GMB QC

Two categories, deliberately kept apart.

---

## 1. Correctness QC — these must pass

Structural and sequence invariants. A violation is a **defect in the annotation**, not a
biological judgement. **A build that fails any of these is not a handover candidate.**

| check | requirement | where |
|---|---|---|
| cDNA sequence mismatches | **0** | `finalise/fasta_qc_report.json` |
| CDS sequence mismatches | **0** | same |
| protein sequence mismatches | **0** | same |
| unintended internal stop codons | **0** | same |
| gene-boundary violations | **0** | gene span == union of its transcripts |
| UTR invariant violations | **0** | `finalise/utr_qc_report.json` |
| one canonical transcript per gene | **PASS** | `finalise/canonical/` |
| cross-seqid chimeras | **0** | no transcript whose parent gene is on another sequence |

`gmb-build --validate-fasta` exits non-zero when sequence QC fails.

> **A non-zero exit here is a QC verdict, not a crash.** The annotation is fully written, and
> the run manifest is deliberately written *before* the QC exit so a failing run still leaves a
> complete record of what produced it.

### Why each one is checked

Each corresponds to a defect that actually shipped and was not caught by anything else.

- **Sequence mismatches.** FASTA was once captured mid-pipeline while the GFF3 continued to be
  trimmed afterwards: **315 of 5,786 cDNAs (5.4%)** silently disagreed with the annotation, and
  QC reported `pass: true` because it only checked ID coverage. Sequences are now regenerated
  from the final GFF3, and a mismatch fails the verdict.
- **Internal stops.** Seven transcripts once carried internal stops because the same
  `transcript_id` appeared on two sequences and the loader fused them into chimeras — three of
  them across *opposite strands* — with a fabricated 40–96 kb intron. IDs are now namespaced by
  seqid.
- **Gene-boundary violations.** Gene records were not recomputed after validation dropped
  transcripts, leaving gene spans far wider than their own transcripts (one 492 kb "gene" whose
  transcripts were all under 8 kb). Affected 10.7% of P. falciparum genes, 22.1% of a fungal
  region and 27% of T. gondii. It also **inflated locus-detection metrics** in every report,
  because `gmb-compare` pairs genes by span overlap.
- **UTR invariants.** UTR trimming once kept the wrong end, leaving exonic bases classified as
  neither CDS nor UTR: 156 of 874 junctions (17.9%).

---

## 2. Biological quality diagnostics — reported, not enforced

These describe how *good* the annotation is. They are judgement calls and depend on the
organism, so they are reported rather than failed on.

| diagnostic | how to read it |
|---|---|
| canonical GT-AG intron fraction | compare with the reference where available. A eukaryote with normal splicing should be near 99% |
| gene count | against the expected gene count for the clade |
| single vs multi-CDS-exon distribution | compare with the reference; under- and over-splitting are both common failure modes |
| very long genes / introns | implausible spans are often chimeras |
| source composition | which evidence produced the selected models |
| merge / split indicators | one GMB gene spanning several reference genes, or vice versa |
| protein support distribution | strong / positional-only / none |

---

## The distinction that matters

> **QC PASS does not mean high biological accuracy.**

All four gene sets in the cross-species validation passed **every** hard QC check. Their
quality was not remotely equal:

| gene set | hard QC | biological quality |
|---|---|---|
| P. falciparum frozen, no long read | PASS | good — recommended |
| P. falciparum frozen + long read | PASS | good, but ~5.8% non-canonical introns from a bad input |
| Z. tritici baseline | PASS | usable; over-predicts by 48% |
| **Z. tritici frozen policy** | **PASS** | **worse than the baseline it was meant to improve** |

The last row is the point. It passed every correctness check while being **measurably worse**
than doing nothing (CDS exact 4,895 → 4,611; improvements : regressions = 249 : 555). Marking
a build collaborator-ready because FASTA QC passed would have shipped the worst of the four.

Correctness QC tells you the annotation is **internally consistent**. It cannot tell you the
gene models are **right**.

---

## Evaluating against a reference

`gmb-compare` is **evaluation tooling and not part of the production path**. Run it when a
reference exists.

```bash
gmb-compare --query "$OUT/finalise/consensus.gff3" \
            --reference "$REFERENCE" --reference-fasta "$GENOME" \
            --evaluation-mode protein_coding \
            --evidence-attribution "$OUT/build/evidence_attribution.tsv" \
            --output-dir "$OUT/comparison"
```

Add `--seqname-map` only when the reference uses different sequence names.

### Which metric to quote

| metric | use |
|---|---|
| **CDS exact** | **the headline.** Reference genes whose coding structure is reproduced exactly |
| **CDS exact — multi-exon** | **the most informative single number** — where evidence integration does real work |
| CDS exact — single-exon | mostly measures whether an unspliced call was kept |
| CDS intron chain | correct splice junctions, CDS ends may differ |
| Exact Match | identical intron chain **and** ≥0.8 reciprocal exonic overlap |
| locus detection / Missed | whether a gene was found at all |
| merges / splits | boundary quality |

> **Exact Match is not coding accuracy.** It includes UTR extent, and GMB emits far fewer and
> shorter UTRs than a curated reference. In one P. falciparum build, **1,206 reference genes
> (22.7%) had a byte-identical CDS and a matching intron chain yet failed Exact Match on UTR
> extent alone**; Exact Match was 20.3% where the reference had no UTR against 7.6% where it
> did — a 2.7× swing from annotation extent, not from accuracy.

> **Locus detection is only meaningful with correct gene records.** `gmb-compare` pairs genes
> by span overlap, so a gene wider than its own transcripts "detects" reference genes it has
> nothing near. Any locus-detection figure from a build predating the gene-boundary fix is
> inflated. CDS exact and Exact Match are unaffected.

### Complementary gene-content QC

BUSCO and OMArk on the final protein set measure gene *content*; the reference comparison
measures gene *structure*. Neither substitutes for the other. Record the lineage and its
release — BUSCO scores are not comparable across lineage versions.
