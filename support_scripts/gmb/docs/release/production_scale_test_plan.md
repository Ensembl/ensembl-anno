# First production-scale comparison — *Zymoseptoria tritici* (GCA_000219625.1)

Purpose: decide, on the whole genome, whether the experimental evidence-integrated fungal
configuration should become the default `fungi` behaviour. The decision criteria (§6) are fixed
here, before the run. Background: [`evidence_integration.md`](evidence_integration.md) (17-window
evaluation), [`handover_review.md`](handover_review.md) (readiness). How to run GMB in general:
[`../../HANDOVER.md`](../../HANDOVER.md).

---

## 1. Runs

| # | run | configuration | long reads | purpose |
|---|---|---|---|---|
| R1 | **Helixer alone** | none (remapped Helixer GFF3) | — | the backbone GMB must improve on |
| R2 | **Current fungal GMB** | `--preset fungi` | no | production baseline (configuration A) |
| R3 | **Experimental, no long reads** | `--preset fungi --config configs/fungi_evidence_integration.experimental.yaml` | no | effect of the selection changes (configuration B) |
| R4 | **Experimental + long reads** | as R3 | collapsed | B as intended for production |
| R5 | Current fungal GMB + long reads | as R2 | collapsed | separates the long-read contribution (R5 − R2) from the selection changes (R4 − R5) |

Optional R6: R4 with protein validation as it would ship. The 17-window evaluation found that the
fungi validation weights (`diamond_weight` 0.7, `min_score` 0.7, penalty 5.0) prune most
alternates and cost 15 exact CDS on held-out core windows (`handover_review.md` §2). Run it only
to measure that, not as a candidate default.

**Identical inputs in every run:** the same genome FASTA, remapped Helixer GFF3, OrthoDB and
UniProt genBlastG alignments, Scallop and StringTie assemblies, the same collapsed long-read
track (R4, R5), the same `--gene-prefix`, and the same GMB commit. Protein validation off
(R1–R5).

## 2. Inputs

```bash
EV=/path/Fungi_release_tests/input_data/GCA_000219625.1      # evidence (read-only)
REF=/path/Fungi_release_tests/reference_data                   # evaluation only
GMB=/path/ensembl-anno/support_scripts/gmb                     # checkout at the tested commit
OUT=/work/ztritici_gmb_$(git -C $GMB rev-parse --short HEAD)
mkdir -p $OUT/inputs
# genome: headers trimmed to the sequence name (1..21)
gzip -dc $REF/Zymoseptoria_tritici_GCA_000219625.1.fa.bgz | awk '/^>/{print $1;next}{print}' > $OUT/inputs/genome.fa
gzip -dc $REF/Zymoseptoria_tritici_GCA_000219625.1.gff3.bgz > $OUT/inputs/reference.gff3
# Helixer names sequences by GenBank accession: rename to the genome's names
python $GMB/tools/remap_helixer.py --input ${EV}_helixer.gff3 \
  --assembly-report ${EV}_MYCGR_v2.0_assembly_report.txt --output $OUT/inputs/helixer.gff3
# long reads: 2.05 M per-read alignments -> consensus transcript models
gmb-longread-consensus --input ${EV}_minimap2_annotation.gtf --preset fungi \
  --output-dir $OUT/inputs/longread
LR=$OUT/inputs/longread/minimap2_consensus.gtf
printf '%s\n' 1 2 3 4 5 6 7 8 9 10 11 12 13 > $OUT/inputs/core.txt
printf '%s\n' 14 15 16 17 18 19 20 21 > $OUT/inputs/accessory.txt
```

## 3. Commands

```bash
G=$OUT/inputs/genome.fa
EVID=(--helixer $OUT/inputs/helixer.gff3
      --scallop ${EV}_scallop_annotation.gtf --stringtie ${EV}_stringtie_annotation.gtf
      --orthodb ${EV}_orthodb_annotation.gtf --uniprot ${EV}_genblast_annotation.gtf)
X=$GMB/configs/fungi_evidence_integration.experimental.yaml

run() {   # run <name> [extra gmb args...]
  local n=$1; shift
  gmb-preflight --preset fungi "$@" --genome $G "${EVID[@]}" --output-dir $OUT/$n/preflight
  /usr/bin/time -v gmb-build --preset fungi "$@" --genome $G "${EVID[@]}" \
    --gene-prefix ZTGMB --output-dir $OUT/$n/build --validate-fasta 2> $OUT/$n/build.time
  gmb-finalise --build-dir $OUT/$n/build --genome $G --output-dir $OUT/$n/finalise
  git -C $GMB rev-parse HEAD > $OUT/$n/gmb_commit.txt
}
run R2_baseline
run R3_experimental        --config $X
run R4_experimental_lr     --config $X --minimap2 $LR
run R5_baseline_lr         --minimap2 $LR
```

`gmb-preflight` exits 1 on FAIL: stop and fix the input rather than adding `--allow-fail`.
`--minimap2` with the **raw** per-read
file must produce `longread_collapsed` in preflight and a `WARNING` in the build log; if you see
either, the wrong file was passed.

Provenance captured automatically per run: `build/resolved_config.yaml` (+ SHA-256),
`build/run_manifest.json` (versions, input hashes, runtime), `finalise/handover_manifest.json`
(output inventory, configuration source). Wall time and peak RSS: `build.time` (GNU time
`Elapsed` and `Maximum resident set size`; on macOS use `/usr/bin/time -l`).
Expected scale: an earlier whole-genome *Z. tritici* build took 29 h with protein validation on
(pre-2.0.0 code). Without it, each of the 17 evaluation windows (14.9 Mb in all, 38% of the
genome) built in 2–15 s with < 0.5 GB, but genome-scale cost has not been measured. Allow
generous limits (≥ 24 h, ≥ 16 GB) for the first run and record what it used. The expected
preflight warnings are 0.4% / 1.7% unstranded rows in Scallop / StringTie, OrthoDB IDs
reused across sequences, and very long OrthoDB spans and introns (the protein filter drops
those).

## 4. Comparison (ensembl-genes `annotation-qc`; nothing is added to GMB)

For every run, three evaluations against the reference, each genome-wide and for core and
accessory chromosomes separately:

```bash
cmpr() {  # cmpr <run> <query gff3> <label> [extra annotation-qc args]
  local r=$1 q=$2 l=$3; shift 3
  for scope in all core accessory; do
    reg=(); [ $scope != all ] && reg=(--regions-file $OUT/inputs/$scope.txt)
    annotation-qc pairwise-compare --query $q --reference $OUT/inputs/reference.gff3 \
      --genome $G --evaluation-mode protein_coding --reference-transcript-biotypes protein_coding \
      "${reg[@]}" "$@" --outdir $OUT/cmp/$r/$l.$scope
  done
}
cmpr R1_helixer $OUT/inputs/helixer.gff3 all_isoforms
for r in R2_baseline R3_experimental R4_experimental_lr R5_baseline_lr; do
  F=$OUT/$r/finalise
  cmpr $r $F/canonical/consensus.canonical_annotated.gff3 canonical --query-transcript-selection canonical
  cmpr $r $F/consensus.gff3 all_isoforms --evidence-attribution $OUT/$r/build/evidence_attribution.tsv
done
```

Canonical is the handover view and the basis of the decision. The full gene set
(`all_isoforms`) rewards extra isoforms, because any isoform may match. Report it, but never
decide on it alone.

## 5. Metrics to report

From `comparison_summary.json` (canonical and all-isoform; all / core / accessory):

| metric | `annotation-qc` field |
|---|---|
| exact CDS (coordinate-exact) | `sensitivity_cds.cds_coordinate_exact_count` / `total_reference_genes` |
| partial CDS agreement | `reference_classification_cds` `Partial_Match`; CDS structural exact (`sensitivity_cds.cds_exact_match_count`) |
| exact CDS intron chain | `sensitivity_cds.cds_intron_chain_recovered` / `cds_intron_chain_evaluated` |
| reference loci recovered / missed | `sensitivity.locus_detected_count`, `sensitivity.missed_count` |
| incorrect splits and merges | `split_merge.gene_split_count`, `gene_merge_count` |
| novel genes | `specificity.novel_consensus_count`, `novel_categories` |
| one-to-one exact-CDS F1 | `cds_exact_one_to_one.f1` |
| CDS-intron precision / recall | `intron_support.cds.query_introns_supported_rate` / `reference_introns_recovered_rate` |

Per gene, from `comparison_details.tsv` of two runs: reference genes gaining and losing an exact
canonical CDS (R3 vs R2, R4 vs R5, each vs R1). Report gains and losses, not only the net.

From GMB's own outputs (reference-free):

| metric | source |
|---|---|
| genes, transcripts, isoforms per gene | `finalise/consensus.gff3` |
| transcript-only loci | `evidence_attribution.tsv`: genes with `gene_has_backbone_model == False` |
| antisense transcript-only coding predictions | as above **and** `opposite_strand_cds_overlap == True` |
| support of new loci | for transcript-only loci: `evidence_sources` (≥ 2 sources / long reads), `protein_alignments_compatible > 0`, `protein_support_strength`, `orf_label` (complete vs partial) |
| promotions over the backbone | `selection_reason == transcript_junction_support` (R3, R4) |
| isoform exclusions | `summary.json` `isoform_cds_overlap_excluded`; `build/excluded_candidates.tsv` |
| canonical ≠ build primary | `canonical_selection_summary.json` `canonical_not_build_primary_by_tier` |
| QC | `finalise/fasta_qc_report.json`, `utr_qc_report.json` (must pass) |
| runtime, peak memory | `build.time`, `run_manifest.json` |

```python
import pandas as pd
a = pd.read_csv("build/evidence_attribution.tsv", sep="\t")
g = a.groupby("gene_id").agg(bb=("gene_has_backbone_model", "first"),
                             anti=("opposite_strand_cds_overlap", "any"),
                             prot=("protein_alignments_compatible", "max"),
                             multi=("evidence_sources", lambda s: s.str.contains(",").any()))
new = g[~g.bb]
print("transcript-only loci", len(new), "| antisense", int(new.anti.sum()),
      "| spliced-protein support", int((new.prot > 0).sum()), "| multi-source", int(new.multi.sum()))
```

Annotation-qc's own novel categories (`Novel_no_reference_overlap`,
`Overlaps_evaluated_reference_opposite_strand`, and `Strand_Mismatch` classes) give the
reference-based view of the same loci.

## 6. Decision criteria for promoting the experimental configuration

Fixed in advance. Canonical output, **core chromosomes**, R3 compared with R2 (and R4 with R5 for
the long-read version). Promote only if **all** hold:

1. **Exact CDS:** R3 ≥ R2 + 1.0 percentage point of core reference genes, and R3 ≥ R1. (The
   17-window held-out estimate was +2.0 points over R2 and +1.2 over R1.)
2. **Regressions bounded:** reference genes losing an exact canonical CDS relative to R1 ≤ ⅓ of
   those gaining one (held-out: 5 lost vs 28 gained).
3. **Exact intron chains:** R3 ≥ R2.
4. **Structure:** canonical splits and merges ≤ R2. Full-gene-set splits + merges ≤ 1.2 × R1 (the
   overlay's isoform rule should remove most of R2's isoform-driven merges).
5. **Loci:** missed reference genes ≤ R2; total gene count ≤ R2 × 1.02.
6. **Accessory chromosomes:** exact CDS not below R2 by more than 1% of accessory reference
   genes.
7. **QC and cost:** finalisation QC passes; wall time ≤ 1.25 × R2; peak RSS within the job limit.

If 1–7 hold, set the four overlay keys in `gmb/configs/fungi.yaml` (and regenerate the golden
fixture and its tests). If they fail, keep configuration A as the default and record which
criterion failed.

**Long reads** are recommended whenever R4 satisfies 1–7 against R5 *and* R4 ≥ R3 on criteria 1
and 3. That makes collapsing the long-read track a pipeline requirement, upstream of GMB.

**Not a criterion, but must be reported:** the number of transcript-only antisense coding
predictions (R2–R5) and their support. A decision on filtering them (e.g. a coding-potential
gate) is separate and needs its own validation (`../known_issues.md`).

## 7. Expected outputs to keep

Per run: `finalise/` (handover), `build/{evidence_attribution.tsv,excluded_candidates.tsv,
resolved_config.yaml,resolved_config_sha256,run_manifest.json,summary.json}`, `preflight/`,
`build.time`, `gmb_commit.txt`; `cmp/<run>/<view>.<scope>/comparison_summary.{json,tsv}`,
`comparison_details.tsv`, `gene_splits.tsv`, `gene_merges.tsv`. Do not keep `build/` FASTA or
GFF3 copies: they duplicate `finalise/`.
