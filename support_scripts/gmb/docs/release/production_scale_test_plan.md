# First production-scale comparison — *Zymoseptoria tritici* (GCA_000219625.1)

Purpose: confirm on the whole genome that the latest fungal settings (`fungi_latest_settings`,
formerly the "experimental evidence-integrated configuration") should replace the historical
`fungi` preset. The decision criteria (§6) are fixed here, before the run. This plan belongs to
the annotation lead: the person running GMB produces R4 by following
[`../../RUN_FUNGI.md`](../../RUN_FUNGI.md) and returns it; the other arms are comparison runs.
Background: [`evidence_integration.md`](evidence_integration.md) (17-window evaluation),
[`handover_review.md`](handover_review.md).

---

## 1. Runs

| # | run | configuration | long reads | purpose |
|---|---|---|---|---|
| R1 | **Helixer alone** | none (remapped Helixer GFF3) | — | the backbone GMB must improve on |
| R2 | **Historical fungal GMB** | `--preset fungi` | no | previous baseline (configuration A) |
| R3 | **Latest settings, no long reads** | `--preset fungi_latest_settings` | no | effect of the selection changes (configuration B) |
| R4 | **Latest settings + long reads** | as R3 | per-read, collapsed automatically | **the production annotation** (RUN_FUNGI.md) |
| R5 | Historical fungal GMB + long reads | as R2 | per-read, collapsed automatically | separates the long-read contribution (R5 − R2) from the selection changes (R4 − R5) |

Optional R6: R4 with protein validation as it would ship. The 17-window evaluation found that the
fungi validation weights (`diamond_weight` 0.7, `min_score` 0.7, penalty 5.0) prune most
alternates and cost 15 exact CDS on held-out core windows (`handover_review.md` §2). Run it only
to measure that, not as a candidate default.

**Identical inputs in every run:** the same genome FASTA, remapped Helixer GFF3, OrthoDB and
UniProt genBlastG alignments, Scallop and StringTie assemblies, the same per-read long-read file
(R4, R5; `gmb-run` collapses it with the same `fungi` consensus settings, recorded in each run's
`longread_consensus/longread_preparation.json`), the same `--gene-prefix`, and the same GMB
commit. Protein validation off (R1–R5).

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
LR=${EV}_minimap2_annotation.gtf          # 2.05 M per-read alignments; gmb-run collapses them
printf '%s\n' 1 2 3 4 5 6 7 8 9 10 11 12 13 > $OUT/inputs/core.txt
printf '%s\n' 14 15 16 17 18 19 20 21 > $OUT/inputs/accessory.txt
```

## 3. Commands

```bash
G=$OUT/inputs/genome.fa
EVID=(--genome $G --helixer $OUT/inputs/helixer.gff3
      --scallop ${EV}_scallop_annotation.gtf --stringtie ${EV}_stringtie_annotation.gtf
      --orthodb ${EV}_orthodb_annotation.gtf --uniprot ${EV}_genblast_annotation.gtf)

run() {   # run <name> <preset> [extra gmb-run args...]
  local n=$1 p=$2; shift 2
  /usr/bin/time -v gmb-run --preset $p "${EVID[@]}" "$@" \
    --gene-prefix ZTGMB --output-dir $OUT/$n 2> $OUT/$n.time
  echo "$n exit $?"
  git -C $GMB rev-parse HEAD > $OUT/$n/gmb_commit.txt
}
run R2_baseline        fungi
run R3_latest          fungi_latest_settings
run R4_latest_lr       fungi_latest_settings --long-read-alignments $LR
run R5_baseline_lr     fungi                 --long-read-alignments $LR
```

Each run's exit code must be 0 (`RUN_FUNGI.md` §6 lists the others). R4 and R5 each collapse the
same long-read file; the consensus GTFs (`R4_latest_lr/longread_consensus/minimap2_consensus.gtf`,
`R5_baseline_lr/…`) must be identical, which `longread_preparation.json` records as SHA-256.

Provenance captured automatically per run: `run_summary.json`, `build/resolved_config.yaml`
(+ SHA-256), `build/run_manifest.json` (versions, input hashes, runtime),
`longread_consensus/longread_preparation.json`, `finalise/handover_manifest.json` (output
inventory, configuration source). Wall time and peak RSS: `<run>.time` (GNU time `Elapsed` and
`Maximum resident set size`; on macOS use `/usr/bin/time -l`). Long-read preparation for the
whole genome took 34 s and ~0.5 GB in testing.
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
for r in R2_baseline R3_latest R4_latest_lr R5_baseline_lr; do
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

## 6. Decision criteria for confirming the latest fungal settings

Fixed in advance. Canonical output, **core chromosomes**, R3 compared with R2 (and R4 with R5 for
the long-read version). Promote only if **all** hold:

1. **Exact CDS:** R3 ≥ R2 + 1.0 percentage point of core reference genes, and R3 ≥ R1. (The
   17-window held-out estimate was +2.0 points over R2 and +1.2 over R1.)
2. **Regressions bounded:** reference genes losing an exact canonical CDS relative to R1 ≤ ⅓ of
   those gaining one (held-out: 5 lost vs 28 gained).
3. **Exact intron chains:** R3 ≥ R2.
4. **Structure:** canonical splits and merges ≤ R2. Full-gene-set splits + merges ≤ 1.2 × R1 (the
   latest settings' isoform rule should remove most of R2's isoform-driven merges).
5. **Loci:** missed reference genes ≤ R2; total gene count ≤ R2 × 1.02.
6. **Accessory chromosomes:** exact CDS not below R2 by more than 1% of accessory reference
   genes.
7. **QC and cost:** finalisation QC passes; wall time ≤ 1.25 × R2; peak RSS within the job limit.

If 1–7 hold, `fungi_latest_settings` is confirmed for production (and may later be folded into
`fungi` itself, with the golden fixture and its tests regenerated). If they fail, record which
criterion failed and decide whether to revert production runs to `fungi`.

**Long reads** are confirmed whenever R4 satisfies 1–7 against R5 *and* R4 ≥ R3 on criteria 1
and 3: per-read long-read alignments then become a standard input to `gmb-run`.

**Not a criterion, but must be reported:** the number of transcript-only antisense coding
predictions (R2–R5) and their support. A decision on filtering them (e.g. a coding-potential
gate) is separate and needs its own validation (`../known_issues.md`).

## 7. Expected outputs to keep

Per run: `finalise/` (handover), `run_summary.json`, `logs/`,
`build/{evidence_attribution.tsv,excluded_candidates.tsv,resolved_config.yaml,resolved_config_sha256,
run_manifest.json,summary.json}`, `preflight/`, `longread_consensus/longread_preparation.json`,
`<run>.time`, `gmb_commit.txt`; `cmp/<run>/<view>.<scope>/comparison_summary.{json,tsv}`,
`comparison_details.tsv`, `gene_splits.tsv`, `gene_merges.tsv`. Do not keep `build/` FASTA or
GFF3 copies: they duplicate `finalise/`.
