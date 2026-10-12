# Z. tritici evidence-integration evaluation (17 windows)

Reproduces every number in [`docs/release/evidence_integration.md`](../../docs/release/evidence_integration.md)
and [`docs/release/handover_review.md`](../../docs/release/handover_review.md) from the original
input data. Not part of the installed `gmb` package. Large intermediate files are written only
to `$WORK`; [`results/`](results/) keeps the small summary tables from the reference run.

## Requirements

| item | version used for `results/` |
|---|---|
| GMB | this repository, `gmb/evidence-tuning` (commit in `results/provenance.txt`) |
| Python env for GMB | Python 3.11, pandas 2.x, pyranges 0.1.4, numpy, pyyaml (`results/provenance.txt`) |
| `annotation-qc` | ensembl-genes, separate Python ≥ 3.12 environment (`results/provenance.txt`) |
| DIAMOND, Psauron, a UniProt DIAMOND DB | only for the optional `U_all_pv` arm |
| original inputs (read only) | `Fungi_release_tests/input_data/GCA_000219625.1_*` (Helixer GFF3, assembly report, Scallop, StringTie, OrthoDB and UniProt genBlastG, per-read Minimap2) and `reference_data/Zymoseptoria_tritici_GCA_000219625.1.{fa,gff3}.bgz` |
| disk / time | 2.3 GB in `$WORK`; on an 8-core laptop ~13 min for input preparation and the 272 jobs without protein validation (`JOBS=6`), plus ~6 min for `U_all_pv` with `JOBS=2` (Psauron is multi-threaded, so that arm runs only when named: `JOBS=2 bash run_matrix.sh U_all_pv`) |

## Windows

[`windows.tsv`](windows.tsv) was fixed before any result was seen. *Development* windows (tuning
allowed) are the nine used in release preparation: core 1:2–3 Mb, 3:1–2 Mb, 5:0.5–1.5 Mb, 7:1–2 Mb,
9:1–2 Mb, 12:0.3–1.3 Mb and accessory chromosomes 14, 18, 21 whole. *Held-out* windows
(evaluated only after every setting was fixed) are on chromosomes not used for tuning: core
2:1–2 Mb, 4:1–2 Mb, 6:1–2 Mb, 8:1–2 Mb, 10:0.3–1.3 Mb, 13:0.1–1.1 Mb and accessory 16, 19
whole. Windows were chosen by hand, not sampled at random: one 1 Mb block per core
chromosome, mostly away from the chromosome start, and whole accessory chromosomes (they are
small). The development set came from release preparation; the held-out set was taken from the
remaining core chromosomes and two unused accessory chromosomes before any held-out result was
computed. Together: 14.9 Mb, 38% of the genome; 3,739 core and 411 accessory reference
protein-coding genes.

A model is included in a window only when its whole transcript/gene lies inside it
(`subset_windows.py`); the genome passed to GMB is the whole chromosome, so coordinates are
unchanged. Reference genes crossing a window boundary therefore count as missed or partial
equally in every arm.

## Run

```bash
export WORK=/scratch/ztritici_eval                     # any empty directory
export INPUT_DATA=/path/Fungi_release_tests/input_data
export REFERENCE_DATA=/path/Fungi_release_tests/reference_data
export GMB_BIN=$(dirname $(which gmb-build))           # env with GMB installed
export AQC=/path/to/ensembl-genes-env/bin/annotation-qc
export ENSEMBL_GENES_DIR=/path/to/ensembl-genes         # optional, recorded in provenance
# export PV_CONFIG=/path/pv_local.yaml                  # optional: enables U_all_pv
bash prepare_inputs.sh      # decompress, remap Helixer, subset windows, collapse long reads
bash run_matrix.sh          # every arm in arms.tsv x every window (or name arms as arguments)
bash make_results.sh        # all summary tables -> results/
```

`results/` was produced this way from commit be54b99 plus the handover-review changes
(`results/provenance.txt`); regenerate it after any code change.

`run_arm.sh` builds, finalises and evaluates the **finalised** outputs three ways with
`annotation-qc pairwise-compare --evaluation-mode protein_coding
--reference-transcript-biotypes protein_coding --region <window>`: `cmp_consensus` (all
isoforms, `fin/consensus.gff3`), `cmp_primary` (the build primary `.t1` of each gene) and
`cmp_canonical` (`fin/canonical/consensus.canonical_annotated.gff3`,
`--query-transcript-selection canonical`, the handover view).

Arms ([`arms.tsv`](arms.tsv)): `A_helixer` (raw Helixer); `B_…`–`E_…` the baseline fungi preset
with the tracks named; `U_…` the same with the experimental overlay; `*_lr` add the collapsed
long-read track; `U_all_pv` adds protein validation.

## Tables

```bash
P=$GMB_BIN/python
$P -I summarise.py $WORK/runs windows.tsv heldout A_helixer E_all U_all E_lr U_all_lr --out $WORK/final_heldout.tsv
$P introns.py $WORK/runs h2,h4,h6,h8,h10,h13 A_helixer:consensus E_all:canonical U_all:canonical U_all_lr:canonical
$P refgene_diff.py $WORK/runs A_helixer U_all canonical h2,h4,h6,h8,h10,h13      # core gains / losses
$P model_changes.py $WORK/runs $WORK/win $WORK/w/ref.gff3 h2,h4,h6,h8,h10,h13 E_all U_all U_all_lr
$P geneset_quality.py $WORK/runs h2,h4,h6,h8,h10,h13,h16,h19 A_helixer E_all U_all E_lr U_all_lr
$P canonical_overrides.py $WORK windows.tsv E_all U_all E_lr U_all_lr U_all_pv
$P promotion_audit.py $WORK U_all windows.tsv
```

`summarise.py` pools windows of one set and class. R (reference genes) and Q (predicted genes)
are disjoint across windows, so counts add up; percentages use pooled R; F1 is pooled
`2 TP / (R + Q)` with TP = one-to-one exact-CDS pairs.

`diagnosis/` holds the exploratory scripts behind the diagnosis sections of the report
(per-structure audit of `select_isoforms`, competing-structure features, missed genes,
new-locus protein coverage, CDS-disjoint isoforms, single-locus inspection). Each documents its
own arguments.
