# GMB quickstart

> **Fungal genome? Follow [`../RUN_FUNGI.md`](../RUN_FUNGI.md) instead.** It is the supported
> fungal procedure. This page covers other clades and running the stages individually.

---

## 1. Install

```bash
mamba create -n gmb -c conda-forge python=3.12 'pandas>=2,<3' 'pyranges<=0.1.4' pyyaml numpy
mamba activate gmb
pip install /path/to/ensembl-anno/support_scripts/gmb
gmb-run --help
```

`pyranges` must be `<=0.1.4` — the API changed after that and it is the most common install
failure. GMB does **not** install Helixer, Tiberius, Scallop, StringTie, Minimap2, DIAMOND or
Psauron.

## 2. Pick a preset

| your situation | preset |
|---|---|
| fungi | `fungi_latest_settings` ([`../RUN_FUNGI.md`](../RUN_FUNGI.md)) |
| Apicomplexa with a general-purpose ab initio backbone | `apicomplexa` |
| new or unknown clade | `standard` |

Details and evidence: `presets.md`. A new clade: `creating_a_preset.md`.

## 3. Run

```bash
gmb-run --preset apicomplexa \
  --genome genome.fa --helixer backbone.gff3 \
  --scallop scallop.gtf --stringtie stringtie.gtf --orthodb orthodb.gtf \
  --gene-prefix XXGMB --output-dir out
```

Use `--tiberius` or `--backbone` instead of `--helixer` for another predictor. Omit any evidence
you do not have. Every evidence file must use the **same sequence names as the genome FASTA**.

**Long reads.** Per-read alignments (`--long-read-alignments`) are collapsed automatically, but
only for presets with a known consensus preset (the fungal ones). For other clades name one,
e.g. `--longread-consensus-preset pfalciparum_pure`, after checking its thresholds suit the
genome (`longread_consensus.md`), or supply already-collapsed models with
`--long-read-models`.

Exit codes (`0` ok, `2` arguments, `3` long-read preparation, `4` preflight, `5` build,
`6` sequence QC, `7` finalise), outputs and logs are as described in
[`../RUN_FUNGI.md`](../RUN_FUNGI.md) §4–6.

## 4. Running the stages individually (advanced)

`gmb-run` calls these; run them directly only for debugging or expert use. Give preflight and
build the **same** `--preset`, `--config` and evidence. `--minimap2` takes **collapsed** long-read
models only (`gmb-longread-consensus`, `longread_consensus.md`).

```bash
gmb-preflight --preset apicomplexa --genome genome.fa --helixer backbone.gff3 \
  --scallop scallop.gtf --orthodb orthodb.gtf --output-dir out/preflight      # exits 1 on FAIL
gmb-build --preset apicomplexa --genome genome.fa --helixer backbone.gff3 \
  --scallop scallop.gtf --orthodb orthodb.gtf \
  --gene-prefix XXGMB --output-dir out/build --validate-fasta
gmb-finalise --build-dir out/build --genome genome.fa --output-dir out/finalise
```

`gmb-finalise` needs no preset: it reads the build's `resolved_config.yaml`.

> **A non-zero `gmb-build` exit under `--validate-fasta` is a QC verdict, not a crash.** The
> annotation is fully written; read `build/fasta_qc_report.json`.

## 5. What you get

```
out/finalise/canonical/consensus.canonical_annotated.gff3   <- the annotation
out/finalise/consensus.gff3                                 all isoforms
out/finalise/{cdna,cds,prot}.fa                             sequences
out/finalise/{fasta_qc_report,utr_qc_report}.json           QC evidence
out/build/evidence_attribution.tsv                          why each model won
out/build/run_manifest.json                                 full provenance
out/run_summary.json                                        gmb-run status
```

**Hand over `finalise/`, not `build/`.** Full list: `output_contract.md`.

> **QC PASS does not mean high biological accuracy.** It means the annotation is internally
> consistent. See `qc.md`.

## 6. Evaluation against a reference (annotation lead, not part of a GMB run)

```bash
annotation-qc pairwise-compare --query out/finalise/consensus.gff3 \
  --reference reference.gff3 --genome genome.fa \
  --evaluation-mode protein_coding --reference-transcript-biotypes protein_coding \
  --outdir comparison
```

`annotation-qc` is in the `ensembl-genes` repository (its own Python ≥ 3.12 environment; GMB does
not depend on it). Metric guidance: `qc.md`.

---

## Next

| I want to… | read |
|---|---|
| know exactly what to supply | `input_contract.md` |
| know exactly what I get | `output_contract.md` |
| understand presets | `presets.md` |
| configure a new clade | `creating_a_preset.md` |
| wire GMB into a pipeline | `pipeline_integration.md` |
| understand QC | `qc.md` |
| fix something | `troubleshooting.md` |
| reproduce a run | `reproducibility.md` |
