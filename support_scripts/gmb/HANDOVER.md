# Gene Model Builder — production handover (fungal genomes)

GMB is ready to annotate fungal genomes at production scale. **To run it, follow
[RUN_FUNGI.md](RUN_FUNGI.md)**: install, supply the evidence, run one command (`gmb-run` with
`--preset fungi_latest_settings`), return the finalised files. This page sets out who is
responsible for what, what to return, and how to put GMB into an automated pipeline.

## Responsibilities

| | responsible for |
|---|---|
| **The person running GMB** | preparing suitable genome and evidence inputs (correct sequence names, stranded transcript evidence; RUN_FUNGI.md §2) · running `gmb-run --preset fungi_latest_settings` · checking that the run finished and passed structural QC (RUN_FUNGI.md §5) · returning the finalised annotation and the metadata listed below |
| **GMB** | collapsing per-read long-read alignments into consensus transcript models · validating every prepared input before building (preflight) · integrating the evidence and selecting gene models · finalisation: regenerating sequences from the final GFF3, structural QC, one canonical transcript per gene · recording configuration, input checksums and versions |
| **The annotation lead** | comparisons with reference annotations and with Helixer alone · biological assessment of annotation quality · deciding whether further tuning or development is needed |

The person running GMB does **not** compare annotations, choose between configurations,
collapse long reads by hand or benchmark anything.

## What to return

Per genome, the archive described in RUN_FUNGI.md §4:

- `finalise/` — the annotation (`canonical/consensus.canonical_annotated.gff3`, all isoforms with
  one canonical transcript per gene tagged `Ensembl_canonical`), `consensus.gff3`,
  `cdna.fa` / `cds.fa` / `prot.fa`, QC reports and `handover_manifest.json`;
- `run_summary.json` and `preflight/preflight_report.txt`;
- `build/evidence_attribution.tsv`, `build/resolved_config.yaml` (+ `_sha256`),
  `build/run_manifest.json`;
- `longread_consensus/longread_preparation.json` when long reads were used;
- the GMB commit installed and the evidence file names.

`build/` and `longread_consensus/` are intermediate. The GFF3 and FASTA files in `build/` are
pre-finalisation copies and must not be shared as results.

## Pipeline integration

- **One step.** Call `gmb-run` (or the Python function `gmb.run_gene_model_builder`, which
  `gmb-run` wraps; both apply the same policy). It is non-interactive, writes everything under
  `--output-dir`, logs to `logs/gmb-run.log` and records the outcome in `run_summary.json`
  (`ok`, `failed_stage`, `failure_reason`, output paths).
- **Exit codes:** 0 success · 2 invalid arguments · 3 long-read preparation · 4 preflight ·
  5 build · 6 sequence QC · 7 finalisation · 1 unexpected error.
- **Templates:** `examples/submit_gmb.slurm` (Slurm, self-contained),
  `examples/run_gene_model_builder.sh` (environment variables → `gmb-run`),
  `examples/run_gene_model_builder.py` (Python API).
- **Long reads.** Pass per-read alignments with `--long-read-alignments`; GMB collapses them
  with `gmb-longread-consensus --preset fungi` into `longread_consensus/` and reuses that
  consensus on a rerun only when the input checksum, settings and GMB version are unchanged.
  Pass already-collapsed models with `--long-read-models`. Per-read data given as models is
  refused.
- **Reruns** into the same directory replace `preflight/`, `build/` and `finalise/`.
- **Resources:** single-threaded; request 24 h and 16 GB for a first run of a ~40 Mb genome
  (RUN_FUNGI.md §3).
- **Installation:** a normal (non-editable) `pip install` contains everything, including
  `fungi_latest_settings` and the long-read consensus presets; run from any directory.

The individual stage commands (`gmb-preflight`, `gmb-build`, `gmb-finalise`,
`gmb-longread-consensus`) stay available for debugging and expert use
([README.md](README.md#commands)).

## Status and limitations

`fungi_latest_settings` was validated on 17 windows covering 38% of the *Z. tritici* genome,
where it improved exact CDS agreement over both Helixer alone and the historical `fungi` preset.
The first production-scale annotations are its genome-wide test, assessed by the annotation lead
([docs/release/production_scale_test_plan.md](docs/release/production_scale_test_plan.md)).
Known limitations relevant to running it at scale are listed in
[docs/known_issues.md](docs/known_issues.md); the scientific background is in
[docs/release/](docs/release/). Neither is needed to run GMB.

The historical fungal baseline (`--preset fungi`) is kept to reproduce earlier results;
`gmb-run` prints a note if it is used.
