# Gene Model Builder (GMB)

**Consensus gene annotation for eukaryotic genomes.** GMB is the
evidence-integration and final gene-model-selection stage of an annotation
pipeline: it takes an ab initio backbone, assembled transcript models and protein
alignments, and emits a finished annotation with per-model evidence attribution,
QC and provenance.

It does not generate evidence and it does not install the tools that do.

> **Annotating a fungal genome? Follow [RUN_FUNGI.md](RUN_FUNGI.md).** It is the complete,
> supported procedure: install, inputs, one command, outputs and troubleshooting.

---

## One command

```bash
gmb-run --preset fungi_latest_settings \
  --genome genome.fa --helixer helixer.gff3 \
  --orthodb orthodb.gtf --uniprot uniprot.gtf \
  --scallop scallop.gtf --stringtie stringtie.gtf \
  --long-read-alignments minimap2_reads.gtf \
  --gene-prefix XXGMB --output-dir out
```

`gmb-run` runs the whole workflow, non-interactively:

```
[long-read consensus]  ->  preflight  ->  build  ->  finalise
 only when per-read        validate the    integrate        regenerate sequences from
 long-read alignments      prepared        evidence,        the final GFF3, QC, canonical
 are supplied              inputs          select models    transcripts, handover
```

Hand over `out/finalise/`. Exit codes, logs and outputs: [RUN_FUNGI.md](RUN_FUNGI.md).
Per-read long-read alignments are collapsed into consensus transcript models automatically
(`out/longread_consensus/`); already-collapsed models go in `--long-read-models` instead.

## Upstream evidence preparation vs GMB

| stage | owner | output GMB consumes |
|---|---|---|
| genome preparation (soft-masking, sequence naming) | upstream pipeline | genome FASTA |
| short-read alignment + assembly (HISAT2/STAR → StringTie, Scallop) | upstream | transcript GTF, **stranded** |
| long-read alignment (Minimap2) with strand from the `ts` tag (`transcriptomic_annotation/transcript_strand.py`) | upstream | per-read GTF, **stranded** |
| protein-to-genome alignment (genBlastG / OrthoDB, UniProt) | upstream | alignment GTF |
| ab initio backbone prediction (Helixer, Tiberius, any predictor) | upstream | GFF3/GTF with CDS |
| long-read consensus | **GMB** (inside `gmb-run`) | `longread_consensus/` |
| input validation | **GMB** `gmb-preflight` | `preflight/` |
| evidence filtering, ORF/CDS/UTR inference, scoring, model selection | **GMB** `gmb-build` | `build/` |
| FASTA regeneration, sequence QC, canonical selection, handover manifest | **GMB** `gmb-finalise` | `finalise/` |
| comparison with a reference annotation, BUSCO/OMArk | evaluation, outside GMB | — |

GMB never re-derives strand from reads: a transcript's strand is whatever the evidence GTF
says, so strand must be correct upstream.

---

## What inputs do I need?

| input | required | flag |
|---|---|---|
| genome FASTA | **yes** | `--genome` |
| ab initio backbone (any predictor; exactly one) | **yes** in practice | `--helixer`, `--tiberius`, or generic `--backbone` (+ optional `--backbone-label`) |
| short-read transcript models | recommended | `--scallop`, `--stringtie` |
| protein-to-genome alignments | recommended | `--orthodb`, `--uniprot`, `--genblast` |
| long-read alignments, one record per read | **optional** | `gmb-run --long-read-alignments` (collapsed automatically) |
| long-read transcript models, already collapsed | **optional** | `gmb-run --long-read-models`, `gmb-build --minimap2` |
| DIAMOND protein DB | optional | config |

Every evidence file must use the **same sequence names as the genome FASTA**.

> **A reference annotation is never an input.** `gmb-preflight`, `gmb-build` and
> `gmb-finalise` expose no option that accepts one. References are for evaluation
> only, after the annotation exists, using `annotation-qc` from the `ensembl-genes`
> repository (see [docs/qc.md](docs/qc.md#evaluating-against-a-reference)).

Details: **[docs/input_contract.md](docs/input_contract.md)**

## What does it output?

```
out/finalise/canonical/consensus.canonical_annotated.gff3   <- the annotation
out/finalise/consensus.gff3                                 all isoforms
out/finalise/{cdna,cds,prot}.fa                             sequences
out/finalise/{fasta_qc_report,utr_qc_report}.json           QC evidence
out/build/evidence_attribution.tsv                          why each model won
out/build/resolved_config.yaml                              exact configuration
out/build/run_manifest.json                                 full provenance
out/longread_consensus/longread_preparation.json            long-read consensus provenance
out/run_summary.json                                        status of the whole run
```

**`finalise/` is the handover; `build/` is intermediate.** Sequences are
regenerated from the final GFF3 in `finalise/`, so annotation and sequence cannot
drift. Details: **[docs/output_contract.md](docs/output_contract.md)**

## Which preset do I use?

| preset | use it for | `backbone_intron_rescue` |
|---|---|---|
| **`fungi_latest_settings`** | **fungi — the recommended fungal settings** ([RUN_FUNGI.md](RUN_FUNGI.md)) | `off` |
| `fungi` | historical fungal baseline, kept to reproduce earlier results | `off` |
| `apicomplexa` | Apicomplexa with a general-purpose ab initio backbone | `auto` |
| `standard` | a new or unknown clade | `off` |

`fungi_latest_settings` is `fungi` plus four selection settings
(`gmb/configs/fungi_latest_settings.yaml`; evidence in
[docs/release/evidence_integration.md](docs/release/evidence_integration.md)). The presets differ
because the **evidence states** differ, not because the clades do. The evidence for each choice
is in **[docs/presets.md](docs/presets.md)**.

For a new clade, start from `configs/new_clade_template.yaml` and follow
**[docs/creating_a_preset.md](docs/creating_a_preset.md)**.

## What is optional?

- **Long-read evidence.** With none supplied there is no long-read preparation, no special
  case and no penalty.
- **Protein validation** (DIAMOND + Psauron). It only *scores* models that already
  exist; it can change which isoform wins, never which structures are possible.
- **Evaluation against a reference.** Not part of GMB; see `ensembl-genes` `annotation-qc`.
- **Every biological policy.** All default off.

## What does QC guarantee?

Every handover annotation satisfies, verifiably:

```
0 cDNA / CDS / protein sequence mismatches
0 unintended internal stop codons
0 gene-boundary violations      (gene span == union of its transcripts)
0 UTR invariant violations
0 cross-seqid chimeras
1 canonical transcript per gene
```

> **QC PASS does not mean high biological accuracy.** It means the annotation is
> internally consistent. In one validation, all four gene sets passed every hard
> check — and one of them was *worse than doing nothing*.

Details: **[docs/qc.md](docs/qc.md)**

## Where do I look when something fails?

| question | file |
|---|---|
| did the run finish, and where did it stop? | `run_summary.json`, `logs/gmb-run.log` |
| were the long reads collapsed, reused or skipped? | `logs/gmb-run.log`, `longread_consensus/longread_preparation.json` |
| were my inputs fit to build from? | `preflight/preflight_report.txt` |
| what did the build do? | `logs/build.log` |
| did it pass QC? | `finalise/fasta_qc_report.json` |
| why was this model chosen? | `build/evidence_attribution.tsv` |
| what settings applied? | `build/resolved_config.yaml` |
| what produced this run? | `build/run_manifest.json` |

Symptom-first guide: **[docs/troubleshooting.md](docs/troubleshooting.md)**

---

## Install

```bash
mamba create -n gmb -c conda-forge python=3.12 'pandas>=2,<3' 'pyranges<=0.1.4' \
    pyyaml numpy
mamba activate gmb
pip install /path/to/ensembl-anno/support_scripts/gmb     # -e for development
gmb-run --help
```

Python ≥ 3.10 (CI: 3.10, 3.11, 3.13). Runtime dependencies are only pandas, pyranges,
PyYAML and numpy. `pyranges` must be `<=0.1.4` — the API changed after that and it is
the most common install failure. Development tools: `pip install -e ".[dev]"`.

**GMB does not install** Helixer, Tiberius, Scallop, StringTie, Minimap2, DIAMOND
or Psauron. DIAMOND and Psauron are *invoked* by GMB if you give it their paths;
they generally need their own environments.

## Commands

| command | purpose | |
|---|---|---|
| **`gmb-run`** | **the whole workflow: long-read consensus, preflight, build, finalise** | **recommended** |
| `gmb-preflight` | validate an evidence bundle before building | advanced / individual stages |
| `gmb-build` | integrate evidence and select gene models | advanced |
| `gmb-finalise` | regenerate FASTA, QC, canonical, handover | advanced |
| `gmb-longread-consensus` | collapse per-read long-read alignments into models (run by `gmb-run`) | expert / debugging ([longread_consensus.md](docs/longread_consensus.md)) |
| `gmb-canonical-selection`, `gmb-interpro-*` | specialised tools (also run by `gmb-finalise`) | no |
| `gmb-compare` | **removed** — prints where reference comparison moved (`ensembl-genes` `annotation-qc`) | no |

Running the stages individually: give `gmb-preflight` and `gmb-build` the same `--preset`,
`--config` and evidence, and `--minimap2` only a collapsed long-read track; `gmb-finalise` needs
only `--build-dir`, `--genome` and `--output-dir`.

## Python API

```python
from gmb import run_gene_model_builder

result = run_gene_model_builder(
    genome="/data/genome.fa",
    backbone="/data/helixer.gff3", backbone_kind="helixer",
    short_read=["/data/scallop.gtf", "/data/stringtie.gtf"],
    protein_alignment=["/data/orthodb.gtf", "/data/genblast.gtf"],
    long_read_alignments="/data/minimap2_reads.gtf",     # optional, per-read
    preset="fungi_latest_settings", output_dir="/work/gmb", gene_prefix="XXGMB",
)
if result.ok:
    ship(result.handover_dir)
else:
    print(result.failed_stage, result.failure_reason)
```

`run_gene_model_builder` is the **only** supported Python entry point; `gmb-run` is a thin
wrapper around it, so both apply the same long-read preparation policy. Everything
under `gmb.pipeline`, `gmb.preflight` and `gmb.provenance` is internal.

## Two design rules

**Evidence roles, not tool names.** Selection logic never tests a literal tool
name. Each source resolves to a role (`backbone`, `short_read_transcriptomic`,
`long_read_transcriptomic`, `protein_alignment`) and everything — ranking,
retention, corroboration, **and numeric weights** — acts on the role. An assembler
called `IsoQuant` or `AssemblerX` behaves exactly like `StringTie` once listed
under `shortread_labels`.

**Correctness is not policy.** UTR correction, FASTA regeneration, duplicate
collapse, gene-boundary recomputation, cross-seqid ID namespacing and sequence QC
are **always on and not configurable**. Biological tuning —
structural corroboration, protein-support mode, long-read handling, backbone
intron rescue, weights — is **explicit and off by default**.

## Documentation

### START HERE

| | |
|---|---|
| **[RUN_FUNGI.md](RUN_FUNGI.md)** | **annotate a fungal genome: the supported procedure** |
| [HANDOVER.md](HANDOVER.md) | production handover: responsibilities, what to return, pipeline integration |
| [docs/pipeline_integration.md](docs/pipeline_integration.md) | wiring GMB into a workflow manager |
| [docs/qc.md](docs/qc.md) | what QC guarantees, and what it does not |

### REFERENCE

Consult when you need the detail.

| | |
|---|---|
| [quickstart.md](docs/quickstart.md) | other clades; running the stages individually |
| [configuration.md](docs/configuration.md) | presets, overrides, adding an evidence source, new clades |
| [input_contract.md](docs/input_contract.md) | formats, coordinates, IDs, failure behaviour |
| [output_contract.md](docs/output_contract.md) | every output file and which ones to hand over |
| [presets.md](docs/presets.md) | the preset catalogue and the evidence behind each |
| [creating_a_preset.md](docs/creating_a_preset.md) | the 8-step procedure for a validated clade preset |
| [architecture.md](docs/architecture.md) | module map and production boundary |
| [reproducibility.md](docs/reproducibility.md) | run manifests and reproducing a run |
| [troubleshooting.md](docs/troubleshooting.md) | symptom-first fixes |
| [known_issues.md](docs/known_issues.md) | open issues by severity, and what was fixed |
| [release/release_readiness.md](docs/release/release_readiness.md) | 2.0.0 release state, changes, outstanding work |
| [release/z_tritici_validation.md](docs/release/z_tritici_validation.md) | targeted fungal validation |

Specialised: [longread_consensus.md](docs/longread_consensus.md),
[canonical_selection.md](docs/canonical_selection.md),
[interpro_resolver.md](docs/interpro_resolver.md).

Examples: `examples/submit_gmb.slurm` (Slurm), `examples/run_gene_model_builder.sh`
(environment-variable wrapper), `examples/run_gene_model_builder.py` (Python API).

### Development and scientific evidence — not needed to run GMB

[docs/release/](docs/release/): release readiness, the *Z. tritici* validation, the
evidence-integration investigation behind `fungi_latest_settings`, the handover review and the
production-scale comparison plan. `validation/` reproduces the evaluation and holds a fixture
smoke test.

## Tests

```bash
pytest tests/ -q                       # ~20 s; includes the 500 kb Z. tritici fixture builds
pytest tests/ -q -m "not integration"  # unit tests only
```

Three slow tests (`RUN_SLOW_INTEGRATION=1`) need a full `z_tritici/` data directory that is
not shipped with the repository; they skip otherwise.

`tests/test_production_contract.py` pins the production contract: role-based
weight resolution using **invented** tool names, the rescue applicability gate,
long-read optionality, preflight splice classification, and the resolved values of
every shipped preset.
