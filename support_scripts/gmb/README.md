# Gene Model Builder (GMB)

**Consensus gene annotation for eukaryotic genomes.** GMB is the
evidence-integration and final gene-model-selection stage of an annotation
pipeline: it takes an ab initio backbone, assembled transcript models and protein
alignments, and emits a finished annotation with per-model evidence attribution,
QC and provenance.

It does not generate evidence and it does not install the tools that do.

> **Release status:** see [docs/release/release_readiness.md](docs/release/release_readiness.md)
> and [docs/known_issues.md](docs/known_issues.md) before a production run.

---

## The 60-second version

```bash
# 1. validate the evidence bundle (exits 1 on FAIL)
gmb-preflight --preset fungi --genome genome.fa --helixer backbone.gff3 \
  --scallop scallop.gtf --stringtie stringtie.gtf --orthodb orthodb.gtf \
  --output-dir out/preflight

# 2. integrate evidence and select gene models
gmb-build --preset fungi --genome genome.fa --helixer backbone.gff3 \
  --scallop scallop.gtf --stringtie stringtie.gtf --orthodb orthodb.gtf \
  --gene-prefix XXGMB --output-dir out/build --validate-fasta

# 3. regenerate FASTA from the final GFF3, QC, canonical, handover
gmb-finalise --build-dir out/build --genome genome.fa \
  --output-dir out/finalise
```

Hand over `out/finalise/`. Full walk-through: **[docs/quickstart.md](docs/quickstart.md)**.

Always pass the **same `--preset` to preflight and build**. If omitted, both fall back to
`fungi` for backward compatibility; the Python API defaults to `standard`.

---

## Upstream evidence preparation vs GMB

| stage | owner | output GMB consumes |
|---|---|---|
| genome preparation (soft-masking, sequence naming) | upstream pipeline | genome FASTA |
| short-read alignment + assembly (HISAT2/STAR → StringTie, Scallop) | upstream | transcript GTF, **stranded** |
| long-read alignment (Minimap2) → strand from the `ts` tag (`transcriptomic_annotation/transcript_strand.py`) → **collapse to transcript models** (`gmb-longread-consensus` or an isoform collapser) | upstream | transcript GTF |
| protein-to-genome alignment (genBlastG / OrthoDB, UniProt) | upstream | alignment GTF |
| ab initio backbone prediction (Helixer, Tiberius, any predictor) | upstream | GFF3/GTF with CDS |
| input validation | **GMB** `gmb-preflight` | `preflight_report.json` |
| evidence filtering, ORF/CDS/UTR inference, scoring, model selection | **GMB** `gmb-build` | `build/` |
| FASTA regeneration, sequence QC, canonical selection, handover manifest | **GMB** `gmb-finalise` | `finalise/` |
| comparison with a reference annotation, BUSCO/OMArk | evaluation, outside GMB | — |

GMB never re-derives strand from reads: a transcript's strand is whatever the evidence GTF
says, so strand must be correct upstream. Raw per-read long-read alignments are **not** valid
`--minimap2` input; collapse them first.

## Fungal example (Zymoseptoria tritici, GCA_000219625.1)

```bash
G=genome.fa                        # sequence names 1..21, identical in every evidence file
EV=input_data/GCA_000219625.1      # Ensembl anno evidence outputs
# Helixer names sequences by GenBank accession (CM001196.1 ...); rename to the genome's
# names first. Preflight fails the track otherwise.
python tools/remap_helixer.py --input ${EV}_helixer.gff3 \
  --assembly-report ${EV}_MYCGR_v2.0_assembly_report.txt --output helixer.gff3
gmb-preflight --preset fungi --genome $G --helixer helixer.gff3 \
  --scallop ${EV}_scallop_annotation.gtf --stringtie ${EV}_stringtie_annotation.gtf \
  --orthodb ${EV}_orthodb_annotation.gtf --uniprot ${EV}_genblast_annotation.gtf \
  --output-dir out/preflight
gmb-build --preset fungi --genome $G --helixer helixer.gff3 \
  --scallop ${EV}_scallop_annotation.gtf --stringtie ${EV}_stringtie_annotation.gtf \
  --orthodb ${EV}_orthodb_annotation.gtf --uniprot ${EV}_genblast_annotation.gtf \
  --gene-prefix ZTGMB --output-dir out/build --validate-fasta
gmb-finalise --build-dir out/build --genome $G --output-dir out/finalise
```

The `genblast_annotation.gtf` from Ensembl anno is the genBlastG alignment of **UniProt**
proteins; pass it as `--uniprot` to match the validated baseline (labels only affect
attribution). Do not pass the raw per-read `minimap2_annotation.gtf` (2 M reads for this
genome; preflight warns) — see the upstream table above. Rename the backbone upstream as
shown rather than with `gmb-build --assembly-report`: preflight and the Python API have no
seqname-mapping option, so only a renamed file can be validated and run the same way
everywhere. The `fungi` preset applies a 20 kb transcript-span filter and
`locus_clustering: transcript_linked` (2.0.0). The evidence behind these choices:
[docs/release/z_tritici_validation.md](docs/release/z_tritici_validation.md).

---

## What inputs do I need?

| input | required | flag |
|---|---|---|
| genome FASTA | **yes** | `--genome` |
| ab initio backbone (any predictor; exactly one) | **yes** in practice | `--helixer`, `--tiberius`, or generic `--backbone` (+ optional `--backbone-label`) |
| short-read transcript models | recommended | `--scallop`, `--stringtie` |
| protein-to-genome alignments | recommended | `--orthodb`, `--uniprot`, `--genblast` |
| long-read transcript models | **optional** | `--minimap2` |
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
```

**`finalise/` is the handover; `build/` is intermediate.** Sequences are
regenerated from the final GFF3 in `finalise/`, so annotation and sequence cannot
drift. Details: **[docs/output_contract.md](docs/output_contract.md)**

## Which preset do I use?

| preset | use it for | `backbone_intron_rescue` |
|---|---|---|
| `standard` | a new or unknown clade | `off` |
| `apicomplexa` | Apicomplexa with a general-purpose ab initio backbone | `auto` |
| `fungi` | fungi with a strong Helixer-like backbone | `off` |

The presets differ because the **evidence states** differ, not because the clades
do — the same policy that improved P. falciparum made Z. tritici measurably worse.
The evidence for each choice is in **[docs/presets.md](docs/presets.md)**.

For a new clade, start from `configs/new_clade_template.yaml` and follow
**[docs/creating_a_preset.md](docs/creating_a_preset.md)**.

## What is optional?

- **Long-read evidence.** With none supplied there is no special case, no penalty,
  and the long-read guard cannot fire.
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
| were my inputs fit to build from? | `preflight/preflight_report.txt` |
| what did the build do? | `build/gmb.log` |
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
pip install -e /path/to/ensembl-anno/support_scripts/gmb
gmb-build --help
```

Python ≥ 3.10 (CI: 3.10, 3.11, 3.13). Runtime dependencies are only pandas, pyranges,
PyYAML and numpy. `pyranges` must be `<=0.1.4` — the API changed after that and it is
the most common install failure. Development tools: `pip install -e ".[dev]"`.

**GMB does not install** Helixer, Tiberius, Scallop, StringTie, Minimap2, DIAMOND
or Psauron. DIAMOND and Psauron are *invoked* by GMB if you give it their paths;
they generally need their own environments.

## Commands

| command | purpose | production path |
|---|---|---|
| `gmb-preflight` | validate an evidence bundle before building | **yes** |
| `gmb-build` | integrate evidence and select gene models | **yes** |
| `gmb-finalise` | regenerate FASTA, QC, canonical, handover | **yes** |
| `gmb-longread-consensus` | collapse long-read alignments into models | optional, upstream |
| `gmb-canonical-selection`, `gmb-interpro-*` | specialised tools (also run by `gmb-finalise`) | no |
| `gmb-compare` | **removed** — prints where reference comparison moved (`ensembl-genes` `annotation-qc`) | no |

## Python API

```python
from gmb import run_gene_model_builder

result = run_gene_model_builder(
    genome="/data/genome.fa",
    backbone="/data/backbone.gff3", backbone_kind="helixer",
    short_read=["/data/scallop.gtf", "/data/stringtie.gtf"],
    protein_alignment=["/data/orthodb.gtf"],
    preset="fungi", output_dir="/work/gmb", gene_prefix="XXGMB",
)
if result.ok:
    ship(result.handover_dir)
```

`run_gene_model_builder` is the **only** supported Python entry point. Everything
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

Four documents cover the normal job of running GMB. You should not need anything else.

| | |
|---|---|
| 1. [quickstart.md](docs/quickstart.md) | run it end to end |
| 2. [pipeline_integration.md](docs/pipeline_integration.md) | wire it into an existing pipeline — what to hand it, what comes back |
| 3. [configuration.md](docs/configuration.md) | presets, overrides, adding an evidence source, new clades |
| 4. [qc.md](docs/qc.md) | what QC guarantees, and what it does not |

### REFERENCE

Consult when you need the detail.

| | |
|---|---|
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

Examples: `examples/run_gene_model_builder.sh`, `.py`, `submit_gmb.slurm`.

### Development evidence — not needed to run GMB

For the validation behind the current production design, see
`diagnostics/gmb_productionisation/final_report.md` in the annotation team's working
evidence. That directory sits **alongside this repository, not inside it**; ask the team for
it if you need it. You do not need it to use this module.

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
