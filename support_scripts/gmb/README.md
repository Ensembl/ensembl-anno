# Gene Model Builder (GMB)

**Consensus gene annotation for eukaryotic genomes.** GMB is the
evidence-integration and final gene-model-selection stage of an annotation
pipeline: it takes an ab initio backbone, assembled transcript models and protein
alignments, and emits a finished annotation with per-model evidence attribution,
QC and provenance.

It does not generate evidence and it does not install the tools that do.

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

---

## What inputs do I need?

| input | required | flag |
|---|---|---|
| genome FASTA | **yes** | `--genome` |
| ab initio backbone (Helixer or Tiberius) | **yes** in practice | `--helixer` / `--tiberius` |
| short-read transcript models | recommended | `--scallop`, `--stringtie` |
| protein-to-genome alignments | recommended | `--orthodb`, `--uniprot`, `--genblast` |
| long-read transcript models | **optional** | `--minimap2` |
| DIAMOND protein DB | optional | config |

Every evidence file must use the **same sequence names as the genome FASTA**.

> **A reference annotation is never an input.** `gmb-preflight`, `gmb-build` and
> `gmb-finalise` expose no option that accepts one. References are for evaluation
> only, through `gmb-compare`, after the annotation exists.

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
- **Evaluation** (`gmb-compare`). Not part of the production path.
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
mamba create -n gmb -c conda-forge python=3.12 pandas 'pyranges<=0.1.4' \
    biopython pyyaml matplotlib numpy
mamba activate gmb
pip install -e /path/to/ensembl-anno/support_scripts/gmb
gmb-build --help
```

Python ≥ 3.10. `pyranges` must be `<=0.1.4` — the API changed after that and it is
the most common install failure.

**GMB does not install** Helixer, Tiberius, Scallop, StringTie, Minimap2, DIAMOND
or Psauron. DIAMOND and Psauron are *invoked* by GMB if you give it their paths;
they generally need their own environments.

## Commands

| command | purpose | production path |
|---|---|---|
| `gmb-preflight` | validate an evidence bundle before building | **yes** |
| `gmb-build` | integrate evidence and select gene models | **yes** |
| `gmb-finalise` | regenerate FASTA, QC, canonical, handover | **yes** |
| `gmb-compare` | evaluate against a reference | evaluation only |
| `gmb-longread-consensus` | collapse long-read alignments into models | optional, upstream |
| `gmb-canonical-selection`, `gmb-interpro-*`, `gmb-visualize` | specialised tools | no |

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

Specialised: [longread_consensus.md](docs/longread_consensus.md),
[canonical_selection.md](docs/canonical_selection.md),
[interpro_resolver.md](docs/interpro_resolver.md).

Examples: `examples/run_gene_model_builder.sh`, `.py`, `submit_gmb.slurm`.

### Development evidence — not needed to run GMB

For the validation behind the current production design, see
`diagnostics/gmb_productionisation/final_report.md`. Everything else under `diagnostics/` is
an audit trail: you do not need to read any of it to use this module.

## Tests

```bash
pytest tests/ -q
```

`tests/test_production_contract.py` pins the production contract: role-based
weight resolution using **invented** tool names, the rescue applicability gate,
long-read optionality, preflight splice classification, and the resolved values of
every shipped preset.
