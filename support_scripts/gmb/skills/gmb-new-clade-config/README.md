# gmb-new-clade-config

An agent skill that builds a **Gene Model Builder (GMB) configuration for a taxonomic group
that has no validated preset** — plants, oomycetes, nematodes, protists, anything — from that
group's own annotation evidence.

It measures your evidence files, decides where the neutral `standard` preset is inappropriate,
writes a small config that differs from `standard` only where a measurement says it should,
runs GMB with it, and reports every change with its reason. It never copies the `fungi` or
`apicomplexa` presets, and it never lets a reference annotation influence the configuration.

Agent behaviour is defined in [`SKILL.md`](SKILL.md). For using GMB itself, see
[`../../README.md`](../../README.md) and [`../../docs/quickstart.md`](../../docs/quickstart.md).

## What you give it

| Input | Required? |
|---|---|
| Clade name | yes |
| Genome FASTA | yes |
| Ab initio / backbone annotation (GFF3 or GTF) | yes |
| Short-read transcript annotation(s) | strongly recommended (up to 2) |
| Long-read transcript annotation | optional (1) |
| Protein-to-genome alignments | optional, recommended (up to 3) |
| Assembly report / seqname map | optional — only if sequence names differ |
| Trusted reference annotation | optional, **evaluation only — GFF3 only** |
| A second representative genome | required for `CLADE_VALIDATED` |

Tools can be called anything — you state each file's **role**. See
[`example_manifest.yaml`](example_manifest.yaml). Evidence files must be uncompressed, with
`.gtf` or `.gff3` extensions matching their content. The reference must be **GFF3**; a GTF
reference is refused rather than risk wrong metrics.

An optional machine overlay may hold **only** tool and database paths and hardware settings
(DIAMOND/Psauron locations, worker counts). Anything that could change which gene model wins is
refused there and belongs in the clade config, with a reason.

## What it produces

- `configs/<clade>.yaml` — the new config, as a delta from `standard`
- `gmb_clade_config_report.md` — evidence state, every setting changed and why, what was
  deliberately left alone, QC, evaluation, and the validation status
- `config_decisions.tsv` — one row per setting considered
- `resolved_config.yaml`, preflight reports, a slot map, and the finished GMB runs
  (`standard` and candidate) with hard-QC summaries

## What "validated" means

| status | meaning |
|---|---|
| `PROVISIONAL` | derived from the evidence and passes GMB's hard QC; not checked against a reference |
| `SINGLE_GENOME_VALIDATED` | also performs acceptably against a trusted reference on **one** genome |
| `CLADE_VALIDATED` | frozen after genome A, then run **unchanged** on an independent genome B, and acceptable there too |

Tuning and testing on the same genome is never "clade validation". With one genome, the skill
stops at `SINGLE_GENOME_VALIDATED` and tells you which genome to test next.

Before any reference is used, the skill freezes the config you wrote **and** the resolved config
the build actually ran with. Evaluation is refused if either has changed since.

## Example request

```
Create a GMB configuration for nematodes.

Genome:
    genome.fa

Backbone:
    predictor.gff3

Short-read transcript evidence:
    assembler_a.gtf
    assembler_b.gtf

Protein evidence:
    orthodb.gtf
    uniprot.gtf

Long-read:
    long_reads.gtf

Reference annotation for evaluation only:
    reference.gff3

Start from the standard GMB preset, assess the evidence using gmb-preflight, create a
provisional configuration, run the validation workflow, and explain every setting that
differs from standard.
```

## Installing the skill

It lives with GMB so it versions with the code it describes. To make it available to Claude
Code, link it into a skills directory, e.g.:

```bash
mkdir -p .claude/skills
ln -s "$PWD/support_scripts/gmb/skills/gmb-new-clade-config" .claude/skills/gmb-new-clade-config
```

The helper script needs the same Python environment as GMB (it uses PyYAML, and reads GMB's
rescue-gate thresholds from the installed package). Its tests are in `tests/`:

```bash
python -m pytest support_scripts/gmb/skills/gmb-new-clade-config/tests
```
