# Gene Model Builder — handover

## What it does

GMB is the **evidence-integration and final gene-model-selection stage** of an annotation
pipeline. Give it a genome, an ab initio backbone, assembled transcript models and protein
alignments; it chooses the best-supported transcript structure at each locus and emits a
finished annotation with sequences, per-model evidence attribution, QC and full provenance. It
does not generate evidence and does not install the tools that do. A reference annotation is
**evaluation-only** and can never influence model selection.

## Start here

1. [`README.md`](README.md) — what it is, inputs, outputs, a copy-paste example
2. [`docs/quickstart.md`](docs/quickstart.md) — a first run, end to end
3. [`docs/release/release_readiness.md`](docs/release/release_readiness.md) — state of the
   2.0.0 release, what changed, and what to do before the first production-scale run
4. [`docs/known_issues.md`](docs/known_issues.md) — open issues by severity

## Normal workflow

```
gmb-preflight   →   gmb-build   →   gmb-finalise
 validate the        integrate        regenerate FASTA from the final
 evidence bundle     evidence and     GFF3, QC, pick canonical, write
 (exits 1 on FAIL)   select models    the handover
```

For *Z. tritici*, rename the Helixer GFF3 to the genome's sequence names first
(`tools/remap_helixer.py`; see README). Evaluation against a reference is **not** part of GMB: it lives in the `ensembl-genes`
repository as `annotation-qc pairwise-compare` (see `docs/qc.md`). The old `gmb-compare`
command only prints that pointer. Always pass the same `--preset` to `gmb-preflight` and
`gmb-build`.

## Presets

| preset | use it when |
|---|---|
| `standard` | new or unvalidated clade — neutral, no policy enabled, no clade assumption |
| `fungi` | fungi with a strong Helixer-like backbone — 2.0.0 adds a working 20 kb read-through filter, `transcript_linked` loci and detached-isoform removal (validated on 22% of Z. tritici; genome-wide confirmation is the production-scale test) |
| `apicomplexa` | Apicomplexa with a general-purpose backbone that under-calls introns |
| `configs/new_clade_template.yaml` | starting a clade of your own — copy and fill in `CHANGE_ME` |

Neither clade preset is universally optimal for every genome in its clade. What they encode is
an **evidence state**, not a taxonomy. See [`docs/configuration.md`](docs/configuration.md).

## Inputs

| input | flag | required |
|---|---|---|
| genome FASTA | `--genome` | **yes** |
| ab initio backbone (exactly one) | `--helixer`, `--tiberius` or generic `--backbone` | **yes** in practice |
| short-read transcript models | `--scallop`, `--stringtie` | recommended |
| protein-to-genome alignments | `--orthodb`, `--uniprot`, `--genblast` | recommended |
| long-read transcript models | `--minimap2` | optional |
| preset / overlays | `--preset`, `--config` | preset yes |
| reference annotation | *(no flag exists)* | **never — evaluation only** |

Flag names are slots, not behaviour: an assembler called anything works once listed under the
right evidence role. Full contract: [`docs/input_contract.md`](docs/input_contract.md).

## Outputs

**Hand over `finalise/`.** `build/` is intermediate.

| file | what it is |
|---|---|
| `finalise/canonical/consensus.canonical_annotated.gff3` | primary annotation, canonical tagged |
| `finalise/consensus.gff3` | all isoforms |
| `finalise/cdna.fa`, `cds.fa`, `prot.fa` | sequences, regenerated from the final GFF3 |
| `finalise/canonical/canonical_transcripts.tsv` | which transcript was chosen per gene |
| `finalise/fasta_qc_report.json`, `utr_qc_report.json` | QC evidence |
| `finalise/handover_manifest.json` | output inventory + which config finalisation used |
| `build/evidence_attribution.tsv` | why each model won |
| `build/resolved_config.yaml` | the exact configuration used |
| `build/run_manifest.json` | versions, input hashes, resolved policy, runtime |

Full list: [`docs/output_contract.md`](docs/output_contract.md).

## Current validation

All three whole-genome, each reproducing its expected annotation with **exact structural
equivalence** — no structure lost, gained or altered:

| genome | sequences | structures unchanged | hard QC | build time |
|---|---|---|---|---|
| **P. falciparum** | 14 | 4,402 / 4,402 | **PASS** | 1.6 h |
| **T. gondii** | 2,263 | 6,902 / 6,902 | **PASS** | 47.9 h |
| **Z. tritici** | 21 | 16,408 / 16,408 | **PASS** | 29 h |

Those builds predate the 2.0.0 release-preparation changes (UTR end support restricted to the
same sequence, performance fixes); see `docs/release/release_readiness.md`.

Tests: **862 passed, 0 failed, 3 skipped** (the 3 need a full `z_tritici/` data directory
that is not in the repository).

Release decision (2026-10-09): **GO** for a production-scale fungal test; **NO-GO** for a stable
tag until that test confirms the 2.0.0 fungal behaviour — see
[`docs/release/release_readiness.md`](docs/release/release_readiness.md).

## Known limitations

1. **Runtime.** T. gondii (2,263 sequences) took **47.9 h** and Z. tritici 29 h. Two costs in
   the isoform-selection loop were removed in 2.0.0 — a second, discarded selection pass and a
   per-model genome-wide UTR end-support scan (1.32 s per model on Z. tritici) — but the
   end-to-end gain has not been measured. Keep generous limits for the first production run
   and record the timing. **It does not hang**: those runs completed with exit code 0.
2. **Alternate isoforms** cause most remaining gene splits and merges, and all of GMB's
   measured CDS gain over Helixer on Z. tritici; the isoform policy for fungi is an open
   decision (`docs/known_issues.md`, P1).
3. **Clade breadth** — two Apicomplexa, one fungus. `fungi.yaml` was itself tuned against
   Z. tritici, so the preset is not independent of that organism; only the selection policy was
   tested independently.
4. **Z. tritici over-predicts** (16,137 genes against 10,931 in the reference). Inherited from
   the validated baseline. Not all of the excess is error: genome-wide, 81.8% of core Helixer
   genes absent from the reference have both protein and RNA-seq support; accessory ones are
   mostly weakly supported (`docs/release/z_tritici_validation.md`).

## If something goes wrong

- QC failed, or you are unsure what a QC pass means → [`docs/qc.md`](docs/qc.md)
- Anything else → [`docs/troubleshooting.md`](docs/troubleshooting.md)

Two things worth knowing in advance:

- A non-zero `gmb-build` exit under `--validate-fasta` is a **QC verdict, not a crash** — the
  annotation is fully written, and so is the run manifest.
- **A QC pass is not a quality verdict.** It proves internal consistency, not biological
  accuracy. Every Z. tritici arm passed all seven hard checks, including one that was
  biologically worse.

## Development evidence

For the validation behind the current production design, see
`diagnostics/gmb_productionisation/final_report.md`.

That directory sits **alongside this repository, not inside it** (it is working evidence, not
shipped code), so ask the annotation team for it if you do not have it. Everything else under
`diagnostics/` is an audit trail — you do not need to read any of it to use this module.
