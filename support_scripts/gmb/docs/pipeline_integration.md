# Integrating GMB into an annotation pipeline

For a production developer whose pipeline already produces the evidence and needs GMB as the
**evidence-integration / final gene-model selection stage**.

---

## Where GMB sits

```
                         Genome FASTA
                              |
        +---------------------+---------------------+---------------------+
        |                     |                     |                     |
        v                     v                     v                     v
  backbone predictor    RNA-seq align +      long-read align +      protein-to-genome
  (Helixer/Tiberius)    assemble (Scallop,   collapse (optional)    alignment (OrthoDB,
                        StringTie, ...)                             UniProt, GenBlast)
        |                     |                     |                     |
        |  backbone_gtf       | short_read_gtf[]    | long_read_gtf[]     | protein_alignment_gtf[]
        +---------------------+----------+----------+---------------------+
                                         |
                                         v
                              +---------------------+
                              |   GMB PREFLIGHT     |   gate: stop on FAIL
                              +---------------------+
                                         |
                                         v
                              +---------------------+
                              |     GMB BUILD       |
                              +---------------------+
                                         |
                                         v
                              +---------------------+
                              |   GMB FINALISE      |
                              +---------------------+
                                         |
                                         v
                                     HANDOVER
                         finalise/ + attribution + manifests

                    ( gmb-compare — EVALUATION ONLY, off the production path )
```

---

## What your pipeline must supply

| pipeline artifact | GMB parameter | required | notes |
|---|---|---|---|
| `genome_fasta` | `--genome` | **yes** | defines the seqids every other file must use |
| `backbone_gtf` | `--helixer` *or* `--tiberius` | **yes** in practice | exactly one; mutually exclusive |
| `short_read_gtf[]` | `--scallop`, `--stringtie` | recommended | 0..N |
| `long_read_gtf[]` | `--minimap2` | optional | collapsed consensus, not raw alignments |
| `protein_alignment_gtf[]` | `--orthodb`, `--uniprot`, `--genblast` | recommended | support/veto only |
| `preset` | `--preset` | **yes** | `standard` \| `apicomplexa` \| `fungi` |
| `config` | `--config` | optional | repeatable overlays; last wins |
| `output_directory` | `--output-dir` | **yes** | one per stage |

**Flag names do not determine behaviour — evidence roles do.** An IsoQuant GTF passed to
`--minimap2` is handled correctly provided `scoring.longread_label` names its label. See
`input_contract.md`.

Your pipeline supplies **only these artifacts and a preset**. Nothing else about GMB needs to
be understood to run it safely.

---

## Calling it

### Shell

```bash
set -euo pipefail

gmb-preflight --preset "$PRESET" ${CONFIG:+--config "$CONFIG"} \
  --genome "$GENOME" --helixer "$BACKBONE" \
  --scallop "$SCALLOP" --stringtie "$STRINGTIE" --orthodb "$PROTEINS" \
  --output-dir "$OUT/preflight"          # exits 1 on FAIL

gmb-build --preset "$PRESET" ${CONFIG:+--config "$CONFIG"} \
  --genome "$GENOME" --helixer "$BACKBONE" \
  --scallop "$SCALLOP" --stringtie "$STRINGTIE" --orthodb "$PROTEINS" \
  --gene-prefix "$PREFIX" --output-dir "$OUT/build" --validate-fasta

# No --preset/--config: finalise runs under the build's own resolved_config.yaml.
gmb-finalise --build-dir "$OUT/build" --genome "$GENOME" \
  --output-dir "$OUT/finalise"
```

A complete runnable version is `examples/run_gene_model_builder.sh`.

### Python

For an orchestration layer that would rather not shell out:

```python
from gmb import run_gene_model_builder

result = run_gene_model_builder(
    genome="/data/genome.fa",
    backbone="/data/backbone.gff3",
    backbone_kind="helixer",          # or "tiberius"
    short_read=["/data/scallop.gtf", "/data/stringtie.gtf"],
    protein_alignment=["/data/orthodb.gtf"],
    long_read=None,                   # optional
    preset="fungi",
    config=["/data/my_overlay.yaml"],  # optional
    output_dir="/work/gmb",
    gene_prefix="XXGMB",
)

result.preflight_verdict   # "PASS" | "WARN" | "FAIL"
result.qc_passed           # bool
result.handover_dir        # the directory to hand over
result.annotation          # canonical annotated GFF3
result.proteins            # prot.fa
result.run_manifest        # run_manifest.json
```

`run_gene_model_builder` is the **only** supported Python entry point. An orchestration layer
should never import from `gmb.pipeline.scoring`, `gmb.pipeline.builder` or any other internal
module — those are free to change.

See `examples/run_gene_model_builder.py`.

---

## Exit codes and gating

| stage | exit 0 | non-zero |
|---|---|---|
| `gmb-preflight` | PASS or WARN | **1** = FAIL (override with `--allow-fail`); 2 = could not run |
| `gmb-build` | build completed and QC passed | **1** = QC failed — *outputs are still complete* |
| `gmb-finalise` | handover written | error |

> **A non-zero `gmb-build` exit under `--validate-fasta` is a QC verdict, not a crash.** The
> annotation is fully written. Read `build/fasta_qc_report.json` and
> `build/run_manifest.json` — the manifest is deliberately written *before* the QC exit, so a
> failing run still leaves a complete record of what produced it.

Recommended gating:

```
preflight FAIL   -> stop; the evidence bundle is unfit
build QC FAIL    -> stop; do not hand over (see qc.md)
finalise error   -> stop
```

---

## What to hand on

The `finalise/` directory plus `evidence_attribution.tsv`, `resolved_config.yaml` and
`run_manifest.json` from `build/`. Full list: `output_contract.md`.

**Never hand over `build/` sequence files** — `finalise/` is where they are regenerated from
the final GFF3.

---

## Resource planning

| | P. falciparum (23 Mb, 14 seqs) | Z. tritici (40 Mb, 21 seqs) |
|---|---|---|
| `gmb-preflight` | < 1 min | ~2 min |
| `gmb-build` wall | 2.1–3.5 h | 10.8–12.0 h |
| `gmb-build` peak RSS | 0.9–1.4 GB | 2.3–2.5 GB |
| `gmb-finalise` | seconds | seconds |
| `gmb-compare` peak RSS | ~0.6 GB | **4.2 GB** |

Measured with up to five builds sharing an 8-core/16 GB host, so wall times are pessimistic.
**GMB's hot path is single-threaded**: `--cpus-per-task 2` is enough. Scale memory from
**evidence volume**, not genome size — Z. tritici is heavier because its OrthoDB track is
732 MB / 1.69 M transcripts. See `examples/submit_gmb.slurm`.
