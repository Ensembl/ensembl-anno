# Integrating GMB into an annotation pipeline

For a production developer whose pipeline already produces the evidence and needs GMB as the
**evidence-integration / final gene-model selection stage**. For fungi, the procedure and the
recommended settings are in [`../RUN_FUNGI.md`](../RUN_FUNGI.md); this page covers the
interface.

---

## Where GMB sits

```
                         Genome FASTA
                              |
        +---------------------+---------------------+---------------------+
        v                     v                     v                     v
  backbone predictor    RNA-seq align +      long-read align        protein-to-genome
  (Helixer/Tiberius)    assemble (Scallop,   (Minimap2, strand      alignment (OrthoDB,
                        StringTie, ...)      from ts tag)           UniProt, genBlastG)
        |                     |                     |                     |
        +---------------------+----------+----------+---------------------+
                                         v
     +----------------------------------------------------------------------+
     |  gmb-run  (one step; = gmb.run_gene_model_builder)                    |
     |    long-read consensus  ->  preflight  ->  build  ->  finalise          |
     |    (per-read input only)    (gate)                   (handover)        |
     +----------------------------------------------------------------------+
                                         v
                       <output-dir>/finalise/  + run_summary.json

      ( reference comparison — EVALUATION ONLY, outside GMB:
        ensembl-genes `annotation-qc pairwise-compare` )
```

---

## What your pipeline must supply

| pipeline artifact | `gmb-run` option | required | notes |
|---|---|---|---|
| `genome_fasta` | `--genome` | **yes** | defines the sequence names every other file must use |
| `backbone` | `--helixer`, `--tiberius` or `--backbone` | **yes** | exactly one |
| `short_read_gtf` | `--scallop`, `--stringtie` | recommended | stranded |
| `protein_alignment_gtf` | `--orthodb`, `--uniprot`, `--genblast` | recommended | support/veto only |
| `long_read_alignments_gtf` | `--long-read-alignments` | optional | **per-read**, stranded; GMB collapses it |
| `long_read_models_gtf` | `--long-read-models` | optional | already-collapsed models (alternative to the above) |
| `preset` | `--preset` | **yes** | fungi: `fungi_latest_settings` |
| `output_directory` | `--output-dir` | **yes** | one per genome |

**Flag names do not determine behaviour — evidence roles do.** See `input_contract.md`.

---

## Calling it

### Shell / workflow manager

```bash
gmb-run --preset fungi_latest_settings \
  --genome "$GENOME" --helixer "$HELIXER" \
  --orthodb "$ORTHODB" --uniprot "$UNIPROT" \
  --scallop "$SCALLOP" --stringtie "$STRINGTIE" \
  ${LONG_READS:+--long-read-alignments "$LONG_READS"} \
  --gene-prefix "$PREFIX" --output-dir "$OUT"
```

Non-interactive. Templates: `examples/run_gene_model_builder.sh` (environment variables),
`examples/submit_gmb.slurm` (Slurm, self-contained).

### Python

```python
from gmb import run_gene_model_builder

result = run_gene_model_builder(
    genome="/data/genome.fa",
    backbone="/data/helixer.gff3", backbone_kind="helixer",
    short_read=["/data/scallop.gtf", "/data/stringtie.gtf"],
    protein_alignment=["/data/orthodb.gtf", "/data/genblast.gtf"],
    long_read_alignments="/data/minimap2_reads.gtf",   # optional, per-read
    preset="fungi_latest_settings",
    output_dir="/work/gmb", gene_prefix="XXGMB",
)
result.ok                    # True only when every stage succeeded and QC passed
result.failed_stage          # None | "longread_consensus" | "preflight" | "build" | "sequence_qc" | "finalise"
result.failure_reason        # human-readable
result.handover_dir          # <output_dir>/finalise
result.annotation            # canonical-tagged GFF3
result.long_read_models      # the consensus GTF that preflight validated and the build used
result.longread_preparation  # action ran / reused / skipped, input sha256, settings, versions
```

`gmb-run` is a thin wrapper around this function, so both apply the same policy.
`run_gene_model_builder` is the **only** supported Python entry point; never import from
`gmb.pipeline` or other internal modules. See `examples/run_gene_model_builder.py`.

---

## Exit codes and gating

`gmb-run` returns: `0` success · `2` invalid arguments · `3` long-read preparation failed ·
`4` preflight failed · `5` build failed · `6` sequence QC failed · `7` finalisation failed ·
`1` unexpected error. Every stage after a failure is skipped, and `run_summary.json` records
`failed_stage` and `failure_reason`. Treat any non-zero code as "do not hand over".

Long-read preparation fails, before any build, when the per-read file is not a GTF, has more
than 5% records without a strand, uses sequence names absent from the genome, is already a GMB
consensus, or when `gmb-longread-consensus` fails or yields no models. Preflight fails a
per-read file given as `--long-read-models` and a long-read track whose splice sites are not
canonical on the declared strand.

If you run the stages individually instead: `gmb-preflight` exits 1 on FAIL;
`gmb-build --validate-fasta` exits 1 when QC fails (the outputs are complete, but must not be
handed over); `gmb-finalise` exits non-zero on error.

---

## What to hand on

`finalise/`, `run_summary.json`, `preflight/preflight_report.txt` and, from `build/`,
`evidence_attribution.tsv`, `resolved_config.yaml` and `run_manifest.json`; with long reads also
`longread_consensus/longread_preparation.json` (`../RUN_FUNGI.md` §4, `output_contract.md`).
**Never hand over `build/` sequence files.**

---

## Resource planning

| | P. falciparum (23 Mb, 14 seqs) | Z. tritici (40 Mb, 21 seqs) |
|---|---|---|
| `gmb-preflight` | < 1 min | ~2 min |
| `gmb-build` wall | 2.1–3.5 h | 10.8–12.0 h |
| `gmb-build` peak RSS | 0.9–1.4 GB | 2.3–2.5 GB |
| `gmb-finalise` | seconds | seconds |

Measured with GMB 1.x/2.0 builds without long reads and up to five builds sharing an
8-core/16 GB host, so wall times are pessimistic. **GMB is single-threaded**:
`--cpus-per-task 2` is enough. Scale memory from **evidence volume**, not genome size —
*Z. tritici* is heavier because its OrthoDB track is 732 MB / 1.69 M alignments. Long-read
consensus streams the per-read file once for validation and once to split it by sequence, then
processes one sequence at a time; the split copy is deleted afterwards.
See `examples/submit_gmb.slurm`.
