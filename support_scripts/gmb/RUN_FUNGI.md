# Running GMB on a fungal genome

Everything you need to annotate a fungal genome with Gene Model Builder (GMB): install it, supply
the evidence, run **one command** with the latest fungal settings, and hand back the finalised,
QC-checked annotation. You do not need any other GMB document to do this.

## Who does what

| you (running GMB) | GMB (automatically) | the annotation lead (afterwards) |
|---|---|---|
| prepare the genome and evidence files (§2) | collapse per-read long-read alignments into transcript models | compare the annotation with references |
| run `gmb-run` (§3) | validate every input before building (preflight) | assess biological annotation quality |
| check that it finished and passed QC (§5) | integrate the evidence and select gene models | decide on further tuning or development |
| return the finalised files and metadata (§4) | finalise, regenerate sequences, run structural QC, record provenance | |

You do **not** need to collapse long reads, choose between configurations, or compare GMB with
Helixer or a reference annotation.

---

## 1. Install

GMB needs Python ≥ 3.10 and four libraries. On a workstation or cluster:

```bash
mamba create -n gmb -c conda-forge python=3.12 'pandas>=2,<3' 'pyranges<=0.1.4' pyyaml numpy
mamba activate gmb
git clone -b feature/gene_model_builder https://github.com/Ensembl/ensembl-anno.git
pip install ./ensembl-anno/support_scripts/gmb
```

If the machine has no internet access during installation, add `--no-build-isolation` to the
`pip install` (setuptools must then be present in the environment).

Check that it works:

```bash
gmb-run --help
gmb-build --list-presets          # must list fungi_latest_settings
git -C ensembl-anno rev-parse HEAD   # record this commit with your results
```

`pyranges` must be `<= 0.1.4`; a newer version is the most common install failure. GMB does not
install Helixer, aligners or assemblers: it only consumes their output.

## 2. Inputs

| input | `gmb-run` option | format | Ensembl anno file | needed |
|---|---|---|---|---|
| genome | `--genome` | FASTA | genome FASTA | **required** |
| Helixer gene predictions | `--helixer` | GFF3 with gene / mRNA / exon / CDS | `<acc>_helixer.gff3`, renamed (below) | **required** |
| protein-to-genome alignments | `--orthodb`, `--uniprot` | GTF (genBlastG) | `<acc>_orthodb_annotation.gtf`, `<acc>_genblast_annotation.gtf` (UniProt) | strongly recommended |
| short-read transcript assemblies | `--scallop`, `--stringtie` | GTF, **stranded** | `<acc>_scallop_annotation.gtf`, `<acc>_stringtie_annotation.gtf` | strongly recommended |
| long-read alignments, one record per read | `--long-read-alignments` | GTF, **stranded** | `<acc>_minimap2_annotation.gtf` | optional; use it when available |

**Sequence names must match.** Every file must use the genome FASTA's sequence names (the first
word of each `>` header). Helixer output for an NCBI assembly uses GenBank accessions (e.g.
`CM001196.1`) instead; rename it with the tool shipped with GMB, using the assembly report:

```bash
python ensembl-anno/support_scripts/gmb/tools/remap_helixer.py \
  --input <acc>_helixer.gff3 --assembly-report <acc>_assembly_report.txt \
  --output helixer.gff3
```

**Strand must already be set.** GMB uses the strand written in each GTF and never infers it. The
Ensembl anno pipeline sets long-read strand from the Minimap2 `ts` tag. GMB stops if more than
5% of long-read records have no strand, and preflight fails a long-read track whose splice sites
are not canonical on the declared strand (the sign of wrongly stranded reads).

**Long reads are supplied as they come from the aligner** (one GTF transcript per read). GMB
collapses them into consensus transcript models itself. If you only have long-read transcript
models that were already collapsed by another tool, pass them with `--long-read-models` instead;
never pass per-read alignments with that option (GMB stops if you do).

## 3. Run

```bash
EV=/data/GCA_000219625.1          # prefix of the Ensembl anno evidence files
gmb-run \
  --preset fungi_latest_settings \
  --genome /data/genome.fa \
  --helixer /data/helixer.gff3 \
  --orthodb  ${EV}_orthodb_annotation.gtf \
  --uniprot  ${EV}_genblast_annotation.gtf \
  --scallop  ${EV}_scallop_annotation.gtf \
  --stringtie ${EV}_stringtie_annotation.gtf \
  --long-read-alignments ${EV}_minimap2_annotation.gtf \
  --gene-prefix ZTRITICI \
  --output-dir /work/gmb_ztritici
```

**Without long reads**, leave out the `--long-read-alignments` line; nothing else changes. Leave
out any other evidence you do not have in the same way; only `--genome` and `--helixer` are
required.

`gmb-run` collapses the long reads (when given), validates all inputs, builds the gene models,
finalises them and runs QC, writing everything under `--output-dir`. It never asks questions,
so it runs unattended in a batch job. `--preset fungi_latest_settings` selects the latest fungal
settings; use it for every fungal genome. Progress goes to the screen and to
`<output-dir>/logs/gmb-run.log`.

**On a Slurm cluster** use the bundled template (edit the `#SBATCH` lines for your site):

```bash
sbatch --export=ALL,GENOME=/data/genome.fa,HELIXER=/data/helixer.gff3,\
ORTHODB=${EV}_orthodb_annotation.gtf,UNIPROT=${EV}_genblast_annotation.gtf,\
SCALLOP=${EV}_scallop_annotation.gtf,STRINGTIE=${EV}_stringtie_annotation.gtf,\
LONG_READ_ALIGNMENTS=${EV}_minimap2_annotation.gtf,\
GENE_PREFIX=ZTRITICI,OUT=/work/gmb_ztritici \
  ensembl-anno/support_scripts/gmb/examples/submit_gmb.slurm
```

**Resources.** GMB is single-threaded. A 40 Mb genome with the full Ensembl anno evidence
(a 730 MB protein alignment file) has taken 11–12 h and ~2.5 GB of memory; request 24 h and
16 GB for the first run of a new genome. Collapsing ~2 M long reads adds minutes, not hours.

**Running again.** Rerunning into the same `--output-dir` replaces the previous results. The
long-read consensus is reused only if the long-read file, its settings and the GMB version are
all unchanged (the log says `REUSED`), and is rebuilt otherwise.

## 4. Outputs

```
<output-dir>/
  finalise/                                     <- THE RESULT
    canonical/consensus.canonical_annotated.gff3   annotation: all genes and isoforms; one
                                                   transcript per gene tagged Ensembl_canonical
    consensus.gff3                                 the same gene set without canonical tags
    cdna.fa  cds.fa  prot.fa                       sequences regenerated from the final GFF3
    canonical/canonical_transcripts.tsv            the canonical transcript of each gene
    fasta_qc_report.json  utr_qc_report.json       structural QC (must pass)
    handover_manifest.json                         file inventory with checksums
  build/       evidence_attribution.tsv  resolved_config.yaml  run_manifest.json  ...
  preflight/   preflight_report.txt  preflight_report.json
  longread_consensus/   minimap2_consensus.gtf  longread_preparation.json  ...   (with long reads)
  logs/        gmb-run.log  longread_consensus.log  preflight.log  build.log  finalise.log
  run_summary.json      overall status and the paths above
```

`build/` and `longread_consensus/` are intermediate; `finalise/` holds the results. Do not use
the GFF3 or FASTA files in `build/`.

**Return for assessment**:

```bash
cd <output-dir>
tar czf gmb_results.tar.gz finalise run_summary.json preflight/preflight_report.txt \
  build/evidence_attribution.tsv build/resolved_config.yaml build/resolved_config_sha256 \
  build/run_manifest.json $( [ -d longread_consensus ] && echo longread_consensus/longread_preparation.json )
```

That is the annotation and sequences, the QC reports, why each model was chosen, the exact
configuration, the GMB version and input checksums, and the long-read preparation record. Add the
GMB commit you recorded at installation and the names of the evidence files you used.

## 5. Check that it worked

`gmb-run` exits 0 and its last lines say `GMB run complete` and list the files to share. Then:

```bash
python -c "import json; s=json.load(open('<output-dir>/run_summary.json')); print('ok:', s['ok'], '| QC passed:', s['qc_passed'])"
```

Both must be `True`. The preflight verdict may be `WARN`; that is normal (on *Z. tritici*, for
example, a few unstranded short-read rows and very long protein alignments are reported and
handled). Biological accuracy is assessed afterwards by the annotation lead, not here.

## 6. If it fails

`gmb-run` stops at the first failing stage, explains why in its last lines and in
`run_summary.json` (`failed_stage`, `failure_reason`), and returns:

| exit | stage | look at | usual cause |
|---|---|---|---|
| 2 | arguments | the error message | a missing file, a typo in an option, an unknown preset |
| 3 | long-read preparation | `logs/gmb-run.log`, `logs/longread_consensus.log` | long-read file has no strand, uses other sequence names, is not a GTF, or is already a GMB consensus (then use `--long-read-models`) |
| 4 | preflight | `preflight/preflight_report.txt` (the FAIL lines) | sequence names that do not match the genome (rename Helixer, §2), an empty or unreadable file, wrongly stranded long reads, per-read alignments given as `--long-read-models` |
| 5 | build | `logs/build.log` | out of memory or time (raise the job limits), a corrupt input |
| 6 | sequence QC | `build/fasta_qc_report.json` | report it with the files below; do not edit outputs |
| 7 | finalisation | `logs/finalise.log` | report it |

To report a problem, send `run_summary.json`, the `logs/` directory and
`preflight/preflight_report.txt`.

---

*Developers and advanced users:* the individual commands (`gmb-preflight`, `gmb-build`,
`gmb-finalise`, `gmb-longread-consensus`), the Python API and the scientific background are
described in [README.md](README.md) and [docs/](docs/). The historical fungal baseline
(`--preset fungi`) is kept only to reproduce earlier results.
