#!/usr/bin/env bash
# =============================================================================
# Gene Model Builder — SLURM template (fungal genomes: see RUN_FUNGI.md)
#
# Runs `gmb-run` with the latest fungal settings. Long-read consensus,
# preflight, build, finalise and QC all happen inside the one job. Self-contained
# (sbatch copies job scripts to a spool directory, so it calls no sibling
# script). Nothing here is site-specific.
#
# Submit with every input passed through --export:
#
#   sbatch \
#     --export=ALL,\
# GENOME=/data/genome.fa,HELIXER=/data/helixer.gff3,\
# ORTHODB=/data/orthodb.gtf,UNIPROT=/data/genblast.gtf,\
# SCALLOP=/data/scallop.gtf,STRINGTIE=/data/stringtie.gtf,\
# LONG_READ_ALIGNMENTS=/data/minimap2_annotation.gtf,\
# GENE_PREFIX=XXGMB,OUT=/work/gmb \
#     examples/submit_gmb.slurm
#
# Only GENOME, HELIXER and OUT are required; leave out evidence you do not have.
# =============================================================================
#SBATCH --job-name=gmb
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=24:00:00
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.err

# --- Resource guidance -------------------------------------------------------
# MEASURED on two genomes (GMB 1.x/2.0 builds, without long reads), with up to
# five builds sharing one 8-core/16 GB host (wall times pessimistic, memory
# realistic):
#
#                         P. falciparum        Z. tritici
#   genome                23 MB, 14 seqs       40 MB, 21 seqs
#   largest evidence      23 MB long-read      732 MB OrthoDB (1.69 M transcripts)
#   gmb-build wall        2.1 - 3.5 h          10.8 - 12.0 h
#   gmb-build peak RSS    0.9 - 1.4 GB         2.3 - 2.5 GB
#   gmb-finalise          3 - 6 s, < 0.5 GB    6 s, ~0.4 GB
#
# These are EXAMPLES, not requirements.
#
#   --cpus-per-task 2   GMB's hot path (scoring and isoform selection) is
#                       single-threaded. More cores buy nothing for GMB itself;
#                       the second absorbs I/O and any DIAMOND subprocess.
#   --mem 16G           8G is enough for a compact genome with modest evidence.
#                       Scale from EVIDENCE VOLUME, not genome size: Z. tritici
#                       is the heavier job because of its protein track, not its
#                       genome. Long-read consensus holds one sequence's reads
#                       at a time.
#   --time 24:00:00     generous for a compact genome on a dedicated node.
#
# After the first run, size from your own measurements:
#   sacct -j <jobid> --format=JobID,Elapsed,MaxRSS,MaxVMSize,State
# GMB also records runtime and peak RSS in build/run_manifest.json.
# -----------------------------------------------------------------------------

set -euo pipefail

: "${GENOME:?}"; : "${HELIXER:?}"; : "${OUT:?}"
: "${PRESET:=fungi_latest_settings}"; : "${GENE_PREFIX:=GMB}"

# --- environment -------------------------------------------------------------
# Activate the environment GMB is installed into. Pick ONE:
#   conda/mamba:  source /shared/conda/etc/profile.d/conda.sh && conda activate gmb
#   venv:         module load python/3.12 && source /shared/venvs/gmb/bin/activate
#   site-install: nothing to do

mkdir -p "$OUT/logs"

{
  echo "date=$(date -Is)"
  echo "hostname=$(hostname)"
  echo "slurm_job_id=${SLURM_JOB_ID:-NA}"
  echo "cpus_per_task=${SLURM_CPUS_PER_TASK:-NA}"
  echo "mem_request=${SLURM_MEM_PER_NODE:-${SLURM_MEM_PER_CPU:-NA}}"
} > "$OUT/logs/slurm_environment.txt"

ARGS=( --preset "$PRESET" --genome "$GENOME" --helixer "$HELIXER"
       --gene-prefix "$GENE_PREFIX" --output-dir "$OUT" )
[ -n "${ORTHODB:-}"   ] && ARGS+=( --orthodb   "$ORTHODB"   )
[ -n "${UNIPROT:-}"   ] && ARGS+=( --uniprot   "$UNIPROT"   )
[ -n "${GENBLAST:-}"  ] && ARGS+=( --genblast  "$GENBLAST"  )
[ -n "${SCALLOP:-}"   ] && ARGS+=( --scallop   "$SCALLOP"   )
[ -n "${STRINGTIE:-}" ] && ARGS+=( --stringtie "$STRINGTIE" )
[ -n "${LONG_READ_ALIGNMENTS:-}" ] && ARGS+=( --long-read-alignments "$LONG_READ_ALIGNMENTS" )
[ -n "${LONG_READ_MODELS:-}"     ] && ARGS+=( --long-read-models     "$LONG_READ_MODELS" )

# Exit status: 0 ok, 2 arguments, 3 long-read preparation, 4 preflight,
# 5 build, 6 sequence QC, 7 finalise (RUN_FUNGI.md section 6).
set +e
gmb-run "${ARGS[@]}"
rc=$?
set -e

# GMB writes its own provenance ($OUT/run_summary.json, build/run_manifest.json);
# this adds the scheduler's view of the same run.
if command -v sacct >/dev/null 2>&1 && [ -n "${SLURM_JOB_ID:-}" ]; then
  sacct -j "$SLURM_JOB_ID" \
    --format=JobID,JobName,Elapsed,MaxRSS,MaxVMSize,ReqMem,State \
    > "$OUT/logs/sacct.txt" 2>/dev/null || true
fi

exit "$rc"
