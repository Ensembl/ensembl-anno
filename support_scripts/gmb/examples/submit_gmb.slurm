#!/usr/bin/env bash
# =============================================================================
# Gene Model Builder — SLURM template
#
# Wraps examples/run_gene_model_builder.sh. Nothing here is site-specific.
#
# Submit with every input passed through --export:
#
#   sbatch \
#     --export=ALL,\
# GENOME=/data/genome.fa,\
# BACKBONE=/data/backbone.gff3,BACKBONE_KIND=helixer,\
# SCALLOP=/data/scallop.gtf,STRINGTIE=/data/stringtie.gtf,\
# ORTHODB=/data/orthodb.gtf,\
# PRESET=fungi,GENE_PREFIX=XXGMB,OUT=/work/gmb \
#     examples/submit_gmb.slurm
#
# MINIMAP2, UNIPROT, GENBLAST and CONFIG are optional.
# =============================================================================
#SBATCH --job-name=gmb
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=24:00:00
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.err

# --- Resource guidance -------------------------------------------------------
# MEASURED on two genomes, with up to five builds sharing one 8-core/16 GB host
# (so wall times are pessimistic and memory figures realistic):
#
#                         P. falciparum        Z. tritici
#   genome                23 MB, 14 seqs       40 MB, 21 seqs
#   largest evidence      23 MB long-read      732 MB OrthoDB (1.69 M transcripts)
#   gmb-build wall        2.1 - 3.5 h          10.8 - 12.0 h
#   gmb-build peak RSS    0.9 - 1.4 GB         2.3 - 2.5 GB
#   gmb-finalise          3 - 6 s, < 0.5 GB    6 s, ~0.4 GB
#   gmb-compare peak RSS  ~0.6 GB              4.2 GB   <-- highest of any stage
#
# These are EXAMPLES, not requirements.
#
#   --cpus-per-task 2   GMB's hot path (scoring and isoform selection) is
#                       single-threaded. More cores buy nothing for GMB itself;
#                       the second absorbs I/O and any DIAMOND subprocess.
#   --mem 16G           8G is enough for a compact genome with modest evidence.
#                       Scale from EVIDENCE VOLUME, not genome size: Z. tritici
#                       is the heavier job because of its protein track, not its
#                       genome. If you also run gmb-compare with a large protein
#                       track, size for that stage -- it peaked highest.
#   --time 24:00:00     generous for a compact genome on a dedicated node.
#
# After the first run, size from your own measurements:
#   sacct -j <jobid> --format=JobID,Elapsed,MaxRSS,MaxVMSize,State
# GMB also records runtime and peak RSS in build/run_manifest.json.
# -----------------------------------------------------------------------------

set -euo pipefail

: "${GENOME:?}"; : "${BACKBONE:?}"; : "${OUT:?}"
: "${PRESET:=standard}"; : "${GENE_PREFIX:=GMB}"; : "${BACKBONE_KIND:=helixer}"

# --- environment -------------------------------------------------------------
# Activate the environment GMB is installed into. Pick ONE:
#   conda/mamba:  source /shared/conda/etc/profile.d/conda.sh && conda activate gmb
#   venv:         module load python/3.12 && source /shared/venvs/gmb/bin/activate
#   site-install: nothing to do
#
# DIAMOND and Psauron are NOT installed by GMB and generally need their own
# environments; give GMB their absolute paths in a config overlay instead of
# trying to co-install them.

mkdir -p "$OUT/logs"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

{
  echo "date=$(date -Is)"
  echo "hostname=$(hostname)"
  echo "slurm_job_id=${SLURM_JOB_ID:-NA}"
  echo "cpus_per_task=${SLURM_CPUS_PER_TASK:-NA}"
  echo "mem_request=${SLURM_MEM_PER_NODE:-${SLURM_MEM_PER_CPU:-NA}}"
} > "$OUT/logs/slurm_environment.txt"

"$HERE/run_gene_model_builder.sh"
rc=$?

# GMB writes its own provenance to $OUT/build/run_manifest.json; this adds the
# scheduler's view of the same run.
if command -v sacct >/dev/null 2>&1 && [ -n "${SLURM_JOB_ID:-}" ]; then
  sacct -j "$SLURM_JOB_ID" \
    --format=JobID,JobName,Elapsed,MaxRSS,MaxVMSize,ReqMem,State \
    > "$OUT/logs/sacct.txt" 2>/dev/null || true
fi

exit "$rc"
