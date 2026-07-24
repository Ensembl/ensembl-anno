#!/bin/bash
# Cluster-ready long-read consensus filtering for a full genome.
#
# Why this is a SLURM array job rather than one local process:
#
# Stage 1 (split-by-seqname) is a single sequential pass over the raw file
# and cannot be parallelized (it's I/O-bound and O(1) memory already; on
# GCA_000002765.3's 14.5GB raw file it took ~90s).
#
# Stage 2 (per-seqname consensus building) is where the real cost is: on
# GCA_000002765.3, seqname "1" (1.09M reads, the SMALLEST chromosome) took
# ~78 minutes even after three rounds of correctness/performance fixes (see
# the run report). The full raw file has ~33x as many reads spread across
# 14 chromosomes, with the largest (chr14) alone carrying ~4.7x as many
# reads as chr1. Naively running Stage 2 for all 14 seqnames sequentially
# on one machine would plausibly take the better part of a day or more --
# impractical and, on a resource-constrained workstation, risks contending
# with other processes for the wall-clock rather than memory (Stage 2's
# memory footprint per seqname is bounded and safe; it is CPU time that
# does not fit a local/interactive session).
#
# Stage 2 is naturally embarrassingly parallel across seqnames (each
# seqname's split file is processed completely independently), so this
# submits it as a SLURM array job: one task per chromosome, all running
# concurrently. Wall-clock time becomes roughly that of the single largest
# chromosome (chr14) rather than the sum of all 14.
#
# Usage:
#   sbatch --array=1-14 longread_consensus_cluster.slurm.sh \
#       <raw_minimap2.gtf> <output_dir> [scallop.gtf] [stringtie.gtf]
#
# Adjust #SBATCH directives (partition, time, mem) for your cluster.

#SBATCH --job-name=lr_consensus
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=12:00:00
#SBATCH --output=lr_consensus_%A_%a.log

set -euo pipefail

RAW_GTF="$1"
OUTPUT_DIR="$2"
SCALLOP="${3:-}"
STRINGTIE="${4:-}"

SEQNAMES=(1 2 3 4 5 6 7 8 9 10 11 12 13 14)
SEQNAME="${SEQNAMES[$((SLURM_ARRAY_TASK_ID - 1))]}"

mkdir -p "$OUTPUT_DIR"

# Stage 1 (split-by-seqname) is cheap and re-run per task deliberately: it
# writes only the ONE seqname this task needs (--seqname flag), so the
# repeated I/O-bound pass costs ~90s x 14 array tasks rather than needing a
# separate shared Stage-1 job with its own dependency wiring. If your
# cluster's shared filesystem makes 14 concurrent full-file reads
# expensive, run Stage 1 once up front instead (drop --seqname, i.e. plain
# `gmb.cli.longread_consensus --input ... --output-dir <shared_split_dir>`
# with only Stage 1's split_by_seqname invoked), then point each array task
# at the pre-split per-seqname file directly.

EXTRA_ARGS=()
if [ -n "$SCALLOP" ]; then
    EXTRA_ARGS+=(--scallop "$SCALLOP")
fi
if [ -n "$STRINGTIE" ]; then
    EXTRA_ARGS+=(--stringtie "$STRINGTIE")
fi

python3 -m gmb.cli.longread_consensus \
    --input "$RAW_GTF" \
    --output-dir "${OUTPUT_DIR}/seqname_${SEQNAME}" \
    --seqname "$SEQNAME" \
    "${EXTRA_ARGS[@]}"

# After all array tasks complete, merge per-seqname outputs:
#   cat "${OUTPUT_DIR}"/seqname_*/minimap2_consensus.gtf > "${OUTPUT_DIR}/minimap2_consensus.gtf"
# Gene/transcript IDs are namespaced by seqname already, so a plain
# concatenation is safe -- no cross-seqname ID collisions.
