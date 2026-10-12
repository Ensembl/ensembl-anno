#!/usr/bin/env bash
# =============================================================================
# Gene Model Builder — environment-variable wrapper around `gmb-run`
#
# gmb-run does everything: automatic long-read consensus (when per-read
# alignments are given), preflight, build, finalise and QC. This wrapper only
# turns environment variables into its options, for workflow managers and the
# Slurm template (submit_gmb.slurm). Nothing here is machine-specific.
#
#   REQUIRED  GENOME HELIXER OUT
#   OPTIONAL  ORTHODB UNIPROT GENBLAST SCALLOP STRINGTIE
#             LONG_READ_ALIGNMENTS   per-read long-read GTF (collapsed automatically)
#             LONG_READ_MODELS       already-collapsed long-read models (rare)
#             GENE_PREFIX (default GMB)
#             PRESET (default fungi_latest_settings: the recommended fungal settings)
#
# Example (fungal genome; see RUN_FUNGI.md):
#   GENOME=/data/genome.fa HELIXER=/data/helixer.gff3 \
#   ORTHODB=/data/orthodb.gtf UNIPROT=/data/genblast.gtf \
#   SCALLOP=/data/scallop.gtf STRINGTIE=/data/stringtie.gtf \
#   LONG_READ_ALIGNMENTS=/data/minimap2_annotation.gtf \
#   GENE_PREFIX=XXGMB OUT=/work/gmb ./run_gene_model_builder.sh
#
# Exit status is gmb-run's: 0 ok, 2 arguments, 3 long-read preparation,
# 4 preflight, 5 build, 6 sequence QC, 7 finalise.
# =============================================================================
set -euo pipefail

: "${GENOME:?set GENOME to the genome FASTA}"
: "${HELIXER:?set HELIXER to the Helixer GFF3 (sequence names as in GENOME)}"
: "${OUT:?set OUT to the output directory}"
: "${PRESET:=fungi_latest_settings}"
: "${GENE_PREFIX:=GMB}"

ARGS=( --preset "$PRESET" --genome "$GENOME" --helixer "$HELIXER"
       --gene-prefix "$GENE_PREFIX" --output-dir "$OUT" )
[ -n "${ORTHODB:-}"   ] && ARGS+=( --orthodb   "$ORTHODB"   )
[ -n "${UNIPROT:-}"   ] && ARGS+=( --uniprot   "$UNIPROT"   )
[ -n "${GENBLAST:-}"  ] && ARGS+=( --genblast  "$GENBLAST"  )
[ -n "${SCALLOP:-}"   ] && ARGS+=( --scallop   "$SCALLOP"   )
[ -n "${STRINGTIE:-}" ] && ARGS+=( --stringtie "$STRINGTIE" )
[ -n "${LONG_READ_ALIGNMENTS:-}" ] && ARGS+=( --long-read-alignments "$LONG_READ_ALIGNMENTS" )
[ -n "${LONG_READ_MODELS:-}"     ] && ARGS+=( --long-read-models     "$LONG_READ_MODELS" )

exec gmb-run "${ARGS[@]}"
