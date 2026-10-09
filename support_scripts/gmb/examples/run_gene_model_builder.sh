#!/usr/bin/env bash
# =============================================================================
# Gene Model Builder — production wrapper
#
# Runs preflight -> build -> finalise and prints the handover directory.
# Every path is an environment variable; nothing here is machine-specific.
#
#   REQUIRED  GENOME BACKBONE OUT
#   OPTIONAL  BACKBONE_KIND PRESET CONFIG GENE_PREFIX
#             SCALLOP STRINGTIE MINIMAP2 ORTHODB UNIPROT GENBLAST
#
# Example:
#   GENOME=/data/genome.fa \
#   BACKBONE=/data/backbone.gff3 BACKBONE_KIND=helixer \
#   SCALLOP=/data/scallop.gtf STRINGTIE=/data/stringtie.gtf \
#   ORTHODB=/data/orthodb.gtf \
#   PRESET=fungi GENE_PREFIX=XXGMB OUT=/work/gmb \
#   ./run_gene_model_builder.sh
# =============================================================================
set -euo pipefail

: "${GENOME:?set GENOME to the genome FASTA}"
: "${BACKBONE:?set BACKBONE to the ab initio backbone annotation}"
: "${OUT:?set OUT to the output directory}"
: "${BACKBONE_KIND:=helixer}"     # helixer | tiberius | backbone (any predictor)
: "${PRESET:=standard}"           # standard | apicomplexa | fungi
: "${GENE_PREFIX:=GMB}"
: "${CONFIG:=}"                   # space-separated YAML overlays

mkdir -p "$OUT"/logs

# --- assemble the evidence arguments once, used by preflight and build -------
EVIDENCE=( "--${BACKBONE_KIND}" "$BACKBONE" )
[ -n "${SCALLOP:-}"   ] && EVIDENCE+=( --scallop   "$SCALLOP"   )
[ -n "${STRINGTIE:-}" ] && EVIDENCE+=( --stringtie "$STRINGTIE" )
[ -n "${MINIMAP2:-}"  ] && EVIDENCE+=( --minimap2  "$MINIMAP2"  )
[ -n "${ORTHODB:-}"   ] && EVIDENCE+=( --orthodb   "$ORTHODB"   )
[ -n "${UNIPROT:-}"   ] && EVIDENCE+=( --uniprot   "$UNIPROT"   )
[ -n "${GENBLAST:-}"  ] && EVIDENCE+=( --genblast  "$GENBLAST"  )

CONFIG_ARGS=()
for c in $CONFIG; do CONFIG_ARGS+=( --config "$c" ); done

# --- 1. preflight ------------------------------------------------------------
# Validates the bundle before the expensive build: seqid compatibility, ID
# collisions, strand completeness, per-track splice quality, resolved roles and
# weights, and which optional policies will actually fire.
# Exits 1 on FAIL. Add --allow-fail to proceed anyway (the verdict is still
# recorded in preflight_report.json and in the run manifest).
echo "[1/3] preflight"
gmb-preflight \
  --preset "$PRESET" "${CONFIG_ARGS[@]+"${CONFIG_ARGS[@]}"}" \
  --genome "$GENOME" "${EVIDENCE[@]}" \
  --output-dir "$OUT/preflight"

# --- 2. build ----------------------------------------------------------------
# NOTE: with --validate-fasta, a non-zero exit is a QC VERDICT, not a crash --
# the annotation is fully written. Inspect build/fasta_qc_report.json.
echo "[2/3] build"
set +e
gmb-build \
  --preset "$PRESET" "${CONFIG_ARGS[@]+"${CONFIG_ARGS[@]}"}" \
  --genome "$GENOME" "${EVIDENCE[@]}" \
  --gene-prefix "$GENE_PREFIX" \
  --output-dir "$OUT/build" \
  --validate-fasta
build_rc=$?
set -e
if [ "$build_rc" -ne 0 ]; then
  echo "gmb-build exited $build_rc — QC did not pass." >&2
  echo "Outputs are complete; see $OUT/build/fasta_qc_report.json" >&2
  echo "Provenance for the failing run: $OUT/build/run_manifest.json" >&2
  echo "DO NOT hand over this build (docs/qc.md)." >&2
  exit "$build_rc"
fi

# --- 3. finalise -------------------------------------------------------------
# Regenerates cDNA/CDS/protein FROM THE FINAL GFF3, re-runs QC, selects one
# canonical transcript per gene, writes the handover manifest.
# No --preset/--config needed: finalise automatically runs under the build's
# own resolved_config.yaml, so it cannot finalise under settings that never
# applied to the build.
echo "[3/3] finalise"
gmb-finalise \
  --build-dir "$OUT/build" \
  --genome "$GENOME" \
  --output-dir "$OUT/finalise"

cat <<SUMMARY

Done. Hand over:

  $OUT/finalise/canonical/consensus.canonical_annotated.gff3   primary annotation
  $OUT/finalise/consensus.gff3                                 all isoforms
  $OUT/finalise/{cdna,cds,prot}.fa                             sequences
  $OUT/finalise/canonical/canonical_transcripts.tsv            canonical decisions
  $OUT/finalise/{fasta_qc_report,utr_qc_report}.json           QC evidence
  $OUT/finalise/handover_manifest.json                         output inventory

  $OUT/build/evidence_attribution.tsv                          why each model won
  $OUT/build/resolved_config.yaml                              exact configuration
  $OUT/build/run_manifest.json                                 full provenance

Do NOT hand over the FASTA files in $OUT/build — finalise/ is where they are
regenerated from the final GFF3. See docs/output_contract.md.
SUMMARY
