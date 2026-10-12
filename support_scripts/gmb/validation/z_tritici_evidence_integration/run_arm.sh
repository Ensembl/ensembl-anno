#!/usr/bin/env bash
# Build + finalise one arm in one window and evaluate the FINALISED outputs three ways:
#   cmp_consensus  all isoforms       fin/consensus.gff3
#   cmp_primary    build primary .t1  fin/primary.gff3 (derived from fin/consensus.gff3)
#   cmp_canonical  handover canonical fin/canonical/consensus.canonical_annotated.gff3
# Usage: run_arm.sh <arm> <tracks> <configs> <window> <region>
#   tracks : comma list of helixer,scallop,stringtie,protein,longread ("-" = raw Helixer, no GMB)
#   configs: comma list of overlay,pv or "-"
source "$(dirname "${BASH_SOURCE[0]}")/env.sh"
ARM=$1; TRACKS=$2; CONFIGS=$3; WN=$4; REGION=$5
D="$WORK/win/$WN"; O="$WORK/runs/$ARM/$WN"; mkdir -p "$O"
cmp_one() {  # <query gff3> <outdir> [extra annotation-qc args]
  local q=$1 out=$2; shift 2
  "$AQC" pairwise-compare --query "$q" --reference "$WORK/w/ref.gff3" --genome "$WORK/w/genome.fa" \
    --evaluation-mode protein_coding --reference-transcript-biotypes protein_coding \
    --region "$REGION" --outdir "$out" "$@" > "$out.log" 2>&1
}
if [ "$TRACKS" = "-" ]; then
  cmp_one "$D/helixer.gff3" "$O/cmp_consensus"
  ln -sfn cmp_consensus "$O/cmp_primary"; ln -sfn cmp_consensus "$O/cmp_canonical"
  echo "done $ARM $WN"; exit 0
fi
ARGS=(--genome "$D/genome.fa")
[[ ",$TRACKS," == *",helixer,"* ]]   && ARGS+=(--helixer "$D/helixer.gff3")
[[ ",$TRACKS," == *",scallop,"* ]]   && ARGS+=(--scallop "$D/scallop.gtf")
[[ ",$TRACKS," == *",stringtie,"* ]] && ARGS+=(--stringtie "$D/stringtie.gtf")
[[ ",$TRACKS," == *",protein,"* ]]   && ARGS+=(--orthodb "$D/orthodb.gtf" --uniprot "$D/genblast.gtf")
[[ ",$TRACKS," == *",longread,"* ]]  && ARGS+=(--minimap2 "$D/minimap2_consensus.gtf")
[[ ",$CONFIGS," == *",overlay,"* ]]  && ARGS+=(--config "$OVERLAY")
[[ ",$CONFIGS," == *",pv,"* ]]       && ARGS+=(--config "${PV_CONFIG:?pv arm needs PV_CONFIG}")
rm -rf "$O/build" "$O/fin"
TIME=(/usr/bin/time -l); [ "$(uname)" = Linux ] && TIME=(/usr/bin/time -v)
"${TIME[@]}" "$GMB_BIN/gmb-build" --preset fungi "${ARGS[@]}" --gene-prefix ZT \
  --output-dir "$O/build" --no-log-file > "$O/build.log" 2>&1
"$GMB_BIN/gmb-finalise" --build-dir "$O/build" --genome "$D/genome.fa" --output-dir "$O/fin" \
  > "$O/fin.log" 2>&1
"$PYTHON" -I "$HERE/primary_only.py" "$O/fin/consensus.gff3" "$O/fin/primary.gff3"
cmp_one "$O/fin/consensus.gff3" "$O/cmp_consensus"
cmp_one "$O/fin/primary.gff3" "$O/cmp_primary"
cmp_one "$O/fin/canonical/consensus.canonical_annotated.gff3" "$O/cmp_canonical" \
  --query-transcript-selection canonical
echo "done $ARM $WN"
