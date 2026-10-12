#!/usr/bin/env bash
# Prepare evaluation inputs from the ORIGINAL data, which is only read:
#   $WORK/w/   genome.fa (headers trimmed), ref.gff3, helixer_wg.gff3 (sequence names remapped)
#   $WORK/win/<window>/   whole-model subsets of every track + whole-chromosome genome
#                         + minimap2_consensus.gtf (per-read long reads collapsed, fungi preset)
source "$(dirname "${BASH_SOURCE[0]}")/env.sh"
: "${INPUT_DATA:?}" "${REFERENCE_DATA:?}"
mkdir -p "$WORK/w"
gzip -dc "$REFERENCE_DATA/Zymoseptoria_tritici_GCA_000219625.1.gff3.bgz" > "$WORK/w/ref.gff3"
gzip -dc "$REFERENCE_DATA/Zymoseptoria_tritici_GCA_000219625.1.fa.bgz" \
  | awk '/^>/{print $1; next}{print}' > "$WORK/w/genome.fa"
"$PYTHON" "$GMB_DIR/tools/remap_helixer.py" --input "$INPUT_DATA/GCA_000219625.1_helixer.gff3" \
  --assembly-report "$INPUT_DATA/GCA_000219625.1_MYCGR_v2.0_assembly_report.txt" \
  --output "$WORK/w/helixer_wg.gff3"
"$PYTHON" -I "$HERE/subset_windows.py" "$INPUT_DATA" "$WORK/w/helixer_wg.gff3" "$WORK/w/genome.fa" \
  "$HERE/windows.tsv" "$WORK/win"
mkdir -p "$WORK/lr"
for w in $(awk 'NR>1{print $1}' "$HERE/windows.tsv"); do
  "$GMB_BIN/gmb-longread-consensus" --input "$WORK/win/$w/minimap2.gtf" --preset fungi \
    --output-dir "$WORK/lr/$w" --no-log-file > "$WORK/lr/$w.log" 2>&1 || { cat "$WORK/lr/$w.log"; exit 1; }
  cp "$WORK/lr/$w/minimap2_consensus.gtf" "$WORK/win/$w/minimap2_consensus.gtf"
done
echo "inputs ready in $WORK"
