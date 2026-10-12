#!/usr/bin/env bash
# Write every summary table of the evaluation into results/ (small text/TSV files).
source "$(dirname "${BASH_SOURCE[0]}")/env.sh"
R="$HERE/results"; mkdir -p "$R"; cd "$HERE"
P="$PYTHON"; RUNS="$WORK/runs"
HC=h2,h4,h6,h8,h10,h13; HA=h16,h19; DC=d1,d3,d5,d7,d9,d12; DA=d14,d18,d21
MAIN=(A_helixer E_all U_all E_lr U_all_lr)
[ -d "$RUNS/U_all_pv" ] && PVARM=(U_all_pv) || PVARM=()

for set in heldout dev; do
  "$P" -I summarise.py "$RUNS" windows.tsv $set "${MAIN[@]}" "${PVARM[@]}" --out "$R/final_$set.tsv" > "$R/final_$set.txt"
  "$P" -I summarise.py "$RUNS" windows.tsv $set A_helixer B_hx C_hx_prot D1_hx_scallop D2_hx_stringtie D_hx_sr D3_hx_lr E_all E_lr \
      --out "$R/ablation_baseline_$set.tsv" > "$R/ablation_baseline_$set.txt"
  "$P" -I summarise.py "$RUNS" windows.tsv $set A_helixer U_B_hx U_C_hx_prot U_D_hx_sr U_D3_hx_lr U_D_hx_sr_lr U_all U_all_lr \
      --out "$R/ablation_updated_$set.tsv" > "$R/ablation_updated_$set.txt"
done
{ for s in "held-out core:$HC" "dev core:$DC"; do
    echo "## ${s%%:*}"; "$P" introns.py "$RUNS" "${s##*:}" A_helixer:consensus E_all:canonical U_all:canonical E_lr:canonical U_all_lr:canonical; echo
  done; } > "$R/introns.txt"
{ for s in "held-out core:$HC" "held-out accessory:$HA" "dev core:$DC" "dev accessory:$DA"; do
    for pair in A_helixer:E_all A_helixer:U_all E_all:U_all A_helixer:U_all_lr E_lr:U_all_lr U_all:U_all_lr; do
      echo "## ${s%%:*}, canonical: ${pair%%:*} -> ${pair##*:} (rows: before, columns: after)"
      "$P" refgene_diff.py "$RUNS" "${pair%%:*}" "${pair##*:}" canonical "${s##*:}" 2>/dev/null | sed -n '1,6p'; echo
    done
  done; } > "$R/transitions.txt"
{ echo "## held-out core"; "$P" model_changes.py "$RUNS" "$WORK/win" "$WORK/w/ref.gff3" $HC E_all U_all E_lr U_all_lr
  echo "## held-out accessory"; "$P" model_changes.py "$RUNS" "$WORK/win" "$WORK/w/ref.gff3" $HA E_all U_all E_lr U_all_lr
} > "$R/model_changes.txt"
"$P" geneset_quality.py "$RUNS" $HC,$HA A_helixer E_all U_all E_lr U_all_lr > "$R/geneset_quality_heldout.txt" 2>/dev/null
"$P" canonical_overrides.py "$WORK" windows.tsv E_all U_all E_lr U_all_lr "${PVARM[@]}" > "$R/canonical_overrides.txt"
{ "$P" -W ignore promotion_audit.py "$WORK" U_all windows.tsv; "$P" -W ignore promotion_audit.py "$WORK" U_all_lr windows.tsv; } > "$R/promotions.txt"
"$P" longread_support.py "$WORK" $HC A_helixer E_all U_all > "$R/longread_support_heldout.txt"
cp "$WORK/provenance.txt" "$R/provenance.txt" 2>/dev/null || "$HERE/provenance.sh" > "$R/provenance.txt"
echo "results written to $R"
