#!/usr/bin/env bash
# Run every arm of arms.tsv (or the arms named as arguments) in every window, in parallel.
# Protein-validation arms (configs containing pv) run only when named.
source "$(dirname "${BASH_SOURCE[0]}")/env.sh"
[ -n "$AQC" ] || { echo "annotation-qc not found; set AQC" >&2; exit 1; }
want=" $* "
mkdir -p "$WORK/runs"
{
  tail -n +2 "$HERE/arms.tsv" | while IFS=$'\t' read -r arm tracks configs _desc; do
    [ $# -gt 0 ] && [[ "$want" != *" $arm "* ]] && continue
    # Protein-validation arms run only when named (and PV_CONFIG is set): Psauron is
    # multi-threaded, so run them with a low JOBS, e.g. JOBS=2 run_matrix.sh U_all_pv.
    [[ ",$configs," == *",pv,"* ]] && { [ -z "${PV_CONFIG:-}" ] || [ $# -eq 0 ]; } && continue
    tail -n +2 "$HERE/windows.tsv" | while IFS=$'\t' read -r w _set chrom s e _cls; do
      printf '%s\t%s\t%s\t%s\t%s\n' "$arm" "$tracks" "$configs" "$w" "$chrom:$s-$e"
    done
  done
} > "$WORK/jobs.tsv"
echo "$(wc -l < "$WORK/jobs.tsv") jobs"
tr '\t' '\n' < "$WORK/jobs.tsv" | xargs -P "$JOBS" -n 5 bash "$HERE/run_arm.sh" | tee "$WORK/runs/matrix.log"
"$HERE/provenance.sh" > "$WORK/provenance.txt"
