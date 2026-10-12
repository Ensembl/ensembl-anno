#!/usr/bin/env bash
# End-to-end smoke test on the bundled 500 kb Z. tritici fixture (chr1:1-500,000):
#
#   1. the recommended workflow, exactly as RUN_FUNGI.md documents it:
#      gmb-run --preset fungi_latest_settings ... --long-read-alignments <per-read GTF>
#      run twice (determinism; the second run must REUSE the long-read consensus)
#   2. the same without long reads
#   3. the historical `fungi` preset through the individual stages, optionally compared
#      byte-for-byte with a baseline GMB checkout
#
# Every handover is checked independently (check_handover.py). Run it from any directory.
#
# Usage: smoke_fixture.sh <outdir> [<baseline ensembl-anno checkout>]
#   GMB_BIN   directory with gmb-* and python (default: from PATH)
# Exit 1 on any failure.
set -euo pipefail
OUT=$1; BASELINE=${2:-}
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
GMB_DIR="$(cd "$HERE/.." && pwd)"
GMB_BIN="${GMB_BIN:-$(dirname "$(command -v gmb-run)")}"
PY="$GMB_BIN/python"; [ -x "$PY" ] || PY="$(command -v python3)"
FX="$GMB_DIR/tests/fixtures/z_tritici_region1"
EVID=(--genome "$FX/genome.fa" --helixer "$FX/helixer_remapped.gff3"
      --orthodb "$FX/orthodb_geneset.gtf" --uniprot "$FX/uniprot_geneset.gtf"
      --scallop "$FX/scallop_geneset.gtf" --stringtie "$FX/stringtie_geneset.gtf")
mkdir -p "$OUT"; OUT="$(cd "$OUT" && pwd)"; fail=0
check() { "$PY" "$HERE/check_handover.py" "$1/finalise" "$1/build" || fail=1; }

echo "--- 1. recommended workflow with per-read long reads (gmb-run)"
for n in latest_lr latest_lr_rerun; do
  rm -rf "$OUT/$n"; [ $n = latest_lr_rerun ] && cp -R "$OUT/latest_lr" "$OUT/$n"
  "$GMB_BIN/gmb-run" --preset fungi_latest_settings "${EVID[@]}" \
    --long-read-alignments "$FX/minimap2_reads.gtf" \
    --gene-prefix ZT --output-dir "$OUT/$n" > "$OUT/$n.log" 2>&1 \
    || { echo "FAIL: gmb-run $n exited $?"; fail=1; }
  check "$OUT/$n"
done
grep -q "Long-read consensus: RAN" "$OUT/latest_lr.log" || { echo "FAIL: consensus not run"; fail=1; }
grep -q "Long-read consensus: REUSED" "$OUT/latest_lr_rerun.log" \
  && echo "rerun reused the long-read consensus" || { echo "FAIL: consensus not reused"; fail=1; }
for f in finalise/consensus.gff3 finalise/canonical/consensus.canonical_annotated.gff3 \
         finalise/prot.fa build/evidence_attribution.tsv longread_consensus/minimap2_consensus.gtf; do
  cmp -s "$OUT/latest_lr/$f" "$OUT/latest_lr_rerun/$f" \
    && echo "identical on rerun: $f" || { echo "DIFFERS on rerun: $f"; fail=1; }
done

echo "--- 2. recommended workflow without long reads"
rm -rf "$OUT/latest"
"$GMB_BIN/gmb-run" --preset fungi_latest_settings "${EVID[@]}" --gene-prefix ZT \
  --output-dir "$OUT/latest" > "$OUT/latest.log" 2>&1 || { echo "FAIL: gmb-run latest"; fail=1; }
check "$OUT/latest"
[ -d "$OUT/latest/longread_consensus" ] && { echo "FAIL: consensus without long reads"; fail=1; }

echo "--- provenance"
"$PY" - "$OUT" <<'PY' || fail=1
import json, sys, yaml, pathlib
out = pathlib.Path(sys.argv[1]); ok = True
for run, lr in (("latest_lr", True), ("latest", False)):
    s = json.load(open(out / run / "run_summary.json"))
    c = yaml.safe_load(open(out / run / "build" / "resolved_config.yaml"))
    got = (c["preset"], c["scoring"]["primary_selection"], c["scoring"]["isoform_cds_overlap"],
           c["scoring"]["min_cds_bp"], c["canonical_selection"]["prefer_build_primary"])
    print(f"{run}: ok={s['ok']} stages={s['stages_run']} settings={got}")
    ok &= s["ok"] and got == ("fungi_latest_settings", "junction_supported", "drop", 150, True)
    if lr:
        p = json.load(open(out / run / "longread_consensus" / "longread_preparation.json"))
        r = p["consensus_settings"]["resolved"]
        print(f"  long reads: {p['consensus_models']} models from {p['input_scan']['transcript_rows']} "
              f"reads; preset={p['consensus_settings']['preset']} span={r['max_span_length']} "
              f"intron={r['max_intron_length']}")
        ok &= p["consensus_settings"]["preset"] == "fungi" and s["long_read_models"] == p["consensus_gtf"]
sys.exit(0 if ok else 1)
PY

echo "--- 3. historical fungi preset, individual stages"
run_stages() {  # run_stages <name> <runner...>
  local name=$1; shift
  local d="$OUT/$name"; rm -rf "$d"; mkdir -p "$d"
  "$@" preflight --preset fungi "${EVID[@]}" --output-dir "$d/preflight" --no-log-file \
    > "$d/preflight.log" 2>&1
  "$@" build --preset fungi "${EVID[@]}" --gene-prefix ZT --output-dir "$d/build" \
    --validate-fasta --no-log-file > "$d/build.log" 2>&1
  "$@" finalise --build-dir "$d/build" --genome "$FX/genome.fa" --output-dir "$d/finalise" \
    > "$d/finalise.log" 2>&1
  check "$d"
}
gmb() { "$GMB_BIN/gmb-$1" "${@:2}"; }
run_stages fungi_stages gmb

if [ -n "$BASELINE" ]; then
  echo "--- historical fungi preset vs baseline checkout $(git -C "$BASELINE" rev-parse --short HEAD)"
  # Call each command's main() directly: older checkouts lack `python -m` guards.
  base() { (cd "$BASELINE/support_scripts/gmb" && "$PY" -c \
    'import importlib, sys; m = importlib.import_module("gmb.cli." + sys.argv[1]); sys.argv = sys.argv[1:]; sys.exit(m.main())' "$@"); }
  run_stages baseline_checkout base
  for f in finalise/consensus.gff3 finalise/canonical/consensus.canonical_annotated.gff3 \
           finalise/cdna.fa finalise/cds.fa finalise/prot.fa build/evidence_attribution.tsv; do
    cmp -s "$OUT/fungi_stages/$f" "$OUT/baseline_checkout/$f" \
      && echo "identical: $f" || { echo "DIFFERS: $f"; fail=1; }
  done
fi

[ $fail = 0 ] && echo "SMOKE TEST PASSED" || { echo "SMOKE TEST FAILED"; exit 1; }
