#!/usr/bin/env bash
# End-to-end smoke test of both fungal configurations on the bundled 500 kb Z. tritici fixture:
# preflight -> build -> finalise, independent output checks, a determinism re-run, and
# (optionally) byte-equivalence of configuration A with a baseline GMB checkout.
#
# Usage: smoke_fixture.sh <outdir> [<baseline ensembl-anno checkout>]
#   GMB_BIN   directory with gmb-* and python (default: from PATH)
# Exit 1 on any failure.
set -euo pipefail
OUT=$1; BASELINE=${2:-}
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
GMB_DIR="$(cd "$HERE/.." && pwd)"
GMB_BIN="${GMB_BIN:-$(dirname "$(command -v gmb-build)")}"
PY="$GMB_BIN/python"; [ -x "$PY" ] || PY="$(command -v python3)"
FX="$GMB_DIR/tests/fixtures/z_tritici_region1"
OVERLAY="$GMB_DIR/configs/fungi_evidence_integration.experimental.yaml"
EVID=(--genome "$FX/genome.fa" --helixer "$FX/helixer_remapped.gff3"
      --scallop "$FX/scallop_geneset.gtf" --stringtie "$FX/stringtie_geneset.gtf"
      --orthodb "$FX/orthodb_geneset.gtf" --uniprot "$FX/uniprot_geneset.gtf")
mkdir -p "$OUT"; fail=0

run_config() {  # run_config <name> <runner prefix...> -- [extra gmb args]
  local name=$1; shift
  local -a runner=(); while [ "$1" != "--" ]; do runner+=("$1"); shift; done; shift
  local d="$OUT/$name"; rm -rf "$d"; mkdir -p "$d"
  "${runner[@]}" preflight --preset fungi "$@" "${EVID[@]}" --output-dir "$d/preflight" \
    --no-log-file > "$d/preflight.log" 2>&1   # fixture: WARN only (exit 0); FAIL stops the test
  "${runner[@]}" build --preset fungi "$@" "${EVID[@]}" --gene-prefix ZT \
    --output-dir "$d/build" --validate-fasta --no-log-file > "$d/build.log" 2>&1
  "${runner[@]}" finalise --build-dir "$d/build" --genome "$FX/genome.fa" \
    --output-dir "$d/finalise" > "$d/finalise.log" 2>&1
  "$PY" "$HERE/check_handover.py" "$d/finalise" "$d/build" || fail=1
}
gmb() { "$GMB_BIN/gmb-$1" "${@:2}"; }

run_config A_baseline gmb --
run_config B_experimental gmb -- --config "$OVERLAY"
run_config B_experimental_rerun gmb -- --config "$OVERLAY"

echo "--- configuration provenance"
"$PY" - "$OUT" <<'PY' || fail=1
import sys, yaml, pathlib
out = pathlib.Path(sys.argv[1]); ok = True
for name, want in (("A_baseline", ("score", "off", 0, False)),
                   ("B_experimental", ("junction_supported", "drop", 150, True))):
    c = yaml.safe_load(open(out / name / "build" / "resolved_config.yaml"))
    got = (c["scoring"]["primary_selection"], c["scoring"]["isoform_cds_overlap"],
           c["scoring"]["min_cds_bp"], c["canonical_selection"]["prefer_build_primary"])
    print(f"{name}: primary_selection={got[0]} isoform_cds_overlap={got[1]} "
          f"min_cds_bp={got[2]} prefer_build_primary={got[3]}")
    ok &= got == want
    man = (out / name / "finalise" / "handover_manifest.json").read_text()
    ok &= "resolved_config" in man
sys.exit(0 if ok else 1)
PY

echo "--- determinism (B run twice)"
for f in finalise/consensus.gff3 finalise/canonical/consensus.canonical_annotated.gff3 \
         finalise/prot.fa build/evidence_attribution.tsv build/excluded_candidates.tsv; do
  cmp -s "$OUT/B_experimental/$f" "$OUT/B_experimental_rerun/$f" \
    && echo "identical: $f" || { echo "DIFFERS: $f"; fail=1; }
done

if [ -n "$BASELINE" ]; then
  echo "--- configuration A vs baseline checkout $(git -C "$BASELINE" rev-parse --short HEAD)"
  # Call each command's main() directly: before this release `python -m gmb.cli.finalise`
  # exited 0 without running.
  base() { (cd "$BASELINE/support_scripts/gmb" && "$PY" -c \
    'import importlib, sys; m = importlib.import_module("gmb.cli." + sys.argv[1]); sys.argv = sys.argv[1:]; sys.exit(m.main())' "$@"); }
  run_config baseline_checkout base --
  for f in finalise/consensus.gff3 finalise/canonical/consensus.canonical_annotated.gff3 \
           finalise/cdna.fa finalise/cds.fa finalise/prot.fa; do
    cmp -s "$OUT/A_baseline/$f" "$OUT/baseline_checkout/$f" \
      && echo "identical: $f" || { echo "DIFFERS: $f"; fail=1; }
  done
  # Attribution and canonical tables only gain trailing columns.
  "$PY" - "$OUT" <<'PY' || fail=1
import sys, pandas as pd
out = sys.argv[1]; ok = True
for rel in ("build/evidence_attribution.tsv", "finalise/canonical/canonical_transcripts.tsv"):
    a = pd.read_csv(f"{out}/A_baseline/{rel}", sep="\t")
    b = pd.read_csv(f"{out}/baseline_checkout/{rel}", sep="\t")
    same = list(a.columns[:len(b.columns)]) == list(b.columns) and a[b.columns].equals(b)
    print(f"{rel}: baseline columns {'identical' if same else 'DIFFER'}; "
          f"appended {list(a.columns[len(b.columns):])}")
    ok &= same
sys.exit(0 if ok else 1)
PY
fi

[ $fail = 0 ] && echo "SMOKE TEST PASSED" || { echo "SMOKE TEST FAILED"; exit 1; }
