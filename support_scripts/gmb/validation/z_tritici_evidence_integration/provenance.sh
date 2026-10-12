#!/usr/bin/env bash
# Record code and tool versions for a run.
source "$(dirname "${BASH_SOURCE[0]}")/env.sh"
echo "date: $(date -u +%Y-%m-%dT%H:%M:%SZ)"
echo "gmb_commit: $(git -C "$GMB_DIR" rev-parse HEAD)$(git -C "$GMB_DIR" diff --quiet || echo ' (+uncommitted changes)')"
echo "gmb_branch: $(git -C "$GMB_DIR" rev-parse --abbrev-ref HEAD)"
"$PYTHON" - <<'PY'
import importlib.metadata as m, platform
print("python:", platform.python_version())
for p in ("gene-model-builder", "pandas", "pyranges", "numpy", "pyyaml"):
    try: print(f"{p}: {m.version(p)}")
    except m.PackageNotFoundError: print(f"{p}: not installed")
PY
AQC_PY="$(head -1 "$AQC" | sed 's/^#!//')"
echo "annotation-qc: $AQC"
[ -x "$AQC_PY" ] && "$AQC_PY" -c "import importlib.metadata as m; print('ensembl-genes:', m.version('ensembl-genes'))" 2>/dev/null || true
[ -n "${ENSEMBL_GENES_DIR:-}" ] && echo "ensembl_genes_commit: $(git -C "$ENSEMBL_GENES_DIR" rev-parse HEAD)"
LATEST_FILE="$("$PYTHON" -c "from importlib.resources import files; print(files('gmb.configs').joinpath('fungi_latest_settings.yaml'))")"
echo "fungi_latest_settings_sha256: $(shasum -a 256 "$LATEST_FILE" | cut -d' ' -f1)"
echo "windows_sha256: $(shasum -a 256 "$HERE/windows.tsv" | cut -d' ' -f1)"
