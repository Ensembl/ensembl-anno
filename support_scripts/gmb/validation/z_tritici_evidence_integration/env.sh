# Shared settings for the Z. tritici evidence-integration evaluation. Sourced by the
# other scripts; override any variable in the environment before running them.
#   WORK            (required) working directory for prepared inputs and runs
#   INPUT_DATA      directory with the GCA_000219625.1_* evidence files (prepare step)
#   REFERENCE_DATA  directory with Zymoseptoria_tritici_GCA_000219625.1.{fa,gff3}.bgz
#   GMB_BIN         directory holding gmb-build / gmb-finalise / gmb-longread-consensus
#   AQC             annotation-qc executable (ensembl-genes)
#   JOBS            parallel jobs (default 6)
#   PV_CONFIG       optional protein-validation overlay; enables the *_pv arms
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
GMB_DIR="$(cd "$HERE/../.." && pwd)"
: "${WORK:?set WORK to a working directory}"
GMB_BIN="${GMB_BIN:-$(dirname "$(command -v gmb-build)")}"
AQC="${AQC:-$(command -v annotation-qc || true)}"
PYTHON="${PYTHON:-$GMB_BIN/python}"
[ -x "$PYTHON" ] || PYTHON="$(command -v python3)"
JOBS="${JOBS:-6}"
# "overlay" arms (arms.tsv) run the latest fungal settings, the packaged preset
# fungi_latest_settings (= fungi + gmb/configs/fungi_latest_settings.yaml; formerly the
# overlay configs/fungi_evidence_integration.experimental.yaml, identical builds).
LATEST_PRESET=fungi_latest_settings
export WORK GMB_BIN AQC PYTHON GMB_DIR LATEST_PRESET HERE JOBS
