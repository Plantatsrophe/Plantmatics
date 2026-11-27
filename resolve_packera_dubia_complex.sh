#!/usr/bin/env bash
# ==============================================================================
# Script: resolve_packera_dubia_complex.sh
# Description:
#   Convenience launcher for the six-stage Packera dubia species-complex
#   discovery, congeneric misidentification purge, purified intraspecific
#   delimitation, unconfounded 15 km spatial/monographer sensitivity audits,
#   and supervised candidate-lineage diagnosability workflow inside a single
#   timestamped parent folder runs/<prefix>_YYYYMMDD_<hex6>/.
#   Forwards all CLI flags (--prefix, --parent-run-id, --resume, --force,
#   --threads, --gpu-batch-size, --dry-run, --skip-discovery, --skip-harvest,
#   --summarize-only, --focal-taxon) to scripts/pipeline/resolve_focal_species_complex.py.
# ==============================================================================

set -euo pipefail

PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export PYTHONPATH="${PROJECT_ROOT}${PYTHONPATH:+:${PYTHONPATH}}"
cd "${PROJECT_ROOT}"

if [[ -x "${PROJECT_ROOT}/.venv/bin/python" ]]; then
  PYTHON_BIN="${PROJECT_ROOT}/.venv/bin/python"
else
  PYTHON_BIN="python3"
fi

exec "${PYTHON_BIN}" "${PROJECT_ROOT}/scripts/pipeline/resolve_focal_species_complex.py" "$@"
