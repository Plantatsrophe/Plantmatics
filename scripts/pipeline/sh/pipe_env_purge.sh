#!/usr/bin/env bash
# scripts/pipeline/sh/pipe_env_purge.sh - Early purge routines, run env/seed setup, and metadata finalization.

_pipe_handle_early_purges() {
# Resolve Python interpreter early if --purge-runs-only was requested (preserves Step 01 & Step 02 segmentation)
if [[ "${PURGE_RUNS_ONLY}" == "true" && "${PURGE_SEGMENTS}" != "true" ]]; then
  if [[ -x "${PROJECT_ROOT}/.venv/bin/python" ]]; then
    _PURGE_PY="${PROJECT_ROOT}/.venv/bin/python"
  else
    _PURGE_PY="python3"
  fi
  echo "[Purge] Purging Step 03 to 07 run outputs and Tier B SQLite tables (preserving Step 01 & Step 02 segmentation)..."
  _PURGE_WS="${PLANTMATICS_WORKSPACE_ROOT:-${PACKERA_WORKSPACE_ROOT:-${PROJECT_ROOT}}}"
  "${_PURGE_PY}" -c "
import sys
from pathlib import Path
sys.path.insert(0, '${PROJECT_ROOT}')
from scripts.core.db import purge_runs_only
db_override = '${CUSTOM_DB_PATH}' or None
res = purge_runs_only(workspace_root=Path('${_PURGE_WS}'), db_path=db_override, refresh_env_buffers=True, vacuum=True)
rf = res['removed_files']
pr = res['purged_db_rows']
print(f'[Purge] Removed {rf} run/output files and {pr} Tier B SQLite rows (Step 01 vouchers & Step 02 segments preserved).')
"
  _NON_PURGE_RUN_ARGS=0
  _SKIP_NEXT_DB_VAL=false
  for _arg in "${ORIG_ARGS[@]}"; do
    if [[ "${_SKIP_NEXT_DB_VAL}" == "true" ]]; then
      _SKIP_NEXT_DB_VAL=false
      continue
    fi
    case "${_arg}" in
      --purge-runs-only|--purge-run-outputs|--db=*)
        ;;
      --db)
        _SKIP_NEXT_DB_VAL=true
        ;;
      *)
        _NON_PURGE_RUN_ARGS=$((_NON_PURGE_RUN_ARGS + 1))
        ;;
    esac
  done
  if [[ "${_NON_PURGE_RUN_ARGS}" -eq 0 ]]; then
    exit 0
  fi
fi

# Resolve Python interpreter early if --purge-segments was requested
if [[ "${PURGE_SEGMENTS}" == "true" ]]; then
  if [[ -x "${PROJECT_ROOT}/.venv/bin/python" ]]; then
    _PURGE_PY="${PROJECT_ROOT}/.venv/bin/python"
  else
    _PURGE_PY="python3"
  fi
  echo "[Purge] Purging all previous Step 02 to 07 segments, masks, contours, run directories, and SQLite segmentation tables..."
  _PURGE_WS="${PLANTMATICS_WORKSPACE_ROOT:-${PACKERA_WORKSPACE_ROOT:-${PROJECT_ROOT}}}"
  "${_PURGE_PY}" -c "
import sys
from pathlib import Path
sys.path.insert(0, '${PROJECT_ROOT}')
from scripts.core.db import purge_segments_and_runs
db_override = '${CUSTOM_DB_PATH}' or None
res = purge_segments_and_runs(workspace_root=Path('${_PURGE_WS}'), db_path=db_override, purge_runs=True, vacuum=True)
rf = res['removed_files']
pr = res['purged_db_rows']
print(f'[Purge] Removed {rf} artifact files and {pr} SQLite rows (Step 01 vouchers & raw images preserved).')
"
  # If invoked standalone (e.g. ./run_pipeline.sh --purge-segments [--db <path>]), exit cleanly after purge
  _NON_PURGE_ARGS=0
  for _arg in "${ORIG_ARGS[@]}"; do
    case "${_arg}" in
      --purge-segments|--purge-old-data|--purge-runs)
        ;;
      *)
        _NON_PURGE_ARGS=$((_NON_PURGE_ARGS + 1))
        ;;
    esac
  done
  if [[ "${_NON_PURGE_ARGS}" -eq 0 ]]; then
    exit 0
  fi
fi
}

_pipe_init_run_env_and_seed() {
if [[ -n "${CUSTOM_RUN_ID}" ]]; then
  RUN_ID="${CUSTOM_RUN_ID}"
else
  RUN_ID="${RUN_ID:-run_$(date +%Y%m%d_%H%M%S)}"
fi
RUN_ID="${RUN_ID//\\//}"
while [[ "${RUN_ID}" == ./* ]]; do
  RUN_ID="${RUN_ID#./}"
done
while [[ "${RUN_ID}" == /* ]]; do
  RUN_ID="${RUN_ID#/}"
done
RUN_ID="${RUN_ID#runs/}"
while [[ "${RUN_ID}" == */ ]]; do
  RUN_ID="${RUN_ID%/}"
done
if [[ -z "${RUN_ID}" ]]; then
  RUN_ID="run_$(date +%Y%m%d_%H%M%S)"
fi

RUN_DIR="runs/${RUN_ID}"
export RUN_ID
export PLANTMATICS_RUN_DIR="${PROJECT_ROOT}/${RUN_DIR}"
export PACKERA_RUN_DIR="${PLANTMATICS_RUN_DIR}"
if [[ ${#DEMO_FLAG[@]} -gt 0 ]]; then
  export PLANTMATICS_DEMO_MODE=1
  export PACKERA_DEMO_MODE=1
  export PLANTMATICS_DB_PATH="${CUSTOM_DB_PATH:-${PROJECT_ROOT}/${RUN_DIR}/data/demo_ephemeral.sqlite}"
  export PACKERA_DB_PATH="${PLANTMATICS_DB_PATH}"
else
  export PLANTMATICS_DB_PATH="${CUSTOM_DB_PATH:-${PROJECT_ROOT}/data/plantmatics.sqlite}"
  export PACKERA_DB_PATH="${PLANTMATICS_DB_PATH}"
fi

# 1. Create standardized run-scoped directory hierarchy
mkdir -p \
  "$(dirname "${RUN_DIR}")" \
  "${RUN_DIR}/data/tables" \
  "${RUN_DIR}/outputs/figures" \
  "${RUN_DIR}/outputs/reports" \
  "${RUN_DIR}/logs"

# Maintain runs/latest symlink pointing to the most recent RUN_DIR (unless disabled for parallel workers)
mkdir -p "runs"
if [[ "${PLANTMATICS_NO_LATEST_SYMLINK:-0}" != "1" ]]; then
  ln -sfn "${RUN_ID}" "runs/latest"
fi

export PYTHONUNBUFFERED=1
LOG_FILE="${RUN_DIR}/logs/pipeline.log"
exec > >(tee -a "${LOG_FILE}") 2>&1

START_TIME="$(date -u +"%Y-%m-%dT%H:%M:%SZ")"
GIT_SHA="$(git rev-parse HEAD 2>/dev/null || echo 'unknown')"

echo "=============================================================================="
echo "Plantmatics Pipeline: Run-Scoped Execution"
echo "  Run ID          : ${RUN_ID}"
echo "  Run Directory   : ${RUN_DIR}"
echo "  SQLite Database : ${PLANTMATICS_DB_PATH}"
echo "  Start Timestamp : ${START_TIME}"
echo "  Git Commit SHA  : ${GIT_SHA}"
echo "=============================================================================="

# Resolve Python interpreters (.venv for primary orchestration, LeafMachine3/.venv or .venv_LM3 for LeafMachine3)
if [[ -x "${PROJECT_ROOT}/.venv/bin/python" ]]; then
  PYTHON_MAIN="${PROJECT_ROOT}/.venv/bin/python"
else
  PYTHON_MAIN="python3"
fi

if [[ -x "${PROJECT_ROOT}/LeafMachine3/.venv/bin/python" ]]; then
  PYTHON_LM3="${PROJECT_ROOT}/LeafMachine3/.venv/bin/python"
elif [[ -x "${PROJECT_ROOT}/.venv_LM3/bin/python" ]]; then
  PYTHON_LM3="${PROJECT_ROOT}/.venv_LM3/bin/python"
else
  PYTHON_LM3="${PYTHON_MAIN}"
fi

# Ensure bundled nvidia-* CUDA/cuDNN shared library directories from PYTHON_LM3 are on LD_LIBRARY_PATH for onnxruntime-gpu
LM3_NVIDIA_LD_DIRS="$("${PYTHON_LM3}" -c "
from pathlib import Path
try:
    import nvidia
    roots = [Path(p) for p in getattr(nvidia, '__path__', ()) or ()]
    print(':'.join(sorted({str(d) for root in roots for d in root.glob('*/lib') if d.is_dir()})))
except Exception:
    pass
" 2>/dev/null || true)"
if [[ -n "${LM3_NVIDIA_LD_DIRS}" ]]; then
  export LD_LIBRARY_PATH="${LM3_NVIDIA_LD_DIRS}${LD_LIBRARY_PATH:+:${LD_LIBRARY_PATH}}"
fi

# Determine whether a non-default geographic scope or custom taxon harvest was requested
CUSTOM_GEO_SCOPE=false
if [[ ( -n "${GEOGRAPHIC_SCOPE}" && "${GEOGRAPHIC_SCOPE}" != "eastern_na" ) || "${INCLUDE_WESTERN}" == "true" || ${#COUNTRIES_RAW[@]} -gt 0 || ${#CUSTOM_TAXA_RAW[@]} -gt 0 ]]; then
  CUSTOM_GEO_SCOPE=true
fi
if [[ "${CUSTOM_GEO_SCOPE}" == "false" && -n "${FOCAL_TAXON}" && "${SKIP_HARVEST}" == "false" && "${START_STEP}" -le 1 && -s "${PROJECT_ROOT}/data/tables/curated_vouchers.csv" ]]; then
  if ! "${PYTHON_MAIN}" -c "
import sys
import pandas as pd
from scripts.core.harvester import standardize_packera_taxon
csv_p, focal_raw = sys.argv[1], sys.argv[2]
target = standardize_packera_taxon(focal_raw)
if not target:
    sys.exit(1)
df = pd.read_csv(csv_p, usecols=lambda c: c in ('scientificName', 'species_raw'), low_memory=False)
present = set()
for col in ('scientificName', 'species_raw'):
    if col in df.columns:
        present.update(df[col].dropna().map(standardize_packera_taxon).dropna())
sys.exit(0 if target in present else 1)
" "${PROJECT_ROOT}/data/tables/curated_vouchers.csv" "${FOCAL_TAXON}" 2>/dev/null; then
    CUSTOM_GEO_SCOPE=true
  fi
fi

# Resolve upstream Step 01-02 seed store (explicit --triage-source run data/ when deriving a filtered run, else canonical root data/)
SEED_DATA_DIR="${PROJECT_ROOT}/data"
HAS_EXPLICIT_SOURCE_RUN=false
if [[ "${EXCLUDE_MISIDENTIFIED}" == "true" && -n "${TRIAGE_SOURCE}" && "${TRIAGE_SOURCE^^}" != "ALL" && "${TRIAGE_SOURCE}" != "*" && "${TRIAGE_SOURCE}" != *,* ]]; then
  RUN_DATA_REAL="$(mkdir -p "${RUN_DIR}/data" && cd "${RUN_DIR}/data" && pwd -P)"
  for cand_src in "${TRIAGE_SOURCE}/data" "${PROJECT_ROOT}/runs/${TRIAGE_SOURCE}/data" "${PROJECT_ROOT}/${TRIAGE_SOURCE}/data"; do
    if [[ -d "${cand_src}" && -s "${cand_src}/tables/curated_vouchers.csv" ]]; then
      CAND_REAL="$(cd "${cand_src}" && pwd -P)"
      if [[ "${CAND_REAL}" != "${RUN_DATA_REAL}" ]]; then
        SEED_DATA_DIR="${CAND_REAL}"
        HAS_EXPLICIT_SOURCE_RUN=true
        break
      fi
    fi
  done
fi

_dir_has_files() {
  [[ -d "$1" ]] && [[ -n "$(find "$1" -maxdepth 2 -type f -print -quit 2>/dev/null)" ]]
}

# Determine whether to reuse upstream/canonical Step 01-02 assets
REUSE_CANONICAL_DATA=false
HAS_UPSTREAM_VOUCHERS=false
if [[ -s "${RUN_DIR}/data/tables/curated_vouchers.csv" || -s "${SEED_DATA_DIR}/tables/curated_vouchers.csv" || -s "${PROJECT_ROOT}/data/tables/curated_vouchers.csv" ]]; then
  HAS_UPSTREAM_VOUCHERS=true
fi
HAS_UPSTREAM_STEP02=false
if [[ -s "${SEED_DATA_DIR}/tables/extracted_leaf_manifest.csv" || -s "${PROJECT_ROOT}/data/tables/extracted_leaf_manifest.csv" ]]; then
  HAS_UPSTREAM_STEP02=true
elif [[ "${START_STEP}" -ge 3 ]] && ( _dir_has_files "${SEED_DATA_DIR}/contours" || _dir_has_files "${SEED_DATA_DIR}/cropped_patches/capitula" || _dir_has_files "${PROJECT_ROOT}/data/contours" || _dir_has_files "${PROJECT_ROOT}/data/cropped_patches/capitula" ); then
  HAS_UPSTREAM_STEP02=true
fi

if [[ "${HAS_UPSTREAM_VOUCHERS}" == "true" && "${HAS_UPSTREAM_STEP02}" == "true" ]]; then
  if [[ "${RESEGMENT}" == "false" && ${#LIMIT_ARG[@]} -eq 0 ]]; then
    if [[ "${START_STEP}" -ge 3 ]] || [[ -z "${FORCE_FLAG}" && ( "${CUSTOM_GEO_SCOPE}" == "false" || "${SKIP_HARVEST}" == "true" || "${HAS_EXPLICIT_SOURCE_RUN}" == "true" || "${GEOGRAPHIC_SCOPE:-eastern_na}" == "eastern_na" ) ]]; then
      REUSE_CANONICAL_DATA=true
    fi
  fi
fi

if [[ "${REUSE_CANONICAL_DATA}" == "true" ]]; then
  echo "[Preflight] Seeding Step 01 to 02 tables, masks, and contours from ${SEED_DATA_DIR}..."
  OVERWRITE_SEED_TABLES=false
  if [[ "${START_STEP}" -le 2 ]] || [[ "${START_STEP}" -eq 3 && "${HAS_EXPLICIT_SOURCE_RUN}" == "true" ]]; then
    OVERWRITE_SEED_TABLES=true
  fi
  for tbl in curated_vouchers.csv extracted_leaf_manifest.csv failed_qc_vouchers.csv voucher_reproductive_metrics.csv voucher_phenological_states.csv; do
    if [[ "${OVERWRITE_SEED_TABLES}" == "true" || ! -s "${RUN_DIR}/data/tables/${tbl}" ]]; then
      SRC_TBL=""
      if [[ -s "${SEED_DATA_DIR}/tables/${tbl}" ]]; then
        SRC_TBL="${SEED_DATA_DIR}/tables/${tbl}"
      elif [[ -s "${PROJECT_ROOT}/data/tables/${tbl}" ]]; then
        SRC_TBL="${PROJECT_ROOT}/data/tables/${tbl}"
      fi
      if [[ -n "${SRC_TBL}" ]]; then
        cp -f "${SRC_TBL}" "${RUN_DIR}/data/tables/${tbl}.tmp"
        mv -f "${RUN_DIR}/data/tables/${tbl}.tmp" "${RUN_DIR}/data/tables/${tbl}"
      fi
    fi
  done
  for asset_sub in masks contours; do
    SRC_ASSET="${PROJECT_ROOT}/data/${asset_sub}"
    if _dir_has_files "${SEED_DATA_DIR}/${asset_sub}"; then
      SRC_ASSET="${SEED_DATA_DIR}/${asset_sub}"
    fi
    DST_ASSET="${RUN_DIR}/data/${asset_sub}"
    if [[ -L "${DST_ASSET}" ]] || ! _dir_has_files "${DST_ASSET}"; then
      rm -rf "${DST_ASSET}"
      ln -sfnr "${SRC_ASSET}" "${DST_ASSET}"
    fi
  done
  if [[ ${#DEMO_FLAG[@]} -eq 0 ]]; then
    "${PYTHON_MAIN}" -c "
from pathlib import Path
from scripts.core.db import get_db_connection, ingest_existing_csvs_to_db
ws = Path('${PROJECT_ROOT}').resolve()
conn = get_db_connection(init_schema=True, workspace_root=ws)
try:
    v_n = conn.execute('SELECT COUNT(*) FROM vouchers;').fetchone()[0]
    el_n = conn.execute('SELECT COUNT(*) FROM extracted_leaves;').fetchone()[0]
finally:
    conn.close()
if v_n == 0 or el_n == 0:
    counts = ingest_existing_csvs_to_db(ws)
    print(f'[Preflight: SQLite Bootstrap] Synchronized canonical assets into SQLite: {counts}')
"
  fi
else
  for asset_sub in masks contours cropped_patches; do
    if [[ -L "${RUN_DIR}/data/${asset_sub}" ]]; then
      rm -f "${RUN_DIR}/data/${asset_sub}"
    elif [[ "${RESEGMENT}" == "true" && -d "${RUN_DIR}/data/${asset_sub}" ]]; then
      rm -rf "${RUN_DIR}/data/${asset_sub}"
    fi
  done
  mkdir -p \
    "${RUN_DIR}/data/masks" \
    "${RUN_DIR}/data/contours" \
    "${RUN_DIR}/data/cropped_patches/capitula"
fi

# Initialize run_metadata.json at start of run
"${PYTHON_MAIN}" -c "
import sys
from pathlib import Path
from scripts.core.config import write_run_metadata
write_run_metadata(
    run_dir=Path('${RUN_DIR}').resolve(),
    run_id='${RUN_ID}',
    start_time='${START_TIME}',
    cli_args=sys.argv[1:],
)
" -- "${ORIG_ARGS[@]}"
}

_pipe_finalize_run_metadata() {
# Finalize run_metadata.json and SQLite pipeline_runs with end timestamp and table row counts
END_TIME="$(date -u +"%Y-%m-%dT%H:%M:%SZ")"
"${PYTHON_MAIN}" -c "
import os
import sys
from pathlib import Path
from scripts.core.config import write_run_metadata
from scripts.core.db import ensure_pipeline_run, get_db_connection
rdir = Path('${RUN_DIR}').resolve()
write_run_metadata(
    run_dir=rdir,
    run_id='${RUN_ID}',
    start_time='${START_TIME}',
    end_time='${END_TIME}',
    cli_args=sys.argv[1:],
    status='completed',
    exit_code=0,
)
if os.environ.get('PLANTMATICS_DEMO_MODE') != '1' and os.environ.get('PACKERA_DEMO_MODE') != '1':
    conn = get_db_connection(init_schema=True, workspace_root=Path('${PROJECT_ROOT}').resolve())
    try:
        ensure_pipeline_run(
            conn,
            run_id='${RUN_ID}',
            run_dir=rdir,
            start_timestamp='${START_TIME}',
            end_timestamp='${END_TIME}',
            git_commit_sha='${GIT_SHA}',
            cli_arguments=sys.argv[1:],
        )
    finally:
        conn.close()
" -- "${ORIG_ARGS[@]}"

echo "=============================================================================="
echo "Pipeline Run Complete: ${RUN_DIR}"
echo "  Metadata JSON : ${RUN_DIR}/outputs/reports/run_metadata.json"
echo "  Figures Dir   : ${RUN_DIR}/outputs/figures/"
echo "  Tables Dir    : ${RUN_DIR}/data/tables/"
echo "  Latest Link   : runs/latest -> ${RUN_ID}"
echo "=============================================================================="
}
