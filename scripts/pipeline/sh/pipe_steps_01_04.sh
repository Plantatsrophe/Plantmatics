#!/usr/bin/env bash
# scripts/pipeline/sh/pipe_steps_01_04.sh - Stage runners for Steps 01 through 04.

_pipe_run_step01_harvest() {
if [[ "${START_STEP}" -le 1 && "${REUSE_CANONICAL_DATA}" == "false" ]]; then
  if [[ "${SKIP_HARVEST}" == "false" && ( -n "${FORCE_FLAG} $([[ "${RESEGMENT}" == "true" ]] && echo "--force")" || ! -s "${RUN_DIR}/data/tables/curated_vouchers.csv" ) ]]; then
    echo "[Step 01/07] Running Voucher Harvester (with high-res image downloads & quality audit)..."
    HARVEST_TAXON_ARGS=()
    if [[ ${#CUSTOM_TAXA_RAW[@]} -gt 0 ]]; then
      HARVEST_TAXON_ARGS=(--taxon "${CUSTOM_TAXA_RAW[@]}" "${EXPAND_SYNONYMS_FLAG}")
    elif [[ -n "${FOCAL_TAXON}" && "${INTRASPECIFIC_MODE}" == "true" ]]; then
      HARVEST_TAXON_ARGS=(--taxon "${FOCAL_TAXON}" "${EXPAND_SYNONYMS_FLAG}")
    else
      HARVEST_TAXON_ARGS=("${EXPAND_SYNONYMS_FLAG}")
    fi
    HARVEST_GEO_ARGS=()
    if [[ -n "${GEOGRAPHIC_SCOPE}" ]]; then
      HARVEST_GEO_ARGS+=(--geographic-scope "${GEOGRAPHIC_SCOPE}")
    fi
    if [[ "${INCLUDE_WESTERN}" == "true" ]]; then
      HARVEST_GEO_ARGS+=(--include-western)
    fi
    if [[ ${#COUNTRIES_RAW[@]} -gt 0 ]]; then
      HARVEST_GEO_ARGS+=(--countries "${COUNTRIES_RAW[@]}")
    fi
    if [[ "${HARVEST_FULL_CONTINENTAL_RANGE:-}" == "true" ]]; then
      HARVEST_GEO_ARGS+=(--harvest-full-continental-range)
    elif [[ "${HARVEST_FULL_CONTINENTAL_RANGE:-}" == "false" ]]; then
      HARVEST_GEO_ARGS+=(--no-harvest-full-continental-range)
    fi
    HARVEST_SPATIAL_ARGS=()
    if [[ -n "${SPATIAL_THIN_KM}" && "${SPATIAL_THIN_KM}" != "0" && "${SPATIAL_THIN_KM}" != "0.0" ]]; then
      HARVEST_SPATIAL_ARGS+=(--spatial-thin-km "${SPATIAL_THIN_KM}")
    fi
    if [[ -n "${MAX_PER_STATE}" ]]; then
      HARVEST_SPATIAL_ARGS+=(--max-per-state "${MAX_PER_STATE}")
    fi
    if [[ -n "${MAX_PER_REGION}" ]]; then
      HARVEST_SPATIAL_ARGS+=(--max-per-region "${MAX_PER_REGION}")
    fi
    if [[ "${USDA_COUNTY_FILTER}" == "true" ]]; then
      HARVEST_SPATIAL_ARGS+=(--usda-county-filter)
    elif [[ "${USDA_COUNTY_FILTER}" == "false" ]]; then
      HARVEST_SPATIAL_ARGS+=(--no-usda-county-filter)
    fi
    if [[ -n "${USDA_COUNTY_BUFFER_KM}" ]]; then
      HARVEST_SPATIAL_ARGS+=(--usda-county-buffer-km "${USDA_COUNTY_BUFFER_KM}")
    fi
    if [[ "${USDA_EXEMPT_GOLD}" == "true" ]]; then
      HARVEST_SPATIAL_ARGS+=(--usda-exempt-gold)
    elif [[ "${USDA_EXEMPT_GOLD}" == "false" ]]; then
      HARVEST_SPATIAL_ARGS+=(--no-usda-exempt-gold)
    fi
    if [[ -n "${THREADS}" ]]; then
      HARVEST_CONC=$(( THREADS > 32 ? 32 : (THREADS < 1 ? 1 : THREADS) ))
      HARVEST_SPATIAL_ARGS+=(--concurrency "${HARVEST_CONC}")
    fi
    "${PYTHON_MAIN}" scripts/data_prep/01_voucher_harvester.py \
      --out-dir "data/raw_vouchers" \
      --output-csv "${RUN_DIR}/data/tables/curated_vouchers.csv" \
      --output-plot "${RUN_DIR}/outputs/figures/fig01_cohort_biogeography_and_authority.png" \
      "${HARVEST_TAXON_ARGS[@]}" \
      "${HARVEST_GEO_ARGS[@]}" \
      "${HARVEST_SPATIAL_ARGS[@]}" \
      --download-images \
      ${FORCE_FLAG} $([[ "${RESEGMENT}" == "true" ]] && echo "--force")
    if [[ "${CUSTOM_GEO_SCOPE}" != "true" && -z "${FOCAL_TAXON}" && -z "${EXCLUDE_TAXON_CSV}" && -z "${DETERMINER_TIER}" && "${GOLD_ONLY}" != "true" && "${EXCLUDE_MISIDENTIFIED}" != "true" && ( -z "${SPATIAL_THIN_KM}" || "${SPATIAL_THIN_KM}" == "0" || "${SPATIAL_THIN_KM}" == "0.0" ) && -z "${MAX_PER_STATE}" && -z "${MAX_PER_REGION}" && -s "${RUN_DIR}/data/tables/curated_vouchers.csv" ]]; then
      mkdir -p "${PROJECT_ROOT}/data/tables"
      cp -f "${RUN_DIR}/data/tables/curated_vouchers.csv" "${PROJECT_ROOT}/data/tables/curated_vouchers.csv.tmp"
      mv -f "${PROJECT_ROOT}/data/tables/curated_vouchers.csv.tmp" "${PROJECT_ROOT}/data/tables/curated_vouchers.csv"
    fi
  fi
fi

# Ensure ${RUN_DIR}/data/tables/curated_vouchers.csv exists prior to Step 02+
if [[ ! -s "${RUN_DIR}/data/tables/curated_vouchers.csv" ]]; then
  if [[ -s "${SEED_DATA_DIR}/tables/curated_vouchers.csv" ]]; then
    echo "[Preflight] Copying ${SEED_DATA_DIR}/tables/curated_vouchers.csv -> ${RUN_DIR}/data/tables/curated_vouchers.csv"
    cp -f "${SEED_DATA_DIR}/tables/curated_vouchers.csv" "${RUN_DIR}/data/tables/curated_vouchers.csv.tmp"
    mv -f "${RUN_DIR}/data/tables/curated_vouchers.csv.tmp" "${RUN_DIR}/data/tables/curated_vouchers.csv"
  elif [[ -s "data/tables/curated_vouchers.csv" ]]; then
    echo "[Preflight] Copying shared data/tables/curated_vouchers.csv -> ${RUN_DIR}/data/tables/curated_vouchers.csv"
    cp -f "data/tables/curated_vouchers.csv" "${RUN_DIR}/data/tables/curated_vouchers.csv.tmp"
    mv -f "data/tables/curated_vouchers.csv.tmp" "${RUN_DIR}/data/tables/curated_vouchers.csv"
  else
    echo "[ERROR] Missing curated_vouchers.csv in ${RUN_DIR}/data/tables/, ${SEED_DATA_DIR}/tables/, and data/tables/." >&2
    exit 1
  fi
fi
}

_pipe_run_step02_segment() {
if [[ "${START_STEP}" -le 2 && "${REUSE_CANONICAL_DATA}" == "false" ]]; then
  echo "[Step 02/07] Running LeafMachine3 Segmentation & Organ Extraction..."
  find "${RUN_DIR}/data/cropped_patches" -type l -delete 2>/dev/null || true
  "${PYTHON_LM3}" scripts/pipeline/02_segment_and_extract.py \
    --vouchers "${RUN_DIR}/data/tables/curated_vouchers.csv" \
    --output-dir "${RUN_DIR}/data" \
    ${FORCE_FLAG} $([[ "${RESEGMENT}" == "true" ]] && echo "--force") \
    "${LIMIT_ARG[@]}" \
    "${COMPUTE_ARGS[@]}"

  # Sync newly extracted masks, contours, cropped patches, and Step 02 manifests back to canonical root data/
  echo "[Step 02/07] Syncing newly extracted Step 02 assets back to canonical root data/ store..."
  mkdir -p \
    "${PROJECT_ROOT}/data/masks" \
    "${PROJECT_ROOT}/data/contours" \
    "${PROJECT_ROOT}/data/ect" \
    "${PROJECT_ROOT}/data/cropped_patches/capitula" \
    "${PROJECT_ROOT}/data/tables"
  cp -alf "${RUN_DIR}/data/masks/." "${PROJECT_ROOT}/data/masks/" 2>/dev/null || cp -rf "${RUN_DIR}/data/masks/." "${PROJECT_ROOT}/data/masks/"
  cp -alf "${RUN_DIR}/data/contours/." "${PROJECT_ROOT}/data/contours/" 2>/dev/null || cp -rf "${RUN_DIR}/data/contours/." "${PROJECT_ROOT}/data/contours/"
  if [[ -d "${RUN_DIR}/data/ect" ]]; then
    cp -alf "${RUN_DIR}/data/ect/." "${PROJECT_ROOT}/data/ect/" 2>/dev/null || cp -rf "${RUN_DIR}/data/ect/." "${PROJECT_ROOT}/data/ect/" 2>/dev/null || true
  fi
  cp -alf "${RUN_DIR}/data/cropped_patches/." "${PROJECT_ROOT}/data/cropped_patches/" 2>/dev/null || cp -rf "${RUN_DIR}/data/cropped_patches/." "${PROJECT_ROOT}/data/cropped_patches/"
  if [[ "${CUSTOM_GEO_SCOPE}" != "true" && -z "${FOCAL_TAXON}" && -z "${EXCLUDE_TAXON_CSV}" && -z "${DETERMINER_TIER}" && "${GOLD_ONLY}" != "true" && "${EXCLUDE_MISIDENTIFIED}" != "true" && ${#LIMIT_ARG[@]} -eq 0 ]]; then
    "${PYTHON_MAIN}" -c "
import math
import os
from pathlib import Path
run_tables = Path('${RUN_DIR}/data/tables')
root_tables = Path('${PROJECT_ROOT}/data/tables')
for tbl in ('extracted_leaf_manifest.csv', 'failed_qc_vouchers.csv', 'voucher_reproductive_metrics.csv', 'voucher_phenological_states.csv'):
    src = run_tables / tbl
    if not src.exists():
        continue
    dst = root_tables / tbl
    txt = src.read_text(encoding='utf-8')
    if 'manifest' in tbl:
        txt = txt.replace('/${RUN_DIR}/data/', '/data/')
    tmp = dst.with_suffix(dst.suffix + '.tmp')
    with open(tmp, 'w', encoding='utf-8', newline='') as f:
        f.write(txt)
        f.flush()
        os.fsync(f.fileno())
    tmp.replace(dst)
"
  fi
  run_preflight_subset_filter
fi
}

_pipe_run_step03_fourier() {
if [[ "${START_STEP}" -le 3 ]]; then
  echo "[Step 03/07] Running Elliptic Fourier Extractor..."
  Rscript scripts/morphometrics/03_fourier_extractor.R \
    --input "${RUN_DIR}/data/contours" \
    --contours-dir "${RUN_DIR}/data/contours" \
    --manifest "${RUN_DIR}/data/tables/extracted_leaf_manifest.csv" \
    --vouchers "${RUN_DIR}/data/tables/curated_vouchers.csv" \
    --output "${RUN_DIR}/data/tables/leaf_efa_harmonics.csv" \
    --output-individual "${RUN_DIR}/data/tables/leaf_efa_harmonics_individual.csv" \
    --figures-dir "${RUN_DIR}/outputs/figures"
fi
}

_pipe_run_step04_gmm() {
if [[ "${START_STEP}" -le 4 ]]; then
  echo "[Step 04/07] Running GMM Clustering & Passive CDA..."
  Rscript scripts/morphometrics/04_gmm_morphotools.R \
    --vouchers "${RUN_DIR}/data/tables/curated_vouchers.csv" \
    --harmonics "${RUN_DIR}/data/tables/leaf_efa_harmonics.csv" \
    --reproductive "${RUN_DIR}/data/tables/voucher_reproductive_metrics.csv" \
    --allometry-report "${RUN_DIR}/outputs/reports/allometry_significance_audit.csv" \
    --allometry-pairwise "${RUN_DIR}/outputs/reports/allometry_pairwise_trajectories.csv" \
    --allometry-plot "${RUN_DIR}/outputs/figures/fig04a_foliar_allometry_audit.png" \
    --output-batch-audit "${RUN_DIR}/outputs/reports/institutional_batch_effect_audit.csv" \
    --output-flags "${RUN_DIR}/data/tables/morphometric_flags.csv" \
    --output-plot "${RUN_DIR}/outputs/figures/fig04b_gmm_cda_morphometrics.png" \
    --ecoregion-plot "${RUN_DIR}/outputs/figures/fig04c_ecoregion_and_batch_audit.png" \
    --output-report "${RUN_DIR}/outputs/reports/gmm_bayes_factors_summary.csv" \
    "${FOCAL_ARGS[@]}" \
    "${INTRASPECIFIC_ARGS[@]}"

  echo "[Step 04/07] Regenerating post-QC Step 02 publication figures..."
  "${PYTHON_MAIN}" scripts/pipeline/02_segment_and_extract.py \
    --vouchers "${RUN_DIR}/data/tables/curated_vouchers.csv" \
    --output-dir "${RUN_DIR}/data" \
    --regenerate-figures-only \
    "${COMPUTE_ARGS[@]}"
fi
}
