#!/usr/bin/env bash
# ==============================================================================
# Script: run_pipeline.sh
# Project: Plantmatics - Packera Species Delimitation & Morphometrics Pipeline
# Affiliation: University of North Carolina at Chapel Hill Herbarium (NCU)
# Description: End-to-end run-scoped orchestrator (Steps 01 to 07).
# ==============================================================================

set -euo pipefail

PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export PYTHONPATH="${PROJECT_ROOT}${PYTHONPATH:+:${PYTHONPATH}}"
cd "${PROJECT_ROOT}"

ORIG_ARGS=("$@")
source "${PROJECT_ROOT}/scripts/pipeline/sh/pipe_cli_args.sh"
source "${PROJECT_ROOT}/scripts/pipeline/sh/pipe_env_purge.sh"
source "${PROJECT_ROOT}/scripts/pipeline/sh/pipe_steps_01_04.sh"

_pipe_init_defaults
while [[ $# -gt 0 ]]; do
  case "$1" in
    -h|--help) print_usage; exit 0 ;;
    -j|--threads|--workers) THREADS="$2"; shift 2 ;;
    --gpu-batch-size|--batch-size) GPU_BATCH_SIZE="$2"; shift 2 ;;
    --purge-runs-only|--purge-run-outputs) PURGE_RUNS_ONLY=true; shift ;;
    --gold-only) GOLD_ONLY=true; DETERMINER_TIER="Tier_1_Gold"; shift ;;
    --determiner-tier) DETERMINER_TIER="$2"; shift 2 ;;
    --determiner-tier=*) DETERMINER_TIER="${1#*=}"; shift ;;
    --exclude-misidentified) EXCLUDE_MISIDENTIFIED=true; shift ;;
    --exclude-misidentified=*) EXCLUDE_MISIDENTIFIED="${1#*=}"; shift ;;
    --triage-source) _pipe_append_triage_source "$2"; shift 2 ;;
    --triage-source=*) _pipe_append_triage_source "${1#*=}"; shift ;;
    --geographic-scope) GEOGRAPHIC_SCOPE="$2"; shift 2 ;;
    --include-western|--no-exclude-western) INCLUDE_WESTERN=true; shift ;;
    --harvest-full-continental-range) HARVEST_FULL_CONTINENTAL_RANGE="true"; shift ;;
    --no-harvest-full-continental-range) HARVEST_FULL_CONTINENTAL_RANGE="false"; shift ;;
    --skip-env-download) SKIP_ENV_DOWNLOAD=true; shift ;;
    --force-env-download) FORCE_ENV_DOWNLOAD=true; shift ;;
    --skip-nwi-network) SKIP_NWI_NETWORK=true; export PLANTMATICS_SKIP_NWI_NETWORK=1 PACKERA_SKIP_NWI_NETWORK=1; shift ;;
    --skip-cytology-download) SKIP_CYTOLOGY_DOWNLOAD=true; shift ;;
    --force-cytology-download) FORCE_CYTOLOGY_DOWNLOAD=true; shift ;;
    --download-public-cytology) DOWNLOAD_PUBLIC_CYTOLOGY=1; shift ;;
    *)
      _pipe_consume_extra_arg "$@"
      shift "${_PIPE_SHIFT}"
      ;;
  esac
done
_pipe_finalize_cli_args
_pipe_handle_early_purges
_pipe_init_run_env_and_seed

if [[ -n "${LINEAGE_DIAG_SOURCE}" ]]; then
  echo "[Stage 06/06] Running Supervised Candidate-Lineage Diagnosability & Key-Trait Resolution..."
  Rscript scripts/analysis/step07_lineage_diagnosability.R \
    --source-run "${LINEAGE_DIAG_SOURCE}" \
    --run-dir "${RUN_DIR}" \
    --run-id "${RUN_ID}" \
    ${FOCAL_ARGS[@]+"${FOCAL_ARGS[@]}"} \
    ${DEMO_FLAG[@]+"${DEMO_FLAG[@]}"}
  _pipe_finalize_run_metadata
  exit 0
fi

_pipe_run_step01_harvest

# Preflight subset filter: select_cascading_outgroup_anchors and Misidentification_Reassign pruning
run_preflight_subset_filter() {
  local _root_data="${PROJECT_ROOT}/data"
  source "${PROJECT_ROOT}/scripts/pipeline/sh/pipe_preflight_cohort_a.sh"
  _pipe_run_preflight_subset_filter "${_root_data}"
}

run_preflight_subset_filter

run_preflight_env_layers() {
  if [[ "${#DEMO_FLAG[@]}" -gt 0 || "${SKIP_ENV_DOWNLOAD}" == "true" || "${PLANTMATICS_DEMO_MODE:-0}" == "1" || "${PACKERA_DEMO_MODE:-0}" == "1" || -n "${PYTEST_CURRENT_TEST:-}" ]]; then
    return 0
  fi
  local env_extra_args=()
  if [[ -s "${RUN_DIR}/data/tables/curated_vouchers.csv" ]]; then
    env_extra_args+=("--vouchers" "${RUN_DIR}/data/tables/curated_vouchers.csv")
  fi
  if [[ -n "${CUSTOM_DB_PATH}" ]]; then
    env_extra_args+=("--db" "${CUSTOM_DB_PATH}")
  fi
  if [[ "${FORCE_ENV_DOWNLOAD}" == "true" ]]; then
    env_extra_args+=("--force")
  fi
  "${PYTHON_MAIN}" scripts/data_prep/download_environmental_layers.py \
    --env-dir "${PROJECT_ROOT}/data/environmental" \
    "${env_extra_args[@]}"
  local geochem_vouchers=""
  if [[ -s "${RUN_DIR}/data/tables/curated_vouchers.csv" ]]; then
    geochem_vouchers="${RUN_DIR}/data/tables/curated_vouchers.csv"
  elif [[ -s "${PROJECT_ROOT}/data/tables/curated_vouchers.csv" ]]; then
    geochem_vouchers="${PROJECT_ROOT}/data/tables/curated_vouchers.csv"
  fi
  if [[ -n "${geochem_vouchers}" ]]; then
    local geochem_args=("--env-dir" "${PROJECT_ROOT}/data/environmental" "--vouchers" "${geochem_vouchers}")
    if [[ -n "${CUSTOM_DB_PATH}" ]]; then
      geochem_args+=("--db" "${CUSTOM_DB_PATH}")
    fi
    if [[ "${FORCE_ENV_DOWNLOAD}" == "true" ]]; then
      geochem_args+=("--force")
    fi
    if [[ "${SKIP_NWI_NETWORK:-false}" == "true" || "${PLANTMATICS_SKIP_NWI_NETWORK:-0}" == "1" || "${PACKERA_SKIP_NWI_NETWORK:-0}" == "1" ]]; then
      geochem_args+=("--skip-nwi-network")
    fi
    "${PYTHON_MAIN}" scripts/data_prep/download_geochem_wetland_layers.py \
      "${geochem_args[@]}" || echo "[Preflight] Geochem/wetland layer provisioning warning (continuing)."
  fi
}

run_preflight_env_layers

run_preflight_public_cytology() {
  if [[ "${#DEMO_FLAG[@]}" -gt 0 || "${SKIP_CYTOLOGY_DOWNLOAD}" == "true" ]]; then
    return 0
  fi
  local cyto_dir="${PLANTMATICS_CYTOLOGY_DIR:-${PACKERA_CYTOLOGY_DIR:-${PROJECT_ROOT}/data/cytology}}"
  local cyto_rec_csv="${cyto_dir}/public_packera_cytology_records.csv"
  local cyto_extra_args=()
  if [[ -n "${CUSTOM_DB_PATH}" ]]; then
    cyto_extra_args+=("--db" "${CUSTOM_DB_PATH}")
  fi
  if [[ "${FORCE_CYTOLOGY_DOWNLOAD}" == "true" || "${DOWNLOAD_PUBLIC_CYTOLOGY}" -eq 1 ]]; then
    echo "[Preflight] Downloading public Packera cytological records into ${cyto_dir}..."
    "${PYTHON_MAIN}" scripts/data_prep/download_public_cytology.py \
      --cytology-dir "${cyto_dir}" \
      --force \
      "${cyto_extra_args[@]}" || echo "[Preflight] Public cytology download skipped cleanly after network warning."
  elif [[ ! -s "${cyto_rec_csv}" ]]; then
    echo "[Preflight] Staging initial public Packera cytological records into ${cyto_dir}..."
    "${PYTHON_MAIN}" scripts/data_prep/download_public_cytology.py \
      --cytology-dir "${cyto_dir}" \
      "${cyto_extra_args[@]}" || echo "[Preflight] Initial public cytology download skipped cleanly after network warning."
  else
    "${PYTHON_MAIN}" scripts/data_prep/download_public_cytology.py \
      --cytology-dir "${cyto_dir}" \
      --sync-db-only \
      "${cyto_extra_args[@]}" || echo "[Preflight] Public cytology SQLite sync warning (continuing)."
  fi
}

run_preflight_public_cytology

run_preflight_taxonomic_backbone() {
  if [[ "${#DEMO_FLAG[@]}" -gt 0 || "${SKIP_BACKBONE_DOWNLOAD}" == "true" || "${PLANTMATICS_DEMO_MODE:-0}" == "1" || "${PACKERA_DEMO_MODE:-0}" == "1" || -n "${PYTEST_CURRENT_TEST:-}" ]]; then
    return 0
  fi
  local bb_script="${PROJECT_ROOT}/scripts/data_prep/download_taxonomic_backbone.py"
  if [[ ! -f "${bb_script}" ]]; then
    return 0
  fi
  local tax_dir="${PLANTMATICS_TAXONOMY_DIR:-${PACKERA_TAXONOMY_DIR:-${PROJECT_ROOT}/data/taxonomy}}"
  local bb_args=()
  if [[ -n "${CUSTOM_DB_PATH}" ]] && grep -q -- '--db' "${bb_script}"; then
    bb_args+=("--db" "${CUSTOM_DB_PATH}")
  fi
  local has_tax_cache=0
  if [[ -d "${tax_dir}" ]] && find "${tax_dir}" -maxdepth 1 -type f \( -name '*.csv' -o -name '*.json' \) -print -quit | grep -q .; then
    has_tax_cache=1
  fi
  if [[ "${FORCE_BACKBONE_DOWNLOAD}" == "true" || "${DOWNLOAD_TAXONOMIC_BACKBONE}" -eq 1 ]]; then
    echo "[Preflight] Downloading dynamic taxonomic backbone profiles into ${tax_dir}..."
    if grep -q -- '--force' "${bb_script}"; then
      bb_args+=("--force")
    fi
    "${PYTHON_MAIN}" "${bb_script}" "${bb_args[@]}" || echo "[Preflight] Taxonomic backbone download skipped cleanly after network warning."
  elif [[ "${has_tax_cache}" -eq 0 ]]; then
    echo "[Preflight] Staging initial taxonomic backbone profiles into ${tax_dir}..."
    "${PYTHON_MAIN}" "${bb_script}" "${bb_args[@]}" || echo "[Preflight] Initial taxonomic backbone download skipped cleanly after network warning."
  else
    if grep -q -- '--sync-db-only' "${bb_script}"; then
      bb_args+=("--sync-db-only")
    elif grep -q -- '--sync-only' "${bb_script}"; then
      bb_args+=("--sync-only")
    else
      return 0
    fi
    "${PYTHON_MAIN}" "${bb_script}" "${bb_args[@]}" || echo "[Preflight] Taxonomic backbone SQLite sync warning (continuing)."
  fi
}

run_preflight_taxonomic_backbone

# Steps 02 to 04: 02_segment_and_extract.py, 03_fourier_extractor.R, 04_gmm_morphotools.R
_pipe_run_step02_segment
_pipe_run_step03_fourier
_pipe_run_step04_gmm

# ------------------------------------------------------------------------------
# Step 05: Multimodal Spatial Random Forests, SoilGrids/SSURGO & Hopkins Phenology
# ------------------------------------------------------------------------------
if [[ "${START_STEP}" -le 5 ]]; then
  echo "[Step 05/07] Running Multimodal Spatial RF & Niche Identity Tests..."
  Rscript scripts/analysis/05_multimodal_spatial_rf.R \
    --vouchers "${RUN_DIR}/data/tables/curated_vouchers.csv" \
    --pheno-states "${RUN_DIR}/data/tables/voucher_phenological_states.csv" \
    --env-dir "data/environmental" \
    --ssurgo-table "${RUN_DIR}/data/tables/ssurgo_edaphic_validation.csv" \
    --contingency-table "${RUN_DIR}/outputs/reports/micro_edaphic_outcrop_contingency.csv" \
    --pheno-anomalies "${RUN_DIR}/data/tables/phenological_anomalies.csv" \
    --pheno-summary "${RUN_DIR}/outputs/reports/phenological_anomaly_summary.csv" \
    --pheno-plot "${RUN_DIR}/outputs/figures/fig05b_hopkins_circular_phenology.png" \
    --output-flags "${RUN_DIR}/data/tables/multimodal_conflict_flags.csv" \
    --output-summary "${RUN_DIR}/outputs/reports/multimodal_spatial_rf_summary.csv" \
    --output-plot "${RUN_DIR}/outputs/figures/fig05a_ecological_niche_and_rf_importance.png" \
    --sdm-plot "${RUN_DIR}/outputs/figures/fig05c_sdm_suitability_and_distribution_map.png" \
    "${FOCAL_ARGS[@]}" \
    "${INTRASPECIFIC_ARGS[@]}" \
    "${DEMO_FLAG[@]}"
fi

# ------------------------------------------------------------------------------
# Step 06: Independent Genomic Ingestion & Phase 2 Phylogenomics Stream
# ------------------------------------------------------------------------------
GENOMIC_STEP07_ARGS=()
if [[ -s "${RUN_DIR}/outputs/genomics/genomic_assignments.csv" ]]; then
  GENOMIC_STEP07_ARGS=(--genomic-dir "${RUN_DIR}/outputs/genomics")
fi
GENOMICS_DIR_EFF="${PLANTMATICS_GENOMICS_DIR:-${PACKERA_GENOMICS_DIR:-${PROJECT_ROOT}/data/genomics}}"
if [[ "${START_STEP}" -le 6 ]]; then
  if [[ "${DOWNLOAD_PUBLIC_GENOMICS}" -eq 1 ]]; then
    echo "[Step 06/07] Downloading public Packera sequences and phylogenomic trees into ${GENOMICS_DIR_EFF}..."
    "${PYTHON_MAIN}" "${PROJECT_ROOT}/main.py" download-public-genomics --genomics-dir "${GENOMICS_DIR_EFF}"
  fi
  HAS_GENO=0
  if [[ -d "${GENOMICS_DIR_EFF}" ]]; then
    GENOMIC_CHECK_OUT="$("${PYTHON_MAIN}" scripts/data_prep/06_ingest_genetic_data.py --genomics-dir "${GENOMICS_DIR_EFF}" --check-only 2>/dev/null || true)"
    if [[ "${GENOMIC_CHECK_OUT}" == *"HAS_GENOMIC_INPUTS=1"* ]]; then
      HAS_GENO=1
    else
      DYN_V_ARGS=()
      if [[ -s "${RUN_DIR}/data/tables/curated_vouchers.csv" ]]; then
        DYN_V_ARGS=(--vouchers "${RUN_DIR}/data/tables/curated_vouchers.csv")
      fi
      GENOMIC_DYN_OUT="$("${PYTHON_MAIN}" scripts/data_prep/06_ingest_genetic_data.py --genomics-dir "${GENOMICS_DIR_EFF}" ${DYN_V_ARGS[@]+"${DYN_V_ARGS[@]}"} ${FOCAL_ARGS[@]+"${FOCAL_ARGS[@]}"} --check-dynamic 2>/dev/null || true)"
      if [[ "${GENOMIC_DYN_OUT}" == *"HAS_DYNAMIC_GENOMIC_INPUTS=1"* ]]; then
        HAS_GENO=1
      fi
    fi
  elif [[ -n "${FOCAL_TAXON}" ]]; then
    HAS_GENO=1
  fi
  if [[ "${HAS_GENO}" -eq 1 ]]; then
    echo "[Step 06/07] Running Independent Genomic Ingestion & Phylogenomic Integration..."
    GENO_IN_DB_ARGS=()
    GENO_R_DB_ARGS=()
    if [[ -n "${PLANTMATICS_DB_PATH:-}" ]]; then
      GENO_IN_DB_ARGS=(--db-path "${PLANTMATICS_DB_PATH}")
      GENO_R_DB_ARGS=(--db "${PLANTMATICS_DB_PATH}")
    fi
    "${PYTHON_MAIN}" scripts/data_prep/06_ingest_genetic_data.py \
      --genomics-dir "${GENOMICS_DIR_EFF}" \
      --out-dir "${RUN_DIR}/outputs/reports" \
      --vouchers "${RUN_DIR}/data/tables/curated_vouchers.csv" \
      --genomic-tier auto \
      ${FOCAL_ARGS[@]+"${FOCAL_ARGS[@]}"} \
      ${INTRASPECIFIC_ARGS[@]+"${INTRASPECIFIC_ARGS[@]}"} \
      ${GENO_IN_DB_ARGS[@]+"${GENO_IN_DB_ARGS[@]}"} \
      ${DEMO_FLAG[@]+"${DEMO_FLAG[@]}"}
    Rscript scripts/analysis/06_genomic_phylogenomics.R \
      --genomics-dir "${GENOMICS_DIR_EFF}" \
      --parsed-dir "${RUN_DIR}/outputs/reports" \
      --out-dir "${RUN_DIR}/outputs/genomics" \
      --tables-dir "${RUN_DIR}/data/tables" \
      --vouchers "${RUN_DIR}/data/tables/curated_vouchers.csv" \
      --genomic-tier auto \
      --figure-a "${RUN_DIR}/outputs/figures/fig06a_genomic_pcoa_tree_concordance.png" \
      --figure-b "${RUN_DIR}/outputs/figures/fig06b_reticulation_gdi_landscape.png" \
      --figure-c "${RUN_DIR}/outputs/figures/fig06c_cytonuclear_paralog_cytotype_audit.png" \
      --run-id "${RUN_ID}" \
      ${FOCAL_ARGS[@]+"${FOCAL_ARGS[@]}"} \
      ${INTRASPECIFIC_ARGS[@]+"${INTRASPECIFIC_ARGS[@]}"} \
      ${GENO_R_DB_ARGS[@]+"${GENO_R_DB_ARGS[@]}"} \
      ${DEMO_FLAG[@]+"${DEMO_FLAG[@]}"}
    GENOMIC_STEP07_ARGS=(--genomic-dir "${RUN_DIR}/outputs/genomics")
  else
    echo "[Step 06/07] No genomic files staged in data/genomics/ - skipping standalone Step 06"
  fi
fi

# ------------------------------------------------------------------------------
# Step 07: Multi-Evidence Taxonomic Synthesis, Complex Delimitation & Dashboard
# ------------------------------------------------------------------------------
if [[ "${START_STEP}" -le 7 ]]; then
  echo "[Step 07/07] Running Integrative Taxonomic Synthesis & Complex Delimitation..."
  if [[ -n "${FOCAL_TAXON}" ]]; then
    FOCAL_TAXON_SLUG="$(printf '%s' "${FOCAL_TAXON}" | sed -E 's/[^A-Za-z0-9]+/_/g; s/^_+|_+$//g')"
    STEP07_FIGURE="${RUN_DIR}/outputs/figures/Figure_Integrative_${FOCAL_TAXON_SLUG}_Revision.png"
    STEP07_REPORT="${RUN_DIR}/outputs/reports/${FOCAL_TAXON_SLUG}_Taxonomic_Revision_Summary.md"
  else
    STEP07_FIGURE="${RUN_DIR}/outputs/figures/fig07a_integrative_taxonomic_synthesis.png"
    STEP07_REPORT="${RUN_DIR}/outputs/reports/Packera_Taxonomic_Revision_Summary.md"
  fi
  Rscript scripts/analysis/07_triage_dashboard_synthesis.R \
    --vouchers "${RUN_DIR}/data/tables/curated_vouchers.csv" \
    --harmonics "${RUN_DIR}/data/tables/leaf_efa_harmonics.csv" \
    --morphometrics "${RUN_DIR}/data/tables/morphometric_flags.csv" \
    --multimodal-flags "${RUN_DIR}/data/tables/multimodal_conflict_flags.csv" \
    --gmm-summary "${RUN_DIR}/outputs/reports/gmm_bayes_factors_summary.csv" \
    --niche-summary "${RUN_DIR}/outputs/reports/multimodal_spatial_rf_summary.csv" \
    --output-queue "${RUN_DIR}/data/tables/triage_queue.csv" \
    --output-plot "${STEP07_FIGURE}" \
    --output-report "${STEP07_REPORT}" \
    --output-html "${RUN_DIR}/outputs/reports/Triage_Dashboard_Interactive.html" \
    --min-affinity "${MIN_AFFINITY}" \
    --output-complexes "${RUN_DIR}/outputs/reports/empirical_taxon_complexes.csv" \
    --output-affinity "${RUN_DIR}/outputs/reports/pairwise_taxon_affinity_matrix.csv" \
    --output-neighbors "${RUN_DIR}/outputs/reports/focal_auto_neighbors.txt" \
    --output-complex-plot "${RUN_DIR}/outputs/figures/fig07b_taxon_complex_delimitation_network.png" \
    ${GENOMIC_STEP07_ARGS[@]+"${GENOMIC_STEP07_ARGS[@]}"} \
    "${FOCAL_ARGS[@]}" \
    "${INTRASPECIFIC_ARGS[@]}" \
    "${DEMO_FLAG[@]}"
fi

_pipe_finalize_run_metadata
