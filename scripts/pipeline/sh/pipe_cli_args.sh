#!/usr/bin/env bash
# scripts/pipeline/sh/pipe_cli_args.sh - Default vars, usage text, CLI option parser, and compute env export.

_pipe_init_defaults() {
  CUSTOM_RUN_ID=""
  CUSTOM_DB_PATH=""
  SKIP_HARVEST=false
  FORCE_FLAG=""
  RESEGMENT=false
  PURGE_SEGMENTS=false
  PURGE_RUNS_ONLY=false
  START_STEP=1
  START_STEP_EXPLICIT=false
  FOCAL_TAXON=""
  CUSTOM_TAXA_RAW=()
  EXCLUDE_TAXA_RAW=()
  INTRASPECIFIC_MODE="false"
  GOLD_ANCHORS=30
  GOLD_ONLY=false
  DETERMINER_TIER=""
  EXCLUDE_MISIDENTIFIED=false
  TRIAGE_SOURCE="ALL"
  GEOGRAPHIC_SCOPE=""
  INCLUDE_WESTERN=false
  COUNTRIES_RAW=()
  EXPAND_SYNONYMS_FLAG="--expand-synonyms"
  DISCOVER_COMPLEXES=false
  AUTO_NEIGHBORS=false
  MIN_AFFINITY="0.03"
  SPATIAL_THIN_KM="0"
  MAX_PER_STATE=""
  MAX_PER_REGION=""
  NO_SPATIAL_WEIGHTS=false
  USDA_COUNTY_FILTER=""
  USDA_COUNTY_BUFFER_KM=""
  USDA_EXEMPT_GOLD=""
  SKIP_ENV_DOWNLOAD=false
  FORCE_ENV_DOWNLOAD=false
  DOWNLOAD_PUBLIC_GENOMICS=0
  SKIP_CYTOLOGY_DOWNLOAD=false
  FORCE_CYTOLOGY_DOWNLOAD=false
  DOWNLOAD_PUBLIC_CYTOLOGY=0
  SKIP_BACKBONE_DOWNLOAD=false
  FORCE_BACKBONE_DOWNLOAD=false
  DOWNLOAD_TAXONOMIC_BACKBONE=0
  LIMIT_ARG=()
  DEMO_FLAG=()
  FOCAL_ARGS=()
  INTRASPECIFIC_ARGS=()
  THREADS=""
  GPU_BATCH_SIZE=""
  WAVE_SIZE=""
  PREFETCH_WORKERS=""
  PREP_WORKERS=""
  POST_WORKERS=""
  COMPUTE_ARGS=()
  SKIP_NWI_NETWORK=false
  HARVEST_FULL_CONTINENTAL_RANGE=""
}

print_usage() {
  cat <<'EOF'
Usage: ./run_pipeline.sh [OPTIONS]

Run & Compute Controls:
  -h, --help                         Show this help message and exit
  --run-id <ID>                      Custom identifier for isolated run directory under runs/<ID>/
  --db <PATH>                        Override path to SQLite WAL database (default: data/plantmatics.sqlite)
  --start-step <1..7>                Start or resume execution at step <1..7> (step >=2 skips harvest automatically)
  --skip-harvest                     Skip Step 01 network harvest; seed cohort from SQLite / curated_vouchers.csv
  --force                            Force re-execution and overwrite cached artifacts (with --start-step 2, re-segments)
  --resegment                        Force Step 02 LeafMachine3 ONNX re-segmentation (--force-reprocess)
  --purge-segments, --purge-old-data Purge all previous Step 02 to 07 segments, masks, contours, runs, and SQLite tables (preserves Step 01 vouchers & raw images)
  --purge-runs-only                  Purge only Step 03 to 07 run outputs and Tier B SQLite tables (preserves Step 01 & Step 02)
  -j, --threads, --workers <INT>     CPU worker thread limit across Python & R stages (exports PACKERA_MAX_WORKERS)
  --gpu-batch-size, --batch-size <N> GPU inference batch size for Step 02 LeafMachine3 (YOLO26x / UNet++ v9 / ONNX Runtime)
  --wave-size <N>                    Cross-voucher GPU wave size in Step 02 (default: 16 or --gpu-batch-size)
  --prefetch-workers <N>             Override Stage A disk image prefetch worker thread count in Step 02
  --prep-workers <N>                 Override Stage B2/B3 CPU crop prep & mask decode worker count in Step 02
  --post-workers <N>                 Override Stage C contour/ECT/rosette post-processing worker count in Step 02
  --limit <N>                        Cap maximum number of vouchers processed in Step 02 segmentation
  --demo                             Run in ephemeral test isolation mode (PACKERA_DEMO_MODE=1)

Taxonomic Scope & Complex Discovery:
  --taxon "<Sp1>" "<Sp2>" ...        Target one or more Packera binomials (omit for genus-wide discovery)
  --focal-taxon "<Binomial>"         Target a focal Packera species for targeted synthesis
  --exclude-taxon "<Binomial>"       Exclude one or more specified taxa from the active run
  --intraspecific-mode               Single-taxon intraspecific morphotype/ecotype delimitation mode
  --expand-synonyms / --no-expand-synonyms
                                     Toggle GBIF synonym/basionym expansion (default: --expand-synonyms)
  --discover-complexes               Enable focal auto-neighbor complex scoping in unified Step 07 synthesis
  --auto-neighbors                   Append high-affinity neighbor taxa from Step 07 complex network output
  --min-affinity <FLOAT>             Minimum multimodal affinity weight threshold in Step 07 (default: 0.03)

Authority Tiering & Curation Filters:
  --gold-only                        Restrict active cohort to Nomenclatural Types and Tier_1_Gold vouchers
  --determiner-tier <TIER>           Restrict cohort to determiner authority tier (gold|silver|bronze or Tier_1_Gold)
  --gold-anchors <N>                 Target Tier 1 Gold reference anchors per outgroup taxon (default: 30)
  --exclude-misidentified            Prune vouchers flagged as Misidentification_Reassign in prior runs
  --triage-source <SOURCE>           Source run ID(s) or ALL for --exclude-misidentified (default: ALL)

Geographic Scope & Spatial Sampling:
  --geographic-scope <SCOPE>         Geographic scope: eastern_na (default), north_america, or global
  --include-western                  Include western North American vouchers west of -106.65W
  --harvest-full-continental-range / --no-harvest-full-continental-range
                                     Toggle full US + CA continental range harvesting for verified Eastern NA species
  --countries "<ISO2> ..."           ISO-2 country codes for GBIF harvesting (default: US)
  --spatial-thin-km <FLOAT>          Minimum spatial grid thinning distance in km per species (default: 0)
  --max-per-state <INT>              Cap maximum vouchers retained per (species, stateProvince)
  --max-per-region <INT>             Cap maximum vouchers retained per EPA Level III Ecoregion
  --no-spatial-weights               Disable inverse spatial-density weighting during reference model fitting
  --usda-county-filter / --no-usda-county-filter
                                     Toggle USDA PLANTS county range filter (1-ring neighbors + geodesic buffer)
  --usda-county-buffer-km <FLOAT>    Geodesic buffer in km around USDA PLANTS documented counties (default: 50.0)
  --usda-exempt-gold / --no-usda-exempt-gold
                                     Toggle Tier_1_Gold exemption for USDA county range filtering (default: true)
  --skip-env-download                Skip automatic preflight verification/download of SOLUS100, SGMC, and 3DEP layers
  --force-env-download               Force re-download of SOLUS100, SGMC, and 3DEP environmental layers
  --skip-nwi-network                 Reuse cached NWI wetland points without REST queries
  --download-public-genomics         Download public Packera sequences and trees into data/genomics/ before Step 06
  --skip-cytology-download           Skip automatic preflight download/sync of public Packera cytology records
  --force-cytology-download          Force fresh download of public Packera cytology records into data/cytology/
  --download-public-cytology         Download public Packera cytology records into data/cytology/ during preflight
  --skip-backbone-download           Skip automatic preflight download/sync of taxonomic backbone profiles
  --force-backbone-download          Force fresh download of taxonomic backbone profiles into data/taxonomy/
  --download-taxonomic-backbone      Download dynamic taxonomic backbone profiles into data/taxonomy/ during preflight
EOF
}

_pipe_append_triage_source() {
  local _val="$1"
  if [[ "${TRIAGE_SOURCE}" == "ALL" || -z "${TRIAGE_SOURCE}" ]]; then
    TRIAGE_SOURCE="${_val}"
  else
    TRIAGE_SOURCE="${TRIAGE_SOURCE},${_val}"
  fi
}

_pipe_consume_extra_arg() {
  _PIPE_SHIFT=1
  case "$1" in
    --threads=*|--workers=*) THREADS="${1#*=}" ;;
    --gpu-batch-size=*|--batch-size=*) GPU_BATCH_SIZE="${1#*=}" ;;
    --wave-size) WAVE_SIZE="$2"; _PIPE_SHIFT=2 ;;
    --wave-size=*) WAVE_SIZE="${1#*=}" ;;
    --prefetch-workers) PREFETCH_WORKERS="$2"; _PIPE_SHIFT=2 ;;
    --prefetch-workers=*) PREFETCH_WORKERS="${1#*=}" ;;
    --prep-workers) PREP_WORKERS="$2"; _PIPE_SHIFT=2 ;;
    --prep-workers=*) PREP_WORKERS="${1#*=}" ;;
    --post-workers) POST_WORKERS="$2"; _PIPE_SHIFT=2 ;;
    --post-workers=*) POST_WORKERS="${1#*=}" ;;
    --run-id) CUSTOM_RUN_ID="$2"; _PIPE_SHIFT=2 ;;
    --run-id=*) CUSTOM_RUN_ID="${1#*=}" ;;
    --db) CUSTOM_DB_PATH="$2"; _PIPE_SHIFT=2 ;;
    --db=*) CUSTOM_DB_PATH="${1#*=}" ;;
    --start-step) START_STEP="$2"; START_STEP_EXPLICIT=true; _PIPE_SHIFT=2 ;;
    --start-step=*) START_STEP="${1#*=}"; START_STEP_EXPLICIT=true ;;
    --skip-harvest) SKIP_HARVEST=true ;;
    --force|--overwrite) FORCE_FLAG="--force" ;;
    --resegment|--force-segment) RESEGMENT=true ;;
    --purge-segments|--purge-old-data|--purge-runs) PURGE_SEGMENTS=true; RESEGMENT=true ;;
    --limit) LIMIT_ARG=("--limit" "$2"); _PIPE_SHIFT=2 ;;
    --demo) DEMO_FLAG=("--demo") ;;
    --focal-taxon) FOCAL_TAXON="$2"; _PIPE_SHIFT=2 ;;
    --focal-taxon=*) FOCAL_TAXON="${1#*=}" ;;
    --exclude-taxon) EXCLUDE_TAXA_RAW+=("$2"); _PIPE_SHIFT=2 ;;
    --exclude-taxon=*) EXCLUDE_TAXA_RAW+=("${1#*=}") ;;
    --intraspecific-mode) INTRASPECIFIC_MODE="true" ;;
    --gold-anchors) GOLD_ANCHORS="$2"; _PIPE_SHIFT=2 ;;
    --gold-anchors=*) GOLD_ANCHORS="${1#*=}" ;;
    --taxon|--taxa)
      _PIPE_SHIFT=1; shift
      while [[ $# -gt 0 && "$1" != -* ]]; do
        CUSTOM_TAXA_RAW+=("$1"); _PIPE_SHIFT=$((_PIPE_SHIFT + 1)); shift
      done
      ;;
    --taxon=*|--taxa=*) CUSTOM_TAXA_RAW+=("${1#*=}") ;;
    --geographic-scope=*) GEOGRAPHIC_SCOPE="${1#*=}" ;;
    --countries)
      _PIPE_SHIFT=1; shift
      while [[ $# -gt 0 && "$1" != -* ]]; do
        COUNTRIES_RAW+=("$1"); _PIPE_SHIFT=$((_PIPE_SHIFT + 1)); shift
      done
      ;;
    --countries=*) COUNTRIES_RAW+=("${1#*=}") ;;
    --expand-synonyms) EXPAND_SYNONYMS_FLAG="--expand-synonyms" ;;
    --no-expand-synonyms) EXPAND_SYNONYMS_FLAG="--no-expand-synonyms" ;;
    --discover-complexes) DISCOVER_COMPLEXES=true ;;
    --auto-neighbors) AUTO_NEIGHBORS=true; DISCOVER_COMPLEXES=true ;;
    --min-affinity) MIN_AFFINITY="$2"; _PIPE_SHIFT=2 ;;
    --min-affinity=*) MIN_AFFINITY="${1#*=}" ;;
    --spatial-thin-km) SPATIAL_THIN_KM="$2"; _PIPE_SHIFT=2 ;;
    --spatial-thin-km=*) SPATIAL_THIN_KM="${1#*=}" ;;
    --max-per-state) MAX_PER_STATE="$2"; _PIPE_SHIFT=2 ;;
    --max-per-state=*) MAX_PER_STATE="${1#*=}" ;;
    --max-per-region) MAX_PER_REGION="$2"; _PIPE_SHIFT=2 ;;
    --max-per-region=*) MAX_PER_REGION="${1#*=}" ;;
    --no-spatial-weights) NO_SPATIAL_WEIGHTS=true ;;
    --usda-county-filter) USDA_COUNTY_FILTER="true" ;;
    --no-usda-county-filter) USDA_COUNTY_FILTER="false" ;;
    --usda-county-buffer-km) USDA_COUNTY_BUFFER_KM="$2"; _PIPE_SHIFT=2 ;;
    --usda-county-buffer-km=*) USDA_COUNTY_BUFFER_KM="${1#*=}" ;;
    --usda-exempt-gold) USDA_EXEMPT_GOLD="true" ;;
    --no-usda-exempt-gold) USDA_EXEMPT_GOLD="false" ;;
    --download-public-genomics) DOWNLOAD_PUBLIC_GENOMICS=1 ;;
    --skip-backbone-download) SKIP_BACKBONE_DOWNLOAD=true ;;
    --force-backbone-download) FORCE_BACKBONE_DOWNLOAD=true ;;
    --download-taxonomic-backbone) DOWNLOAD_TAXONOMIC_BACKBONE=1 ;;
    *) ;;
  esac
}

_pipe_finalize_cli_args() {
  CUSTOM_TAXA_CSV="$(IFS=,; echo "${CUSTOM_TAXA_RAW[*]:-}")"
  EXCLUDE_TAXON_CSV="$(IFS=,; echo "${EXCLUDE_TAXA_RAW[*]:-}")"

  case "${DETERMINER_TIER,,}" in
    gold|tier_1_gold)
      DETERMINER_TIER="Tier_1_Gold"
      GOLD_ONLY=true
      ;;
    silver|tier_2_silver)
      DETERMINER_TIER="Tier_2_Silver"
      ;;
    bronze|tier_3_bronze)
      DETERMINER_TIER="Tier_3_Bronze"
      ;;
  esac

  if [[ "${START_STEP}" -ge 2 ]]; then
    SKIP_HARVEST=true
  fi
  if [[ "${START_STEP_EXPLICIT}" == "true" && "${START_STEP}" -eq 2 && -n "${FORCE_FLAG} $([[ "${RESEGMENT}" == "true" ]] && echo "--force")" ]]; then
    RESEGMENT=true
  fi
  if [[ "${DISCOVER_COMPLEXES}" == "true" && -n "${FOCAL_TAXON}" ]]; then
    AUTO_NEIGHBORS=true
  fi
  if [[ -n "${THREADS}" ]]; then
    export PLANTMATICS_MAX_WORKERS="${THREADS}"
    export PACKERA_MAX_WORKERS="${THREADS}"
    COMPUTE_ARGS+=("--threads" "${THREADS}")
  fi
  if [[ -n "${GPU_BATCH_SIZE}" ]]; then
    export PLANTMATICS_GPU_BATCH_SIZE="${GPU_BATCH_SIZE}"
    export PACKERA_GPU_BATCH_SIZE="${GPU_BATCH_SIZE}"
    COMPUTE_ARGS+=("--gpu-batch-size" "${GPU_BATCH_SIZE}")
  fi
  if [[ -n "${WAVE_SIZE}" ]]; then
    export PLANTMATICS_WAVE_SIZE="${WAVE_SIZE}"
    export PACKERA_WAVE_SIZE="${WAVE_SIZE}"
    COMPUTE_ARGS+=("--wave-size" "${WAVE_SIZE}")
  fi
  if [[ -n "${PREFETCH_WORKERS}" ]]; then
    export PLANTMATICS_PREFETCH_WORKERS="${PREFETCH_WORKERS}"
    export PACKERA_PREFETCH_WORKERS="${PREFETCH_WORKERS}"
    COMPUTE_ARGS+=("--prefetch-workers" "${PREFETCH_WORKERS}")
  fi
  if [[ -n "${PREP_WORKERS}" ]]; then
    export PLANTMATICS_PREP_WORKERS="${PREP_WORKERS}"
    export PACKERA_PREP_WORKERS="${PREP_WORKERS}"
    COMPUTE_ARGS+=("--prep-workers" "${PREP_WORKERS}")
  fi
  if [[ -n "${POST_WORKERS}" ]]; then
    export PLANTMATICS_POST_WORKERS="${POST_WORKERS}"
    export PACKERA_POST_WORKERS="${POST_WORKERS}"
    COMPUTE_ARGS+=("--post-workers" "${POST_WORKERS}")
  fi

  if [[ -n "${FOCAL_TAXON}" ]]; then
    FOCAL_ARGS=("--focal-taxon" "${FOCAL_TAXON}")
  fi
  if [[ "${INTRASPECIFIC_MODE}" == "true" ]]; then
    INTRASPECIFIC_ARGS=("--intraspecific-mode")
  fi
}
