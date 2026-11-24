"""
===============================================================================
Module: taxonomy_authority.py
Project: Packera Species Delimitation & Morphometrics Pipeline
Affiliation: University of North Carolina at Chapel Hill Herbarium (NCU)

Description:
    Backward-compatible facade for Step 01 Darwin Core exsiccatae normalization,
    duplicate stratification, spatial balancing/thinning, 3-layer dynamic
    determiner authority engine, ICN nomenclatural type classification, Eastern
    NA cohort gating, and cascading outgroup anchor selection.
===============================================================================
"""

from __future__ import annotations

import logging
import math
from pathlib import Path
import re
from typing import Optional
import unicodedata

try:
    import numpy as np
except ImportError:
    np = None  # type: ignore[assignment]
try:
    import pandas as pd
except ImportError:
    pd = None  # type: ignore[assignment]
try:
    from PIL import Image
except ImportError:
    Image = None  # type: ignore[assignment]

if Image is not None:
    Image.MAX_IMAGE_PIXELS = None

PROJECT_ROOT = Path(__file__).resolve().parents[2]

from scripts.core.config import (
    CIS_MISSISSIPPI_STATES_PROVINCES,
    DEFAULT_CIS_MISSISSIPPI_LON,
    DEFAULT_COUNTRIES,
    DEFAULT_GEOGRAPHIC_SCOPE,
    DEFAULT_HOME_HERBARIUM,
    DEFAULT_MIN_CIS_MISSISSIPPI_FRAC,
    EXCLUDED_WESTERN_STATES,
    GREAT_PLAINS_STRADDLING_STATES,
    MAJOR_HERBARIA_CODES,
    PRIMARY_TYPE_STATUSES,
    SECONDARY_TYPE_STATUSES,
    SPECIALIST_PATTERNS,
    VALID_GEOGRAPHIC_SCOPES,
    WESTERN_LONGITUDE_THRESHOLD,
    normalize_countries_list,
)
from scripts.core.logger import setup_logging
from scripts.core.step01_concept_authority_utils import (
    CORE_PACKERA_GOLD_REGEXES,
    CORE_PACKERA_GOLD_SURNAMES,
    REGIONAL_AND_HISTORICAL_SILVER_REGEXES,
    REGIONAL_AND_HISTORICAL_SILVER_SURNAMES,
    apply_concept_and_segregate_recency_caps,
)

EXPORT_COLUMNS = [
    "catalogNumber",
    "institutionCode",
    "scientificName",
    "species_raw",
    "recordedBy",
    "recordNumber",
    "identifiedBy",
    "determiner_tier",
    "type_status",
    "county",
    "stateProvince",
    "decimalLatitude",
    "decimalLongitude",
    "coordinateUncertainty",
    "year",
    "month",
    "day",
    "eventDate",
    "regional_group",
    "image_path",
    "exsiccatae_key",
    "is_primary_duplicate",
    "duplicate_count",
]


def setup_logger(log_file_path: Optional[Path] = None, verbose: bool = False) -> logging.Logger:
    """Configures or retrieves a formatted logger for voucher harvesting."""
    return setup_logging(log_file=log_file_path, verbose=verbose, name="VoucherHarvester")


def sanitize_filename(name: str) -> str:
    """Sanitizes arbitrary strings into safe, valid filesystem filenames across platforms."""
    clean = re.sub(r'[\\/*?:"<>|\s]+', "_", str(name).strip())
    clean = clean.strip("._")
    return clean if clean else "voucher_unknown"


from scripts.core.step01_exsiccatae_fingerprint import (  # noqa: E402
    _as_bool_flag,
    _is_extralimital_status,
    _resolve_row_standardized_packera,
    generate_exsiccatae_key,
    get_image_resolution_score,
    haversine_distance_meters,
    normalize_collector,
    normalize_event_date,
    normalize_record_number,
    resolve_voucher_image_path,
)
from scripts.core.step01_duplicate_stratifier import stratify_duplicates  # noqa: E402
from scripts.core.step01_exsiccatae_reconcile import (  # noqa: E402
    reconcile_cross_taxon_exsiccatae_with_db,
)
from scripts.core.step01_spatial_thinning import (  # noqa: E402
    spatially_balance_and_thin_vouchers,
)
from scripts.core.step01_determiner_surnames import (  # noqa: E402
    CORE_PACKERA_MONOGRAPHERS,
    GENERIC_DETERMINER_TOKENS,
    PACKERA_SPECIALISTS,
    REGIONAL_FLORA_SPECIALISTS,
    _AMBIGUOUS_COMMON_SURNAMES,
    _COLLECTOR_REFERENCE_REGEX,
    _CORE_QUALIFIED_SPECIALIST_REGEXES,
    _FLORA_CITATION_REGEX,
    _HOLOTYPE_TIER_STATUSES,
    _INVALID_INSTITUTION_CODES,
    _ISOTYPE_TIER_STATUSES,
    _MISSING_IDENTIFIED_BY_TOKENS,
    _NEGATED_TYPE_RE,
    _NON_PERSON_DETERMINATION_TOKENS,
    _PARATYPE_TIER_STATUSES,
    _QUALIFIED_SPECIALIST_REGEXES,
    _REGIONAL_QUALIFIED_SPECIALIST_REGEXES,
    _SYNTYPE_TIER_STATUSES,
    _fold_ascii,
    _has_core_monographer_match,
    _has_non_person_determination_token,
    _has_regional_or_historical_specialist_match,
    _has_specialist_match,
    _is_explicit_collector_reference,
    _is_missing_identified_by,
    _is_verified_type_status,
    _should_fallback_to_collector,
    classify_nomenclatural_type_rank,
    classify_nomenclatural_type_subrank,
    normalize_determiner_surnames,
)
from scripts.core.step01_corpus_authority_discovery import (  # noqa: E402
    discover_corpus_authorities,
)
from scripts.core.step01_corpus_authority_tiers import (  # noqa: E402
    apply_corpus_authority_tiers,
    parse_determiner_tier,
)
from scripts.core.step01_eastern_cohort_scope import (  # noqa: E402
    _clean_upper_state_province,
    evaluate_eastern_na_cohort_filter,
    is_cis_mississippi_record,
    is_excluded_western_region,
    resolve_geographic_scope_settings,
)
from scripts.core.step01_outgroup_anchors import (  # noqa: E402
    select_cascading_outgroup_anchors,
)
