"""
Unified botanical voucher ingestion, Darwin Core metadata normalization,
determiner authority evaluation, and asynchronous media acquisition engine
for the Packera Species Delimitation & Morphometrics Pipeline (NCU).
"""

from __future__ import annotations

import asyncio
import concurrent.futures
import datetime
import json
import logging
import math
import os
import re
import sys
import tempfile
import time
import urllib.request
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Set, Tuple, Union

try:
    import aiohttp
except ImportError:
    aiohttp = None  # type: ignore[assignment]
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


class _MissingPygbifSubmodule:
    """Mock-friendly stub when optional `pygbif` is not installed in the current environment."""

    _is_missing_pygbif = True

    def __init__(self, module_name: str) -> None:
        self._module_name = module_name

    def _raise(self, *_args: Any, **_kwargs: Any) -> Any:
        raise RuntimeError(
            f"Optional dependency '{self._module_name}' (pygbif) is required for live GBIF API queries. "
            "Install via `pip install pygbif`."
        )

    search = _raise
    name_backbone = _raise
    name_usage = _raise
    name_lookup = _raise


try:
    import pygbif.occurrences as occ
    import pygbif.species as species
except ImportError:
    occ = _MissingPygbifSubmodule("pygbif.occurrences")  # type: ignore[assignment]
    species = _MissingPygbifSubmodule("pygbif.species")  # type: ignore[assignment]

_AIOHTTP_CLIENT_ERRORS: Tuple[type, ...] = (aiohttp.ClientError,) if aiohttp is not None else ()
try:
    from shapely.geometry import Point
except ImportError:
    Point = None  # type: ignore[assignment]
try:
    from tqdm.asyncio import tqdm as async_tqdm
except ImportError:
    async_tqdm = None  # type: ignore[assignment]

# Allow decompression of high-resolution botanical herbarium sheets (>89 MP)
if Image is not None:
    Image.MAX_IMAGE_PIXELS = None

PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from scripts.core.config import (
    CIS_MISSISSIPPI_STATES_PROVINCES, DEFAULT_CIS_MISSISSIPPI_LON, DEFAULT_GEOGRAPHIC_SCOPE,
    DEFAULT_HOME_HERBARIUM, DEFAULT_INCLUDE_UNGEOREFERENCED_STATE_VOUCHERS,
    DEFAULT_MIN_CIS_MISSISSIPPI_FRAC, DEFAULT_MIN_FILE_SIZE_KB, DEFAULT_MIN_MEGAPIXELS,
    DEFAULT_MIN_SHARPNESS_LAPLACIAN, DEFAULT_OUTPUT_CSV, DEFAULT_RAW_DIR,
    DEFAULT_TARGET_TAXA, DEFAULT_USDA_COUNTY_BUFFER_KM, DEFAULT_USDA_COUNTY_FILTER,
    DEFAULT_USDA_EXEMPT_GOLD, DEFAULT_WORKSPACE, WESTERN_LONGITUDE_THRESHOLD,
    _DEFAULT_CONFIG, resolve_target_genus,
)
from scripts.core.coord_hygiene_utils import evaluate_coordinate_centroid_quarantine
from scripts.core.db import (
    get_db_connection, is_canonical_or_explicit_db_context, load_taxon_synonyms,
    resolve_db_path, upsert_taxon_synonyms, upsert_vouchers_records,
)
from scripts.core.epa_ecoregions import (
    _US_STATE_ABBR_TO_FULL, _normalize_us_state_name, infer_regional_group,
)
from scripts.core.step01_biogeography_utils import (
    normalize_packera_orthography, resolve_eastern_resident_whitelist, resolve_synonym_genera,
)
from scripts.core.step01_concept_authority_utils import (
    extract_identification_year, extract_typified_basionym, has_open_nomenclature_qualifier,
    is_interspecific_hybrid_formula, parse_voucher_taxonomic_concept,
)
from scripts.core.taxonomy_authority import (
    EXPORT_COLUMNS, _INVALID_INSTITUTION_CODES, apply_corpus_authority_tiers,
    classify_nomenclatural_type_rank, discover_corpus_authorities,
    evaluate_eastern_na_cohort_filter, is_cis_mississippi_record, is_excluded_western_region,
    parse_determiner_tier, reconcile_cross_taxon_exsiccatae_with_db,
    resolve_geographic_scope_settings, resolve_voucher_image_path,
    sanitize_filename, stratify_duplicates,
)
from scripts.core.usda_plants_ranges import filter_vouchers_by_usda_county_range

# Single shared mutable state and caches owned by scripts.core.harvester
_GBIF_SYNONYM_STATUSES: Set[str] = {"SYNONYM", "HOMOTYPIC_SYNONYM", "HETEROTYPIC_SYNONYM"}
_GBIF_REJECTED_SYNONYM_STATUSES: Set[str] = {"MISAPPLIED", "DOUBTFUL", "PROPARTE_SYNONYM"}
_DYNAMIC_SYNONYM_TO_CANONICAL: Dict[str, str] = {}
_DYNAMIC_TAXON_AUTHOR_STRINGS: Dict[str, Set[str]] = {}
_GBIF_SYNONYM_CACHE: Dict[str, Tuple[Optional[str], List[str], List[str]]] = {}
_GBIF_BACKBONE_SYNONYM_FLAG: Dict[str, bool] = {}
_GBIF_GENUS_KEY_CACHE: Dict[str, Optional[int]] = {}
_GBIF_GENUS_SPECIES_CACHE: Dict[int, List[str]] = {}
_GBIF_NETWORK_UNAVAILABLE: bool = False

from scripts.core.step01_harvester_filters import (
    _CA_PROVINCE_ALIASES_TO_KEY, _CULTIVATED_STATUS_TOKENS, _CULTIVATED_TEXT_RE,
    _FATAL_GBIF_GEO_ISSUES, _NA_STATE_PROVINCE_BBOXES, _OPEN_NOMENCLATURE_EPITHETS,
    _SYNONYM_GENERA, _extract_canonical_binomial, _extract_infraspecific_taxon,
    _is_cultivated_occurrence, _normalize_ca_province_bbox_key,
)

for _cfg_canon, _cfg_syns in (_DEFAULT_CONFIG.taxa.synonyms or {}).items():
    _canon_norm = _extract_canonical_binomial(_cfg_canon, {"Packera"}) or _cfg_canon.strip()
    if _canon_norm:
        _DYNAMIC_SYNONYM_TO_CANONICAL[_canon_norm.lower()] = _canon_norm
        for _cs in _cfg_syns or []:
            _cs_low = re.sub(r"\s+", " ", str(_cs).strip().lower())
            if _cs_low:
                _DYNAMIC_SYNONYM_TO_CANONICAL[_cs_low] = _canon_norm
            _cs_bin = _extract_canonical_binomial(_cs, _SYNONYM_GENERA)
            if _cs_bin:
                _DYNAMIC_SYNONYM_TO_CANONICAL[_cs_bin.lower()] = _canon_norm

from scripts.core.step01_harvester_synonyms import (
    _call_and_normalize_gbif_backbone, _get_packera_genus_key_gbif,
    _hydrate_synonyms_from_sqlite, _is_gbif_species_mocked,
    _persist_dynamic_synonyms_for_taxon, _register_dynamic_synonyms,
    discover_genus_packera_species_gbif,
)
from scripts.core.step01_harvester_taxon_std import (
    _query_gbif_backbone_and_synonyms, expand_taxon_synonyms_gbif, standardize_packera_taxon,
)
from scripts.core.step01_harvester_prescreen import prescreen_eastern_na_residents_gbif
from scripts.core.step01_harvester_media import (
    download_all_voucher_images, download_single_image, extract_high_res_image_url,
    optimize_herbarium_image_url, validate_image_quality,
)
from scripts.core.step01_harvester_webp_rescue import (
    download_and_validate_voucher_media, migrate_existing_vouchers_to_webp,
)
from scripts.core.step01_harvester_export import (
    _CONFLICT_STOP_TOKENS, _NEGATED_TAXON_RE, _detect_verbatim_identification_conflict,
    export_curated_table, print_and_log_summary,
)
from scripts.core.step01_harvester_record_ingest import ingest_gbif_occurrence_record
from scripts.core.step01_harvester_pagination import execute_gbif_harvest_loop

_hydrate_synonyms_from_sqlite()


class VoucherHarvester:
    """Consolidated botanical voucher harvester and curator for the genus `Packera`."""

    def __init__(
        self,
        taxa: Optional[List[str]] = None,
        max_uncertainty_meters: float = 5000.0, max_records_per_taxon: int = 5000,
        exclude_western: bool = True, min_megapixels: float = DEFAULT_MIN_MEGAPIXELS,
        min_file_size_kb: float = DEFAULT_MIN_FILE_SIZE_KB, check_sharpness: bool = False,
        min_sharpness: float = DEFAULT_MIN_SHARPNESS_LAPLACIAN, concurrency: int = 15,
        output_csv: Path = DEFAULT_OUTPUT_CSV, raw_dir: Path = DEFAULT_RAW_DIR,
        workspace_dir: Path = DEFAULT_WORKSPACE, logger: Optional[logging.Logger] = None,
        force: bool = False, geographic_scope: str = DEFAULT_GEOGRAPHIC_SCOPE,
        countries: Optional[Union[str, List[str]]] = None, include_western: bool = False,
        western_longitude_threshold: float = WESTERN_LONGITUDE_THRESHOLD,
        expand_synonyms: bool = False, cfg_synonyms: Optional[Dict[str, List[str]]] = None,
        home_herbarium: str = DEFAULT_HOME_HERBARIUM,
        db_path: Optional[Union[str, Path]] = None,
        usda_county_filter: bool = DEFAULT_USDA_COUNTY_FILTER,
        usda_county_buffer_km: float = DEFAULT_USDA_COUNTY_BUFFER_KM,
        usda_exempt_gold: bool = DEFAULT_USDA_EXEMPT_GOLD,
        usda_cache_dir: Optional[Path] = None,
        usda_county_geojson_path: Optional[Path] = None,
        harvest_full_continental_range: bool = False,
        min_eastern_residency_records: int = 10,
        cis_mississippi_longitude_threshold: float = DEFAULT_CIS_MISSISSIPPI_LON,
        min_cis_mississippi_fraction: float = DEFAULT_MIN_CIS_MISSISSIPPI_FRAC,
        include_ungeoreferenced_state_vouchers: bool = DEFAULT_INCLUDE_UNGEOREFERENCED_STATE_VOUCHERS,
        webp_quality: int = 92, **_unused_kwargs: Any,
    ):
        self.logger = logger or logging.getLogger("VoucherHarvester")
        self.db_path = db_path
        if self.db_path is not None:
            _hydrate_synonyms_from_sqlite(self.db_path)
        raw_input_taxa = [str(t).strip() for t in (taxa if taxa is not None else DEFAULT_TARGET_TAXA) if str(t).strip()]
        self.expand_synonyms = bool(expand_synonyms)
        self.cfg_synonyms = cfg_synonyms if cfg_synonyms is not None else dict(_DEFAULT_CONFIG.taxa.synonyms or {})

        explicit_genus_only = len(raw_input_taxa) == 1 and raw_input_taxa[0].lower() == "packera"
        scope, eff_excl_west, eff_countries = resolve_geographic_scope_settings(
            geographic_scope=geographic_scope, exclude_western=exclude_western,
            include_western=include_western, countries=countries,
        )
        if not raw_input_taxa:
            self.is_genus_wide = True
            self.logger.info("Discovering accepted Packera species across the genus via GBIF Backbone...")
            discovered_species = discover_genus_packera_species_gbif()
            base_taxa = discovered_species if discovered_species else ["Packera"]
            self.logger.info(f"Discovered {len(base_taxa)} accepted Packera base taxa.")
            if eff_excl_west and len(base_taxa) > 1:
                base_taxa = prescreen_eastern_na_residents_gbif(
                    base_taxa=base_taxa, eff_countries=eff_countries,
                    cis_mississippi_longitude_threshold=cis_mississippi_longitude_threshold,
                    min_eastern_residency_records=min_eastern_residency_records,
                    min_cis_mississippi_fraction=min_cis_mississippi_fraction,
                    db_path=self.db_path, workspace_dir=workspace_dir, logger=self.logger,
                )
        elif explicit_genus_only:
            self.is_genus_wide = True
            base_taxa = ["Packera"]
        else:
            self.is_genus_wide = False
            base_taxa = raw_input_taxa

        if self.is_genus_wide:
            self._protected_packera_taxa: Set[str] = {
                b.lower() for tx in base_taxa for b in [_extract_canonical_binomial(tx, {"Packera"})] if b
            }
        else:
            self._protected_packera_taxa = set()
            if not self.expand_synonyms:
                self._protected_packera_taxa = {
                    b.lower() for tx in raw_input_taxa for b in [_extract_canonical_binomial(tx, {"Packera"})] if b
                }

        if self.expand_synonyms and not (len(base_taxa) == 1 and base_taxa[0].lower() == "packera"):
            self.logger.info(
                f"Resolving GBIF Backbone synonyms (unaccepted Packera + Senecio/Cineraria basionyms) "
                f"across {len(base_taxa)} base taxa..."
            )
            if len(base_taxa) > 1 and not _is_gbif_species_mocked() and not _GBIF_NETWORK_UNAVAILABLE:
                uncached_taxa = [tx for tx in base_taxa if tx.strip().lower() not in _GBIF_SYNONYM_CACHE]
                if uncached_taxa:
                    with concurrent.futures.ThreadPoolExecutor(max_workers=min(12, len(uncached_taxa))) as syn_pool:
                        list(syn_pool.map(_query_gbif_backbone_and_synonyms, uncached_taxa))
            if not self.is_genus_wide:
                # Protect explicitly co-requested ACCEPTED Packera species while allowing unaccepted
                # synonyms (e.g. Packera tomentosa -> Packera dubia) to merge into their accepted species.
                for tx in base_taxa:
                    b_tx = _extract_canonical_binomial(tx, {"Packera"})
                    if b_tx and not _GBIF_BACKBONE_SYNONYM_FLAG.get(b_tx.lower(), False):
                        self._protected_packera_taxa.add(b_tx.lower())
            expanded_taxa: List[str] = []
            seen_t: Set[str] = set()
            for idx_tx, tx in enumerate(base_taxa, start=1):
                for syn_tx in expand_taxon_synonyms_gbif(
                    tx, cfg_synonyms=self.cfg_synonyms, protected_packera_taxa=self._protected_packera_taxa,
                ):
                    low_tx = syn_tx.lower()
                    if low_tx not in seen_t:
                        seen_t.add(low_tx)
                        expanded_taxa.append(syn_tx)
                if len(base_taxa) >= 15 and (idx_tx % 25 == 0 or idx_tx == len(base_taxa)):
                    self.logger.info(
                        f"  [Synonym Expansion] Resolved {idx_tx}/{len(base_taxa)} base taxa "
                        f"-> {len(expanded_taxa)} total query binomials..."
                    )
            self.taxa = expanded_taxa if expanded_taxa else base_taxa
        else:
            self.taxa = base_taxa

        scope, eff_excl_west, eff_countries = resolve_geographic_scope_settings(
            geographic_scope=geographic_scope, exclude_western=exclude_western,
            include_western=include_western, countries=countries,
        )
        self.geographic_scope = scope
        self.exclude_western = eff_excl_west
        self.countries = eff_countries
        self.western_longitude_threshold = float(western_longitude_threshold)
        self.home_herbarium = str(home_herbarium or DEFAULT_HOME_HERBARIUM).strip().upper() or "NCU"
        self.max_uncertainty_meters = max_uncertainty_meters
        self.max_records_per_taxon = max_records_per_taxon
        self.min_megapixels = min_megapixels
        self.min_file_size_kb = min_file_size_kb
        self.check_sharpness = check_sharpness
        self.min_sharpness = min_sharpness
        self.concurrency = concurrency
        self.output_csv = Path(output_csv)
        self.raw_dir = Path(raw_dir)
        self.workspace_dir = Path(workspace_dir)
        self.logger = logger or logging.getLogger("VoucherHarvester")
        self.force = force
        self.usda_county_filter = bool(usda_county_filter)
        self.usda_county_buffer_km = float(usda_county_buffer_km)
        self.usda_exempt_gold = bool(usda_exempt_gold)
        self.usda_cache_dir = Path(usda_cache_dir) if usda_cache_dir is not None else None
        self.usda_county_geojson_path = Path(usda_county_geojson_path) if usda_county_geojson_path is not None else None
        self.harvest_full_continental_range = bool(harvest_full_continental_range)
        self.min_eastern_residency_records = max(1, int(min_eastern_residency_records))
        self.cis_mississippi_longitude_threshold = float(cis_mississippi_longitude_threshold)
        self.min_cis_mississippi_fraction = float(min_cis_mississippi_fraction)
        self.include_ungeoreferenced_state_vouchers = bool(include_ungeoreferenced_state_vouchers)
        self.webp_quality = max(1, min(100, int(webp_quality)))
        self.last_usda_filter_stats: Optional[Dict[str, Any]] = None
        self._secondary_duplicate_pool: Optional[pd.DataFrame] = None
        self._pre_thin_curated_pool: Optional[pd.DataFrame] = None

    def harvest(self) -> pd.DataFrame:
        """Harvests and normalizes Darwin Core occurrence records across configured taxa."""
        return execute_gbif_harvest_loop(self)

    def download_and_validate_media(self, df: pd.DataFrame) -> Tuple[pd.DataFrame, Dict[str, Any]]:
        """Asynchronously downloads voucher specimen sheets, rescuing failed primaries via duplicate fallbacks."""
        return download_and_validate_voucher_media(self, df)

    def export(self, df: pd.DataFrame) -> Path:
        """Atomically persists the curated vouchers table to CSV and SQLite."""
        return export_curated_table(
            df, self.output_csv, logger=self.logger,
            home_herbarium=self.home_herbarium, db_path=self.db_path,
        )

    def run(
        self,
        download_images: bool = False,
        force: Optional[bool] = None,
    ) -> Tuple[pd.DataFrame, Optional[Dict[str, Any]]]:
        """Executes the end-to-end voucher ingestion and curation workflow."""
        if force is not None:
            self.force = force

        migrate_existing_vouchers_to_webp(
            self.raw_dir, db_path=self.db_path, curated_csv=self.output_csv,
            quality=self.webp_quality, logger=self.logger,
        )
        self.logger.info("Starting Packera Voucher Ingestion & Authority Stratification Pipeline...")
        self.logger.info(f"Target Taxa ({len(self.taxa)}): {self.taxa}")
        self.logger.info(
            f"Geographic Scope: {self.geographic_scope} (countries={self.countries or 'GLOBAL'}, "
            f"exclude_western={self.exclude_western})"
        )
        self.logger.info(f"Max Coordinate Uncertainty Threshold: {self.max_uncertainty_meters} m")
        self.logger.info(f"Force Overwrite Mode: {self.force}")

        df_curated = self.harvest()
        download_stats = None
        if download_images and not df_curated.empty:
            df_curated, download_stats = self.download_and_validate_media(df_curated)

        self.export(df_curated)
        print_and_log_summary(
            df_curated, download_stats, self.logger,
            usda_filter_stats=self.last_usda_filter_stats,
        )
        self.logger.info("Pipeline execution completed successfully.")
        return df_curated, download_stats


def harvest_taxa_occurrences(
    taxa_list: Optional[List[str]] = None,
    max_uncertainty_meters: float = 5000.0, max_records_per_taxon: int = 1000,
    exclude_western: bool = True, logger: Optional[logging.Logger] = None,
    geographic_scope: str = DEFAULT_GEOGRAPHIC_SCOPE,
    countries: Optional[Union[str, List[str]]] = None, include_western: bool = False,
    expand_synonyms: bool = False, usda_county_filter: bool = DEFAULT_USDA_COUNTY_FILTER,
    usda_county_buffer_km: float = DEFAULT_USDA_COUNTY_BUFFER_KM,
    usda_exempt_gold: bool = DEFAULT_USDA_EXEMPT_GOLD,
) -> pd.DataFrame:
    """
    Functional entry point for occurrence harvesting (defaults to `expand_synonyms=False`
    for deterministic programmatic/unit-test queries; pass `expand_synonyms=True` to enable
    dynamic GBIF Backbone synonym expansion).
    """
    harvester = VoucherHarvester(
        taxa=taxa_list, max_uncertainty_meters=max_uncertainty_meters,
        max_records_per_taxon=max_records_per_taxon, exclude_western=exclude_western,
        logger=logger, geographic_scope=geographic_scope, countries=countries,
        include_western=include_western, expand_synonyms=expand_synonyms,
        usda_county_filter=usda_county_filter, usda_county_buffer_km=usda_county_buffer_km,
        usda_exempt_gold=usda_exempt_gold,
    )
    return harvester.harvest()


def main() -> None:
    """CLI execution wrapper delegating to data prep runner."""
    import importlib
    cli_mod = importlib.import_module("scripts.data_prep.01_voucher_harvester")
    cli_mod.main()


if __name__ == "__main__":
    main()
