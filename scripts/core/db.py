"""
===============================================================================
Module: db.py
Project: Packera Species Delimitation & Morphometrics Pipeline
Affiliation: University of North Carolina at Chapel Hill Herbarium (NCU)

Description:
    Public facade for the central transactional SQLite (WAL mode) + lazy DuckDB
    Two-Tier Storage relational layer (`data/plantmatics.sqlite`, with automatic
    migration from legacy `data/packera_pipeline.sqlite`). Relies strictly on
    Python's standard library at import time for `.venv_LM3` ABI compatibility.
===============================================================================
"""

from __future__ import annotations

from scripts.core.db_cohort_and_env_cache import (
    _normalize_semicolon_list_str, normalize_run_id_token, populate_run_active_cohort,
    query_environmental_bbox, query_misidentified_vouchers_from_db,
    upsert_environmental_cache_records,
)
from scripts.core.db_cohort_consolidation import consolidate_false_multi_plant_splits_db
from scripts.core.db_column_constants import (
    DEFAULT_SQLITE_REL_PATH, DWC_VOUCHER_COLUMNS, ENVIRONMENTAL_BOOL_COLS,
    ENVIRONMENTAL_CACHE_COLUMNS, ENVIRONMENTAL_TEXT_COLS, HARMONIC_COLS,
    LEGACY_SQLITE_REL_PATH, PC_COLS, PROJECT_ROOT, SYM_HARMONIC_COLS,
    WC_BIO_COLUMNS, _MIGRATED_LEGACY_SQLITE_PATHS, _is_missing, _to_bool_int,
    _to_float, _to_int, _to_str, format_lat_lon_key, logger,
)
from scripts.core.db_conn_pragmas import (
    _apply_sqlite_pragmas, _ensure_dynamic_table_migrated,
    _migrate_legacy_sqlite_file, _migrate_legacy_sqlite_if_needed,
    get_db_connection, get_duckdb_connection, get_table_columns,
    is_canonical_or_explicit_db_context, resolve_db_path, upsert_rows_generic,
)
from scripts.core.db_csv_ingest_sync import ingest_existing_csvs_to_db
from scripts.core.db_genomics_cytology_persist import (
    persist_genetic_lineage_assignments, upsert_cytological_records,
    upsert_genomic_samples, upsert_natureserve_subnational_ranges,
    upsert_taxon_backbone_profiles, upsert_taxon_cytology_summary,
)
from scripts.core.db_morph_and_rf_persist import (
    _LEAF_EFA_BOOL_COLS, _LEAF_EFA_FLOAT_COLS, _LEAF_EFA_INT_COLS,
    _LEAF_EFA_STR_COLS, _OTU_MORPH_BOOL_COLS, _OTU_MORPH_FLOAT_COLS,
    _OTU_MORPH_INT_COLS, _OTU_MORPH_STR_COLS, _SRF_BOOL_COLS, _SRF_FLOAT_COLS,
    _SRF_INT_COLS, _SRF_LIST_STR_COLS, _SRF_STR_COLS, persist_leaf_efa_harmonics,
    persist_otu_morphometrics, persist_spatial_rf_predictions,
)
from scripts.core.db_purge_maintenance import purge_runs_only, purge_segments_and_runs
from scripts.core.db_run_tables_ingest import (
    _read_csv_dicts, ingest_existing_run_tables_into_db,
)
from scripts.core.db_schema import (
    MASTER_VIEW_DDL, PACKERA_SCHEMA_USER_VERSION, PLANTMATICS_SCHEMA_USER_VERSION,
    SCHEMA_DDL, _CANONICAL_EXTRACTED_LEAVES_COLS, _CANONICAL_PLANT_INDIVIDUALS_COLS,
    _CANONICAL_VOUCHER_REPRO_COLS, _DYNAMIC_TABLE_MIGRATION_SPECS,
    _MASTER_VIEW_REQUIRED_TABLES, _VOUCHERS_MIGRATION_COLUMNS,
    _build_master_view_ddl, _build_schema_ddl, _migrate_dynamic_schema_columns,
    _migrate_vouchers_schema_columns, _rebuild_tier_a_step02_tables_if_needed,
    _table_needs_canonical_rebuild, _tier_a_step02_needs_rebuild,
    init_db_schema, migrate_database_schema,
)
from scripts.core.db_step02_contours_tx import (
    _resolve_contour_file, read_contour_csv_xy,
    record_voucher_segmentation_transaction,
)
from scripts.core.db_step02_loader import (
    _PidAwareDict, ensure_pipeline_run, load_step02_records_from_db,
)
from scripts.core.db_vouchers_synonyms import (
    _CSV_BOOL_VOUCHER_COLUMNS, _format_voucher_csv_value,
    _infer_species_standardized_fallback, get_completed_qc_vouchers,
    load_taxon_synonyms, reconcile_curated_vouchers_csv_with_db,
    upsert_taxon_synonyms, upsert_vouchers_records,
)

__all__ = (
    "logger", "PROJECT_ROOT", "DEFAULT_SQLITE_REL_PATH", "LEGACY_SQLITE_REL_PATH",
    "_MIGRATED_LEGACY_SQLITE_PATHS", "DWC_VOUCHER_COLUMNS", "_is_missing", "_to_str", "_to_float",
    "_to_int", "_to_bool_int", "format_lat_lon_key", "_migrate_legacy_sqlite_file",
    "_migrate_legacy_sqlite_if_needed", "resolve_db_path", "is_canonical_or_explicit_db_context",
    "_apply_sqlite_pragmas", "get_db_connection", "get_duckdb_connection", "get_table_columns",
    "upsert_rows_generic", "_ensure_dynamic_table_migrated", "upsert_taxon_synonyms",
    "load_taxon_synonyms", "_infer_species_standardized_fallback", "upsert_vouchers_records",
    "_CSV_BOOL_VOUCHER_COLUMNS", "_format_voucher_csv_value",
    "reconcile_curated_vouchers_csv_with_db", "get_completed_qc_vouchers", "read_contour_csv_xy",
    "_resolve_contour_file", "record_voucher_segmentation_transaction", "_PidAwareDict",
    "load_step02_records_from_db", "ensure_pipeline_run", "consolidate_false_multi_plant_splits_db",
    "populate_run_active_cohort", "normalize_run_id_token", "query_misidentified_vouchers_from_db",
    "_normalize_semicolon_list_str", "upsert_environmental_cache_records",
    "query_environmental_bbox", "_LEAF_EFA_BOOL_COLS", "_LEAF_EFA_INT_COLS", "_LEAF_EFA_STR_COLS",
    "_LEAF_EFA_FLOAT_COLS", "persist_leaf_efa_harmonics", "_OTU_MORPH_BOOL_COLS",
    "_OTU_MORPH_INT_COLS", "_OTU_MORPH_STR_COLS", "_OTU_MORPH_FLOAT_COLS",
    "persist_otu_morphometrics", "_SRF_BOOL_COLS", "_SRF_INT_COLS", "_SRF_STR_COLS",
    "_SRF_LIST_STR_COLS", "_SRF_FLOAT_COLS", "persist_spatial_rf_predictions", "_read_csv_dicts",
    "ingest_existing_run_tables_into_db", "upsert_genomic_samples", "upsert_cytological_records",
    "upsert_taxon_cytology_summary", "upsert_taxon_backbone_profiles",
    "upsert_natureserve_subnational_ranges", "persist_genetic_lineage_assignments",
    "ingest_existing_csvs_to_db", "purge_segments_and_runs", "purge_runs_only", "HARMONIC_COLS",
    "SYM_HARMONIC_COLS", "PC_COLS", "WC_BIO_COLUMNS", "ENVIRONMENTAL_CACHE_COLUMNS",
    "ENVIRONMENTAL_TEXT_COLS", "ENVIRONMENTAL_BOOL_COLS", "_build_master_view_ddl",
    "_build_schema_ddl", "MASTER_VIEW_DDL", "SCHEMA_DDL", "_VOUCHERS_MIGRATION_COLUMNS",
    "_DYNAMIC_TABLE_MIGRATION_SPECS", "_MASTER_VIEW_REQUIRED_TABLES",
    "_migrate_vouchers_schema_columns", "PLANTMATICS_SCHEMA_USER_VERSION",
    "PACKERA_SCHEMA_USER_VERSION", "_CANONICAL_PLANT_INDIVIDUALS_COLS",
    "_CANONICAL_EXTRACTED_LEAVES_COLS", "_CANONICAL_VOUCHER_REPRO_COLS",
    "_table_needs_canonical_rebuild", "_tier_a_step02_needs_rebuild",
    "_rebuild_tier_a_step02_tables_if_needed", "_migrate_dynamic_schema_columns",
    "migrate_database_schema", "init_db_schema",
)
