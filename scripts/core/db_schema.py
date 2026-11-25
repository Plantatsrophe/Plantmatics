"""
===============================================================================
Module: db_schema.py
Project: Packera Species Delimitation & Morphometrics Pipeline
Affiliation: University of North Carolina at Chapel Hill Herbarium (NCU)

Description:
    Public facade re-exporting the Two-Tier Relational Schema DDL (Tier A +
    Tier B + R*Tree + Master View) and idempotent dynamic migration helpers.
===============================================================================
"""

from __future__ import annotations

from scripts.core.db_schema_ddl_tier_b import SCHEMA_DDL, _build_schema_ddl
from scripts.core.db_schema_env_cols import (
    ENVIRONMENTAL_BOOL_COLS,
    ENVIRONMENTAL_CACHE_COLUMNS,
    ENVIRONMENTAL_TEXT_COLS,
    HARMONIC_COLS,
    PC_COLS,
    SYM_HARMONIC_COLS,
    WC_BIO_COLUMNS,
)
from scripts.core.db_schema_master_view_b import MASTER_VIEW_DDL, _build_master_view_ddl
from scripts.core.db_schema_migration_specs import (
    PACKERA_SCHEMA_USER_VERSION,
    PLANTMATICS_SCHEMA_USER_VERSION,
    _DYNAMIC_TABLE_MIGRATION_SPECS,
    _MASTER_VIEW_REQUIRED_TABLES,
    _VOUCHERS_MIGRATION_COLUMNS,
)
from scripts.core.db_schema_migrations import (
    _migrate_dynamic_schema_columns,
    _migrate_vouchers_schema_columns,
    init_db_schema,
    migrate_database_schema,
)
from scripts.core.db_schema_step02_rebuild import (
    _CANONICAL_EXTRACTED_LEAVES_COLS,
    _CANONICAL_PLANT_INDIVIDUALS_COLS,
    _CANONICAL_VOUCHER_REPRO_COLS,
    _rebuild_tier_a_step02_tables_if_needed,
    _table_needs_canonical_rebuild,
    _tier_a_step02_needs_rebuild,
)

__all__ = (
    "HARMONIC_COLS", "SYM_HARMONIC_COLS", "PC_COLS", "WC_BIO_COLUMNS",
    "ENVIRONMENTAL_CACHE_COLUMNS", "ENVIRONMENTAL_TEXT_COLS", "ENVIRONMENTAL_BOOL_COLS",
    "_build_master_view_ddl", "_build_schema_ddl", "MASTER_VIEW_DDL", "SCHEMA_DDL",
    "_VOUCHERS_MIGRATION_COLUMNS", "_DYNAMIC_TABLE_MIGRATION_SPECS", "_MASTER_VIEW_REQUIRED_TABLES",
    "_migrate_vouchers_schema_columns", "PLANTMATICS_SCHEMA_USER_VERSION",
    "PACKERA_SCHEMA_USER_VERSION", "_CANONICAL_PLANT_INDIVIDUALS_COLS",
    "_CANONICAL_EXTRACTED_LEAVES_COLS", "_CANONICAL_VOUCHER_REPRO_COLS",
    "_table_needs_canonical_rebuild", "_tier_a_step02_needs_rebuild",
    "_rebuild_tier_a_step02_tables_if_needed", "_migrate_dynamic_schema_columns",
    "migrate_database_schema", "init_db_schema",
)
