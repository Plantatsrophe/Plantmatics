"""
SQLite WAL connection factory, legacy path migration, DuckDB attachment, and generic batch UPSERT.
"""

from __future__ import annotations

import os
from pathlib import Path
import sqlite3
from typing import Any, List, Mapping, Optional, Sequence, Set, Tuple, Union

from scripts.core.db_column_constants import (
    DEFAULT_SQLITE_REL_PATH,
    LEGACY_SQLITE_REL_PATH,
    PROJECT_ROOT,
    _MIGRATED_LEGACY_SQLITE_PATHS,
    _is_missing,
    logger,
)
from scripts.core.db_schema import (
    _CANONICAL_EXTRACTED_LEAVES_COLS,
    _CANONICAL_PLANT_INDIVIDUALS_COLS,
    _CANONICAL_VOUCHER_REPRO_COLS,
    _DYNAMIC_TABLE_MIGRATION_SPECS,
    _migrate_dynamic_schema_columns,
    _table_needs_canonical_rebuild,
    init_db_schema,
)

def _migrate_legacy_sqlite_file(target_path: Path) -> Path:
    """
    Atomically migrate a legacy `packera_pipeline.sqlite` database (and its `-wal`/`-shm`
    sidecars) to `plantmatics.sqlite` in the same parent directory when `target_path`
    does not yet exist.
    """
    resolved_target = Path(target_path).resolve()
    if resolved_target.name != "plantmatics.sqlite":
        return resolved_target
    legacy_path = resolved_target.parent / "packera_pipeline.sqlite"
    if not resolved_target.exists() and legacy_path.exists():
        resolved_target.parent.mkdir(parents=True, exist_ok=True)
        os.replace(legacy_path, resolved_target)
        for sfx in ("-wal", "-shm"):
            legacy_sidecar = resolved_target.parent / f"packera_pipeline.sqlite{sfx}"
            if legacy_sidecar.exists():
                os.replace(legacy_sidecar, Path(f"{resolved_target}{sfx}"))
        key = str(resolved_target)
        if key not in _MIGRATED_LEGACY_SQLITE_PATHS:
            _MIGRATED_LEGACY_SQLITE_PATHS.add(key)
            logger.info("Migrated legacy SQLite database %s -> %s", legacy_path, resolved_target)
    return resolved_target

def _migrate_legacy_sqlite_if_needed(workspace_root: Path = PROJECT_ROOT) -> Path:
    """Resolve canonical `data/plantmatics.sqlite`, migrating legacy `data/packera_pipeline.sqlite` if needed."""
    canonical = (Path(workspace_root).resolve() / DEFAULT_SQLITE_REL_PATH).resolve()
    return _migrate_legacy_sqlite_file(canonical)

def resolve_db_path(
    db_path: Optional[Union[str, Path]] = None,
    workspace_root: Path = PROJECT_ROOT,
    demo: bool = False,
    run_dir: Optional[Union[str, Path]] = None,
) -> Union[str, Path]:
    """
    Resolve the SQLite database path while strictly isolating `--demo` runs from
    the canonical `data/plantmatics.sqlite` database (migrating legacy
    `data/packera_pipeline.sqlite` when applicable in non-demo contexts).
    """
    canonical = (Path(workspace_root).resolve() / DEFAULT_SQLITE_REL_PATH).resolve()
    legacy_canonical = (Path(workspace_root).resolve() / LEGACY_SQLITE_REL_PATH).resolve()

    if demo:
        if db_path is not None and str(db_path).strip() == ":memory:":
            return ":memory:"
        if db_path is not None and str(db_path).strip():
            cand = Path(str(db_path).strip()).resolve()
            if cand not in {canonical, legacy_canonical}:
                return cand
        eff_run = (
            (str(run_dir).strip() if run_dir is not None else "")
            or os.environ.get("PLANTMATICS_RUN_DIR", "").strip()
            or os.environ.get("PACKERA_RUN_DIR", "").strip()
        )
        if eff_run:
            rdir = Path(eff_run).resolve()
            if rdir.exists():
                return rdir / "demo_ephemeral.sqlite"
        return ":memory:"

    if db_path is not None and str(db_path).strip():
        raw = str(db_path).strip()
        if raw == ":memory:":
            return ":memory:"
        p = Path(raw)
        resolved = p.resolve() if p.is_absolute() else (Path(workspace_root).resolve() / p).resolve()
        if resolved.name == "plantmatics.sqlite":
            return _migrate_legacy_sqlite_file(resolved)
        return resolved

    env_db = (
        os.environ.get("PLANTMATICS_DB_PATH", "").strip()
        or os.environ.get("PACKERA_DB_PATH", "").strip()
    )
    if env_db:
        if env_db == ":memory:":
            return ":memory:"
        p = Path(env_db)
        resolved = p.resolve() if p.is_absolute() else (Path(workspace_root).resolve() / p).resolve()
        if resolved.name == "plantmatics.sqlite":
            return _migrate_legacy_sqlite_file(resolved)
        return resolved

    return _migrate_legacy_sqlite_if_needed(Path(workspace_root))

def is_canonical_or_explicit_db_context(
    target_path: Optional[Union[str, Path]] = None,
    workspace_root: Path = PROJECT_ROOT,
) -> bool:
    """
    Return True when `PLANTMATICS_DB_PATH` or `PACKERA_DB_PATH` is explicitly configured
    OR when `target_path` lives inside `workspace_root` (either canonical `data/` or
    `runs/<run_id>/`), and False when a standalone unit test writes into an isolated
    `/tmp` directory without setting either DB path environment variable.
    """
    if (
        os.environ.get("PLANTMATICS_DB_PATH", "").strip()
        or os.environ.get("PACKERA_DB_PATH", "").strip()
    ):
        return True
    if target_path is None:
        return True
    try:
        tp = Path(target_path).resolve()
        wr = Path(workspace_root).resolve()
        return tp == wr or wr in tp.parents
    except OSError:
        return False

def _apply_sqlite_pragmas(conn: sqlite3.Connection, read_only: bool = False) -> None:
    """Apply high-throughput workstation SQLite PRAGMAs (256 MiB page cache, 2 GiB mmap, RAM temp store)."""
    conn.execute("PRAGMA busy_timeout = 30000;")
    conn.execute("PRAGMA foreign_keys = ON;")
    conn.execute("PRAGMA cache_size = -262144;")
    conn.execute("PRAGMA temp_store = MEMORY;")
    try:
        conn.execute("PRAGMA mmap_size = 2147483648;")
    except sqlite3.Error:
        pass
    if not read_only:
        conn.execute("PRAGMA journal_mode = WAL;")
        conn.execute("PRAGMA synchronous = NORMAL;")
        try:
            conn.execute("PRAGMA wal_autocheckpoint = 4000;")
        except sqlite3.Error:
            pass

def get_db_connection(
    db_path: Optional[Union[str, Path]] = None,
    read_only: bool = False,
    demo: bool = False,
    init_schema: bool = False,
    workspace_root: Path = PROJECT_ROOT,
    check_same_thread: bool = False,
) -> sqlite3.Connection:
    """
    Open a SQLite connection configured with WAL mode, NORMAL synchronous durability,
    strict foreign-key enforcement, 256 MiB page cache, 2 GiB mmap, and a 30-second busy timeout.
    """
    resolved = resolve_db_path(db_path=db_path, workspace_root=workspace_root, demo=demo)

    if read_only and resolved != ":memory:":
        p = Path(resolved)
        if not p.exists():
            raise FileNotFoundError(f"SQLite database not found for read-only connection: {p}")
        conn = sqlite3.connect(
            f"file:{p.as_posix()}?mode=ro",
            uri=True,
            timeout=30.0,
            check_same_thread=check_same_thread,
        )
        conn.row_factory = sqlite3.Row
        _apply_sqlite_pragmas(conn, read_only=True)
        return conn

    if resolved != ":memory:":
        Path(resolved).parent.mkdir(parents=True, exist_ok=True)
        conn = sqlite3.connect(str(resolved), timeout=30.0, check_same_thread=check_same_thread)
    else:
        conn = sqlite3.connect(":memory:", timeout=30.0, check_same_thread=check_same_thread)

    conn.row_factory = sqlite3.Row
    _apply_sqlite_pragmas(conn, read_only=False)

    if init_schema:
        init_db_schema(conn)
    return conn

def get_duckdb_connection(
    db_path: Optional[Union[str, Path]] = None,
    read_only: bool = True,
    workspace_root: Path = PROJECT_ROOT,
) -> Any:
    """
    Lazily attach the SQLite database inside an in-memory DuckDB analytical session
    when `duckdb` is installed in `.venv`.
    """
    import duckdb  # Lazy import so .venv_LM3 never requires duckdb

    resolved = resolve_db_path(db_path=db_path, workspace_root=workspace_root, demo=False)
    dconn = duckdb.connect(database=":memory:")
    n_threads = max(1, min(32, os.cpu_count() or 8))
    try:
        dconn.execute(f"PRAGMA threads={n_threads};")
    except Exception:
        pass
    try:
        dconn.execute("LOAD sqlite;")
    except Exception:
        dconn.execute("INSTALL sqlite; LOAD sqlite;")
    ro_clause = " (TYPE SQLITE, READ_ONLY)" if read_only else " (TYPE SQLITE)"
    escaped_path = str(resolved).replace("'", "''")
    dconn.execute(f"ATTACH '{escaped_path}' AS packera{ro_clause};")
    dconn.execute("USE packera;")
    return dconn

def get_table_columns(conn: sqlite3.Connection, table_name: str) -> List[str]:
    cur = conn.execute(f"PRAGMA table_info({table_name});")
    return [str(r[1]) for r in cur.fetchall()]

def upsert_rows_generic(
    conn: sqlite3.Connection,
    table_name: str,
    rows: Sequence[Mapping[str, Any]],
    pk_cols: Sequence[str],
    preserve_existing_on_null: bool = True,
) -> int:
    """
    Parameterized batch UPSERT (`ON CONFLICT(...) DO UPDATE`) for arbitrary rows.
    """
    if not rows:
        return 0
    valid_cols = get_table_columns(conn, table_name)
    valid_set = set(valid_cols)

    present_cols: List[str] = []
    seen: Set[str] = set()
    for pk in pk_cols:
        if pk in valid_set and pk not in seen:
            present_cols.append(pk)
            seen.add(pk)
    for r in rows:
        for k in r.keys():
            if k in valid_set and k not in seen:
                present_cols.append(k)
                seen.add(k)

    non_pk_cols = [c for c in present_cols if c not in pk_cols]
    col_list_sql = ", ".join(present_cols)
    placeholders_sql = ", ".join("?" for _ in present_cols)
    pk_list_sql = ", ".join(pk_cols)

    if non_pk_cols:
        if preserve_existing_on_null:
            updates_sql = ", ".join(
                f"{c} = COALESCE(excluded.{c}, {table_name}.{c})" for c in non_pk_cols
            )
        else:
            updates_sql = ", ".join(f"{c} = excluded.{c}" for c in non_pk_cols)
        sql = (
            f"INSERT INTO {table_name} ({col_list_sql}) VALUES ({placeholders_sql}) "
            f"ON CONFLICT ({pk_list_sql}) DO UPDATE SET {updates_sql};"
        )
    else:
        sql = f"INSERT OR IGNORE INTO {table_name} ({col_list_sql}) VALUES ({placeholders_sql});"

    payload: List[Tuple[Any, ...]] = []
    for r in rows:
        vals: List[Any] = []
        for c in present_cols:
            v = r.get(c)
            if _is_missing(v):
                vals.append(None)
            elif isinstance(v, bool):
                vals.append(1 if v else 0)
            else:
                vals.append(v)
        payload.append(tuple(vals))

    with conn:
        conn.executemany(sql, payload)
    return len(payload)

def _ensure_dynamic_table_migrated(conn: sqlite3.Connection, table_name: str) -> None:
    canon_map = {
        "plant_individuals": _CANONICAL_PLANT_INDIVIDUALS_COLS,
        "extracted_leaves": _CANONICAL_EXTRACTED_LEAVES_COLS,
        "voucher_reproductive_phenology": _CANONICAL_VOUCHER_REPRO_COLS,
    }
    if table_name in canon_map and _table_needs_canonical_rebuild(conn, table_name, canon_map[table_name]):
        _migrate_dynamic_schema_columns(conn)
        return
    specs = _DYNAMIC_TABLE_MIGRATION_SPECS.get(table_name)
    if not specs:
        return
    existing = set(get_table_columns(conn, table_name))
    if existing and any(col not in existing for col, _ in specs):
        _migrate_dynamic_schema_columns(conn)

__all__ = (
    "_migrate_legacy_sqlite_file", "_migrate_legacy_sqlite_if_needed", "resolve_db_path",
    "is_canonical_or_explicit_db_context", "_apply_sqlite_pragmas", "get_db_connection",
    "get_duckdb_connection", "get_table_columns", "upsert_rows_generic",
    "_ensure_dynamic_table_migrated",
)
