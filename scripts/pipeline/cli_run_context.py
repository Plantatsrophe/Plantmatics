"""Shared CLI run-context helpers, env-var exports, preflight guards, and `main` patch bridge."""

from __future__ import annotations

import argparse
import logging
import os
from pathlib import Path
import sys
from typing import Any, Dict, Optional, Tuple, Union

PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from scripts.core.config import (
    PipelineConfig,
    get_active_run_dir,
    initialize_run_directory,
    resolve_canonical_run_id,
)
from scripts.core.db import (
    _read_csv_dicts,
    ensure_pipeline_run,
    get_db_connection,
    is_canonical_or_explicit_db_context,
    populate_run_active_cohort,
)
from scripts.core.logger import setup_logging

logger = setup_logging(name="PackeraPipeline")


def _main_attr(name: str, default: Any) -> Any:
    """Resolve attribute on root `main` module first so unittest.mock.patch('main.*') intercepts calls."""
    main_mod = sys.modules.get("main")
    if main_mod is not None and hasattr(main_mod, name):
        return getattr(main_mod, name)
    return default


def _apply_explicit_db_arg(args: Optional[argparse.Namespace], run_dir: Optional[Path] = None) -> None:
    """Export optional `--db`, `--demo`, `--threads`, and `--gpu-batch-size` CLI overrides to environment."""
    if args is None:
        return
    threads_arg = getattr(args, "threads", None)
    if threads_arg is not None and int(threads_arg) >= 1:
        os.environ["PLANTMATICS_MAX_WORKERS"] = str(int(threads_arg))
        os.environ["PACKERA_MAX_WORKERS"] = str(int(threads_arg))
    gpu_bs_arg = getattr(args, "gpu_batch_size", None)
    if gpu_bs_arg is not None and int(gpu_bs_arg) >= 1:
        os.environ["PLANTMATICS_GPU_BATCH_SIZE"] = str(int(gpu_bs_arg))
        os.environ["PACKERA_GPU_BATCH_SIZE"] = str(int(gpu_bs_arg))
    db_arg = getattr(args, "db", None)
    if db_arg is not None and str(db_arg).strip():
        raw = str(db_arg).strip()
        resolved_db = raw if raw == ":memory:" else str(Path(raw).resolve())
        os.environ["PLANTMATICS_DB_PATH"] = resolved_db
        os.environ["PACKERA_DB_PATH"] = resolved_db
    elif getattr(args, "demo", False):
        os.environ["PLANTMATICS_DEMO_MODE"] = "1"
        os.environ["PACKERA_DEMO_MODE"] = "1"
        if run_dir is not None:
            demo_db = str((run_dir / "data" / "demo_ephemeral.sqlite").resolve())
            os.environ["PLANTMATICS_DB_PATH"] = demo_db
            os.environ["PACKERA_DB_PATH"] = demo_db
        else:
            fallback_db = os.environ.get("PLANTMATICS_DB_PATH") or os.environ.get("PACKERA_DB_PATH") or ":memory:"
            os.environ["PLANTMATICS_DB_PATH"] = fallback_db
            os.environ["PACKERA_DB_PATH"] = fallback_db
    if getattr(args, "demo", False):
        os.environ["PLANTMATICS_DEMO_MODE"] = "1"
        os.environ["PACKERA_DEMO_MODE"] = "1"
    det_tier_arg = getattr(args, "determiner_tier", None)
    gold_only_arg = bool(getattr(args, "gold_only", False))
    if det_tier_arg is not None and str(det_tier_arg).strip():
        raw_tier = str(det_tier_arg).strip()
        tier_norm = {
            "gold": "Tier_1_Gold", "tier_1_gold": "Tier_1_Gold",
            "silver": "Tier_2_Silver", "tier_2_silver": "Tier_2_Silver",
            "bronze": "Tier_3_Bronze", "tier_3_bronze": "Tier_3_Bronze",
        }.get(raw_tier.lower(), raw_tier)
        os.environ["PLANTMATICS_DETERMINER_TIER"] = tier_norm
        os.environ["PACKERA_DETERMINER_TIER"] = tier_norm
        if tier_norm == "Tier_1_Gold":
            gold_only_arg = True
    if gold_only_arg:
        os.environ["PLANTMATICS_GOLD_ONLY"] = "1"
        os.environ["PACKERA_GOLD_ONLY"] = "1"
        if not (det_tier_arg is not None and str(det_tier_arg).strip()):
            os.environ["PLANTMATICS_DETERMINER_TIER"] = "Tier_1_Gold"
            os.environ["PACKERA_DETERMINER_TIER"] = "Tier_1_Gold"
    if getattr(args, "exclude_misidentified", False):
        os.environ["PLANTMATICS_EXCLUDE_MISIDENTIFIED"] = "1"
        os.environ["PACKERA_EXCLUDE_MISIDENTIFIED"] = "1"
    triage_src_arg = getattr(args, "triage_source", None)
    if triage_src_arg is not None:
        triage_src_str = (
            ",".join(str(x).strip() for x in triage_src_arg if str(x).strip())
            if isinstance(triage_src_arg, (list, tuple))
            else str(triage_src_arg).strip()
        )
        if triage_src_str:
            os.environ["PLANTMATICS_TRIAGE_SOURCE"] = triage_src_str
            os.environ["PACKERA_TRIAGE_SOURCE"] = triage_src_str
    if getattr(args, "auto_neighbors", False):
        os.environ["PLANTMATICS_AUTO_NEIGHBORS"] = "1"
        os.environ["PACKERA_AUTO_NEIGHBORS"] = "1"
    if getattr(args, "skip_nwi_network", False):
        os.environ["PLANTMATICS_SKIP_NWI_NETWORK"] = "1"
        os.environ["PACKERA_SKIP_NWI_NETWORK"] = "1"


def _db_has_extracted_contours(contours_dir: Path) -> bool:
    """Return True if SQLite `extracted_leaves` already holds 120-vertex contours for this workspace/DB."""
    _proj_root = _main_attr("PROJECT_ROOT", PROJECT_ROOT)
    if not is_canonical_or_explicit_db_context(contours_dir, workspace_root=_proj_root):
        return False
    try:
        conn = get_db_connection(read_only=True, workspace_root=_proj_root)
        try:
            cur = conn.execute(
                "SELECT COUNT(*) FROM extracted_leaves WHERE contour_xy_json IS NOT NULL AND LENGTH(contour_xy_json) > 2;"
            )
            return int((cur.fetchone() or [0])[0] or 0) > 0
        finally:
            conn.close()
    except Exception:
        return False


_COHORT_SYNC_CACHE: Dict[str, Tuple[str, int, int, int, int]] = {}


def _sync_run_cohort_in_db(
    run_dir: Path,
    args: Optional[argparse.Namespace] = None,
    cfg: Optional[Union[PipelineConfig, Dict[str, Any]]] = None,
    end_timestamp: Optional[str] = None,
) -> None:
    """Register `pipeline_runs` and synchronize `run_active_cohort` in SQLite when active."""
    if os.environ.get("PLANTMATICS_DEMO_MODE") == "1" or os.environ.get("PACKERA_DEMO_MODE") == "1":
        return
    _proj_root = _main_attr("PROJECT_ROOT", PROJECT_ROOT)
    _logger = _main_attr("logger", logger)
    if not is_canonical_or_explicit_db_context(run_dir, workspace_root=_proj_root):
        return
    focal_taxon, intraspecific_mode = (
        _resolve_focal_and_intraspecific(args, cfg or {})
        if args is not None
        else (None, None)
    )
    run_id = resolve_canonical_run_id(run_dir=run_dir, workspace_root=_proj_root)
    try:
        conn = get_db_connection(init_schema=True, workspace_root=_proj_root)
        try:
            ensure_pipeline_run(
                conn,
                run_id=run_id,
                run_dir=run_dir,
                end_timestamp=end_timestamp,
                focal_taxon=focal_taxon,
                intraspecific_mode=intraspecific_mode,
            )
            vouchers_csv = run_dir / "data" / "tables" / "curated_vouchers.csv"
            if vouchers_csv.is_file():
                vst = vouchers_csv.stat()
                if vst.st_size > 0:
                    pi_row = conn.execute("SELECT COUNT(*), COALESCE(MAX(rowid), 0) FROM plant_individuals;").fetchone() or (0, 0)
                    has_cohort = conn.execute(
                        "SELECT 1 FROM run_active_cohort WHERE run_id = ? LIMIT 1;",
                        (run_id,),
                    ).fetchone() is not None
                    fp = (str(vouchers_csv.resolve()), int(vst.st_mtime_ns), int(vst.st_size), int(pi_row[0] or 0), int(pi_row[1] or 0))
                    if has_cohort and _COHORT_SYNC_CACHE.get(run_id) == fp:
                        return
                    v_rows = _read_csv_dicts(vouchers_csv)
                    if v_rows:
                        populate_run_active_cohort(
                            conn,
                            run_id=run_id,
                            voucher_rows=v_rows,
                            clear_existing=True,
                        )
                        _COHORT_SYNC_CACHE[run_id] = fp
        finally:
            conn.close()
    except Exception as exc:
        _logger.debug("SQLite cohort registration skipped for %s: %s", run_dir, exc)


def _attach_run_log_handler(run_dir: Path) -> None:
    """Attach a FileHandler writing to `${RUN_DIR}/logs/pipeline.log` if not already attached."""
    log_path = (run_dir / "logs" / "pipeline.log").resolve()
    log_path.parent.mkdir(parents=True, exist_ok=True)
    _logger = _main_attr("logger", logger)
    for h in _logger.handlers:
        if isinstance(h, logging.FileHandler) and Path(getattr(h, "baseFilename", "")).resolve() == log_path:
            return
    fh = logging.FileHandler(log_path, encoding="utf-8")
    fh.setLevel(logging.DEBUG)
    fh.setFormatter(logging.Formatter("[%(asctime)s] [%(levelname)s] %(message)s", datefmt="%Y-%m-%d %H:%M:%S"))
    _logger.addHandler(fh)


def _resolve_active_run_context(
    args: argparse.Namespace,
    cfg: Union[PipelineConfig, Dict[str, Any]],
) -> Tuple[Optional[Path], Union[PipelineConfig, Dict[str, Any]]]:
    """
    If `--run-id` is passed on `args` or `PACKERA_RUN_DIR` is set in `os.environ`,
    initialize the run directory, attach `${RUN_DIR}/logs/pipeline.log`, and return
    `(run_dir, run_scoped_cfg)`. Otherwise return `(None, cfg)` for backward-compatible
    direct unit-test invocations.
    """
    _proj_root = _main_attr("PROJECT_ROOT", PROJECT_ROOT)
    _apply_db = _main_attr("_apply_explicit_db_arg", _apply_explicit_db_arg)
    _apply_db(args)
    explicit_run_id = getattr(args, "run_id", None)
    active_run_dir = get_active_run_dir(run_id=explicit_run_id, workspace_root=_proj_root)
    if active_run_dir is None:
        return None, cfg

    run_dir = initialize_run_directory(run_id=explicit_run_id, workspace_root=_proj_root)
    _apply_db(args, run_dir=run_dir)
    _main_attr("_attach_run_log_handler", _attach_run_log_handler)(run_dir)
    if isinstance(cfg, PipelineConfig):
        cfg = cfg.with_run_dir(run_dir)
    _main_attr("_sync_run_cohort_in_db", _sync_run_cohort_in_db)(run_dir, args=args, cfg=cfg)
    return run_dir, cfg


def _resolve_focal_and_intraspecific(
    args: argparse.Namespace,
    cfg: Union[PipelineConfig, Dict[str, Any]],
) -> Tuple[Optional[str], bool]:
    """Resolve optional focal taxon and intraspecific mode from CLI args or configuration."""
    taxa_cfg = cfg.get("taxa", {}) if hasattr(cfg, "get") else {}
    if not isinstance(taxa_cfg, dict):
        taxa_cfg = {}
    focal_raw = (
        getattr(args, "focal_taxon", None)
        or taxa_cfg.get("focal_taxon")
        or (cfg.get("focal_taxon") if hasattr(cfg, "get") else None)
    )
    focal_clean = str(focal_raw).strip() if focal_raw else None
    intraspecific = bool(
        getattr(args, "intraspecific_mode", False)
        or taxa_cfg.get("intraspecific_mode", False)
        or (cfg.get("intraspecific_mode", False) if hasattr(cfg, "get") else False)
    )
    return (focal_clean or None), intraspecific



# =============================================================================
# Defensive Preflight Sanity Checks
# =============================================================================

def check_gpu_availability(warn_only: bool = True) -> bool:
    """Checks whether PyTorch detects an active CUDA GPU device."""
    _logger = _main_attr("logger", logger)
    try:
        import torch
        if torch.cuda.is_available():
            device_name = torch.cuda.get_device_name(0)
            _logger.info(f"GPU Preflight: CUDA device detected: {device_name}")
            return True
        else:
            msg = "GPU Preflight: No CUDA-capable device detected by PyTorch."
            if warn_only:
                _logger.warning(f"{msg} Operations will run on CPU, which may be slower.")
            else:
                _logger.error(msg)
            return False
    except ImportError:
        _logger.warning("GPU Preflight: PyTorch not found in current environment.")
        return False


def verify_file_exists(file_path: Path, description: str, hint_cmd: Optional[str] = None) -> None:
    """Verifies that a required input file exists, failing early with clear instructions."""
    if not file_path.exists():
        err = f"Missing required {description} at: {file_path}"
        if hint_cmd:
            err += f"\n  -> Upstream phase skipped or incomplete! Please run: `{hint_cmd}`"
        _main_attr("logger", logger).error(err)
        sys.exit(1)


def verify_dir_has_files(dir_path: Path, pattern: str, description: str, hint_cmd: Optional[str] = None) -> None:
    """Verifies that a required directory exists and is non-empty."""
    if not dir_path.exists() or not any(dir_path.glob(pattern)):
        err = f"Missing or empty {description} directory at: {dir_path} (pattern: {pattern})"
        if hint_cmd:
            err += f"\n  -> Upstream phase skipped or incomplete! Please run: `{hint_cmd}`"
        _main_attr("logger", logger).error(err)
        sys.exit(1)
