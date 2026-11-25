"""CLI execution handlers for Phase 1 (`harvest`) and Phase 2 (`segment`)."""

from __future__ import annotations

import argparse
from datetime import datetime
from pathlib import Path
import subprocess
import sys
from typing import List

from scripts.core.artifact_manager import ensure_model_weights
from scripts.core.config import (
    PipelineConfig,
    ensure_run_curated_vouchers_snapshot,
    write_run_metadata,
)
from scripts.pipeline.cli_run_context import (
    PROJECT_ROOT,
    _db_has_extracted_contours,
    _main_attr,
    _resolve_active_run_context,
    _sync_run_cohort_in_db,
    check_gpu_availability,
    logger,
    verify_dir_has_files,
    verify_file_exists,
)


def run_harvest(args: argparse.Namespace, cfg: PipelineConfig) -> None:
    """Executes Phase 1: Voucher harvesting and authority tier stratification."""
    _logger = _main_attr("logger", logger)
    _proj_root = _main_attr("PROJECT_ROOT", PROJECT_ROOT)
    _subprocess = _main_attr("subprocess", subprocess)
    _verify_file_exists = _main_attr("verify_file_exists", verify_file_exists)
    _resolve_ctx = _main_attr("_resolve_active_run_context", _resolve_active_run_context)
    _sync_cohort = _main_attr("_sync_run_cohort_in_db", _sync_run_cohort_in_db)
    _logger.info("=== Starting Phase 1: Voucher Ingestion & Authority Stratification ===")
    start_ts = datetime.now().astimezone().isoformat()
    run_dir, cfg = _resolve_ctx(args, cfg)

    script_path = _proj_root / "scripts" / "data_prep" / "01_voucher_harvester.py"
    if not script_path.exists():
        script_path = _proj_root / "scripts" / "pipeline" / "01_voucher_harvester.py"
    _verify_file_exists(script_path, "harvesting script")

    out_dir = Path(args.out_dir or cfg.get("paths", {}).get("raw_vouchers_dir", "data/raw_vouchers"))
    default_csv = (
        run_dir / "data" / "tables" / "curated_vouchers.csv"
        if run_dir is not None
        else cfg.get("paths", {}).get("curated_vouchers_csv", "data/tables/curated_vouchers.csv")
    )
    output_csv = Path(args.output_csv or default_csv)
    out_dir.mkdir(parents=True, exist_ok=True)
    output_csv.parent.mkdir(parents=True, exist_ok=True)

    harvest_conc = args.concurrency
    if harvest_conc is None:
        t_req = getattr(args, "threads", None)
        if t_req is not None and int(t_req) >= 1:
            harvest_conc = min(32, max(1, int(t_req)))
        else:
            harvest_conc = cfg["harvesting"]["download_concurrency"]
    cmd = [
        sys.executable, str(script_path),
        "--out-dir", str(out_dir),
        "--output-csv", str(output_csv),
        "--max-records", str(args.max_records or cfg["harvesting"]["max_records_per_taxon"]),
        "--max-uncertainty", str(args.max_uncertainty or cfg["harvesting"]["max_uncertainty_meters"]),
        "--min-megapixels", str(args.min_megapixels or cfg["harvesting"]["min_megapixels"]),
        "--min-file-size-kb", str(args.min_file_size_kb or cfg["harvesting"]["min_file_size_kb"]),
        "--concurrency", str(harvest_conc),
    ]
    cfg_target_species: List[str] = []
    if isinstance(cfg, PipelineConfig):
        cfg_target_species = list(cfg.taxa.target_species or [])
    elif isinstance(cfg, dict):
        cfg_target_species = list(
            (cfg.get("taxa") or {}).get("target_species")
            or (cfg.get("taxonomy") or {}).get("target_taxa")
            or []
        )
    explicit_taxa = list(args.taxa) if args.taxa else cfg_target_species
    if not explicit_taxa and getattr(args, "intraspecific_mode", False) and getattr(args, "focal_taxon", None):
        focal_clean = str(args.focal_taxon).strip()
        if focal_clean:
            explicit_taxa = [focal_clean]
    if explicit_taxa:
        cmd.extend(["--taxon", *explicit_taxa])
    if getattr(args, "geographic_scope", None):
        cmd.extend(["--geographic-scope", str(args.geographic_scope)])
    if getattr(args, "include_western", False):
        cmd.append("--include-western")
    if getattr(args, "countries", None):
        cmd.extend(["--countries", *[str(c) for c in args.countries]])
    if getattr(args, "expand_synonyms", None) is False:
        cmd.append("--no-expand-synonyms")
    if getattr(args, "usda_county_filter", None) is True:
        cmd.append("--usda-county-filter")
    elif getattr(args, "usda_county_filter", None) is False:
        cmd.append("--no-usda-county-filter")
    if getattr(args, "usda_county_buffer_km", None) is not None:
        cmd.extend(["--usda-county-buffer-km", str(args.usda_county_buffer_km)])
    if getattr(args, "usda_exempt_gold", None) is True:
        cmd.append("--usda-exempt-gold")
    elif getattr(args, "usda_exempt_gold", None) is False:
        cmd.append("--no-usda-exempt-gold")
    if getattr(args, "harvest_full_continental_range", None) is True:
        cmd.append("--harvest-full-continental-range")
    elif getattr(args, "harvest_full_continental_range", None) is False:
        cmd.append("--no-harvest-full-continental-range")

    if args.download_images:
        cmd.append("--download-images")
    if getattr(args, "force", False):
        cmd.append("--force")
    if args.verbose:
        cmd.append("--verbose")

    _logger.info(f"Running command: {' '.join(cmd)}")
    _subprocess.run(cmd, check=True)
    _verify_file_exists(output_csv, "harvested vouchers CSV output")
    if run_dir is not None:
        _sync_cohort(run_dir, args=args, cfg=cfg)
        write_run_metadata(run_dir, start_time=start_ts)
    _logger.info("Phase 1 completed successfully.")


def get_lm3_python_executable() -> Path:
    """Resolves the Python interpreter executable for LeafMachine3 execution.

    Checks for:
      - Unix/Ubuntu: LeafMachine3/.venv/bin/python and .venv_LM3/bin/python
      - Windows: LeafMachine3/.venv/Scripts/python.exe and .venv_LM3/Scripts/python.exe

    Returns:
        Path: Path to the LM3 virtual environment Python interpreter if found,
        otherwise falls back to sys.executable and logs a warning advising setup
        via setup_leafmachine3.sh.
    """
    _proj_root = _main_attr("PROJECT_ROOT", PROJECT_ROOT)
    _logger = _main_attr("logger", logger)
    venv_lm3_primary = _proj_root / "LeafMachine3" / ".venv"
    venv_lm3_alt = _proj_root / ".venv_LM3"
    candidates = [
        venv_lm3_primary / "bin" / "python",
        venv_lm3_alt / "bin" / "python",
        venv_lm3_primary / "Scripts" / "python.exe",
        venv_lm3_alt / "Scripts" / "python.exe",
    ]
    for candidate in candidates:
        if candidate.exists():
            return candidate

    _logger.warning(
        "LeafMachine3 dedicated virtual environment (LeafMachine3/.venv or .venv_LM3) not found at %s. "
        "Falling back to current active interpreter (%s). "
        "Please run 'bash setup_leafmachine3.sh' to configure the dedicated LM3 environment.",
        venv_lm3_primary,
        sys.executable,
    )
    return Path(sys.executable)


def run_segment(args: argparse.Namespace, cfg: PipelineConfig) -> None:
    """Executes Phase 2: LeafMachine3 (YOLO26x / UNet++ v9 / ONNX Runtime) segmentation, upright mask/ECT export, and 120-vertex contour extraction."""
    _logger = _main_attr("logger", logger)
    _proj_root = _main_attr("PROJECT_ROOT", PROJECT_ROOT)
    _subprocess = _main_attr("subprocess", subprocess)
    _verify_file_exists = _main_attr("verify_file_exists", verify_file_exists)
    _verify_dir_has_files = _main_attr("verify_dir_has_files", verify_dir_has_files)
    _ensure_weights = _main_attr("ensure_model_weights", ensure_model_weights)
    _check_gpu = _main_attr("check_gpu_availability", check_gpu_availability)
    _get_lm3_py = _main_attr("get_lm3_python_executable", get_lm3_python_executable)
    _resolve_ctx = _main_attr("_resolve_active_run_context", _resolve_active_run_context)
    _has_contours = _main_attr("_db_has_extracted_contours", _db_has_extracted_contours)
    _sync_cohort = _main_attr("_sync_run_cohort_in_db", _sync_run_cohort_in_db)
    _logger.info("=== Starting Phase 2: Segmentation & Geometric Leaf Extraction ===")
    start_ts = datetime.now().astimezone().isoformat()
    run_dir, cfg = _resolve_ctx(args, cfg)

    script_path = _proj_root / "scripts" / "pipeline" / "02_segment_and_extract.py"
    _verify_file_exists(script_path, "segmentation script")

    if run_dir is not None and not args.vouchers:
        vouchers_csv = ensure_run_curated_vouchers_snapshot(run_dir, workspace_root=_proj_root)
    else:
        vouchers_csv = Path(args.vouchers or cfg.get("paths", {}).get("curated_vouchers_csv", "data/tables/curated_vouchers.csv"))
    _verify_file_exists(vouchers_csv, "curated vouchers metadata table", "python main.py harvest")

    if args.weights:
        model_weights = Path(args.weights)
    else:
        seg_weights_cand = None
        if isinstance(cfg, PipelineConfig) and getattr(cfg, "segmentation", None) is not None:
            seg_weights_cand = Path(cfg.segmentation.model_weights)
        elif isinstance(cfg, dict) and cfg.get("segmentation", {}).get("model_weights"):
            cand_p = Path(cfg["segmentation"]["model_weights"])
            seg_weights_cand = cand_p if cand_p.is_absolute() else (_proj_root / cand_p)

        if seg_weights_cand is not None and not seg_weights_cand.exists():
            try:
                from scripts.vision.configure_leafmachine3 import resolve_lm3_weights
                resolved_seg = resolve_lm3_weights(seg_weights_cand.name, explicit_path=seg_weights_cand)
                if resolved_seg.exists():
                    seg_weights_cand = resolved_seg
            except Exception:
                pass

        if seg_weights_cand is not None and seg_weights_cand.exists():
            model_weights = seg_weights_cand
        else:
            try:
                model_weights = _ensure_weights(cfg)
            except Exception as exc:
                _logger.warning(f"Automated weights acquisition could not complete: {exc}")
                default_seg_weights = (
                    "LeafMachine3/leafmachine3/models/"
                    "lm3_leaf_segmenter__yolo26x_seg_1024.onnx"
                )
                model_weights = Path(
                    cfg.get("paths", {}).get("model_weights", default_seg_weights)
                )
    _verify_file_exists(model_weights, "LeafMachine3 model weights checkpoint")

    device = args.device or cfg.get("segmentation", {}).get("device", "cuda")
    if device == "cuda":
        _check_gpu(warn_only=True)

    lm3_python = _get_lm3_py()
    try:
        display_path = lm3_python.relative_to(_proj_root).as_posix()
    except ValueError:
        display_path = str(lm3_python)
    _logger.info(f"Executing LeafMachine3 segmentation via: {display_path}")

    default_out_dir = (run_dir / "data") if run_dir is not None else (_proj_root / "data")
    out_dir = Path(args.output_dir or default_out_dir)
    cmd = [
        str(lm3_python), str(script_path),
        "--vouchers", str(vouchers_csv),
        "--model-weights", str(model_weights),
        "--output-dir", str(out_dir),
        "--device", device,
        "--score-thresh", str(args.score_thresh or cfg["segmentation"]["score_thresh"]),
    ]
    if args.limit:
        cmd.extend(["--limit", str(args.limit)])
    if getattr(args, "threads", None) is not None:
        cmd.extend(["--threads", str(args.threads)])
    if getattr(args, "gpu_batch_size", None) is not None:
        cmd.extend(["--gpu-batch-size", str(args.gpu_batch_size)])
    if getattr(args, "wave_size", None) is not None:
        cmd.extend(["--wave-size", str(args.wave_size)])
    if getattr(args, "prefetch_workers", None) is not None:
        cmd.extend(["--prefetch-workers", str(args.prefetch_workers)])
    if getattr(args, "prep_workers", None) is not None:
        cmd.extend(["--prep-workers", str(args.prep_workers)])
    if getattr(args, "post_workers", None) is not None:
        cmd.extend(["--post-workers", str(args.post_workers)])
    if getattr(args, "purge_segments", False):
        cmd.append("--purge-segments")
    if getattr(args, "force", False):
        cmd.append("--force")

    _logger.info(f"Running command: {' '.join(cmd)}")
    try:
        _subprocess.run(cmd, check=True)
    except _subprocess.CalledProcessError as exc:
        _logger.error(f"LeafMachine3 segmentation failed with exit code {exc.returncode}")
        sys.exit(exc.returncode)

    contours_dir = (run_dir / "data" / "contours") if run_dir is not None else Path(cfg["paths"]["contours_dir"])
    if not _has_contours(contours_dir):
        _verify_dir_has_files(contours_dir, "*.csv", "extracted leaf contours")
    if run_dir is not None:
        _sync_cohort(run_dir, args=args, cfg=cfg)
        write_run_metadata(run_dir, start_time=start_ts)
    _logger.info("Phase 2 completed successfully.")
