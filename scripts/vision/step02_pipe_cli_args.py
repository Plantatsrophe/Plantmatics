#!/usr/bin/env python3
"""Step 02 pipeline runner orchestration (wave/sequential dispatch) and CLI argument helpers."""

from __future__ import annotations

import argparse
import logging
import os
import sys
from pathlib import Path
from typing import Any, Callable, Dict, List, Optional, Sequence, Tuple

import pandas as pd

from scripts.core.config import PipelineConfig
from scripts.core.db import purge_segments_and_runs, record_voucher_segmentation_transaction
from scripts.vision.lm3_inference_engines import torch_cuda_available
from scripts.vision.step02_figures import (
    generate_step02_publication_figures,
    regenerate_step02_figures_from_tables,
)
from scripts.vision.step02_pipe_manifest_flush import (
    load_resumption_state,
    retire_skipped_voucher,
    should_skip_voucher,
)
from scripts.vision.step02_pipe_wave_runner import execute_wave_engine

PROJECT_ROOT = Path(__file__).resolve().parents[2]
logger = logging.getLogger("SegmentAndExtract")


def _run_sequential_vouchers(
    pipe: Any,
    voucher_rows: Sequence[Tuple[Any, Any]],
    state: Dict[str, Any],
    manifest_path: Path,
    total_vouchers: int,
    all_extracted_records: List[Dict[str, Any]],
    all_failed_records: List[Dict[str, Any]],
    project_root: Path,
) -> Tuple[int, int]:
    """Process vouchers sequentially via `pipe.process_voucher(...)` and commit per-voucher SQLite transactions."""
    processed_count = skipped_count = 0
    db_conn = state["db_conn"]

    for _, row in voucher_rows:
        cat_num = str(row["catalogNumber"]).strip()
        raw_img_path = Path(str(row["image_path"]).strip())
        if not raw_img_path.is_absolute():
            raw_img_path = project_root / raw_img_path
        if not raw_img_path.exists():
            for ext in (".webp", ".jpg", ".jpeg", ".png"):
                cand_p = raw_img_path.with_suffix(ext)
                if cand_p.exists():
                    raw_img_path = cand_p
                    break

        if should_skip_voucher(pipe, cat_num, state, manifest_path, total_vouchers):
            skipped_count = retire_skipped_voucher(
                pipe, cat_num, state, all_extracted_records, all_failed_records,
                processed_count, skipped_count, total_vouchers,
            )
            continue

        processed_count += 1
        taxon = str(row.get("scientificName", row.get("species_raw", ""))).strip()
        try:
            extracted, failure = pipe.process_voucher(cat_num, raw_img_path, taxon=taxon)
        except Exception as exc:
            logger.warning(f"Unexpected error processing voucher {cat_num} ({raw_img_path}): {exc}")
            extracted, failure = [], {
                "catalogNumber": cat_num, "image_path": str(raw_img_path),
                "failure_reason": "SEGMENTATION_RUNTIME_ERROR",
                "details": f"Runtime exception ({type(exc).__name__}) during voucher extraction: {exc}",
            }

        if extracted:
            state["completed_on_disk"].add(cat_num)
        all_extracted_records.extend(extracted)
        if failure:
            all_failed_records.append(failure)

        if db_conn is not None:
            try:
                record_voucher_segmentation_transaction(
                    conn=db_conn, catalog_number=cat_num, extracted_records=extracted, failure_record=failure,
                    repro_record=pipe.voucher_repro_records.get(cat_num), pheno_record=pipe.voucher_pheno_records.get(cat_num),
                    repro_by_pid=pipe.voucher_repro_by_pid.get(cat_num), pheno_by_pid=pipe.voucher_pheno_by_pid.get(cat_num),
                    crowns_by_pid=pipe.voucher_crowns_by_pid.get(cat_num),
                    rosette_centroids_by_pid=pipe.voucher_rosette_centroids_by_pid.get(cat_num),
                    by_plant_petiole=pipe.voucher_petiole_by_pid.get(cat_num),
                    contour_coords_by_leaf=pipe.voucher_contour_coords_by_leaf.get(cat_num),
                    output_dir=pipe.output_dir, image_path=str(raw_img_path),
                )
            except Exception as db_tx_exc:
                logger.warning(f"Per-voucher SQLite transaction failed for {cat_num}: {db_tx_exc}")

        if processed_count % 25 == 0 or (processed_count + skipped_count) == total_vouchers:
            logger.info(
                f"Progress: [Processed: {processed_count} | Skipped: {skipped_count} | Total: {total_vouchers}] "
                f"({len(all_extracted_records)} leaves extracted, {len(all_failed_records)} failed QC)."
            )
    return processed_count, skipped_count


def run_segment_pipeline(
    pipe: Any,
    limit: Optional[int] = None,
    facade_mod: Optional[Any] = None,
    base_is_completed_fn: Optional[Callable[..., Any]] = None,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Execute full Step 02 extraction pipeline across curated vouchers with resumption support."""
    logger.info("==================================================================")
    logger.info("Starting Step 02: Direct Segmentation & Botanical Leaf Extraction")
    logger.info(f"Vouchers Metadata: {pipe.vouchers_csv}")
    logger.info(f"LeafMachine3 Model: {pipe.model_weights}")
    logger.info(f"Output Directory: {pipe.output_dir}")
    logger.info(f"Device: {pipe.device}")
    logger.info(f"Force Overwrite: {pipe.force}")
    logger.info(
        "Segmentation Mode: Read-Only LeafMachine3 Ingestion Adapter "
        "(100% LM3 masks + native QC & organismal indexing saved to SQLite; EFA filtering in Step 03)"
    )
    logger.info("==================================================================")

    if not pipe.vouchers_csv.exists():
        raise FileNotFoundError(f"Curated vouchers table not found at {pipe.vouchers_csv}")

    df_vouchers = pd.read_csv(pipe.vouchers_csv)
    logger.info(f"Loaded {len(df_vouchers)} voucher records.")
    if "catalogNumber" not in df_vouchers.columns or "image_path" not in df_vouchers.columns:
        raise KeyError("Vouchers CSV must contain 'catalogNumber' and 'image_path' columns.")
    if limit and limit > 0:
        df_vouchers = df_vouchers.head(limit)
        logger.info(f"Applied execution limit: processing first {limit} vouchers.")

    manifest_path = pipe.tables_dir / "extracted_leaf_manifest.csv"
    failed_qc_path = pipe.tables_dir / "failed_qc_vouchers.csv"
    repro_manifest_path = pipe.tables_dir / "voucher_reproductive_metrics.csv"
    pheno_manifest_path = pipe.tables_dir / "voucher_phenological_states.csv"

    state = load_resumption_state(
        pipe, df_vouchers, manifest_path, failed_qc_path, repro_manifest_path, pheno_manifest_path,
        base_is_completed_fn=base_is_completed_fn,
    )
    all_extracted_records: List[Dict[str, Any]] = []
    all_failed_records: List[Dict[str, Any]] = []
    total_vouchers = len(df_vouchers)
    voucher_rows = list(df_vouchers.iterrows())

    try:
        runner_fn = execute_wave_engine if state["use_wave_engine"] else _run_sequential_vouchers
        processed_count, skipped_count = runner_fn(
            pipe, voucher_rows, state, manifest_path, total_vouchers,
            all_extracted_records, all_failed_records, PROJECT_ROOT,
        )
    finally:
        db_conn = state.get("db_conn")
        if db_conn is not None:
            try:
                db_conn.close()
            except Exception:
                pass

    logger.info(f"Extraction Summary: [Processed: {processed_count} | Skipped: {skipped_count} | Total: {total_vouchers}]")
    extracted_df, failed_df, repro_df, pheno_df = pipe._flush_manifests(
        all_extracted_records, all_failed_records, manifest_path, failed_qc_path, repro_manifest_path, pheno_manifest_path,
    )
    logger.info("==================================================================")
    logger.info("Extraction completed!")
    logger.info(f"  Extracted Leaves: {len(extracted_df)} -> {manifest_path}")
    logger.info(f"  Failed QC Vouchers: {len(failed_df)} -> {failed_qc_path}")
    logger.info(f"  Reproductive Metrics: {len(repro_df)} -> {repro_manifest_path}")
    logger.info(f"  Phenological States: {len(pheno_df)} -> {pheno_manifest_path}")
    logger.info("==================================================================")

    if len(extracted_df) > 0:
        try:
            env_run_dir = os.environ.get("PLANTMATICS_RUN_DIR", "").strip() or os.environ.get("PACKERA_RUN_DIR", "").strip()
            if env_run_dir:
                figures_dir = Path(env_run_dir).resolve() / "outputs" / "figures"
            elif pipe.output_dir.name == "data":
                figures_dir = pipe.output_dir.parent / "outputs" / "figures"
            else:
                figures_dir = pipe.output_dir / "outputs" / "figures"
            facade = facade_mod or sys.modules.get("scripts.pipeline.02_segment_and_extract")
            fig_fn = getattr(facade, "generate_step02_publication_figures", generate_step02_publication_figures)
            fig_fn(extracted_df=extracted_df, failed_df=failed_df, repro_df=repro_df, pheno_df=pheno_df, figures_dir=figures_dir)
        except Exception as fig_exc:
            logger.warning(f"Step 02 publication figure generation skipped due to non-fatal error: {fig_exc}")

    pipe.shutdown()
    return extracted_df, failed_df


def build_step02_arg_parser(
    cfg: Optional[PipelineConfig] = None, cuda_available: Optional[bool] = None
) -> argparse.ArgumentParser:
    """Construct the canonical Step 02 CLI argument parser."""
    cfg_use = cfg or PipelineConfig.from_yaml()
    has_cuda = bool(torch_cuda_available() if cuda_available is None else cuda_available)
    parser = argparse.ArgumentParser(
        description="Direct LeafMachine3 ONNX segmentation, upright leaf/ECT extraction, and 120-vertex contour export."
    )
    parser.add_argument("--vouchers", type=Path, default=cfg_use.paths.curated_vouchers_csv, help="Path to curated vouchers CSV containing catalogNumber and image_path.")
    parser.add_argument("--model-weights", type=Path, default=cfg_use.segmentation.model_weights, help="Path to LeafMachine3 YOLO26x leaf segmenter ONNX model checkpoint.")
    parser.add_argument("--output-dir", type=Path, default=cfg_use.paths.workspace_root / "data", help="Root output directory for masks, contours, and manifest tables.")
    parser.add_argument("--lm3-dir", dest="lm3_project_dir", type=Path, default=None, help="Optional LeafMachine3 project output directory for cached stage-1 detections and ruler calibrations.")
    parser.add_argument("--device", type=str, default=cfg_use.segmentation.device if has_cuda else "cpu", choices=["cuda", "cpu"], help="Compute device for inference (default: cuda if available, else cpu).")
    parser.add_argument("--score-thresh", type=float, default=cfg_use.segmentation.score_thresh, help=f"Detection confidence score threshold (default: {cfg_use.segmentation.score_thresh}).")
    parser.add_argument("-j", "--threads", "--workers", dest="threads", type=int, default=None, help="Maximum CPU worker threads for image prefetching, tensor prep, and post-processing (default: auto-scaled for 32-core/64-thread CPUs).")
    parser.add_argument("--gpu-batch-size", dest="gpu_batch_size", type=int, default=None, help="GPU inference batch size for Stage 2 LeafMachine3 leaf crops and Stage 1 wave scaling (default: 10 on CUDA, 8 on CPU).")
    parser.add_argument("--batch-size", dest="gpu_batch_size", type=int, default=argparse.SUPPRESS, help=argparse.SUPPRESS)
    parser.add_argument("--wave-size", dest="wave_size", type=int, default=None, help="Number of herbarium sheets processed per cross-voucher GPU wave (default: 16 or --gpu-batch-size).")
    parser.add_argument("--prefetch-workers", dest="prefetch_workers", type=int, default=None, help="Optional override for Stage A disk image prefetch worker thread count.")
    parser.add_argument("--prep-workers", dest="prep_workers", type=int, default=None, help="Optional override for Stage B2/B3 CPU crop prep and prototype mask decode worker thread count.")
    parser.add_argument("--post-workers", dest="post_workers", type=int, default=None, help="Optional override for Stage C contour/ECT/rosette post-processing worker thread count.")
    parser.add_argument("--purge-segments", "--purge-old-data", dest="purge_segments", action="store_true", default=False, help="Purge all previous Step 02 to 07 segments, masks, contours, run directories, and SQLite segmentation tables before running.")
    parser.add_argument("--force", dest="force", action="store_true", default=False, help="Force re-segmentation and overwrite existing contour CSVs and masks.")
    parser.add_argument("--overwrite", dest="force", action="store_true", default=argparse.SUPPRESS, help=argparse.SUPPRESS)
    parser.add_argument("--limit", type=int, default=None, help="Optional limit on number of vouchers to process.")
    parser.add_argument("--db", type=Path, default=None, help="Optional SQLite database path override (default: data/plantmatics.sqlite).")
    parser.add_argument("--regenerate-figures-only", dest="regenerate_figures_only", action="store_true", default=False, help="Reload Step 02 to 04 tables and regenerate post-QC Step 02 publication figures without running segmentation.")
    parser.add_argument("--save-debug-crops", dest="save_debug_crops", action="store_true", default=bool(getattr(cfg_use.segmentation, "save_debug_crops", False)), help="Save per-capitulum RGB debug crop JPEGs under data/cropped_patches/capitula/.")
    return parser


def parse_step02_args(argv: Optional[Sequence[str]] = None, facade_mod: Optional[Any] = None) -> argparse.Namespace:
    """Parse Step 02 CLI arguments using the facade's `PipelineConfig` and `torch_cuda_available`."""
    facade = facade_mod or sys.modules.get("scripts.pipeline.02_segment_and_extract")
    cfg_cls = getattr(facade, "PipelineConfig", PipelineConfig)
    cuda_fn = getattr(facade, "torch_cuda_available", torch_cuda_available)
    return build_step02_arg_parser(cfg=cfg_cls.from_yaml(), cuda_available=bool(cuda_fn())).parse_args(argv)


def run_step02_cli_main(facade_mod: Optional[Any] = None) -> None:
    """Execute the Step 02 CLI entrypoint via facade-resolved hooks."""
    facade = facade_mod or sys.modules.get("scripts.pipeline.02_segment_and_extract")
    parse_fn = getattr(facade, "parse_args", parse_step02_args)
    purge_fn = getattr(facade, "purge_segments_and_runs", purge_segments_and_runs)
    regen_fn = getattr(facade, "regenerate_step02_figures_from_tables", regenerate_step02_figures_from_tables)
    pipe_cls = getattr(facade, "SegmentAndExtractPipeline", None)
    if pipe_cls is None:
        import importlib

        pipe_cls = importlib.import_module("scripts.pipeline.02_segment_and_extract").SegmentAndExtractPipeline

    args = parse_fn()
    if getattr(args, "purge_segments", False):
        purge_fn(
            workspace_root=PROJECT_ROOT, db_path=getattr(args, "db", None),
            output_dir=getattr(args, "output_dir", None), purge_runs=False, vacuum=True,
        )
    if getattr(args, "regenerate_figures_only", False):
        regen_fn(output_dir=args.output_dir, vouchers_csv=args.vouchers)
        return
    pipeline = pipe_cls(
        vouchers_csv=args.vouchers, model_weights=args.model_weights, output_dir=args.output_dir,
        device=args.device, score_thresh=args.score_thresh,
        force=bool(args.force or getattr(args, "purge_segments", False)),
        lm3_project_dir=getattr(args, "lm3_project_dir", None), db_path=getattr(args, "db", None),
        threads=getattr(args, "threads", None), gpu_batch_size=getattr(args, "gpu_batch_size", None),
        wave_size=getattr(args, "wave_size", None), prefetch_workers=getattr(args, "prefetch_workers", None),
        prep_workers=getattr(args, "prep_workers", None), post_workers=getattr(args, "post_workers", None),
        save_debug_crops=bool(getattr(args, "save_debug_crops", False)),
    )
    pipeline.run(limit=args.limit)
