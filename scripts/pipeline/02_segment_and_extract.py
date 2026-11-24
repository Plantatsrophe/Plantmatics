#!/usr/bin/env python3
"""
===============================================================================
Script: 02_segment_and_extract.py
Project: Packera dubia Species Delimitation & Morphometrics Pipeline
Affiliation: University of North Carolina at Chapel Hill Herbarium (NCU)

Description:
    Production Runner for Step 02:
    Direct voucher image segmentation and leaf extraction pipeline for
    downstream Elliptic Fourier Analysis (EFA).

    Key Capabilities:
      1. Direct Ingestion: Reads curated_vouchers.csv directly without intermediate
         staging or symlink copying.
      2. LeafMachine3 ONNX Segmentation: Direct inference using canonical LeafMachine3
         YOLO26x plant/leaf/pose engines (`models/lm3/lm3_leaf_segmenter__yolo26x_seg_1024.onnx`)
         and UNet++ v9 specimen segmentation.
      3. LeafMachine3 Scale Calibration: Resolves metric scale (`pixels_per_mm`,
         `scale_source`, `is_calibrated_mm`) from LeafMachine3's Module 7 Ruler
         Lattice Solver and optional Module 1 Megapixel Conversion Prior.
      4. In-Sheet Organ & Rosette Crown Clustering: Groups organs into plant
         individuals (`plant_individual_id`), classifies anatomical `leaf_zone`
         (`basal` vs. `cauline`), and estimates per-individual rosette crowns.
      5. Read-Only LeafMachine3 Ingestion & Native QC Persistence:
           - Persists segmented LeafMachine3 leaves with native bilateral symmetry
             (`bilateral_dice`, `si_a`, `lm3_gates_pass`), unbridged focal pinna
             component counts (`n_lamina_components`), and scale provenance.
           - Defers closed-contour EFA QC gating to Step 03 while routing sheets
             with zero extractable leaves to data/tables/failed_qc_vouchers.csv.
      6. Standardized Upright Asset Export: Saves upright base-to-tip binary
         lamina masks to data/masks/{catalogNumber}_leaf{id}.png, 120-vertex
         tip-up contours to data/contours/{catalogNumber}_leaf{id}.csv, and
         32 x 360 Euler Characteristic Transform matrices to data/ect/.

Usage:
    python scripts/pipeline/02_segment_and_extract.py \
        --vouchers data/tables/curated_vouchers.csv \
        --model-weights models/lm3/lm3_leaf_segmenter__yolo26x_seg_1024.onnx \
        --output-dir data/ \
        --device cuda
===============================================================================
"""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Set, Tuple, Union

import numpy as np
import pandas as pd

# Ensure project root is in sys.path
PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from scripts.core.config import (  # noqa: E402
    PipelineConfig,
    resolve_cpu_threads,
    resolve_gpu_batch_size,
)
from scripts.core.db import purge_segments_and_runs  # noqa: E402
from scripts.vision.capitulum_phenology_classifier import (  # noqa: E402
    classify_single_capitulum,
    classify_voucher_phenology,
)
from scripts.vision.leaf_geometry_utils import (  # noqa: E402
    approximate_petiole_and_scape_metrics,
    cluster_plant_individuals,
    compute_lm3_shape_metrics,
    estimate_rosette_crown_point,
    export_standardized_contour,
    save_leaf_mask_png,
)
from scripts.vision.lm3_compat import DBSCAN  # noqa: E402
from scripts.vision.lm3_inference_engines import (  # noqa: E402
    DetectedInstance,
    LM3LeafSegmenterEngine,
    LM3PlantDetector,
    LM3PoseLandmarkEngine,
    LM3ScaleCalibrator,
    LM3SpecimenSegmenter,
    _VoucherFinalizeBundle,
    _extract_native_lm3_petiole_anchor,
    _reduce_crop_predictions_native_lm3,
    compute_capitulum_metrics,
    compute_leaf_ect_descriptor,
    crop_leaf_bbox_with_padding,
    evaluate_leaf_damage_and_fold_qc,
    is_stock_lm3_weights,
    orient_and_score_leaf_lm3,
    reproject_crop_mask_to_sheet,
    torch_cuda_available,
)
from scripts.vision.organ_crown_metrics import classify_leaf_zone_anatomical  # noqa: E402
from scripts.vision.step02_figures import (  # noqa: E402
    _atomic_save_matplotlib_fig,
    generate_step02_publication_figures,
    regenerate_step02_figures_from_tables,
)
from scripts.vision.step02_photometry_corymb_utils import (  # noqa: E402
    estimate_sheet_paper_cie_l_ref,
)
from scripts.vision.step02_pipe_cli_args import (  # noqa: E402
    build_step02_arg_parser,
    parse_step02_args,
    run_segment_pipeline,
    run_step02_cli_main,
)
from scripts.vision.step02_pipe_init_checkpoint import (  # noqa: E402
    apply_finalize_bundle,
    atomic_save_csv,
    check_voucher_completed,
    index_completed_vouchers,
    init_segment_pipeline_state,
    purge_voucher_disk_artifacts,
)
from scripts.vision.step02_pipe_leaf_readout_b import process_voucher_impl  # noqa: E402
from scripts.vision.step02_pipe_manifest_flush import flush_pipeline_manifests  # noqa: E402
from scripts.vision.step02_pipe_stage1_candidates import (  # noqa: E402
    _resolve_scale_source,
    _unpack_reduced_lm3_crop_preds,
    unpack_stage1_candidates,
)

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
logger = logging.getLogger("SegmentAndExtract")


def _resolve_engine_namespace() -> Dict[str, Any]:
    """Resolve LeafMachine3 engine symbols from the facade module so `mock.patch` intercepts cleanly."""
    facade = sys.modules.get("scripts.pipeline.02_segment_and_extract")
    mod_globals = globals()

    def _pick(name: str) -> Any:
        if facade is not None and hasattr(facade, name):
            return getattr(facade, name)
        return mod_globals[name]

    return {
        "LM3PlantDetector": _pick("LM3PlantDetector"),
        "LM3LeafSegmenterEngine": _pick("LM3LeafSegmenterEngine"),
        "LM3PoseLandmarkEngine": _pick("LM3PoseLandmarkEngine"),
        "LM3ScaleCalibrator": _pick("LM3ScaleCalibrator"),
        "LM3SpecimenSegmenter": _pick("LM3SpecimenSegmenter"),
        "torch_cuda_available": _pick("torch_cuda_available"),
    }


class SegmentAndExtractPipeline:
    """Unified LeafMachine3 ONNX detection, specimen cleanup, leaf segmentation, pose landmarking, and extraction runner."""

    def __init__(
        self,
        vouchers_csv: Union[str, Path],
        model_weights: Union[str, Path],
        output_dir: Union[str, Path],
        device: str = "cuda",
        score_thresh: float = 0.35,
        force: bool = False,
        config: Optional[PipelineConfig] = None,
        lm3_project_dir: Optional[Union[str, Path]] = None,
        stage1_weights: Optional[Union[str, Path]] = None,
        ruler_confidence_thresh: float = 0.50,
        db_path: Optional[Union[str, Path]] = None,
        threads: Optional[int] = None,
        gpu_batch_size: Optional[int] = None,
        wave_size: Optional[int] = None,
        prefetch_workers: Optional[int] = None,
        prep_workers: Optional[int] = None,
        post_workers: Optional[int] = None,
        save_debug_crops: Optional[bool] = None,
    ) -> None:
        init_segment_pipeline_state(
            self,
            _resolve_engine_namespace(),
            vouchers_csv=vouchers_csv,
            model_weights=model_weights,
            output_dir=output_dir,
            device=device,
            score_thresh=score_thresh,
            force=force,
            config=config,
            lm3_project_dir=lm3_project_dir,
            stage1_weights=stage1_weights,
            ruler_confidence_thresh=ruler_confidence_thresh,
            db_path=db_path,
            threads=threads,
            gpu_batch_size=gpu_batch_size,
            wave_size=wave_size,
            prefetch_workers=prefetch_workers,
            prep_workers=prep_workers,
            post_workers=post_workers,
            save_debug_crops=save_debug_crops,
        )

    def _apply_finalize_bundle(self, bundle: _VoucherFinalizeBundle) -> None:
        """Merge a completed per-voucher bundle into pipeline state on the coordinator thread."""
        apply_finalize_bundle(self, bundle)

    def shutdown(self) -> None:
        """Clean up background thread pools."""
        pass

    def _index_completed_vouchers(self) -> Set[str]:
        """Build an O(1) lookup set of catalogNumbers that already have mask or contour files."""
        return index_completed_vouchers(self)

    def is_voucher_completed(self, catalog_number: str) -> bool:
        """Verify whether non-empty contour CSV or binary mask artifacts exist for `catalog_number`."""
        return check_voucher_completed(self, catalog_number)

    def _purge_voucher_disk_artifacts(self, catalog_number: str) -> None:
        """Purge prior mask, contour, ECT, and capitula crop files for `catalog_number` before `--force`."""
        purge_voucher_disk_artifacts(self, catalog_number)

    @staticmethod
    def _atomic_save_csv(df: pd.DataFrame, target_path: Path) -> None:
        """Atomically persist DataFrame to CSV via temporary file replacement."""
        atomic_save_csv(df, target_path)

    def _unpack_stage1_candidates(
        self,
        catalog_number: str,
        img_bgr: np.ndarray,
        stage1_detections: Sequence[Any],
        pixels_per_mm: Optional[float],
        cf_source: str = "measured_from_ruler",
    ) -> Tuple[
        List[Tuple[Tuple[int, int, int, int], float, int, np.ndarray, Tuple[int, int, int, int]]],
        List[Tuple[Tuple[int, int, int, int], Optional[np.ndarray], float, int]],
        List[Tuple[int, int, int, int]],
        List[Tuple[float, float]],
        bool,
        int,
        bool,
        float,
        float,
    ]:
        """Thread-safe CPU Stage 1 -> Stage 2 candidate crop extractor."""
        return unpack_stage1_candidates(
            self,
            catalog_number=catalog_number,
            img_bgr=img_bgr,
            stage1_detections=stage1_detections,
            pixels_per_mm=pixels_per_mm,
            cf_source=cf_source,
        )

    def process_voucher(
        self,
        catalog_number: str,
        image_path: Path,
        taxon: Optional[str] = None,
        preloaded_img: Optional[np.ndarray] = None,
        precomputed_s1_letterbox: Optional[np.ndarray] = None,
        precomputed_s1_detections: Optional[List[Any]] = None,
        precomputed_s2_payload: Optional[Tuple[Any, ...]] = None,
        precomputed_single_stage_detections: Optional[Union[List[Any], Tuple[List[Any], bool]]] = None,
        _return_bundle: bool = False,
        precomputed_spec_mask: Optional[np.ndarray] = None,
    ) -> Any:
        """Process a single herbarium voucher specimen sheet using LeafMachine3 ONNX engines."""
        return process_voucher_impl(
            self,
            catalog_number=catalog_number,
            image_path=image_path,
            taxon=taxon,
            preloaded_img=preloaded_img,
            precomputed_s1_letterbox=precomputed_s1_letterbox,
            precomputed_s1_detections=precomputed_s1_detections,
            precomputed_s2_payload=precomputed_s2_payload,
            precomputed_single_stage_detections=precomputed_single_stage_detections,
            _return_bundle=_return_bundle,
            precomputed_spec_mask=precomputed_spec_mask,
        )

    def _flush_manifests(
        self,
        all_extracted_records: List[Dict[str, Any]],
        all_failed_records: List[Dict[str, Any]],
        manifest_path: Path,
        failed_qc_path: Path,
        repro_manifest_path: Path,
        pheno_manifest_path: Path,
    ) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
        """Build and atomically persist all 4 canonical Step 02 output tables."""
        return flush_pipeline_manifests(
            self,
            all_extracted_records=all_extracted_records,
            all_failed_records=all_failed_records,
            manifest_path=manifest_path,
            failed_qc_path=failed_qc_path,
            repro_manifest_path=repro_manifest_path,
            pheno_manifest_path=pheno_manifest_path,
        )

    def run(self, limit: Optional[int] = None) -> Tuple[pd.DataFrame, pd.DataFrame]:
        """Execute full extraction pipeline across curated vouchers with resumption support."""
        return run_segment_pipeline(
            self,
            limit=limit,
            facade_mod=sys.modules.get(__name__),
            base_is_completed_fn=SegmentAndExtractPipeline.is_voucher_completed,
        )


Step02SegmentAndExtract = SegmentAndExtractPipeline


def parse_args() -> argparse.Namespace:
    """Parse CLI arguments for Step 02 direct segmentation and botanical leaf extraction."""
    return parse_step02_args(facade_mod=sys.modules.get(__name__))


def main() -> None:
    """Run the Step 02 CLI entrypoint."""
    run_step02_cli_main(facade_mod=sys.modules.get(__name__))


__all__ = [
    "PROJECT_ROOT",
    "PipelineConfig",
    "DBSCAN",
    "logger",
    "DetectedInstance",
    "LM3PlantDetector",
    "LM3LeafSegmenterEngine",
    "LM3PoseLandmarkEngine",
    "LM3ScaleCalibrator",
    "LM3SpecimenSegmenter",
    "_VoucherFinalizeBundle",
    "_unpack_reduced_lm3_crop_preds",
    "_resolve_scale_source",
    "cluster_plant_individuals",
    "compute_lm3_shape_metrics",
    "estimate_rosette_crown_point",
    "approximate_petiole_and_scape_metrics",
    "classify_leaf_zone_anatomical",
    "classify_single_capitulum",
    "classify_voucher_phenology",
    "export_standardized_contour",
    "save_leaf_mask_png",
    "compute_capitulum_metrics",
    "compute_leaf_ect_descriptor",
    "crop_leaf_bbox_with_padding",
    "reproject_crop_mask_to_sheet",
    "evaluate_leaf_damage_and_fold_qc",
    "is_stock_lm3_weights",
    "orient_and_score_leaf_lm3",
    "estimate_sheet_paper_cie_l_ref",
    "torch_cuda_available",
    "generate_step02_publication_figures",
    "regenerate_step02_figures_from_tables",
    "_atomic_save_matplotlib_fig",
    "purge_segments_and_runs",
    "SegmentAndExtractPipeline",
    "Step02SegmentAndExtract",
    "build_step02_arg_parser",
    "parse_args",
    "main",
]

if __name__ == "__main__":
    main()
