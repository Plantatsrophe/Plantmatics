#!/usr/bin/env python3
"""
Unit and integration tests for Step 02 Stage 1 to Stage 2 crop reprojection,
reproductive routing, multi-plant clustering, isotropic contour aspect ratio,
high-DPI blade filtering, and Stage 1 cache handling.
"""

from __future__ import annotations

import importlib
import shutil
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import cv2
import numpy as np
import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

step02 = importlib.import_module("scripts.pipeline.02_segment_and_extract")
SegmentAndExtractPipeline = step02.SegmentAndExtractPipeline
export_standardized_contour = step02.export_standardized_contour


class TestRoutingAndManifestSchema(unittest.TestCase):
    """Test suite for the 5-step routing decision tree and diagnostic manifest schema."""

    def setUp(self) -> None:
        self.temp_dir = Path(tempfile.mkdtemp(prefix="test_routing_"))
        self.vouchers_csv = self.temp_dir / "vouchers.csv"
        self.output_dir = self.temp_dir / "data"
        self.img_path = self.temp_dir / "NCU123.jpg"
        self.img_path.write_bytes(b"\xff\xd8\xff" + b"0" * 200)

        pd.DataFrame([
            {"catalogNumber": "NCU123", "image_path": str(self.img_path), "scientificName": "Packera paupercula"}
        ]).to_csv(self.vouchers_csv, index=False)

        self.pipeline = SegmentAndExtractPipeline(
            vouchers_csv=self.vouchers_csv,
            model_weights=Path("models/dummy.pth"),
            output_dir=self.output_dir,
            device="cpu",
            force=True,
        )

    def tearDown(self) -> None:
        self.pipeline.shutdown()
        if self.temp_dir.exists():
            shutil.rmtree(self.temp_dir, ignore_errors=True)

    def test_stage1_to_stage2_crop_reprojection(self) -> None:
        """Verify Stage 1 (LM3PlantDetector 11-class) -> 12%-padded crop -> Stage 2 LM3LeafSegmenterEngine reprojection."""
        from unittest.mock import patch

        sheet_img = np.full((1000, 800, 3), 220, dtype=np.uint8)
        # Stage 1 box for leaf_whole (class 0): (ymin=200, xmin=100, ymax=400, xmax=300) -> 200x200 box
        # 12% padding -> pad_h=24, pad_w=24 -> crop_box = (176, 76, 424, 324) -> size (248, 248)
        s1_bbox = (200, 100, 400, 300)
        crop_bgr, crop_box = step02.crop_leaf_bbox_with_padding(sheet_img, s1_bbox, pad_ratio=0.12)
        self.assertEqual(crop_box, (176, 76, 424, 324))
        self.assertEqual(crop_bgr.shape[:2], (248, 248))

        # Build a pristine elliptic mask inside the 248x248 crop
        crop_mask = np.zeros((248, 248), dtype=np.uint8)
        cv2.ellipse(crop_mask, (124, 124), (55, 95), 0, 0, 360, 255, -1)

        sheet_mask, full_bbox = step02.reproject_crop_mask_to_sheet(crop_mask, crop_box, (1000, 800))
        self.assertEqual(sheet_mask.shape, (1000, 800))
        self.assertEqual(np.count_nonzero(sheet_mask), np.count_nonzero(crop_mask))
        # Verify reprojected mask is positioned at crop_box offset
        np.testing.assert_array_equal(sheet_mask[176:424, 76:324], crop_mask)

        # Now test full two-stage pipeline execution with mixed 11-class Stage 1 detections:
        # 0: leaf_whole, 1: leaf_partial, 2: leaflet (ignored), 8: specimen (ignored), 10: wood (ignored)
        stage1_dets = [
            ((200, 100, 400, 300), 0.94, 0),  # leaf_whole -> routed to Stage 2
            ((550, 400, 750, 600), 0.89, 1),  # leaf_partial -> routed to Stage 2 when segment_partial_leaves=True
            ((100, 50, 150, 90), 0.80, 2),    # leaflet -> ignored
            ((10, 10, 990, 790), 0.99, 8),    # specimen -> ignored
            ((800, 100, 950, 140), 0.77, 10), # wood -> ignored
        ]

        def fake_stage2_predict(crop_img: np.ndarray):
            ch, cw = crop_img.shape[:2]
            m = np.zeros((ch, cw), dtype=np.uint8)
            cv2.ellipse(m, (cw // 2, ch // 2), (int(cw * 0.24), int(ch * 0.40)), 0, 0, 360, 255, -1)
            return [((10, 10, ch - 10, cw - 10), m, 0.96, 0)]

        prev_seg_partial = getattr(self.pipeline.config.segmentation, "segment_partial_leaves", False)
        object.__setattr__(self.pipeline.config.segmentation, "segment_partial_leaves", True)
        try:
            with patch("cv2.imread", return_value=sheet_img):
                with patch.object(self.pipeline.stage1_detector, "predict", return_value=stage1_dets):
                    with patch.object(self.pipeline.model_engine, "predict", side_effect=fake_stage2_predict) as mock_s2:
                        records, failure = self.pipeline.process_voucher("NCU123", self.img_path)
        finally:
            object.__setattr__(self.pipeline.config.segmentation, "segment_partial_leaves", prev_seg_partial)

        self.assertIsNone(failure)
        # Stage 2 LM3LeafSegmenterEngine must be called ONLY for class 0 (leaf_whole) and class 1 (leaf_partial)
        self.assertEqual(mock_s2.call_count, 2)
        self.assertEqual(len(records), 2)
        for rec in records:
            self.assertTrue(bool(rec["lm3_gates_pass"]))
            self.assertTrue(Path(rec["contour_path"]).exists())
            self.assertTrue(Path(rec["mask_path"]).exists())

    def test_reproductive_routing_and_capitula_patch_saving(self) -> None:
        """Verify Stage 1 reproductive detections (classes 2..6) route to phenology classifier & save capitula crops."""
        from unittest.mock import patch

        sheet_img = np.full((1000, 800, 3), 230, dtype=np.uint8)
        # Paint golden-yellow corollas inside the Flower_ONE (class 4) region (50..170, 300..400)
        # BGR for golden yellow (RGB ~ 240, 200, 20 -> BGR = 20, 200, 240)
        sheet_img[60:120, 310:390] = (20, 200, 240)
        # Paint dark green erect involucre at lower half (120..165, 325..375)
        sheet_img[120:165, 325:375] = (30, 90, 40)

        stage1_dets = [
            ((300, 200, 550, 420), 0.95, 0),  # 0: Leaf_WHOLE
            ((50, 300, 170, 400), 0.92, 4),   # 4: Flower_ONE (H=120, W=100, AR=1.2)
            ((60, 450, 170, 550), 0.88, 2),   # 2: Seed_Fruit_ONE (H=110, W=100, AR=1.1)
            ((70, 600, 160, 680), 0.84, 6),   # 6: Bud (H=90, W=80, AR=1.125)
        ]

        def fake_stage2_predict(crop_img: np.ndarray):
            ch, cw = crop_img.shape[:2]
            m = np.zeros((ch, cw), dtype=np.uint8)
            cv2.ellipse(m, (cw // 2, ch // 2), (int(cw * 0.25), int(ch * 0.40)), 0, 0, 360, 255, -1)
            return [((10, 10, ch - 10, cw - 10), m, 0.95, 0)]

        with patch("cv2.imread", return_value=sheet_img):
            with patch.object(self.pipeline.stage1_detector, "predict", return_value=stage1_dets):
                with patch.object(self.pipeline.model_engine, "predict", side_effect=fake_stage2_predict):
                    extracted_df, failed_df = self.pipeline.run()

        self.assertEqual(len(extracted_df), 1)
        self.assertEqual(len(failed_df), 0)

        # Default save_debug_crops=False skips writing debug capitula JPEGs while preserving phenology tally
        cap_patches_default = sorted(self.pipeline.capitula_patches_dir.glob("NCU123_cap*.jpg"))
        self.assertEqual(cap_patches_default, [], "Debug capitula crops must be skipped when save_debug_crops=False.")

        # Verify voucher_phenological_states.csv was written and contains total_capitula == 3
        pheno_csv = self.pipeline.tables_dir / "voucher_phenological_states.csv"
        self.assertTrue(pheno_csv.exists())
        pheno_df = pd.read_csv(pheno_csv)
        self.assertEqual(len(pheno_df), 1)
        self.assertEqual(int(pheno_df.iloc[0]["total_capitula"]), 3)
        self.assertIn(
            pheno_df.iloc[0]["voucher_phenological_state"],
            {"anthesis", "fruit", "bud"},
        )

        # Enabling save_debug_crops=True exports all 3 capitula patch JPEGs
        self.pipeline.save_debug_crops = True
        self.pipeline.force = True
        with patch("cv2.imread", return_value=sheet_img):
            with patch.object(self.pipeline.stage1_detector, "predict", return_value=stage1_dets):
                with patch.object(self.pipeline.model_engine, "predict", side_effect=fake_stage2_predict):
                    extracted_dbg, failed_dbg = self.pipeline.run()

        self.assertEqual(len(extracted_dbg), 1)
        self.assertEqual(len(failed_dbg), 0)
        cap_patches = sorted(self.pipeline.capitula_patches_dir.glob("NCU123_cap*.jpg"))
        self.assertEqual(len(cap_patches), 3, "All 3 reproductive detections must be saved as capitula crops when save_debug_crops=True.")

    def test_multi_plant_clustering_and_tier1_mask_generation(self) -> None:
        """Verify multi-plant DBSCAN clustering and Tier 1 binary mask saving per plant_individual_id."""
        from unittest.mock import patch

        sheet_img = np.full((2000, 2000, 3), 210, dtype=np.uint8)
        # Two spatially separated plants on the same sheet:
        # Plant 0 near top-left (100..300, 100..280) and Plant 1 near bottom-right (1500..1700, 1500..1680)
        stage1_dets = [
            ((100, 100, 320, 280), 0.95, 0),
            ((120, 140, 340, 310), 0.93, 0),
            ((1500, 1500, 1720, 1680), 0.91, 0),
        ]

        def fake_stage2_predict(crop_img: np.ndarray):
            ch, cw = crop_img.shape[:2]
            m = np.zeros((ch, cw), dtype=np.uint8)
            cv2.ellipse(m, (cw // 2, ch // 2), (int(cw * 0.25), int(ch * 0.40)), 0, 0, 360, 255, -1)
            return [((10, 10, ch - 10, cw - 10), m, 0.95, 0)]

        with patch("cv2.imread", return_value=sheet_img):
            with patch.object(self.pipeline.stage1_detector, "predict", return_value=stage1_dets):
                with patch.object(self.pipeline.model_engine, "predict", side_effect=fake_stage2_predict):
                    records, failure = self.pipeline.process_voucher("NCU123", self.img_path)

        self.assertIsNone(failure)
        self.assertEqual(len(records), 3)

        # Verify 2 distinct plant_individual_id clusters were formed
        plant_ids = sorted({r["plant_individual_id"] for r in records})
        self.assertEqual(plant_ids, [0, 1])

        # Verify binary leaf masks exist directly in data/masks/
        for leaf_id in [1, 2, 3]:
            leaf_mask = self.pipeline.masks_dir / f"NCU123_leaf{leaf_id}.png"
            self.assertTrue(leaf_mask.exists(), f"Binary leaf mask {leaf_mask} must exist.")

    def test_export_standardized_contour_isotropic_aspect_ratio_preservation(self) -> None:
        """
        Verify an oblong tip-up leaf with L/W = 3.0 exported via export_standardized_contour()
        retains (max(y) - min(y)) / (max(x) - min(x)) == 3.0 within 1e-3.
        """
        mask = np.zeros((500, 500), dtype=np.uint8)
        # Oblong tip-up leaf with semi-axes a=150 (L=300 along y), b=50 (W=100 along x) -> L/W = 3.0
        cv2.ellipse(mask, (250, 250), (50, 150), 0.0, 0, 360, 255, -1)

        csv_path = export_standardized_contour(
            mask,
            catalog_number="TEST_OBLONG_LW3",
            leaf_id=1,
            output_dir=self.temp_dir,
            num_points=120,
            petiole_attachment_pt=(250, 400),
        )
        self.assertIsNotNone(csv_path)
        df = pd.read_csv(csv_path)

        span_x = float(df["x"].max() - df["x"].min())
        span_y = float(df["y"].max() - df["y"].min())
        lw_ratio = max(span_x, span_y) / max(min(span_x, span_y), 1e-9)
        self.assertAlmostEqual(
            lw_ratio,
            3.0,
            delta=1e-3,
            msg=f"Isotropically normalized contour must preserve L/W = 3.0 within 1e-3 (got {lw_ratio:.6f}).",
        )

    def test_stage1_high_dpi_filtering_and_single_primary_blade_selection(self) -> None:
        """Verify on high-DPI herbarium sheets (>2000px), Stage 1 forwards all Class 0 (`Leaf_WHOLE`) boxes to Stage 2 while skipping Class 1 (`Leaf_PARTIAL`)."""
        sheet_img = np.full((3000, 2400, 3), 230, dtype=np.uint8)
        self.assertTrue(max(sheet_img.shape[:2]) > 2000)

        stage1_dets = [
            ((1500, 500, 1900, 900), 0.90, 1),   # leaf_partial -> skipped
            ((400, 800, 800, 1200), 0.92, 0),    # upper-sheet Class 0 leaf_whole -> forwarded
            ((1600, 600, 1680, 680), 0.88, 0),   # 80x80 Class 0 box -> forwarded to Stage 2
            ((1800, 800, 2400, 1400), 0.95, 0),  # lower-sheet Class 0 leaf_whole -> forwarded
        ]

        def fake_high_dpi_stage2(crop_img: np.ndarray):
            ch, cw = crop_img.shape[:2]
            # Prediction 1: Class 1 (petiole/stem fragment)
            m1 = np.zeros((ch, cw), dtype=np.uint8)
            cv2.rectangle(m1, (10, 10), (min(cw - 1, 30), min(ch - 1, 200)), 255, -1)

            # Prediction 2: Class 0, tiny noise fragment
            m2 = np.zeros((ch, cw), dtype=np.uint8)
            cv2.circle(m2, ( min(cw - 1, 50), min(ch - 1, 50)), 5, 255, -1)

            # Prediction 3: Class 0, primary blade in crop
            m3 = np.zeros((ch, cw), dtype=np.uint8)
            cv2.ellipse(m3, (cw // 2, ch // 2), (max(8, cw // 3), max(8, ch // 3)), 0, 0, 360, 255, -1)

            # Prediction 4: Class 0, smaller blade fragment (lower score)
            m4 = np.zeros((ch, cw), dtype=np.uint8)
            cv2.ellipse(m4, (cw // 2, ch // 2), (max(5, cw // 5), max(5, ch // 5)), 0, 0, 360, 255, -1)

            return [
                ((10, 10, min(ch, 200), min(cw, 30)), m1, 0.90, 1),
                ((10, 10, min(ch, 30), min(cw, 30)), m2, 0.85, 0),
                ((10, 10, max(20, ch - 10), max(20, cw - 10)), m3, 0.95, 0),
                ((20, 20, max(30, ch - 20), max(30, cw - 20)), m4, 0.70, 0),
            ]

        prev_seg_partial = getattr(self.pipeline.config.segmentation, "segment_partial_leaves", False)
        object.__setattr__(self.pipeline.config.segmentation, "segment_partial_leaves", False)
        try:
            with patch("cv2.imread", return_value=sheet_img):
                with patch.object(self.pipeline.stage1_detector, "predict", return_value=stage1_dets):
                    with patch.object(self.pipeline.model_engine, "predict", side_effect=fake_high_dpi_stage2) as mock_s2:
                        records, failure = self.pipeline.process_voucher("NCU_HIGHDPI", self.img_path)
        finally:
            object.__setattr__(self.pipeline.config.segmentation, "segment_partial_leaves", prev_seg_partial)

        self.assertIsNone(failure)
        # Stage 2 LM3LeafSegmenterEngine is called for all 3 Class 0 (`Leaf_WHOLE`) boxes (Class 1 skipped)
        self.assertEqual(mock_s2.call_count, 3)
        self.assertEqual(len(records), 3)
        self.assertTrue(all(bool(r["lm3_gates_pass"]) for r in records))

    def test_stage1_use_cache_flag(self) -> None:
        """Verify LM3PlantDetector.predict respects use_cache flag."""
        with tempfile.TemporaryDirectory() as tmp_dir:
            proj_dir = Path(tmp_dir)
            labels_dir = proj_dir / "Plant_Components" / "labels"
            labels_dir.mkdir(parents=True)
            label_file = labels_dir / "VOUCHER1.txt"
            label_file.write_text("0 0.5 0.5 0.2 0.2\n")

            detector = step02.LM3PlantDetector(
                weights_path="dummy.pt",
                lm3_project_dir=proj_dir,
            )
            # When use_cache=True, loads precomputed labels
            dets_cached = detector.predict(np.zeros((100, 100, 3), dtype=np.uint8), catalog_number="VOUCHER1", use_cache=True)
            self.assertEqual(len(dets_cached), 1)

            # When use_cache=False, bypasses precomputed labels (and returns [] since model is None)
            dets_uncached = detector.predict(np.zeros((100, 100, 3), dtype=np.uint8), catalog_number="VOUCHER1", use_cache=False)
            self.assertEqual(len(dets_uncached), 0)


if __name__ == "__main__":
    unittest.main()
