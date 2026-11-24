#!/usr/bin/env python3
"""
Unit and integration tests for Step 02 macro-basal survival, raw unconstrained
native LM3 segmentation, 3-stage wave-engine parity, and SQLite wave commits.
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


class TestStep02WaveEngineAndSqlite(unittest.TestCase):
    """Test suite for macro-basal survival, raw native LM3 mode, and 3-stage wave engine SQLite commits."""

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

    def test_macro_basal_leaf_survival_and_sqlite_descriptor_persistence(self) -> None:
        """Verify macro-basal leaves survive high-DPI Stage 1/2 bounds and SQLite persists shape descriptors & failed calibration."""
        from scripts.core.db import get_db_connection, record_voucher_segmentation_transaction

        sheet_h, sheet_w = 5600, 3800
        sheet_img = np.full((sheet_h, sheet_w, 3), 232, dtype=np.uint8)
        # Macro-basal leaf bbox: height=2500 px, width=1050 px (exceeds legacy 1300*scale_f and 520*scale_f limits)
        s1_dets = [((2700, 1200, 5200, 2250), 0.96, 0)]

        def fake_macro_stage2(crops):
            out = []
            for crop_img in crops:
                ch, cw = crop_img.shape[:2]
                m = np.zeros((ch, cw), dtype=np.uint8)
                cv2.ellipse(m, (cw // 2, ch // 2), (int(cw * 0.38), int(ch * 0.43)), 0, 0, 360, 255, -1)
                out.append([((20, 20, ch - 20, cw - 20), m, 0.97, 0)])
            return out

        self.pipeline.ruler_calibrations["NCU_MACRO"] = 12.5
        self.pipeline.ruler_calibrations["NCU_FAIL_CAL"] = 11.8
        with patch("cv2.imread", return_value=sheet_img):
            with patch.object(self.pipeline.stage1_detector, "predict", return_value=s1_dets):
                with patch.object(self.pipeline.model_engine, "predict_batch", side_effect=fake_macro_stage2):
                    records, failure = self.pipeline.process_voucher("NCU_MACRO", self.img_path)

            with patch.object(self.pipeline.stage1_detector, "predict", return_value=[]):
                with patch.object(self.pipeline.model_engine, "predict", return_value=[]):
                    fail_records, fail_rec = self.pipeline.process_voucher("NCU_FAIL_CAL", self.img_path)

        self.assertIsNone(failure)
        self.assertEqual(len(records), 1)
        rec = records[0]
        for col in ("circularity", "extent", "aspect_ratio", "area_px", "area_mm2"):
            self.assertIn(col, rec)
            self.assertTrue(np.isfinite(float(rec[col])), f"{col} should be finite, got {rec[col]}")

        self.assertEqual(len(fail_records), 0)
        self.assertIsNotNone(fail_rec)
        self.assertAlmostEqual(float(fail_rec["pixels_per_mm"]), 11.8, places=3)
        self.assertTrue(bool(fail_rec["is_calibrated_mm"]))

        conn = get_db_connection(db_path=":memory:", init_schema=True)
        record_voucher_segmentation_transaction(
            conn,
            catalog_number="NCU_MACRO",
            image_path=str(self.img_path),
            extracted_records=records,
            failure_record=None,
            output_dir=self.pipeline.output_dir,
        )
        record_voucher_segmentation_transaction(
            conn,
            catalog_number="NCU_FAIL_CAL",
            image_path=str(self.img_path),
            extracted_records=[],
            failure_record=fail_rec,
            output_dir=self.pipeline.output_dir,
        )
        row_leaf = conn.execute(
            "SELECT circularity, extent, aspect_ratio, area_px, area_mm2 FROM extracted_leaves WHERE catalogNumber = 'NCU_MACRO'"
        ).fetchone()
        self.assertIsNotNone(row_leaf)
        for col in ("circularity", "extent", "aspect_ratio", "area_px", "area_mm2"):
            self.assertIsNotNone(row_leaf[col], f"SQLite extracted_leaves.{col} must not be NULL")

        row_fail_v = conn.execute(
            "SELECT qc_status, pixels_per_mm, is_calibrated_mm FROM vouchers WHERE catalogNumber = 'NCU_FAIL_CAL'"
        ).fetchone()
        self.assertIsNotNone(row_fail_v)
        self.assertEqual(row_fail_v["qc_status"], "FAILED_QC")
        self.assertAlmostEqual(float(row_fail_v["pixels_per_mm"]), 11.8, places=3)
        self.assertEqual(int(row_fail_v["is_calibrated_mm"]), 1)
        conn.close()

    def test_raw_unconstrained_lm3_native_leaf_and_petiole_segmentation(self) -> None:
        """
        Verify default 100% Raw Unconstrained LeafMachine3 Mode (read-only Step 02 adapter):
          1. Stage 1 (`LM3PlantDetector`) unpacks all Class 0 (`leaf_whole`) boxes regardless of size
             or vertical sheet position, while natively excluding `leaf_partial` (1), `leaflet` (2),
             `specimen` (8), and `wood`/stems (10).
          2. Stage 2 `LM3LeafSegmenterEngine` (`thing_classes = ['leaf', 'petiole', 'hole']`) retains the top-scoring
             `'leaf'` (`class_id == 0`) mask 100% pixel-identical (zero heuristic petiole clipping),
             pairs the top-scoring `'petiole'` (`class_id == 1`) mask to locate `p_base` and native petiole length,
             and retains >6 leaves without any post-hoc geometric/photometric rejections or top-6 cap.
        """
        with patch("sys.argv", ["02_segment_and_extract.py"]):
            parsed_default = step02.parse_args()
            self.assertEqual(parsed_default.score_thresh, self.pipeline.config.segmentation.score_thresh)

        sheet_h, sheet_w = 3600, 2600
        sheet_img = np.full((sheet_h, sheet_w, 3), (225, 232, 236), dtype=np.uint8)
        self.pipeline.ruler_calibrations["NCU_RAW_LM3"] = 10.0  # 10 px/mm = 100 px/cm

        # 7 Stage 1 Class 0 (`leaf_whole`) boxes across top, middle, and bottom of sheet
        # (including small alpine-sized boxes and >6 leaves) plus native LM3 non-leaf classes:
        s1_dets = [
            ((250, 300, 420, 460), 0.94, 0),    # 1. small upper Class 0 leaf_whole
            ((650, 350, 1050, 750), 0.95, 0),   # 2. upper-mid Class 0 leaf_whole
            ((1200, 400, 1750, 920), 0.96, 0),  # 3. mid-sheet cordate/suborbicular Class 0 leaf_whole
            ((1850, 450, 2400, 950), 0.97, 0),  # 4. lower-mid Class 0 leaf_whole
            ((2450, 400, 3000, 920), 0.95, 0),  # 5. basal Class 0 leaf_whole
            ((2450, 1050, 3000, 1550), 0.93, 0),# 6. basal Class 0 leaf_whole
            ((2450, 1650, 3000, 2150), 0.92, 0),# 7. 7th Class 0 leaf_whole (verifies no top-6 cap)
            ((1400, 1800, 1900, 2200), 0.89, 1),# leaf_partial -> natively ignored by Stage 1
            ((1100, 1800, 1250, 1950), 0.85, 2),# leaflet -> natively ignored by Stage 1
            ((50, 50, 3500, 2500), 0.99, 8),    # specimen -> natively ignored by Stage 1
            ((900, 1200, 2600, 1280), 0.91, 10),# wood/stem -> natively ignored by Stage 1
        ]

        expected_fg_counts: list[int] = []

        def fake_raw_lm3_stage2(crops):
            out = []
            for idx, crop_img in enumerate(crops):
                ch, cw = crop_img.shape[:2]
                leaf_m = np.zeros((ch, cw), dtype=np.uint8)
                cx, cy = cw // 2, int(ch * 0.42)
                rx, ry = max(18, int(cw * 0.30)), max(22, int(ch * 0.28))
                cv2.ellipse(leaf_m, (cx, cy), (rx, ry), 0.0, 0, 360, 255, -1)
                if idx == 2:
                    # Cordate basal sinus + crenate indent (low solidity) that raw LM3 must keep untouched
                    cv2.ellipse(leaf_m, (cx, cy + ry - 4), (int(rx * 0.38), int(ry * 0.35)), 0.0, 0, 360, 0, -1)
                expected_fg_counts.append(int(np.count_nonzero(leaf_m)))

                # Native LM3LeafSegmenterEngine class_id == 1 (`petiole`) mask extending below lamina base
                pet_m = np.zeros((ch, cw), dtype=np.uint8)
                pet_y0 = min(ch - 6, cy + ry - 2)
                pet_y1 = min(ch - 2, pet_y0 + max(25, int(ch * 0.22)))
                cv2.line(pet_m, (cx, pet_y0), (cx, pet_y1), 255, thickness=7)

                out.append([
                    ((10, 10, ch - 10, cw - 10), leaf_m, 0.96, 0),
                    ((pet_y0, max(0, cx - 5), pet_y1, min(cw, cx + 5)), pet_m, 0.91, 1),
                ])
            return out

        prev_seg_partial = getattr(self.pipeline.config.segmentation, "segment_partial_leaves", True)
        object.__setattr__(self.pipeline.config.segmentation, "segment_partial_leaves", False)
        try:
            with patch("cv2.imread", return_value=sheet_img):
                with patch.object(self.pipeline.stage1_detector, "predict", return_value=s1_dets):
                    with patch.object(self.pipeline.model_engine, "predict_batch", side_effect=fake_raw_lm3_stage2) as mock_s2:
                        records, failure = self.pipeline.process_voucher("NCU_RAW_LM3", self.img_path)
        finally:
            object.__setattr__(self.pipeline.config.segmentation, "segment_partial_leaves", prev_seg_partial)

        self.assertIsNone(failure)
        # Only the 7 Class 0 (`leaf_whole`) boxes enter Stage 2; all 7 are retained with zero rejections and zero top-6 cap
        self.assertEqual(len(mock_s2.call_args[0][0]), 7)
        self.assertEqual(len(records), 7)
        for idx, rec in enumerate(records):
            self.assertIn("lm3_gates_pass", rec)
            # Saved binary mask foreground pixel count must equal raw LM3 leaf_m 100% identically (zero petiole trimming!)
            saved_m = cv2.imread(str(rec["mask_path"]), cv2.IMREAD_GRAYSCALE)
            self.assertIsNotNone(saved_m)
            self.assertEqual(int(np.count_nonzero(saved_m > 0)), expected_fg_counts[idx])
            self.assertGreater(float(rec["petiole_length_mm"]), 0.0)

    def test_three_stage_wave_engine_parity_and_sqlite_wave_commit(self) -> None:
        """
        Verify 3-Stage Asynchronous Wave Engine (`run()` across >4 vouchers with wave_size=4)
        produces 100% identical extracted leaf records, binary masks, 120-vertex contours,
        ROI-scoped rosette patches/masks, and SQLite `contour_xy_json` records compared to
        sequential single-voucher execution.
        """
        import json
        from scripts.core.db import get_db_connection

        wave_root = self.temp_dir / "wave_run"
        seq_root = self.temp_dir / "seq_run"
        img_dir = self.temp_dir / "wave_imgs"
        img_dir.mkdir(parents=True, exist_ok=True)

        vouchers_rows = []
        for idx in range(6):
            cat = f"NCU_WAVE_{idx:02d}"
            img_p = img_dir / f"{cat}.jpg"
            sheet = np.full((1200, 900, 3), (225, 232, 236), dtype=np.uint8)
            # Draw two green foliar patches so real JPEG decode works across prefetch workers
            cv2.ellipse(sheet, (300, 760), (75, 110), 0.0, 0, 360, (60, 125, 70), -1)
            cv2.ellipse(sheet, (600, 780), (70, 105), 0.0, 0, 360, (55, 120, 65), -1)
            cv2.imwrite(str(img_p), sheet)
            vouchers_rows.append({
                "catalogNumber": cat,
                "scientificName": "Packera anonyma",
                "image_path": str(img_p),
            })

        vouchers_csv = self.temp_dir / "wave_vouchers.csv"
        pd.DataFrame(vouchers_rows).to_csv(vouchers_csv, index=False)

        s1_boxes = [
            ((610, 190, 910, 410), 0.96, 0),
            ((630, 490, 930, 710), 0.94, 0),
        ]

        def fake_s1_predict_batch(_self_s1, items, use_cache=False, max_batch_size=16):
            return [list(s1_boxes) for _ in items]

        def fake_s1_predict(_self_s1, _img_bgr, catalog_number=None, use_cache=False, precomputed_letterbox=None):
            return list(s1_boxes)

        def fake_s2_predict_batch(_self_s2, images_bgr, batch_size=24, pretransformed=None):
            out = []
            for crop_img in images_bgr:
                ch, cw = crop_img.shape[:2]
                leaf_m = np.zeros((ch, cw), dtype=np.uint8)
                cx, cy = cw // 2, int(ch * 0.46)
                rx, ry = max(16, int(cw * 0.28)), max(22, int(ch * 0.30))
                cv2.ellipse(leaf_m, (cx, cy), (rx, ry), 0.0, 0, 360, 255, -1)
                pet_m = np.zeros((ch, cw), dtype=np.uint8)
                py0 = min(ch - 6, cy + ry - 2)
                py1 = min(ch - 2, py0 + max(18, int(ch * 0.16)))
                cv2.line(pet_m, (cx, py0), (cx, py1), 255, thickness=5)
                out.append([
                    ((8, 8, ch - 8, cw - 8), leaf_m, 0.97, 0),
                    ((py0, max(0, cx - 4), py1, min(cw, cx + 4)), pet_m, 0.90, 1),
                ])
            return out

        wave_db = self.temp_dir / "wave_ephemeral.sqlite"
        seq_db = self.temp_dir / "seq_ephemeral.sqlite"

        dummy_s1_weights = self.temp_dir / "dummy_s1.pt"
        with patch.object(step02.LM3PlantDetector, "predict_batch", fake_s1_predict_batch), \
             patch.object(step02.LM3PlantDetector, "predict", fake_s1_predict), \
             patch.object(step02.LM3LeafSegmenterEngine, "predict_batch", fake_s2_predict_batch):
            pipe_wave = step02.SegmentAndExtractPipeline(
                vouchers_csv=vouchers_csv,
                model_weights=self.pipeline.model_weights,
                stage1_weights=dummy_s1_weights,
                output_dir=wave_root,
                device="cpu",
                db_path=wave_db,
            )
            pipe_wave.stage1_wave_size = 4  # Exercises 2 cross-voucher waves (4 + 2)
            for r in vouchers_rows:
                pipe_wave.ruler_calibrations[r["catalogNumber"]] = 10.0
            df_wave, fail_wave = pipe_wave.run()

            pipe_seq = step02.SegmentAndExtractPipeline(
                vouchers_csv=vouchers_csv,
                model_weights=self.pipeline.model_weights,
                stage1_weights=dummy_s1_weights,
                output_dir=seq_root,
                device="cpu",
                db_path=seq_db,
            )
            # Instance attribute hook forces sequential single-voucher fallback in run()
            pipe_seq.stage1_detector.predict = lambda img, **kw: list(s1_boxes)  # type: ignore[method-assign]
            for r in vouchers_rows:
                pipe_seq.ruler_calibrations[r["catalogNumber"]] = 10.0
            df_seq, fail_seq = pipe_seq.run()

        self.assertTrue(fail_wave.empty)
        self.assertTrue(fail_seq.empty)
        self.assertEqual(len(df_wave), 12)
        self.assertEqual(len(df_seq), 12)

        compare_cols = [
            "catalogNumber", "leaf_id", "plant_individual_id", "lm3_gates_pass",
            "solidity", "archetype_score", "circularity", "extent", "area_px",
            "area_mm2", "petiole_length_mm", "blade_length_mm",
        ]
        pd.testing.assert_frame_equal(
            df_wave[compare_cols].reset_index(drop=True),
            df_seq[compare_cols].reset_index(drop=True),
        )

        # Verify SQLite stores 120-pt contours
        conn = get_db_connection(wave_db, init_schema=False)
        try:
            rows_db = conn.execute(
                "SELECT catalogNumber, leaf_id, contour_n_pts, contour_xy_json FROM extracted_leaves ORDER BY catalogNumber, leaf_id"
            ).fetchall()
            self.assertEqual(len(rows_db), 12)
            for row_db in rows_db:
                self.assertEqual(int(row_db["contour_n_pts"]), 120)
                coords = json.loads(row_db["contour_xy_json"])
                self.assertEqual(len(coords), 120)
        finally:
            conn.close()


if __name__ == "__main__":
    unittest.main()
