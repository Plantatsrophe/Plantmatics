#!/usr/bin/env python3
"""Step 02 Stage 1 rescue, working-frame unscaling, 6-pt pose, letterbox, and CLI minibatch tests."""

from __future__ import annotations

import importlib
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


class TestDissectedLeafSegmentationEnhancements(unittest.TestCase):
    """Tests for high-DPI dissected leaf geometry, rachis bridging, and petiole trimming."""

    def test_stage1_rescues_leaf_partial_when_leaf_whole_absent(self) -> None:
        """When Stage 1 finds 0 Class 0 (leaf_whole) boxes, Class 1 (leaf_partial) boxes are rescued for LM3LeafSegmenterEngine."""
        with tempfile.TemporaryDirectory() as tmp:
            tmp_p = Path(tmp)
            dummy_weights = tmp_p / "dummy.pt"
            dummy_weights.write_bytes(b"0")
            v_csv = tmp_p / "v.csv"
            pd.DataFrame([{"catalogNumber": "TEST_GLABELLA_1", "scientificName": "Packera glabella"}]).to_csv(
                v_csv, index=False
            )
            with patch.object(step02, "LM3LeafSegmenterEngine"), patch.object(
                step02, "LM3PlantDetector"
            ):
                pipe = step02.SegmentAndExtractPipeline(
                    vouchers_csv=v_csv,
                    model_weights=dummy_weights,
                    stage1_weights=dummy_weights,
                    output_dir=tmp_p / "out",
                    device="cpu",
                    db_path=tmp_p / "test.sqlite",
                )
            sheet_img = np.zeros((2000, 1500, 3), dtype=np.uint8)
            # Case A: Only Class 1 (leaf_partial) detected on sheet -> rescued into candidate_s1_crops
            s1_only_partial = [
                ((100, 200, 600, 550), 0.72, 1),
            ]
            unpacked_a = pipe._unpack_stage1_candidates("TEST_GLABELLA_1", sheet_img, s1_only_partial, 10.0)
            candidate_crops_a = unpacked_a[0]
            self.assertEqual(len(candidate_crops_a), 1)
            self.assertEqual(candidate_crops_a[0][2], 1)

            # Case B: Both Class 0 and Class 1 detected on sheet
            # With default segment_partial_leaves=True, both Class 0 and Class 1 are included;
            # when segment_partial_leaves=False is explicitly set, only Class 0 is included.
            s1_mixed = [
                ((100, 200, 600, 550), 0.85, 0),
                ((700, 400, 1000, 600), 0.65, 1),
            ]
            unpacked_b_default = pipe._unpack_stage1_candidates("TEST_GLABELLA_1", sheet_img, s1_mixed, 10.0)
            candidate_crops_b_default = unpacked_b_default[0]
            if getattr(pipe.config.segmentation, "segment_partial_leaves", False):
                self.assertEqual(len(candidate_crops_b_default), 2)
                self.assertEqual({c[2] for c in candidate_crops_b_default}, {0, 1})

            object.__setattr__(pipe.config.segmentation, "segment_partial_leaves", False)
            unpacked_b = pipe._unpack_stage1_candidates("TEST_GLABELLA_1", sheet_img, s1_mixed, 10.0)
            candidate_crops_b = unpacked_b[0]
            self.assertEqual(len(candidate_crops_b), 1)
            self.assertEqual(candidate_crops_b[0][2], 0)

    def test_lm3_working_frame_coordinate_unscaling(self) -> None:
        """
        Verify LeafMachine3 working-frame coordinate unscaling (`compute_lm3_work_scale`
        and `unscale_working_bbox_to_sheet`) converts from a 3,200 px working frame back
        to native sheet dimensions (3753x5634).
        """
        from scripts.vision.lm3_data_loader import (
            compute_lm3_work_scale,
            unscale_working_bbox_to_sheet,
        )

        sheet_h, sheet_w = 5634, 3753
        work_scale = compute_lm3_work_scale(sheet_h, sheet_w, max_working_dim=3200)
        expected_scale = 3200.0 / 5634.0
        self.assertAlmostEqual(work_scale, expected_scale, places=6)

        # Sheets already <= 3200 px on longest side must retain work_scale == 1.0
        self.assertAlmostEqual(compute_lm3_work_scale(2400, 1800, max_working_dim=3200), 1.0, places=6)

        # Working-frame bbox -> native sheet bbox (ymin, xmin, ymax, xmax)
        work_bbox = (1600, 1000, 2400, 1500)
        sheet_bbox = unscale_working_bbox_to_sheet(
            work_bbox,
            work_scale,
            sheet_shape=(sheet_h, sheet_w),
        )
        self.assertEqual(len(sheet_bbox), 4)
        ymin, xmin, ymax, xmax = sheet_bbox
        self.assertAlmostEqual(ymin, int(round(1600 / expected_scale)), delta=1)
        self.assertAlmostEqual(xmin, int(round(1000 / expected_scale)), delta=1)
        self.assertAlmostEqual(ymax, int(round(2400 / expected_scale)), delta=1)
        self.assertAlmostEqual(xmax, int(round(1500 / expected_scale)), delta=1)
        self.assertGreaterEqual(ymin, 0)
        self.assertGreaterEqual(xmin, 0)
        self.assertLessEqual(ymax, sheet_h)
        self.assertLessEqual(xmax, sheet_w)

    def test_lm3_unetpp_specimen_mask_rachis_bridging_and_6pt_petiole_pose(self) -> None:
        """
        Verify LM3 6-point petiole pose skeleton (`petiole_1..6`, `lamina_base`) handling in
        `estimate_rosette_crown_point` and `approximate_petiole_and_scape_metrics`.
        """
        sheet_h, sheet_w = 1000, 800
        sheet_img = np.full((sheet_h, sheet_w, 3), 230, dtype=np.uint8)
        true_crown_xy = (400.0, 820.0)

        leaf_records = []
        for idx, angle_deg in enumerate([-25.0, 0.0, 25.0], start=1):
            rad = np.deg2rad(angle_deg)
            dx = np.sin(rad)
            dy = -np.cos(rad)
            # Lamina base 90 px outward from crown; tip 220 px outward from crown
            bx = true_crown_xy[0] + 90.0 * dx
            by = true_crown_xy[1] + 90.0 * dy
            tx = true_crown_xy[0] + 220.0 * dx
            ty = true_crown_xy[1] + 220.0 * dy

            kp31 = np.zeros((31, 3), dtype=np.float32)
            kp31[0] = [tx, ty, 0.96]  # lamina_tip
            kp31[1] = [bx, by, 0.96]  # lamina_base
            for p_i in range(6):
                # petiole_1 (near lamina_base) -> petiole_6 (near rosette crown)
                frac = (p_i + 1) / 6.0
                px = bx + frac * (true_crown_xy[0] - bx)
                py = by + frac * (true_crown_xy[1] - by)
                kp31[2 + p_i] = [px, py, 0.94]

            m_sheet = np.zeros((sheet_h, sheet_w), dtype=np.uint8)
            cx_i, cy_i = int(round((bx + tx) * 0.5)), int(round((by + ty) * 0.5))
            cv2.circle(m_sheet, (cx_i, cy_i), 35, 255, -1)

            leaf_records.append({
                "catalogNumber": "NCU_POSE31",
                "plant_individual_id": 0,
                "leaf_id": idx,
                "bbox": (cy_i - 40, cx_i - 40, cy_i + 40, cx_i + 40),
                "mask": m_sheet,
                "p_base": (int(round(bx)), int(round(by))),
                "p_apex": (int(round(tx)), int(round(ty))),
                "keypoints_31": kp31,
                "pixels_per_mm": 10.0,
                "is_calibrated_mm": True,
            })

        crown_pt = step02.estimate_rosette_crown_point(
            leaf_records,
            sheet_shape=(sheet_h, sheet_w),
        )
        self.assertIsNotNone(crown_pt)
        dist_crown = float(np.hypot(crown_pt[0] - true_crown_xy[0], crown_pt[1] - true_crown_xy[1]))
        self.assertLessEqual(
            dist_crown,
            35.0,
            f"Estimated crown {crown_pt} with 6-pt petiole pose skeleton should converge near {true_crown_xy}.",
        )

        pet_metrics = step02.approximate_petiole_and_scape_metrics(
            leaf_records=leaf_records,
            reproductive_boxes=[],
            sheet_bgr=sheet_img,
            pixels_per_mm=10.0,
        )
        self.assertIn("petiole_length_px", pet_metrics)
        self.assertGreater(float(pet_metrics["petiole_length_px"]), 50.0)
        self.assertGreater(float(pet_metrics["petiole_length_mm"]), 5.0)

    def test_canonical_tip_up_contour_ect_and_scale_calibrator(self) -> None:
        """
        Verify canonical tip-up contour rotation in export_standardized_contour,
        LM3ScaleCalibrator megapixel prior, orient_and_score_leaf_lm3, and
        compute_leaf_ect_descriptor.
        """
        from scripts.vision.leaf_geometry_utils import export_standardized_contour
        from scripts.vision.lm3_inference_engines import (
            LM3ScaleCalibrator,
            compute_leaf_ect_descriptor,
            orient_and_score_leaf_lm3,
        )

        # 2. Canonical tip-up contour rotation on an oblique 135-degree leaf
        oblique_mask = np.zeros((400, 400), dtype=np.uint8)
        cv2.ellipse(oblique_mask, (200, 200), (95, 38), 45.0, 0, 360, 255, -1)
        kp31 = np.zeros((31, 3), dtype=np.float32)
        # Base at top-left (133, 133), tip at bottom-right (267, 267) -> pointing down-right (+135 deg from tip-up)
        kp31[0] = [267.0, 267.0, 0.98]  # lamina_tip
        kp31[1] = [133.0, 133.0, 0.98]  # lamina_base
        for i in range(23):
            frac = i / 22.0
            kp31[8 + i] = [133.0 + frac * 134.0, 133.0 + frac * 134.0, 0.95]

        with tempfile.TemporaryDirectory() as tmp:
            tmp_p = Path(tmp)
            c_path = export_standardized_contour(
                oblique_mask,
                catalog_number="TEST_OBLIQUE_135",
                leaf_id=1,
                output_dir=tmp_p,
                num_points=120,
                keypoints_31=kp31,
            )
            self.assertIsNotNone(c_path)
            cdf = pd.read_csv(c_path)
            # Point 0 (lamina_base) must sit near the bottom (max y), and point 60 (lamina_tip) near the top (min y)
            y_min_c, y_max_c = float(cdf["y"].min()), float(cdf["y"].max())
            span_y = max(y_max_c - y_min_c, 1e-6)
            span_x = float(cdf["x"].max() - cdf["x"].min())
            self.assertGreater((float(cdf.loc[0, "y"]) - y_min_c) / span_y, 0.85)
            self.assertLess((float(cdf.loc[60, "y"]) - y_min_c) / span_y, 0.15)
            self.assertGreater(span_y, 1.8 * span_x)

            # 3. Module 18 ECT descriptor export
            ect_file = tmp_p / "TEST_OBLIQUE_135_leaf1_ect.npy"
            ect_arr = compute_leaf_ect_descriptor(
                oblique_mask,
                output_path=ect_file,
                num_dirs=32,
                num_thresh=32,
            )
            self.assertIsNotNone(ect_arr)
            self.assertTrue(ect_file.exists())
            self.assertEqual(ect_arr.shape, (32, 32))
            self.assertEqual(ect_arr.dtype, np.float32)

        # 4. LM3 native leaf orientation and bilateral archetype scoring
        orient_res = orient_and_score_leaf_lm3(oblique_mask, keypoints_31=kp31)
        self.assertTrue(orient_res["success"])
        self.assertGreater(abs(float(orient_res["angle_cw"])), 100.0)
        self.assertGreaterEqual(float(orient_res["archetype_score"]), 0.70)

        # 5. LM3ScaleCalibrator megapixel prior fallback when no ruler is present
        calibrator = LM3ScaleCalibrator(device="cpu")
        calibrator.archival_session = None
        sheet_dummy = np.zeros((3000, 2000, 3), dtype=np.uint8)  # 6.0 MP
        cal_res = calibrator.calibrate_sheet(sheet_dummy, catalog_number="NCU_MP_FALLBACK")
        expected_px_per_mm = 27.112084 * np.sqrt(6.0) / 10.0
        self.assertAlmostEqual(float(cal_res["pixels_per_mm"]), expected_px_per_mm, places=3)
        self.assertEqual(cal_res["cf_source"], "predicted_from_megapixels")

    def test_native_lm3_stride32_letterbox_tight_crops_and_cordate_survival(self) -> None:
        """
        Verify Phases 1-4 native LeafMachine3 alignment:
          1. Dynamic stride-32 letterboxing (`auto=True, stride=32`) matches `ultra_rep.letterbox`
             instead of forcing a square canvas.
          2. Default `crop_leaf_bbox_with_padding` uses tight `pad_ratio=0.0` matching `detector_common.py`.
          3. `LM3LeafSegmenterEngine.segment_leaf_crop` unions multiple `Petiole` and `Hole` instances (`>= 64 px`).
          4. Stage 1 `_deduplicate_sheet_boxes` uses `overlap=0.95` (`inter / min_area >= 0.95`) so adjacent
             rosette leaves with moderate overlap (`IoU ~ 0.55`) both survive.
          5. Broad cordate/orbicular leaves (`W/L > 0.78, solidity < 0.84`) survive high-DPI QC.
        """
        from scripts.vision.lm3_inference_engines import (
            LM3LeafSegmenterEngine,
            LM3PlantDetector,
            _letterbox_bgr,
            _letterbox_gain_pad_hw,
            crop_leaf_bbox_with_padding,
        )

        # 1. Dynamic stride-32 letterboxing on a 600x200 crop at img_size=1024
        dummy_crop = np.full((600, 200, 3), 180, dtype=np.uint8)
        lb_bgr, gain, (pad_x, pad_y) = _letterbox_bgr(dummy_crop, 1024)
        self.assertEqual(lb_bgr.shape[0], 1024)
        self.assertLess(lb_bgr.shape[1], 1024)
        self.assertEqual(lb_bgr.shape[1] % 32, 0)
        gain_hw, (pad_hw_x, pad_hw_y) = _letterbox_gain_pad_hw(lb_bgr.shape[:2], dummy_crop.shape[:2])
        self.assertAlmostEqual(gain_hw, gain, places=6)
        self.assertAlmostEqual(pad_hw_x, pad_x, places=6)
        self.assertAlmostEqual(pad_hw_y, pad_y, places=6)

        # 2. Default tight crop (`pad_ratio=0.0`)
        sheet_img = np.full((1000, 800, 3), 235, dtype=np.uint8)
        crop_tight, box_tight = crop_leaf_bbox_with_padding(sheet_img, (200, 150, 500, 450))
        self.assertEqual(box_tight, (200, 150, 500, 450))
        self.assertEqual(crop_tight.shape[:2], (300, 300))

        # 3. Multi-instance Petiole and Hole unioning in LM3LeafSegmenterEngine.segment_leaf_crop
        seg_eng = LM3LeafSegmenterEngine(device="cpu")
        lam_m = np.zeros((200, 200), dtype=np.uint8)
        cv2.circle(lam_m, (100, 90), 55, 255, -1)
        pet_1 = np.zeros((200, 200), dtype=np.uint8)
        cv2.rectangle(pet_1, (96, 145), (104, 168), 255, -1)  # 9 * 24 = 216 px >= 64
        pet_2 = np.zeros((200, 200), dtype=np.uint8)
        cv2.rectangle(pet_2, (96, 165), (104, 190), 255, -1)  # overlaps pet_1 across rows 165..168
        hole_1 = np.zeros((200, 200), dtype=np.uint8)
        cv2.circle(hole_1, (85, 85), 6, 255, -1)  # ~113 px >= 64
        hole_2 = np.zeros((200, 200), dtype=np.uint8)
        cv2.circle(hole_2, (115, 95), 6, 255, -1)  # ~113 px >= 64
        fake_preds = [
            ((35, 45, 145, 155), lam_m, 0.94, 0),
            ((145, 96, 168, 104), pet_1, 0.85, 1),
            ((165, 96, 190, 104), pet_2, 0.80, 1),
            ((79, 79, 91, 91), hole_1, 0.82, 2),
            ((89, 109, 101, 121), hole_2, 0.79, 2),
        ]
        with patch.object(seg_eng, "predict", return_value=fake_preds):
            lam_out, pet_out, hole_out, sc_out = seg_eng.segment_leaf_crop(np.zeros((200, 200, 3), dtype=np.uint8))
        self.assertIsNotNone(lam_out)
        self.assertIsNotNone(pet_out)
        self.assertIsNotNone(hole_out)
        self.assertGreater(int(np.count_nonzero(pet_out)), int(np.count_nonzero(pet_1)))
        self.assertGreater(int(np.count_nonzero(hole_out)), int(np.count_nonzero(hole_1)))
        self.assertAlmostEqual(sc_out, 0.94, places=3)

        # 4. Stage 1 overlap=0.95 deduplication keeps moderately overlapping rosette leaves
        boxes_s1 = [
            ((500, 500, 900, 900), 0.95, 0),
            ((550, 550, 950, 950), 0.92, 0),  # ~62% containment (< 0.95) -> both kept!
            ((510, 510, 890, 890), 0.88, 0),  # 100% contained (>= 0.95) -> suppressed!
        ]
        det_eng = LM3PlantDetector(device="cpu")
        deduped = det_eng._deduplicate_sheet_boxes(boxes_s1)
        self.assertEqual(len(deduped), 2)

        # Confidence-ordered deduplication keeps high-confidence tight box and drops lower-confidence oversized box
        boxes_nested_conf = [
            ((480, 480, 920, 920), 0.45, 0),
            ((520, 520, 880, 880), 0.96, 0),
        ]
        deduped_conf = det_eng._deduplicate_sheet_boxes(boxes_nested_conf)
        self.assertEqual(len(deduped_conf), 1)
        self.assertAlmostEqual(deduped_conf[0][1], 0.96, places=3)
        self.assertEqual(deduped_conf[0][0], (520, 520, 880, 880))

    def test_step02_cli_compute_flags_and_back_to_back_minibatch_halving(self) -> None:
        """
        Verify CLI flags (--threads/--workers, --gpu-batch-size/--batch-size, --wave-size,
        --prefetch-workers, --prep-workers, --post-workers) and _run_onnx_minibatches_back_to_back
        automatic mini-batch halving on transient memory errors.
        """
        import main as main_cli
        from scripts.vision.lm3_inference_engines import _run_onnx_minibatches_back_to_back

        with patch(
            "sys.argv",
            [
                "02_segment_and_extract.py",
                "--workers",
                "32",
                "--batch-size",
                "16",
                "--wave-size",
                "24",
                "--prefetch-workers",
                "16",
                "--prep-workers",
                "24",
                "--post-workers",
                "32",
            ],
        ):
            parsed_s2 = step02.parse_args()
            self.assertEqual(parsed_s2.threads, 32)
            self.assertEqual(parsed_s2.gpu_batch_size, 16)
            self.assertEqual(parsed_s2.wave_size, 24)
            self.assertEqual(parsed_s2.prefetch_workers, 16)
            self.assertEqual(parsed_s2.prep_workers, 24)
            self.assertEqual(parsed_s2.post_workers, 32)

        main_parser = main_cli.build_parser()
        parsed_main = main_parser.parse_args(
            [
                "segment",
                "--workers",
                "64",
                "--gpu-batch-size",
                "16",
                "--wave-size",
                "32",
                "--prefetch-workers",
                "24",
                "--prep-workers",
                "32",
                "--post-workers",
                "48",
            ]
        )
        self.assertEqual(parsed_main.threads, 64)
        self.assertEqual(parsed_main.gpu_batch_size, 16)
        self.assertEqual(parsed_main.wave_size, 32)
        self.assertEqual(parsed_main.prefetch_workers, 24)
        self.assertEqual(parsed_main.prep_workers, 32)
        self.assertEqual(parsed_main.post_workers, 48)

        # Verify _run_onnx_minibatches_back_to_back automatically halves sub_bs when B>2 raises an OOM error
        class _FakeBFCExhaustedSession:
            def __init__(self) -> None:
                self.batch_sizes_seen: list[int] = []

            def run(self, _out_names, feed_dict):
                x = feed_dict["images"]
                b = int(x.shape[0])
                self.batch_sizes_seen.append(b)
                if b > 2:
                    raise RuntimeError("BFCArena::AllocateRawInternal Failed to allocate memory")
                return [np.zeros((b, 8, 10), dtype=np.float32)]

        fake_sess = _FakeBFCExhaustedSession()
        tensors = [np.zeros((3, 8, 8), dtype=np.float32) for _ in range(6)]
        records = _run_onnx_minibatches_back_to_back(fake_sess, "images", tensors, initial_sub_bs=8)
        # Initial attempt at B=6 fails -> halves to 3 (fails) -> halves to 1 (succeeds for all 6 items)
        total_covered = sum(end_i - start_i for start_i, end_i, _ in records)
        self.assertEqual(total_covered, 6)
        self.assertGreaterEqual(len(fake_sess.batch_sizes_seen), 3)


if __name__ == "__main__":
    unittest.main()
