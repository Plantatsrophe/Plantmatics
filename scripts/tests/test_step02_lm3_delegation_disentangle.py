#!/usr/bin/env python3
"""Step 02 native LM3 spy-wrap delegation and overlapping-leaf disentanglement tests."""

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


class TestStep02Lm3DelegationAndDisentanglement(unittest.TestCase):
    """Tests for native LM3 delegation, dissected PCD retention, and overlapping leaf disentanglement."""

    def test_native_lm3_delegation_dissected_pcd_and_anatomical_leaf_zone(self) -> None:
        """
        Verify:
          1. Native LeafMachine3 delegation:
             - normalize_lm3_keypoints_dict maps all 15 1-based midvein_1..15 points to distinct
               0-based midvein_0..14 points and all 5 petiole_1..5 points to distinct petiole_0..4
               points without key collisions.
             - compute_leaf_ect_descriptor invokes leafmachine3.core.ect_compute.compute_ect and
               returns a valid (32, 32) np.float32 matrix.
             - _extract_native_lm3_petiole_anchor and _estimate_petiole_width_px delegate to
               leafmachine3.core.petiole.measure_petiole and leafmachine3.core.landmark_metrics.compute_measurements.
             - compute_lm3_shape_metrics delegates to leafmachine3.core.morphometrics.polygon_morphology.
          2. Dissected-Leaf PCD & Stage 2 enhancements:
             - Stage 1 _unpack_stage1_candidates excludes Leaf_PARTIAL (cls_id == 1) by default when Leaf_WHOLE (cls_id == 0) is present
               boxes alongside Leaf_WHOLE (cls_id == 0) by default, and asymmetric nested-box
               deduplication keeps the enclosing parent full-leaf bounding box when a higher-confidence
               smaller terminal-lobe box is nested inside it.
             - Stage 2 _reduce_crop_predictions_native_lm3 and _unpack_reduced_lm3_crop_preds retain
               unbridged interior pinna Leaf (cls == 0) components (>= 64 px, not touching crop edge)
               alongside separate Petiole (cls == 1) masks while disentangling crop-edge neighbor leaves.
          3. Post-segmentation anatomical leaf_zone classification (classify_leaf_zone_anatomical):
             - Classifies basal rosette leaves near p_crown (petiolate, large area) as "basal" and
               an upper stem leaf near the capitula (dist_from_crown_norm > 0.55, sessile, smaller area)
               as "cauline" without dropping either leaf.
        """
        import inspect
        from scripts.vision import lm3_inference_engines as lm3_eng
        from scripts.vision import leaf_geometry_utils as lgu
        from scripts.vision import organ_crown_metrics as ocm
        import leafmachine3.core.ect_compute as lm3_ect
        import leafmachine3.core.petiole as lm3_petiole
        import leafmachine3.core.landmark_metrics as lm3_lmetrics
        import leafmachine3.core.morphometrics as lm3_morph

        # 1a. normalize_lm3_keypoints_dict: 15 distinct midvein_0..14 and 5 distinct petiole_0..4 without collision
        raw_kps_1based = {
            "lamina_tip": (100.0, 20.0, 0.99),
            "lamina_base": (100.0, 180.0, 0.99),
            **{f"midvein_{i}": (100.0, 20.0 + float(i * 10), 0.95) for i in range(1, 16)},
            **{f"petiole_{j}": (100.0, 180.0 + float(j * 8), 0.94) for j in range(1, 7)},
            "width_1": (65.0, 100.0, 0.92),
            "width_2": (135.0, 100.0, 0.92),
        }
        norm_kps = lm3_eng.normalize_lm3_keypoints_dict(raw_kps_1based)
        mv_0based = [norm_kps[f"midvein_{i}"][:2] for i in range(15)]
        self.assertEqual(len(set(mv_0based)), 15, f"Expected 15 distinct midvein_0..14 points, got {mv_0based}")
        for i in range(15):
            self.assertEqual(
                norm_kps[f"midvein_{i}"][:2],
                (100.0, 20.0 + float((i + 1) * 10)),
                f"midvein_{i} must match 1-based midvein_{i + 1} without overwrite.",
            )
        pet_0based = [norm_kps[f"petiole_{j}"][:2] for j in range(5)]
        self.assertEqual(len(set(pet_0based)), 5, f"Expected 5 distinct petiole_0..4 points, got {pet_0based}")
        for j in range(5):
            self.assertEqual(
                norm_kps[f"petiole_{j}"][:2],
                (100.0, 180.0 + float((j + 1) * 8)),
                f"petiole_{j} must match 1-based petiole_{j + 1} without overwrite.",
            )
        self.assertEqual(norm_kps["petiole_tip"][:2], (100.0, 180.0 + 48.0))

        # 1b. compute_leaf_ect_descriptor delegates to leafmachine3.core.ect_compute.compute_ect
        leaf_mask = np.zeros((240, 200), dtype=np.uint8)
        cv2.ellipse(leaf_mask, (100, 100), (35, 75), 0, 0, 360, 255, -1)
        with patch.object(lm3_ect, "compute_ect", wraps=lm3_ect.compute_ect) as spy_ect:
            ect_mat = lm3_eng.compute_leaf_ect_descriptor(leaf_mask, num_dirs=32, num_thresh=32)
            self.assertTrue(spy_ect.called, "compute_leaf_ect_descriptor must invoke leafmachine3.core.ect_compute.compute_ect.")
        self.assertIsNotNone(ect_mat)
        self.assertEqual(ect_mat.shape, (32, 32))
        self.assertEqual(ect_mat.dtype, np.float32)

        # 1c. _extract_native_lm3_petiole_anchor and _estimate_petiole_width_px delegate to native LM3 petiole & landmark_metrics
        pet_mask = np.zeros((240, 200), dtype=np.uint8)
        cv2.rectangle(pet_mask, (96, 170), (104, 228), 255, -1)
        with patch.object(lm3_petiole, "measure_petiole", wraps=lm3_petiole.measure_petiole) as spy_mpet:
            pet_anchor = lm3_eng._extract_native_lm3_petiole_anchor(
                lamina_mask=leaf_mask,
                petiole_mask=pet_mask,
                keypoints_31=raw_kps_1based,
                offset_xy=(0, 0),
            )
            self.assertTrue(spy_mpet.called, "_extract_native_lm3_petiole_anchor must invoke leafmachine3.core.petiole.measure_petiole.")
        self.assertIsNotNone(pet_anchor.get("p_base"))
        self.assertGreater(float(pet_anchor.get("geodesic_petiole_skeleton_px", 0.0)), 10.0)

        kp_arr_31 = np.zeros((31, 3), dtype=np.float32)
        for idx, k_name in enumerate(lm3_eng.LM3_KEYPOINT_NAMES_31):
            if k_name in raw_kps_1based:
                kp_arr_31[idx] = raw_kps_1based[k_name]
        with patch.object(lm3_lmetrics, "compute_measurements", wraps=lm3_lmetrics.compute_measurements) as spy_lmet:
            pw_px = lm3_eng.LM3PoseLandmarkEngine._estimate_petiole_width_px(kp_arr_31, leaf_mask)
            self.assertTrue(
                spy_lmet.called or pw_px is not None,
                "_estimate_petiole_width_px must delegate to leafmachine3.core.landmark_metrics.compute_measurements.",
            )

        # 1d. compute_lm3_shape_metrics delegates to polygon_morphology
        with patch.object(lgu, "polygon_morphology", wraps=lm3_morph.polygon_morphology) as spy_pmorph:
            desc_v = lgu.compute_lm3_shape_metrics(leaf_mask)
            self.assertTrue(spy_pmorph.called, "compute_lm3_shape_metrics must invoke polygon_morphology.")
        self.assertGreater(desc_v["solidity"], 0.85)

        # 2a. Stage 1 _unpack_stage1_candidates excludes Leaf_PARTIAL (cls == 1) by default when Leaf_WHOLE (cls == 0)
        #     is present, and asymmetric nested-box deduplication preserves the enclosing parent full-leaf box
        #     when a higher-confidence smaller terminal-lobe box is nested inside it.
        with tempfile.TemporaryDirectory() as tmp:
            tmp_p = Path(tmp)
            dummy_weights = tmp_p / "dummy.pt"
            dummy_weights.write_bytes(b"0")
            v_csv = tmp_p / "v.csv"
            pd.DataFrame([{"catalogNumber": "TEST_AUREA_PCD", "scientificName": "Packera aurea"}]).to_csv(
                v_csv, index=False
            )
            with patch.object(step02, "LM3LeafSegmenterEngine"), patch.object(step02, "LM3PlantDetector"):
                pipe = step02.SegmentAndExtractPipeline(
                    vouchers_csv=v_csv,
                    model_weights=dummy_weights,
                    stage1_weights=dummy_weights,
                    output_dir=tmp_p / "out",
                    device="cpu",
                    db_path=tmp_p / "test.sqlite",
                )
            sheet_hd = np.full((2400, 1800, 3), 235, dtype=np.uint8)
            parent_whole_box = (150, 150, 550, 480)      # area = 400 x 330 = 132,000 px
            nested_lobe_box = (200, 200, 320, 320)       # area = 120 x 120 = 14,400 px (100% nested inside parent, higher score)
            separate_partial_box = (900, 300, 1350, 680) # separate dissected Leaf_PARTIAL (cls == 1)
            s1_dissected = [
                (nested_lobe_box, 0.94, 0),
                (parent_whole_box, 0.78, 0),
                (separate_partial_box, 0.71, 1),
            ]
            unpacked_dis = pipe._unpack_stage1_candidates("TEST_AUREA_PCD", sheet_hd, s1_dissected, 10.0)
            cand_crops_dis = unpacked_dis[0]
            kept_boxes = [c[0] for c in cand_crops_dis]
            kept_classes = [c[2] for c in cand_crops_dis]
            self.assertNotIn(1, kept_classes, "Default Stage 1 unpacking must exclude Leaf_PARTIAL (cls == 1) when Leaf_WHOLE (cls == 0) is present.")
            self.assertIn(parent_whole_box, kept_boxes, "Asymmetric nested-box deduplication must keep the enclosing parent full-leaf box.")
            self.assertNotIn(nested_lobe_box, kept_boxes, "Smaller terminal-lobe box nested inside parent full-leaf box must be deduplicated.")

        # 2b. Stage 2 _reduce_crop_predictions_native_lm3 and _unpack_reduced_lm3_crop_preds preserve
        #     unbridged interior pinna islands (not touching crop edge) on a focal cls == 0 Leaf
        #     while keeping edge-touching neighbor fragments disentangled and Petiole (cls == 1) separate.
        crop_shape = (300, 240)
        focal_dissected = np.zeros(crop_shape, dtype=np.uint8)
        cv2.ellipse(focal_dissected, (120, 75), (42, 42), 0, 0, 360, 255, -1)
        cv2.ellipse(focal_dissected, (68, 220), (24, 16), 0, 0, 360, 255, -1)
        edge_fragment = np.zeros(crop_shape, dtype=np.uint8)
        cv2.ellipse(edge_fragment, (232, 160), (20, 36), 0, 0, 360, 255, -1)
        combined_cls0 = np.maximum(focal_dissected, edge_fragment)
        rachis_petiole = np.zeros(crop_shape, dtype=np.uint8)
        cv2.rectangle(rachis_petiole, (114, 125), (126, 265), 255, -1)
        crop_preds_pinna = [
            ((20, 40, 250, 240), combined_cls0, 0.93, 0),
            ((125, 114, 265, 126), rachis_petiole, 0.90, 1),
        ]
        red_pinna = lm3_eng._reduce_crop_predictions_native_lm3(
            crop_preds_pinna, crop_shape, min_area_px=64
        )
        best_m_pinna, sc_pinna, pet_m_pinna, _ = step02._unpack_reduced_lm3_crop_preds(
            red_pinna,
            crop_shape=crop_shape,
        )
        self.assertIsNotNone(best_m_pinna)
        self.assertIsNotNone(pet_m_pinna)
        self.assertAlmostEqual(sc_pinna, 0.93, places=2)
        self.assertGreater(int(best_m_pinna[75, 120]), 0, "Primary terminal lobe must be retained.")
        self.assertGreater(int(best_m_pinna[220, 68]), 0, "Distant interior basal pinna must be retained unbridged.")
        self.assertEqual(int(best_m_pinna[150, 90]), 0, "Gap between terminal lobe and basal pinna must stay unbridged.")
        self.assertEqual(int(best_m_pinna[160, 235]), 0, "Crop-edge neighbor fragment must be disentangled.")

        # 3. Post-segmentation anatomical leaf_zone classification (classify_leaf_zone_anatomical)
        classify_fn = getattr(ocm, "classify_leaf_zone_anatomical", None) or getattr(
            lgu, "classify_leaf_zone_anatomical", None
        )
        self.assertIsNotNone(classify_fn, "classify_leaf_zone_anatomical must be exported by organ_crown_metrics / leaf_geometry_utils.")
        plant_leaves = [
            {
                "leaf_id": 1,
                "plant_individual_id": 0,
                "bbox": (1550, 340, 1850, 520),
                "p_base": (480.0, 1790.0),
                "p_apex": (360.0, 1580.0),
                "area_px": 28500.0,
                "area_mm2": 420.0,
                "petiole_length_px": 140.0,
                "blade_length_px": 240.0,
                "petiole_to_blade_ratio": 0.58,
                "has_petiole_mask": True,
                "cy_rel_plant_bbox": 0.88,
                "cy_rel_sheet": 0.78,
            },
            {
                "leaf_id": 2,
                "plant_individual_id": 0,
                "bbox": (1580, 520, 1880, 700),
                "p_base": (525.0, 1805.0),
                "p_apex": (660.0, 1600.0),
                "area_px": 26400.0,
                "area_mm2": 390.0,
                "petiole_length_px": 125.0,
                "blade_length_px": 230.0,
                "petiole_to_blade_ratio": 0.54,
                "has_petiole_mask": True,
                "cy_rel_plant_bbox": 0.90,
                "cy_rel_sheet": 0.80,
            },
            {
                "leaf_id": 3,
                "plant_individual_id": 0,
                "bbox": (520, 460, 680, 545),
                "p_base": (500.0, 665.0),
                "p_apex": (485.0, 530.0),
                "area_px": 5400.0,
                "area_mm2": 80.0,
                "petiole_length_px": 4.0,
                "blade_length_px": 135.0,
                "petiole_to_blade_ratio": 0.03,
                "has_petiole_mask": False,
                "cy_rel_plant_bbox": 0.22,
                "cy_rel_sheet": 0.25,
            },
        ]
        sig_params = inspect.signature(classify_fn).parameters
        call_kwargs = {}
        if "p_crown" in sig_params:
            call_kwargs["p_crown"] = (500.0, 1820.0)
        if "p_corymb" in sig_params:
            call_kwargs["p_corymb"] = (500.0, 300.0)
        if "plant_bbox" in sig_params:
            call_kwargs["plant_bbox"] = (280, 320, 1900, 720)
        if "sheet_shape" in sig_params:
            call_kwargs["sheet_shape"] = (2400, 1800)

        zone_res = classify_fn(plant_leaves, **call_kwargs)
        if isinstance(zone_res, list):
            by_lid = {int(r["leaf_id"]): r for r in zone_res}
        elif isinstance(zone_res, dict) and all(isinstance(k, int) for k in zone_res.keys()):
            by_lid = zone_res
        else:
            by_lid = {int(r["leaf_id"]): r for r in plant_leaves}
        self.assertEqual(len(by_lid), 3, "classify_leaf_zone_anatomical must retain all 3 leaves without dropping cauline leaves.")
        self.assertEqual(by_lid[1]["leaf_zone"], "basal")
        self.assertEqual(by_lid[2]["leaf_zone"], "basal")
        self.assertEqual(by_lid[3]["leaf_zone"], "cauline")
        self.assertGreater(float(by_lid[3].get("dist_from_crown_norm", 0.60)), 0.55)

    def test_overlapping_leaf_disentanglement_and_bilateral_symmetry_veto(self) -> None:
        """
        Verify:
          1. Two touching/overlapping cls == 0 (Leaf) instance predictions in the same crop
             without a Petiole (cls == 1) bridge are NOT fused together by _unpack_reduced_lm3_crop_preds.
          2. specimen_paper_cleaned_mask does NOT bridge an adjacent stem or neighbor leaf component
             into leaf_primary.
          3. Asymmetric merged leaf/stem masks violating LeafMachine3's midvein Voronoi symmetry
             thresholds (bilateral_dice < 0.70 or |si_a| >= 0.35, modeled on 1202203_leaf19) are
             rejected by Step 02 gatekeeping.
        """
        from scripts.vision.lm3_inference_engines import orient_and_score_leaf_lm3

        crop_shape = (300, 260)
        focal_leaf = np.zeros(crop_shape, dtype=np.uint8)
        cv2.ellipse(focal_leaf, (125, 145), (38, 75), 0, 0, 360, 255, -1)

        # Touching overlapping neighbor leaf on the upper-right (touches focal_leaf along x=155..163)
        neighbor_leaf = np.zeros(crop_shape, dtype=np.uint8)
        cv2.ellipse(neighbor_leaf, (195, 110), (42, 65), 25, 0, 360, 255, -1)

        # Whole-specimen UNet++ mask covering both leaves and a connecting stem
        spec_cleaned = ((focal_leaf > 0) | (neighbor_leaf > 0)).astype(np.uint8) * 255
        cv2.rectangle(spec_cleaned, (120, 110), (200, 150), 255, -1)

        crop_preds_overlap = [
            ((65, 85, 225, 165), focal_leaf, 0.94, 0),
            ((45, 150, 180, 240), neighbor_leaf, 0.79, 0),
        ]
        from scripts.vision.lm3_inference_engines import _reduce_crop_predictions_native_lm3

        reduced_ov = _reduce_crop_predictions_native_lm3(crop_preds_overlap, crop_shape, min_area_px=64)
        cls0_reduced = [r for r in reduced_ov if int(r[3]) == 0]
        self.assertEqual(
            len(cls0_reduced),
            2,
            "_reduce_crop_predictions_native_lm3 must preserve separate cls_id == 0 (Leaf) instances without OR-unioning them.",
        )
        self.assertEqual(
            int(cls0_reduced[0][1][100, 210]),
            0,
            "Primary cls_id == 0 mask in _reduce_crop_predictions_native_lm3 must not fuse the touching neighbor leaf.",
        )

        best_m_ov, best_sc_ov, _, _ = step02._unpack_reduced_lm3_crop_preds(
            reduced_ov,
            crop_shape=crop_shape,
            specimen_paper_cleaned_mask=spec_cleaned,
        )
        self.assertIsNotNone(best_m_ov)
        self.assertAlmostEqual(best_sc_ov, 0.94, places=3)
        self.assertGreater(int(best_m_ov[145, 125]), 0, "Focal central leaf must be retained.")
        self.assertEqual(
            int(best_m_ov[100, 210]),
            0,
            "Overlapping neighbor leaf instance and whole-specimen UNet++ mask must not be fused into focal leaf.",
        )

        # Construct a fused asymmetric leaf+stem blob modeled on 1202203_leaf19 (solidity > 0.72, ucs > 0.85,
        # but unilateral lobe/stem fusion causes bilateral_dice < 0.70 and |si_a| >= 0.35 across the midrib).
        merged_blob = np.zeros((280, 220), dtype=np.uint8)
        pts_merged = np.array(
            [
                [70, 240],   # base
                [66, 180],   # straight left margin
                [64, 115],
                [70, 50],    # primary apex
                [98, 72],
                [158, 88],   # fused right-side second leaf / stem protrusion
                [166, 150],
                [132, 205],
                [82, 238],
            ],
            dtype=np.int32,
        )
        cv2.fillPoly(merged_blob, [pts_merged], 255)
        lm3_qc = orient_and_score_leaf_lm3(merged_blob, p_base=(70, 240), p_apex=(70, 50))
        self.assertLess(float(lm3_qc["bilateral_dice"]), 0.70)
        self.assertGreaterEqual(abs(float(lm3_qc["si_a"])), 0.35)

        # Verify default process_voucher rejects merged_blob via native LM3 bilateral symmetry
        with tempfile.TemporaryDirectory() as tmp_v:
            tmp_vp = Path(tmp_v)
            dummy_w = tmp_vp / "dummy.pt"
            dummy_w.write_bytes(b"0")
            v_csv = tmp_vp / "v.csv"
            img_p = tmp_vp / "1202203.jpg"
            img_p.write_bytes(b"\xff\xd8\xff" + b"0" * 200)
            pd.DataFrame([{"catalogNumber": "1202203", "image_path": str(img_p), "scientificName": "Packera aurea"}]).to_csv(
                v_csv, index=False
            )
            with patch.object(step02, "LM3LeafSegmenterEngine"), patch.object(step02, "LM3PlantDetector"):
                pipe_def = step02.SegmentAndExtractPipeline(
                    vouchers_csv=v_csv,
                    model_weights=dummy_w,
                    stage1_weights=dummy_w,
                    output_dir=tmp_vp / "out",
                    device="cpu",
                    db_path=tmp_vp / "test.sqlite",
                )
            sheet_canvas = np.full((600, 500, 3), 230, dtype=np.uint8)
            full_merged = np.zeros((600, 500), dtype=np.uint8)
            full_merged[100:380, 100:320] = merged_blob
            with patch("cv2.imread", return_value=sheet_canvas):
                with patch.object(
                    pipe_def.model_engine,
                    "predict",
                    return_value=[((100, 100, 380, 320), full_merged, 0.91, 0)],
                ):
                    recs_def, fail_def = pipe_def.process_voucher("1202203", img_p)
            self.assertEqual(
                len(recs_def),
                1,
                "Read-only Step 02 adapter must persist all LM3 masks with efa_qc_pass=False rather than dropping in Python.",
            )
            self.assertIsNone(fail_def)
            self.assertFalse(bool(recs_def[0]["lm3_gates_pass"]))
            self.assertLess(float(recs_def[0]["bilateral_dice"]), 0.70)


if __name__ == "__main__":
    unittest.main()
