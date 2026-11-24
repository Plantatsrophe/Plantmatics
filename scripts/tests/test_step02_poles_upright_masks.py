#!/usr/bin/env python3
"""Step 02 organ-pole ray tracing, upright mask persistence, reniform L/W, and 32x360 ECT tests."""

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
SegmentAndExtractPipeline = step02.SegmentAndExtractPipeline
compute_lm3_shape_metrics = step02.compute_lm3_shape_metrics


class TestStep02PolesAndUprightMasks(unittest.TestCase):
    """Tests for 4-case organ-pole vascular ray tracing and upright mask/ECT persistence."""

    def test_four_organ_pole_cases_and_vascular_ray_tracing(self) -> None:
        """
        Verify classify_leaf_zone_anatomical across all 4 biological organ-pole cases:
          - Case 1 (Roots + Capitula): Inverted mount (NCU00385199: roots at top y=197,
            capitula at bottom y=1572) preserves top root crown without overwriting and
            classifies top rosette leaf as basal and lower stem leaves as cauline.
          - Case 2 (Roots Only, No Capitula): Erect long-petiolate basal tuft (NCU00438219:
            roots at bottom y=2027, no capitula, elevated lamina blades at y=1580 with
            downward petiole rays tracing into root crown) classifies all tuft leaves as basal.
          - Case 3 (Capitula Only, No Roots): Clipped upper flowering stem (capitula at top,
            no roots, spaced sessile stem leaves without basal convergence hub) classifies
            internodal leaves as cauline.
          - Case 4 (Neither Roots Nor Capitula): Multi-plant sheet (NCSC00103660: pid=0
            flowering stem on right, pid=1 detached sterile basal rosette on upper-left)
            does not leak pid=0 capitula onto pid=1 and classifies all pid=1 rosette leaves
            as basal regardless of sheet y-position.
        """
        from scripts.vision.organ_crown_metrics import classify_leaf_zone_anatomical

        # Case 1: Inverted mount (roots at top y=197, capitula at bottom y=1572)
        inverted_leaves = [
            {
                "leaf_id": 2,
                "plant_individual_id": 0,
                "p_base_sheet": (723.0, 349.0),
                "p_apex_sheet": (640.0, 480.0),
                "p_petiole_tip": (760.0, 245.0),
                "blade_length_px": 155.0,
                "has_petiole_mask": True,
            },
            {
                "leaf_id": 4,
                "plant_individual_id": 0,
                "p_base_sheet": (1085.0, 1061.0),
                "p_apex_sheet": (1165.0, 1185.0),
                "blade_length_px": 148.0,
                "has_petiole_mask": False,
                "petiole_to_blade_ratio": 0.02,
            },
        ]
        res_case1 = classify_leaf_zone_anatomical(
            inverted_leaves,
            crowns_by_pid={0: (781.0, 206.0)},
            rosette_centroids_by_pid={0: (980.0, 920.0)},
            root_boxes=[((160, 740, 235, 820), 0.88, 7)],
            capitula_detections=[(1520, 1280, 1620, 1360)],
        )
        self.assertEqual(res_case1[2]["leaf_zone"], "basal")
        self.assertEqual(res_case1[4]["leaf_zone"], "cauline")

        # Case 2: Erect long-petiolate basal tuft with Roots only and no Capitula (NCU00438219)
        tuft_leaves = [
            {
                "leaf_id": 1,
                "plant_individual_id": 0,
                "p_base_sheet": (855.0, 1910.0),
                "p_apex_sheet": (845.0, 1730.0),
                "p_petiole_tip": (860.0, 1995.0),
                "ray_dir": (0.04, 0.999),
                "blade_length_px": 180.0,
                "has_petiole_mask": True,
            },
            {
                "leaf_id": 2,
                "plant_individual_id": 0,
                "p_base_sheet": (820.0, 1580.0),
                "p_apex_sheet": (790.0, 1370.0),
                "p_petiole_tip": (832.0, 1665.0),
                "ray_dir": (0.09, 0.996),
                "blade_length_px": 210.0,
                "has_petiole_mask": True,
            },
            {
                "leaf_id": 3,
                "plant_individual_id": 0,
                "p_base_sheet": (910.0, 1595.0),
                "p_apex_sheet": (945.0, 1390.0),
                "p_petiole_tip": (895.0, 1680.0),
                "ray_dir": (-0.09, 0.996),
                "blade_length_px": 205.0,
                "has_petiole_mask": True,
            },
        ]
        res_case2 = classify_leaf_zone_anatomical(
            tuft_leaves,
            crowns_by_pid={0: (863.0, 2027.0)},
            root_boxes=[((1980, 810, 2080, 915), 0.91, 7)],
            capitula_detections=[],
        )
        self.assertEqual(res_case2[1]["leaf_zone"], "basal")
        self.assertEqual(res_case2[2]["leaf_zone"], "basal")
        self.assertEqual(res_case2[3]["leaf_zone"], "basal")

        # Case 3: Capitula only, no Roots, clipped upper flowering stem without basal rosette hub
        clipped_top_leaves = [
            {
                "leaf_id": 10,
                "plant_individual_id": 0,
                "p_base_sheet": (600.0, 1150.0),
                "p_apex_sheet": (510.0, 1060.0),
                "blade_length_px": 125.0,
                "has_petiole_mask": False,
                "petiole_to_blade_ratio": 0.0,
            },
            {
                "leaf_id": 11,
                "plant_individual_id": 0,
                "p_base_sheet": (605.0, 780.0),
                "p_apex_sheet": (700.0, 690.0),
                "blade_length_px": 115.0,
                "has_petiole_mask": False,
                "petiole_to_blade_ratio": 0.0,
            },
        ]
        res_case3 = classify_leaf_zone_anatomical(
            clipped_top_leaves,
            root_boxes=[],
            capitula_detections=[(380, 560, 460, 650)],
        )
        self.assertEqual(res_case3[10]["leaf_zone"], "cauline")
        self.assertEqual(res_case3[11]["leaf_zone"], "cauline")

        # Case 4: Multi-plant sheet where pid=1 is a detached sterile basal rosette (NCSC00103660)
        multi_plant_leaves = [
            # pid=0: flowering plant on right side (x ~ 1350)
            {
                "leaf_id": 1,
                "plant_individual_id": 0,
                "p_base_sheet": (1340.0, 1820.0),
                "p_apex_sheet": (1240.0, 1690.0),
                "p_petiole_tip": (1365.0, 1870.0),
                "blade_length_px": 165.0,
                "has_petiole_mask": True,
            },
            {
                "leaf_id": 2,
                "plant_individual_id": 0,
                "p_base_sheet": (1360.0, 920.0),
                "p_apex_sheet": (1450.0, 810.0),
                "blade_length_px": 140.0,
                "has_petiole_mask": False,
                "petiole_to_blade_ratio": 0.02,
            },
            # pid=1: detached sterile basal rosette on upper-left (x ~ 440, y ~ 725)
            {
                "leaf_id": 3,
                "plant_individual_id": 1,
                "p_base_sheet": (380.0, 740.0),
                "p_apex_sheet": (220.0, 780.0),
                "p_petiole_tip": (435.0, 728.0),
                "blade_length_px": 165.0,
                "has_petiole_mask": True,
            },
            {
                "leaf_id": 4,
                "plant_individual_id": 1,
                "p_base_sheet": (500.0, 755.0),
                "p_apex_sheet": (650.0, 810.0),
                "p_petiole_tip": (448.0, 730.0),
                "blade_length_px": 160.0,
                "has_petiole_mask": True,
            },
            {
                "leaf_id": 5,
                "plant_individual_id": 1,
                "p_base_sheet": (442.0, 650.0),
                "p_apex_sheet": (440.0, 480.0),
                "p_petiole_tip": (442.0, 720.0),
                "blade_length_px": 170.0,
                "has_petiole_mask": True,
            },
        ]
        res_case4 = classify_leaf_zone_anatomical(
            multi_plant_leaves,
            crowns_by_pid={0: (1370.0, 1890.0), 1: (442.0, 725.0)},
            root_boxes=[((1870, 1330, 1960, 1420), 0.90, 7)],
            capitula_detections=[(380, 1310, 470, 1410)],
        )
        self.assertEqual(res_case4[1]["leaf_zone"], "basal")
        self.assertEqual(res_case4[2]["leaf_zone"], "cauline")
        self.assertEqual(res_case4[3]["leaf_zone"], "basal")
        self.assertEqual(res_case4[4]["leaf_zone"], "basal")
        self.assertEqual(res_case4[5]["leaf_zone"], "basal")

    def test_step02_upright_mask_persistence_reniform_lw_crown_tiebreak_and_32x360_ect(self) -> None:
        """
        Verify Step 2 botanical systematics & LM3 downstream handoff invariants:
          1. Oblique / inverted leaves persist upright base-to-tip binary masks in `data/masks/*.png`
             with `mask_is_upright == True` and `32 x 360` ECT `.npy` matrices matching `LM3_settings.yaml`.
          2. Reniform / cordate wider-than-long blades (`W > L`, e.g., Packera aurea) preserve true
             anatomical `aspect_ratio < 1.0` (`blade_length / blade_width`) without `90 deg` sideways rotation.
          3. Downward-radiating rosette / inverted-mount leaves (`NCU00385199` pattern) tie-break `lamina_base`
             toward the root/crown hub rather than flipping `180 deg` via `argmax(y)`.
          4. `RulerCalibrationMap` keeps `predicted_from_megapixels` as `is_calibrated_mm == False`
             (`scale_source == "mp_estimated"`) in both CSV manifests and SQLite `vouchers` / `extracted_leaves`.
          5. Pre-oriented `momocs_outline_xy` passed to `export_standardized_contour` is never double-rotated.
        """
        from scripts.core.db import get_db_connection, record_voucher_segmentation_transaction
        from scripts.vision.leaf_geometry_utils import export_standardized_contour
        from scripts.vision.lm3_data_loader import load_ruler_calibrations

        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_p = Path(tmpdir)
            vouchers_csv = tmp_p / "vouchers.csv"
            img_obl = tmp_p / "NCU_OBLIQUE.jpg"
            img_ren = tmp_p / "NCU_RENIFORM.jpg"
            img_inv = tmp_p / "NCU_INVERTED.jpg"
            for ip in (img_obl, img_ren, img_inv):
                cv2.imwrite(str(ip), np.full((1200, 1000, 3), 235, dtype=np.uint8))

            pd.DataFrame([
                {"catalogNumber": "NCU_OBLIQUE", "scientificName": "Packera anonyma", "image_path": str(img_obl)},
                {"catalogNumber": "NCU_RENIFORM", "scientificName": "Packera aurea", "image_path": str(img_ren)},
                {"catalogNumber": "NCU_INVERTED", "scientificName": "Packera dubia", "image_path": str(img_inv)},
            ]).to_csv(vouchers_csv, index=False)

            # Set up LM3 ruler CSV with one empirical ruler and one megapixel prior
            ruler_dir = tmp_p / "lm3_proj" / "Data" / "Ruler"
            ruler_dir.mkdir(parents=True, exist_ok=True)
            pd.DataFrame([
                {
                    "filename": "NCU_OBLIQUE__R__10-10-100-500",
                    "ruler_success": True,
                    "conversion_mean": 100.0,
                    "ruler_class_confidence": 95.0,
                    "conversion_determination": "measured_from_ruler",
                },
                {
                    "filename": "NCU_RENIFORM__R__10-10-100-500",
                    "ruler_success": True,
                    "conversion_mean": 105.0,
                    "ruler_class_confidence": 88.0,
                    "conversion_determination": "predicted_from_megapixels",
                },
            ]).to_csv(ruler_dir / "ruler_summary.csv", index=False)

            calibs = load_ruler_calibrations(tmp_p / "lm3_proj", min_confidence=0.50)
            self.assertEqual(calibs.sources["NCU_OBLIQUE"], "measured_from_ruler")
            self.assertEqual(calibs.sources["NCU_RENIFORM"], "predicted_from_megapixels")

            pipe = SegmentAndExtractPipeline(
                vouchers_csv=vouchers_csv,
                model_weights=Path("models/dummy.pth"),
                output_dir=tmp_p / "out_data",
                lm3_project_dir=tmp_p / "lm3_proj",
                device="cpu",
                force=True,
            )
            try:
                # 1. Oblique 135-deg leaf (base at top-left, tip at bottom-right) -> saved mask must be upright (H > 1.8 * W)
                obl_crop = np.zeros((420, 420), dtype=np.uint8)
                cv2.ellipse(obl_crop, (210, 210), (110, 42), 45.0, 0, 360, 255, -1)
                obl_pet = np.zeros((420, 420), dtype=np.uint8)
                cv2.line(obl_pet, (132, 132), (92, 92), 255, thickness=6)

                with patch.object(pipe.stage1_detector, "predict", return_value=[((200, 200, 620, 620), 0.96, 0)]):
                    with patch.object(
                        pipe.model_engine,
                        "predict_batch",
                        return_value=[[((10, 10, 410, 410), obl_crop, 0.96, 0), ((85, 85, 140, 140), obl_pet, 0.90, 1)]],
                    ):
                        recs_obl, fail_obl = pipe.process_voucher("NCU_OBLIQUE", img_obl)

                self.assertIsNone(fail_obl)
                self.assertEqual(len(recs_obl), 1)
                r_obl = recs_obl[0]
                self.assertTrue(bool(r_obl["is_calibrated_mm"]))
                self.assertEqual(r_obl["scale_source"], "ruler_calibrated")
                self.assertTrue(bool(r_obl["mask_is_upright"]))
                saved_obl_mask = cv2.imread(str(r_obl["mask_path"]), cv2.IMREAD_GRAYSCALE)
                self.assertIsNotNone(saved_obl_mask)
                ys_o, xs_o = np.where(saved_obl_mask > 0)
                h_obl = float(ys_o.max() - ys_o.min())
                w_obl = float(xs_o.max() - xs_o.min())
                self.assertGreater(h_obl, 1.8 * w_obl, "Saved PNG mask in data/masks/ must be rotated upright base-to-tip.")
                ect_obl = np.load(str(r_obl["ect_path"]))
                self.assertEqual(ect_obl.shape, (32, 360), "Step 02 ECT .npy must honor LM3_settings.yaml 32x360 grid.")

                # 2. Reniform wider-than-long leaf (L ~ 120, W ~ 190 -> L/W < 0.72) + mp_estimated provenance
                ren_crop = np.zeros((380, 380), dtype=np.uint8)
                cv2.ellipse(ren_crop, (190, 160), (95, 62), 0.0, 0, 360, 255, -1)
                cv2.ellipse(ren_crop, (190, 220), (34, 24), 0.0, 0, 360, 0, -1)
                ren_pet = np.zeros((380, 380), dtype=np.uint8)
                cv2.line(ren_pet, (190, 222), (190, 315), 255, thickness=6)

                with patch.object(pipe.stage1_detector, "predict", return_value=[((300, 300, 680, 680), 0.97, 0)]):
                    with patch.object(
                        pipe.model_engine,
                        "predict_batch",
                        return_value=[[((10, 10, 370, 370), ren_crop, 0.96, 0), ((220, 184, 320, 196), ren_pet, 0.91, 1)]],
                    ):
                        recs_ren, fail_ren = pipe.process_voucher("NCU_RENIFORM", img_ren)

                self.assertIsNone(fail_ren)
                r_ren = recs_ren[0]
                self.assertLess(float(r_ren["aspect_ratio"]), 0.80, "Reniform W > L blade must preserve aspect_ratio < 1.0.")
                self.assertFalse(bool(r_ren["is_calibrated_mm"]), "predicted_from_megapixels must not set is_calibrated_mm=True.")
                self.assertEqual(r_ren["scale_source"], "mp_estimated")
                self.assertTrue(pd.isna(r_ren["area_mm2"]))

                # Verify SQLite persistence keeps is_calibrated_mm == 0 for mp_estimated
                db_file = tmp_p / "test_step02_prov.sqlite"
                conn = get_db_connection(db_file, init_schema=True)
                try:
                    record_voucher_segmentation_transaction(
                        conn,
                        catalog_number="NCU_RENIFORM",
                        image_path=str(img_ren),
                        extracted_records=recs_ren,
                        failure_record=None,
                        output_dir=pipe.output_dir,
                    )
                    v_cal = conn.execute("SELECT is_calibrated_mm FROM vouchers WHERE catalogNumber = 'NCU_RENIFORM';").fetchone()[0]
                    el_cal, el_ss = conn.execute(
                        "SELECT is_calibrated_mm, scale_source FROM extracted_leaves WHERE catalogNumber = 'NCU_RENIFORM';"
                    ).fetchone()
                    self.assertEqual(int(v_cal), 0)
                    self.assertEqual(int(el_cal), 0)
                    self.assertEqual(str(el_ss), "mp_estimated")
                finally:
                    conn.close()

                # 3. Inverted herbarium mount (root crown at top y=150, downward-pointing blade lacking petiole_mask)
                #    Draw ovate blade wider near top base (y ~ 310) and tapering to acute apex at bottom (y ~ 480)
                inv_crop = np.zeros((320, 320), dtype=np.uint8)
                pts_inv = np.array([[160, 40], [215, 95], [195, 190], [160, 280], [125, 190], [105, 95]], dtype=np.int32)
                cv2.fillPoly(inv_crop, [pts_inv], 255)
                s1_inv = [
                    ((120, 460, 190, 540), 0.93, 7),   # Root tuft at top of sheet (inverted mount)
                    ((240, 340, 560, 660), 0.95, 0),   # Downward-pointing rosette leaf below top root crown
                ]
                with patch.object(pipe.stage1_detector, "predict", return_value=s1_inv):
                    with patch.object(
                        pipe.model_engine,
                        "predict_batch",
                        return_value=[[((10, 10, 310, 310), inv_crop, 0.95, 0)]],
                    ):
                        recs_inv, fail_inv = pipe.process_voucher("NCU_INVERTED", img_inv)

                self.assertIsNone(fail_inv)
                # Sheet-space midrib points downward (~ +90 deg), so Step 02 rotates ~180 deg to make upright mask tip-up
                self.assertGreater(float(recs_inv[0]["midrib_angle_deg"]), 45.0)
                inv_saved = cv2.imread(str(recs_inv[0]["mask_path"]), cv2.IMREAD_GRAYSCALE)
                # After 180-deg upright rotation, broader lamina base sits in lower half and acute apex sits at top
                h_i = inv_saved.shape[0]
                top_w = int(np.count_nonzero(inv_saved[: h_i // 3, :]))
                bot_w = int(np.count_nonzero(inv_saved[(2 * h_i) // 3 :, :]))
                self.assertGreater(bot_w, top_w, "Crown-hub tie-breaker must orient broader basal half to bottom (+Y) and acute apex to top (-Y).")

                # 4. Double-rotation guard when pre-oriented momocs_outline_xy is supplied alongside sheet keypoints
                th = np.linspace(0.0, 2.0 * np.pi, 120, endpoint=False)
                pre_oriented_xy = np.column_stack([200.0 - 35.0 * np.sin(th), 200.0 + 95.0 * np.cos(th)])
                kp_oblique = {
                    "lamina_base": (120.0, 120.0, 0.95),
                    "lamina_tip": (280.0, 280.0, 0.95),
                }
                c_pre = export_standardized_contour(
                    None,
                    catalog_number="PRE_ORI",
                    leaf_id=1,
                    output_dir=tmp_p / "pre_ori_out",
                    momocs_outline_xy=pre_oriented_xy,
                    keypoints_31=kp_oblique,
                )
                df_pre = pd.read_csv(c_pre)
                span_y_pre = float(df_pre["y"].max() - df_pre["y"].min())
                span_x_pre = float(df_pre["x"].max() - df_pre["x"].min())
                self.assertGreater(span_y_pre, 2.2 * span_x_pre, "Pre-oriented momocs_outline_xy must not be double-rotated by sheet keypoints.")

                # 5. Unbridged focal multi-pinna retention in _reduce_crop_predictions_native_lm3 + 32x360 multi-component ECT
                from scripts.vision import lm3_inference_engines as lm3_eng_chk
                pinna_crop_shape = (320, 240)
                p_inst_u8 = np.zeros(pinna_crop_shape, dtype=np.uint8)
                cv2.ellipse(p_inst_u8, (120, 80), (42, 48), 0, 0, 360, 255, -1)   # terminal lobe
                cv2.ellipse(p_inst_u8, (82, 168), (24, 14), 0, 0, 360, 255, -1)   # detached left pinna
                cv2.ellipse(p_inst_u8, (158, 168), (24, 14), 0, 0, 360, 255, -1)  # detached right pinna
                cv2.rectangle(p_inst_u8, (0, 270), (35, 319), 255, -1)            # border-clipped neighbor bleed -> rejected
                pet_stalk_u8 = np.zeros(pinna_crop_shape, dtype=np.uint8)
                cv2.rectangle(pet_stalk_u8, (110, 122), (130, 260), 255, -1)
                reduced_pinna = lm3_eng_chk._reduce_crop_predictions_native_lm3(
                    [
                        ((0, 0, 320, 240), p_inst_u8, 0.94, 0),
                        ((122, 110, 260, 130), pet_stalk_u8, 0.91, 1),
                    ],
                    pinna_crop_shape,
                    min_area_px=64,
                )
                lam_pinna_u8 = [m for _b, m, _s, cid in reduced_pinna if cid == 0][0]
                n_cc_p, _, stats_p, _ = cv2.connectedComponentsWithStats((lam_pinna_u8 > 0).astype(np.uint8), connectivity=8)
                sig_p = sum(1 for cid in range(1, n_cc_p) if int(stats_p[cid, cv2.CC_STAT_AREA]) >= 64)
                self.assertEqual(sig_p, 3, "Unbridged terminal lobe + 2 petiole-attached lateral pinnae must be retained (rejecting border neighbor bleed).")
                self.assertEqual(int(np.count_nonzero(lam_pinna_u8[130:150, 112:128])), 0, "Zero synthetic petiole bridge pixels may be added to lamina mask.")
                term_only_u8 = np.zeros(pinna_crop_shape, dtype=np.uint8)
                term_only_u8[:135, :] = p_inst_u8[:135, :]
                ect_term_only = lm3_eng_chk.compute_leaf_ect_descriptor(term_only_u8, num_dirs=360, num_thresh=32, already_oriented=True)
                ect_multi_pinna = lm3_eng_chk.compute_leaf_ect_descriptor(lam_pinna_u8, num_dirs=360, num_thresh=32, already_oriented=True)
                self.assertGreater(float(np.max(ect_multi_pinna)), float(np.max(ect_term_only)), "Multi-pinna 32x360 ECT must capture higher multi-component Euler characteristic than terminal-lobe-only.")
                mc_metrics = compute_lm3_shape_metrics(lam_pinna_u8, angle_cw=0.0)
                self.assertLess(mc_metrics["solidity"], 0.82, "Multi-pinna solidity must reflect inter-pinna sinuses across all islands.")
            finally:
                pipe.shutdown()


if __name__ == "__main__":
    unittest.main()
