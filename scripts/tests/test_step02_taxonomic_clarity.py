#!/usr/bin/env python3
"""
Step 02 read-only LeafMachine3 taxonomic clarity & photometry test suite:
  1. Wool-preserving CIELAB photometry, q85/basal-zone deltas, midrib safety, and abaxial flags.
  2. Native LM3 landmark angle/curvature/trace forwarding and read-only Packera margin/base traits.
  3. 180-deg obovate polarity swap, 2D stacked multi-root clustering, zero petiole inflation, and cordate blade_length_mm.
  4. Imperial ruler rejection, ochre/top-down capitulum gating, radiate/fruiting phenology, bbox rescore, and multi-pid DB round-trip.
"""

from __future__ import annotations

import importlib, math, sqlite3, sys, tempfile, unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

import cv2, numpy as np, pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from scripts.core.db import get_db_connection, load_step02_records_from_db, record_voucher_segmentation_transaction
from scripts.vision.capitulum_phenology_classifier import HEAD_ANTHESIS, HEAD_FRUIT, VOUCHER_ANTHESIS, VOUCHER_FRUIT, classify_single_capitulum, classify_voucher_phenology
from scripts.vision.lm3_inference_engines import DetectedInstance, LM3PoseResult, LM3ScaleCalibrator, compute_capitulum_metrics, orient_and_score_leaf_lm3
from scripts.vision.organ_crown_metrics import approximate_petiole_and_scape_metrics, cluster_plant_individuals
from scripts.vision.rescore_step02_qc_and_zones import rescore_tier_a_leaves
from scripts.vision.step02_leaf_traits_utils import evaluate_rosette_abaxial_leaf_flags, extract_packera_leaf_margin_and_base_traits
from scripts.vision.step02_photometry_corymb_utils import compute_calibrated_blade_photometry

SegmentAndExtractPipeline = importlib.import_module("scripts.pipeline.02_segment_and_extract").SegmentAndExtractPipeline


def _make_cordate_mask(shape: tuple[int, int] = (220, 180)) -> np.ndarray:
    """Build a proportional upright bilateral cordate mask (for (220, 180): apex y=22, sinus y=152, auricles y=194)."""
    h, w = shape
    sx, sy, cx = w / 180.0, h / 220.0, w // 2
    mask = np.zeros(shape, dtype=np.uint8)
    cv2.ellipse(mask, (cx, int(round(92 * sy))), (max(2, int(round(52 * sx))), max(2, int(round(70 * sy)))), 0, 0, 360, 255, -1)
    r_aur = max(2, int(round(28 * min(sx, sy))))
    for dx in (-26, 26):
        cv2.circle(mask, (cx + int(round(dx * sx)), int(round(166 * sy))), r_aur, 255, -1)
    notch = np.array([[cx, int(round(152 * sy))], [cx - int(round(18 * sx)), h - 1], [cx + int(round(18 * sx)), h - 1]], dtype=np.int32)
    cv2.fillPoly(mask, [notch], 0)
    return mask


def _make_margin_mask(mode: str) -> np.ndarray:
    """Build upright 280x180 masks with entire, rounded crenate, or sharp serrate margins."""
    mask = np.zeros((280, 180), dtype=np.uint8)
    for y in range(20, 261):
        t = (y - 20.0) / 240.0
        hw = 46.0 * (math.sin(math.pi * t) ** 0.55)
        if 0.14 <= t <= 0.86:
            if mode == "serrate":
                hw += 12.0 * ((1.0 - abs(2.0 * (((y - 20) % 16) / 16.0) - 1.0)) ** 1.2)
            elif mode == "crenate":
                hw += 8.0 * (math.sin(math.pi * (((y - 20) % 24) / 24.0)) ** 0.45)
        iw = max(1, int(round(hw)))
        mask[y, 90 - iw : 90 + iw + 1] = 255
    return mask


def _make_full_kps(cx: float = 90.0, y_tip: float = 22.0, y_base: float = 152.0, conf: float = 0.95, bow_px: float = 5.0) -> dict[str, tuple[float, float, float]]:
    """Build a complete 31-keypoint LM3 dict including apex/base angle arms and midvein_1..15."""
    ym = 0.5 * (y_tip + y_base)
    return {
        "lamina_tip": (cx, y_tip, conf), "apex_1": (cx - 18.0, y_tip + 22.0, conf), "apex_2": (cx, y_tip, conf), "apex_3": (cx + 18.0, y_tip + 22.0, conf),
        "lamina_base": (cx, y_base, conf), "base_1": (cx - 24.0, y_base - 18.0, conf), "base_2": (cx, y_base, conf), "base_3": (cx + 24.0, y_base - 18.0, conf),
        "width_1": (cx - 48.0, ym, conf), "width_2": (cx + 48.0, ym, conf),
        **{f"midvein_{i}": (cx + bow_px * math.sin(math.pi * i / 16.0), y_tip + (i / 16.0) * (y_base - y_tip), conf) for i in range(1, 16)},
    }


def _fibrous_gray_bgr(xx: np.ndarray, mean_val: float, amp: float, lo: int, hi: int) -> np.ndarray:
    """Synthesize high-frequency fibrous wool luminance (DoG MAD >> 1.8) in BGR."""
    lum = np.clip(np.round(mean_val + amp * np.sin(2.0 * np.pi * xx / 8.0)), lo, hi).astype(np.uint8)
    return np.stack([lum, lum, lum], axis=2)


def _extract_pid_maps(res_pid: object, cat: str) -> tuple[dict, dict]:
    """Resolve `{pid: rec}` maps for `cat` from `load_step02_records_from_db(..., include_by_pid=True)`."""
    r_raw, p_raw = res_pid[4], res_pid[5]  # type: ignore[index]
    r_map = {int(k[1]): v for k, v in r_raw.items() if isinstance(k, tuple) and k[0] == cat} if any(isinstance(k, tuple) for k in r_raw) else r_raw.get(cat, r_raw)
    p_map = {int(k[1]): v for k, v in p_raw.items() if isinstance(k, tuple) and k[0] == cat} if any(isinstance(k, tuple) for k in p_raw) else p_raw.get(cat, p_raw)
    return r_map, p_map


class TestStep02TaxonomicClarity(unittest.TestCase):
    """Focused verification of Step 02 read-only LM3 taxonomic clarity upgrades."""

    def test_wool_preserving_photometry_midrib_and_abaxial_flags(self) -> None:
        """Verify wool-safe strap exclusion, q85/basal deltas, narrow/midrib safety, and abaxial flags."""
        paper_bgr = (226, 228, 230)
        sheet = np.full((220, 220, 3), paper_bgr, dtype=np.uint8)
        mask = np.zeros((220, 220), dtype=np.uint8)
        cv2.circle(mask, (110, 110), 60, 255, -1)
        sheet_glab = sheet.copy()
        sheet_glab[mask > 0] = (46, 90, 56)

        # Smooth rim-to-rim cream linen strap (DoG MAD <= 1.8, ar_rect >= 2.8) vs fibrous snow-white wool (DoG MAD > 1.8)
        sheet_strap, sheet_wool, sheet_patch = sheet_glab.copy(), sheet_glab.copy(), sheet_glab.copy()
        sheet_strap[104:116, 45:175] = (224, 234, 236)
        yy, xx = np.mgrid[:220, :220]
        w_reg = (mask > 0) & (yy < 134)
        sheet_wool[w_reg] = _fibrous_gray_bgr(xx, 238.5, 16.5, 222, 255)[w_reg]
        p_reg = (mask > 0) & (((xx - 90) ** 2 + (yy - 90) ** 2) <= (30 ** 2))
        sheet_patch[p_reg] = _fibrous_gray_bgr(xx, 212.0, 16.0, 190, 236)[p_reg]
        m_glab = compute_calibrated_blade_photometry(sheet_glab, mask, pixels_per_mm=10.0)
        m_strap = compute_calibrated_blade_photometry(sheet_strap, mask, pixels_per_mm=10.0)
        m_wool = compute_calibrated_blade_photometry(sheet_wool, mask, pixels_per_mm=10.0)
        m_patch = compute_calibrated_blade_photometry(sheet_patch, mask, pixels_per_mm=10.0)
        self.assertAlmostEqual(m_glab["blade_cie_l_norm"], m_strap["blade_cie_l_norm"], delta=0.03)
        self.assertLess(m_strap["blade_cie_l_q85_sheet_delta"], -25.0)
        self.assertTrue(m_wool["blade_cie_l_norm"] > 0.85 and m_wool["blade_cie_l_sheet_delta"] > -8.0)
        self.assertTrue(m_patch["blade_cie_l_sheet_delta"] < -35.0 and m_patch["blade_cie_l_q85_sheet_delta"] > -15.0)

        # Persistent basal wool tuft near lamina_base captured by basal_zone_cie_l_delta
        tall_m = np.zeros((220, 140), dtype=np.uint8)
        cv2.ellipse(tall_m, (70, 110), (32, 74), 0, 0, 360, 255, -1)
        tall_glab = np.full((220, 140, 3), paper_bgr, dtype=np.uint8)
        tall_glab[tall_m > 0] = (46, 90, 56)
        tall_tuft, (yt, xt) = tall_glab.copy(), np.mgrid[:220, :140]
        t_reg = (tall_m > 0) & (yt >= 142)
        tall_tuft[t_reg] = _fibrous_gray_bgr(xt, 210.0, 16.0, 188, 234)[t_reg]
        kps_tall = _make_full_kps(cx=70.0, y_tip=36.0, y_base=184.0, bow_px=0.0)
        m_no_tuft = compute_calibrated_blade_photometry(tall_glab, tall_m, pixels_per_mm=10.0, keypoints_31=kps_tall)
        m_tuft = compute_calibrated_blade_photometry(tall_tuft, tall_m, pixels_per_mm=10.0, keypoints_31=kps_tall)
        self.assertGreater(m_tuft["basal_zone_cie_l_delta"], m_tuft["blade_cie_l_sheet_delta"] + 20.0)
        self.assertGreater(m_tuft["basal_zone_cie_l_delta"], m_no_tuft["basal_zone_cie_l_delta"] + 25.0)

        # Narrow 6-px blade stays finite; wider blade excludes pale midrib corridor via keypoints_31
        nar_m, nar_s = np.zeros((120, 80), dtype=np.uint8), np.full((120, 80, 3), paper_bgr, dtype=np.uint8)
        nar_m[15:105, 37:43], nar_s[15:105, 37:43] = 255, (46, 90, 56)
        m_nar = compute_calibrated_blade_photometry(nar_s, nar_m, pixels_per_mm=10.0, keypoints_31=_make_full_kps(40.0, 18.0, 102.0, bow_px=0.0))
        self.assertTrue(len(m_nar) == 7 and all(np.isfinite(m_nar[k]) for k in m_nar), f"Narrow blade failed: {m_nar}")
        rib_m, rib_s = np.zeros((200, 160), dtype=np.uint8), np.full((200, 160, 3), paper_bgr, dtype=np.uint8)
        rib_m[20:180, 72:88], rib_s[20:180, 72:88], rib_s[20:180, 79:82] = 255, (46, 90, 56), (180, 218, 224)
        m_rib_no = compute_calibrated_blade_photometry(rib_s, rib_m, pixels_per_mm=15.0)
        m_rib_kp = compute_calibrated_blade_photometry(rib_s, rib_m, pixels_per_mm=15.0, keypoints_31=_make_full_kps(80.0, 22.0, 178.0, bow_px=0.0))
        self.assertLess(m_rib_kp["blade_cie_l_q85_sheet_delta"], m_rib_no["blade_cie_l_q85_sheet_delta"] - 15.0)

        # Within-rosette abaxial pressed-leaf flagging across 4 basal leaves + 1 cauline control
        recs = [
            {"catalogNumber": "ABAX_1", "plant_individual_id": 0, "leaf_id": i, "leaf_zone": "basal", "blade_cie_l_norm": 0.43, "blade_cie_l_sheet_delta": -45.0 + i, "blade_cie_a_delta": -3.0, "blade_cie_b_delta": 8.0}
            for i in (1, 2, 3)
        ] + [
            {"catalogNumber": "ABAX_1", "plant_individual_id": 0, "leaf_id": 4, "leaf_zone": "basal", "blade_cie_l_norm": 0.72, "blade_cie_l_sheet_delta": -22.0, "blade_cie_a_delta": -1.0, "blade_cie_b_delta": 4.0},
            {"catalogNumber": "ABAX_1", "plant_individual_id": 0, "leaf_id": 5, "leaf_zone": "cauline", "blade_cie_l_norm": 0.74, "blade_cie_l_sheet_delta": -20.0, "blade_cie_a_delta": -1.0, "blade_cie_b_delta": 4.0},
        ]
        evaluate_rosette_abaxial_leaf_flags(recs)
        self.assertTrue(recs[3]["is_abaxial_pressed_leaf"] and math.isclose(recs[3]["abaxial_cie_l_sheet_delta"], -22.0, abs_tol=1e-3))
        self.assertTrue(np.isnan(recs[3]["blade_cie_l_norm"]) and np.isnan(recs[3]["blade_cie_l_sheet_delta"]))
        for idx in (0, 1, 2, 4):
            self.assertTrue(not recs[idx]["is_abaxial_pressed_leaf"] and np.isfinite(recs[idx]["blade_cie_l_sheet_delta"]))

    def test_lm3_landmark_angles_and_packera_margin_base_traits(self) -> None:
        """Verify LM3 angle/curvature/trace forwarding and read-only Packera margin, base, and lobe traits."""
        ovate = np.zeros((220, 180), dtype=np.uint8)
        cv2.ellipse(ovate, (90, 110), (44, 78), 0, 0, 360, 255, -1)
        ori_res = orient_and_score_leaf_lm3(ovate, keypoints_31=_make_full_kps(90.0, 32.0, 188.0, bow_px=6.0))
        for key in ("base_angle_deg", "apex_angle_deg", "lamina_curvature", "lamina_trace_length_px"):
            self.assertTrue(key in ori_res and np.isfinite(float(ori_res[key])), f"Invalid {key}: {ori_res.get(key)}")

        cordate = _make_cordate_mask()
        cord_orig = cordate.copy()
        t_cord = extract_packera_leaf_margin_and_base_traits(cordate, {"lamina_tip": (90.0, 22.0), "lamina_base": (90.0, 152.0)})
        self.assertTrue(np.array_equal(cordate, cord_orig) and t_cord["basal_cordate_index"] > 0.05)

        cuneate = np.zeros((220, 180), dtype=np.uint8)
        cv2.fillPoly(cuneate, [np.array([[90, 25], [35, 72], [90, 195], [145, 72]], dtype=np.int32)], 255)
        self.assertLess(extract_packera_leaf_margin_and_base_traits(cuneate, {"lamina_tip": (90.0, 25.0), "lamina_base": (90.0, 195.0)})["basal_cordate_index"], -0.05)

        t_ent = extract_packera_leaf_margin_and_base_traits(_make_margin_mask("entire"))
        t_cre = extract_packera_leaf_margin_and_base_traits(_make_margin_mask("crenate"))
        t_ser = extract_packera_leaf_margin_and_base_traits(_make_margin_mask("serrate"))
        self.assertTrue(t_ent["margin_tooth_acuteness"] == 0.0 and t_ent["margin_serration_index"] < 0.02)
        self.assertTrue(t_ser["margin_serration_index"] > 0.0 and t_ser["margin_tooth_acuteness"] > 1.0)
        self.assertTrue(0.0 < t_cre["margin_tooth_acuteness"] < 0.6 and t_cre["margin_tooth_acuteness"] < t_ser["margin_tooth_acuteness"])

        lyrate = np.zeros((260, 200), dtype=np.uint8)
        cv2.ellipse(lyrate, (100, 68), (54, 46), 0, 0, 360, 255, -1)
        cv2.rectangle(lyrate, (94, 108), (106, 232), 255, -1)
        cv2.ellipse(lyrate, (100, 156), (42, 18), 0, 0, 360, 255, -1)
        cv2.ellipse(lyrate, (100, 210), (36, 16), 0, 0, 360, 255, -1)
        self.assertTrue(0.2 < extract_packera_leaf_margin_and_base_traits(lyrate)["terminal_lobe_area_frac"] < 0.85)
        self.assertEqual(extract_packera_leaf_margin_and_base_traits(ovate)["terminal_lobe_area_frac"], 1.0)

        tiny = np.zeros((20, 20), dtype=np.uint8)
        tiny[5:10, 5:10] = 255
        self.assertTrue(all(np.isnan(v) for v in extract_packera_leaf_margin_and_base_traits(tiny).values()))

    def test_obovate_polarity_stacked_roots_petiole_and_cordate_blade_mm(self) -> None:
        """Verify 180-deg obovate swap, 2D stacked root split, zero petiole inflation, and cordate upright blade_length_mm."""
        obovate = np.zeros((260, 200), dtype=np.uint8)
        cv2.fillPoly(obovate, [np.array([[100, 34], [52, 68], [86, 186], [100, 192], [114, 186], [148, 68]], dtype=np.int32)], 255)
        pet_m = np.zeros((260, 200), dtype=np.uint8)
        cv2.rectangle(pet_m, (96, 188), (104, 244), 255, -1)
        inv_kps = {
            "lamina_tip": (100.0, 190.0, 0.78), "lamina_base": (100.0, 40.0, 0.82),
            **{f"midvein_{i}": (100.0, 190.0 - i * (150.0 / 16.0), 0.80) for i in range(1, 16)},
        }
        ori_ob = orient_and_score_leaf_lm3(obovate, keypoints_31=inv_kps, petiole_mask=pet_m, crown_hint_xy=(100.0, 255.0))
        self.assertLess(abs(((float(ori_ob["angle_cw"]) + 180.0) % 360.0) - 180.0), 20.0)
        om = ori_ob["oriented_mask"]
        self.assertGreater(int(np.count_nonzero(om[: om.shape[0] // 2])), int(np.count_nonzero(om[om.shape[0] // 2 :])) * 1.2)

        # Vertically stacked multi-root sheet connected by a touching specimen_mask bridge
        d_m = np.ones((12, 12), dtype=np.uint8)
        stacked = [DetectedInstance("STK", i + 1, b, d_m, 0.95 - 0.01 * i, 0) for i, b in enumerate([(1020, 420, 1120, 500), (1030, 500, 1130, 580), (2220, 420, 2320, 500), (2230, 500, 2330, 580)])]
        spec_m = np.zeros((3000, 2000), dtype=np.uint8)
        cv2.line(spec_m, (500, 1050), (500, 2400), 255, thickness=24)
        cl_out = cluster_plant_individuals(stacked, 2000, 3000, specimen_mask=spec_m, root_boxes=[(1160, 450, 1240, 550), (2360, 450, 2440, 550)])
        self.assertEqual([inst.plant_individual_id for inst in cl_out], [0, 0, 1, 1])

        # Zero tip_to_crown_px inflation when petiole_1..6 pose exists
        pet_kps = {"lamina_tip": (500.0, 1480.0, 0.95), "lamina_base": (500.0, 1600.0, 0.95), **{f"petiole_{j}": (500.0, 1600.0 + float(j * 10), 0.92) for j in range(1, 7)}}
        pet_res = approximate_petiole_and_scape_metrics(
            basal_leaves=[{"leaf_id": 1, "plant_individual_id": 0, "blade_length_px": 120.0, "keypoints_31": pet_kps}],
            p_crown=(500.0, 1780.0), pixels_per_mm=10.0,
        )
        self.assertAlmostEqual(float(pet_res["petiole_length_px"]), 60.0, delta=1.5)
        self.assertAlmostEqual(float(pet_res["petiole_length_mm"]), 6.0, delta=0.15)

        # Cordate upright blade_length_mm preserved from inst.blade_length_px / pixels_per_mm
        cord_m = _make_cordate_mask((220, 180))
        cord_kps = LM3PoseResult(_make_full_kps(cx=90.0, y_tip=22.0, y_base=152.0, bow_px=0.0))
        with tempfile.TemporaryDirectory() as tmp_d:
            tmp_p = Path(tmp_d)
            img_p, v_csv = tmp_p / "CORD_1.jpg", tmp_p / "v.csv"
            img_p.write_bytes(b"\xff\xd8\xff" + b"0" * 128)
            pd.DataFrame([{"catalogNumber": "CORD_1", "image_path": str(img_p)}]).to_csv(v_csv, index=False)
            pipe = SegmentAndExtractPipeline(vouchers_csv=v_csv, model_weights=Path("models/dummy.pth"), output_dir=tmp_p / "data", device="cpu", force=True, save_debug_crops=False)
            try:
                pipe._is_dummy_test_weights = False
                cal_ret = {"pixels_per_mm": 10.0, "is_calibrated_mm": True, "cf_source": "measured_from_ruler", "scale_source": "ruler_calibrated"}
                with patch("cv2.imread", return_value=np.full((800, 600, 3), 235, dtype=np.uint8)), \
                     patch.object(pipe.scale_calibrator, "calibrate_sheet", return_value=cal_ret), \
                     patch.object(pipe.stage1_detector, "predict", return_value=[((200, 200, 420, 380), 0.96, 0)]), \
                     patch.object(pipe.model_engine, "predict", return_value=[((10, 10, 210, 170), cord_m, 0.95, 0)]), \
                     patch.object(pipe.pose_engine, "predict", return_value=cord_kps), \
                     patch.object(pipe.pose_engine, "predict_landmarks", return_value=cord_kps), \
                     patch.object(pipe.pose_engine, "predict_batch", return_value=[cord_kps]):
                    bundle = pipe.process_voucher("CORD_1", img_p, _return_bundle=True)
            finally:
                pipe.shutdown()
            self.assertAlmostEqual(float(bundle.extracted_records[0]["blade_length_mm"]), 17.2, delta=0.8)

    def test_imperial_ruler_capitulum_and_radiate_fruiting_phenology(self) -> None:
        """Verify imperial ruler rejection, ochre/top-down capitulum gating, and radiate/fruiting phenology."""
        cal = LM3ScaleCalibrator(device="cpu", providers=[], use_cf_predicted_by_mp=False)
        cal.archival_session = object()
        cal.ruler_lattice = SimpleNamespace(predict=lambda _rgb: SimpleNamespace(px_per_cm=120.0, confidence=0.95))
        blank = np.full((1200, 1600, 3), 240, dtype=np.uint8)
        for imp_lbl in ("INCH", "STD_IN", "IMPERIAL"):
            cal.ruler_classifier = SimpleNamespace(predict=lambda _c, lbl=imp_lbl: {"label": lbl, "conf": 0.90})
            with patch.object(cal, "_detect_ruler_boxes", return_value=[((40, 100, 180, 900), 0.95)]):
                self.assertFalse(cal.calibrate_sheet(blank, catalog_number="NCU_IMP", use_cf_predicted_by_mp=False)["is_calibrated_mm"])

        # Lateral class_id=4 ochre exclusion vs class_id=6 retention, plus top-down radial disc NaN gating
        ochre_bgr = tuple(int(c) for c in cv2.cvtColor(np.uint8([[[14, 64, 180]]]), cv2.COLOR_HSV2BGR)[0, 0])
        lat_crop = np.full((120, 100, 3), 238, dtype=np.uint8)
        cv2.rectangle(lat_crop, (35, 65), (65, 105), (35, 95, 30), -1)
        cv2.rectangle(lat_crop, (25, 22), (75, 65), ochre_bgr, -1)
        m_lat = compute_capitulum_metrics(None, (0, 0, 120, 100), pixels_per_mm=10.0, crop_bgr=lat_crop, class_id=4)
        m_bud = compute_capitulum_metrics(None, (0, 0, 120, 100), pixels_per_mm=10.0, crop_bgr=lat_crop, class_id=6)
        self.assertTrue(m_lat is not None and bool(m_lat.get("is_lateral_profile", True)) and float(m_lat["involucre_height_px"]) < 55.0)
        self.assertGreater(float(m_bud["involucre_height_px"]), float(m_lat["involucre_height_px"]) + 15.0)

        top_crop = np.full((120, 120, 3), 238, dtype=np.uint8)
        yel_bgr = tuple(int(c) for c in cv2.cvtColor(np.uint8([[[25, 220, 235]]]), cv2.COLOR_HSV2BGR)[0, 0])
        cv2.circle(top_crop, (60, 60), 44, yel_bgr, -1)
        cv2.circle(top_crop, (60, 60), 24, (35, 95, 30), -1)
        m_top = compute_capitulum_metrics(None, (0, 0, 120, 120), pixels_per_mm=10.0, crop_bgr=top_crop, class_id=4)
        self.assertTrue(m_top is not None and not bool(m_top.get("is_lateral_profile", True)))
        self.assertTrue(all(np.isnan(m_top[k]) for k in ("involucre_height_px", "involucre_width_px", "AR_inv")))

        # Radiate vs discoid ligules and fruiting corymb ochre guard
        rad_rgb = np.full((120, 120, 3), 240, dtype=np.uint8)
        cv2.rectangle(rad_rgb, (44, 58), (76, 114), (35, 95, 30), -1)
        cv2.circle(rad_rgb, (60, 44), 14, (235, 175, 20), -1)
        for ang in (0, 30, 60, 120, 150, 180):
            cv2.ellipse(rad_rgb, (60, 44), (50, 7), ang, 0, 360, (235, 175, 20), -1)
        disc_rgb = np.full((100, 100, 3), 240, dtype=np.uint8)
        cv2.rectangle(disc_rgb, (24, 46), (76, 98), (35, 95, 30), -1)
        cv2.ellipse(disc_rgb, (50, 34), (23, 19), 0, 0, 360, (235, 175, 20), -1)
        c_rad, c_disc = classify_single_capitulum(rad_rgb, class_id=4), classify_single_capitulum(disc_rgb, class_id=4)
        self.assertTrue(c_rad["has_radiate_ligules"] and not c_disc["has_radiate_ligules"])
        self.assertAlmostEqual(float(classify_voucher_phenology("RAD_1", [c_rad, c_rad, c_disc])["radiate_head_fraction"]), 2.0 / 3.0, delta=0.02)

        f_heads = [{"head_state": HEAD_FRUIT, "has_pappus_plume": True, "has_fresh_yellow_corollas": False, "has_aged_ochre_corollas": False, "has_radiate_ligules": False, "white_fibrous_ratio": 0.42} for _ in range(4)]
        ochre_fp = {"head_state": HEAD_ANTHESIS, "has_pappus_plume": False, "has_fresh_yellow_corollas": False, "has_aged_ochre_corollas": True, "has_radiate_ligules": False, "white_fibrous_ratio": 0.02}
        self.assertEqual(classify_voucher_phenology("FR_GUARD", f_heads + [ochre_fp])["voucher_phenological_state"], VOUCHER_FRUIT)
        self.assertEqual(classify_voucher_phenology("FR_FRESH", f_heads + [{**ochre_fp, "has_fresh_yellow_corollas": True, "has_aged_ochre_corollas": False, "yellow_ratio": 0.25}])["voucher_phenological_state"], VOUCHER_ANTHESIS)

    def test_bbox_photometry_rescore_and_multipid_sqlite_roundtrip(self) -> None:
        """Verify bbox-backed CPU rescore (photometry + abaxial + margin traits) and multi-pid SQLite round-trip."""
        with tempfile.TemporaryDirectory() as tmp_d:
            data_dir = Path(tmp_d) / "data"
            for sub in ("raw_vouchers", "masks", "contours", "tables"):
                (data_dir / sub).mkdir(parents=True, exist_ok=True)
            sheet = np.full((300, 400, 3), (234, 238, 240), dtype=np.uint8)
            img_p, db_file = data_dir / "raw_vouchers" / "BBOX_1.png", data_dir / "plantmatics.sqlite"
            cord_sm, oval_sm = _make_cordate_mask((100, 80)), np.zeros((100, 80), dtype=np.uint8)
            cv2.ellipse(oval_sm, (40, 50), (24, 38), 0, 0, 360, 255, -1)
            boxes = [(20, 40, 100, 140), (110, 40, 190, 140), (200, 40, 280, 140), (290, 40, 370, 140)]
            m_paths = []
            for idx, ((x1, y1, x2, y2), m_i, bgr_i) in enumerate(zip(boxes, [cord_sm, oval_sm, oval_sm, oval_sm], [(46, 90, 56)] * 3 + [(200, 202, 200)]), start=1):
                sheet[y1:y2, x1:x2][m_i > 0] = bgr_i
                mp = data_dir / "masks" / f"BBOX_1_leaf{idx}.png"
                cv2.imwrite(str(mp), m_i)
                m_paths.append(str(mp))
            cv2.imwrite(str(img_p), sheet)

            conn = get_db_connection(db_path=db_file, init_schema=True)
            try:
                with conn:
                    conn.execute("INSERT INTO vouchers (catalogNumber, scientificName, species_standardized, determiner_tier, qc_status, image_path, pixels_per_mm, is_calibrated_mm) VALUES ('BBOX_1', 'Packera aurea', 'Packera aurea', 'Tier_1_Gold', 'PASS', ?, 10.0, 1);", (str(img_p),))
                    conn.execute("INSERT INTO plant_individuals (catalogNumber, plant_individual_id, crown_x, crown_y, plant_height_px, is_bent_stem) VALUES ('BBOX_1', 0, 200.0, 260.0, 500.0, 0);")
                    conn.execute("INSERT INTO voucher_reproductive_phenology (catalogNumber, plant_individual_id, capitula_count) VALUES ('BBOX_1', 0, 4);")
                    conn.executemany(
                        "INSERT INTO extracted_leaves (catalogNumber, plant_individual_id, leaf_id, leaf_zone, leaf_zone_confidence, dist_from_crown_norm, solidity, circularity, aspect_ratio, blade_length_px, petiole_length_px, petiole_to_blade_ratio, area_px, midrib_angle_deg, pixels_per_mm, is_calibrated_mm, scale_source, bilateral_dice, si_a, lm3_gates_pass, n_lamina_components, cy_rel_plant_bbox, mask_path, mask_is_upright, crop_path, contour_xy_json, bbox_x1, bbox_y1, bbox_x2, bbox_y2) VALUES ('BBOX_1', 0, ?, 'basal', 0.92, 0.05, 0.90, 0.80, 1.5, 80.0, 50.0, 0.62, 2200.0, 0.0, 10.0, 1, 'ruler_calibrated', 0.94, 0.02, 1, 1, 0.95, ?, 1, NULL, '[[40,88],[16,50],[40,12],[64,50]]', ?, ?, ?, ?);",
                        [(idx, mp, x1, y1, x2, y2) for idx, ((x1, y1, x2, y2), mp) in enumerate(zip(boxes, m_paths), start=1)],
                    )
            finally:
                conn.close()

            man_csv = data_dir / "tables" / "extracted_leaf_manifest.csv"
            pd.DataFrame([{"catalogNumber": "BBOX_1", "plant_individual_id": 0, "leaf_id": i, "leaf_zone": "basal"} for i in range(1, 5)]).to_csv(man_csv, index=False)
            stats = rescore_tier_a_leaves(db_path=db_file, data_dir=data_dir, manifest_csv=man_csv, recompute_photometry=True, workers=1)
            self.assertGreaterEqual(stats["photometry_recomputed"], 4)
            with sqlite3.connect(str(db_file)) as c_r:
                c_r.row_factory = sqlite3.Row
                rows = {int(r["leaf_id"]): dict(r) for r in c_r.execute("SELECT * FROM extracted_leaves WHERE catalogNumber='BBOX_1';").fetchall()}
            self.assertTrue(float(rows[1]["basal_cordate_index"]) > 0.05 and float(rows[1]["blade_cie_l_sheet_delta"]) < -30.0)
            self.assertTrue(int(rows[4]["is_abaxial_pressed_leaf"]) == 1 and rows[4]["blade_cie_l_norm"] is None and float(rows[4]["abaxial_cie_l_sheet_delta"]) > -25.0)
            man_out = pd.read_csv(man_csv).set_index("leaf_id")
            self.assertTrue(bool(man_out.loc[4, "is_abaxial_pressed_leaf"]) and float(man_out.loc[1, "basal_cordate_index"]) > 0.05)

            # Multi-pid SQLite round-trip covering corymb_head_count_est, radiate_head_fraction, and all 16 new columns
            new16 = {
                "blade_cie_l_q85_sheet_delta": -18.5, "basal_zone_cie_l_delta": -12.4, "is_abaxial_pressed_leaf": False, "abaxial_cie_l_sheet_delta": None,
                "base_angle_deg": 108.5, "apex_angle_deg": 74.2, "lamina_curvature": 9.4, "lamina_trace_length_px": 142.6,
                "basal_cordate_index": 0.19, "margin_serration_index": 0.042, "margin_tooth_acuteness": 1.35, "terminal_lobe_area_frac": 0.62,
                "bbox_x1": 25, "bbox_y1": 45, "bbox_x2": 185, "bbox_y2": 245,
            }
            conn2 = get_db_connection(db_path=db_file, init_schema=True)
            try:
                with conn2:
                    conn2.execute("INSERT INTO vouchers (catalogNumber, scientificName, species_standardized, determiner_tier) VALUES ('MPID_1', 'Packera aurea', 'Packera aurea', 'Tier_1_Gold');")
                ext_recs = [
                    {"catalogNumber": "MPID_1", "plant_individual_id": 0, "leaf_id": 1, "solidity": 0.91, "pixels_per_mm": 10.0, "is_calibrated_mm": True, "contour_xy_json": "[[10,20],[20,10],[30,20]]", **new16},
                    {"catalogNumber": "MPID_1", "plant_individual_id": 1, "leaf_id": 2, "solidity": 0.88, "pixels_per_mm": 10.0, "is_calibrated_mm": True, "contour_xy_json": "[[12,22],[22,12],[32,22]]",
                     **{**new16, "is_abaxial_pressed_leaf": True, "abaxial_cie_l_sheet_delta": -19.2, "bbox_x1": 220, "bbox_x2": 360}},
                ]
                repro_pid = {0: {"capitula_count": 7, "corymb_head_count_est": 11, "AR_inv": 1.45}, 1: {"capitula_count": 4, "corymb_head_count_est": 6, "AR_inv": 1.18}}
                pheno_pid = {
                    0: {"total_capitula": 7, "n_anthesis": 5, "n_fruiting": 2, "voucher_phenological_state": VOUCHER_ANTHESIS, "radiate_head_fraction": 0.85},
                    1: {"total_capitula": 4, "n_anthesis": 0, "n_fruiting": 4, "voucher_phenological_state": VOUCHER_FRUIT, "radiate_head_fraction": 0.0},
                }
                record_voucher_segmentation_transaction(conn2, "MPID_1", ext_recs, repro_record=repro_pid[0], pheno_record=pheno_pid[0], repro_by_pid=repro_pid, pheno_by_pid=pheno_pid, output_dir=data_dir)
                res4 = load_step02_records_from_db(conn2, ["MPID_1"])
                self.assertTrue(len(res4) == 4 and getattr(res4[2], "by_pid", None) is not None and getattr(res4[3], "by_pid", None) is not None)
                res_pid = load_step02_records_from_db(conn2, ["MPID_1"], include_by_pid=True)
                ext_ld = {int(r["plant_individual_id"]): r for r in res_pid[0]["MPID_1"]}
                self.assertEqual(set(ext_ld.keys()), {0, 1})
                for col, exp_v in new16.items():
                    got_v = ext_ld[0].get(col)
                    if exp_v is None:
                        self.assertIsNone(got_v)
                    else:
                        self.assertEqual(bool(got_v), exp_v) if isinstance(exp_v, bool) else self.assertAlmostEqual(float(got_v), float(exp_v), places=3)
                self.assertTrue(bool(ext_ld[1]["is_abaxial_pressed_leaf"]) and math.isclose(float(ext_ld[1]["abaxial_cie_l_sheet_delta"]), -19.2, abs_tol=1e-3))
                r_map, p_map = _extract_pid_maps(res_pid, "MPID_1")
                self.assertEqual((int(r_map[0]["corymb_head_count_est"]), int(r_map[1]["corymb_head_count_est"])), (11, 6))
                self.assertTrue(math.isclose(float(p_map[0]["radiate_head_fraction"]), 0.85, abs_tol=1e-3) and math.isclose(float(p_map[1]["radiate_head_fraction"]), 0.0, abs_tol=1e-3))
            finally:
                conn2.close()


if __name__ == "__main__":
    unittest.main()
