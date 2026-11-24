#!/usr/bin/env python3
"""Step 02 calibrated blade photometry, bolting scape consolidation, void QC, and damage/fold QC tests."""

from __future__ import annotations

import importlib
import sys
import unittest
from pathlib import Path

import cv2
import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

step02 = importlib.import_module("scripts.pipeline.02_segment_and_extract")
DetectedInstance = step02.DetectedInstance
cluster_plant_individuals = step02.cluster_plant_individuals


class TestStep02PhotometryAndVoidQc(unittest.TestCase):
    """Tests for blade photometry, bolting shoot clustering, enclosed void QC, and fold/damage QC."""

    def test_eroded_blade_cie_l_tomentum_and_paper_strap_masking(self) -> None:
        """Verify compute_calibrated_blade_photometry separates floccose white tomentum from glabrous blades and masks white paper mounting straps, petioles, and holes."""
        from scripts.vision.step02_photometry_corymb_utils import (
            compute_calibrated_blade_photometry,
        )

        rng = np.random.default_rng(42)
        sheet_bgr = np.full((240, 320, 3), (228, 236, 242), dtype=np.uint8)
        mask = np.zeros((120, 120), dtype=np.uint8)
        cv2.circle(mask, (60, 60), 46, 255, -1)

        sheet_glab = sheet_bgr.copy()
        sub_g = sheet_glab[60:180, 40:160].copy()
        glab_tex = np.clip(68.0 + rng.normal(0.0, 2.0, size=(120, 120, 1)), 40, 95).astype(np.uint8)
        leaf_glab_bgr = np.concatenate([glab_tex - 18, glab_tex + 8, glab_tex + 4], axis=2)
        sub_g[mask > 0] = leaf_glab_bgr[mask > 0]
        sheet_glab[60:180, 40:160] = sub_g

        sheet_glab_strap = sheet_glab.copy()
        sheet_glab_strap[114:124, 40:160] = (232, 238, 243)

        pet_m = np.zeros((120, 120), dtype=np.uint8)
        pet_m[90:108, 52:68] = 255
        hole_m = np.zeros((120, 120), dtype=np.uint8)
        cv2.circle(hole_m, (60, 42), 8, 255, -1)
        sheet_glab_excl = sheet_glab.copy()
        sub_ex = sheet_glab_excl[60:180, 40:160].copy()
        sub_ex[(pet_m > 0) | (hole_m > 0)] = (238, 240, 244)
        sheet_glab_excl[60:180, 40:160] = sub_ex

        sheet_tom = sheet_bgr.copy()
        sub_t = sheet_tom[60:180, 40:160].copy()
        wool_grid = rng.normal(0.0, 26.0, size=(30, 30)).astype(np.float32)
        wool_up = cv2.resize(wool_grid, (120, 120), interpolation=cv2.INTER_CUBIC)
        tom_lum = np.clip(168.0 + wool_up, 95.0, 222.0).astype(np.uint8)
        leaf_tom_bgr = np.stack([tom_lum, tom_lum, tom_lum], axis=2)
        sub_t[mask > 0] = leaf_tom_bgr[mask > 0]
        sheet_tom[60:180, 40:160] = sub_t

        paper_l = step02.estimate_sheet_paper_cie_l_ref(sheet_glab)
        self.assertGreater(paper_l, 85.0)

        m_glab = compute_calibrated_blade_photometry(
            sheet_glab, mask, mask_origin=(40, 60), paper_l_ref=paper_l, pixels_per_mm=10.0
        )
        m_strap = compute_calibrated_blade_photometry(
            sheet_glab_strap, mask, mask_origin=(40, 60), paper_l_ref=paper_l, pixels_per_mm=10.0
        )
        m_excl = compute_calibrated_blade_photometry(
            sheet_glab_excl,
            mask,
            mask_origin=(40, 60),
            paper_l_ref=paper_l,
            petiole_mask=pet_m,
            hole_mask=hole_m,
            pixels_per_mm=10.0,
        )
        m_tom = compute_calibrated_blade_photometry(
            sheet_tom, mask, mask_origin=(40, 60), paper_l_ref=paper_l, pixels_per_mm=10.0
        )

        self.assertGreater(m_tom["blade_cie_l_norm"], m_glab["blade_cie_l_norm"] + 0.25)
        self.assertGreater(m_tom["blade_cie_l_sd"], m_glab["blade_cie_l_sd"] * 2.5)
        self.assertAlmostEqual(m_glab["blade_cie_l_norm"], m_strap["blade_cie_l_norm"], delta=0.02)
        self.assertAlmostEqual(m_glab["blade_cie_l_norm"], m_excl["blade_cie_l_norm"], delta=0.02)

    def test_bolting_scape_specimen_mask_and_splayed_leaf_consolidation(self) -> None:
        """Verify cluster_plant_individuals reunites bolting stems, mask-bridged shoots, and 1-leaf strays while separating side-by-side rooted plants."""
        dummy_mask = np.ones((12, 12), dtype=np.uint8)
        # Single bolting plant on a 2000x3000 sheet: basal rosette at y~2550, upper cauline leaves at y~1250 (dy=1300 > 0.15*3000=450),
        # plus 1 splayed outer basal leaf 270 px laterally from the rosette (<= 0.12*3000=360 px).
        b1 = DetectedInstance(catalog_number="BOLT1", leaf_id=1, bbox=(2500, 860, 2620, 960), mask=dummy_mask, score=0.95, class_id=0)
        b2 = DetectedInstance(catalog_number="BOLT1", leaf_id=2, bbox=(2520, 930, 2640, 1030), mask=dummy_mask, score=0.94, class_id=0)
        c1 = DetectedInstance(catalog_number="BOLT1", leaf_id=3, bbox=(1200, 890, 1300, 970), mask=dummy_mask, score=0.92, class_id=0)
        c2 = DetectedInstance(catalog_number="BOLT1", leaf_id=4, bbox=(1120, 910, 1210, 980), mask=dummy_mask, score=0.91, class_id=0)
        splay = DetectedInstance(catalog_number="BOLT1", leaf_id=5, bbox=(2530, 1200, 2630, 1300), mask=dummy_mask, score=0.89, class_id=0)

        out_single = cluster_plant_individuals([b1, b2, c1, c2, splay], sheet_width=2000, sheet_height=3000)
        self.assertEqual([inst.plant_individual_id for inst in out_single], [0, 0, 0, 0, 0])

        # Diagonal V-bent plant bridged by Specimen_Segmenter mask across > 0.22 * sheet_width
        v1 = DetectedInstance(catalog_number="VBENT", leaf_id=1, bbox=(2400, 300, 2520, 420), mask=dummy_mask, score=0.95, class_id=0)
        v2 = DetectedInstance(catalog_number="VBENT", leaf_id=2, bbox=(2420, 360, 2540, 480), mask=dummy_mask, score=0.93, class_id=0)
        v3 = DetectedInstance(catalog_number="VBENT", leaf_id=3, bbox=(900, 1400, 1020, 1520), mask=dummy_mask, score=0.92, class_id=0)
        v4 = DetectedInstance(catalog_number="VBENT", leaf_id=4, bbox=(940, 1450, 1060, 1560), mask=dummy_mask, score=0.90, class_id=0)
        spec_mask = np.zeros((3000, 2000), dtype=np.uint8)
        cv2.line(spec_mask, (390, 2460), (1460, 960), 255, thickness=18)
        out_bridged = cluster_plant_individuals([v1, v2, v3, v4], sheet_width=2000, sheet_height=3000, specimen_mask=spec_mask)
        self.assertEqual([inst.plant_individual_id for inst in out_bridged], [0, 0, 0, 0])

        # Two genuinely separate side-by-side rooted rosettes remain distinct [0, 0, 1, 1]
        p0_a = DetectedInstance(catalog_number="TWO_ROOT", leaf_id=1, bbox=(2300, 250, 2420, 360), mask=dummy_mask, score=0.95, class_id=0)
        p0_b = DetectedInstance(catalog_number="TWO_ROOT", leaf_id=2, bbox=(2320, 310, 2440, 420), mask=dummy_mask, score=0.94, class_id=0)
        p1_a = DetectedInstance(catalog_number="TWO_ROOT", leaf_id=3, bbox=(2300, 1480, 2420, 1590), mask=dummy_mask, score=0.93, class_id=0)
        p1_b = DetectedInstance(catalog_number="TWO_ROOT", leaf_id=4, bbox=(2320, 1540, 2440, 1650), mask=dummy_mask, score=0.92, class_id=0)
        roots = [(2460, 280, 2580, 380), (2460, 1510, 2580, 1610)]
        out_two = cluster_plant_individuals([p0_a, p0_b, p1_a, p1_b], sheet_width=2000, sheet_height=3000, root_boxes=roots)
        self.assertEqual([inst.plant_individual_id for inst in out_two], [0, 0, 1, 1])

    def test_unfilled_dissected_masks_paper_pockets_and_void_qc(self) -> None:
        """
        Verify:
          1. `_reduce_crop_predictions_native_lm3` and `_unpack_reduced_lm3_crop_preds` preserve unfilled
             `cls == 0` (`Leaf`) interior `0`-pixel sinus voids and intersect interior mounting-paper
             pockets with `specimen_paper_cleaned_mask`.
          2. Read-only `measure_enclosed_interlobe_void_frac` leaves input masks 100% pixel-identical,
             quantifies enclosed interlobe void fraction (`> 0.020`), excludes known `hole_mask` overlaps,
             and drives void-aware `compute_lm3_shape_metrics` solidity plus `ENCLOSED_INTERLOBE_VOID` QC.
        """
        from scripts.vision.leaf_geometry_utils import compute_lm3_shape_metrics
        from scripts.vision.lm3_inference_engines import (
            _reduce_crop_predictions_native_lm3,
            evaluate_leaf_damage_and_fold_qc,
        )
        from scripts.vision.step02_leaf_traits_utils import (
            measure_enclosed_interlobe_void_frac,
        )

        # 1A. Dissected leaf mask with enclosed sinus pocket (no cls==2 Hole prediction)
        raw_ring = np.zeros((120, 120), dtype=np.uint8)
        cv2.circle(raw_ring, (60, 60), 42, 1, -1)
        cv2.circle(raw_ring, (60, 60), 16, 0, -1)
        crop_det_tuple = [((18, 18, 102, 102), raw_ring, 0.93, 0)]

        red_dets = _reduce_crop_predictions_native_lm3(crop_det_tuple, (120, 120))
        self.assertEqual(len(red_dets), 1)
        _red_bbox, red_lam, red_sc, red_cls = red_dets[0]
        self.assertEqual(red_cls, 0)
        self.assertAlmostEqual(red_sc, 0.93, places=2)
        self.assertEqual(int(red_lam[60, 60]), 0, "Unfilled cls==0 Leaf mask must preserve interior sinus zero-pocket.")

        built_lam, built_sc, built_pet, built_hole = step02._unpack_reduced_lm3_crop_preds(red_dets, (120, 120))
        self.assertIsNotNone(built_lam)
        self.assertAlmostEqual(built_sc, 0.93, places=2)
        self.assertIsNone(built_pet)
        self.assertIsNone(built_hole)
        self.assertEqual(int(built_lam[60, 60]), 0, "_unpack_reduced_lm3_crop_preds must not flood-fill outer Leaf polygons.")

        # 1B. Interior mounting-paper pocket intersected via specimen_paper_cleaned_mask
        yy, xx = np.ogrid[:140, :140]
        solid_bool = (((yy - 70) ** 2) / (52.0 ** 2) + ((xx - 70) ** 2) / (38.0 ** 2)) <= 1.0
        solid_u8 = solid_bool.astype(np.uint8) * 255
        pocket_reg = ((yy - 52) ** 2 + (xx - 70) ** 2) <= (9.0 ** 2)
        spec_clean = solid_u8.copy()
        spec_clean[pocket_reg] = 0

        red_solid = _reduce_crop_predictions_native_lm3([((18, 32, 122, 108), solid_u8, 0.92, 0)], (140, 140))
        unpacked_lam, _, _, _ = step02._unpack_reduced_lm3_crop_preds(
            red_solid,
            (140, 140),
            specimen_paper_cleaned_mask=spec_clean,
        )
        self.assertIsNotNone(unpacked_lam)
        self.assertEqual(int(unpacked_lam[52, 70]), 0, "Interior pocket absent from specimen_paper_cleaned_mask must be carved out.")
        self.assertGreater(int(unpacked_lam[90, 70]), 0, "Intact lamina interior must be preserved.")

        # 2. Read-only enclosed interlobe void measurement, hole_mask exclusion, and ENCLOSED_INTERLOBE_VOID QC
        broad_mask = np.zeros((140, 140), dtype=np.uint8)
        cv2.ellipse(broad_mask, (70, 70), (44, 54), 0, 0, 360, 255, -1)
        broad_mask[52:88, 54:86] = 0
        orig_copy = broad_mask.copy()

        evf_broad = measure_enclosed_interlobe_void_frac(broad_mask)
        self.assertTrue(np.array_equal(broad_mask, orig_copy), "measure_enclosed_interlobe_void_frac must be strictly read-only.")
        self.assertGreater(evf_broad, 0.05)

        hole_m = np.zeros_like(broad_mask)
        hole_m[52:88, 54:86] = 255
        evf_with_hole = measure_enclosed_interlobe_void_frac(broad_mask, hole_mask=hole_m)
        self.assertTrue(np.array_equal(broad_mask, orig_copy))
        self.assertEqual(evf_with_hole, 0.0, "Known hole_mask overlaps (>30%) must be excluded from enclosed void fraction.")

        filled_ref = broad_mask.copy()
        cv2.ellipse(filled_ref, (70, 70), (44, 54), 0, 0, 360, 255, -1)
        m_void = compute_lm3_shape_metrics(broad_mask)
        m_fill = compute_lm3_shape_metrics(filled_ref)
        self.assertLess(m_void["solidity"], m_fill["solidity"] - 0.08, "Void-aware solidity must subtract interior interlobe voids.")

        qc_broad = evaluate_leaf_damage_and_fold_qc(
            oriented_mask_u8=broad_mask,
            solidity=m_void["solidity"],
            bilateral_dice=0.92,
            si_a=0.02,
            lm3_gates_pass=True,
            enclosed_void_frac=evf_broad,
        )
        self.assertFalse(qc_broad["efa_qc_pass"])
        self.assertIn("ENCLOSED_INTERLOBE_VOID", qc_broad["qc_flag_reason"])

    def test_folded_damaged_leaf_qc_gates_and_bilateral_dissection_symmetry(self) -> None:
        """
        Verify `evaluate_leaf_damage_and_fold_qc` accurately separates:
          - `LEAF_PARTIAL_FRAGMENT`, `CROP_BORDER_CLIPPED`, `SEVERE_HERBIVORY_HOLES`
          - `FOLDED_MIDVEIN_TRACE` and `APICAL_TIP_FOLD` on 31-keypoint pose predictions
          - `LONGITUDINAL_MARGIN_CREASE` on half-creased ("taco") leaves
          - `UNILATERAL_ASYMMETRY_OR_TEAR` on one-sided torn blades with `solidity < 0.72`
          - `PASS` (`bilateral_lobe_symmetry == True`) on genuine bilaterally pinnatifid Packera leaves.
        """
        from scripts.vision.lm3_inference_engines import evaluate_leaf_damage_and_fold_qc

        ovate = np.zeros((160, 120), dtype=np.uint8)
        cv2.ellipse(ovate, (60, 80), (32, 60), 0, 0, 360, 255, -1)

        r_part = evaluate_leaf_damage_and_fold_qc(ovate, detector_class_id=1, segment_partial_leaves=False)
        self.assertFalse(r_part["efa_qc_pass"])
        self.assertIn("LEAF_PARTIAL_FRAGMENT", r_part["qc_flag_reason"])

        r_clip = evaluate_leaf_damage_and_fold_qc(ovate, is_crop_clipped=True)
        self.assertFalse(r_clip["efa_qc_pass"])
        self.assertIn("CROP_BORDER_CLIPPED", r_clip["qc_flag_reason"])

        r_hole = evaluate_leaf_damage_and_fold_qc(ovate, hole_area_frac=0.055)
        self.assertFalse(r_hole["efa_qc_pass"])
        self.assertIn("SEVERE_HERBIVORY_HOLES", r_hole["qc_flag_reason"])

        creased = ovate.copy()
        creased[:, :58] = 0
        r_crease = evaluate_leaf_damage_and_fold_qc(
            creased, solidity=0.93, bilateral_dice=0.82, si_a=0.12, lm3_gates_pass=True
        )
        self.assertFalse(r_crease["efa_qc_pass"])
        self.assertIn("LONGITUDINAL_MARGIN_CREASE", r_crease["qc_flag_reason"])

        bent_kps = {
            "lamina_tip": (60.0, 22.0, 0.92),
            "midvein_3": (66.0, 48.0, 0.90),
            "midvein_7": (96.0, 80.0, 0.91),
            "midvein_11": (68.0, 112.0, 0.89),
            "lamina_base": (60.0, 138.0, 0.93),
        }
        r_bent = evaluate_leaf_damage_and_fold_qc(
            ovate, keypoints_31=bent_kps, is_geometric_fallback_pose=False
        )
        self.assertFalse(r_bent["efa_qc_pass"])
        self.assertIn("FOLDED_MIDVEIN_TRACE", r_bent["qc_flag_reason"])

        tucked_kps = {
            "lamina_tip": (60.0, 52.0, 0.88),
            "midvein_7": (60.0, 95.0, 0.90),
            "lamina_base": (60.0, 138.0, 0.92),
        }
        r_tip = evaluate_leaf_damage_and_fold_qc(
            ovate,
            bilateral_dice=0.76,
            si_a=0.18,
            oriented_points={"lamina_tip": (60.0, 52.0), "lamina_base": (60.0, 138.0)},
            keypoints_31=tucked_kps,
            is_geometric_fallback_pose=False,
        )
        self.assertFalse(r_tip["efa_qc_pass"])
        self.assertIn("APICAL_TIP_FOLD", r_tip["qc_flag_reason"])

        torn_one_side = np.zeros((180, 140), dtype=np.uint8)
        for y_r in range(25, 155):
            frac = (y_r - 25.0) / 130.0
            hw_r = int(round(8.0 + 34.0 * np.sin(np.pi * frac) * (0.45 + 0.55 * (np.sin(3.0 * np.pi * frac) ** 2))))
            hw_l = 6
            torn_one_side[y_r, 70 - hw_l : 70 + hw_r] = 255
        r_torn = evaluate_leaf_damage_and_fold_qc(
            torn_one_side,
            solidity=0.61,
            bilateral_dice=0.60,
            si_a=0.24,
            lm3_gates_pass=False,
            oriented_points={"lamina_tip": (70.0, 25.0), "lamina_base": (70.0, 154.0)},
        )
        self.assertFalse(r_torn["bilateral_lobe_symmetry"])
        self.assertFalse(r_torn["efa_qc_pass"])
        self.assertIn("UNILATERAL_ASYMMETRY_OR_TEAR", r_torn["qc_flag_reason"])

        pinnatifid = np.zeros((180, 140), dtype=np.uint8)
        for y_r in range(25, 155):
            frac = (y_r - 25.0) / 130.0
            hw = int(round(7.0 + 36.0 * (0.35 + 0.65 * (np.sin(4.0 * np.pi * frac) ** 2))))
            pinnatifid[y_r, 70 - hw : 70 + hw] = 255
        r_pinn = evaluate_leaf_damage_and_fold_qc(
            pinnatifid,
            solidity=0.56,
            bilateral_dice=0.88,
            si_a=0.04,
            lm3_gates_pass=False,
            oriented_points={"lamina_tip": (70.0, 25.0), "lamina_base": (70.0, 154.0)},
        )
        self.assertTrue(r_pinn["bilateral_lobe_symmetry"])
        self.assertTrue(r_pinn["efa_qc_pass"], f"Expected bilateral pinnatifid leaf to PASS, got {r_pinn}")
        self.assertEqual(r_pinn["qc_flag_reason"], "PASS")


if __name__ == "__main__":
    unittest.main()
