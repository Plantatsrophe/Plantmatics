#!/usr/bin/env python3
"""Step 02 flowering-shoot basal-vs-cauline zonation, CPU rescorer, and missing-crop photometry preservation tests."""

from __future__ import annotations

import importlib
import sys
import tempfile
import unittest
from pathlib import Path

import cv2
import numpy as np
import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

step02 = importlib.import_module("scripts.pipeline.02_segment_and_extract")


class TestStep02ShootZonationAndRescorer(unittest.TestCase):
    """Tests for flowering-shoot zonation tuning and read-only Step 02 CPU rescorer persistence."""

    def test_flowering_shoot_basal_vs_cauline_tuning_and_cpu_rescorer(self) -> None:
        """
        Verify:
          1. Subpetiolate lower-cauline stem leaves (`s_stem ~ 0.30`) on flowering shoots are classified
             `cauline` rather than `basal` even when angled toward the crown.
          2. Unrooted upper flowering stem fragments with sessile alternate leaves stay `cauline` (Case 4B).
          3. `rescore_step02_qc_and_zones.py` migrates SQLite, refreshes QC gates,
             measures read-only enclosed interlobe voids, and updates `leaf_zone` across SQLite and `extracted_leaf_manifest.csv`.
        """
        from scripts.core.db import get_db_connection
        from scripts.core.db_schema import PACKERA_SCHEMA_USER_VERSION
        from scripts.vision.organ_crown_metrics import classify_leaf_zone_anatomical
        from scripts.vision.rescore_step02_qc_and_zones import rescore_tier_a_leaves

        leaves_flowering = [
            {
                "leaf_id": 1,
                "plant_individual_id": 0,
                "p_petiole_tip": (495.0, 2190.0),
                "p_base_sheet": (460.0, 2120.0),
                "p_apex_sheet": (420.0, 2030.0),
                "blade_length_px": 98.0,
                "petiole_length_px": 80.0,
                "petiole_to_blade_ratio": 0.82,
                "has_petiole_mask": True,
            },
            {
                "leaf_id": 2,
                "plant_individual_id": 0,
                "p_petiole_tip": (506.0, 2185.0),
                "p_base_sheet": (548.0, 2115.0),
                "p_apex_sheet": (590.0, 2025.0),
                "blade_length_px": 99.0,
                "petiole_length_px": 81.0,
                "petiole_to_blade_ratio": 0.82,
                "has_petiole_mask": True,
            },
            {
                "leaf_id": 3,
                "plant_individual_id": 0,
                "p_petiole_tip": (502.0, 1660.0),
                "p_base_sheet": (530.0, 1635.0),
                "p_apex_sheet": (585.0, 1585.0),
                "blade_length_px": 74.0,
                "petiole_length_px": 22.0,
                "petiole_to_blade_ratio": 0.30,
                "has_petiole_mask": True,
            },
        ]
        zones_fl = classify_leaf_zone_anatomical(
            leaves_flowering,
            crowns_by_pid={0: (500.0, 2200.0)},
            capitula_detections=[(380.0, 480.0, 420.0, 520.0)],
        )
        self.assertEqual(zones_fl[1]["leaf_zone"], "basal")
        self.assertEqual(zones_fl[2]["leaf_zone"], "basal")
        self.assertEqual(zones_fl[3]["leaf_zone"], "cauline")

        leaves_cut_stem = [
            {
                "leaf_id": 10,
                "plant_individual_id": 0,
                "p_petiole_tip": (400.0, 1280.0),
                "p_base_sheet": (400.0, 1280.0),
                "p_apex_sheet": (455.0, 1210.0),
                "blade_length_px": 89.0,
                "petiole_length_px": 0.0,
                "petiole_to_blade_ratio": 0.0,
                "has_petiole_mask": False,
            },
            {
                "leaf_id": 11,
                "plant_individual_id": 0,
                "p_petiole_tip": (398.0, 920.0),
                "p_base_sheet": (398.0, 920.0),
                "p_apex_sheet": (345.0, 855.0),
                "blade_length_px": 84.0,
                "petiole_length_px": 0.0,
                "petiole_to_blade_ratio": 0.0,
                "has_petiole_mask": False,
            },
            {
                "leaf_id": 12,
                "plant_individual_id": 0,
                "p_petiole_tip": (402.0, 610.0),
                "p_base_sheet": (402.0, 610.0),
                "p_apex_sheet": (450.0, 550.0),
                "blade_length_px": 77.0,
                "petiole_length_px": 0.0,
                "petiole_to_blade_ratio": 0.0,
                "has_petiole_mask": False,
            },
        ]
        zones_cut = classify_leaf_zone_anatomical(
            leaves_cut_stem,
            crowns_by_pid={},
            capitula_detections=[(240.0, 380.0, 280.0, 420.0)],
        )
        self.assertEqual([zones_cut[k]["leaf_zone"] for k in (10, 11, 12)], ["cauline", "cauline", "cauline"])

        with tempfile.TemporaryDirectory() as tmp_d:
            tmp_p = Path(tmp_d)
            data_dir = tmp_p / "data"
            masks_dir = data_dir / "masks"
            contours_dir = data_dir / "contours"
            tables_dir = data_dir / "tables"
            for d in (masks_dir, contours_dir, tables_dir):
                d.mkdir(parents=True, exist_ok=True)

            db_file = data_dir / "plantmatics.sqlite"
            conn = get_db_connection(db_path=db_file, init_schema=True)
            try:
                uv = int(conn.execute("PRAGMA user_version;").fetchone()[0])
                self.assertEqual(uv, PACKERA_SCHEMA_USER_VERSION)
                self.assertEqual(uv, 20261127)

                m1 = np.zeros((150, 110), dtype=np.uint8)
                cv2.ellipse(m1, (55, 75), (30, 54), 0, 0, 360, 255, -1)
                m2 = m1.copy()
                m2[:, :53] = 0
                m3 = np.zeros((130, 96), dtype=np.uint8)
                cv2.ellipse(m3, (48, 65), (22, 46), 0, 0, 360, 255, -1)

                p_m1 = masks_dir / "RESC_1_leaf1.png"
                p_m2 = masks_dir / "RESC_1_leaf2.png"
                p_m3 = masks_dir / "RESC_1_leaf3.png"
                cv2.imwrite(str(p_m1), m1)
                cv2.imwrite(str(p_m2), m2)
                cv2.imwrite(str(p_m3), m3)

                dummy_cjson = "[[55.0,129.0],[25.0,75.0],[55.0,21.0],[85.0,75.0]]"
                with conn:
                    conn.execute(
                        "INSERT INTO vouchers (catalogNumber, scientificName, species_standardized, determiner_tier, qc_status) VALUES ('RESC_1', 'Packera dubia', 'Packera dubia', 'Tier_1_Gold', 'PASS');"
                    )
                    conn.execute(
                        "INSERT INTO plant_individuals (catalogNumber, plant_individual_id, crown_x, crown_y, plant_height_px, is_bent_stem) VALUES ('RESC_1', 0, 500.0, 2200.0, 1600.0, 0);"
                    )
                    conn.execute(
                        "INSERT INTO voucher_reproductive_phenology (catalogNumber, plant_individual_id, capitula_count) VALUES ('RESC_1', 0, 6);"
                    )
                    rows_ins = [
                        ("RESC_1", 0, 1, "basal", 0.92, 0.05, 0.94, 0.85, 1.80, 95.0, 75.0, 0.79, 4500.0, 0.94, 0.03, 1, 1, 0.95, str(p_m1), dummy_cjson),
                        ("RESC_1", 0, 2, "basal", 0.90, 0.06, 0.93, 0.80, 3.60, 95.0, 70.0, 0.74, 2250.0, 0.82, 0.12, 1, 1, 0.93, str(p_m2), dummy_cjson),
                        ("RESC_1", 0, 3, "basal", 0.70, 0.20, 0.92, 0.82, 2.09, 80.0, 16.0, 0.20, 2800.0, 0.93, 0.04, 1, 1, 0.62, str(p_m3), dummy_cjson),
                    ]
                    conn.executemany(
                        """
                        INSERT INTO extracted_leaves (
                            catalogNumber, plant_individual_id, leaf_id, leaf_zone, leaf_zone_confidence,
                            dist_from_crown_norm, solidity, circularity, aspect_ratio, blade_length_px,
                            petiole_length_px, petiole_to_blade_ratio, area_px, bilateral_dice, si_a,
                            lm3_gates_pass, n_lamina_components, cy_rel_plant_bbox, mask_path, contour_xy_json
                        ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?);
                        """,
                        rows_ins,
                    )
                pd.DataFrame(
                    [
                        {"catalogNumber": "RESC_1", "plant_individual_id": 0, "leaf_id": 1, "leaf_zone": "basal", "solidity": 0.94},
                        {"catalogNumber": "RESC_1", "plant_individual_id": 0, "leaf_id": 2, "leaf_zone": "basal", "solidity": 0.93},
                        {"catalogNumber": "RESC_1", "plant_individual_id": 0, "leaf_id": 3, "leaf_zone": "basal", "solidity": 0.92},
                    ]
                ).to_csv(tables_dir / "extracted_leaf_manifest.csv", index=False)
            finally:
                conn.close()

            stats = rescore_tier_a_leaves(
                db_path=db_file,
                data_dir=data_dir,
                manifest_csv=tables_dir / "extracted_leaf_manifest.csv",
            )
            self.assertEqual(stats["total"], 3)
            self.assertEqual(stats["efa_qc_pass"], 2)
            self.assertEqual(stats["flagged"], 1)
            self.assertEqual(stats["basal"], 2)
            self.assertEqual(stats["cauline"], 1)

            man_after = pd.read_csv(tables_dir / "extracted_leaf_manifest.csv").set_index("leaf_id")
            self.assertTrue(bool(man_after.loc[1, "efa_qc_pass"]))
            self.assertFalse(bool(man_after.loc[2, "efa_qc_pass"]))
            self.assertIn("LONGITUDINAL_MARGIN_CREASE", str(man_after.loc[2, "qc_flag_reason"]))
            self.assertEqual(str(man_after.loc[1, "leaf_zone"]), "basal")
            self.assertEqual(str(man_after.loc[3, "leaf_zone"]), "cauline")

    def test_rescore_recompute_photometry_preserves_existing_values_when_crop_missing(self) -> None:
        import sqlite3
        from scripts.core.db import get_db_connection
        from scripts.vision.rescore_step02_qc_and_zones import rescore_tier_a_leaves

        with tempfile.TemporaryDirectory() as tmp_d:
            data_dir = Path(tmp_d) / "data"
            for sub in ("raw_vouchers", "masks", "contours", "tables"):
                (data_dir / sub).mkdir(parents=True, exist_ok=True)
            img_path = data_dir / "raw_vouchers" / "PHOT_1.png"
            mask_path = data_dir / "masks" / "PHOT_1_leaf1.png"
            cv2.imwrite(str(img_path), np.full((200, 300, 3), 240, dtype=np.uint8))
            mask = np.zeros((80, 60), dtype=np.uint8)
            cv2.ellipse(mask, (30, 40), (16, 26), 0, 0, 360, 255, -1)
            cv2.imwrite(str(mask_path), mask)

            db_file = data_dir / "plantmatics.sqlite"
            conn = get_db_connection(db_path=db_file, init_schema=True)
            try:
                with conn:
                    conn.execute(
                        "INSERT INTO vouchers (catalogNumber, scientificName, species_standardized, determiner_tier, qc_status, image_path, is_calibrated_mm) VALUES ('PHOT_1', 'Packera dubia', 'Packera dubia', 'Tier_1_Gold', 'PASS', ?, 1);",
                        (str(img_path),),
                    )
                    conn.execute("INSERT INTO plant_individuals (catalogNumber, plant_individual_id, crown_x, crown_y, plant_height_px, is_bent_stem) VALUES ('PHOT_1', 0, 150.0, 170.0, 140.0, 0);")
                    conn.execute("INSERT INTO voucher_reproductive_phenology (catalogNumber, plant_individual_id, capitula_count) VALUES ('PHOT_1', 0, 2);")
                    conn.execute(
                        """
                        INSERT INTO extracted_leaves (
                            catalogNumber, plant_individual_id, leaf_id, leaf_zone, leaf_zone_confidence,
                            dist_from_crown_norm, solidity, circularity, aspect_ratio, blade_length_px,
                            petiole_length_px, petiole_to_blade_ratio, area_px, bilateral_dice, si_a,
                            lm3_gates_pass, n_lamina_components, cy_rel_plant_bbox, mask_path, crop_path,
                            contour_xy_json, blade_cie_l_norm, blade_cie_l_sd
                        ) VALUES ('PHOT_1', 0, 1, 'basal', 0.90, 0.05, 0.94, 0.84, 1.62, 52.0, 30.0, 0.58, 1300.0, 0.95, 0.02, 1, 1, 0.85, ?, NULL, '[[30,66],[14,40],[30,14],[46,40]]', 0.712, 0.045);
                        """,
                        (str(mask_path),),
                    )
            finally:
                conn.close()

            man_csv = data_dir / "tables" / "extracted_leaf_manifest.csv"
            pd.DataFrame([{"catalogNumber": "PHOT_1", "plant_individual_id": 0, "leaf_id": 1, "leaf_zone": "basal", "solidity": 0.94, "blade_cie_l_norm": 0.712, "blade_cie_l_sd": 0.045}]).to_csv(man_csv, index=False)
            stats = rescore_tier_a_leaves(db_path=db_file, data_dir=data_dir, manifest_csv=man_csv, recompute_photometry=True, workers=2)
            self.assertEqual(stats["photometry_recomputed"], 0)
            with sqlite3.connect(str(db_file)) as c_chk:
                l_norm, l_sd = c_chk.execute("SELECT blade_cie_l_norm, blade_cie_l_sd FROM extracted_leaves WHERE catalogNumber='PHOT_1';").fetchone()
            self.assertAlmostEqual(float(l_norm), 0.712, places=5)
            self.assertAlmostEqual(float(l_sd), 0.045, places=5)
            man_out = pd.read_csv(man_csv)
            self.assertAlmostEqual(float(man_out.loc[0, "blade_cie_l_norm"]), 0.712, places=5)
            self.assertAlmostEqual(float(man_out.loc[0, "blade_cie_l_sd"]), 0.045, places=5)


if __name__ == "__main__":
    unittest.main()
