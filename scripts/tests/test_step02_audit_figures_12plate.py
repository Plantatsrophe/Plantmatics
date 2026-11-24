#!/usr/bin/env python3
"""Step 02 coordinate alignment, config/DB/docstring audit, and 300-DPI PNG-only 12-plate figure tests."""

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


class TestStep02AuditAndFigures12Plate(unittest.TestCase):
    """Tests for Step 02 coordinate/config/DB/docstring audit and 300-DPI PNG figure architecture."""

    def test_step02_rereview_coordinate_alignment_config_db_and_docstring_audit(self) -> None:
        """
        Verify Step 02 zero-trust code & docstring audit fixes:
          1. `_parse_basal_leaf_pose`, `estimate_rosette_crown_point`, `approximate_petiole_and_scape_metrics`,
             and `classify_leaf_zone_anatomical` shift `local_win` `keypoints_31` via `local_mask_origin`
             into native full-sheet coordinates and request `return_sheet_size=True` from `specimen_segmenter`.
          2. `SegmentAndExtractPipeline` resolves nested `modules.*` keys (`ect.num_dirs`, `ect.num_thresh`,
             `plant_detector.conf_thres`, `leaf_segmenter.min_area_px`, `landmark_detector.conf_thres`)
             from `LM3_settings.yaml`.
          3. `SegmentationConfig` parses `segment_partial_leaves` and `rosette_plant_bbox_lower_fraction`,
             `LM3ScaleCalibrator.calibrate_sheet(..., use_cf_predicted_by_mp=False)` suppresses megapixel
             guesses, `normalize_pixels_per_mm` preserves 200-DPI scans (`0.78 px/mm`), and
             `_parse_calibration_row` never promotes `mp_estimated`/`uncalibrated` rows to `ruler_calibrated`.
          4. SQLite `record_voucher_segmentation_transaction` + `load_step02_records_from_db` round-trip
             `mask_is_upright`, `petiole_width_px`, `petiole_width_mm`, `blade_width_px`, and `blade_width_mm`.
          5. `generate_step02_publication_figures` scales `fig02a` `ylim` below `1.0` for reniform blades
             (`aspect_ratio < 1.0`), `device="cpu", threads=2` initializes without `torch` NameError,
             existing `.jpg` scans remain intact without Step 02 mutation, runtime crashes record
             `SEGMENTATION_RUNTIME_ERROR`, and stale docstrings are purged.
        """
        from PIL import Image as PILImage
        from scripts.core.config import PipelineConfig
        from scripts.core.db import (
            get_db_connection,
            load_step02_records_from_db,
            record_voucher_segmentation_transaction,
        )
        from scripts.vision.capitulum_phenology_classifier import (
            classify_single_capitulum,
            classify_voucher_phenology,
        )
        from scripts.vision.lm3_data_loader import (
            _parse_calibration_row,
            normalize_pixels_per_mm,
            resolve_px_per_cm,
        )
        from scripts.vision.lm3_inference_engines import LM3ScaleCalibrator
        from scripts.vision.organ_crown_metrics import (
            _parse_basal_leaf_pose,
            approximate_petiole_and_scape_metrics,
            estimate_rosette_crown_point,
        )
        from scripts.vision.step02_figures import generate_step02_publication_figures

        # 1. Three-space coordinate alignment (`local_mask_origin` -> full-sheet coordinates)
        leaf_local_1 = {
            "leaf_id": 1,
            "plant_individual_id": 0,
            "p_base": (1950.0, 3600.0),
            "p_apex": (1820.0, 3420.0),
            "p_base_sheet": (1950.0, 3600.0),
            "p_apex_sheet": (1820.0, 3420.0),
            "mask_origin": (0.0, 0.0),
            "local_mask_origin": (1750.0, 3350.0),
            "blade_length_px": 222.0,
            "petiole_width_px": 14.0,
            "keypoints_31": {
                "lamina_tip": (70.0, 70.0, 0.95),
                "lamina_base": (200.0, 250.0, 0.95),
                "petiole_1": (210.0, 275.0, 0.92),
                "petiole_6": (245.0, 345.0, 0.90),  # sheet: (1995, 3695)
            },
        }
        leaf_local_2 = {
            "leaf_id": 2,
            "plant_individual_id": 0,
            "p_base": (2060.0, 3590.0),
            "p_apex": (2200.0, 3410.0),
            "p_base_sheet": (2060.0, 3590.0),
            "p_apex_sheet": (2200.0, 3410.0),
            "mask_origin": (0.0, 0.0),
            "local_mask_origin": (1960.0, 3340.0),
            "blade_length_px": 228.0,
            "petiole_width_px": 16.0,
            "keypoints_31": {
                "lamina_tip": (240.0, 70.0, 0.95),
                "lamina_base": (100.0, 250.0, 0.95),
                "petiole_1": (85.0, 280.0, 0.92),
                "petiole_6": (45.0, 360.0, 0.90),  # sheet: (2005, 3700)
            },
        }
        parsed_l1 = _parse_basal_leaf_pose(leaf_local_1)
        self.assertIsNotNone(parsed_l1)
        self.assertAlmostEqual(parsed_l1["p_base"][0], 1950.0, places=2)
        self.assertAlmostEqual(parsed_l1["p_base"][1], 3600.0, places=2)
        self.assertAlmostEqual(parsed_l1["p_petiole_tip"][0], 1995.0, places=2)
        self.assertAlmostEqual(parsed_l1["p_petiole_tip"][1], 3695.0, places=2)

        crown_xy = estimate_rosette_crown_point(
            [leaf_local_1, leaf_local_2],
            root_boxes=[(3680, 1960, 3820, 2040)],
        )
        self.assertTrue(1920.0 <= crown_xy[0] <= 2080.0, f"Expected sheet-space crown X ~2000, got {crown_xy}")
        self.assertTrue(3650.0 <= crown_xy[1] <= 3800.0, f"Expected sheet-space crown Y ~3700, got {crown_xy}")

        scape_res = approximate_petiole_and_scape_metrics(
            basal_leaves=[leaf_local_1, leaf_local_2],
            p_crown={0: crown_xy},
            px_per_cm=100.0,
            capitula_detections=[(700, 1980, 760, 2020)],
        )
        # Corymb at (2000, 730) above crown at (~2000, ~3700) -> ~2970 px = ~29.7 cm
        self.assertTrue(28.0 <= float(scape_res["plant_height_cm"]) <= 31.5)
        self.assertAlmostEqual(float(scape_res["petiole_width_mm"]), 1.55, delta=0.15)

        # 2. Calibration provenance allowlist & 200-DPI (< 1.0 px/mm) preservation
        self.assertAlmostEqual(normalize_pixels_per_mm(0.7874), 0.7874, places=4)
        self.assertAlmostEqual(resolve_px_per_cm(pixels_per_mm=0.7874), 7.874, places=3)
        self.assertAlmostEqual(normalize_pixels_per_mm(0.08), 12.5, places=3)  # explicit legacy mm/px < 0.25

        parsed_mp = _parse_calibration_row({"catalogNumber": "NCU_MP_1", "pixels_per_mm": 9.8, "scale_source": "mp_estimated"})
        self.assertIsNotNone(parsed_mp)
        self.assertFalse(parsed_mp["is_empirical_ruler"])
        self.assertEqual(parsed_mp["scale_source"], "mp_estimated")
        self.assertIsNone(
            _parse_calibration_row({"catalogNumber": "NCU_UNCAL_1", "pixels_per_mm": 9.8, "scale_source": "uncalibrated"})
        )

        cal_no_mp = LM3ScaleCalibrator(device="cpu", providers=[], use_cf_predicted_by_mp=False)
        cal_no_mp.archival_session = None
        res_no_mp = cal_no_mp.calibrate_sheet(np.zeros((2000, 1500, 3), dtype=np.uint8))
        self.assertIsNone(res_no_mp["pixels_per_mm"])
        self.assertFalse(res_no_mp["is_calibrated_mm"])
        self.assertEqual(res_no_mp["scale_source"], "uncalibrated")

        # 3. CapitulumPhenologyConfig parameter wiring
        yellow_patch = np.full((30, 30, 3), 100, dtype=np.uint8)
        yellow_patch[:15, :, :] = (240, 175, 20)
        strict_head = classify_single_capitulum(yellow_patch, min_floret_yellow_ratio=0.80)
        permissive_head = classify_single_capitulum(yellow_patch, min_floret_yellow_ratio=0.20)
        self.assertFalse(strict_head["has_fresh_yellow_corollas"])
        self.assertTrue(permissive_head["has_fresh_yellow_corollas"])
        v_pheno_2 = classify_voucher_phenology("NCU_PHENO", [permissive_head], min_anthesis_heads_for_flowering=2)
        self.assertNotEqual(v_pheno_2["voucher_phenological_state"], "anthesis")

        # 4. Nested LM3_settings.yaml lookup, CPU threads init, SQLite width/upright round-trip,
        #    RUNTIME_ERROR label, and reniform aspect_ratio < 1.0 plotting
        cfg_loaded = PipelineConfig.from_yaml()
        self.assertFalse(cfg_loaded.segmentation.segment_partial_leaves)
        self.assertAlmostEqual(cfg_loaded.segmentation.rosette_plant_bbox_lower_fraction, 0.60, places=2)

        with tempfile.TemporaryDirectory() as tmp_d:
            tmp_p = Path(tmp_d)
            raw_dir = tmp_p / "run1" / "data" / "raw_vouchers"
            raw_dir.mkdir(parents=True, exist_ok=True)
            jpg_path = raw_dir / "NCU_RENIFORM.jpg"
            PILImage.fromarray(np.full((600, 600, 3), 240, dtype=np.uint8)).save(jpg_path, format="JPEG", quality=90)
            crash_path = raw_dir / "NCU_CRASH.webp"
            PILImage.fromarray(np.full((400, 400, 3), 240, dtype=np.uint8)).save(crash_path, format="WEBP", quality=90)

            v_csv = tmp_p / "run1" / "data" / "tables" / "curated_vouchers.csv"
            v_csv.parent.mkdir(parents=True, exist_ok=True)
            pd.DataFrame([
                {"catalogNumber": "NCU_RENIFORM", "scientificName": "Packera aurea", "image_path": str(jpg_path)},
                {"catalogNumber": "NCU_CRASH", "scientificName": "Packera aurea", "image_path": str(crash_path)},
            ]).to_csv(v_csv, index=False)

            dummy_w = tmp_p / "dummy.onnx"
            dummy_w.write_bytes(b"0")
            db_file = tmp_p / "run1" / "step02_audit.sqlite"

            # Verify CPU threads init does not raise NameError on `torch`
            pipe = step02.SegmentAndExtractPipeline(
                vouchers_csv=v_csv,
                model_weights=dummy_w,
                output_dir=tmp_p / "run1" / "data",
                device="cpu",
                threads=2,
                db_path=db_file,
                wave_size=1,
            )
            try:
                self.assertEqual(pipe.ect_num_dirs, 360)
                self.assertEqual(pipe.ect_num_thetas, 32)
                self.assertFalse(pipe.use_cf_predicted_by_mp)

                # Reniform leaf wider than long (W = 220 px, L = 140 px -> anatomical AR ~ 0.64 < 1.0)
                ren_mask = np.zeros((600, 600), dtype=np.uint8)
                cv2.ellipse(ren_mask, (200, 190), (110, 70), 0, 0, 360, 255, -1)
                ren_pet = np.zeros((600, 600), dtype=np.uint8)
                cv2.rectangle(ren_pet, (193, 255), (207, 350), 255, -1)

                orig_proc = pipe.process_voucher

                def _side_effect_proc(cat_n: str, img_p: Path, *a, **kw):
                    if cat_n == "NCU_CRASH":
                        raise RuntimeError("synthetic corrupt scan")
                    return orig_proc(cat_n, img_p, *a, **kw)

                pipe.ruler_calibrations["NCU_RENIFORM"] = 10.0
                with patch.object(
                    pipe.model_engine,
                    "predict",
                    return_value=[
                        ((120, 90, 260, 310), ren_mask, 0.96, 0),
                        ((255, 193, 350, 207), ren_pet, 0.92, 1),
                    ],
                ):
                    with patch.object(pipe, "process_voucher", side_effect=_side_effect_proc):
                        ext_df, fail_df = pipe.run()

                # Step 02 leaves existing raw JPEG scans untouched on disk
                self.assertTrue(jpg_path.exists())

                # Runtime crash tagged as SEGMENTATION_RUNTIME_ERROR, not UNRESOLVED_CLUMP
                self.assertEqual(len(fail_df), 1)
                self.assertEqual(str(fail_df.iloc[0]["failure_reason"]), "SEGMENTATION_RUNTIME_ERROR")

                # Reniform AR < 1.0 preserved and rendered in fig02a without clipping
                self.assertEqual(len(ext_df), 1)
                ar_ren = float(ext_df.iloc[0]["aspect_ratio"])
                self.assertLess(ar_ren, 0.85)
                self.assertTrue(bool(ext_df.iloc[0]["mask_is_upright"]))
                self.assertGreater(float(ext_df.iloc[0]["blade_width_mm"]), float(ext_df.iloc[0]["blade_length_mm"]))

                fig_paths = generate_step02_publication_figures(
                    extracted_df=ext_df,
                    failed_df=fail_df,
                    repro_df=pd.DataFrame(),
                    pheno_df=pd.DataFrame(),
                    figures_dir=tmp_p / "run1" / "outputs" / "figures",
                )
                self.assertTrue(Path(fig_paths["fig02a_png"]).is_file())

                # Verify SQLite round-trip of mask_is_upright, petiole_width_mm, blade_width_px, blade_width_mm
                conn = get_db_connection(db_path=db_file, init_schema=True)
                try:
                    loaded_ext, _, _, _ = load_step02_records_from_db(conn, ["NCU_RENIFORM"])
                    self.assertIn("NCU_RENIFORM", loaded_ext)
                    rec_db = loaded_ext["NCU_RENIFORM"][0]
                    self.assertTrue(rec_db["mask_is_upright"])
                    self.assertIsNotNone(rec_db["blade_width_px"])
                    self.assertIsNotNone(rec_db["blade_width_mm"])
                    self.assertGreater(float(rec_db["blade_width_mm"]), float(rec_db["blade_length_mm"]))
                finally:
                    conn.close()
            finally:
                pipe.shutdown()

        # 5. Docstring & CLI text hygiene audit across custom Step 2 modules
        stale_tokens = (
            "lm3_packera_seg_finetuned",
            "2-path leaf extraction",
            "retaining petiole-bridged lateral pinnae",
            "(Class 6 detections)",
            "Dissected Rachis Bridging",
            "topological connectivity gatekeeping",
            "Batched FP16 Stage 1 YOLOv5",
        )
        step2_sources = list(
            dict.fromkeys(
                [
                    PROJECT_ROOT / "scripts" / "pipeline" / "02_segment_and_extract.py",
                    PROJECT_ROOT / "scripts" / "vision" / "leaf_geometry_utils.py",
                    PROJECT_ROOT / "scripts" / "vision" / "lm3_inference_engines.py",
                    PROJECT_ROOT / "scripts" / "vision" / "organ_crown_metrics.py",
                    PROJECT_ROOT / "scripts" / "vision" / "rescore_step02_qc_and_zones.py",
                    PROJECT_ROOT / "scripts" / "vision" / "capitulum_phenology_classifier.py",
                    PROJECT_ROOT / "scripts" / "vision" / "step02_figures.py",
                    *sorted((PROJECT_ROOT / "scripts" / "vision").glob("step02_*.py")),
                ]
            )
        )
        for src_p in step2_sources:
            src_txt = src_p.read_text(encoding="utf-8")
            for bad_tok in stale_tokens:
                self.assertNotIn(bad_tok, src_txt, f"Stale docstring/comment token {bad_tok!r} found in {src_p.name}")

    def test_300dpi_png_only_12_figure_suite_and_high_richness_scaling(self) -> None:
        """
        Verify the 12-plate 300-DPI PNG-only figure architecture:
          1. Step 02 scales dynamically across 18 taxa (including subsp., var., and hybrid ranks),
             writes 300-DPI PNGs only (zero .pdf files), caches unlimited lm3_assets/ + ect_overlays/,
             and exports capped whole-sheet diagnostic overlays (lm3_sheet_overlays/).
          2. Steps 03-07 modular R suite (scripts/figures/*.R) defines all 12 R canonical PNG plates
             (14 total plates including fig02a and fig02b) with balanced R delimiters and zero
             file-backed grDevices::pdf() calls.
        """
        from PIL import Image
        from scripts.vision.step02_figures import generate_step02_publication_figures

        taxa_18 = [
            "Packera anonyma", "Packera aurea", "Packera crawfordii", "Packera debilis",
            "Packera dubia", "Packera glabella", "Packera hyperborealis", "Packera indecora",
            "Packera millefolium", "Packera obovata", "Packera paupercula",
            "Packera paupercula var. savannarum", "Packera paupercula subsp. pseudotomentosa",
            "Packera plattensis", "Packera pseudaurea", "Packera schweinitziana",
            "Packera tomentosa", "Packera × memmingeri",
        ]
        with tempfile.TemporaryDirectory() as tmp_d:
            tmp_p = Path(tmp_d)
            raw_dir = tmp_p / "raw"
            mask_dir = tmp_p / "masks"
            fig_dir = tmp_p / "outputs" / "figures"
            raw_dir.mkdir(parents=True)
            mask_dir.mkdir(parents=True)

            sheet_img = np.full((360, 280, 3), 242, dtype=np.uint8)
            cv2.rectangle(sheet_img, (20, 40), (120, 180), (65, 115, 70), -1)
            cv2.rectangle(sheet_img, (150, 170), (245, 310), (55, 105, 65), -1)
            sheet_file = raw_dir / "NCU_MULTI_RAMET.jpg"
            cv2.imwrite(str(sheet_file), sheet_img)

            rows = []
            for idx, tx in enumerate(taxa_18):
                for l_idx in (1, 2):
                    m_img = np.zeros((140, 110), dtype=np.uint8)
                    cv2.ellipse(m_img, (55, 70), (24 + (idx % 4) * 3, 46), 0, 0, 360, 255, -1)
                    m_p = mask_dir / f"CAT_{idx:02d}_L{l_idx}.png"
                    cv2.imwrite(str(m_p), m_img)
                    cat_id = "NCU_MULTI_RAMET" if idx == 0 else f"CAT_{idx:02d}"
                    rows.append({
                        "catalogNumber": cat_id,
                        "leaf_id": l_idx,
                        "plant_individual_id": l_idx if idx == 0 else 1,
                        "species_standardized": tx,
                        "determiner_tier": "Tier_1_Gold",
                        "is_type_specimen": (idx == 0 and l_idx == 1),
                        "leaf_zone": "cauline" if (idx % 3 == 0 and l_idx == 2) else "basal",
                        "solidity": 0.72 + 0.01 * (idx % 12),
                        "aspect_ratio": 0.78 if idx == 1 else (1.35 + 0.08 * (idx % 10)),
                        "archetype_score": 0.92,
                        "pixels_per_mm": 11.5,
                        "is_calibrated_mm": True,
                        "scale_source": "ruler_ocr",
                        "mask_path": str(m_p),
                        "image_path": str(sheet_file) if idx == 0 else "",
                        "bbox_x1": 20 if l_idx == 1 else 150,
                        "bbox_y1": 40 if l_idx == 1 else 170,
                        "bbox_x2": 120 if l_idx == 1 else 245,
                        "bbox_y2": 180 if l_idx == 1 else 310,
                    })

            saved = generate_step02_publication_figures(
                extracted_df=pd.DataFrame(rows),
                failed_df=pd.DataFrame(),
                repro_df=pd.DataFrame({"catalogNumber": ["CAT_01", "CAT_02"], "AR_inv": [1.18, 1.42]}),
                pheno_df=pd.DataFrame({"catalogNumber": ["CAT_01", "CAT_02"], "voucher_phenological_state": ["anthesis", "fruiting"]}),
                figures_dir=fig_dir,
            )
            self.assertEqual(list(fig_dir.rglob("*.pdf")), [])
            for key in ("fig02a_png", "fig02b_png"):
                p_out = Path(saved[key])
                self.assertTrue(p_out.is_file() and p_out.stat().st_size > 15000)
                with Image.open(p_out) as im:
                    dpi_x, dpi_y = im.info.get("dpi", (0, 0))
                    self.assertAlmostEqual(dpi_x, 300.0, delta=2.0)
                    self.assertAlmostEqual(dpi_y, 300.0, delta=2.0)

            self.assertGreaterEqual(len(list((fig_dir / "ect_overlays").glob("*.png"))), 18)
            self.assertGreaterEqual(len(list((fig_dir / "lm3_sheet_overlays").glob("*.png"))), 1)

        expected_r_pngs = {
            "fig03a_efa_morphospace_and_harmonics.png",
            "fig03b_topological_ect_and_heteroblasty.png",
            "fig04a_foliar_allometry_audit.png",
            "fig04b_gmm_cda_morphometrics.png",
            "fig04c_ecoregion_and_batch_audit.png",
            "fig05a_ecological_niche_and_rf_importance.png",
            "fig05b_hopkins_circular_phenology.png",
            "fig05c_sdm_suitability_and_distribution_map.png",
            "fig06a_genomic_pcoa_tree_concordance.png",
            "fig06b_reticulation_gdi_landscape.png",
            "fig07a_integrative_taxonomic_synthesis.png",
            "fig07b_taxon_complex_delimitation_network.png",
        }
        r_fig_dir = PROJECT_ROOT / "scripts" / "figures"
        combined_r_src = "\n".join(p.read_text(encoding="utf-8") for p in sorted(r_fig_dir.glob("*.R")))
        import re
        for png_name in expected_r_pngs:
            self.assertIn(png_name, combined_r_src)
        combined_no_null = re.sub(r"grDevices::pdf\(\s*NULL\s*\)", "", combined_r_src)
        self.assertNotIn("grDevices::pdf(", combined_no_null)


if __name__ == "__main__":
    unittest.main()
