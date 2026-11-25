#!/usr/bin/env python3
"""Warm-R tests for Step 03 basal gating, SQLite synonym resolution, 32x360 ECT handoff, and void QC."""

from __future__ import annotations

import os
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np
import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from scripts.tests._r_session import rscript_available, run_warm_r_cmd


def _read_db_utils_bundle() -> str:
    paths = [PROJECT_ROOT / "scripts/core/db_utils.R", *sorted((PROJECT_ROOT / "scripts/core").glob("db_r_*.R"))]
    return "\n".join(p.read_text(encoding="utf-8") for p in paths if p.exists())


def _read_dissected_utils_bundle() -> str:
    paths = [PROJECT_ROOT / "scripts/morphometrics/dissected_leaf_utils.R", *sorted((PROJECT_ROOT / "scripts/morphometrics").glob("dissect_*.R"))]
    return "\n".join(p.read_text(encoding="utf-8") for p in paths if p.exists())


class TestStep03VoidQCAndSQLite(unittest.TestCase):
    """Warm-R tests for Step 03 basal gating, SQLite synonym resolution, 32x360 ECT handoff, and void QC."""

    def setUp(self) -> None:
        self.r_script_path = PROJECT_ROOT.joinpath("scripts", "morphometrics", "03_fourier_extractor.R")

    def test_basal_positional_gate_lowest_lamina_fallback_ribbon_filter_and_otu_medians(self):
        """Verify basal positional homology filter, lowest-lamina fallback, and OTU median/q25 dissection summaries."""
        if not rscript_available():
            self.skipTest("Rscript not installed.")
        t = np.linspace(0, 2 * np.pi, 120, endpoint=False)

        def write_xy_csv(path: Path, cat: str, lid: int, xy: np.ndarray):
            pd.DataFrame({
                "catalogNumber": cat,
                "leaf_id": lid,
                "point_index": np.arange(len(xy)),
                "x": xy[:, 0] + 250.0,
                "y": xy[:, 1] + 250.0,
            }).to_csv(path, index=False)

        entire_xy = np.column_stack([85.0 * np.cos(t), 52.0 * np.sin(t)])
        u_t = np.cos(t)
        lobe_env = (0.25 + 0.75 * ((u_t + 1.0) / 2.0) ** 0.6) * (0.35 + 0.65 * np.abs(np.sin(4.5 * np.pi * (u_t + 1.0) / 2.0)))
        pinn_xy = np.column_stack([100.0 * u_t, 42.0 * np.sin(t) * lobe_env])

        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_p = Path(tmpdir)
            contours_dir = tmp_p / "contours"
            contours_dir.mkdir(parents=True, exist_ok=True)
            vouchers_csv = tmp_p / "curated_vouchers.csv"
            manifest_csv = tmp_p / "extracted_leaf_manifest.csv"

            vouchers = [
                {"catalogNumber": "V_MIX", "taxon": "Packera paupercula", "scientificName": "Packera paupercula", "species_raw": "Packera paupercula", "determiner_tier": "Tier_1_Gold"},
                {"catalogNumber": "V_DET", "taxon": "Packera paupercula", "scientificName": "Packera paupercula", "species_raw": "Packera paupercula", "determiner_tier": "Tier_1_Gold"},
                {"catalogNumber": "V_REF1", "taxon": "Packera obovata", "scientificName": "Packera obovata", "species_raw": "Packera obovata", "determiner_tier": "Tier_1_Gold"},
                {"catalogNumber": "V_REF2", "taxon": "Packera obovata", "scientificName": "Packera obovata", "species_raw": "Packera obovata", "determiner_tier": "Tier_1_Gold"},
            ]
            pd.DataFrame(vouchers).to_csv(vouchers_csv, index=False)

            leaf_specs = [
                # V_MIX: leaf 1 is upper cauline (cy=0.18 -> excluded); leaf 2 entire basal; leaves 3-4 pinnatifid basal
                ("V_MIX", 1, entire_xy, 0.18, 0.94, False),
                ("V_MIX", 2, entire_xy, 0.72, 0.92, False),
                ("V_MIX", 3, pinn_xy, 0.78, 0.44, True),
                ("V_MIX", 4, pinn_xy * 1.03, 0.82, 0.46, True),
                # V_DET: top-mounted detached rosette (all cy < 0.40, rescued by max - 0.15 fallback)
                ("V_DET", 1, entire_xy * 0.98, 0.22, 0.91, False),
                ("V_DET", 2, pinn_xy * 0.97, 0.31, 0.45, True),
                # Reference vouchers
                ("V_REF1", 1, entire_xy * 1.02, 0.76, 0.93, False),
                ("V_REF2", 1, entire_xy * 0.96, 0.79, 0.92, False),
            ]
            manifest_rows = []
            for cat, lid, xy_arr, cy_p, sol_v, is_diss in leaf_specs:
                cpath = contours_dir / f"{cat}_leaf{lid}.csv"
                write_xy_csv(cpath, cat, lid, xy_arr)
                manifest_rows.append({
                    "catalogNumber": cat, "plant_individual_id": 0, "leaf_id": lid, "pixels_per_mm": 10.0,
                    "is_calibrated_mm": True, "solidity": sol_v, "archetype_score": 0.93, "is_dissected": is_diss,
                    "cy_rel_plant_bbox": cy_p, "cy_rel_sheet": cy_p, "area_mm2": 200.0, "contour_path": str(cpath),
                })
            pd.DataFrame(manifest_rows).to_csv(manifest_csv, index=False)

            out_spec = tmp_p / "leaf_efa_harmonics.csv"
            out_ind = tmp_p / "leaf_efa_harmonics_ind.csv"
            proc = run_warm_r_cmd(
                ["Rscript", str(self.r_script_path), "--contours-dir", str(contours_dir), "--masks-dir", str(tmp_p / "masks"), "--vouchers", str(vouchers_csv), "--manifest", str(manifest_csv), "--output", str(out_spec), "--output-individual", str(out_ind)],
                stage_key="step03", cwd=PROJECT_ROOT, skip_figures=True,
            )
            self.assertEqual(proc.returncode, 0, f"Step 03 basal/dissection test failed:\n{proc.stderr}\n{proc.stdout}")
            ind_df = pd.read_csv(out_ind)
            spec_df = pd.read_csv(out_spec)

            vmix_rows = ind_df.loc[ind_df["catalogNumber"] == "V_MIX"]
            if "leaf_zone" in vmix_rows.columns and len(vmix_rows) == 4:
                vmix_basal_lids = sorted(
                    vmix_rows.loc[vmix_rows["leaf_zone"].astype(str).str.lower() == "basal", "leaf_id"]
                    .astype(int)
                    .tolist()
                )
                vmix_caul_lids = sorted(
                    vmix_rows.loc[vmix_rows["leaf_zone"].astype(str).str.lower() == "cauline", "leaf_id"]
                    .astype(int)
                    .tolist()
                )
                self.assertEqual(vmix_basal_lids, [2, 3, 4], f"Expected basal leaves [2, 3, 4], got {vmix_basal_lids}")
                self.assertEqual(vmix_caul_lids, [1], f"Expected cauline leaf [1], got {vmix_caul_lids}")
            else:
                vmix_lids = sorted(vmix_rows["leaf_id"].astype(int).tolist())
                self.assertEqual(vmix_lids, [2, 3, 4], f"Expected basal leaves [2, 3, 4], got {vmix_lids}")
            vdet_lids = sorted(ind_df.loc[ind_df["catalogNumber"] == "V_DET", "leaf_id"].astype(int).tolist())
            self.assertEqual(vdet_lids, [1, 2], f"Expected detached rosette leaves rescued by fallback, got {vdet_lids}")

            vmix_otu = spec_df.loc[spec_df["catalogNumber"] == "V_MIX"].iloc[0]
            self.assertEqual(int(vmix_otu["basal_leaf_count"]), 3)
            # True median of [0.92, 0.44, 0.46] is 0.46 (NOT first-leaf 0.92)
            self.assertAlmostEqual(float(vmix_otu["solidity"]), 0.46, places=2)
            self.assertLess(float(vmix_otu["q25_solidity"]), 0.50)
            self.assertAlmostEqual(float(vmix_otu["dissected_leaf_fraction"]), 2.0 / 3.0, places=3)
            self.assertGreater(float(vmix_otu["dissection_index"]), 1.35)

    def test_dynamic_sqlite_synonym_resolution_and_tier_b_cohort_hygiene(self):
        """Verify dynamic SQLite taxon_synonyms lookup in R and Tier B orphan/synonym refresh in populate_run_active_cohort."""
        from scripts.core.db import get_db_connection, init_db_schema, populate_run_active_cohort

        with tempfile.TemporaryDirectory() as tmp:
            eph_db = Path(tmp) / "ephemeral_synonym_test.sqlite"
            conn = get_db_connection(eph_db)
            init_db_schema(conn)
            with conn:
                conn.executemany(
                    "INSERT INTO taxon_synonyms (raw_name_lower, canonical_packera) VALUES (?, ?);",
                    [
                        ("packera synalpha", "Packera canonbeta"),
                        ("packera synalpha c.auth", "Packera canonbeta"),
                        ("senecio synalpha", "Packera canonbeta"),
                    ],
                )
                conn.executemany(
                    "INSERT INTO vouchers (catalogNumber, scientificName, species_standardized, determiner_tier, qc_status) VALUES (?, ?, ?, ?, 'PASS');",
                    [
                        ("CAT_KEEP", "Packera synalpha C.Auth", "Packera canonbeta", "Tier_1_Gold"),
                        ("CAT_DROP", "Packera gammataxon", "Packera gammataxon", "Tier_2_Silver"),
                    ],
                )
                conn.executemany(
                    "INSERT INTO plant_individuals (catalogNumber, plant_individual_id) VALUES (?, ?);",
                    [("CAT_KEEP", 1), ("CAT_DROP", 1)],
                )
            populate_run_active_cohort(
                conn,
                run_id="test_run",
                voucher_rows=[
                    {"catalogNumber": "CAT_KEEP", "species_standardized": "Packera synalpha", "determiner_tier": "Tier_1_Gold"},
                    {"catalogNumber": "CAT_DROP", "species_standardized": "Packera gammataxon", "determiner_tier": "Tier_2_Silver"},
                ],
            )
            with conn:
                conn.executemany(
                    "INSERT INTO otu_morphometrics (run_id, catalogNumber, plant_individual_id, species_standardized) VALUES (?, ?, ?, ?);",
                    [
                        ("test_run", "CAT_KEEP", 1, "Packera synalpha"),
                        ("test_run", "CAT_DROP", 1, "Packera gammataxon"),
                    ],
                )
            # Re-populate narrowed cohort with canonical standardized species
            populate_run_active_cohort(
                conn,
                run_id="test_run",
                voucher_rows=[
                    {"catalogNumber": "CAT_KEEP", "species_standardized": "Packera canonbeta", "determiner_tier": "Tier_1_Gold"},
                ],
            )
            rows_om = conn.execute(
                "SELECT catalogNumber, plant_individual_id, species_standardized FROM otu_morphometrics WHERE run_id = 'test_run';"
            ).fetchall()
            conn.close()
            self.assertEqual([(r[0], int(r[1]), r[2]) for r in rows_om], [("CAT_KEEP", 1, "Packera canonbeta")])

            if rscript_available():
                env = dict(os.environ)
                env["PLANTMATICS_DB_PATH"] = str(eph_db)
                r_expr = (
                    "if (grDevices::dev.cur() == 1L) grDevices::pdf(NULL); "
                    "source('scripts/morphometrics/pub_figures_utils.R'); "
                    "res <- standardize_taxon_name(c('Packera synalpha', 'Packera synalpha C.Auth', 'P. synalpha', 'Senecio synalpha', 'Packera canonbeta var. minor')); "
                    "cat(paste(res, collapse='|'))"
                )
                proc = run_warm_r_cmd(["Rscript", "-e", r_expr], stage_key='step03', cwd=PROJECT_ROOT, env=env, skip_figures=False)
                self.assertEqual(proc.returncode, 0, f"R synonym test failed:\n{proc.stderr}\n{proc.stdout}")
                self.assertEqual(
                    proc.stdout.strip().splitlines()[-1],
                    "|".join(["Packera canonbeta"] * 5),
                )

    def test_step03_rigid_tipup_dissected_homology_reniform_lw_and_32x360_ect(self) -> None:
        """
        Verify Step 02 -> Step 03 rigid (-dy, -dx) contour handoff, reniform L/W < 1.0 preservation,
        mp_estimated centroid_size gating, and dynamic 32x360 .npy ECT ingestion.
        """
        # 1. Verify fourier_parallel_utils.R and dissected_leaf_utils.R static contracts
        fpu_content = (PROJECT_ROOT / "scripts" / "morphometrics" / "fourier_parallel_utils.R").read_text(encoding="utf-8")
        dlu_content = _read_dissected_utils_bundle()
        dbu_content = _read_db_utils_bundle()
        self.assertIn("is_step02_tipup", fpu_content)
        self.assertIn("cbind(-dy + cx, -dx + cy)", fpu_content)
        self.assertNotIn("floor(nrow(c_xy) / 2L) + 1L", fpu_content)
        self.assertIn("524288L", dlu_content)
        self.assertIn(".resample_ect_directions_r", dlu_content)
        self.assertIn("mp_estimated", dbu_content)

        # Mathematical verification: asymmetrically pinnatifid Step 02 tip-up leaf (extra lobes on left margin)
        # preserves exact 0.000-deg midrib alignment under (-dy + cx, -dx + cy) whereas vertex-61 drifts off-axis
        t = np.linspace(0.0, 2.0 * np.pi, 120, endpoint=False)
        # Clockwise in image space (x_right, y_down): t=0 base (250, 340), t=pi/2 left (x<250), t=pi tip (250, 160), t=3pi/2 right (x>250)
        y_img = 250.0 + 90.0 * np.cos(t)
        left_mask = np.sin(t) > 0.0
        lobe_asym = np.where(left_mask, 1.0 + 0.42 * np.sin(5.0 * t) ** 2, 1.0)
        x_img = 250.0 - 38.0 * np.sin(t) * lobe_asym
        # Resample by equal perimeter arc-length (just like Step 02 resample_closed) so asymmetric left perimeter shifts index 60 off tip!
        xy_raw = np.column_stack([x_img, y_img])
        closed = np.vstack([xy_raw, xy_raw[0:1]])
        seg = np.hypot(np.diff(closed[:, 0]), np.diff(closed[:, 1]))
        cum = np.concatenate([[0.0], np.cumsum(seg)])
        targets = np.linspace(0.0, cum[-1], 120, endpoint=False)
        coo_step02 = np.column_stack([np.interp(targets, cum, closed[:, 0]), np.interp(targets, cum, closed[:, 1])])

        # Old vertex-60 heuristic drifts significantly off the vertical midrib x=250:
        old_drift_deg = abs(np.degrees(np.arctan2(coo_step02[60, 0] - coo_step02[0, 0], coo_step02[0, 1] - coo_step02[60, 1])))
        self.assertGreater(old_drift_deg, 3.0)

        # New Step 03 rigid transform (-dy + cx, -dx + cy) maps true midrib (x_img == 250) onto constant Y == cy_m (0.000-deg drift)
        cx_m, cy_m = float(coo_step02[:, 0].mean()), float(coo_step02[:, 1].mean())
        dx_m, dy_m = coo_step02[:, 0] - cx_m, coo_step02[:, 1] - cy_m
        momocs_xy = np.column_stack([-dy_m + cx_m, -dx_m + cy_m])
        # Clockwise in Cartesian Momocs space (shoelace signed area < 0)
        shoelace = 0.5 * np.sum(momocs_xy[:, 0] * np.roll(momocs_xy[:, 1], -1) - np.roll(momocs_xy[:, 0], -1) * momocs_xy[:, 1])
        self.assertLess(shoelace, 0.0)
        s_03_chk = (PROJECT_ROOT / "scripts" / "morphometrics" / "03_fourier_extractor.R").read_text(encoding="utf-8")
        s_db_chk = _read_db_utils_bundle()
        self.assertIn("allow_multi_component_ect = TRUE", s_03_chk)
        self.assertIn("allow_multi_component_ect = TRUE", s_db_chk)

        if not rscript_available():
            self.skipTest("Rscript not installed.")

        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_p = Path(tmpdir)
            contours_dir = tmp_p / "contours"
            ect_dir = tmp_p / "ect"
            contours_dir.mkdir(parents=True, exist_ok=True)
            ect_dir.mkdir(parents=True, exist_ok=True)
            vouchers_csv = tmp_p / "curated_vouchers.csv"
            manifest_csv = tmp_p / "extracted_leaf_manifest.csv"

            vouchers = []
            manifest_rows = []
            for idx in range(6):
                cat = f"STEP23_HO_{idx:02d}"
                is_reniform = idx >= 3
                sp = "Packera aurea" if is_reniform else "Packera paupercula"
                if is_reniform:
                    # Reniform wider-than-long tip-up blade: L = 110 (y span), W = 160 (x span) -> L/W = 0.6875 < 1.0
                    y_c = 250.0 + 55.0 * np.cos(t)
                    x_c = 250.0 - 80.0 * np.sin(t)
                    lw_true = 110.0 / 160.0
                    sol_v = 0.93
                else:
                    y_c = coo_step02[:, 1]
                    x_c = coo_step02[:, 0]
                    lw_true = 2.15
                    sol_v = 0.64

                cpath = contours_dir / f"{cat}_leaf1.csv"
                pd.DataFrame({
                    "catalogNumber": cat,
                    "leaf_id": 1,
                    "point_index": np.arange(120),
                    "x": x_c,
                    "y": y_c,
                }).to_csv(cpath, index=False)

                epath = ect_dir / f"{cat}_leaf1_ect.npy"
                theta_grid = np.linspace(0.0, 2.0 * np.pi, 360, endpoint=False, dtype=np.float32)
                rad_grid = np.linspace(-1.0, 1.0, 32, dtype=np.float32)[:, None]
                ect_32x360 = (
                    np.exp(-(rad_grid ** 2) * (2.0 + (0.0 if is_reniform else 1.5 * np.cos(4.0 * theta_grid)[None, :])))
                ).astype(np.float32)
                np.save(str(epath), ect_32x360)

                is_mp = idx == 5
                vouchers.append({
                    "catalogNumber": cat,
                    "species_raw": sp,
                    "scientificName": sp,
                    "determiner_tier": "Tier_1_Gold",
                })
                manifest_rows.append({
                    "catalogNumber": cat, "plant_individual_id": 0, "leaf_id": 1, "pixels_per_mm": 10.5 if is_mp else 10.0,
                    "is_calibrated_mm": (not is_mp), "scale_source": "mp_estimated" if is_mp else "ruler_calibrated",
                    "aspect_ratio": round(lw_true, 4), "solidity": sol_v, "archetype_score": 0.94, "contour_path": str(cpath), "ect_path": str(epath),
                })

            pd.DataFrame(vouchers).to_csv(vouchers_csv, index=False)
            pd.DataFrame(manifest_rows).to_csv(manifest_csv, index=False)

            out_spec = tmp_p / "leaf_efa_harmonics.csv"
            out_ind = tmp_p / "leaf_efa_harmonics_ind.csv"
            proc = run_warm_r_cmd(
                ["Rscript", str(self.r_script_path), "--contours-dir", str(contours_dir), "--masks-dir", str(tmp_p / "masks"), "--vouchers", str(vouchers_csv), "--manifest", str(manifest_csv), "--output", str(out_spec), "--output-ind", str(out_ind), "--figures-dir", str(tmp_p / "figs"), "--harmonics", "12", "--num-pcs", "5"],
                stage_key="step03", cwd=PROJECT_ROOT, skip_figures=True,
            )
            self.assertEqual(proc.returncode, 0, f"Step 03 handoff Rscript test failed:\n{proc.stderr}\n{proc.stdout}")
            ind_df = pd.read_csv(out_ind)
            ren_rows = ind_df[ind_df["catalogNumber"].isin(["STEP23_HO_03", "STEP23_HO_04", "STEP23_HO_05"])]
            self.assertTrue((ren_rows["aspect_ratio"].astype(float) < 0.80).all(), "Reniform Packera aurea leaves must retain aspect_ratio < 1.0 in Step 03.")
            mp_row = ind_df[ind_df["catalogNumber"] == "STEP23_HO_05"].iloc[0]
            self.assertTrue(pd.isna(mp_row["centroid_size"]), "mp_estimated leaf must keep centroid_size == NA in Step 03.")
            self.assertTrue(ind_df["ECT_PC1"].notna().all(), "Dynamic 32x360 .npy ECT reader must populate ECT_PC1 for all leaves.")

    def test_step03_enclosed_void_efa_withholding_bilateral_qc_and_stray_pid_rezoning(self) -> None:
        """
        Verify Step 03 R contracts across `db_utils.R`, `dissected_leaf_utils.R`, and `03_fourier_extractor.R`:
          1. Schema version `20261121L` and 7 new `extracted_leaves` QC columns (`enclosed_void_frac`,
             `bilateral_lobe_symmetry`, `hole_area_frac`, `is_crop_clipped`, `detector_class_id`, `efa_qc_pass`, `qc_flag_reason`).
          2. `03_fourier_extractor.R` (`non_jordan_mask`) and `db_utils.R` withhold single-loop EFA harmonics when
             `enclosed_void_frac > 0.020` (`ENCLOSED_INTERLOBE_VOID`) while retaining the row for hole-aware `32x360` ECT descriptors.
          3. `evaluate_step03_lm3_efa_qc` enforces `bilateral_lobe_symmetry` (`is_verified_diss`) before granting
             `solidity < 0.72` asymmetry exemptions, and rejects `hole_area_frac` (`haf > 0.035`) and `is_crop_clipped`.
          4. `dissected_leaf_utils.R` tightens `d_crown[ii]` fallback (`<= 0.18` basal vs `> 0.24` cauline, plus `> 0.28`
             cauline-only stem guard) and `03_fourier_extractor.R` re-zones merged upper-stem stray PIDs (`rezoned_strays_s03`) to `cauline`.
        """
        dbu_txt = _read_db_utils_bundle()
        dlu_txt = _read_dissected_utils_bundle()
        s03_txt = (PROJECT_ROOT / "scripts" / "morphometrics" / "03_fourier_extractor.R").read_text(encoding="utf-8")

        self.assertIn("PLANTMATICS_SCHEMA_USER_VERSION <- 20261127L", dbu_txt)
        self.assertIn("PACKERA_SCHEMA_USER_VERSION <- PLANTMATICS_SCHEMA_USER_VERSION", dbu_txt)
        for col_name in ("detector_class_id", "is_crop_clipped", "hole_area_frac",
                         "enclosed_void_frac", "bilateral_lobe_symmetry", "efa_qc_pass", "qc_flag_reason"):
            self.assertIn(col_name, dbu_txt)

        self.assertIn("is_enc_void", dbu_txt)
        self.assertIn("ENCLOSED_INTERLOBE_VOID", dbu_txt)
        self.assertIn("is_verified_diss", dbu_txt)
        self.assertIn("haf > 0.035", dbu_txt)

        self.assertIn("d_crown[ii] <= 0.18", dlu_txt)
        self.assertIn("d_crown[ii] > 0.24", dlu_txt)
        self.assertIn("> 0.28", dlu_txt)
        self.assertIn("crown_verified_cauline_only", dlu_txt)

        self.assertIn("non_jordan_mask", s03_txt)
        self.assertIn("ENCLOSED_INTERLOBE_VOID", s03_txt)
        self.assertIn("enclosed_void_frac", s03_txt)
        self.assertIn("consolidate_fragmented_sheet_individuals", s03_txt)
        self.assertIn("min_prim_bas_cy - 0.12", s03_txt + "\n" + dlu_txt)
        self.assertIn("rezoned_strays_s03", s03_txt + "\n" + dlu_txt)

        if not rscript_available():
            return

        r_expr = (
            "if (grDevices::dev.cur() == 1L) grDevices::pdf(NULL); "
            "source('scripts/core/db_utils.R'); "
            "df_q <- data.frame("
            "  solidity = c(0.60, 0.60, 0.88, 0.88), "
            "  bilateral_dice = c(0.62, 0.62, 0.91, 0.91), "
            "  si_a = c(0.22, 0.22, 0.03, 0.03), "
            "  lm3_gates_pass = c(FALSE, FALSE, TRUE, TRUE), "
            "  bilateral_lobe_symmetry = c(FALSE, TRUE, TRUE, TRUE), "
            "  enclosed_void_frac = c(0.0, 0.0, 0.045, 0.0), "
            "  hole_area_frac = c(0.0, 0.0, 0.0, 0.050), "
            "  stringsAsFactors = FALSE"
            "); "
            "res_efa <- evaluate_step03_lm3_efa_qc(df_q, allow_multi_component_ect = FALSE); "
            "res_ect <- evaluate_step03_lm3_efa_qc(df_q, allow_multi_component_ect = TRUE); "
            "stopifnot(identical(as.logical(res_efa), c(FALSE, TRUE, FALSE, FALSE))); "
            "stopifnot(identical(as.logical(res_ect), c(FALSE, TRUE, TRUE, FALSE))); "
            "cat('STEP03_QC_GATES_OK\n')"
        )
        proc = run_warm_r_cmd(["Rscript", "-e", r_expr], stage_key="step03", cwd=PROJECT_ROOT)
        self.assertEqual(proc.returncode, 0, f"R QC gate check failed:\n{proc.stderr}\n{proc.stdout}")
        self.assertIn("STEP03_QC_GATES_OK", proc.stdout)


if __name__ == "__main__":
    unittest.main()
