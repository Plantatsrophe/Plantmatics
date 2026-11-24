#!/usr/bin/env python3
"""Batched warm-R contract suite for Step 03 taxonomic clarity invariants."""

from __future__ import annotations

import json, re, sqlite3, sys, tempfile, unittest
from pathlib import Path

import numpy as np
import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from scripts.core.db_schema import PLANTMATICS_SCHEMA_USER_VERSION, init_db_schema
from scripts.tests._r_session import rscript_available, run_warm_r_cmd

RUN_ID = "run_s03_clarity"


def _contour_xy(rx: float, ry: float, lobe: float = 0.0, asym: float = 0.0, shear: float = 0.0) -> np.ndarray:
    t = np.linspace(0.0, 2.0 * np.pi, 120, endpoint=False)
    mod = 1.0 + asym * np.cos(t) + lobe * np.cos(4.0 * t)
    x = 250.0 - rx * np.sin(t) * mod + shear * rx * (1.0 - np.cos(t))
    y = 250.0 + ry * np.cos(t) + shear * ry * np.sin(t)
    return np.column_stack([x, y])


def _col_diff(a: np.ndarray, b: np.ndarray) -> float:
    return min(float(np.max(np.abs(a - b))), float(np.max(np.abs(a + b))))


def _seed_sqlite_clarity_cohort(tmp_p: Path) -> Path:
    db_path, sql_dir = tmp_p / "step03_clarity.sqlite", tmp_p / "sql_run"
    for sub in ("contours", "tables", "figures"):
        (sql_dir / sub).mkdir(parents=True, exist_ok=True)
    conn = sqlite3.connect(str(db_path))
    init_db_schema(conn)
    conn.execute("INSERT INTO pipeline_runs (run_id) VALUES (?);", (RUN_ID,))
    specs = [
        ("V_ABAX_1", "Packera aurea", 64.0, 38.0, 0.04, 0.10, 0.10, False), ("V_AUR_2", "Packera aurea", 60.0, 41.0, 0.08, 0.04, 0.04, False),
        ("V_AUR_3", "Packera aurea", 56.0, 44.0, 0.02, -0.03, 0.02, False), ("V_DUB_1", "Packera dubia", 42.0, 72.0, 0.06, 0.03, 0.05, False),
        ("V_DUB_2", "Packera dubia", 45.0, 68.0, 0.10, -0.04, 0.03, False), ("V_DUB_3", "Packera dubia", 40.0, 75.0, 0.05, 0.02, 0.04, False),
        ("V_VOID_SQL", "Packera dubia", 38.0, 70.0, 0.22, 0.05, 0.06, True),
    ]
    cfg_keys = (
        "leaf_id", "is_abaxial_pressed_leaf", "blade_cie_l_norm", "blade_cie_l_sd", "blade_cie_l_sheet_delta",
        "blade_cie_l_q85_sheet_delta", "basal_zone_cie_l_delta", "blade_cie_a_delta", "blade_cie_b_delta",
        "abaxial_cie_l_sheet_delta", "base_angle_deg", "apex_angle_deg", "lamina_curvature", "lamina_trace_length_px",
        "basal_cordate_index", "margin_serration_index", "margin_tooth_acuteness", "terminal_lobe_area_frac",
        "bbox_x1", "bbox_y1", "bbox_x2", "bbox_y2",
    )
    base_l = {
        "plant_individual_id": 0, "circularity": 0.84, "petiole_length_mm": 12.0, "petiole_width_mm": 1.2,
        "blade_length_mm": 26.0, "blade_width_mm": 16.0, "petiole_to_blade_ratio": 0.46, "pixels_per_mm": 10.0,
        "is_calibrated_mm": 1, "scale_source": "ruler_calibrated", "mask_is_upright": 1, "leaf_zone": "basal",
        "leaf_zone_confidence": 0.95, "dist_from_crown_norm": 0.10, "cy_rel_plant_bbox": 0.82, "archetype_score": 0.95,
        "bilateral_dice": 0.90, "si_a": 0.05, "lm3_gates_pass": 1, "n_lamina_components": 1, "efa_qc_pass": 1, "qc_flag_reason": "PASS",
    }
    v_rows, p_rows, rac_rows, leaf_dicts = [], [], [], []
    for cat, sp, rx, ry, lobe, asym, shear, is_void in specs:
        v_rows.append((cat, sp, sp, "Tier_1_Gold", "", 1, 10.0, 1))
        p_rows.append((cat, 0))
        rac_rows.append((RUN_ID, cat, 0, sp, "Tier_1_Gold", "", 0, 1.0, 1.0))
        c_json = json.dumps([[round(float(px), 2), round(float(py), 2)] for px, py in _contour_xy(rx, ry, lobe, asym, shear)])
        cfgs = (
            [(1, 0, 0.42, 4.1, -19.0, -14.0, -16.5, 2.4, 11.2, None, 112.0, 78.0, 0.14, 420.0, 0.28, 0.34, 0.52, 0.41, 100, 120, 220, 260),
             (2, 1, 0.81, 2.8, -3.5, -2.0, -4.0, 0.5, 3.1, -3.5, 116.0, 82.0, 0.18, 460.0, 0.32, 0.38, 0.56, 0.45, 130, 130, 250, 270)]
            if cat == "V_ABAX_1" else
            [(1, 0, 0.46, 3.5, -15.0, -11.0, -13.0, 1.8, 9.0, None, 98.0, 72.0, 0.11, 380.0, 0.18, 0.26, 0.44, 0.36, 110, 110, 230, 250)]
        )
        for cfg in cfgs:
            row = {**base_l, **dict(zip(cfg_keys, cfg)), "catalogNumber": cat, "shape_id": f"{cat}_p0_leaf{cfg[0]}",
                   "solidity": 0.64 if is_void else 0.90, "aspect_ratio": round(ry / rx, 3), "area_px": float(np.pi * rx * ry),
                   "area_mm2": round(float(np.pi * rx * ry) / 100.0, 2), "contour_xy_json": c_json, "enclosed_void_frac": 0.045 if is_void else 0.0}
            leaf_dicts.append(row)
    conn.executemany("INSERT INTO vouchers (catalogNumber, scientificName, species_standardized, determiner_tier, type_status, is_primary_duplicate, pixels_per_mm, is_calibrated_mm) VALUES (?, ?, ?, ?, ?, ?, ?, ?);", v_rows)
    conn.executemany("INSERT INTO plant_individuals (catalogNumber, plant_individual_id) VALUES (?, ?);", p_rows)
    conn.executemany("INSERT INTO run_active_cohort (run_id, catalogNumber, plant_individual_id, species_standardized, determiner_tier, type_status, is_type_specimen, authority_weight, spatial_balance_weight) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?);", rac_rows)
    l_cols = list(leaf_dicts[0].keys())
    conn.executemany(f"INSERT INTO extracted_leaves ({', '.join(l_cols)}) VALUES ({', '.join(['?'] * len(l_cols))});", [tuple(r[c] for c in l_cols) for r in leaf_dicts])
    conn.commit()
    conn.close()
    pd.DataFrame([{"catalogNumber": r[0], "scientificName": r[1], "species_standardized": r[2], "determiner_tier": r[3],
                   "type_status": r[4], "is_primary_duplicate": 1, "pixels_per_mm": r[6], "is_calibrated_mm": r[7]} for r in v_rows]).to_csv(sql_dir / "tables" / "curated_vouchers.csv", index=False)
    return db_path


def _seed_ect_cohort(base_dir: Path, rep_g1: int = 1, include_void: bool = True, include_bronze: bool = False) -> None:
    cdir = base_dir / "contours"
    cdir.mkdir(parents=True, exist_ok=True)
    v_rows, m_rows = [], []
    items = [
        ("G1", "Tier_1_Gold", 64.0, 38.0, 0.00, 0.90, 0.0, 1), ("G2", "Tier_1_Gold", 60.0, 41.0, 0.06, 0.88, 0.0, 1),
        ("G3", "Tier_1_Gold", 56.0, 44.0, 0.12, 0.85, 0.0, 1), ("G4", "Tier_1_Gold", 52.0, 48.0, 0.04, 0.89, 0.0, 1),
        ("G5", "Tier_1_Gold", 48.0, 52.0, 0.10, 0.86, 0.0, 1), ("G6", "Tier_1_Gold", 44.0, 56.0, 0.16, 0.83, 0.0, 1),
    ]
    if include_void:
        items.append(("G_VOID", "Tier_1_Gold", 28.0, 86.0, 0.28, 0.62, 0.045, 2))
    if include_bronze:
        items.append(("B_EXTREME", "Tier_3_Bronze", 18.0, 96.0, 0.34, 0.58, 0.0, 1))
    for cat, tier, rx, ry, lobe, sol, void_f, n_comp in items:
        v_rows.append({"catalogNumber": cat, "scientificName": "Packera aurea", "species_standardized": "Packera aurea",
                       "determiner_tier": tier, "type_status": "", "is_primary_duplicate": 1})
        pts = _contour_xy(rx, ry, lobe=lobe)
        for lid in range(1, (rep_g1 if cat == "G1" else 1) + 1):
            cp = cdir / f"{cat}_leaf{lid}.csv"
            pd.DataFrame({"x": pts[:, 0], "y": pts[:, 1]}).to_csv(cp, index=False)
            m_rows.append({
                "catalogNumber": cat, "plant_individual_id": 0, "leaf_id": lid, "leaf_zone": "basal", "leaf_zone_confidence": 0.95,
                "mask_is_upright": 1, "pixels_per_mm": 10.0, "is_calibrated_mm": 1, "scale_source": "ruler_calibrated", "solidity": sol,
                "archetype_score": 0.95, "lm3_gates_pass": 1, "efa_qc_pass": 1, "is_crop_clipped": 0, "n_lamina_components": n_comp,
                "enclosed_void_frac": void_f, "qc_flag_reason": "PASS", "contour_path": str(cp),
            })
    pd.DataFrame(v_rows).to_csv(base_dir / "vouchers.csv", index=False)
    pd.DataFrame(m_rows).to_csv(base_dir / "manifest.csv", index=False)


class TestStep03TaxonomicClarity(unittest.TestCase):
    """Warm-R contract suite for Step 03 contour orientation, ECT, traits, and plasticity."""

    @classmethod
    def setUpClass(cls) -> None:
        if not rscript_available():
            raise unittest.SkipTest("Rscript not installed.")
        cls._tmp_dir = tempfile.TemporaryDirectory()
        cls.tmp_p = Path(cls._tmp_dir.name)
        cls.db_path = _seed_sqlite_clarity_cohort(cls.tmp_p)
        for tag, rep, v_on, b_on in [("ect_base", 1, True, False), ("ect_rep6", 6, True, False), ("ect_bronze", 1, True, True), ("ect_novoid", 1, False, False)]:
            _seed_ect_cohort(cls.tmp_p / tag, rep_g1=rep, include_void=v_on, include_bronze=b_on)
        (cls.tmp_p / "clarity_s03.yaml").write_text(
            "morphometrics:\n  harmonics: 6\n  num_pcs: 5\n  fourier:\n    min_harmonics: 6\n    max_harmonics: 8\n    min_export_pcs: 5\n    max_export_pcs: 6\n",
            encoding="utf-8",
        )
        cls.json_out, r_script = cls.tmp_p / "clarity_results.json", cls.tmp_p / "run_clarity.R"
        r_script.write_text(
            f"""
            if (grDevices::dev.cur() == 1L) grDevices::pdf(NULL)
            for (.f in c(
              "scripts/morphometrics/step03_config_scale_utils.R", "scripts/morphometrics/step03_ect_traits_utils.R",
              "scripts/morphometrics/dissected_leaf_utils.R", "scripts/morphometrics/fourier_parallel_utils.R",
              "scripts/morphometrics/plasticity_mitigation_utils.R", "scripts/morphometrics/authority_weighting_utils.R",
              "scripts/morphometrics/pub_figures_utils.R", "scripts/core/db_utils.R", "scripts/morphometrics/03_fourier_extractor.R"
            )) source(.f)

            # 1. Vertex-1 CW preservation & horizontal X-axis orientation
            s_area <- function(m) 0.5 * sum(m[, 1L] * c(m[-1L, 2L], m[1L, 2L]) - c(m[-1L, 1L], m[1L, 1L]) * m[, 2L])
            t1 <- seq(0, 2 * pi, length.out = 121L)[-121L]; lx <- -60 * cos(t1); ly <- 24 * sin(t1) * (1 + 0.15 * cos(t1)); phi <- 0.65
            cw_obl <- cbind(250 + lx * cos(phi) - ly * sin(phi), 250 + lx * sin(phi) + ly * cos(phi))
            ccw_obl <- cw_obl[c(1L, seq(120L, 2L)), , drop = FALSE]; cw_rev <- cbind(250 - lx, 250 - ly)
            out_ccw <- ensure_clockwise_keep_start_r(ccw_obl); out_cw <- ensure_clockwise_keep_start_r(cw_obl)
            al_obl <- align_contour_horizontal_xax(cw_obl, is_upright = FALSE); al_rev <- align_contour_horizontal_xax(cw_rev, is_upright = FALSE)

            # 2. Cordate-reniform outline reconstruction & k=1 harmonic power
            h_cord <- compute_efourier_coords_pure(cbind(-0.40 * cos(t1) + 0.18 * cos(2 * t1), 0.72 * sin(t1) + 0.18 * sin(2 * t1)), nb_h = 8L, norm = FALSE, start = FALSE)
            names(h_cord) <- as.vector(outer(c("A", "B", "C", "D"), seq_len(8L), paste0))
            k1_pow_cord <- 0.5 * sum(h_cord[c("A1", "B1", "C1", "D1")]^2)
            rec_cord <- reconstruct_efa_outline(h_cord); tip_idx <- as.integer(round((nrow(rec_cord) + 1L) / 2L))

            # 3. Sheet individual consolidation via bbox_x1..y2 and all-NA refusal
            df_cons <- data.frame(
              catalogNumber = c(rep("SH_BBOX", 4L), rep("SH_NOCOORD", 3L)), plant_individual_id = c(0L, 0L, 1L, 2L, 0L, 0L, 1L),
              leaf_id = c(1L, 2L, 3L, 4L, 1L, 2L, 3L), leaf_zone = rep("basal", 7L),
              dist_from_crown_norm = c(0.08, 0.09, 0.10, 0.09, 0.08, 0.09, 0.10), cy_rel_plant_bbox = c(0.80, 0.82, 0.79, 0.80, 0.80, 0.82, 0.80),
              pixels_per_mm = rep(10.0, 7L), is_calibrated_mm = rep(TRUE, 7L), scale_source = rep("ruler_calibrated", 7L),
              cx_px = rep(NA_real_, 7L), cy_px = rep(NA_real_, 7L), cx_rel_sheet = rep(NA_real_, 7L), cy_rel_sheet = rep(NA_real_, 7L),
              bbox_x1 = c(480, 500, 680, 1100, NA_real_, NA_real_, NA_real_), bbox_y1 = c(480, 490, 490, 490, NA_real_, NA_real_, NA_real_),
              bbox_x2 = c(520, 540, 720, 1140, NA_real_, NA_real_, NA_real_), bbox_y2 = c(520, 530, 530, 530, NA_real_, NA_real_, NA_real_), stringsAsFactors = FALSE
            )
            cons_out <- consolidate_fragmented_sheet_individuals(df_cons)

            # 4 & 5. SQLite Step 03 run (with figures) + CSV ECT invariance runs
            cfg_y <- "{cls.tmp_p / 'clarity_s03.yaml'}"; sql_dir <- "{cls.tmp_p / 'sql_run'}"
            Sys.setenv(PLANTMATICS_DB_PATH = "{cls.db_path}", RUN_ID = "{RUN_ID}")
            Sys.unsetenv("PLANTMATICS_SKIP_FIGURES"); Sys.unsetenv("PACKERA_SKIP_FIGURES")
            run_fourier_extraction(parse_args_robust(c(
              "--config", cfg_y, "--contours-dir", file.path(sql_dir, "contours"), "--vouchers", file.path(sql_dir, "tables", "curated_vouchers.csv"),
              "--output", file.path(sql_dir, "tables", "leaf_efa_harmonics.csv"), "--output-individual", file.path(sql_dir, "tables", "leaf_efa_harmonics_individual.csv"),
              "--figures-dir", file.path(sql_dir, "figures")
            )))
            Sys.unsetenv("PLANTMATICS_DB_PATH"); Sys.unsetenv("RUN_ID"); Sys.setenv(PLANTMATICS_SKIP_FIGURES = "1")
            for (tag in c("ect_base", "ect_rep6", "ect_bronze", "ect_novoid")) {{
              pdir <- file.path("{cls.tmp_p}", tag)
              run_fourier_extraction(parse_args_robust(c(
                "--config", cfg_y, "--contours-dir", file.path(pdir, "contours"), "--manifest", file.path(pdir, "manifest.csv"),
                "--vouchers", file.path(pdir, "vouchers.csv"), "--output", file.path(pdir, "spec.csv"), "--output-individual", file.path(pdir, "ind.csv")
              )))
            }}

            # 6. Ontogenetic / phenophase plasticity residualization guards
            set.seed(101)
            df_nf <- data.frame(
              species_standardized = rep(c("Packera aurea", "Packera dubia"), each = 12L), authority_weight = rep(1.0, 24L),
              spatial_balance_weight = rep(1.0, 24L), morph_tier_source = rep("Basal_Rosette", 24L), is_calibrated_mm = rep(TRUE, 24L),
              total_capitula = rep(10L, 24L), n_fruiting = rep(0L, 24L), n_senescent = rep(0L, 24L),
              log_centsize = rep(seq(1.2, 2.3, length.out = 12L), 2L), stringsAsFactors = FALSE
            )
            df_nf$solidity <- 0.65 + 0.14 * (df_nf$log_centsize - 1.75) + rnorm(24L, 0, 0.008)
            out_nf <- residualize_ontogenetic_phenophase_plasticity(df_nf, target_cols = c("solidity"))
            fr_all <- rep(c(rep(0.0, 8L), rep(c(0.4, 0.6, 0.8, 1.0), 2L)), 2L); is_a <- rep(c(TRUE, FALSE), each = 16L)
            df_fr <- data.frame(
              species_standardized = ifelse(is_a, "Packera aurea", "Packera dubia"), authority_weight = rep(1.0, 32L),
              spatial_balance_weight = rep(1.0, 32L), morph_tier_source = rep("Basal_Rosette", 32L), is_calibrated_mm = rep(TRUE, 32L),
              total_capitula = rep(10L, 32L), n_fruiting = as.integer(round(fr_all * 10L)), n_senescent = rep(0L, 32L),
              log_centsize = rep(seq(1.2, 2.2, length.out = 16L), 2L),
              blade_cie_b_delta = ifelse(is_a, 5.0, 8.0) + 2.0 * fr_all + rnorm(32L, 0, 0.02),
              dissection_index = 1.5 + ifelse(is_a, 2.5, 0.2) * fr_all + rnorm(32L, 0, 0.02),
              AR_inv = ifelse(is_a, 1.2, 1.8) + 1.5 * fr_all + rnorm(32L, 0, 0.02),
              radiate_head_fraction = ifelse(is_a, 0.2, 0.3) + 0.6 * fr_all + rnorm(32L, 0, 0.015), stringsAsFactors = FALSE
            )
            out_fr <- residualize_ontogenetic_phenophase_plasticity(df_fr, target_cols = c("blade_cie_b_delta", "dissection_index", "AR_inv", "radiate_head_fraction"))

            # 7. Herbarium EB batch centering: centroid-pin exclusion & Cramer's V >= 0.45 guard
            set.seed(202)
            inst_c <- c(rep("NCU", 12L), rep("LSU", 12L), rep("NCU", 20L), rep("LSU", 20L))
            df_cent <- data.frame(
              catalogNumber = paste0("C_", seq_len(64L)), species_standardized = rep(c("Packera aurea", "Packera dubia"), times = 32L),
              institutionCode = inst_c, regional_group = c(rep(c("Piedmont", "Coastal_Plain"), each = 6L), rep(c("Piedmont", "Coastal_Plain"), each = 6L), rep("Piedmont", 20L), rep("Coastal_Plain", 20L)),
              is_centroid_or_inst_coord = c(rep(FALSE, 24L), rep(TRUE, 40L)), authority_weight = rep(1.0, 64L), is_passive = rep(FALSE, 64L),
              PC1 = rnorm(64L, 0.0, 0.1), blade_cie_b_delta = ifelse(inst_c == "NCU", 5.0, 14.0) + rnorm(64L, 0, 0.2), stringsAsFactors = FALSE
            )
            eb_cent <- apply_eb_institution_batch_centering(df_cent, target_cols = c("PC1", "blade_cie_b_delta"), min_inst_n = 5L)
            inst_m <- rep(c("NCU", "LSU", "BRIT"), each = 20L)
            reg_m <- c(rep("R1", 14L), rep("R2", 3L), rep("R3", 3L), rep("R1", 3L), rep("R2", 14L), rep("R3", 3L), rep("R1", 3L), rep("R2", 3L), rep("R3", 14L))
            df_mod <- data.frame(
              catalogNumber = paste0("M_", seq_len(60L)), species_standardized = rep(c("Packera aurea", "Packera dubia"), times = 30L),
              institutionCode = inst_m, regional_group = reg_m, is_centroid_or_inst_coord = rep(FALSE, 60L), authority_weight = rep(1.0, 60L),
              is_passive = rep(FALSE, 60L), blade_cie_b_delta = ifelse(inst_m == "NCU", 4.0, ifelse(inst_m == "LSU", 9.0, 14.0)) + rnorm(60L, 0, 0.2), stringsAsFactors = FALSE
            )
            eb_mod <- apply_eb_institution_batch_centering(df_mod, target_cols = c("blade_cie_b_delta"), min_inst_n = 5L)

            res <- list(
              in_ccw_area = s_area(ccw_obl), out_ccw_area = s_area(out_ccw), ccw_p1_diff = max(abs(out_ccw[1L, ] - ccw_obl[1L, ])),
              in_cw_area = s_area(cw_obl), out_cw_area = s_area(out_cw), cw_p1_diff = max(abs(out_cw[1L, ] - cw_obl[1L, ])),
              obl_p1_dx = al_obl[1L, 1L] - mean(al_obl[, 1L]), obl_p61_dx = al_obl[61L, 1L] - mean(al_obl[, 1L]),
              obl_p1_dy = abs(al_obl[1L, 2L] - mean(al_obl[, 2L])), obl_p61_dy = abs(al_obl[61L, 2L] - mean(al_obl[, 2L])),
              rev_p1_dx = al_rev[1L, 1L] - mean(al_rev[, 1L]), rev_p61_dx = al_rev[61L, 1L] - mean(al_rev[, 1L]),
              cord_y_base = as.numeric(rec_cord[1L, 2L]), cord_y_tip = as.numeric(rec_cord[tip_idx, 2L]),
              cord_span_x = diff(range(rec_cord[, 1L], na.rm = TRUE)), cord_span_y = diff(range(rec_cord[, 2L], na.rm = TRUE)),
              k1_pow_cord = k1_pow_cord, cons_pids = as.integer(cons_out$plant_individual_id),
              nf_sol_diff = max(abs(out_nf$solidity - df_nf$solidity)), fr_het_diff = max(abs(out_fr$dissection_index - df_fr$dissection_index)),
              fr_ar_diff = max(abs(out_fr$AR_inv - df_fr$AR_inv)), fr_rad_diff = max(abs(out_fr$radiate_head_fraction - df_fr$radiate_head_fraction)),
              fr_hom_anthesis_diff = max(abs(out_fr$blade_cie_b_delta[fr_all == 0] - df_fr$blade_cie_b_delta[fr_all == 0])),
              fr_hom_fruit_shift = max(abs(out_fr$blade_cie_b_delta[fr_all > 0] - df_fr$blade_cie_b_delta[fr_all > 0])),
              fr_hom_gap_a = abs(mean(out_fr$blade_cie_b_delta[is_a & fr_all > 0]) - mean(out_fr$blade_cie_b_delta[is_a & fr_all == 0])),
              fr_hom_gap_b = abs(mean(out_fr$blade_cie_b_delta[!is_a & fr_all > 0]) - mean(out_fr$blade_cie_b_delta[!is_a & fr_all == 0])),
              eb_cent_phot_diff = max(abs(eb_cent$blade_cie_b_delta - df_cent$blade_cie_b_delta)),
              eb_cent_pc1_diff = max(abs(eb_cent$PC1 - df_cent$PC1)), eb_mod_phot_diff = max(abs(eb_mod$blade_cie_b_delta - df_mod$blade_cie_b_delta))
            )
            writeLines(jsonlite::toJSON(res, auto_unbox = TRUE, na = "null"), "{cls.json_out}")
            """,
            encoding="utf-8",
        )
        proc = run_warm_r_cmd(["Rscript", str(r_script)], stage_key="step03", cwd=PROJECT_ROOT)
        if proc.returncode != 0:
            raise RuntimeError(f"Step 03 clarity Rscript failed:\n{proc.stderr}\n{proc.stdout}")
        cls.data = json.loads(cls.json_out.read_text(encoding="utf-8"))

    @classmethod
    def tearDownClass(cls) -> None:
        if hasattr(cls, "_tmp_dir"):
            cls._tmp_dir.cleanup()

    def test_01_vertex1_cw_preservation_and_horizontal_xax_alignment(self) -> None:
        d = self.data
        self.assertGreater(float(d["in_ccw_area"]), 0.0)
        for k in ("out_ccw_area", "in_cw_area", "out_cw_area"):
            self.assertLess(float(d[k]), 0.0, msg=k)
        for k in ("ccw_p1_diff", "cw_p1_diff"):
            self.assertLess(float(d[k]), 1e-12, msg=k)
        for k in ("obl_p1_dx", "rev_p1_dx"):
            self.assertLess(float(d[k]), -40.0, msg=k)
        for k in ("obl_p61_dx", "rev_p61_dx"):
            self.assertGreater(float(d[k]), 40.0, msg=k)
        self.assertLess(max(float(d["obl_p1_dy"]), float(d["obl_p61_dy"])), 1e-2)

    def test_02_cordate_outline_reconstruction_and_k1_harmonic_power(self) -> None:
        d = self.data
        self.assertLess(float(d["cord_y_base"]), 0.0)
        self.assertGreater(float(d["cord_y_tip"]), 0.0)
        self.assertLess(float(d["cord_y_base"]), float(d["cord_y_tip"]))
        self.assertGreater(float(d["cord_span_x"]), float(d["cord_span_y"]))
        self.assertGreater(float(d["k1_pow_cord"]), 0.0)
        fig_src = (PROJECT_ROOT / "scripts" / "figures" / "step03_figures.R").read_text(encoding="utf-8")
        self.assertIsNone(re.search(r"\|\|[ \t]*k[ \t]*==[ \t]*1L?", fig_src))
        fig_dir = self.tmp_p / "sql_run" / "figures"
        fig03a = fig_dir / "fig03a_efa_morphospace_and_harmonics.png"
        self.assertTrue(fig03a.exists() and fig03a.stat().st_size > 1000)
        self.assertEqual(list(fig_dir.glob("*.pdf")), [])
        self.assertFalse((PROJECT_ROOT / "Rplots.pdf").exists())
        spec_df = pd.read_csv(self.tmp_p / "sql_run" / "tables" / "leaf_efa_harmonics.csv")
        self.assertGreater(float(0.5 * (spec_df["A1"].dropna() ** 2 + spec_df["D1"].dropna() ** 2).mean()), 0.0)

    def test_03_consolidate_fragmented_sheet_individuals_bbox_and_all_na(self) -> None:
        self.assertEqual(self.data["cons_pids"], [0, 0, 0, 2, 0, 0, 1])

    def test_04_ect_otu_invariance_multi_island_and_bronze_zero_weight(self) -> None:
        dfs = {t: pd.read_csv(self.tmp_p / t / "spec.csv").sort_values("catalogNumber") for t in ("ect_base", "ect_rep6", "ect_bronze", "ect_novoid")}
        s_base, s_rep6, s_brz, s_nov = dfs["ect_base"], dfs["ect_rep6"], dfs["ect_bronze"], dfs["ect_novoid"]
        ect_cols = ["ECT_PC1", "ECT_PC2", "ECT_PC3"]
        for col in ect_cols:
            self.assertLess(_col_diff(s_base[col].to_numpy(dtype=float), s_rep6[col].to_numpy(dtype=float)), 1e-6)
        v_row = s_base[s_base["catalogNumber"] == "G_VOID"].iloc[0]
        self.assertTrue(all(np.isfinite(float(v_row[c])) for c in ect_cols) and pd.isna(v_row["PC1"]))
        for hcol in [c for c in s_base.columns if re.match(r"^[ABCD][0-9]+$", c)]:
            self.assertTrue(pd.isna(v_row[hcol]), f"Expected NA {hcol} for non-Jordan G_VOID")
        g_ids = [f"G{i}" for i in range(1, 7)]
        b_gold, n_gold = s_base[s_base["catalogNumber"].isin(g_ids)], s_nov[s_nov["catalogNumber"].isin(g_ids)]
        self.assertGreater(max(_col_diff(b_gold[c].to_numpy(dtype=float), n_gold[c].to_numpy(dtype=float)) for c in ect_cols), 1e-4)
        brz_gold = s_brz[s_brz["catalogNumber"] != "B_EXTREME"]
        self.assertEqual(float(s_brz[s_brz["catalogNumber"] == "B_EXTREME"].iloc[0]["authority_weight"]), 0.0)
        for col in ect_cols:
            self.assertLess(_col_diff(s_base[col].to_numpy(dtype=float), brz_gold[col].to_numpy(dtype=float)), 1e-6)

    def test_05_adaxial_abaxial_photometry_step02_traits_and_sqlite_harmonics(self) -> None:
        tables_dir = self.tmp_p / "sql_run" / "tables"
        ind_df, spec_df = pd.read_csv(tables_dir / "leaf_efa_harmonics_individual.csv"), pd.read_csv(tables_dir / "leaf_efa_harmonics.csv")
        self.assertTrue({"B1", "C1"}.issubset(ind_df.columns) and {"B1", "C1"}.isdisjoint(spec_df.columns))
        conn = sqlite3.connect(str(self.db_path))
        conn.row_factory = sqlite3.Row
        uv = int(conn.execute("PRAGMA user_version;").fetchone()[0])
        self.assertEqual(uv, PLANTMATICS_SCHEMA_USER_VERSION)
        self.assertEqual(uv, 20261127)
        self.assertTrue({"B1", "C1"}.isdisjoint({r[1] for r in conn.execute("PRAGMA table_info(otu_morphometrics);").fetchall()}))
        lh_abax = dict(conn.execute("SELECT * FROM leaf_efa_harmonics WHERE run_id=? AND catalogNumber='V_ABAX_1' AND leaf_id=1;", (RUN_ID,)).fetchone())
        lh_void = dict(conn.execute("SELECT * FROM leaf_efa_harmonics WHERE run_id=? AND catalogNumber='V_VOID_SQL' AND leaf_id=1;", (RUN_ID,)).fetchone())
        om_row = dict(conn.execute("SELECT * FROM otu_morphometrics WHERE run_id=? AND catalogNumber='V_ABAX_1';", (RUN_ID,)).fetchone())
        vm_row = dict(conn.execute("SELECT * FROM v_integrative_triage_master WHERE run_id=? AND catalogNumber='V_ABAX_1';", (RUN_ID,)).fetchone())
        conn.close()
        self.assertTrue(lh_abax["B1"] is not None and lh_abax["C1"] is not None and np.isfinite(float(lh_abax["B1"])) and np.isfinite(float(lh_abax["C1"])))
        for hcol in ("A1", "B1", "C1", "D1"):
            self.assertIsNone(lh_void[hcol], f"Expected SQL NULL for {hcol} on non-Jordan leaf")
        expected = [("blade_cie_l_sheet_delta", -19.0), ("blade_cie_l_q85_sheet_delta", -14.0), ("abaxial_leaf_fraction", 0.5),
                    ("abaxial_cie_l_sheet_delta", -3.5), ("basal_cordate_index", 0.30), ("base_angle_deg", 114.0)]
        for label, row in (("otu_morphometrics", om_row), ("v_integrative_triage_master", vm_row)):
            for col, val in expected:
                self.assertAlmostEqual(float(row[col]), val, places=4, msg=f"{label}.{col}")
            self.assertEqual(int(row["is_abaxial_pressed_leaf"]), 1, msg=label)
            self.assertTrue(42.0 <= float(row["lamina_trace_length_mm"]) <= 46.0, msg=label)

    def test_06_ontogenetic_phenophase_plasticity_residualization_guards(self) -> None:
        d = self.data
        for k in ("nf_sol_diff", "fr_het_diff", "fr_ar_diff", "fr_rad_diff", "fr_hom_anthesis_diff"):
            self.assertLess(float(d[k]), 1e-12, msg=k)
        self.assertGreater(float(d["fr_hom_fruit_shift"]), 1.0)
        self.assertLess(max(float(d["fr_hom_gap_a"]), float(d["fr_hom_gap_b"])), 0.08)

    def test_07_eb_institution_batch_centering_centroid_and_cramers_v_guards(self) -> None:
        d = self.data
        self.assertGreater(float(d["eb_cent_phot_diff"]), 1.0)
        self.assertLess(float(d["eb_cent_pc1_diff"]), 1e-12)
        self.assertLess(float(d["eb_mod_phot_diff"]), 1e-12)


if __name__ == "__main__":
    unittest.main()
