#!/usr/bin/env python3
"""Warm-R tests for Step 07 Hopkins shrinkage, active-anthesis gating, edaphic z-scores, OOF precision, and 4-stream simplex."""

from __future__ import annotations

import json
import os
import shutil
import sys
import tempfile
import unittest
from pathlib import Path
from typing import Any, Dict

PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from scripts.tests._r_session import run_warm_r_cmd

SKIP_R_COND = (
    os.environ.get("PLANTMATICS_FAST_AUDIT") == "1"
    or os.environ.get("PACKERA_FAST_AUDIT") == "1"
    or shutil.which("Rscript") is None
)

_WARM_R_CODE = r"""
if (grDevices::dev.cur() == 1L) grDevices::pdf(NULL)
args <- commandArgs(trailingOnly = TRUE)
out_json <- args[1L]
suppressPackageStartupMessages({ library(jsonlite) })
source("scripts/analysis/07_triage_dashboard_synthesis.R")

# 1. Anthesis-dominance + spatial_balance weighted Empirical Bayes Hopkins shrinkage & residual coalescing
mk_hop_df <- function(low_wt_doy, full_wt_doy) {
  lats <- c(seq(32.0, 38.0, length.out = 10L), 40.0, 40.0)
  data.frame(
    catalogNumber = sprintf("HOP_%02d", seq_len(12L)),
    species_standardized = rep("Packera dubia", 12L),
    determiner_tier = rep("Tier_1_Gold", 12L),
    authority_weight = rep(2.0, 12L),
    spatial_balance_weight = c(rep(1.0, 10L), 0.25, 1.0),
    anthesis_dominance_weight = c(rep(1.0, 10L), 0.20, 1.0),
    decimalLatitude = lats,
    elevation_m = rep(150.0, 12L),
    bio1_temp_mean = rep(14.5, 12L),
    phenological_state = rep("anthesis", 12L),
    doy = c(100.0 + 4.0 * (lats[1:10] - 35.0), low_wt_doy, full_wt_doy),
    stringsAsFactors = FALSE
  )
}
h_base <- compute_hopkins_adjusted_phenology(mk_hop_df(120.0, 120.0))
h_low <- compute_hopkins_adjusted_phenology(mk_hop_df(210.0, 120.0))
h_full <- compute_hopkins_adjusted_phenology(mk_hop_df(120.0, 210.0))
shift_low_wt <- max(abs(h_base$hopkins_adjusted_doy[1:10] - h_low$hopkins_adjusted_doy[1:10]))
shift_full_wt <- max(abs(h_base$hopkins_adjusted_doy[1:10] - h_full$hopkins_adjusted_doy[1:10]))

# Residual coalescing precedence: delta_doy_resid > hopkins_adjusted_doy_resid > delta_doy
df_coal <- data.frame(
  catalogNumber = c("C_D_RES", "C_HOP_RES", "C_DDOY"),
  species_standardized = rep("Packera dubia", 3L),
  phenological_state = rep("anthesis", 3L),
  doy = c(115, 118, 120),
  delta_doy_resid = c(7.5, NA_real_, NA_real_),
  hopkins_adjusted_doy_resid = c(99.0, -11.2, NA_real_),
  delta_doy = c(88.0, 77.0, 5.4),
  stringsAsFactors = FALSE
)
h_coal <- compute_hopkins_adjusted_phenology(df_coal)

# 2. Active-anthesis gating: Late_Fruiting_Senescent_Excluded, anthesis_dominance_weight < 0.20, unknown/fruit/bud
df_gate <- data.frame(
  catalogNumber = c("G_PASS", "G_LATE_EXCL", "G_LOW_WT", "G_FRUIT", "G_UNK_FALSE", "G_UNK_VER", "G_SENESC"),
  species_standardized = rep("Packera dubia", 7L),
  phenological_state = c("anthesis", "anthesis", "anthesis", "fruit", "unknown", "unknown", "senescent"),
  anthesis_filter_status = c("Active_Anthesis", "Late_Fruiting_Senescent_Excluded", "Active_Anthesis", "Active_Anthesis", "Active_Anthesis", "Active_Anthesis", "Active_Anthesis"),
  anthesis_dominance_weight = c(0.85, 0.90, 0.12, 0.85, 0.85, 0.85, 0.85),
  is_verified_anthesis = c(TRUE, TRUE, TRUE, FALSE, FALSE, TRUE, FALSE),
  doy = c(115, 115, 115, 140, 115, 115, 150),
  stringsAsFactors = FALSE
)
anth_mask <- resolve_active_anthesis_mask(df_gate)
h_gate <- compute_hopkins_adjusted_phenology(df_gate)

# 3. Edaphic predictor hygiene & dynamic edaphic_z_score (banned cartographic/buffer excluded, smooth indices included)
df_ed <- data.frame(
  catalogNumber = c(paste0("ED_CORE_", 1:6), "ED_OUTCROP", "ED_NOISE"),
  species_standardized = rep("Packera dubia", 8L),
  determiner_tier = c(rep("Tier_1_Gold", 7L), "Tier_3_Bronze"),
  is_confident_label_noise = c(rep(FALSE, 7L), TRUE),
  primary_env_predictors = rep(
    "lithic_outcrop_index;wetland_seepage_index;summer_moisture_deficit_idx;maritime_continentality_idx;ssurgo_rock_outcrop_pct;ssurgo_aws_0_50cm;ssurgo_lithic_series_pct;ultramafic_calcareous_gradient;min_bedrock_depth_buffer_cm",
    8L
  ),
  lithic_outcrop_index = c(0.10, 0.12, 0.08, 0.11, 0.09, 0.13, 0.95, 0.99),
  wetland_seepage_index = c(0.15, 0.18, 0.14, 0.16, 0.17, 0.13, 0.92, 0.99),
  summer_moisture_deficit_idx = c(0.20, 0.22, 0.19, 0.21, 0.23, 0.18, 0.95, 0.99),
  maritime_continentality_idx = c(0.40, 0.42, 0.39, 0.41, 0.43, 0.38, 0.96, 0.99),
  ssurgo_rock_outcrop_pct = c(0, 0, 5, 0, 10, 0, 90, 95),
  ssurgo_aws_0_50cm = c(8, 9, 10, 8, 9, 10, 1, 1),
  ssurgo_lithic_series_pct = c(5, 10, 5, 10, 5, 10, 85, 90),
  ultramafic_calcareous_gradient = c(-0.2, -0.1, 0.0, -0.2, 0.1, 0.0, 0.9, 0.9),
  min_bedrock_depth_buffer_cm = c(150, 140, 160, 155, 145, 150, 10, 5),
  fine_soil_ph = c(5.5, 5.6, 5.4, 5.5, 5.6, 5.5, 7.6, 8.8),
  stringsAsFactors = FALSE
)
ed_preds_meta <- .resolve_step07_edaphic_predictors(df_ed)
df_ed_fb <- df_ed; df_ed_fb$primary_env_predictors <- NULL
ed_preds_fb <- .resolve_step07_edaphic_predictors(df_ed_fb)
ez_vals <- .compute_dynamic_edaphic_z_scores(
  df_ed, df_ed$species_standardized, "Packera dubia",
  list("Packera dubia" = 1:6), clean_coh_mask = c(rep(TRUE, 6L), FALSE, FALSE)
)

# 4. Taxon-dispersion tau_ecotype / tau_cryptic scaling & noise isolation
df_tau <- data.frame(
  catalogNumber = sprintf("T_%02d", 1:18),
  species_standardized = c(rep("Packera tight", 8L), rep("Packera broad", 8L), "Packera tight", "Packera tight"),
  determiner_tier = c(rep("Tier_1_Gold", 16L), "Tier_3_Bronze", "Tier_3_Bronze"),
  is_confident_label_noise = c(rep(FALSE, 16L), TRUE, TRUE),
  is_active_anthesis = rep(TRUE, 18L),
  delta_doy_resid = c(-3, -2, -1, 0, 0, 1, 2, 3, -18, -14, -9, -3, 3, 9, 14, 18, 95, -95),
  stringsAsFactors = FALSE
)
taus_1 <- .compute_taxon_pheno_thresholds(df_tau)
df_tau_clean <- df_tau[1:16, , drop = FALSE]
taus_2 <- .compute_taxon_pheno_thresholds(df_tau_clean)

# 5. OOF precision map, Stage 1 cross-modal non-contradiction gate, and 4-stream simplex renormalization
prec_map <- .compute_oof_class_precision_map(
  y_true_vec = c(rep("Packera A", 8L), rep("Packera B", 8L)),
  y_pred_vec = c(rep("Packera A", 8L), rep("Packera A", 5L), rep("Packera B", 3L))
)

df_sx <- data.frame(
  catalogNumber = c(paste0("PA_", 1:5), paste0("PB_", 1:5), "SX_ALL4", "SX_NO_PHENO_GENO", "SX_MORPH_ONLY", "SX_G2_NO_EDAPH", "GOLD_NONCONTRA"),
  species_standardized = c(rep("Packera aurea", 5L), rep("Packera obovata", 5L), rep("Packera aurea", 5L)),
  determiner_tier = c(rep("Tier_1_Gold", 10L), rep("Tier_2_Silver", 4L), "Tier_1_Gold"),
  decimalLatitude = c(rep(35.8, 5L), rep(35.9, 5L), 35.8, 38.5, NA_real_, NA_real_, 35.8),
  decimalLongitude = c(rep(-79.0, 5L), rep(-79.1, 5L), -79.0, -81.5, NA_real_, NA_real_, -79.0),
  stateProvince = c(rep("NC", 12L), NA_character_, NA_character_, "NC"),
  county = c(rep("Orange", 12L), NA_character_, NA_character_, "Orange"),
  usda_county_status = rep("in_range", 15L),
  phenological_state = c(rep("anthesis", 11L), "fruit", "fruit", "anthesis", "anthesis"),
  doy = c(rep(115, 5L), rep(140, 5L), 116, 150, 150, 117, 115),
  hopkins_adjusted_doy = c(rep(115, 5L), rep(140, 5L), 116, NA_real_, NA_real_, 117, 115),
  delta_doy_resid = c(rep(0.5, 10L), 1.0, NA_real_, NA_real_, 2.0, 0.5),
  PC1 = c(rep(-2.0, 5L), rep(2.0, 5L), -2.0, -1.9, -2.0, -1.95, 1.4),
  PC2 = rep(0.0, 15L),
  can1 = c(rep(-2.0, 5L), rep(2.0, 5L), -2.0, -1.9, -2.0, -1.95, 1.4),
  can2 = rep(0.0, 15L),
  cda_predicted_taxon = c(rep("Packera aurea", 5L), rep("Packera obovata", 5L), rep("Packera aurea", 4L), "Packera obovata"),
  cda_posterior_prob = c(rep(0.95, 10L), 0.90, 0.90, 0.90, 0.90, 0.80),
  self_morph_mahal_d2 = c(rep(1.0, 14L), 15.0),
  edaphic_best_fit_taxon = c(rep("Packera aurea", 5L), rep("Packera obovata", 5L), "Packera aurea", "Packera aurea", NA_character_, NA_character_, "Packera aurea"),
  rf_predicted_taxon = c(rep("Packera aurea", 5L), rep("Packera obovata", 5L), "Packera aurea", "Packera aurea", NA_character_, NA_character_, "Packera aurea"),
  rf_predicted_prob = c(rep(0.92, 10L), 0.88, 0.88, NA_real_, NA_real_, 0.92),
  edaphic_niche_overlap_d = c(rep(0.85, 10L), 0.82, 0.82, NA_real_, NA_real_, 0.88),
  edaphic_z_score = c(rep(0.3, 10L), 0.4, 0.4, NA_real_, NA_real_, 0.4),
  geno_linkage_tier = c("G1_Direct_Voucher", "G1_Direct_Voucher", rep(NA_character_, 3L), "G1_Direct_Voucher", "G1_Direct_Voucher", rep(NA_character_, 3L), "G1_Direct_Voucher", NA_character_, NA_character_, "G2_Conditional_Deme_Hull", NA_character_),
  geno_quality_pass = c(TRUE, TRUE, rep(NA, 3L), TRUE, TRUE, rep(NA, 3L), TRUE, NA, NA, TRUE, NA),
  geno_predicted_taxon = c("Packera aurea", "Packera aurea", rep(NA_character_, 3L), "Packera obovata", "Packera obovata", rep(NA_character_, 3L), "Packera aurea", NA_character_, NA_character_, "Packera aurea", NA_character_),
  geno_posterior_prob = c(0.96, 0.95, rep(NA_real_, 3L), 0.96, 0.95, rep(NA_real_, 3L), 0.95, NA_real_, NA_real_, 0.85, NA_real_),
  geno_deme_weight_scale = c(rep(NA_real_, 13L), 0.66, NA_real_),
  morph_tier_source = rep("Basal_Rosette", 15L),
  basal_leaf_count = rep(3L, 15L),
  stringsAsFactors = FALSE
)
dec_sx <- apply_taxonomic_decision_matrix(df_sx, intraspecific_mode = FALSE, opts = list(demo = TRUE))
get_sx <- function(cid) {
  k <- which(dec_sx$catalogNumber == cid)[1L]
  c(m = dec_sx$w_morph_simplex[k], e = dec_sx$w_edaphic_simplex[k], p = dec_sx$w_pheno_simplex[k], g = dec_sx$w_geno_simplex[k])
}

write_json(list(
  shift_low_wt = shift_low_wt, shift_full_wt = shift_full_wt,
  coal_ddoy = as.numeric(h_coal$delta_doy_resid),
  anth_mask = as.logical(anth_mask),
  gate_ddoy_na = is.na(h_gate$delta_doy_resid),
  ed_preds_meta = ed_preds_meta,
  ed_preds_fb = ed_preds_fb,
  ez_core = as.numeric(ez_vals[1L]),
  ez_outcrop = as.numeric(ez_vals[7L]),
  tau_eco_tight = as.numeric(taus_1$tau_ecotype[1L]),
  tau_eco_broad = as.numeric(taus_1$tau_ecotype[9L]),
  tau_cry_tight = as.numeric(taus_1$tau_cryptic[1L]),
  tau_cry_broad = as.numeric(taus_1$tau_cryptic[9L]),
  tau_noise_invariant = isTRUE(all.equal(taus_1$tau_ecotype[1:16], taus_2$tau_ecotype[1:16], tolerance = 1e-8)),
  prec_a = as.numeric(prec_map[["Packera A"]]),
  prec_b = as.numeric(prec_map[["Packera B"]]),
  sx_all4 = as.list(get_sx("SX_ALL4")),
  sx_no_pg = as.list(get_sx("SX_NO_PHENO_GENO")),
  sx_morph_only = as.list(get_sx("SX_MORPH_ONLY")),
  sx_g2_no_ed = as.list(get_sx("SX_G2_NO_EDAPH")),
  gold_noncontra_sub = as.character(dec_sx$taxonomic_treatment_subtype[dec_sx$catalogNumber == "GOLD_NONCONTRA"][1L])
), out_json, auto_unbox = TRUE, na = "null", digits = NA)
"""


@unittest.skipIf(SKIP_R_COND, "Skipping warm-R Hopkins/consensus suite under fast audit or missing Rscript.")
class TestStep07HopkinsAndConsensus(unittest.TestCase):
    """Contract tests for Hopkins weighting, anthesis gating, smooth edaphic z-scores, and 4-stream simplex."""

    res: Dict[str, Any] = {}

    @classmethod
    def setUpClass(cls) -> None:
        with tempfile.TemporaryDirectory(prefix="step07_hop_cons_") as td:
            tmp = Path(td)
            r_p, out_j = tmp / "probe.R", tmp / "out.json"
            r_p.write_text(_WARM_R_CODE, encoding="utf-8")
            proc = run_warm_r_cmd(["Rscript", str(r_p), str(out_j)], stage_key="step07", cwd=PROJECT_ROOT, skip_figures=True)
            if proc.returncode != 0:
                raise AssertionError(f"Warm R probe failed ({proc.returncode}):\n{proc.stdout}\n{proc.stderr}")
            cls.res = json.loads(out_j.read_text(encoding="utf-8"))

    def test_hopkins_weighting_residual_coalescing_and_anthesis_gating(self) -> None:
        r = self.res
        self.assertLess(float(r["shift_low_wt"]), float(r["shift_full_wt"]))
        self.assertAlmostEqual(float(r["coal_ddoy"][0]), 7.5, places=3)
        self.assertAlmostEqual(float(r["coal_ddoy"][1]), -11.2, places=3)
        self.assertAlmostEqual(float(r["coal_ddoy"][2]), 5.4, places=3)
        self.assertEqual(r["anth_mask"], [True, False, False, False, False, True, False])
        self.assertEqual(r["gate_ddoy_na"], [False, True, True, True, True, False, True])

    def test_smooth_edaphic_predictor_hygiene_and_taxon_scaled_pheno_taus(self) -> None:
        r = self.res
        banned = {
            "ssurgo_rock_outcrop_pct",
            "ssurgo_aws_0_50cm",
            "ssurgo_lithic_series_pct",
            "ultramafic_calcareous_gradient",
            "min_bedrock_depth_buffer_cm",
        }
        smooth = {
            "lithic_outcrop_index",
            "wetland_seepage_index",
            "summer_moisture_deficit_idx",
            "maritime_continentality_idx",
        }
        for key in ("ed_preds_meta", "ed_preds_fb"):
            preds = set(r[key])
            self.assertTrue(banned.isdisjoint(preds), f"Banned edaphic predictors leaked into {key}: {preds & banned}")
            self.assertTrue(smooth.issubset(preds), f"Missing smooth predictors in {key}: {smooth - preds}")
        self.assertLess(float(r["ez_core"]), 1.0)
        self.assertGreater(float(r["ez_outcrop"]), 2.5)
        self.assertLess(float(r["tau_eco_tight"]), float(r["tau_eco_broad"]))
        self.assertLess(float(r["tau_cry_tight"]), float(r["tau_cry_broad"]))
        self.assertGreaterEqual(float(r["tau_eco_tight"]), 10.0)
        self.assertLessEqual(float(r["tau_cry_broad"]), 35.0)
        self.assertTrue(bool(r["tau_noise_invariant"]))

    def test_oof_precision_stage1_non_contradiction_and_four_stream_simplex(self) -> None:
        r = self.res
        self.assertLess(float(r["prec_a"]), float(r["prec_b"]))
        self.assertNotEqual(r["gold_noncontra_sub"], "Misidentification_Reassign")
        for key, zero_cols, pos_cols in (
            ("sx_all4", (), ("m", "e", "p", "g")),
            ("sx_no_pg", ("p", "g"), ("m", "e")),
            ("sx_morph_only", ("e", "p", "g"), ("m",)),
            ("sx_g2_no_ed", ("e",), ("m", "p", "g")),
        ):
            sx = r[key]
            self.assertAlmostEqual(sum(float(sx[k]) for k in ("m", "e", "p", "g")), 1.0, places=5)
            for k in zero_cols:
                self.assertAlmostEqual(float(sx[k]), 0.0, places=6)
            for k in pos_cols:
                self.assertGreater(float(sx[k]), 0.0)


if __name__ == "__main__":
    unittest.main()
