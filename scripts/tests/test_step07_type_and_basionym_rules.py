#!/usr/bin/env python3
"""Warm-R tests for Step 07 nomenclatural type Rules A/B/C, heterotypic synonymy, passive type quarantine, and open concepts."""

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

# 1. Multi-taxon type control Rules A, B, C + heterotypic Type_Synonymy_Candidate + passive anchor quarantine
df_mt <- data.frame(
  catalogNumber = c(
    "TYPE_HOMO_DUB", "TYPE_HET_DIV", "TYPE_UNC_COORD", "OPEN_AFF_CONCEPT", "HYB_FORMULA_CONCEPT",
    paste0("DUB_", 1:5), paste0("OBO_", 1:5), "NON_TYPE_RULE_BC_REASSIGN", "NON_TYPE_PROTECTED", "NON_TYPE_WEAK_RF"
  ),
  species_standardized = c(
    "Packera dubia", "Packera dubia", "Packera dubia", "Packera dubia", "Packera dubia",
    rep("Packera dubia", 5L), rep("Packera obovata", 5L), "Packera dubia", "Packera dubia", "Packera dubia"
  ),
  species_raw = c(
    "Packera dubia", "Senecio smallii", "Packera dubia", "Packera aff. dubia", "Packera dubia x obovata",
    rep("Packera dubia", 5L), rep("Packera obovata", 5L), "Packera dubia", "Packera dubia", "Packera dubia"
  ),
  determiner_tier = c(
    "Tier_0_Type_Anchor", "Tier_0_Type_Anchor", "Tier_0_Type_Anchor", "Tier_1_Gold", "Tier_2_Silver",
    rep("Tier_1_Gold", 10L), "Tier_2_Silver", "Tier_2_Silver", "Tier_2_Silver"
  ),
  is_type_specimen = c(TRUE, TRUE, TRUE, FALSE, FALSE, rep(FALSE, 13L)),
  type_status = c("HOLOTYPE", "HOLOTYPE", "ISOTYPE", "none", "none", rep("none", 13L)),
  typified_basionym = c("Senecio tomentosus", "Senecio smallii", "Senecio tomentosus", NA_character_, NA_character_, rep(NA_character_, 13L)),
  concept_synonymy_status = c("Homotypic_Synonym", "Heterotypic_Synonym", "Homotypic_Synonym", "Open_Nomenclature", "Hybrid_Formula", rep("Homotypic_Synonym", 13L)),
  anchor_tier_source = c("Tier_0_Type_Anchor", "Passive_Heterotypic_Type_Audit", "Passive_Type_Uncertain_Coord", "Passive_Open_or_Hybrid_Concept", "Passive_Open_or_Hybrid_Concept", rep("Tier_1_Gold_Core", 10L), rep("Tier_2_Silver_Regional", 3L)),
  is_open_or_hybrid_concept = c(FALSE, FALSE, FALSE, TRUE, TRUE, rep(FALSE, 13L)),
  is_passive = c(FALSE, TRUE, TRUE, TRUE, TRUE, rep(FALSE, 13L)),
  decimalLatitude = c(35.80, 35.81, 35.82, 35.80, 35.81, rep(35.80, 5L), rep(35.81, 5L), 35.81, 35.80, 35.80),
  decimalLongitude = c(-79.00, -79.01, -79.02, -79.00, -79.01, rep(-79.00, 5L), rep(-79.01, 5L), -79.01, -79.00, -79.00),
  stateProvince = rep("NC", 18L), county = rep("Orange", 18L), usda_county_status = rep("in_range", 18L),
  phenological_state = rep("anthesis", 18L),
  doy = c(115, 120, 115, 115, 117, rep(115, 5L), rep(118, 5L), 118, 115, 115),
  PC1 = c(-2.0, 2.1, -2.0, -1.9, 0.0, rep(-2.0, 5L), rep(2.0, 5L), 2.15, -1.95, 1.65),
  PC2 = c(0.0, 0.05, 0.0, 0.1, 0.02, rep(0.0, 10L), 0.05, 0.0, 0.10),
  can1 = c(-2.0, 2.1, -2.0, -1.9, 0.0, rep(-2.0, 5L), rep(2.0, 5L), 2.15, -1.95, 1.65),
  can2 = c(0.0, 0.05, 0.0, 0.1, 0.02, rep(0.0, 10L), 0.05, 0.0, 0.10),
  cda_predicted_taxon = c("Packera dubia", "Packera obovata", "Packera dubia", "Packera dubia", "Packera obovata", rep("Packera dubia", 5L), rep("Packera obovata", 5L), "Packera obovata", "Packera dubia", "Packera obovata"),
  cda_posterior_prob = c(0.95, 0.93, 0.94, 0.90, 0.51, rep(0.95, 10L), 0.89, 0.92, 0.65),
  cda_second_taxon = c("Packera obovata", "Packera dubia", "Packera obovata", "Packera obovata", "Packera dubia", rep("Packera obovata", 5L), rep("Packera dubia", 5L), "Packera dubia", "Packera obovata", "Packera dubia"),
  cda_second_posterior_prob = c(0.05, 0.07, 0.06, 0.10, 0.49, rep(0.05, 10L), 0.11, 0.08, 0.35),
  cda_posterior_margin = c(0.90, 0.86, 0.88, 0.80, 0.02, rep(0.90, 10L), 0.78, 0.84, 0.30),
  gmm_uncertainty = c(0.02, 0.03, 0.02, 0.05, 0.36, rep(0.02, 10L), 0.03, 0.02, 0.04),
  is_typical_morphotype = c(TRUE, FALSE, TRUE, TRUE, FALSE, rep(TRUE, 10L), FALSE, TRUE, FALSE),
  type_morph_mahal_d2 = c(0.5, 26.0, 0.8, 1.2, 8.0, rep(1.0, 10L), 28.0, 1.4, 22.0),
  self_morph_mahal_d2 = c(0.5, 26.0, 0.8, 1.2, 8.0, rep(1.0, 10L), 28.0, 1.4, 22.0),
  nearest_type_morph_taxon = c("Packera dubia", "Packera obovata", "Packera dubia", "Packera dubia", "Packera obovata", rep("Packera dubia", 5L), rep("Packera obovata", 5L), "Packera obovata", "Packera dubia", "Packera obovata"),
  nearest_type_niche_taxon = c("Packera dubia", "Packera obovata", "Packera dubia", "Packera dubia", "Packera dubia", rep("Packera dubia", 5L), rep("Packera obovata", 5L), "Packera obovata", "Packera obovata", NA_character_),
  edaphic_best_fit_taxon = c("Packera dubia", "Packera obovata", "Packera dubia", "Packera dubia", "Packera dubia", rep("Packera dubia", 5L), rep("Packera obovata", 5L), "Packera obovata", "Packera obovata", "Packera obovata"),
  rf_predicted_taxon = c("Packera dubia", "Packera obovata", "Packera dubia", "Packera dubia", "Packera dubia", rep("Packera dubia", 5L), rep("Packera obovata", 5L), "Packera obovata", "Packera obovata", "Packera obovata"),
  rf_predicted_prob = c(0.93, 0.90, 0.92, 0.90, 0.60, rep(0.93, 10L), 0.88, 0.85, 0.35),
  edaphic_z_score = c(0.3, 3.2, 0.4, 0.5, 1.0, rep(0.3, 10L), 3.4, 2.8, 2.6),
  type_edaphic_z_score = c(0.3, 3.2, 0.4, 0.5, 1.0, rep(0.3, 10L), 3.6, 2.8, 1.5),
  gmm_cluster = rep(1L, 18L),
  morph_tier_source = rep("Basal_Rosette", 18L),
  basal_leaf_count = rep(3L, 18L),
  stringsAsFactors = FALSE
)

anch_step07 <- .resolve_step07_anchor_mask(df_mt)
anch_comp <- select_taxon_anchors(df_mt)$anchor_mask
res_mt <- apply_taxonomic_decision_matrix(df_mt, intraspecific_mode = FALSE, opts = list(demo = TRUE))

# 2. Divergent heterotypic synonym type in --intraspecific-mode Stage 3
df_intra <- data.frame(
  catalogNumber = c("INTRA_HOMO_TYPE", "INTRA_HET_TYPE", sprintf("INTRA_CORE_%02d", 1:6), sprintf("INTRA_K2_%02d", 1:6)),
  species_standardized = rep("Packera anonyma", 14L),
  species_raw = c("Packera anonyma", "Senecio smallii", rep("Packera anonyma", 12L)),
  determiner_tier = c("Tier_0_Type_Anchor", "Tier_0_Type_Anchor", rep("Tier_1_Gold", 6L), rep("Tier_2_Silver", 6L)),
  is_type_specimen = c(TRUE, TRUE, rep(FALSE, 12L)),
  type_status = c("HOLOTYPE", "HOLOTYPE", rep("none", 12L)),
  typified_basionym = c("Senecio anonymus", "Senecio smallii", rep(NA_character_, 12L)),
  concept_synonymy_status = c("Homotypic_Synonym", "Heterotypic_Synonym", rep("Homotypic_Synonym", 12L)),
  anchor_tier_source = c("Tier_0_Type_Anchor", "Passive_Heterotypic_Type_Audit", rep("Tier_1_Gold_Core", 6L), rep("Tier_2_Silver_Regional", 6L)),
  is_passive = c(FALSE, TRUE, rep(FALSE, 12L)),
  decimalLatitude = c(35.8, 36.2, seq(35.75, 35.85, length.out = 6L), seq(36.15, 36.45, length.out = 6L)),
  decimalLongitude = c(-79.0, -80.2, rep(-79.0, 6L), seq(-80.1, -80.6, length.out = 6L)),
  recordedBy = c("TypeColl", "Small", paste0("CoreC", 1:6), rep(c("Smith", "Jones", "Brown"), 2L)),
  recordNumber = as.character(1:14),
  institutionCode = rep(c("NCU", "DUKE", "US"), length.out = 14L),
  phenological_state = rep("anthesis", 14L),
  doy = c(115, 160, rep(115, 6L), rep(160, 6L)),
  delta_doy_resid = c(0.0, 42.0, rep(0.5, 6L), rep(42.0, 6L)),
  allochronic_bimodality_delta_bic = c(0.2, 6.5, rep(0.2, 6L), rep(6.5, 6L)),
  PC1 = c(-2.0, 3.5, seq(-2.1, -1.9, length.out = 6L), seq(3.2, 3.8, length.out = 6L)),
  PC2 = c(0.0, 0.2, rep(0.0, 6L), rep(0.2, 6L)),
  can1 = c(-2.0, 3.5, seq(-2.1, -1.9, length.out = 6L), seq(3.2, 3.8, length.out = 6L)),
  can2 = c(0.0, 0.2, rep(0.0, 6L), rep(0.2, 6L)),
  cda_predicted_taxon = rep("Packera anonyma", 14L),
  cda_posterior_prob = c(0.95, 0.42, rep(0.94, 6L), rep(0.45, 6L)),
  self_morph_mahal_d2 = c(0.5, 29.0, rep(1.0, 6L), rep(28.0, 6L)),
  type_morph_mahal_d2 = c(0.5, 29.0, rep(1.0, 6L), rep(28.0, 6L)),
  nearest_type_basionym = c("Senecio anonymus", "Senecio smallii", rep("Senecio anonymus", 6L), rep("Senecio smallii", 6L)),
  nearest_type_basionym_d2 = c(0.5, 0.8, rep(1.0, 6L), rep(1.5, 6L)),
  nearest_type_niche_basionym = c("Senecio anonymus", "Senecio smallii", rep("Senecio anonymus", 6L), rep("Senecio smallii", 6L)),
  nearest_type_basionym_edaphic_z = c(0.3, 0.6, rep(0.4, 6L), rep(1.2, 6L)),
  edaphic_z_score = c(0.3, 3.4, rep(0.4, 6L), rep(3.3, 6L)),
  edaphic_niche_overlap_d = c(0.85, 0.15, rep(0.85, 6L), rep(0.16, 6L)),
  gmm_cluster = c(1L, 2L, rep(1L, 6L), rep(2L, 6L)),
  morphometric_outlier_flag = c(FALSE, TRUE, rep(FALSE, 6L), rep(TRUE, 6L)),
  morph_tier_source = rep("Basal_Rosette", 14L),
  basal_leaf_count = rep(3L, 14L),
  stringsAsFactors = FALSE
)
res_intra <- apply_taxonomic_decision_matrix(df_intra, intraspecific_mode = TRUE, gmm_two_delta_bic = 16.0, opts = list(demo = TRUE))

# 3. Direct evaluate_step07_hydro_edaphic_and_basionym_guards(): concordant heterotypic resurrection vs offline homotypic lock
df_guard <- data.frame(
  catalogNumber = c("G_HET_RES", "G_HOMO_OFFLINE", "G_SERP", "G_FEN"),
  species_standardized = c("Packera anonyma", "Packera aurea", "Packera serpenticola", "Packera aurea"),
  resurrection_basionym = c("Senecio smallii", "Senecio aureus", NA_character_, NA_character_),
  nearest_type_basionym = c("Senecio smallii", "Senecio aureus", NA_character_, NA_character_),
  nearest_type_niche_basionym = c("Senecio smallii", "Senecio aureus", NA_character_, NA_character_),
  nearest_type_basionym_edaphic_z = c(1.1, 0.8, NA_real_, NA_real_),
  heterotypic_synonyms_json = c("[\"Senecio smallii\"]", "[]", NA_character_, NA_character_),
  sgmc_lithology_class = c("granite", "shale", "ultramafic_serpentine", "alluvium"),
  soil_ca_mg_molar_ratio = c(0.8, 1.5, 0.20, 1.4),
  soil_ultramafic_ni_cr_index = c(1.0, 1.0, 4.5, 1.0),
  is_nwi_wetland = c(0L, 0L, 0L, 1L),
  wetland_seepage_index = c(0.05, 0.20, 0.05, 0.92),
  nwpl_wetland_score = c(0.10, 0.75, 0.05, 0.90),
  stringsAsFactors = FALSE
)
res_guard <- evaluate_step07_hydro_edaphic_and_basionym_guards(df_guard, db_conn = NULL)

pick_row <- function(d, cid) {
  k <- which(d$catalogNumber == cid)[1L]
  list(
    call = as.character(d$taxonomic_status_call[k]),
    sub = as.character(d$taxonomic_treatment_subtype[k]),
    rec = as.character(d$recommended_determination[k]),
    act = as.character(d$synthesis_triage_action[k]),
    res_bas = if (is.na(d$resurrection_basionym[k])) NULL else as.character(d$resurrection_basionym[k])
  )
}
idx_of <- function(cid) which(df_mt$catalogNumber == cid)[1L]
write_json(list(
  anch_s7 = list(
    homo = isTRUE(anch_step07[idx_of("TYPE_HOMO_DUB")]),
    het = isTRUE(anch_step07[idx_of("TYPE_HET_DIV")]),
    unc = isTRUE(anch_step07[idx_of("TYPE_UNC_COORD")]),
    open = isTRUE(anch_step07[idx_of("OPEN_AFF_CONCEPT")]),
    hyb = isTRUE(anch_step07[idx_of("HYB_FORMULA_CONCEPT")])
  ),
  anch_comp = list(
    homo = isTRUE(anch_comp[idx_of("TYPE_HOMO_DUB")]),
    het = isTRUE(anch_comp[idx_of("TYPE_HET_DIV")]),
    unc = isTRUE(anch_comp[idx_of("TYPE_UNC_COORD")]),
    open = isTRUE(anch_comp[idx_of("OPEN_AFF_CONCEPT")]),
    hyb = isTRUE(anch_comp[idx_of("HYB_FORMULA_CONCEPT")])
  ),
  mt_homo = pick_row(res_mt, "TYPE_HOMO_DUB"),
  mt_het = pick_row(res_mt, "TYPE_HET_DIV"),
  mt_hyb = pick_row(res_mt, "HYB_FORMULA_CONCEPT"),
  mt_bc = pick_row(res_mt, "NON_TYPE_RULE_BC_REASSIGN"),
  mt_prot = pick_row(res_mt, "NON_TYPE_PROTECTED"),
  mt_weak_rf = pick_row(res_mt, "NON_TYPE_WEAK_RF"),
  intra_homo = pick_row(res_intra, "INTRA_HOMO_TYPE"),
  intra_het = pick_row(res_intra, "INTRA_HET_TYPE"),
  intra_k2_res = pick_row(res_intra, "INTRA_K2_01"),
  guard_het_type = as.character(res_guard$basionym_synonymy_type[1L]),
  guard_het_bas = if (is.na(res_guard$resurrection_basionym[1L])) NULL else as.character(res_guard$resurrection_basionym[1L]),
  guard_het_niche_corrob = isTRUE(res_guard$basionym_niche_corroborated[1L]),
  guard_homo_type = as.character(res_guard$basionym_synonymy_type[2L]),
  guard_homo_status = as.character(res_guard$basionym_guard_status[2L]),
  guard_homo_bas = if (is.na(res_guard$resurrection_basionym[2L])) NULL else as.character(res_guard$resurrection_basionym[2L]),
  guard_regimes = as.character(res_guard$hydro_edaphic_niche_regime),
  guard_wfc_fen = as.numeric(res_guard$wetland_fidelity_concordance[4L])
), out_json, auto_unbox = TRUE, na = "null", null = "null")
"""


@unittest.skipIf(SKIP_R_COND, "Skipping warm-R type/basionym suite under fast audit or missing Rscript.")
class TestStep07TypeAndBasionymRules(unittest.TestCase):
    """Contract tests for Step 07 type control Rules A/B/C, heterotypic synonymy, and basionym resurrection guards."""

    res: Dict[str, Any] = {}

    @classmethod
    def setUpClass(cls) -> None:
        with tempfile.TemporaryDirectory(prefix="step07_type_bas_") as td:
            tmp = Path(td)
            r_p, out_j = tmp / "probe.R", tmp / "out.json"
            r_p.write_text(_WARM_R_CODE, encoding="utf-8")
            proc = run_warm_r_cmd(["Rscript", str(r_p), str(out_j)], stage_key="step07", cwd=PROJECT_ROOT, skip_figures=True)
            if proc.returncode != 0:
                raise AssertionError(f"Warm R probe failed ({proc.returncode}):\n{proc.stdout}\n{proc.stderr}")
            cls.res = json.loads(out_j.read_text(encoding="utf-8"))

    def test_rule_a_homotypic_control_and_heterotypic_synonymy_both_modes(self) -> None:
        r = self.res
        self.assertEqual(r["mt_homo"]["call"], "Confirmed_Type_Control")
        self.assertEqual(r["mt_homo"]["sub"], "Confirmed_Species")
        self.assertEqual(r["mt_het"]["call"], "Type_Synonymy_Candidate")
        self.assertEqual(r["mt_het"]["sub"], "Type_Synonymy_Candidate")
        self.assertEqual(r["mt_het"]["res_bas"], "Senecio smallii")
        self.assertEqual(r["intra_homo"]["call"], "Confirmed_Type_Control")
        self.assertEqual(r["intra_het"]["call"], "Type_Synonymy_Candidate")
        self.assertEqual(r["intra_het"]["sub"], "Type_Synonymy_Candidate")
        self.assertEqual(r["intra_het"]["res_bas"], "Senecio smallii")
        self.assertEqual(r["intra_k2_res"]["call"], "Taxonomic_Resurrection_Candidate")
        self.assertEqual(r["intra_k2_res"]["res_bas"], "Senecio smallii")

    def test_passive_type_and_open_concept_anchor_quarantine_and_hybrid_routing(self) -> None:
        r = self.res
        for key in ("anch_s7", "anch_comp"):
            self.assertTrue(r[key]["homo"], f"Homotypic type must remain anchor in {key}")
            for blocked in ("het", "unc", "open", "hyb"):
                self.assertFalse(r[key][blocked], f"{blocked} must be quarantined from {key}")
        self.assertEqual(r["mt_hyb"]["sub"], "Hybrid_Introgressant")
        self.assertIn(" x ", r["mt_hyb"]["rec"])

    def test_rule_bc_concordance_and_offline_homotypic_vs_heterotypic_basionym_guards(self) -> None:
        r = self.res
        self.assertEqual(r["mt_bc"]["sub"], "Misidentification_Reassign")
        self.assertEqual(r["mt_bc"]["rec"], "Packera obovata")
        self.assertEqual(r["mt_prot"]["sub"], "Confirmed_Species")
        self.assertNotEqual(r["mt_weak_rf"]["sub"], "Misidentification_Reassign")
        self.assertEqual(r["guard_het_type"], "Heterotypic")
        self.assertEqual(r["guard_het_bas"], "Senecio smallii")
        self.assertTrue(bool(r["guard_het_niche_corrob"]))
        self.assertEqual(r["guard_homo_type"], "Homotypic")
        self.assertEqual(r["guard_homo_status"], "Homotypic_Sensu_Stricto_Lock")
        self.assertIsNone(r["guard_homo_bas"])
        self.assertEqual(r["guard_regimes"][2], "Ultramafic_Serpentine_Barren")
        self.assertEqual(r["guard_regimes"][3], "Palustrine_Wetland_Fen")
        self.assertGreaterEqual(float(r["guard_wfc_fen"]), 0.95)


if __name__ == "__main__":
    unittest.main()
