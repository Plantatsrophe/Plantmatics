#!/usr/bin/env python3
"""
Cross-stage publication figure and visualization audit suite (Steps 01-07).

Verifies Python/R 24-hue Tol-Okabe palette parity, ICN botanical typography
(including flush nothospecies multiplication signs and Roman ranks/authorities),
double-ringed homotypic vs. heterotypic type symbols, US+Canada EPSG:5070
basemaps and centroid quarantine, 300-DPI PNG-only rendering across fig01 and
fig03a..fig07b plus Figure_Integrative_Packera_dubia_Revision.png with zero PDF
leaks, Ecotype_* Kruskal-Wallis p suppression, and interactive HTML overlay links.
"""

from __future__ import annotations

import json
import os
import struct
import sys
import tempfile
import unittest
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from scripts.tests._r_session import get_warm_r_session, rscript_available
from scripts.vision.step01_figures import generate_step01_cohort_audit_figure
from scripts.vision.step02_fig_style_assets import (
    HOLLOW_MARKERS_5,
    ICN_ROMAN_TOKENS,
    TOL_OKABE_12,
    TOL_OKABE_24,
    format_icn_taxon_mathtext,
)


def _read_png_meta(png_path: Path) -> Tuple[int, int, Optional[float]]:
    """Parse PNG dimensions (w, h) and DPI from IHDR/pHYs chunks without external deps."""
    raw = png_path.read_bytes()
    if raw[:8] != b"\x89PNG\r\n\x1a\n":
        raise AssertionError(f"Invalid PNG signature in {png_path}")
    w, h = struct.unpack(">II", raw[16:24])
    dpi: Optional[float] = None
    pos = 8
    while pos + 8 <= len(raw):
        length = struct.unpack(">I", raw[pos : pos + 4])[0]
        ctype = raw[pos + 4 : pos + 8]
        data = raw[pos + 8 : pos + 8 + length]
        if ctype == b"pHYs" and len(data) == 9:
            x_ppu, _y_ppu, unit = struct.unpack(">IIB", data)
            if unit == 1 and x_ppu > 0:
                dpi = x_ppu * 0.0254
            break
        pos += 12 + length
    return w, h, dpi


_R_AUDIT_SCRIPT = r"""
suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(jsonlite)
})
args <- commandArgs(trailingOnly = TRUE)
run_dir <- args[1L]
json_out <- args[2L]
fig_dir <- file.path(run_dir, "outputs", "figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

source("scripts/figures/theme_packera_pub.R")
source("scripts/morphometrics/pub_taxon_label_utils.R")
source("scripts/figures/step03_figures.R")
source("scripts/figures/step04_figures.R")
source("scripts/figures/step05_figures.R")
source("scripts/figures/step05_sdm_figures.R")
source("scripts/figures/step06_figures.R")
source("scripts/figures/step06c_figures.R")
source("scripts/figures/step07_figures.R")
source("scripts/figures/step07b_network_figures.R")
source("scripts/morphometrics/pub_synthesis_plates.R")
source("scripts/analysis/step07_reports_html_utils.R")

deparse1_safe <- function(x) paste(deparse(x), collapse = "")
icn_exprs <- list(
  notho_spaced = deparse1_safe(format_icn_label_expr("Packera \u00d7 memmingeri")),
  notho_flush = deparse1_safe(format_icn_label_expr("Packera \u00d7memmingeri")),
  notho_abbr = deparse1_safe(format_icn_label_expr("Packera \u00d7memmingeri", abbreviate_genus = TRUE)),
  sect_roman = deparse1_safe(format_icn_label_expr("Packera sect. Aurei")),
  var_auth = deparse1_safe(format_icn_label_expr("Packera paupercula var. savannarum R.R.Kowal")),
  subsp_sp_nov = deparse1_safe(format_icn_label_expr("Packera dubia subsp. dubia sp. nov.")),
  title_notho = deparse1_safe(format_icn_title_expr("A. ", "Packera \u00d7memmingeri", " Revision", bold_all = TRUE))
)

type_df <- data.frame(
  x = c(0.1, 0.2, 0.3, 0.4),
  y = c(0.1, 0.2, 0.3, 0.4),
  decimalLongitude = c(-79.0, -78.5, -77.5, -76.5),
  decimalLatitude = c(35.9, 36.1, 44.5, 37.2),
  is_type_specimen = c(FALSE, TRUE, TRUE, TRUE),
  taxonomic_treatment_subtype = c("Confirmed_Species", "Confirmed_Type_Control", "Type_Synonymy_Candidate", "Confirmed_Type_Control"),
  taxonomic_status_call = c("Confirmed_Species", "Confirmed_Species", "Type_Synonymy_Candidate", "Confirmed_Species"),
  type_niche_status = c("Clean", "Clean", "Clean", "Passive_Type_Uncertain_Coord"),
  is_centroid_or_inst_coord = c(1L, 0L, 0L, 1L),
  stringsAsFactors = FALSE
)
tiers <- classify_type_plot_tier_r(type_df)
p0 <- ggplot2::ggplot(type_df, ggplot2::aes(x = .data$x, y = .data$y)) + ggplot2::geom_point()
p_no_unc <- add_nomenclatural_type_layer_gg(p0, type_df, "x", "y", show_uncertain_coord = FALSE)
p_with_unc <- add_nomenclatural_type_layer_gg(p0, type_df, "x", "y", show_uncertain_coord = TRUE)
lyr_shapes <- vapply(p_with_unc$layers[-1L], function(l) as.integer(l$aes_params$shape), integer(1))
lyr_cols <- vapply(p_with_unc$layers[-1L], function(l) as.character(l$aes_params$colour), character(1))

occ_filt <- filter_map_occurrences_r(type_df)
bm <- build_na_projected_basemap_gg(c(-82.0, -75.0), c(35.0, 45.5))
bg_fill <- tryCatch(as.character(bm$plot$theme$panel.background$fill), error = function(e) "")

set.seed(42)
n_sp <- 16L
sp_vec <- rep(c("Packera dubia", "Packera anonyma"), each = 8L)
base_df <- data.frame(
  catalogNumber = sprintf("NCU_%04d", seq_len(n_sp)),
  plant_individual_id = rep(1L, n_sp),
  species_standardized = sp_vec,
  scientificName = sp_vec,
  regional_group = rep(c("Piedmont", "Blue_Ridge", "Southeastern_Plains", "Ridge_and_Valley"), length.out = n_sp),
  institutionCode = rep(c("NCU", "DUKE", "US", "NY"), length.out = n_sp),
  determiner_tier = rep(c("Tier_1_Gold", "Tier_2_Silver"), length.out = n_sp),
  authority_weight = rep(1.0, n_sp),
  is_type_specimen = seq_len(n_sp) %in% c(1L, 9L, 10L),
  taxonomic_treatment_subtype = ifelse(seq_len(n_sp) == 1L, "Confirmed_Type_Control", ifelse(seq_len(n_sp) == 9L, "Type_Synonymy_Candidate", "Confirmed_Species")),
  taxonomic_status_call = ifelse(seq_len(n_sp) == 9L, "Type_Synonymy_Candidate", "Confirmed_Species"),
  type_niche_status = ifelse(seq_len(n_sp) == 10L, "Passive_Type_Uncertain_Coord", "Clean"),
  is_centroid_or_inst_coord = as.integer(seq_len(n_sp) %in% c(2L, 10L)),
  decimalLongitude = seq(-82.5, -75.5, length.out = n_sp),
  decimalLatitude = seq(34.5, 44.5, length.out = n_sp),
  PC1 = c(rnorm(8, -1.2, 0.3), rnorm(8, 1.2, 0.3)),
  PC2 = rnorm(n_sp, 0, 0.4),
  PC3 = rnorm(n_sp, 0, 0.3),
  ECT_PC1 = c(rnorm(8, -0.8, 0.25), rnorm(8, 0.8, 0.25)),
  ECT_PC2 = rnorm(n_sp, 0, 0.3),
  log_centsize = seq(2.0, 3.2, length.out = n_sp),
  is_calibrated_mm = rep(TRUE, n_sp),
  can1 = c(rnorm(8, -1.4, 0.3), rnorm(8, 1.4, 0.3)),
  can2 = rnorm(n_sp, 0, 0.35),
  solidity = seq(0.55, 0.92, length.out = n_sp),
  AR_inv = seq(1.1, 2.2, length.out = n_sp),
  ect_lobe_complexity = seq(0.2, 0.8, length.out = n_sp),
  cauline_ect_lobe_complexity = seq(0.3, 0.9, length.out = n_sp),
  blade_cie_l_sheet_delta = rnorm(n_sp, -10, 2),
  abaxial_cie_l_sheet_delta = rnorm(n_sp, -6, 2),
  log1p_capitulum_count = runif(n_sp, 1.0, 2.8),
  lithic_outcrop_index = seq(0.05, 0.90, length.out = n_sp),
  wetland_seepage_index = seq(0.85, 0.10, length.out = n_sp),
  fine_soil_ph = seq(4.8, 6.8, length.out = n_sp),
  bio1_temp_mean = seq(10.0, 17.5, length.out = n_sp),
  doy = seq(105, 165, length.out = n_sp),
  hopkins_adjusted_doy = seq(115, 155, length.out = n_sp),
  hopkins_adjusted_doy_resid = seq(-18, 18, length.out = n_sp),
  anthesis_dominance_weight = rep(0.85, n_sp),
  phenology_state = rep("Peak_Anthesis", n_sp),
  gmm_cluster = rep(c(1L, 2L), each = 8L),
  probabilistic_consensus_score = runif(n_sp, 0.72, 0.96),
  morph_tier_source = rep("Basal_Rosette", n_sp),
  final_taxonomic_decision = rep("Confirmed_Species", n_sp),
  stringsAsFactors = FALSE
)

# Step 03 plates
hmat <- matrix(rnorm(n_sp * 24L, 0, 0.08), nrow = n_sp, ncol = 24L)
colnames(hmat) <- c(paste0("A", 1:6), paste0("B", 1:6), paste0("C", 1:6), paste0("D", 1:6))
render_step03_figures_suite(base_df, base_df, harm_mat = hmat, figures_dir = fig_dir)

# Step 04 plates
render_step04_allometry_png(base_df, r2_size = 0.09, p_size = 0.04, r2_int = 0.02, p_int = 0.45,
                            r2_within = 0.07, p_within = 0.03, scale_tier = "verified_ruler", corrected = TRUE,
                            output_path = file.path(fig_dir, "fig04a_foliar_allometry_audit.png"))
gmm_mock <- list(best_k = 2L, bic_table = data.frame(K = 1:3, BIC = c(-80, -60, -62), Delta_BIC_prev = c(0, 20, -2),
                 Two_Delta_BIC = c(0, 20, 18), optimal_K_raw = 3L, optimal_K_merged = 2L,
                 Delta_ICL_vs_Null = c(0, 14.5, 10.2), Boot_LRT_P_prev = c(NA_real_, 0.005, 0.42)))
render_step04_gmm_cda_png(base_df, cda_res = list(loocv_accuracy = 0.91), gmm_res = gmm_mock,
                          output_path = file.path(fig_dir, "fig04b_gmm_cda_morphometrics.png"))
render_step04_ecoregion_batch_png(base_df, eco_cda_res = list(loocv_accuracy = 0.84), figures_dir = fig_dir,
                                  output_path = file.path(fig_dir, "fig04c_ecoregion_and_batch_audit.png"))

# Step 05 plates (use Ecotype_* on fig05a to verify zero KW p leakage)
eco_df <- base_df
eco_df$analysis_group <- rep(c("Ecotype_1", "Ecotype_2"), each = 8L)
srf_mock <- list(importance = data.frame(variable = c("lithic_outcrop_index", "wetland_seepage_index", "fine_soil_ph", "bio1_temp_mean"),
                 importance = c(0.18, 0.14, 0.09, 0.06), stringsAsFactors = FALSE))
render_step05a_ecological_niche_png(eco_df, srf_res = srf_mock, out_plot_path = file.path(fig_dir, "fig05a_ecological_niche_and_rf_importance.png"))
render_step05b_hopkins_phenology_png(base_df, baseline_info = list(slope = 3.6, selected_model = "M2_Latitude_Elevation"),
                                     out_pheno_path = file.path(fig_dir, "fig05b_hopkins_circular_phenology.png"),
                                     bimod_info = list(allochronic_bimodality_delta_bic = 4.2, ashman_d = 2.4))
render_step05c_sdm_suitability_png(base_df, figures_dir = fig_dir, focal_taxon = "Packera dubia",
                                   out_sdm_path = file.path(fig_dir, "fig05c_sdm_suitability_and_distribution_map.png"))

# Step 06 plates
pass_df <- data.frame(genetic_sample_id = sprintf("S%02d", 1:6), species_standardized = rep(c("Packera dubia", "Packera anonyma"), each = 3L),
                      replicate_pair_id = c("R1", "R1", NA, NA, NA, NA), type_status = c("Holotype", rep("", 5L)),
                      is_putative_novel_or_hybrid_query = rep(0L, 6L), stringsAsFactors = FALSE)
assign_df <- data.frame(genetic_sample_id = pass_df$genetic_sample_id, geno_predicted_taxon = pass_df$species_standardized,
                        PCo1 = c(-0.3, -0.28, -0.25, 0.25, 0.3, 0.32), PCo2 = c(0.1, 0.12, -0.08, 0.05, -0.1, 0.02),
                        inferred_cytotype = c("2x", "2x", "2x", "4x", "4x", "4x"),
                        nquire_delta_logl = c(-6.5, -5.2, -1.0, 6.2, 7.1, 5.4), stringsAsFactors = FALSE)
q_mat <- cbind(c(0.95, 0.92, 0.88, 0.08, 0.05, 0.04), c(0.05, 0.08, 0.12, 0.92, 0.95, 0.96))
rownames(q_mat) <- pass_df$genetic_sample_id
colnames(q_mat) <- c("Packera dubia", "Packera anonyma")
tr_mock <- ape::read.tree(text = "((S01:0.1,S02:0.1)95:0.2,((S03:0.1,S04:0.1)40:0.1,(S05:0.1,S06:0.1)88:0.2)75:0.1);")
render_fig06a_genomic_pcoa_tree_png(file.path(fig_dir, "fig06a_genomic_pcoa_tree_concordance.png"),
                                    pass_df, assign_df, q_mat, c(54.2, 21.8), tr_mock, data.frame(), "", "test_run", focal_taxon = "Packera dubia")
tri_df <- data.frame(hybrid_index_h = c(0.05, 0.5, 0.95), interspecific_ancestry_q12 = c(0.06, 0.88, 0.08),
                     hybrid_generation_class = c("Parental_A", "F1_Hybrid", "Parental_B"), stringsAsFactors = FALSE)
gdi_df <- data.frame(taxon = c("Packera dubia", "Packera anonyma"), coalescent_gdi = c(0.76, 0.81),
                     gdi_ci_low = c(0.65, 0.72), gdi_ci_high = c(0.85, 0.89),
                     gdi_classification = c("Distinct_Species", "Distinct_Species"),
                     bf_method = c("bpp_stepping_stone", "bpp_stepping_stone"), stringsAsFactors = FALSE)
render_fig06b_reticulation_gdi_png(file.path(fig_dir, "fig06b_reticulation_gdi_landscape.png"),
                                   tri_df, data.frame(), gdi_df, data.frame(), data.frame(), focal_taxon = "Packera dubia")
render_fig06c_cytonuclear_paralog_cytotype_png(file.path(fig_dir, "fig06c_cytonuclear_paralog_cytotype_audit.png"),
                                               parsed_dir = file.path(run_dir, "missing_parsed"), out_dir = file.path(run_dir, "missing_out"),
                                               full_audit_df = merge(pass_df, assign_df, by = "genetic_sample_id"))

# Step 07 + focal revision plates + interactive HTML
render_step07a_integrative_synthesis_png(base_df, out_plot_path = file.path(fig_dir, "fig07a_integrative_taxonomic_synthesis.png"), focal_taxon = "Packera dubia")
W <- matrix(c(1.0, 0.35, 0.12, 0.35, 1.0, 0.42, 0.12, 0.42, 1.0), nrow = 3L,
            dimnames = list(c("Packera dubia", "Packera anonyma", "Packera aurea"), c("Packera dubia", "Packera anonyma", "Packera aurea")))
cdf <- data.frame(species_standardized = rownames(W), n_vouchers = c(8L, 8L, 6L), consensus_complex_id = c(1L, 1L, 2L), stringsAsFactors = FALSE)
render_step07b_complex_network_png(cdf, W, out_complex_path = file.path(fig_dir, "fig07b_taxon_complex_delimitation_network.png"))
render_focal_revision_synthesis_plate_png(base_df[base_df$species_standardized == "Packera dubia", , drop = FALSE],
                                          file.path(fig_dir, "Figure_Integrative_Packera_dubia_Revision.png"),
                                          full_df = base_df, focal_taxon = "Packera dubia")
html_path <- file.path(run_dir, "outputs", "reports", "Triage_Dashboard_Interactive.html")
dir.create(dirname(html_path), recursive = TRUE, showWarnings = FALSE)
generate_interactive_triage_html(base_df, html_path, focal_taxon = "Packera dubia")

res <- list(
  pal24 = unname(PACKERA_TOL_OKABE_24),
  pal12 = unname(PACKERA_TOL_OKABE_12),
  icn = icn_exprs,
  tiers = tiers,
  n_layers_no_unc = length(p_no_unc$layers),
  n_layers_with_unc = length(p_with_unc$layers),
  lyr_shapes = lyr_shapes,
  lyr_cols = lyr_cols,
  n_clean = nrow(occ_filt$clean_df),
  n_unc_type = nrow(occ_filt$uncertain_type_df),
  n_quar_non_type = occ_filt$n_quarantined_non_type,
  bm_projected = isTRUE(bm$projected),
  bm_crs = as.integer(bm$crs),
  bm_has_bbox = !is.null(bm$bbox_proj),
  bm_water_bg = bg_fill,
  html_path = html_path
)
jsonlite::write_json(res, json_out, auto_unbox = TRUE)
"""


class TestPubFigureAuditSuite(unittest.TestCase):
    """End-to-end verification of the upgraded 300-DPI PNG figure suite across Steps 01-07."""

    _tmp_dir: tempfile.TemporaryDirectory[str]
    _run_dir: Path
    _fig_dir: Path
    _r_data: Optional[Dict[str, Any]] = None

    @classmethod
    def setUpClass(cls) -> None:
        cls._tmp_dir = tempfile.TemporaryDirectory(prefix="plantmatics_pub_fig_audit_")
        cls._run_dir = Path(cls._tmp_dir.name) / "run_audit"
        cls._fig_dir = cls._run_dir / "outputs" / "figures"
        cls._fig_dir.mkdir(parents=True, exist_ok=True)
        cls._run_step01_fixture()
        if rscript_available():
            json_out = Path(cls._tmp_dir.name) / "r_audit.json"
            proc = get_warm_r_session("step07").run_script(
                _R_AUDIT_SCRIPT,
                cli_args=[str(cls._run_dir), str(json_out)],
                skip_figures=False,
            )
            if proc.returncode != 0 or not json_out.exists():
                raise AssertionError(f"Warm-R figure audit failed (rc={proc.returncode}):\nSTDOUT:\n{proc.stdout}\nSTDERR:\n{proc.stderr}")
            cls._r_data = json.loads(json_out.read_text(encoding="utf-8"))

    @classmethod
    def tearDownClass(cls) -> None:
        cls._tmp_dir.cleanup()

    @classmethod
    def _run_step01_fixture(cls) -> None:
        v_csv = cls._run_dir / "data" / "tables" / "curated_vouchers.csv"
        v_csv.parent.mkdir(parents=True, exist_ok=True)
        pd.DataFrame(
            {
                "catalogNumber": ["NCU_01", "NCU_02", "NCU_03", "NCU_04", "NCU_05"],
                "species_standardized": ["Packera dubia", "Packera ×memmingeri", "Packera dubia", "Packera anonyma", "Packera anonyma"],
                "decimalLongitude": [-79.1, -82.4, -78.2, -77.1, -76.4],
                "decimalLatitude": [35.9, 35.6, 36.2, 44.6, 37.3],
                "determiner_tier": ["Tier_1_Gold", "Tier_2_Silver", "Tier_3_Bronze", "Tier_1_Gold", "Tier_2_Silver"],
                "is_type_specimen": [1, 0, 0, 1, 1],
                "taxonomic_treatment_subtype": ["Confirmed_Type_Control", "Confirmed_Species", "Confirmed_Species", "Type_Synonymy_Candidate", "Confirmed_Type_Control"],
                "type_niche_status": ["Clean", "Clean", "Clean", "Clean", "Passive_Type_Uncertain_Coord"],
                "is_centroid_or_inst_coord": [0, 0, 1, 0, 1],
                "regional_group": ["Piedmont", "Blue_Ridge", "Piedmont", "Northern_Lakes", "Southeastern_Plains"],
                "year": [1965, 1988, 2002, 2011, 2019],
                "sharpness_score": [120.0, 145.0, 98.0, 210.0, 175.0],
            }
        ).to_csv(v_csv, index=False)
        out_p = generate_step01_cohort_audit_figure(v_csv)
        assert out_p is not None and out_p.exists()

    def test_python_palette_icn_and_fig01_300dpi_png(self) -> None:
        self.assertEqual(len(TOL_OKABE_24), 24)
        self.assertEqual(len(set(TOL_OKABE_24)), 24)
        self.assertEqual(TOL_OKABE_12, TOL_OKABE_24[:12])
        self.assertEqual(len(HOLLOW_MARKERS_5), 5)
        for tok in ("sect.", "var.", "subsp.", "sp.", "nov."):
            self.assertIn(tok, ICN_ROMAN_TOKENS)
        self.assertIn("×$\\it{memmingeri}$", format_icn_taxon_mathtext("Packera × memmingeri"))
        self.assertIn("×$\\it{memmingeri}$", format_icn_taxon_mathtext("Packera ×memmingeri", abbreviate_genus=True))
        self.assertIn("var.", format_icn_taxon_mathtext("Packera paupercula var. savannarum"))
        self.assertIn("sect.", format_icn_taxon_mathtext("Packera sect. Aurei"))
        self.assertIn("sp.", format_icn_taxon_mathtext("Packera dubia sp. nov."))
        fig01 = self._fig_dir / "fig01_cohort_biogeography_and_authority.png"
        self.assertTrue(fig01.exists(), "Auto-resolved fig01 PNG missing")
        w, h, dpi = _read_png_meta(fig01)
        self.assertGreater(w, 1200)
        self.assertGreater(h, 800)
        if dpi is not None:
            self.assertGreater(dpi, 290.0)
        self.assertEqual(list(self._run_dir.rglob("*.pdf")), [])

    def test_r_palette_parity_icn_types_basemap_and_all_plates(self) -> None:
        if not rscript_available():
            self.skipTest("Rscript unavailable or fast Python-only audit mode active")
        assert self._r_data is not None
        d = self._r_data
        self.assertEqual(tuple(d["pal24"]), TOL_OKABE_24)
        self.assertEqual(tuple(d["pal12"]), TOL_OKABE_12)
        self.assertIn('paste("\u00d7", italic("memmingeri"))', d["icn"]["notho_spaced"])
        self.assertIn('paste("\u00d7", italic("memmingeri"))', d["icn"]["notho_flush"])
        self.assertIn('italic("P.")', d["icn"]["notho_abbr"])
        self.assertIn('paste("\u00d7", italic("memmingeri"))', d["icn"]["notho_abbr"])
        self.assertIn('"sect."', d["icn"]["sect_roman"])
        self.assertIn('"var."', d["icn"]["var_auth"])
        self.assertIn('"R.R.Kowal"', d["icn"]["var_auth"])
        self.assertIn('"sp."', d["icn"]["subsp_sp_nov"])
        self.assertIn('"nov."', d["icn"]["subsp_sp_nov"])
        self.assertIn('bolditalic("memmingeri")', d["icn"]["title_notho"])

        self.assertEqual(d["tiers"], ["None", "Confirmed_Type_Control", "Type_Synonymy_Candidate", "Passive_Type_Uncertain_Coord"])
        self.assertEqual(d["n_layers_no_unc"], 5)
        self.assertEqual(d["n_layers_with_unc"], 7)
        self.assertEqual(d["lyr_shapes"], [23, 5, 24, 2, 1, 5])
        self.assertEqual(d["lyr_cols"], ["#332288", "#111827", "#882255", "#882255", "#64748B", "#64748B"])

        self.assertEqual(d["n_clean"], 2)
        self.assertEqual(d["n_unc_type"], 1)
        self.assertEqual(d["n_quar_non_type"], 1)
        self.assertTrue(d["bm_projected"])
        self.assertEqual(d["bm_crs"], 5070)
        self.assertTrue(d["bm_has_bbox"])
        self.assertEqual(d["bm_water_bg"], "#EEF4F8")

        expected_pngs: List[str] = [
            "fig01_cohort_biogeography_and_authority.png",
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
            "fig06c_cytonuclear_paralog_cytotype_audit.png",
            "fig07a_integrative_taxonomic_synthesis.png",
            "fig07b_taxon_complex_delimitation_network.png",
            "Figure_Integrative_Packera_dubia_Revision.png",
        ]
        for fname in expected_pngs:
            fpath = self._fig_dir / fname
            self.assertTrue(fpath.exists(), f"Missing rendered plate: {fname}")
            w, h, dpi = _read_png_meta(fpath)
            self.assertGreater(w, 800, f"Unexpectedly narrow PNG: {fname}")
            self.assertGreater(h, 500, f"Unexpectedly short PNG: {fname}")
            if dpi is not None:
                self.assertGreater(dpi, 290.0, f"Expected 300-DPI PNG for {fname}, got {dpi}")

        self.assertEqual(list(self._run_dir.rglob("*.pdf")), [])
        self.assertFalse((PROJECT_ROOT / "Rplots.pdf").exists())

        step05_src = (PROJECT_ROOT / "scripts" / "figures" / "step05_figures.R").read_text(encoding="utf-8")
        self.assertIn("kw_p <- NA_real_", step05_src)
        self.assertIn("^(Ecotype|Morphotype)_", step05_src)

        html_text = Path(d["html_path"]).read_text(encoding="utf-8")
        self.assertIn("_lm3_sheet_overlay.png", html_text)
        self.assertIn(".s-Confirmed_Type_Control", html_text)
        self.assertIn("fig01_cohort_biogeography_and_authority.png", html_text)
        self.assertIn("fig06c_cytonuclear_paralog_cytotype_audit.png", html_text)
        self.assertIn("Figure_Integrative_Packera_dubia_Revision.png", html_text)


if __name__ == "__main__":
    unittest.main()
