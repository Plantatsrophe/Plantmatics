#!/usr/bin/env python3
"""Warm-R & Python tests for Step 07 taxon-complex discovery, 4-stream genomic affinity, eco hygiene, and fig07b PNG."""

from __future__ import annotations

import json
import os
import shutil
import sys
import tempfile
import unittest
from pathlib import Path
from typing import Any, Dict
from unittest import mock

PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

import scripts.core.harvester as _harv
from scripts.core.harvester import standardize_packera_taxon
from scripts.tests._r_session import run_warm_r_cmd

SKIP_R_COND = (
    os.environ.get("PLANTMATICS_FAST_AUDIT") == "1"
    or os.environ.get("PACKERA_FAST_AUDIT") == "1"
    or shutil.which("Rscript") is None
)

_WARM_R_CODE = r"""
if (grDevices::dev.cur() == 1L) grDevices::pdf(NULL)
args <- commandArgs(trailingOnly = TRUE)
out_json <- args[1L]; work_dir <- args[2L]
suppressPackageStartupMessages({ library(jsonlite) })
source("scripts/morphometrics/pub_figures_utils.R")
source("scripts/analysis/07_triage_dashboard_synthesis.R")
source("scripts/figures/step07b_network_figures.R")

# 1. Pure POSIX partition_taxon_complexes() on mixed Packera & P. [K#] labels + RNG preservation
set.seed(123L)
seed_before <- .Random.seed
lbl_mix <- c("Packera dubia", "P. anonyma [K1]", "P. anonyma [K2]", "Packera aurea", "P. obovata [K1]", "Packera paupercula")
W_mix <- matrix(0.05, nrow = 6L, ncol = 6L, dimnames = list(lbl_mix, lbl_mix))
W_mix[1:3, 1:3] <- 0.85; W_mix[4:5, 4:5] <- 0.85; diag(W_mix) <- 0.0
nv_m <- stats::setNames(rep(12L, 6L), lbl_mix); ng_m <- stats::setNames(rep(8L, 6L), lbl_mix)
part_mix <- partition_taxon_complexes(W_mix, lbl_mix, nv_m, ng_m)
seed_after <- .Random.seed
rng_preserved <- identical(seed_before, seed_after)

# 2. 8-species 3-complex discovery + fig07b 300-DPI PNG + zero PDFs
set.seed(42L)
sp8 <- c("Packera dubia", "Packera anonyma", "Packera millefolium",
         "Packera aurea", "Packera obovata", "Packera schweinitziana",
         "Packera paupercula", "Packera plattensis")
centers8 <- list(
  "Packera dubia" = c(-5.0, 0.0, 5.2, 115), "Packera anonyma" = c(-4.6, 0.2, 5.4, 118), "Packera millefolium" = c(-4.8, -0.2, 5.3, 120),
  "Packera aurea" = c(0.0, 5.0, 6.8, 130), "Packera obovata" = c(0.3, 4.7, 6.9, 133), "Packera schweinitziana" = c(-0.2, 5.2, 6.7, 132),
  "Packera paupercula" = c(5.0, -4.0, 7.6, 155), "Packera plattensis" = c(5.2, -3.8, 7.8, 158)
)
rows8 <- do.call(rbind, lapply(sp8, function(sp) {
  cc <- centers8[[sp]]; n_i <- 8L
  data.frame(
    catalogNumber = sprintf("%s_%02d", gsub("[^A-Za-z]", "", sp), seq_len(n_i)),
    species_standardized = rep(sp, n_i),
    determiner_tier = rep("Tier_1_Gold", n_i),
    morph_tier_source = rep("Basal_Rosette", n_i),
    basal_leaf_count = rep(3L, n_i),
    PC1 = rnorm(n_i, cc[1L], 0.25), PC2 = rnorm(n_i, cc[2L], 0.25), PC3 = rnorm(n_i, 0.0, 0.15),
    can1 = rnorm(n_i, cc[1L], 0.25), can2 = rnorm(n_i, cc[2L], 0.25),
    cda_predicted_taxon = rep(sp, n_i), cda_posterior_prob = rep(0.92, n_i),
    edaphic_best_fit_taxon = rep(sp, n_i), rf_predicted_taxon = rep(sp, n_i), rf_predicted_prob = rep(0.90, n_i),
    fine_soil_ph = rnorm(n_i, cc[3L], 0.15),
    lithic_outcrop_index = rnorm(n_i, if (cc[1L] < -2) 0.8 else 0.2, 0.05),
    wetland_seepage_index = rnorm(n_i, if (cc[2L] > 2) 0.8 else 0.1, 0.05),
    summer_moisture_deficit_idx = rnorm(n_i, if (cc[1L] > 2) 0.75 else 0.25, 0.05),
    maritime_continentality_idx = rnorm(n_i, if (cc[1L] > 2) 0.80 else 0.35, 0.05),
    decimalLatitude = rnorm(n_i, 36.0, 0.2), decimalLongitude = rnorm(n_i, -80.0, 0.2),
    phenological_state = rep("anthesis", n_i),
    doy = rnorm(n_i, cc[4L], 2.0), hopkins_adjusted_doy = rnorm(n_i, cc[4L], 2.0),
    gmm_cluster = rep(1L, n_i),
    stringsAsFactors = FALSE
  )
}))
c_dir <- file.path(work_dir, "complex_out")
dir.create(c_dir, recursive = TRUE, showWarnings = FALSE)
opts_c <- list(
  focal_taxon = "Packera dubia", min_affinity = 0.03, demo = TRUE,
  output_complexes = file.path(c_dir, "empirical_taxon_complexes.csv"),
  output_affinity = file.path(c_dir, "pairwise_taxon_affinity_matrix.csv"),
  output_neighbors = file.path(c_dir, "focal_auto_neighbors.txt"),
  output_complex_plot = file.path(c_dir, "fig07b_taxon_complex_delimitation_network.png")
)
comp8 <- run_step07_complex_discovery(rows8, opts_c)
cdf8 <- comp8$complexes_df
cid_map <- stats::setNames(as.integer(cdf8$consensus_complex_id), cdf8$species_standardized)

# 3. Subcluster mode (2 species with K1/K2 -> TRUE preserving base_species) vs multi-species suppression (4 species -> FALSE)
df_sub2 <- rows8[rows8$species_standardized %in% c("Packera dubia", "Packera anonyma"), , drop = FALSE]
df_sub2$gmm_cluster <- rep(c(rep(1L, 5L), rep(2L, 3L)), 2L)
df_sub2 <- rbind(df_sub2, df_sub2[df_sub2$gmm_cluster == 2L, , drop = FALSE])
sub2_res <- build_taxon_morphotype_subclusters(df_sub2, min_vouchers = 5L)
df_sub4 <- rows8[rows8$species_standardized %in% sp8[1:4], , drop = FALSE]
df_sub4$gmm_cluster <- rep(c(1L, 1L, 1L, 1L, 2L, 2L, 2L, 2L), 4L)
sub4_res <- build_taxon_morphotype_subclusters(df_sub4, min_vouchers = 4L)

# 4. 4-stream genomic affinity vs exact 3-stream invariance when genomics absent
df_4s <- rows8[rows8$species_standardized %in% c("Packera dubia", "Packera anonyma", "Packera aurea", "Packera paupercula"), , drop = FALSE]
t4 <- sort(unique(df_4s$species_standardized))
anc4 <- select_taxon_anchors(df_4s)
mr4 <- compute_morphometric_affinity(df_4s, anc4$anchor_df, t4)
er4 <- compute_ecological_affinity(df_4s, t4, anc4$anchor_mask)
pr4 <- compute_hopkins_phenological_affinity(df_4s, t4, anc4$anchor_mask)
kw3 <- compute_modality_kappa_weights(df_4s, anc4$anchor_mask)
gr_none <- compute_genomic_complex_affinity(df_4s, t4, anc4$anchor_mask)
W_3_null <- compute_multimodal_complex_affinity(mr4, er4, pr4, kw3, t4, geno_res = NULL)
W_3_none <- compute_multimodal_complex_affinity(mr4, er4, pr4, kw3, t4, geno_res = gr_none)
max_diff_no_geno <- max(abs(W_3_null - W_3_none))

df_4g <- df_4s
df_4g$geno_quality_pass <- TRUE
df_4g$geno_predicted_taxon <- df_4g$species_standardized
df_4g$coancestry_cluster <- ifelse(df_4g$species_standardized %in% c("Packera dubia", "Packera anonyma"), "C_DUB_ANO", ifelse(df_4g$species_standardized == "Packera aurea", "C_AUR", "C_PAU"))
df_4g$public_backbone_clade <- ifelse(df_4g$species_standardized %in% c("Packera dubia", "Packera anonyma"), "Clade_East_A", "Clade_East_B")
df_4g$coalescent_gdi <- ifelse(df_4g$species_standardized %in% c("Packera dubia", "Packera anonyma"), 0.25, 0.85)
kw4 <- compute_modality_kappa_weights(df_4g, anc4$anchor_mask)
gr_yes <- compute_genomic_complex_affinity(df_4g, t4, anc4$anchor_mask)
W_4_yes <- compute_multimodal_complex_affinity(mr4, er4, pr4, kw4, t4, geno_res = gr_yes)

# 5. Banned cartographic/buffer exclusion + smooth index inclusion in resolve_step07_complex_eco_predictors()
df_eco_hyg <- df_4s
n_eh <- nrow(df_eco_hyg); ph_idx <- seq_len(n_eh)
df_eco_hyg$lithic_outcrop_index <- 0.5 + 0.4 * sin(2 * pi * ph_idx / 4)
df_eco_hyg$wetland_seepage_index <- 0.5 + 0.4 * cos(2 * pi * ph_idx / 4)
df_eco_hyg$summer_moisture_deficit_idx <- 0.5 + 0.4 * sin(2 * pi * ph_idx / 7)
df_eco_hyg$maritime_continentality_idx <- 0.5 + 0.4 * cos(2 * pi * ph_idx / 11)
df_eco_hyg$ssurgo_rock_outcrop_pct <- seq(0, 100, length.out = n_eh)
df_eco_hyg$ssurgo_aws_0_50cm <- seq(1, 15, length.out = n_eh)
df_eco_hyg$ssurgo_lithic_series_pct <- seq(0, 90, length.out = n_eh)
df_eco_hyg$ultramafic_calcareous_gradient <- seq(-1, 1, length.out = n_eh)
df_eco_hyg$min_bedrock_depth_buffer_cm <- seq(5, 200, length.out = n_eh)
eco_hyg_cols <- resolve_step07_complex_eco_predictors(df_eco_hyg)

png_p <- opts_c$output_complex_plot
png_magic_ok <- file.exists(png_p) && identical(as.integer(readBin(png_p, "raw", 4L)), c(137L, 80L, 78L, 71L))
pdf_files <- list.files(c_dir, pattern = "\\.pdf$", recursive = TRUE, ignore.case = TRUE)

write_json(list(
  rng_preserved = rng_preserved,
  mix_labels = as.character(part_mix$complexes_df$complex_label),
  comp8_subcluster = isTRUE(comp8$subcluster_mode),
  cid_map = as.list(cid_map),
  n_complexes_8 = length(unique(cdf8$consensus_complex_id)),
  focal_neighbors = as.character(comp8$selected_neighbors),
  sub2_active = isTRUE(sub2_res$subcluster_mode),
  sub2_has_base_sp = "base_species" %in% names(sub2_res$morph_df),
  sub4_suppressed = !isTRUE(sub4_res$subcluster_mode),
  gr_none_avail = isTRUE(gr_none$available),
  max_diff_no_geno = max_diff_no_geno,
  gr_yes_avail = isTRUE(gr_yes$available),
  gr_yes_n_taxa = as.integer(gr_yes$n_taxa),
  geno_aff_sisters = as.numeric(gr_yes$mat["Packera dubia", "Packera anonyma"]),
  geno_aff_distant = as.numeric(gr_yes$mat["Packera dubia", "Packera paupercula"]),
  w4_sisters = as.numeric(W_4_yes["Packera dubia", "Packera anonyma"]),
  w4_distant = as.numeric(W_4_yes["Packera dubia", "Packera paupercula"]),
  w4_w_geno = as.numeric(attr(W_4_yes, "w_geno")),
  eco_hyg_cols = eco_hyg_cols,
  png_magic_ok = png_magic_ok,
  n_pdfs = length(pdf_files)
), out_json, auto_unbox = TRUE, na = "null")
"""


class TestStep07ComplexDiscovery(unittest.TestCase):
    """Contract tests for taxon standardization, 8-species complex discovery, 4-stream genomic affinity, and fig07b."""

    res: Dict[str, Any] = {}

    @classmethod
    def setUpClass(cls) -> None:
        if SKIP_R_COND:
            return
        with tempfile.TemporaryDirectory(prefix="step07_complex_") as td:
            tmp = Path(td)
            r_p, out_j = tmp / "probe.R", tmp / "out.json"
            r_p.write_text(_WARM_R_CODE, encoding="utf-8")
            proc = run_warm_r_cmd(["Rscript", str(r_p), str(out_j), str(tmp)], stage_key="step07", cwd=PROJECT_ROOT, skip_figures=False)
            if proc.returncode != 0:
                raise AssertionError(f"Warm R probe failed ({proc.returncode}):\n{proc.stdout}\n{proc.stderr}")
            cls.res = json.loads(out_j.read_text(encoding="utf-8"))

    def test_standardize_packera_taxon_eastern_species_and_basionyms(self) -> None:
        direct_cases = {
            "Packera dubia (Spreng.) Trock & Mabb.": "Packera dubia",
            "Packera x memmingeri": "Packera memmingeri",
            "Packera millefolia": "Packera millefolium",
            "Packera paupercula var. appalachiana": "Packera paupercula",
            "Packera serpenticola": "Packera serpenticola",
            "Packera glabella": "Packera glabella",
        }
        basionym_cases = {
            "Senecio tomentosus": "Packera dubia",
            "Senecio anonymus": "Packera anonyma",
            "Senecio smallii": "Packera anonyma",
            "Senecio aureus": "Packera aurea",
            "Senecio obovatus": "Packera obovata",
            "Senecio pauperculus": "Packera paupercula",
            "Senecio millefolium": "Packera millefolium",
            "Senecio plattensis": "Packera plattensis",
            "Senecio schweinitzianus": "Packera schweinitziana",
        }

        def _fake_backbone(name: str = "", rank: str = "species", **kwargs: Any) -> Dict[str, Any]:
            acc = basionym_cases.get(str(name).strip())
            if not acc:
                return {"matchType": "NONE"}
            return {
                "matchType": "EXACT",
                "status": "SYNONYM",
                "scientificName": name,
                "canonicalName": name,
                "species": acc,
                "accepted": acc,
                "usageKey": 1000,
                "acceptedUsageKey": 1000,
            }

        with (
            mock.patch.dict(_harv._DYNAMIC_SYNONYM_TO_CANONICAL, {}, clear=True),
            mock.patch.dict(_harv._GBIF_BACKBONE_SYNONYM_FLAG, {}, clear=True),
            mock.patch.object(_harv.species, "name_backbone", side_effect=_fake_backbone),
            mock.patch.object(_harv.species, "name_usage", return_value={"results": []}),
        ):
            for raw, expected in direct_cases.items():
                self.assertEqual(standardize_packera_taxon(raw), expected, f"Failed standardization for {raw}")
            for raw, expected in basionym_cases.items():
                self.assertEqual(standardize_packera_taxon(raw), expected, f"Failed standardization for {raw}")

    @unittest.skipIf(SKIP_R_COND, "Fast audit or Rscript unavailable")
    def test_posix_partitioning_8_species_3_complexes_and_subcluster_suppression(self) -> None:
        r = self.res
        self.assertTrue(bool(r["rng_preserved"]))
        self.assertTrue(all(lbl.startswith("Complex_") and "Packera " not in lbl and "P. " not in lbl for lbl in r["mix_labels"]))
        self.assertFalse(bool(r["comp8_subcluster"]))
        self.assertEqual(int(r["n_complexes_8"]), 3)
        cids = r["cid_map"]
        self.assertEqual(len({cids["Packera dubia"], cids["Packera anonyma"], cids["Packera millefolium"]}), 1)
        self.assertEqual(len({cids["Packera aurea"], cids["Packera obovata"], cids["Packera schweinitziana"]}), 1)
        self.assertEqual(len({cids["Packera paupercula"], cids["Packera plattensis"]}), 1)
        self.assertEqual(len({cids["Packera dubia"], cids["Packera aurea"], cids["Packera paupercula"]}), 3)
        self.assertIn("Packera anonyma", r["focal_neighbors"])
        self.assertTrue(bool(r["sub2_active"]) and bool(r["sub2_has_base_sp"]))
        self.assertTrue(bool(r["sub4_suppressed"]))

    @unittest.skipIf(SKIP_R_COND, "Fast audit or Rscript unavailable")
    def test_four_stream_genomic_affinity_eco_hygiene_and_fig07b_png(self) -> None:
        r = self.res
        self.assertFalse(bool(r["gr_none_avail"]))
        self.assertAlmostEqual(float(r["max_diff_no_geno"]), 0.0, places=9)
        self.assertTrue(bool(r["gr_yes_avail"]))
        self.assertEqual(int(r["gr_yes_n_taxa"]), 4)
        self.assertGreater(float(r["geno_aff_sisters"]), float(r["geno_aff_distant"]))
        self.assertGreater(float(r["w4_sisters"]), float(r["w4_distant"]))
        self.assertGreater(float(r["w4_w_geno"]), 0.0)
        banned = {
            "ssurgo_rock_outcrop_pct",
            "ssurgo_aws_0_50cm",
            "ssurgo_lithic_series_pct",
            "ultramafic_calcareous_gradient",
            "min_bedrock_depth_buffer_cm",
        }
        eco_set = set(r["eco_hyg_cols"])
        self.assertTrue(banned.isdisjoint(eco_set), f"Banned complex eco predictors leaked: {eco_set & banned}")
        self.assertTrue(
            {"lithic_outcrop_index", "wetland_seepage_index", "summer_moisture_deficit_idx", "maritime_continentality_idx"}.issubset(eco_set)
        )
        self.assertTrue(bool(r["png_magic_ok"]))
        self.assertEqual(int(r["n_pdfs"]), 0)


if __name__ == "__main__":
    unittest.main()
