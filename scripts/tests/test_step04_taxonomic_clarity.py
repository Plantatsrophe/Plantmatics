#!/usr/bin/env python3
"""Batched warm-R and static contract suite for Step 04 taxonomic clarity invariants."""

from __future__ import annotations

import json, math, re, sqlite3, sys, tempfile, unittest
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from scripts.core.db_schema import PLANTMATICS_SCHEMA_USER_VERSION, init_db_schema
from scripts.tests._r_session import rscript_available, run_warm_r_cmd

RUN_MI, RUN_P9 = "run_s04_clarity_mi", "run_s04_clarity_p9"
SKIP_R_REASON = "Skipping warm-R Step 04 taxonomic clarity tests under fast audit or missing Rscript."
SKIP_R_COND = not rscript_available()


def _is_finite(val: object) -> bool:
    if val is None or isinstance(val, bool):
        return False
    try:
        return math.isfinite(float(val))
    except (TypeError, ValueError):
        return False


def _seed_clarity_fixtures(tmp_p: Path) -> dict[str, Path]:
    mi_dir, intra_dir = tmp_p / "mi_run", tmp_p / "intra_run"
    mi_dir.mkdir(parents=True, exist_ok=True); intra_dir.mkdir(parents=True, exist_ok=True)
    db_path = tmp_p / "step04_clarity.sqlite"
    conn = sqlite3.connect(str(db_path)); init_db_schema(conn)
    conn.executemany("INSERT INTO pipeline_runs (run_id) VALUES (?);", [(RUN_MI,), (RUN_P9,)])
    v_rows, p_rows, rac_rows, om_rows = [], [], [], []
    v_lines = ["catalogNumber,scientificName,species_standardized,determiner_tier,type_status,typified_basionym,concept_synonymy_status,institutionCode,regional_group,decimalLatitude,decimalLongitude,is_primary_duplicate,pixels_per_mm,is_calibrated_mm,scale_source"]
    offs = [-0.42, -0.30, -0.18, -0.06, 0.06, 0.18, 0.30, 0.42]
    for idx in range(1, 17):
        is_a, j = idx <= 8, (idx - 1) % 8
        sp, cat, is_t = ("Packera aurea" if is_a else "Packera dubia"), f"MI_{idx:02d}", (1 if j == 0 else 0)
        ts, bas = ("Holotype" if is_t else ""), (("Senecio aureus" if is_a else "Cineraria dubia") if is_t else "")
        cs, inst, rg = ("Homotypic_Synonym" if is_t else "Accepted_Sensu_Stricto"), ("NCU" if idx % 2 else "DUKE"), ("Blue_Ridge" if j < 4 else "Piedmont")
        lat = round((35.50 + 0.002 * j) if (is_a and j < 4) else (35.10 + 0.24 * idx), 4)
        lon = round((-82.50 + 0.002 * j) if (is_a and j < 4) else (-83.20 + 0.24 * idx), 4)
        v_rows.append((cat, sp, sp, "Tier_1_Gold", ts, bas, cs, inst, rg, lat, lon, 1, 10.0, 1))
        p_rows.append((cat, 0)); rac_rows.append((RUN_MI, cat, 0, sp, "Tier_1_Gold", ts, is_t, 1.0, None))
        if idx <= 4:
            rac_rows.append((RUN_P9, cat, 0, sp, "Tier_1_Gold", ts, is_t, 1.0, None))
        e1, e2, e3 = round((-2.4 if is_a else 2.4) + offs[j], 4), round((0.35 if j % 2 else -0.35) + 0.4 * offs[j], 4), round((-0.25 if (j // 2) % 2 else 0.25) + 0.25 * offs[j], 4)
        om_rows.append((RUN_MI, cat, 0, sp, "Basal_Rosette", 10.0, 1, "ruler_calibrated", e1, e2, e3, round(abs(e1) * 0.4 + 0.2, 4), round((0.88 if is_a else 0.62) + 0.04 * offs[j], 4), round((1.15 if is_a else 1.85) - 0.12 * offs[j], 4), round(1.5 + 0.2 * offs[j], 4)))
        v_lines.append(f"{cat},{sp},{sp},Tier_1_Gold,{ts},{bas},{cs},{inst},{rg},{lat},{lon},1,10.0,1,ruler_calibrated")
    conn.executemany("INSERT INTO vouchers (catalogNumber,scientificName,species_standardized,determiner_tier,type_status,typified_basionym,concept_synonymy_status,institutionCode,regional_group,decimalLatitude,decimalLongitude,is_primary_duplicate,pixels_per_mm,is_calibrated_mm) VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?);", v_rows)
    conn.executemany("INSERT INTO plant_individuals (catalogNumber,plant_individual_id) VALUES (?,?);", p_rows)
    conn.executemany("INSERT INTO run_active_cohort (run_id,catalogNumber,plant_individual_id,species_standardized,determiner_tier,type_status,is_type_specimen,authority_weight,spatial_balance_weight) VALUES (?,?,?,?,?,?,?,?,?);", rac_rows)
    conn.executemany("INSERT INTO otu_morphometrics (run_id,catalogNumber,plant_individual_id,species_standardized,morph_tier_source,pixels_per_mm,is_calibrated_mm,scale_source,ECT_PC1,ECT_PC2,ECT_PC3,ect_lobe_complexity,solidity,dissection_index,aspect_ratio) VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?);", om_rows)
    conn.commit(); conn.close()
    v_mi_csv = mi_dir / "vouchers_mi.csv"; v_mi_csv.write_text("\n".join(v_lines) + "\n", encoding="utf-8")

    vi_lines = ["catalogNumber,scientificName,species_standardized,determiner_tier,institutionCode,regional_group,is_centroid_or_inst_coord,decimalLatitude,decimalLongitude,is_primary_duplicate,pixels_per_mm,is_calibrated_mm,scale_source"]
    hi_lines = ["catalogNumber,plant_individual_id,species_standardized,scientificName,determiner_tier,morph_tier_source,is_calibrated_mm,scale_source,solidity,dissection_index,aspect_ratio,A1,D1,A2,D2,A3,D3,A4,D4,PC1,PC2,PC3"]
    z12 = [-0.44, -0.36, -0.28, -0.20, -0.12, -0.04, 0.04, 0.12, 0.20, 0.28, 0.36, 0.44]
    for i in range(24):
        is_br, j = i < 12, i % 12
        cat, rg, inst = f"IN_{i + 1:02d}", ("Blue_Ridge" if is_br else "Piedmont"), ("NCU" if i % 2 else "DUKE")
        lat, lon = round((35.2 if is_br else 35.8) + 0.16 * j, 4), round((-83.0 if is_br else -79.2) + 0.16 * j, 4)
        p1, p2, p3 = round((-2.6 if is_br else 2.6) + z12[j], 4), round((0.40 if j % 2 else -0.40) + 0.3 * z12[j], 4), round((-0.30 if (j // 2) % 2 else 0.30) - 0.2 * z12[j], 4)
        vi_lines.append(f"{cat},Packera dubia,Packera dubia,Tier_1_Gold,{inst},{rg},FALSE,{lat},{lon},1,10.0,1,ruler_calibrated")
        hi_lines.append(f"{cat},0,Packera dubia,Packera dubia,Tier_1_Gold,Basal_Rosette,TRUE,ruler_calibrated,0.84,1.20,1.60,{round(0.4*p1,4)},{round(0.3*p1,4)},{round(0.4*p2,4)},{round(0.3*p2,4)},{round(0.4*p3,4)},{round(0.2*p3,4)},0.01,-0.01,{p1},{p2},{p3}")
    v_in_csv, h_in_csv = intra_dir / "vouchers_intra.csv", intra_dir / "harmonics_intra.csv"
    v_in_csv.write_text("\n".join(vi_lines) + "\n", encoding="utf-8"); h_in_csv.write_text("\n".join(hi_lines) + "\n", encoding="utf-8")
    return {"db_path": db_path, "mi_dir": mi_dir, "v_mi": v_mi_csv, "intra_dir": intra_dir, "v_in": v_in_csv, "h_in": h_in_csv}


_WARM_R_SCRIPT = r"""
if (grDevices::dev.cur() == 1L) grDevices::pdf(NULL)
args <- commandArgs(trailingOnly = TRUE)
out_json <- args[1L]; db_path <- args[2L]; mi_dir <- args[3L]; v_mi <- args[4L]; intra_dir <- args[5L]; v_in <- args[6L]; h_in <- args[7L]
Sys.setenv(PLANTMATICS_SKIP_FIGURES = "1", PACKERA_SKIP_FIGURES = "1")
suppressPackageStartupMessages({ library(jsonlite); library(mclust) })
for (.f in c("scripts/core/db_utils.R", "scripts/figures/theme_packera_pub.R", "scripts/morphometrics/dissected_leaf_utils.R",
             "scripts/morphometrics/pub_figures_utils.R", "scripts/morphometrics/authority_weighting_utils.R",
             "scripts/morphometrics/plasticity_mitigation_utils.R", "scripts/morphometrics/gmm_icl_bootstrap_utils.R",
             "scripts/morphometrics/pub_ecoregion_manova_utils.R", "scripts/morphometrics/step04_allometry_merger_utils.R",
             "scripts/morphometrics/step04_trait_projection_utils.R", "scripts/morphometrics/04_gmm_morphotools.R")) source(.f)
res <- list()
mk_opts <- function(d, v, h, intra = FALSE, ol = FALSE, ft = "") list(vouchers = v, harmonics = h, reproductive = "", outline_only = ol, correct_allometry = "false", allometry_report = file.path(d, "allom.csv"), allometry_pairwise = "", allometry_plot = "", use_clustvarsel = FALSE, output_flags = file.path(d, "flags.csv"), output_plot = "", output_report = file.path(d, "gmm.csv"), output_batch_audit = file.path(d, "batch.csv"), max_k = 3L, num_pcs = 3L, permutations = 49L, seed = 42L, intraspecific_mode = intra, focal_taxon = ft)

# 1. Type-I SS allometry false-positive guard & homogeneous PC1 sign-lock
n1 <- 32L; sp1 <- rep(c("Packera aurea", "Packera dubia"), each = 16L); zw1 <- rep(seq(-0.35, 0.35, length.out = 16L), 2L); is_a1 <- sp1 == "Packera aurea"
orth_w <- function(v) { o <- numeric(n1); for (s in unique(sp1)) { ii <- which(sp1 == s); o[ii] <- stats::resid(stats::lm(v[ii] ~ zw1[ii])) }; o }
set.seed(11L)
df_fp <- data.frame(catalogNumber = sprintf("FP_%02d", 1:n1), plant_individual_id = 0L, species_standardized = sp1, scientificName = sp1, determiner_tier = "Tier_1_Gold", morph_tier_source = "Basal_Rosette", solidity = 0.85, is_dissected = FALSE, is_calibrated_mm = TRUE, scale_source = "ruler_calibrated", decimalLatitude = seq(34.5, 39.2, length.out = n1), decimalLongitude = seq(-83.5, -76.8, length.out = n1), log_centsize = ifelse(is_a1, 2.0, 4.0) + zw1, A1 = ifelse(is_a1, -0.85, 0.85) + orth_w(rnorm(n1, 0, 0.04)), D1 = ifelse(is_a1, -0.60, 0.60) + orth_w(rnorm(n1, 0, 0.04)), A2 = ifelse(is_a1, 0.50, -0.50) + orth_w(rnorm(n1, 0, 0.04)), D2 = ifelse(is_a1, 0.40, -0.40) + orth_w(rnorm(n1, 0, 0.04)), A3 = orth_w(rnorm(n1, 0, 0.04)), D3 = orth_w(rnorm(n1, 0, 0.04)), A4 = orth_w(rnorm(n1, 0, 0.04)), D4 = orth_w(rnorm(n1, 0, 0.04)), stringsAsFactors = FALSE)
df_fp$centroid_size <- 10^df_fp$log_centsize; df_fp$PC1 <- df_fp$A1 + 0.5 * df_fp$D1; df_fp$PC2 <- df_fp$A2 - 0.5 * df_fp$D2
r_fp <- test_and_correct_allometry(df_fp, correct_mode = "auto", permutations = 199L, seed = 42L, r2_threshold = 0.10, output_report = "", output_plot = "", output_pairwise = "")
set.seed(12L); zw_h <- rep(seq(-1.2, 1.2, length.out = 16L), 2L); df_hp <- df_fp; df_hp$log_centsize <- 3.0 + zw_h; df_hp$centroid_size <- 10^df_hp$log_centsize
df_hp$A1 <- 0.65 * zw_h + ifelse(is_a1, -0.90, 0.90) + rnorm(n1, 0, 0.03); df_hp$D1 <- 0.45 * zw_h + ifelse(is_a1, -0.35, 0.35) + rnorm(n1, 0, 0.03)
df_hp$A2 <- 0.35 * zw_h + rnorm(n1, 0, 0.03); df_hp$D2 <- 0.25 * zw_h + rnorm(n1, 0, 0.03); df_hp$A3 <- ifelse(is_a1, 0.20, -0.20) + rnorm(n1, 0, 0.03)
pr_h <- stats::prcomp(as.matrix(df_hp[, c("A1", "D1", "A2", "D2", "A3", "D3", "A4", "D4")])); df_hp$PC1 <- pr_h$x[, 1L]; df_hp$PC2 <- pr_h$x[, 2L]; df_hn <- df_hp; df_hn$PC1 <- -pr_h$x[, 1L]
r_hp <- test_and_correct_allometry(df_hp, correct_mode = "auto", permutations = 199L, seed = 42L, r2_threshold = 0.10, output_report = "", output_plot = "", output_pairwise = "")
r_hn <- test_and_correct_allometry(df_hn, correct_mode = "auto", permutations = 199L, seed = 42L, r2_threshold = 0.10, output_report = "", output_plot = "", output_pairwise = "")
res$c1 <- list(fp_corr = isTRUE(r_fp$corrected), fp_a1_diff = max(abs(r_fp$df$A1 - df_fp$A1)), hp_corr = isTRUE(r_hp$corrected), hn_corr = isTRUE(r_hn$corrected), cor_pos = stats::cor(r_hp$df$PC1, df_hp$PC1), cor_neg = stats::cor(r_hn$df$PC1, df_hn$PC1))

# 2. Uncalibrated voucher protection & dimensionless active retention
z12 <- seq(-0.40, 0.40, length.out = 12L); sp2 <- rep(c("Packera aurea", "Packera dubia"), each = 12L); cal2 <- c(rep(TRUE, 6L), rep(FALSE, 6L), rep(TRUE, 6L), rep(FALSE, 6L))
df_cal <- data.frame(catalogNumber = sprintf("UC_%02d", 1:24), plant_individual_id = 0L, species_standardized = sp2, scientificName = sp2, determiner_tier = "Tier_1_Gold", morph_tier_source = "Basal_Rosette", authority_weight = 1.0, spatial_balance_weight = 1.0, is_calibrated_mm = cal2, scale_source = ifelse(cal2, "ruler_calibrated", rep(c("uncalibrated", "tick_lattice_unboxed"), length.out = 24L)), decimalLatitude = seq(34.5, 38.0, length.out = 24L), decimalLongitude = seq(-83.0, -77.5, length.out = 24L), petiole_length_mm = seq(15.0, 38.0, length.out = 24L), blade_length_mm = seq(25.0, 48.0, length.out = 24L), lamina_trace_length_mm = seq(35.0, 58.0, length.out = 24L), plant_height_cm = seq(22.0, 45.0, length.out = 24L), log_centsize = seq(2.4, 3.6, length.out = 24L), centroid_size = 10^seq(2.4, 3.6, length.out = 24L), area_mm2 = seq(220.0, 680.0, length.out = 24L), PC1 = ifelse(sp2 == "Packera aurea", -2.5, 2.5) + rep(z12, 2L), PC2 = rep(c(-0.3, 0.3), 12L) + 0.4 * rep(z12, 2L), PC3 = rep(c(0.25, 0.25, -0.25, -0.25), 6L), ECT_PC1 = ifelse(sp2 == "Packera aurea", -2.0, 2.0) + 0.5 * rep(z12, 2L), ECT_PC2 = rep(c(0.3, -0.3), 12L), ECT_PC3 = rep(c(-0.2, -0.2, 0.2, 0.2), 6L), solidity = ifelse(sp2 == "Packera aurea", 0.88, 0.64) + 0.03 * rep(z12, 2L), is_dissected = FALSE, basal_cordate_index = ifelse(sp2 == "Packera aurea", 0.35, 0.12) + 0.05 * rep(z12, 2L), AR_inv = ifelse(sp2 == "Packera aurea", 1.25, 1.65) + 0.08 * rep(z12, 2L), radiate_head_fraction = ifelse(sp2 == "Packera aurea", 0.90, 0.20) + 0.05 * rep(z12, 2L), blade_cie_l_q85_sheet_delta = ifelse(sp2 == "Packera aurea", -18.0, -7.0) + 1.5 * rep(z12, 2L), stringsAsFactors = FALSE)
df_msk <- mask_uncalibrated_step04_mm_cols(df_cal)
phys_c <- c("petiole_length_mm", "blade_length_mm", "lamina_trace_length_mm", "plant_height_cm", "log_centsize", "centroid_size", "area_mm2")
diml_c <- c("PC1", "PC2", "PC3", "ECT_PC1", "ECT_PC2", "ECT_PC3", "solidity", "basal_cordate_index", "AR_inv", "radiate_head_fraction", "blade_cie_l_q85_sheet_delta")
prep_u <- prepare_step04_standardized_features(df_msk, c("PC1", "PC2", "PC3"), outline_only = FALSE)
gmm_u <- run_gmm_cluster_analysis(prep_u$df, prep_u$gmm_cols, max_k = 3L, use_clustvarsel = FALSE, intraspecific_mode = FALSE, seed = 42L)
cda_u <- run_canonical_discriminant_analysis(prep_u$df, prep_u$cda_cols, c("Packera aurea", "Packera dubia"), gmm_res = gmm_u, intraspecific_mode = FALSE)
res$c2 <- list(uncal_phys_all_na = all(vapply(phys_c, function(cl) all(is.na(df_msk[[cl]][!cal2])), logical(1L))), cal_phys_all_fin = all(vapply(phys_c, function(cl) all(is.finite(df_msk[[cl]][cal2])), logical(1L))), diml_max_diff = max(vapply(diml_c, function(cl) max(abs(df_msk[[cl]] - df_cal[[cl]])), numeric(1L))), uncal_gmm_ok = all(!is.na(gmm_u$classification[!cal2])), uncal_cda_active = all(as.logical(cda_u$is_active[!cal2])), uncal_auth_pos = all(as.numeric(cda_u$authority_weight[!cal2]) > 0.0))

# 3. Step 02-03 trait integration, near-constant rejection & deterministic collinearity pruning
set.seed(33L); Q3 <- qr.Q(qr(matrix(rnorm(40L * 12L), 40L, 12L))); t_const <- 0.999 + rnorm(40L, 0, 1e-4); t_var <- rep(c(0.55, 0.85), each = 20L) + rnorm(40L, 0, 0.02)
df_f3 <- data.frame(catalogNumber = sprintf("F3_%02d", 1:40), plant_individual_id = 0L, species_standardized = rep(c("Packera aurea", "Packera dubia"), each = 20L), authority_weight = 1.0, spatial_balance_weight = 1.0, morph_tier_source = "Basal_Rosette", is_calibrated_mm = TRUE, scale_source = "ruler_calibrated", PC1 = 2.4 * Q3[, 1L], PC2 = 1.8 * Q3[, 2L], PC3 = 1.2 * Q3[, 3L], ECT_PC1 = 2.5 * Q3[, 4L], ECT_PC2 = 2.0 * Q3[, 5L], ECT_PC3 = 1.6 * Q3[, 6L], basal_cordate_index = 0.25 + 0.70 * Q3[, 7L], blade_cie_l_q85_sheet_delta = -14.0 + 40.0 * Q3[, 8L], radiate_head_fraction = pmin(1.0, pmax(0.0, 0.50 + 2.2 * Q3[, 9L])), terminal_lobe_area_frac = t_const, log_centsize = 3.0 + 2.0 * Q3[, 10L], petiole_length_mm = 25.0 + 40.0 * Q3[, 11L], plant_height_cm = 35.0 + 50.0 * Q3[, 12L], lamina_trace_length_mm = 45.0 + 30.0 * Q3[, 11L], stringsAsFactors = FALSE)
df_f3$terminal_lobe_area_frac_std <- as.numeric(scale(t_const)); df_fv <- data.frame(terminal_lobe_area_frac = t_var, terminal_lobe_area_frac_std = as.numeric(scale(t_var)))
flt_const <- filter_biologically_variable_dissection_scalars(df_f3, "terminal_lobe_area_frac_std"); flt_var <- filter_biologically_variable_dissection_scalars(df_fv, "terminal_lobe_area_frac_std")
df_col <- data.frame(PC1 = Q3[, 1L], PC2 = Q3[, 2L], ECT_PC1 = 3.0 * Q3[, 4L], ECT_PC1_std = as.numeric(scale(Q3[, 4L])) * 0.999, q25_solidity = 0.85 + 0.01 * Q3[, 4L], q25_solidity_std = as.numeric(scale(Q3[, 4L])) * 1.001, basal_cordate_index = 0.25 + 0.70 * Q3[, 7L], basal_cordate_index_std = as.numeric(scale(Q3[, 7L])) * 0.999, blade_cie_l_norm = 0.70 + 0.01 * Q3[, 7L], blade_cie_l_norm_std = as.numeric(scale(Q3[, 7L])) * 1.001)
cands_ab <- c("q25_solidity_std", "ECT_PC1_std", "blade_cie_l_norm_std", "basal_cordate_index_std")
pr_ab <- prune_collinear_dissection_scalars(df_col, base_cols = c("PC1", "PC2"), cand_scalars = cands_ab, group_vec = NULL, max_r = 0.75)
pr_ba <- prune_collinear_dissection_scalars(df_col, base_cols = c("PC1", "PC2"), cand_scalars = rev(cands_ab), group_vec = NULL, max_r = 0.75)
prep_f3 <- prepare_step04_standardized_features(df_f3, c("PC1", "PC2", "PC3"), outline_only = FALSE); prep_ol <- prepare_step04_standardized_features(df_f3, c("PC1", "PC2", "PC3"), outline_only = TRUE)
res$c3 <- list(const_rejected = !("terminal_lobe_area_frac_std" %in% flt_const), var_retained = ("terminal_lobe_area_frac_std" %in% flt_var), prune_ab = as.character(pr_ab), prune_ba = as.character(pr_ba), gmm_cols = as.character(prep_f3$gmm_cols), cda_cols = as.character(prep_f3$cda_cols), outline_cda_cols = as.character(prep_ol$cda_cols))

# 4. Multi-island shared-subspace accuracy (G>=2 & G==1) + all-multi-island SQLite ECT ingestion
df_mi <- df_cal[c(1:12, 13:24, 24L, 12L), , drop = FALSE]; rownames(df_mi) <- NULL; df_mi$catalogNumber <- sprintf("MI2_%02d", 1:26); df_mi$is_calibrated_mm <- TRUE; df_mi$scale_source <- "ruler_calibrated"
df_mi$type_status <- ""; df_mi$is_type_specimen <- FALSE; df_mi$type_status[c(1L, 13L)] <- "Holotype"; df_mi$is_type_specimen[c(1L, 13L)] <- TRUE
df_mi[25:26, c("PC1", "PC2", "PC3")] <- NA_real_; df_mi$ECT_PC1[25:26] <- c(2.45, -1.55); df_mi$ECT_PC2[25:26] <- c(-0.75, 0.75); df_mi$ECT_PC3[25:26] <- c(0.55, -0.55); df_mi$solidity[25:26] <- c(0.56, 0.80)
prep_mi <- prepare_step04_standardized_features(df_mi, c("PC1", "PC2", "PC3"), outline_only = FALSE)
gmm_mi <- run_gmm_cluster_analysis(prep_mi$df, prep_mi$gmm_cols, max_k = 3L, use_clustvarsel = FALSE, intraspecific_mode = FALSE, seed = 42L)
cda_mi <- run_canonical_discriminant_analysis(prep_mi$df, prep_mi$cda_cols, c("Packera aurea", "Packera dubia"), gmm_res = gmm_mi, intraspecific_mode = FALSE)
sh_mi <- project_shared_subspace_step04(prep_mi$df, cda_mi, gmm_mi, c("PC1", "PC2", "PC3"))
tc_mi <- compute_type_specimen_morph_control(prep_mi$df, sh_mi$cda_res$X_proj, sh_mi$cda_res$S_reg, sh_mi$cda_res, sh_mi$gmm_res)
df_g1 <- df_mi[c(13:24, 13:16, 25:26), , drop = FALSE]; rownames(df_g1) <- NULL; df_g1$catalogNumber <- sprintf("G1_%02d", 1:18); df_g1$species_standardized <- "Packera dubia"; df_g1$scientificName <- "Packera dubia"
df_g1$PC1[1:16] <- seq(-0.6, 0.6, length.out = 16L); df_g1$PC2[1:16] <- rep(c(-0.3, 0.3), 8L); df_g1$PC3[1:16] <- rep(c(-0.2, -0.2, 0.2, 0.2), 4L)
df_g1$ECT_PC1[1:16] <- seq(-0.5, 0.5, length.out = 16L); df_g1$ECT_PC2[1:16] <- rep(c(0.25, -0.25), 8L); df_g1$ECT_PC3[1:16] <- rep(c(-0.2, -0.2, 0.2, 0.2), 4L)
df_g1[17:18, c("PC1", "PC2", "PC3")] <- NA_real_; df_g1$ECT_PC1[17:18] <- c(1.6, -1.5); df_g1$ECT_PC2[17:18] <- c(-1.1, 1.2); df_g1$solidity[17:18] <- c(0.50, 0.78)
prep_g1 <- prepare_step04_standardized_features(df_g1, c("PC1", "PC2", "PC3"), outline_only = FALSE)
gmm_g1 <- list(classification = c(rep(1L, 16L), NA_integer_, NA_integer_), uncertainty = c(rep(0.0, 16L), NA_real_, NA_real_), best_k = 1L, selected_features = prep_g1$gmm_cols)
cda_g1 <- run_canonical_discriminant_analysis(prep_g1$df, prep_g1$cda_cols, "Packera dubia", gmm_res = gmm_g1, intraspecific_mode = TRUE)
sh_g1 <- project_shared_subspace_step04(prep_g1$df, cda_g1, gmm_g1, c("PC1", "PC2", "PC3"))
Sys.setenv(PLANTMATICS_DB_PATH = db_path, PACKERA_DB_PATH = db_path, RUN_ID = "run_s04_clarity_mi")
out_mi_sql <- run_gmm_morphotools_pipeline(mk_opts(mi_dir, v_mi, file.path(mi_dir, "missing_harmonics.csv")))
df_p9 <- out_mi_sql[order(out_mi_sql$catalogNumber), , drop = FALSE][1:4, , drop = FALSE]; df_p9$spatial_balance_weight <- c(0.375, 1.625, 0.8125, 1.1875)
ok_p9 <- persist_step04_to_sqlite(df_p9, vouchers_path = v_mi, run_id = "run_s04_clarity_p9")
Sys.unsetenv(c("PLANTMATICS_DB_PATH", "PACKERA_DB_PATH", "RUN_ID"))
res$c4 <- list(has_sh_obj = is.list(sh_mi$cda_res$shared_subspace) && length(sh_mi$cda_res$shared_subspace$shared_cols) >= 2L, mi_pc1_na = all(is.na(sh_mi$cda_res$X_proj[25:26, "PC1"])), mi_type_d2 = as.numeric(tc_mi$type_morph_mahal_d2[25:26]), mi_cong_d2 = as.numeric(tc_mi$min_congener_mahal_d2[25:26]), g1_pc1_na = all(is.na(sh_g1$cda_res$X_proj[17:18, "PC1"])), g1_self_d2 = as.numeric(sh_g1$cda_res$self_morph_mahal_d2[17:18]), g1_self_df = as.integer(sh_g1$cda_res$self_morph_df[17:18]), g1_grp = as.character(cda_g1$grouping_col[1L]), sql_n = nrow(out_mi_sql), sql_can1_fin = all(is.finite(as.numeric(out_mi_sql$can1))), sql_pred_ok = all(!is.na(out_mi_sql$cda_predicted_taxon)), ok_p9 = isTRUE(ok_p9))

# 5. Sequential ICL/bootstrap K=3->2 step-down, post-step-down Hennig ridgeline merge & Kish-N_eff BF direction
run_step3 <- function(Xm) {
  ic <- compute_gmm_icl_profile(Xm, data.frame(K = 1:3, BIC = c(0, 0, 0)), modelNames = "EEI", seed = 42L)
  gm <- ic$Classification_Entropy[2L] + ic$Classification_Entropy[3L] - 0.25; b <- c(-500.0, -500.0 + gm, -500.0 + gm + 0.50)
  set.seed(42L); f3 <- mclust::Mclust(Xm, G = 3L, modelNames = "EEI", verbose = FALSE); p3 <- stats::predict(f3, newdata = Xm)
  apply_gmm_icl_bootstrap_gate(Xm, Xm, data.frame(K = 1:3, BIC = b, BIC_EEI = b), best_k = 3L, best_name = "EEI", cls_sub = as.integer(p3$classification), unc_sub = as.numeric(p3$uncertainty), modelNames = "EEI", is_intra_gmm = TRUE, seed = 42L, z_eval_mat = as.matrix(p3$z))
}
set.seed(501L); X_sep <- cbind(c(rnorm(40L, -4.0, 0.7), rnorm(40L, 4.0, 0.7)), rnorm(80L, 0, 0.7)); gs_sep <- run_step3(X_sep)
set.seed(505L); zq90 <- stats::qnorm(stats::ppoints(90L)); ys180 <- rnorm(180L, 0, 1.0); ys180[1:90] <- ys180[1:90] - mean(ys180[1:90]); ys180[91:180] <- ys180[91:180] - mean(ys180[91:180])
X_sad <- cbind(c(-1.55 + zq90, 1.55 + zq90), ys180); gs_sad <- run_step3(X_sad)
mrg_sad <- merge_plastic_gmm_components_ridgeline(X_sad, X_sad, gs_sad$z_eval_mat, gs_sad$z_eval_mat, intraspecific_mode = TRUE)
bt_k <- data.frame(K = 1:3, BIC = c(-400.0, -350.0, -340.0), stringsAsFactors = FALSE)
bf_half <- .finalize_gmm_bic_table(bt_k, 3L, 3L, "VVI", 100.0, 0.50, 3L, NULL, FALSE, FALSE); bf_full <- .finalize_gmm_bic_table(bt_k, 3L, 3L, "VVI", 100.0, 1.00, 3L, NULL, FALSE, FALSE)
res$c5 <- list(sep_dicl2 = as.numeric(gs_sep$bic_table$Delta_ICL_vs_Null[2L]), sep_dicl3 = as.numeric(gs_sep$bic_table$Delta_ICL_vs_Null[3L]), sep_k = as.integer(gs_sep$best_k), sep_z_ncol = if (is.matrix(gs_sep$z_eval_mat)) ncol(gs_sep$z_eval_mat) else 0L, sad_step_k = as.integer(gs_sad$best_k), sad_z_ncol = if (is.matrix(gs_sad$z_eval_mat)) ncol(gs_sad$z_eval_mat) else 0L, sad_mrg_k = as.integer(mrg_sad$best_k), sad_mrg_any = isTRUE(mrg_sad$merged_any), bf_half = as.numeric(bf_half$BayesFactor_vs_Null[2:3]), bic_half = as.numeric(bf_half$Two_Delta_BIC[2:3]), bf_full = as.numeric(bf_full$BayesFactor_vs_Null[2:3]), bic_full = as.numeric(bf_full$Two_Delta_BIC[2:3]))

# 6. Intraspecific loocv_accuracy separation vs ecoregion_loocv_accuracy
out_in <- run_gmm_morphotools_pipeline(mk_opts(intra_dir, v_in, h_in, intra = TRUE, ol = TRUE, ft = "Packera dubia"))
res$c6 <- list(has_eco_col = "ecoregion_loocv_accuracy" %in% names(out_in), eco_loocv = if ("ecoregion_loocv_accuracy" %in% names(out_in)) as.numeric(out_in$ecoregion_loocv_accuracy[1L]) else NA_real_, intra_loocv_all_na = all(is.na(out_in$loocv_accuracy)))

# 7. Homotypic HOLO_EXTREME immunity vs extreme heterotypic type demotion & MCD outlier zero weight
set.seed(707L); n7 <- 25L
df_t7 <- data.frame(catalogNumber = c(sprintf("A_%02d", 1:9), "HOLO_EXTREME", "DUB_HOLO", sprintf("D_%02d", 1:12), "HET_TYPE_EXTREME", "MCD_GOLD_OUTLIER"), plant_individual_id = 0L, species_standardized = c(rep("Packera aurea", 10L), rep("Packera dubia", 15L)), scientificName = c(rep("Packera aurea", 10L), rep("Packera dubia", 15L)), determiner_tier = "Tier_1_Gold", morph_tier_source = "Basal_Rosette", type_status = c(rep("", 9L), "Holotype", "Holotype", rep("", 12L), "Holotype", ""), is_type_specimen = c(rep(FALSE, 9L), TRUE, TRUE, rep(FALSE, 12L), TRUE, FALSE), typified_basionym = c(rep("", 9L), "Senecio aureus", "Cineraria dubia", rep("", 12L), "Senecio alabamensis", ""), concept_synonymy_status = c(rep("Accepted_Sensu_Stricto", 9L), "Homotypic_Synonym", "Homotypic_Synonym", rep("Accepted_Sensu_Stricto", 12L), "Heterotypic_Synonym", "Accepted_Sensu_Stricto"), solidity = 0.85, is_dissected = FALSE, decimalLatitude = seq(34.5, 38.1, length.out = n7), decimalLongitude = seq(-83.0, -79.4, length.out = n7), PC1 = c(rnorm(9L, -3.0, 0.2), -12.0, 3.0, rnorm(12L, 3.0, 0.2), 9.5, 11.0), PC2 = c(rnorm(9L, 0, 0.2), 8.0, 0.0, rnorm(12L, 0, 0.2), 8.0, 9.0), PC3 = c(rnorm(9L, 0, 0.2), 6.0, 0.0, rnorm(12L, 0, 0.2), 6.0, -7.0), stringsAsFactors = FALSE)
anc_t7 <- select_adaptive_reference_anchors(df_t7, feature_matrix = as.matrix(df_t7[, c("PC1", "PC2", "PC3")]), taxon_col = "species_standardized")
cda_t7 <- run_canonical_discriminant_analysis(df_t7, c("PC1", "PC2", "PC3"), c("Packera aurea", "Packera dubia"), intraspecific_mode = FALSE)
df_t7_tc <- df_t7; df_t7_tc$anchor_tier_source <- cda_t7$anchor_tier_source
tc_t7 <- compute_type_specimen_morph_control(df_t7_tc, cda_t7$X_proj, cda_t7$S_reg, cda_t7, flags = rep(FALSE, n7), triage = rep("CLEAN", n7), reasons = rep("Verified", n7))
aud_t7 <- audit_misidentifications(df_t7, cda_t7)
res$c7 <- list(holo_anc_act = isTRUE(anc_t7$active_mask[10L]), holo_anc_tier = as.character(anc_t7$anchor_tier_source[10L]), holo_aw = as.numeric(cda_t7$authority_weight[10L]), het_anc_act = isTRUE(anc_t7$active_mask[24L]), het_anc_tier = as.character(anc_t7$anchor_tier_source[24L]), het_aw = as.numeric(cda_t7$authority_weight[24L]), het_tc_flag = isTRUE(tc_t7$flags[24L]), het_tc_triage = as.character(tc_t7$triage[24L]), het_tc_reason = as.character(tc_t7$reasons[24L]), het_aud_triage = as.character(aud_t7$triage_priority[24L]), het_aud_reason = as.character(aud_t7$discordance_reason[24L]), het_aud_misid = isTRUE(aud_t7$misidentification_flag[24L]), mcd_act = isTRUE(cda_t7$is_active[25L]), mcd_aw = as.numeric(cda_t7$authority_weight[25L]))

# 8. Open/hybrid concept audit overrides (exact symmetric centroids -> post_p = 0.50 at midpoint)
za <- c(0.0, -0.45, 0.45, -0.35, 0.35, -0.25, 0.25, -0.15, 0.15, 0.0); zb <- c(0.0, 0.30, 0.30, -0.40, -0.40, 0.20, 0.20, -0.10, -0.10, 0.0); zc <- c(0.0, 0.20, -0.20, -0.30, 0.30, 0.35, -0.35, -0.25, 0.25, 0.0)
df_oh <- data.frame(catalogNumber = c(sprintf("OA_%02d", 1:10), sprintf("OD_%02d", 1:10), "HYB_01", "OPEN_AFF_01"), plant_individual_id = 0L, species_standardized = c(rep("Packera aurea", 10L), rep("Packera dubia", 12L)), scientificName = c(rep("Packera aurea", 10L), rep("Packera dubia", 12L)), verbatim_taxon_concept = c(rep("Packera aurea", 10L), rep("Packera dubia", 10L), "Packera aurea x Packera dubia", "Packera aff. dubia"), concept_synonymy_status = c("Homotypic_Synonym", rep("Accepted_Sensu_Stricto", 9L), "Homotypic_Synonym", rep("Accepted_Sensu_Stricto", 9L), "Hybrid_Formula", "Open_Nomenclature"), is_open_or_hybrid_concept = c(rep(FALSE, 20L), TRUE, TRUE), determiner_tier = c(rep("Tier_1_Gold", 21L), "Tier_2_Silver"), type_status = c("Holotype", rep("", 9L), "Holotype", rep("", 11L)), is_type_specimen = c(TRUE, rep(FALSE, 9L), TRUE, rep(FALSE, 11L)), typified_basionym = c("Senecio aureus", rep("", 9L), "Cineraria dubia", rep("", 11L)), morph_tier_source = "Basal_Rosette", solidity = 0.85, is_dissected = FALSE, decimalLatitude = c(seq(34.5, 36.3, length.out = 10L), seq(37.0, 38.8, length.out = 10L), 36.6, 36.8), decimalLongitude = c(seq(-83.0, -81.2, length.out = 10L), seq(-80.5, -78.7, length.out = 10L), -80.9, -80.7), PC1 = c(-3.0 + za, 3.0 - za, 0.0, -3.0), PC2 = c(zb, -zb, 0.0, 0.0), PC3 = c(zc, -zc, 0.0, 0.0), stringsAsFactors = FALSE)
cda_oh <- run_canonical_discriminant_analysis(df_oh, c("PC1", "PC2", "PC3"), c("Packera aurea", "Packera dubia"), intraspecific_mode = FALSE)
aud_oh <- audit_misidentifications(df_oh, cda_oh); df_oh_ov <- df_oh; df_oh_ov$self_morph_mahal_d2 <- cda_oh$self_morph_mahal_d2; df_oh_ov$self_morph_df <- cda_oh$self_morph_df
ov_oh <- apply_open_or_hybrid_concept_audit_overrides(df_oh_ov, cda_oh, flags = rep(FALSE, 22L), triage = rep("CLEAN", 22L), reasons = rep("Verified", 22L), misid_flags = rep(FALSE, 22L), p_audit = 3L)
res$c8 <- list(hyb_post_p = as.numeric(cda_oh$posterior_prob[21L]), hyb_misid = isTRUE(aud_oh$misidentification_flag[21L]), hyb_ov_misid = isTRUE(ov_oh$misid_flags[21L]), hyb_triage = as.character(aud_oh$triage_priority[21L]), hyb_reason = as.character(aud_oh$discordance_reason[21L]), hyb_tier = as.character(aud_oh$anchor_tier_source[21L]), aff_post_p = as.numeric(cda_oh$posterior_prob[22L]), aff_self_d2 = as.numeric(cda_oh$self_morph_mahal_d2[22L]), aff_misid = isTRUE(aud_oh$misidentification_flag[22L]), aff_ov_flag = isTRUE(ov_oh$flags[22L]), aff_ov_misid = isTRUE(ov_oh$misid_flags[22L]), aff_triage = as.character(aud_oh$triage_priority[22L]), aff_reason = as.character(aud_oh$discordance_reason[22L]))

# 9. Trait label contract
res$c9 <- list(ar_inv_lbl = as.character(format_morph_trait_labels("AR_inv")[1L]), ar_inv_std_lbl = as.character(format_morph_trait_labels("AR_inv_std")[1L]))
jsonlite::write_json(res, out_json, auto_unbox = TRUE, na = "null", digits = 10)
"""


class TestStep04TaxonomicClarity(unittest.TestCase):
    """Batched warm-R and static contract suite for Step 04 taxonomic clarity invariants."""

    results: dict = {}

    @classmethod
    def setUpClass(cls) -> None:
        if SKIP_R_COND:
            cls.results = {}
            return
        cls._tmp_dir = tempfile.TemporaryDirectory(prefix="step04_clarity_")
        cls.tmp_p = Path(cls._tmp_dir.name)
        cls.fx = _seed_clarity_fixtures(cls.tmp_p)
        r_path, out_json = cls.tmp_p / "run_step04_clarity.R", cls.tmp_p / "clarity_results.json"
        r_path.write_text(_WARM_R_SCRIPT, encoding="utf-8")
        proc = run_warm_r_cmd(
            ["Rscript", str(r_path), str(out_json), str(cls.fx["db_path"]), str(cls.fx["mi_dir"]),
             str(cls.fx["v_mi"]), str(cls.fx["intra_dir"]), str(cls.fx["v_in"]), str(cls.fx["h_in"])],
            stage_key="step04", cwd=PROJECT_ROOT,
        )
        if proc.returncode != 0 or not out_json.exists():
            raise RuntimeError(f"Batched Step 04 clarity R failed (rc={proc.returncode}):\n{proc.stdout}\n{proc.stderr}")
        cls.results = json.loads(out_json.read_text(encoding="utf-8"))

    @classmethod
    def tearDownClass(cls) -> None:
        if hasattr(cls, "_tmp_dir"):
            cls._tmp_dir.cleanup()

    @unittest.skipIf(SKIP_R_COND, SKIP_R_REASON)
    def test_01_allometry_type1_ss_false_positive_guard_and_pc1_sign_lock(self) -> None:
        c1 = self.results["c1"]
        self.assertFalse(bool(c1["fp_corr"]))
        self.assertEqual(float(c1["fp_a1_diff"]), 0.0)
        self.assertTrue(bool(c1["hp_corr"]) and bool(c1["hn_corr"]))
        self.assertGreater(float(c1["cor_pos"]), 0.5)
        self.assertGreater(float(c1["cor_neg"]), 0.5)

    @unittest.skipIf(SKIP_R_COND, SKIP_R_REASON)
    def test_02_uncalibrated_voucher_mm_masking_and_dimensionless_retention(self) -> None:
        c2 = self.results["c2"]
        self.assertTrue(bool(c2["uncal_phys_all_na"]) and bool(c2["cal_phys_all_fin"]))
        self.assertEqual(float(c2["diml_max_diff"]), 0.0)
        self.assertTrue(bool(c2["uncal_gmm_ok"]) and bool(c2["uncal_cda_active"]) and bool(c2["uncal_auth_pos"]))

    @unittest.skipIf(SKIP_R_COND, SKIP_R_REASON)
    def test_03_step02_03_trait_integration_constant_rejection_and_collinearity(self) -> None:
        c3 = self.results["c3"]
        self.assertTrue(bool(c3["const_rejected"]) and bool(c3["var_retained"]))
        self.assertEqual(set(c3["prune_ab"]), {"ECT_PC1_std", "basal_cordate_index_std"})
        self.assertEqual(set(c3["prune_ba"]), {"ECT_PC1_std", "basal_cordate_index_std"})
        gmm_c, cda_c = list(c3["gmm_cols"]), list(c3["cda_cols"])
        self.assertTrue(all(re.match(r"^(PC[0-9]+|ECT_PC[0-9]+(_std)?)$", col) for col in gmm_c), gmm_c)
        self.assertTrue({"PC1", "PC2", "PC3"}.issubset(gmm_c) and any(re.match(r"^ECT_PC1(_std)?$", col) for col in gmm_c))
        for req in ("basal_cordate_index_std", "blade_cie_l_q85_sheet_delta_std", "radiate_head_fraction_std"):
            self.assertIn(req, cda_c)
        for excl in ("terminal_lobe_area_frac_std", "log_centsize_std", "petiole_length_mm_std", "plant_height_cm_std", "lamina_trace_length_mm_std"):
            self.assertNotIn(excl, cda_c); self.assertNotIn(excl, gmm_c)
        self.assertNotIn("basal_cordate_index_std", c3["outline_cda_cols"])

    @unittest.skipIf(SKIP_R_COND, SKIP_R_REASON)
    def test_04_multi_island_shared_subspace_and_all_ect_sqlite_ingestion(self) -> None:
        c4 = self.results["c4"]
        self.assertTrue(bool(c4["has_sh_obj"]) and bool(c4["mi_pc1_na"]))
        self.assertTrue(all(_is_finite(v) and float(v) > 0.0 for v in c4["mi_type_d2"] + c4["mi_cong_d2"]))
        self.assertGreater(float(c4["mi_cong_d2"][0]), float(c4["mi_type_d2"][0]))
        self.assertEqual(str(c4["g1_grp"]), "Morphotype_K1")
        self.assertTrue(bool(c4["g1_pc1_na"]) and all(_is_finite(v) and float(v) > 0.0 for v in c4["g1_self_d2"]) and all(int(d) >= 2 for d in c4["g1_self_df"]))
        self.assertEqual(int(c4["sql_n"]), 16)
        self.assertTrue(bool(c4["sql_can1_fin"]) and bool(c4["sql_pred_ok"]))

    @unittest.skipIf(SKIP_R_COND, SKIP_R_REASON)
    def test_05_sequential_icl_stepdown_post_ridgeline_and_kish_neff_bf(self) -> None:
        c5 = self.results["c5"]
        self.assertGreater(float(c5["sep_dicl2"]), 0.0)
        self.assertLessEqual(float(c5["sep_dicl3"]), 0.0)
        self.assertEqual((int(c5["sep_k"]), int(c5["sep_z_ncol"]), int(c5["sad_step_k"]), int(c5["sad_z_ncol"])), (2, 2, 2, 2))
        self.assertTrue(bool(c5["sad_mrg_any"]) and int(c5["sad_mrg_k"]) == 1)
        for bf_v, bic_v in zip(c5["bf_half"], c5["bic_half"]):
            self.assertTrue(float(bic_v) > 0.0 and float(bf_v) < float(bic_v))
        for bf_v, bic_v in zip(c5["bf_full"], c5["bic_full"]):
            self.assertAlmostEqual(float(bf_v), float(bic_v), places=4)

    @unittest.skipIf(SKIP_R_COND, SKIP_R_REASON)
    def test_06_intraspecific_ecoregion_loocv_separation(self) -> None:
        c6 = self.results["c6"]
        self.assertTrue(bool(c6["has_eco_col"]) and bool(c6["intra_loocv_all_na"]))
        self.assertTrue(_is_finite(c6["eco_loocv"]) and float(c6["eco_loocv"]) >= 0.85)

    @unittest.skipIf(SKIP_R_COND, SKIP_R_REASON)
    def test_07_homotypic_immunity_heterotypic_demotion_and_mcd_zero_weight(self) -> None:
        c7 = self.results["c7"]
        self.assertTrue(bool(c7["holo_anc_act"]) and str(c7["holo_anc_tier"]) == "Tier_0_Type_Anchor" and float(c7["holo_aw"]) > 0.0)
        self.assertTrue((not bool(c7["het_anc_act"])) and str(c7["het_anc_tier"]) == "Passive_Heterotypic_Type_Audit" and float(c7["het_aw"]) == 0.0)
        self.assertFalse(bool(c7["het_tc_flag"]) or bool(c7["het_aud_misid"]))
        self.assertEqual((str(c7["het_tc_triage"]), str(c7["het_aud_triage"])), ("HIGH", "HIGH"))
        self.assertTrue(str(c7["het_tc_reason"]).startswith("Heterotypic_Type_Divergence_") and str(c7["het_aud_reason"]).startswith("Heterotypic_Type_Divergence_"))
        self.assertTrue((not bool(c7["mcd_act"])) and float(c7["mcd_aw"]) == 0.0)

    @unittest.skipIf(SKIP_R_COND, SKIP_R_REASON)
    def test_08_open_and_hybrid_concept_audit_overrides(self) -> None:
        c8 = self.results["c8"]
        self.assertLess(float(c8["hyb_post_p"]), 0.80)
        self.assertFalse(bool(c8["hyb_misid"]) or bool(c8["hyb_ov_misid"]))
        self.assertEqual((str(c8["hyb_triage"]), str(c8["hyb_tier"])), ("HIGH", "Passive_Open_or_Hybrid_Concept"))
        self.assertTrue(str(c8["hyb_reason"]).startswith("Hybrid_Concept_Intermediate_"))
        self.assertTrue(float(c8["aff_post_p"]) >= 0.90 and float(c8["aff_self_d2"]) > 15.0)
        self.assertTrue(bool(c8["aff_misid"]) and bool(c8["aff_ov_flag"]) and bool(c8["aff_ov_misid"]) and str(c8["aff_triage"]) == "HIGH")
        self.assertTrue(str(c8["aff_reason"]).startswith("Open_Nomenclature_Reassignment_Candidate_to_"))

    @unittest.skipIf(SKIP_R_COND, SKIP_R_REASON)
    def test_09_trait_label_and_sqlite_spatial_balance_persistence(self) -> None:
        c9 = self.results["c9"]
        self.assertEqual((str(c9["ar_inv_lbl"]), str(c9["ar_inv_std_lbl"])), ("Involucre Aspect Ratio (H / W)", "Involucre Aspect Ratio (H / W)"))
        self.assertTrue(bool(self.results["c4"]["ok_p9"]))
        conn = sqlite3.connect(str(self.fx["db_path"]))
        uv = int(conn.execute("PRAGMA user_version;").fetchone()[0])
        mi_rows = conn.execute("SELECT catalogNumber, spatial_balance_weight FROM run_active_cohort WHERE run_id=? ORDER BY catalogNumber;", (RUN_MI,)).fetchall()
        p9_rows = conn.execute("SELECT catalogNumber, spatial_balance_weight FROM run_active_cohort WHERE run_id=? ORDER BY catalogNumber;", (RUN_P9,)).fetchall()
        conn.close()
        self.assertEqual((uv, uv), (PLANTMATICS_SCHEMA_USER_VERSION, 20261127))
        self.assertTrue(len(mi_rows) == 16 and all(w is not None and _is_finite(w) and float(w) > 0.0 for _, w in mi_rows))
        self.assertEqual([r[0] for r in p9_rows], ["MI_01", "MI_02", "MI_03", "MI_04"])
        for (_, got_w), exp_w in zip(p9_rows, (0.375, 1.625, 0.8125, 1.1875)):
            self.assertIsNotNone(got_w); self.assertAlmostEqual(float(got_w), exp_w, places=4)

    def test_10_static_step04_taxonomic_clarity_hygiene_and_exports(self) -> None:
        tp_path = PROJECT_ROOT / "scripts" / "morphometrics" / "step04_trait_projection_utils.R"
        s04_path = PROJECT_ROOT / "scripts" / "morphometrics" / "04_gmm_morphotools.R"
        self.assertLess(len(Path(__file__).resolve().read_text(encoding="utf-8").splitlines()), 380)
        self.assertTrue(tp_path.is_file(), f"Missing module: {tp_path}")
        tp_src, s04_src = tp_path.read_text(encoding="utf-8"), s04_path.read_text(encoding="utf-8")
        self.assertLess(len(tp_src.splitlines()), 500)
        for fn in ("mask_uncalibrated_step04_mm_cols", "prepare_step04_standardized_features", "project_shared_subspace_step04", "apply_open_or_hybrid_concept_audit_overrides"):
            self.assertIsNotNone(re.search(rf"^{re.escape(fn)}\s*<-\s*function\(", tp_src, flags=re.MULTILINE), f"Missing export {fn}")
            self.assertIn(fn, s04_src, f"04_gmm_morphotools.R does not wire {fn}")
        self.assertIn("ecoregion_loocv_accuracy", s04_src)
        self.assertIsNone(re.search(r"audited_df\$loocv_accuracy\s*<-\s*eco_cda_res", s04_src))
        for r_rel in ("scripts/figures/theme_packera_pub.R", "scripts/figures/step04_figures.R", "scripts/morphometrics/step04_trait_projection_utils.R"):
            code_lines = [ln for ln in (PROJECT_ROOT / r_rel).read_text(encoding="utf-8").splitlines() if not ln.lstrip().startswith("#")]
            self.assertIsNone(re.search(r"\\\\[sdw]", "\n".join(code_lines)), f"Non-POSIX regex escape in {r_rel}")


if __name__ == "__main__":
    unittest.main()
