if (grDevices::dev.cur() == 1L) grDevices::pdf(NULL)
# ==============================================================================
# Module: scripts/morphometrics/ecoregion_manova_plot.R
# Multivariate MANOVA audit orchestrator, Pillai subtitle formatter, and Fig 04c hook
# ==============================================================================

format_manova_pillai_subtitle <- function(
  manova_res,
  factor_name = NULL,
  alt_factor = NULL,
  fallback_label = "",
  group_pattern = NULL
) {
  fb_out <- if (is.null(fallback_label)) "" else as.character(fallback_label)[1L]
  gt <- if (!is.null(manova_res) && is.list(manova_res) && !is.null(manova_res$group_table)) {
    manova_res$group_table
  } else if (is.data.frame(manova_res)) {
    manova_res
  } else {
    NULL
  }
  if (!is.null(gt) && nrow(gt) > 0L && "grouping_factor" %in% names(gt)) {
    idx <- integer(0)
    if (!is.null(factor_name) && nzchar(as.character(factor_name)[1L])) {
      idx <- which(gt$grouping_factor == as.character(factor_name)[1L])
    }
    if (length(idx) == 0L && !is.null(alt_factor) && nzchar(as.character(alt_factor)[1L])) {
      idx <- which(gt$grouping_factor == as.character(alt_factor)[1L])
    }
    if (length(idx) == 0L && !is.null(group_pattern) && nzchar(as.character(group_pattern)[1L])) {
      pat_g <- as.character(group_pattern)[1L]
      m_gf <- grepl(pat_g, as.character(gt$grouping_factor), ignore.case = TRUE)
      m_fl <- if ("factor_label" %in% names(gt)) grepl(pat_g, as.character(gt$factor_label), ignore.case = TRUE) else rep(FALSE, nrow(gt))
      idx <- which(m_gf | m_fl)
    }
    if (length(idx) > 0L) {
      r <- gt[idx[1], , drop = FALSE]
      v_val <- suppressWarnings(as.numeric(r$pillai_trace[1]))
      f_val <- suppressWarnings(as.numeric(r$approx_F[1]))
      ndf <- suppressWarnings(as.numeric(r$num_df[1]))
      ddf <- suppressWarnings(as.numeric(r$den_df[1]))
      p_val <- suppressWarnings(as.numeric(r$p_value[1]))
      if (is.finite(v_val) && is.finite(f_val)) {
        ndf_s <- if (is.finite(ndf) && abs(ndf - round(ndf)) < 1e-5) sprintf("%d", as.integer(round(ndf))) else sprintf("%.1f", ndf)
        ddf_s <- if (is.finite(ddf) && abs(ddf - round(ddf)) < 1e-5) sprintf("%d", as.integer(round(ddf))) else sprintf("%.1f", ddf)
        p_s <- if (!is.finite(p_val)) "NA" else if (p_val < 1e-4) sprintf("%.2e", p_val) else sprintf("%.4f", p_val)
        return(sprintf("Omnibus MANOVA: Pillai's Trace V = %.3f, F(%s, %s) = %.2f, p = %s", v_val, ndf_s, ddf_s, f_val, p_s))
      } else if (is.finite(v_val) && !is.finite(p_val)) {
        return(sprintf("Descriptive Multivariate Separation: Pillai's Trace V = %.3f (post-clustering p suppressed)", v_val))
      }
    }
  }
  fb_out
}

run_multivariate_manova_audit <- function(
  df,
  feature_cols,
  gmm_res = NULL,
  output_group_csv = "outputs/reports/manova_group_significance.csv",
  output_pairwise_csv = "outputs/reports/manova_pairwise_posthoc.csv",
  include_gmm_clusters = FALSE
) {
  run_dir <- trimws(Sys.getenv("PLANTMATICS_RUN_DIR", unset = ""))
  if (!nzchar(run_dir)) run_dir <- trimws(Sys.getenv("PACKERA_RUN_DIR", unset = ""))
  if (nzchar(run_dir)) {
    if (missing(output_group_csv) && !is.null(output_group_csv)) {
      output_group_csv <- file.path(run_dir, "outputs/reports/manova_group_significance.csv")
    }
    if (missing(output_pairwise_csv) && !is.null(output_pairwise_csv)) {
      output_pairwise_csv <- file.path(run_dir, "outputs/reports/manova_pairwise_posthoc.csv")
    }
  }

  empty_group_df <- data.frame(
    grouping_factor = character(0),
    factor_label = character(0),
    n_groups = integer(0),
    n_active_vouchers = integer(0),
    n_features = integer(0),
    active_features = character(0),
    pillai_trace = numeric(0),
    wilks_lambda = numeric(0),
    approx_F = numeric(0),
    num_df = numeric(0),
    den_df = numeric(0),
    partial_eta2 = numeric(0),
    p_value = numeric(0),
    q_value = numeric(0),
    stringsAsFactors = FALSE
  )
  empty_pairwise_df <- data.frame(
    grouping_factor = character(0),
    group_1 = character(0),
    group_2 = character(0),
    contrast = character(0),
    n_1 = integer(0),
    n_2 = integer(0),
    n_features = integer(0),
    hotelling_t2 = numeric(0),
    mahalanobis_d = numeric(0),
    pillai_trace = numeric(0),
    wilks_lambda = numeric(0),
    approx_F = numeric(0),
    num_df = numeric(0),
    den_df = numeric(0),
    partial_eta2 = numeric(0),
    p_value = numeric(0),
    q_value = numeric(0),
    stringsAsFactors = FALSE
  )

  atomic_write_manova_csv <- function(tbl, path) {
    if (!is.null(path) && nzchar(path)) {
      dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
      tmp_p <- paste0(path, ".", Sys.getpid(), ".tmp")
      utils::write.csv(tbl, file = tmp_p, row.names = FALSE)
      file.rename(tmp_p, path)
    }
  }

  if (is.null(df) || !is.data.frame(df) || nrow(df) < 4L) {
    atomic_write_manova_csv(empty_group_df, output_group_csv)
    atomic_write_manova_csv(empty_pairwise_df, output_pairwise_csv)
    return(list(
      group_table = empty_group_df,
      pairwise_table = empty_pairwise_df,
      output_group_csv = output_group_csv,
      output_pairwise_csv = output_pairwise_csv
    ))
  }

  n <- nrow(df)
  avail_cols <- .resolve_default_manova_features(df, feature_cols)
  avail_cols <- avail_cols[vapply(avail_cols, function(cn) {
    sum(is.finite(suppressWarnings(as.numeric(df[[cn]])))) >= 3L
  }, logical(1))]

  if (length(avail_cols) == 0L) {
    atomic_write_manova_csv(empty_group_df, output_group_csv)
    atomic_write_manova_csv(empty_pairwise_df, output_pairwise_csv)
    return(list(
      group_table = empty_group_df,
      pairwise_table = empty_pairwise_df,
      output_group_csv = output_group_csv,
      output_pairwise_csv = output_pairwise_csv
    ))
  }

  # Strict active/passive authority weighting (Tier_1_Gold & Tier_2_Silver active; Tier_3_Bronze = 0.0)
  is_type_row <- if ("is_type_specimen" %in% names(df)) as.logical(df$is_type_specimen) else rep(FALSE, n)
  is_type_row[is.na(is_type_row)] <- FALSE
  if ("authority_weight" %in% names(df)) {
    auth_w <- suppressWarnings(as.numeric(df$authority_weight))
    auth_w[is.na(auth_w) | !is.finite(auth_w)] <- 0.0
  } else if ("determiner_tier" %in% names(df)) {
    if (exists("compute_authority_precision_weights", mode = "function")) {
      auth_w <- compute_authority_precision_weights(df)
    } else {
      det_tier <- as.character(df$determiner_tier)
      err_col <- intersect(c("dawid_skene_error_rate", "label_error_rate", "error_rate_k"), names(df))
      prec_w <- if (length(err_col) > 0L) pmax(0.0, pmin(1.0, 1.0 - suppressWarnings(as.numeric(df[[err_col[1]]])))) else rep(1.0, n)
      prec_w[is.na(prec_w) | !is.finite(prec_w)] <- 1.0
      auth_w <- ifelse(det_tier %in% c("Tier_0_Type_Anchor", "Tier_1_Gold", "Tier_2_Silver") | is_type_row, prec_w, 0.0)
    }
  } else {
    auth_w <- rep(1.0, n)
  }
  is_pas_row <- if ("is_passive" %in% names(df)) { p_v <- as.logical(df$is_passive); !is.na(p_v) & p_v } else rep(FALSE, n)
  is_cln_row <- if ("is_confident_label_noise" %in% names(df)) { c_v <- as.logical(df$is_confident_label_noise); !is.na(c_v) & c_v } else rep(FALSE, n)
  auth_w[is_pas_row | is_cln_row] <- 0.0
  has_positive_anchors <- any(auth_w > 0)

  is_t1 <- if ("morph_tier_source" %in% names(df)) {
    !is.na(df$morph_tier_source) & (as.character(df$morph_tier_source) == "Basal_Rosette")

  } else {
    rep(TRUE, n)
  }

  pc_req <- grep("^PC[0-9]+$", avail_cols, value = TRUE)
  ect_req <- grep("^ECT_PC[0-9]+(_std)?$", avail_cols, value = TRUE)
  req_cols <- if (length(pc_req) > 0L) pc_req else if (length(ect_req) > 0L) ect_req else avail_cols
  has_pcs <- if (length(req_cols) > 0L) stats::complete.cases(df[, req_cols, drop = FALSE]) else rep(FALSE, n)

  # 1. species_standardized (>= 2 species)
  sp_raw <- if ("species_standardized" %in% names(df)) {
    as.character(df$species_standardized)
  } else if ("scientificName" %in% names(df)) {
    standardize_taxon_name(df$scientificName)
  } else if ("species_raw" %in% names(df)) {
    standardize_taxon_name(df$species_raw)
  } else {
    rep(NA_character_, n)
  }
  sp_vec <- trimws(as.character(sp_raw))
  sp_vec[is.na(sp_vec) | !nzchar(sp_vec) | tolower(sp_vec) == "unknown"] <- NA_character_

  # 2. Unsupervised Morphotype_K1..Kn from gmm_cluster / gmm_res (>= 2 clusters)
  cl_raw <- if (!is.null(gmm_res) && !is.null(gmm_res$classification) && length(gmm_res$classification) == n) {
    gmm_res$classification
  } else if ("gmm_cluster" %in% names(df)) {
    df$gmm_cluster
  } else {
    rep(NA_character_, n)
  }
  cl_chr <- trimws(as.character(cl_raw))
  cl_num <- suppressWarnings(as.integer(gsub("[^0-9]", "", cl_chr)))
  cl_vec <- ifelse(
    is.na(cl_chr) | !nzchar(cl_chr) | cl_chr == "NA",
    NA_character_,
    ifelse(grepl("^Morphotype_K[0-9]+$", cl_chr), cl_chr, ifelse(is.na(cl_num), NA_character_, paste0("Morphotype_K", cl_num)))
  )

  # 3. Non-NA regional_group US EPA Level III Ecoregions (>= 2 regions with >= 3 vouchers)
  quar_mask <- .is_quarantined_ecoregion_coord(df)
  rg_vec <- if ("regional_group" %in% names(df)) standardize_epa_l3_label_r(df$regional_group) else rep(NA_character_, n)
  if (all(is.na(rg_vec))) {
    lat_col <- intersect(c("decimalLatitude"), names(df))
    lon_col <- intersect(c("decimalLongitude"), names(df))
    if ((length(lat_col) > 0L && length(lon_col) > 0L) || any(c("stateProvince", "county") %in% names(df))) {
      lat_v <- if (length(lat_col) > 0L) suppressWarnings(as.numeric(df[[lat_col[1]]])) else rep(NA_real_, n)
      lon_v <- if (length(lon_col) > 0L) suppressWarnings(as.numeric(df[[lon_col[1]]])) else rep(NA_real_, n)
      if (any(quar_mask)) {
        lat_v[quar_mask] <- NA_real_
        lon_v[quar_mask] <- NA_real_
      }
      st_v <- if ("stateProvince" %in% names(df)) df$stateProvince else NULL
      co_v <- if ("county" %in% names(df)) df$county else NULL
      rg_vec <- tryCatch(infer_epa_level3_ecoregion(lat_v, lon_v, state_province = st_v, county = co_v), error = function(e) rg_vec)
    }
  }
  rg_vec[is.na(rg_vec) | !nzchar(trimws(rg_vec))] <- NA_character_

  sort_group_levels <- function(lvls) {
    if (length(lvls) > 0L && all(grepl("^Morphotype_K[0-9]+$", lvls))) {
      lvls[order(as.integer(sub("^Morphotype_K", "", lvls)))]
    } else {
      sort(lvls)
    }
  }

  factor_specs <- Filter(Negate(is.null), list(
    list(name = "species_standardized", label = "Taxonomic Species (species_standardized)", vec = sp_vec, min_n = 2L),
    if (isTRUE(include_gmm_clusters)) list(name = "gmm_cluster", label = "Unsupervised GMM Morphotype (Morphotype_K1..Kn)", vec = cl_vec, min_n = 2L) else NULL,
    list(name = "regional_group", label = "US EPA Level III Ecoregion (regional_group)", vec = rg_vec, min_n = 3L)
  ))

  group_rows <- list()
  pairwise_rows <- list()

  for (spec in factor_specs) {
    fac_res <- .eval_manova_factor_and_pairwise(
      spec = spec, df = df, avail_cols = avail_cols, n = n,
      has_positive_anchors = has_positive_anchors, auth_w = auth_w,
      is_t1 = is_t1, has_pcs = has_pcs, is_pas_row = is_pas_row,
      is_cln_row = is_cln_row, is_type_row = is_type_row,
      quar_mask = quar_mask, sort_group_levels = sort_group_levels
    )
    if (is.null(fac_res)) next
    if (!is.null(fac_res$group_row)) group_rows[[length(group_rows) + 1L]] <- fac_res$group_row
    if (!is.null(fac_res$pairwise_df)) pairwise_rows[[length(pairwise_rows) + 1L]] <- fac_res$pairwise_df
  }

  group_df <- if (length(group_rows) > 0L) do.call(rbind, group_rows) else empty_group_df
  if (nrow(group_df) > 0L) {
    fin_grp <- is.finite(group_df$p_value)
    if (any(fin_grp)) group_df$q_value[fin_grp] <- signif(stats::p.adjust(group_df$p_value[fin_grp], method = "BH"), 6)
    rownames(group_df) <- NULL
  }

  pairwise_df <- if (length(pairwise_rows) > 0L) do.call(rbind, pairwise_rows) else empty_pairwise_df
  if (nrow(pairwise_df) > 0L) {
    rownames(pairwise_df) <- NULL
  }

  atomic_write_manova_csv(group_df, output_group_csv)
  atomic_write_manova_csv(pairwise_df, output_pairwise_csv)

  list(
    group_table = group_df,
    pairwise_table = pairwise_df,
    output_group_csv = output_group_csv,
    output_pairwise_csv = output_pairwise_csv
  )
}

render_fig04c_supervised_ecoregion_cda <- function(
  audit_df,
  eco_cda_res = NULL,
  gmm_res = NULL,
  figures_dir = NULL,
  output_plot_path = NULL,
  manova_res = NULL,
  batch_csv = NULL,
  ecoregion_plot_path = NULL
) {
  fdir <- if (!is.null(figures_dir) && nzchar(figures_dir)) {
    as.character(figures_dir)
  } else if (exists("resolve_active_figures_dir", mode = "function")) {
    resolve_active_figures_dir(NULL, if (!is.null(output_plot_path)) output_plot_path else ecoregion_plot_path)
  } else {
    "outputs/figures"
  }
  if (!exists("render_step04_ecoregion_batch_png", mode = "function") && file.exists("scripts/figures/step04_figures.R")) {
    if (file.exists("scripts/figures/theme_packera_pub.R")) source("scripts/figures/theme_packera_pub.R")
    source("scripts/figures/step04_figures.R")
  }
  render_step04_ecoregion_batch_png(
    audit_df = audit_df,
    eco_cda_res = eco_cda_res,
    gmm_res = gmm_res,
    figures_dir = fdir,
    manova_res = manova_res,
    batch_csv = batch_csv,
    output_path = ecoregion_plot_path
  )
}
