if (grDevices::dev.cur() == 1L) grDevices::pdf(NULL)
# Step 07 Deme & Cluster Propagation Utilities
.safe_bool_prop <- function(x) if (is.logical(x)) !is.na(x) & x else if (is.numeric(x) || is.integer(x)) !is.na(x) & (x != 0) else tolower(trimws(as.character(x))) %in% c("1", "true", "t", "yes", "y")
.pick_chr_prop <- function(df, cols, default = "") {
  n <- nrow(df); if (n == 0L) return(character(0)); out <- rep(NA_character_, n)
  for (c in cols) if (c %in% names(df)) { v <- as.character(df[[c]]); f <- (is.na(out) | !nzchar(trimws(out))) & !is.na(v) & nzchar(trimws(v)); if (any(f)) out[f] <- v[f] }
  out[is.na(out) | !nzchar(trimws(out))] <- default; out
}
.pick_num_prop <- function(df, cols) {
  n <- nrow(df); if (n == 0L) return(numeric(0)); out <- rep(NA_real_, n)
  for (c in cols) if (c %in% names(df)) { v <- suppressWarnings(as.numeric(df[[c]])); f <- !is.finite(out) & is.finite(v); if (any(f)) out[f] <- v[f] }
  out
}
.pick_bool_prop <- function(df, cols, default = FALSE) { for (c in cols) if (c %in% names(df)) return(.safe_bool_prop(df[[c]])); rep(default, nrow(df)) }

.is_type_prop <- function(df) {
  if (nrow(df) == 0L) return(logical(0))
  ts <- .pick_chr_prop(df, c("typeStatus", "type_status"), "")
  .pick_bool_prop(df, "is_type_specimen", FALSE) | (.pick_chr_prop(df, c("authority_tier", "determiner_tier"), "") == "Tier_0_Type_Anchor") |
    (grepl("type", ts, ignore.case = TRUE) & !grepl("non[- ]?type|not a type", ts, ignore.case = TRUE))
}

.is_step07_quarantined_concept_prop <- function(df) {
  if (nrow(df) == 0L) return(logical(0))
  (trimws(.pick_chr_prop(df, "anchor_tier_source", "")) %in% c("Passive_Heterotypic_Type_Audit", "Passive_Type_Uncertain_Coord", "Passive_Open_or_Hybrid_Concept")) |
    (trimws(.pick_chr_prop(df, "concept_synonymy_status", "")) == "Heterotypic_Synonym") | .pick_bool_prop(df, "is_open_or_hybrid_concept", FALSE)
}

if (!exists(".is_step07_centroid_coord", mode = "function")) {
  .is_step07_centroid_coord <- function(df) {
    n <- nrow(df); if (n == 0L) return(logical(0)); if (!("is_centroid_or_inst_coord" %in% names(df))) return(rep(FALSE, n))
    .safe_bool_prop(df$is_centroid_or_inst_coord) & !(.is_type_prop(df) & !.is_step07_quarantined_concept_prop(df))
  }
}

.init_propagation_cols <- function(df, reset = FALSE) {
  n <- nrow(df)
  for (c in c("genomic_propagation_tier", "tier0_calibrated_taxon")) df[[c]] <- if (reset || !(c %in% names(df))) rep(NA_character_, n) else as.character(df[[c]])
  for (c in c("deme_mahal_d2", "geno_deme_weight_scale", "tier0_calibrated_prob", "cluster_propagated_gdi")) df[[c]] <- if (reset || !(c %in% names(df))) rep(NA_real_, n) else suppressWarnings(as.numeric(df[[c]]))
  df[["tier0_abstention_flag"]] <- if (reset || !("tier0_abstention_flag" %in% names(df))) rep(NA_integer_, n) else suppressWarnings(as.integer(df[["tier0_abstention_flag"]]))
  df
}

.resolve_prop_csv <- function(df, opts, attr_name, opt_names, filename) {
  for (on in opt_names) if (!is.null(opts[[on]]) && is.data.frame(opts[[on]])) return(opts[[on]])
  at <- attr(df, attr_name); if (!is.null(at) && is.data.frame(at)) return(at)
  od <- if (is.character(opts$output_dir) && length(opts$output_dir) == 1L && nzchar(opts$output_dir)) opts$output_dir else ""
  dirs <- unique(Filter(function(x) is.character(x) && length(x) == 1L && nzchar(x), c(opts$genomics_dir, opts$genomic_dir, opts$tables_dir, if (nzchar(od)) file.path(od, c("tables", "reports")) else NULL)))
  for (d in dirs) if (file.exists(p <- file.path(d, filename))) { res <- tryCatch(utils::read.csv(p, stringsAsFactors = FALSE), error = function(e) NULL); if (is.data.frame(res)) return(res) }
  NULL
}

.extract_collector_surname <- function(s) {
  vapply(as.character(s), function(x) {
    if (is.na(x) || !nzchar(trimws(x))) return("")
    x <- trimws(strsplit(x, "[;|&]|[[:space:]]+[Ee][Tt][[:space:]]+|[[:space:]]+[Aa][Nn][Dd][[:space:]]+")[[1L]][1L])
    if (grepl(",", x, fixed = TRUE)) part <- trimws(strsplit(x, ",", fixed = TRUE)[[1L]][1L])
    else { toks <- unlist(strsplit(x, "[[:space:]]+")); toks <- toks[!grepl("^[A-Za-z]\\.?$|^[A-Z]{1,3}\\.$", toks)]; part <- if (length(toks) > 0L) toks[length(toks)] else "" }
    tolower(gsub("[^a-zA-Z0-9]+", "", part))
  }, character(1L), USE.NAMES = FALSE)
}

.norm_record_number <- function(x) {
  vapply(as.character(x), function(s) {
    if (is.na(s)) return(""); raw <- tolower(trimws(s)); if (!nzchar(raw) || grepl("^0+$", raw)) return("")
    cln <- gsub("[^a-z0-9]+", "", raw); if (!nzchar(cln) || cln %in% c("sn", "unnumbered", "none", "na", "null", "unknown", "missing")) "" else cln
  }, character(1L), USE.NAMES = FALSE)
}

.haversine_km_vec <- function(lat1, lon1, lat2, lon2) {
  r <- pi / 180; a <- sin((lat2 - lat1) * r / 2)^2 + cos(lat1 * r) * cos(lat2 * r) * sin((lon2 - lon1) * r / 2)^2
  2.0 * 6371.0088 * asin(pmin(1.0, sqrt(pmax(0.0, a))))
}

.is_basal_mask_prop <- function(df) {
  if ("morph_tier_source" %in% names(df) && any(!is.na(df$morph_tier_source) & nzchar(as.character(df$morph_tier_source)))) return(!is.na(df$morph_tier_source) & as.character(df$morph_tier_source) == "Basal_Rosette")
  if ("basal_leaf_count" %in% names(df) && any(!is.na(blc <- suppressWarnings(as.integer(df$basal_leaf_count))))) return(!is.na(blc) & blc >= 1L)
  rep(TRUE, nrow(df))
}

.reg_cov_inv_prop <- function(X_sub, ridge = 1e-4) {
  p <- ncol(X_sub); if (nrow(X_sub) < 2L) return(diag(1.0, p))
  S <- stats::cov(X_sub); S[!is.finite(S)] <- 0.0; S <- 0.5 * (S + t(S))
  S_reg <- 0.90 * S + 0.10 * diag(pmax(diag(S), 1e-3), nrow = p) + diag(ridge, p)
  tryCatch(solve(S_reg), error = function(e) diag(1.0 / pmax(diag(S_reg), 1e-3), nrow = p))
}

.compute_morph_pool_prop <- function(df, pc_cols = NULL, ridge = 1e-4) {
  if (nrow(df) == 0L) return(NULL)
  if (is.null(pc_cols)) {
    pc_cols <- intersect(c(paste0("PC", 1:6), paste0("ECT_PC", 1:3)), names(df))
    if (length(pc_cols) < 2L) pc_cols <- intersect(c(paste0("morph_PC", 1:6), paste0("ECT_PC", 1:3)), names(df))
    pc_cols <- pc_cols[vapply(pc_cols, function(c) sum(is.finite(suppressWarnings(as.numeric(df[[c]])))) >= 2L, logical(1L))]
  }
  if (length(pc_cols) < 2L) return(NULL)
  X <- matrix(NA_real_, nrow = nrow(df), ncol = length(pc_cols), dimnames = list(NULL, pc_cols))
  for (k in seq_along(pc_cols)) X[, k] <- suppressWarnings(as.numeric(df[[pc_cols[k]]]))
  is_basal <- .is_basal_mask_prop(df); cc <- stats::complete.cases(X); idx <- which(is_basal & cc); if (length(idx) < max(3L, length(pc_cols) + 1L)) idx <- which(cc)
  list(X = X, S_inv = .reg_cov_inv_prop(X[idx, , drop = FALSE], ridge = ridge), pc_cols = pc_cols, is_basal = is_basal)
}

.subspace_mahal_d2_prop <- function(x_row, mu_row, X_pool, min_dims = 2L, S_inv_full = NULL, is_basal = NULL, ridge = 1e-4) {
  if (is.list(X_pool) && !is.null(X_pool$X)) { if (is.null(S_inv_full)) S_inv_full <- X_pool$S_inv; if (is.null(is_basal)) is_basal <- X_pool$is_basal; X_pool <- X_pool$X }
  if (is.null(x_row) || is.null(mu_row) || is.null(X_pool)) return(list(d2 = NA_real_, df = 0L))
  cols_ok <- which(is.finite(x_row) & is.finite(mu_row)); k <- length(cols_ok); if (k < min_dims) return(list(d2 = NA_real_, df = 0L))
  dv <- x_row[cols_ok] - mu_row[cols_ok]
  if (!is.null(S_inv_full) && k == ncol(X_pool) && is.matrix(S_inv_full) && nrow(S_inv_full) == k && ncol(S_inv_full) == k) {
    d2 <- as.numeric(crossprod(dv, S_inv_full %*% dv)); return(list(d2 = if (is.finite(d2)) d2 else NA_real_, df = k))
  }
  X_sub <- X_pool[, cols_ok, drop = FALSE]; cc <- stats::complete.cases(X_sub)
  idx <- if (!is.null(is_basal) && length(is_basal) == nrow(X_sub)) { ib <- which(is_basal & cc); if (length(ib) < max(3L, k + 1L)) which(cc) else ib } else which(cc)
  d2 <- as.numeric(crossprod(dv, .reg_cov_inv_prop(X_sub[idx, , drop = FALSE], ridge = ridge) %*% dv))
  list(d2 = if (is.finite(d2)) d2 else NA_real_, df = k)
}

.build_taxon_centroids <- function(df, morph_X, is_basal, t0) {
  quar <- .is_step07_quarantined_concept_prop(df); clean_type <- .is_type_prop(df) & !quar
  det <- .pick_chr_prop(df, c("authority_tier", "determiner_tier"), "Tier_3_Bronze")
  wt <- ifelse(t0$is_tier0_dna_anchor | clean_type, 4.0, ifelse(grepl("tier_?1|gold", det, ignore.case = TRUE), 2.0,
               ifelse(grepl("tier_?2|silver", det, ignore.case = TRUE), 0.75, ifelse(.pick_bool_prop(df, "is_promoted_bronze", FALSE), 0.50, 0.25))))
  rect_sp <- trimws(t0$rectified_taxon); pure_ev <- t0$evict_from_given_anchor & !t0$is_tier0_dna_anchor
  ok_base <- is_basal & !pure_ev & !quar & !.pick_bool_prop(df, "is_confident_label_noise", FALSE) &
    (!.pick_bool_prop(df, "is_passive", FALSE) | clean_type | t0$is_tier0_dna_anchor) &
    (!.pick_bool_prop(df, c("morphometric_outlier_flag", "is_morph_outlier"), FALSE) | clean_type | t0$is_tier0_dna_anchor)
  mu <- list(); p <- ncol(morph_X)
  for (tx in sort(unique(rect_sp[nzchar(rect_sp) & rect_sp != "Unknown"]))) {
    mu_v <- vapply(seq_len(p), function(j) {
      fin_j <- is.finite(morph_X[, j]); idx_j <- which(ok_base & rect_sp == tx & fin_j)
      if (length(idx_j) < 2L) idx_j <- which(is_basal & !pure_ev & !quar & rect_sp == tx & fin_j)
      if (length(idx_j) >= 2L) sum(morph_X[idx_j, j] * wt[idx_j]) / sum(wt[idx_j]) else NA_real_
    }, numeric(1L))
    if (sum(is.finite(mu_v)) >= 2L) mu[[tx]] <- mu_v
  }
  mu
}

.extract_known_ploidy_vec <- function(tbl) {
  n <- nrow(tbl); out <- rep(NA_integer_, n); if (n == 0L) return(out)
  for (col in c("voucher_ploidy_level", "geno_ploidy_level", "ploidy_level", "inferred_cytotype")) {
    if (!col %in% names(tbl)) next
    v <- suppressWarnings(as.integer(round(as.numeric(sub("x$", "", tolower(trimws(as.character(tbl[[col]]))))))))
    v[is.na(v) | v < 2L | v > 16L] <- NA_integer_; m <- is.na(out) & !is.na(v); if (any(m)) out[m] <- v[m]
  }
  out
}

.fit_gda_predict <- function(X_tr, y_tr, X_ev, classes, lambda = 0.20) {
  d <- ncol(X_tr); k_cnt <- length(classes)
  for (j in seq_len(d)) { v <- X_tr[is.finite(X_tr[, j]), j]; med <- if (length(v) > 0L) stats::median(v) else 0.0; X_tr[!is.finite(X_tr[, j]), j] <- med; X_ev[!is.finite(X_ev[, j]), j] <- med }
  mu_g <- colMeans(X_tr); sd_g <- apply(X_tr, 2L, stats::sd); sd_g[!is.finite(sd_g) | sd_g < 1e-6] <- 1.0
  Z_tr <- sweep(sweep(X_tr, 2L, mu_g, "-"), 2L, sd_g, "/"); Z_ev <- sweep(sweep(X_ev, 2L, mu_g, "-"), 2L, sd_g, "/")
  mu_k <- matrix(0.0, nrow = k_cnt, ncol = d); n_k <- integer(k_cnt); res_list <- vector("list", k_cnt)
  for (idx in seq_len(k_cnt)) { sub <- Z_tr[y_tr == classes[idx], , drop = FALSE]; n_k[idx] <- nrow(sub); mu_k[idx, ] <- colMeans(sub); res_list[[idx]] <- sweep(sub, 2L, mu_k[idx, ], "-") }
  S_w <- crossprod(do.call(rbind, res_list)) / max(1.0, nrow(Z_tr) - k_cnt); S_w[!is.finite(S_w)] <- 0.0; S_w <- 0.5 * (S_w + t(S_w))
  S_lw <- (1.0 - lambda) * S_w + lambda * diag(pmax(diag(S_w), 1e-4), nrow = d) + diag(1e-5, d)
  W_inv <- tryCatch(solve(S_lw), error = function(e) diag(1.0 / pmax(diag(S_lw), 1e-4), nrow = d))
  priors <- pmax(n_k / sum(n_k), 1.0 / (2.0 * k_cnt)); priors <- priors / sum(priors)
  logits <- matrix(0.0, nrow = nrow(Z_ev), ncol = k_cnt); d2_mat <- matrix(0.0, nrow = nrow(Z_ev), ncol = k_cnt)
  for (idx in seq_len(k_cnt)) {
    a_k <- W_inv %*% mu_k[idx, ]; logits[, idx] <- as.numeric(Z_ev %*% a_k - 0.5 * sum(mu_k[idx, ] * a_k) + log(priors[idx]))
    diff_k <- sweep(Z_ev, 2L, mu_k[idx, ], "-"); d2_mat[, idx] <- rowSums((diff_k %*% W_inv) * diff_k)
  }
  list(logits = logits, d2 = d2_mat)
}

.softmax_T <- function(logits, T_val = 1.0) { sc <- logits / max(T_val, 1e-4); mx <- apply(sc, 1L, max); ex <- exp(sweep(sc, 1L, mx, "-")); sweep(ex, 1L, rowSums(ex), "/") }

.extract_deme_anchors_prop <- function(df, ga_df = NULL) {
  extract_anc <- function(tbl, allow_g3, track_row) {
    tr <- .pick_chr_prop(tbl, c("geno_linkage_tier", "linkage_tier"), ""); qcf <- .pick_bool_prop(tbl, c("geno_quality_pass", "qc_pass"), FALSE)
    la <- .pick_num_prop(tbl, "decimalLatitude"); lo <- .pick_num_prop(tbl, "decimalLongitude"); ad <- .pick_bool_prop(tbl, "is_geno_admixed", FALSE)
    pr <- .pick_num_prop(tbl, "geno_posterior_prob"); hg <- trimws(.pick_chr_prop(tbl, "hybrid_generation_class", ""))
    pd <- trimws(.pick_chr_prop(tbl, "geno_predicted_taxon", "")); spv <- trimws(.pick_chr_prop(tbl, "species_standardized", ""))
    ok <- (tr %in% c("G1_Direct_Voucher", "G1_Exsiccatae_Isoduplicate", "G2_Population_Voucher_Proxy", if (allow_g3) "G3_Georeferenced_Unvouchered")) &
      qcf & !.pick_bool_prop(tbl, "is_centroid_or_inst_coord", FALSE) & !.is_step07_quarantined_concept_prop(tbl) &
      is.finite(la) & is.finite(lo) & !ad & ((is.finite(pr) & pr >= 0.80) | (hg == "Divergent_NonHybrid_Lineage"))
    tx <- ifelse(nzchar(pd) & pd != "Unknown", pd, ifelse(hg == "Divergent_NonHybrid_Lineage", spv, "")); ok <- ok & nzchar(tx) & (tx != "Unknown")
    doy_v <- .pick_num_prop(tbl, "hopkins_adjusted_doy"); rstat_v <- tolower(.pick_chr_prop(tbl, c("reproductive_trait_status", "capitulum_phenology_state"), ""))
    anth_v <- if ("is_active_anthesis" %in% names(tbl)) .safe_bool_prop(tbl$is_active_anthesis) else (is.finite(doy_v) & !grepl("vegetative|bud_only|post_fruiting", rstat_v))
    data.frame(
      lat = la[ok], lon = lo[ok], taxon = tx[ok], prob = ifelse(is.finite(pr[ok]), pr[ok], 0.85),
      azonal = .pick_num_prop(tbl, "is_azonal_lithology")[ok], ploidy = .extract_known_ploidy_vec(tbl)[ok],
      ph = .pick_num_prop(tbl, "fine_soil_ph")[ok], doy = ifelse(anth_v[ok] & is.finite(doy_v[ok]), doy_v[ok], NA_real_),
      cat = .pick_chr_prop(tbl, c("catalogNumber", "linked_catalog_number"), "")[ok], sid = .pick_chr_prop(tbl, "genetic_sample_id", "")[ok], gcl = trimws(.pick_chr_prop(tbl, "genetic_cluster", ""))[ok],
      df_row = if (track_row) which(ok) else rep(NA_integer_, sum(ok)), stringsAsFactors = FALSE
    )
  }
  anc_df <- extract_anc(df, FALSE, TRUE)
  if (!is.null(ga_df) && is.data.frame(ga_df) && nrow(ga_df) > 0L) {
    ga_rows <- extract_anc(ga_df, TRUE, FALSE)
    if (nrow(ga_rows) > 0L) {
      df_keys <- paste(anc_df$cat, anc_df$sid, sep = "::")[nzchar(anc_df$cat) | nzchar(anc_df$sid)]
      ga_keys <- paste(ga_rows$cat, ga_rows$sid, sep = "::")
      anc_df <- rbind(anc_df, ga_rows[!((nzchar(ga_rows$cat) | nzchar(ga_rows$sid)) & (ga_keys %in% df_keys)), , drop = FALSE])
    }
  }
  anc_df
}

.resolve_verified_sympatry_rows_prop <- function(df, t0, lat, lon, pure_ev) {
  quar <- .is_step07_quarantined_concept_prop(df); clean_type <- .is_type_prop(df) & !quar
  det <- .pick_chr_prop(df, c("authority_tier", "determiner_tier"), "")
  is_gold <- grepl("tier_?1|gold", det, ignore.case = TRUE) & !quar & !.pick_bool_prop(df, "is_confident_label_noise", FALSE) &
    !.pick_bool_prop(df, c("morphometric_outlier_flag", "is_morph_outlier"), FALSE) & !pure_ev
  which((clean_type | t0$is_tier0_dna_anchor | is_gold) & !quar & !.is_step07_centroid_coord(df) & is.finite(lat) & is.finite(lon) & nzchar(t0$rectified_taxon) & (t0$rectified_taxon != "Unknown") & !pure_ev)
}

.eval_deme_edaphic_gate_prop <- function(cand_sp, sp_i, ed_z_i, ph_i, ph_donor, ph_ref_cand) {
  d_ph <- if (is.finite(ph_i) && is.finite(ph_donor)) abs(ph_i - ph_donor) else Inf
  if (cand_sp == sp_i) {
    if (!is.finite(ed_z_i) || ed_z_i <= 1.75 || d_ph <= 0.60) return(list(pass = TRUE, mult = 1.0))
    if (ed_z_i <= 2.50 || (ed_z_i <= 3.0 && d_ph <= 1.00)) return(list(pass = TRUE, mult = 0.88))
    return(list(pass = FALSE, mult = 0.0))
  }
  if (!is.finite(ph_i) || is.null(ph_ref_cand)) return(list(pass = TRUE, mult = 1.0))
  z_ph <- abs(ph_i - ph_ref_cand[1L]) / ph_ref_cand[2L]
  if (z_ph > 2.5 && (!is.finite(ph_donor) || d_ph > 1.0)) return(list(pass = FALSE, mult = 0.0))
  list(pass = TRUE, mult = if (z_ph > 1.8) 0.88 else 1.0)
}

.eval_deme_anthesis_gate_prop <- function(is_anth_i, doy_i, doy_donor, mu_doy_cand, tau_eco_i) {
  if (!is_anth_i || !is.finite(doy_i)) return(list(pass = TRUE, mult = 1.0))
  doy_ref <- if (is.finite(doy_donor)) doy_donor else if (!is.null(mu_doy_cand) && is.finite(mu_doy_cand)) mu_doy_cand else NA_real_
  if (!is.finite(doy_ref)) return(list(pass = TRUE, mult = 1.0))
  dd <- abs(doy_i - doy_ref) %% 365.25; d_doy <- pmin(dd, 365.25 - dd)
  tau_d <- if (is.finite(tau_eco_i)) pmax(24.0, pmin(28.0, tau_eco_i)) else 26.0
  if (d_doy > tau_d) return(list(pass = FALSE, mult = 0.0))
  list(pass = TRUE, mult = if (d_doy > 18.0) 0.88 else 1.0)
}

.eval_mcl_cluster_gdi_prop <- function(i, mcl_key, mcl_col, src_ok, g_col, pred, rec_ploidy, lat, lon, is_cent = rep(FALSE, length(lat))) {
  if (!nzchar(mcl_key)) return(list(conflict = FALSE, gdi = NA_real_, pl = NULL))
  m_rows <- which(mcl_col == mcl_key)
  if (!is_cent[i] && is.finite(lat[i]) && is.finite(lon[i])) {
    has_c <- is.finite(lat[m_rows]) & is.finite(lon[m_rows]) & !is_cent[m_rows]; keep <- rep(TRUE, length(m_rows))
    if (any(has_c)) keep[has_c] <- .haversine_km_vec(lat[i], lon[i], lat[m_rows[has_c]], lon[m_rows[has_c]]) <= 250.0
    m_rows <- m_rows[keep]
  }
  if (length(m_rows) == 0L) return(list(conflict = FALSE, gdi = NA_real_, pl = NULL))
  mk_pls <- unique(rec_ploidy[m_rows[!is.na(rec_ploidy[m_rows])]])
  if (length(mk_pls) > 1L) return(list(conflict = TRUE, gdi = NA_real_, pl = NULL))
  c_pl <- if (length(mk_pls) == 1L) mk_pls[1L] else NULL
  if (!is.na(rec_ploidy[i]) && !is.null(c_pl) && rec_ploidy[i] != c_pl) return(list(conflict = TRUE, gdi = NA_real_, pl = c_pl))
  a_idx <- m_rows[src_ok[m_rows]]
  if (length(a_idx) >= 2L) {
    tb <- sort(table(pred[a_idx]), decreasing = TRUE)
    if (length(tb) > 0L && (as.numeric(tb[1L]) / length(a_idx)) >= 0.75) return(list(conflict = FALSE, gdi = stats::median(g_col[a_idx[pred[a_idx] == names(tb)[1L]]]), pl = c_pl))
  }
  list(conflict = FALSE, gdi = NA_real_, pl = c_pl)
}
