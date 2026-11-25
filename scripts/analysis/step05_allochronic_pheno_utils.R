# ==============================================================================
# Step 05 Module: Allochronic Phenology, Anthesis Dominance & von Mises Mixture
# ==============================================================================
if (grDevices::dev.cur() == 1L) grDevices::pdf(NULL)

extract_doy_robust <- function(df) {
  doy_vec <- if ("doy" %in% names(df)) suppressWarnings(as.numeric(df$doy)) else rep(NA_real_, nrow(df))
  if ("eventDate" %in% names(df)) {
    ni <- which(is.na(doy_vec) & !is.na(df$eventDate))
    if (length(ni) > 0L) doy_vec[ni] <- as.integer(format(as.Date(df$eventDate[ni], format = "%Y-%m-%d"), "%j"))
  }
  if (any(is.na(doy_vec)) && all(c("year", "month", "day") %in% names(df))) {
    ni <- which(is.na(doy_vec) & !is.na(df$year) & !is.na(df$month) & !is.na(df$day))
    if (length(ni) > 0L) doy_vec[ni] <- as.integer(format(as.Date(paste(as.character(as.integer(df$year[ni])), sprintf("%02d", as.integer(df$month[ni])), sprintf("%02d", as.integer(df$day[ni])), sep = "-"), format = "%Y-%m-%d"), "%j"))
  }
  doy_vec
}

read_config_phenology_bounds <- function(config_path = "config/config.yaml") {
  min_doy <- 30.0; max_doy <- 220.0
  if (exists("load_step05_macroecology_config", mode = "function")) {
    cfg <- tryCatch(load_step05_macroecology_config(config_path), error = function(e) NULL)
    if (is.list(cfg)) {
      if (!is.null(cfg$min_flowering_doy) && is.finite(cfg$min_flowering_doy) && cfg$min_flowering_doy >= 1.0 && cfg$min_flowering_doy <= 366.0) min_doy <- as.numeric(cfg$min_flowering_doy)
      if (!is.null(cfg$max_flowering_doy) && is.finite(cfg$max_flowering_doy) && cfg$max_flowering_doy >= 1.0 && cfg$max_flowering_doy <= 366.0) max_doy <- as.numeric(cfg$max_flowering_doy)
      return(list(min_flowering_doy = min_doy, max_flowering_doy = max_doy))
    }
  }
  if (!is.null(config_path) && nzchar(config_path) && file.exists(config_path)) {
    lines <- tryCatch(readLines(config_path, warn = FALSE), error = function(e) character(0))
    for (pr in list(c("min_flowering_doy", "min"), c("max_flowering_doy", "max"))) {
      ml <- grep(paste0("^[[:space:]]*", pr[1L], "[[:space:]]*:"), lines, value = TRUE)
      if (length(ml) > 0L) {
        val <- suppressWarnings(as.numeric(sub("#.*$", "", sub("^[^:]+:[[:space:]]*", "", ml[1L]))))
        if (is.finite(val) && val >= 1.0 && val <= 366.0) { if (pr[2L] == "min") min_doy <- val else max_doy <- val }
      }
    }
  }
  list(min_flowering_doy = min_doy, max_flowering_doy = max_doy)
}

compute_circular_mean_doy <- function(doy_vec) {
  v <- suppressWarnings(as.numeric(doy_vec)); v <- v[is.finite(v) & v >= 1.0 & v <= 366.0]
  if (length(v) == 0L) return(NA_real_)
  th <- 2.0 * pi * (v / 365.25); sm <- mean(sin(th)); cm <- mean(cos(th))
  if (!is.finite(sm) || !is.finite(cm) || sqrt(sm^2 + cm^2) < 1e-8) return(as.numeric(stats::median(v)))
  ang <- atan2(sm, cm); if (ang <= 0.0) ang <- ang + 2.0 * pi
  c_doy <- ang * 365.25 / (2.0 * pi); as.numeric(if (c_doy <= 0.0) c_doy + 365.25 else c_doy)
}

circular_signed_doy_offset <- function(doy_vec, center_doy) {
  v <- suppressWarnings(as.numeric(doy_vec)); if (length(v) == 0L) return(numeric(0))
  c_v <- suppressWarnings(as.numeric(center_doy))
  if (length(c_v) == 0L || all(!is.finite(c_v))) return(rep(NA_real_, length(v)))
  if (length(c_v) == 1L) c_v <- rep(c_v, length(v)) else if (length(c_v) != length(v)) stop("circular_signed_doy_offset length mismatch")
  out <- rep(NA_real_, length(v)); ok <- is.finite(v) & is.finite(c_v)
  if (any(ok)) { d_th <- 2.0 * pi * ((v[ok] - c_v[ok]) / 365.25); out[ok] <- atan2(sin(d_th), cos(d_th)) * 365.25 / (2.0 * pi) }
  out
}

wrap_doy_365 <- function(doy_vec) {
  w <- suppressWarnings(as.numeric(doy_vec)) %% 365.25
  idx_le <- which(is.finite(w) & w <= 0.0); if (length(idx_le) > 0L) w[idx_le] <- w[idx_le] + 365.25; w
}

model_latitudinal_spring_baseline <- function(flowering_df, circ_mean_doy = NULL, lat_bounds = NULL, obs_weights = NULL) {
  message("Modeling Circular von Mises Latitudinal Anthesis Baseline (AICc Selection over Hopkins 1918 & Empirical Models)...")
  df_in <- as.data.frame(flowering_df, check.names = FALSE, stringsAsFactors = FALSE)
  df_in$decimalLatitude <- if ("decimalLatitude" %in% names(df_in) && any(!is.na(df_in$decimalLatitude))) suppressWarnings(as.numeric(df_in$decimalLatitude))
    else if ("latitude" %in% names(df_in)) suppressWarnings(as.numeric(df_in$latitude)) else rep(NA_real_, nrow(df_in))
  df_in$doy <- if ("doy" %in% names(df_in)) suppressWarnings(as.numeric(df_in$doy)) else rep(NA_real_, nrow(df_in))
  keep_v <- is.finite(df_in$decimalLatitude) & df_in$decimalLatitude >= -90.0 & df_in$decimalLatitude <= 90.0 &
            is.finite(df_in$doy) & df_in$doy >= 1.0 & df_in$doy <= 366.0
  if (!is.null(lat_bounds) && length(lat_bounds) == 2L && all(is.finite(lat_bounds))) keep_v <- keep_v & df_in$decimalLatitude >= min(lat_bounds) & df_in$decimalLatitude <= max(lat_bounds)

  valid <- df_in[keep_v, , drop = FALSE]
  w_valid <- if (!is.null(obs_weights) && length(obs_weights) == nrow(df_in)) suppressWarnings(as.numeric(obs_weights))[keep_v] else NULL
  day_per_rad <- 365.25 / (2.0 * pi)
  ref_circ_doy <- if (!is.null(circ_mean_doy) && length(circ_mean_doy) >= 1L && is.finite(circ_mean_doy[1L])) as.numeric(circ_mean_doy[1L]) else compute_circular_mean_doy(valid$doy)
  cf <- select_circular_hopkins_baseline_model(valid, obs_weights = w_valid)
  n_v <- nrow(valid); mu_0_doy <- wrap_doy_365(cf$mu_0 * day_per_rad)
  b_lat <- as.numeric(cf$b_lat); b_yr <- as.numeric(cf$b_year)
  min_lat <- if (n_v > 0L) min(valid$decimalLatitude, na.rm = TRUE) else NA_real_

  mk_base <- function(int_v, slp_v, s_yr, p_ml, h_exp, r2_v, pv, mod_obj, st_lbl, unif_fl) {
    list(
      intercept = int_v, slope = slp_v, slope_lat = slp_v, slope_year = s_yr,
      b_lon = cf$b_lon, ref_lon = cf$lon_ref, b_year = s_yr, ref_year = cf$year_ref,
      b_elev = cf$b_elev, ref_elev = cf$elev_ref, b_mat = cf$b_mat, ref_mat = cf$mat_ref,
      b_vernal = cf$b_vernal, vernal_ref = cf$vernal_ref, vernal_col = cf$vernal_col,
      min_lat = min_lat, pred_min_lat = p_ml, circ_mean_doy = ref_circ_doy, hopkins_expected = h_exp,
      r_squared = r2_v, p_value = pv, n = as.integer(n_v), n_eff = cf$n_eff, model = mod_obj, df = valid,
      r2 = r2_v, p_val = pv, status = st_lbl, selected_model = cf$selected_model,
      model_selection_table = cf$model_selection_table, aicc_weight = cf$aicc_weight,
      delta_aicc_vs_hopkins1918 = cf$delta_aicc_vs_hopkins1918, rayleigh_z = cf$rayleigh_z,
      rayleigh_p = cf$rayleigh_p, r_bar = cf$r_bar, is_uniform_null = unif_fl, predict_doy_rad = cf$predict_doy_rad
    )
  }

  if (n_v == 0L) {
    message("[WARN] Zero georeferenced anthesis records with valid DOY available for latitudinal baseline.")
    valid$expected_doy <- numeric(0); valid$hopkins_expected <- numeric(0); valid$delta_doy <- numeric(0)
    return(mk_base(ref_circ_doy, 0.0, 0.0, ref_circ_doy, numeric(0), 0.0, NA_real_, NULL, "insufficient_sample", TRUE))
  }

  theta_obs <- doy_to_rad(valid$doy); theta_hat <- cf$predict_doy_rad(valid)
  hopkins_expected <- wrap_doy_365(theta_hat * day_per_rad)
  delta_doy <- wrap_rad(theta_obs - theta_hat) * day_per_rad
  valid$expected_doy <- hopkins_expected; valid$hopkins_expected <- hopkins_expected; valid$delta_doy <- delta_doy
  pred_min_lat <- wrap_doy_365(cf$predict_doy_rad(valid[which.min(valid$decimalLatitude)[1L], , drop = FALSE])[1L] * day_per_rad)
  lin_intercept <- mu_0_doy - b_lat * cf$lat_ref
  if (n_v < 6L) {
    message(sprintf("[WARN] Small anthesis cohort (N=%d < 6); returning conservative circular fallback baseline.", n_v))
    return(mk_base(lin_intercept, b_lat, b_yr, pred_min_lat, hopkins_expected, 0.0, cf$rayleigh_p, cf, "insufficient_sample", cf$is_uniform_null))
  }

  r_null <- cf$r_bar; r_res <- mean(cos(wrap_rad(theta_obs - theta_hat)), na.rm = TRUE)
  r2_val <- if (is.finite(r_null) && is.finite(r_res) && (1.0 - r_null) > 1e-8) max(0.0, min(1.0, (r_res - r_null) / (1.0 - r_null))) else 0.0
  tbl <- cf$model_selection_table; s_i <- which(tbl$Model == cf$selected_model)[1L]; m0_i <- which(tbl$Model == "M0_Constant_Phase")[1L]
  k_sel <- if (!is.na(s_i)) tbl$k_params[s_i] else 2L
  ll_sel <- if (!is.na(s_i)) tbl$LogLik[s_i] else NA_real_; ll_m0 <- if (!is.na(m0_i)) tbl$LogLik[m0_i] else NA_real_

  p_val <- if (is.finite(k_sel) && k_sel > 2L && is.finite(ll_sel) && is.finite(ll_m0)) {
    stats::pchisq(max(0.0, 2.0 * (ll_sel - ll_m0)), df = k_sel - 2L, lower.tail = FALSE)
  } else if (identical(cf$selected_model, "M1_Hopkins_1918_Fixed")) {
    d_obs_c <- circular_signed_doy_offset(valid$doy, mu_0_doy); d_pred_c <- wrap_rad(theta_hat - cf$mu_0) * day_per_rad
    if (length(d_pred_c) >= 3L && is.finite(stats::sd(d_pred_c, na.rm = TRUE)) && stats::sd(d_pred_c, na.rm = TRUE) > 1e-8) {
      tryCatch({ sm <- summary(stats::lm(d_obs_c ~ d_pred_c)); if (nrow(sm$coefficients) >= 2L && ncol(sm$coefficients) >= 4L) as.numeric(sm$coefficients[2L, 4L]) else cf$rayleigh_p }, error = function(e) cf$rayleigh_p)
    } else cf$rayleigh_p
  } else cf$rayleigh_p

  message(sprintf("Circular Hopkins AICc Selected [%s]: slope=%.3f d/deg (b_elev=%.4f, b_mat=%.3f, b_vernal=%.3f, b_yr=%.4f, kappa=%.2f, AICc_w=%.3f, Rayleigh_p=%.3e, R2_circ=%.4f, N=%d)",
                  cf$selected_model, b_lat, cf$b_elev, cf$b_mat, cf$b_vernal, b_yr, cf$kappa, ifelse(is.finite(cf$aicc_weight), cf$aicc_weight, NA_real_), cf$rayleigh_p, r2_val, n_v))
  mk_base(lin_intercept, b_lat, b_yr, pred_min_lat, hopkins_expected, r2_val, p_val, cf, if (isTRUE(cf$is_uniform_null)) "uniform_null" else "ok", cf$is_uniform_null)
}

test_allochronic_bimodality_von_mises <- function(resid_days, weights = NULL) {
  d_raw <- suppressWarnings(as.numeric(resid_days)); ok <- is.finite(d_raw)
  if (!is.null(weights) && length(weights) == length(d_raw)) {
    w_raw <- suppressWarnings(as.numeric(weights)); ok <- ok & is.finite(w_raw) & (w_raw > 0.0); w <- w_raw[ok]
  } else w <- rep(1.0, sum(ok))
  d <- d_raw[ok]; n <- length(d)
  na_out <- list(allochronic_bimodality_delta_bic = NA_real_, ashman_d = NA_real_, mu1_days = NA_real_, mu2_days = NA_real_, pi1 = NA_real_, n = as.integer(n))
  if (n < 8L || sum(w) <= 1e-12) return(na_out)

  day_per_rad <- 365.25 / (2.0 * pi); th <- wrap_rad(2.0 * pi * (d / 365.25))
  w_sum <- sum(w); w_norm <- w / w_sum; n_eff <- (w_sum^2) / sum(w^2); w_eff <- w_norm * n_eff
  if (n_eff < 6.0) return(na_out)

  log_dvm <- function(th_v, mu_v, kap_v) {
    if (!is.finite(kap_v) || kap_v <= 1e-8) return(rep(-log(2.0 * pi), length(th_v)))
    kap_v * cos(th_v - mu_v) - log(2.0 * pi) - (log(max(1e-300, besselI(kap_v, nu = 0.0, expon.scaled = TRUE))) + kap_v)
  }
  c0 <- sum(w_norm * cos(th)); s0 <- sum(w_norm * sin(th))
  mu_k1 <- atan2(s0, c0); r_k1 <- min(1.0 - 1e-5, sqrt(c0^2 + s0^2))
  kap_k1 <- estimate_von_mises_kappa(r_k1, n_eff)
  bic_k1 <- -2.0 * sum(w_eff * log_dvm(th, mu_k1, kap_k1)) + 2.0 * log(n_eff)
  u <- wrap_rad(th - mu_k1)
  unim_out <- list(allochronic_bimodality_delta_bic = round(-3.0 * log(n_eff), 4), ashman_d = 0.0,
                   mu1_days = round(mu_k1 * day_per_rad, 3), mu2_days = round(mu_k1 * day_per_rad, 3), pi1 = 1.0, n = as.integer(n))
  if (stats::sd(u) <= 1e-6) return(unim_out)

  min_nw <- max(if (n_eff >= 16.0) 3.0 else 2.0, 0.08 * n_eff)
  run_em <- function(mu1, mu2, pi1, kap_i) {
    pi1 <- min(0.92, max(0.08, pi1)); k1 <- kap_i; k2 <- kap_i; ll_prev <- -Inf; r1 <- r_k1; r2 <- r_k1
    for (it in seq_len(80L)) {
      lp1 <- log(pi1) + log_dvm(th, mu1, k1); lp2 <- log(1.0 - pi1) + log_dvm(th, mu2, k2)
      lmx <- pmax(lp1, lp2); lse <- lmx + log(exp(lp1 - lmx) + exp(lp2 - lmx)); ll_cur <- sum(w_eff * lse)
      g1 <- exp(lp1 - lse); w1 <- w_eff * g1; w2 <- w_eff * (1.0 - g1); n1 <- sum(w1); n2 <- sum(w2)
      if (!is.finite(n1) || !is.finite(n2) || n1 < min_nw || n2 < min_nw) return(NULL)
      pi1 <- n1 / n_eff
      c1 <- sum(w1 * cos(th)) / n1; s1 <- sum(w1 * sin(th)) / n1; c2 <- sum(w2 * cos(th)) / n2; s2 <- sum(w2 * sin(th)) / n2
      mu1 <- atan2(s1, c1); mu2 <- atan2(s2, c2)
      r1 <- min(1.0 - 1e-4, max(1e-6, sqrt(c1^2 + s1^2))); r2 <- min(1.0 - 1e-4, max(1e-6, sqrt(c2^2 + s2^2)))
      k1 <- min(400.0, estimate_von_mises_kappa(r1, n1)); k2 <- min(400.0, estimate_von_mises_kappa(r2, n2))
      if (is.finite(ll_cur) && abs(ll_cur - ll_prev) < 1e-7) break
      ll_prev <- ll_cur
    }
    lp1 <- log(pi1) + log_dvm(th, mu1, k1); lp2 <- log(1.0 - pi1) + log_dvm(th, mu2, k2); lmx <- pmax(lp1, lp2)
    list(ll = sum(w_eff * (lmx + log(exp(lp1 - lmx) + exp(lp2 - lmx)))), pi1 = pi1, mu1 = mu1, mu2 = mu2, r1 = r1, r2 = r2)
  }

  qs <- stats::quantile(u, probs = c(0.15, 0.25, 0.30, 0.70, 0.75, 0.85), na.rm = TRUE, names = FALSE)
  c_a <- qs[1L]; c_b <- qs[6L]
  for (si in seq_len(6L)) { cl <- abs(u - c_a) <= abs(u - c_b); if (sum(cl) >= 2L && sum(!cl) >= 2L) { c_a <- mean(u[cl]); c_b <- mean(u[!cl]) } }
  k_st <- max(2.0, min(80.0, kap_k1 * 1.5))
  starts <- list(list(m1 = qs[2L], m2 = qs[5L], p1 = 0.50), list(m1 = c_a, m2 = c_b, p1 = 0.50),
                 list(m1 = qs[1L], m2 = qs[4L], p1 = 0.30), list(m1 = qs[3L], m2 = qs[6L], p1 = 0.70))
  best <- NULL
  for (st in starts) {
    fit <- run_em(wrap_rad(mu_k1 + st$m1), wrap_rad(mu_k1 + st$m2), st$p1, k_st)
    if (!is.null(fit) && is.finite(fit$ll) && (is.null(best) || fit$ll > best$ll)) best <- fit
  }
  if (is.null(best)) return(unim_out)

  delta_bic <- bic_k1 - (-2.0 * best$ll + 5.0 * log(n_eff))
  sd1 <- sqrt(max(1e-8, -2.0 * log(max(1e-8, min(1.0 - 1e-6, best$r1)))))
  sd2 <- sqrt(max(1e-8, -2.0 * log(max(1e-8, min(1.0 - 1e-6, best$r2)))))
  ash_d <- sqrt(2.0) * abs(wrap_rad(best$mu1 - best$mu2)) / sqrt(sd1^2 + sd2^2)
  m1_d <- wrap_rad(best$mu1) * day_per_rad; m2_d <- wrap_rad(best$mu2) * day_per_rad; p1_out <- best$pi1
  if (m1_d > m2_d) { tmp <- m1_d; m1_d <- m2_d; m2_d <- tmp; p1_out <- 1.0 - p1_out }
  list(allochronic_bimodality_delta_bic = round(delta_bic, 4), ashman_d = round(ash_d, 4),
       mu1_days = round(m1_d, 3), mu2_days = round(m2_d, 3), pi1 = round(p1_out, 4), n = as.integer(n))
}

compute_phenological_anomalies <- function(df, output_csv = NULL, summary_csv = NULL, demo = FALSE, config_path = "config/config.yaml") {
  message("Computing cohort-adaptive circular phenological anomalies (Delta DOY)...")
  n_r <- nrow(df); df$doy <- extract_doy_robust(df)
  lat_val <- if ("decimalLatitude" %in% names(df) && any(!is.na(df$decimalLatitude))) suppressWarnings(as.numeric(df$decimalLatitude))
    else if ("latitude" %in% names(df)) suppressWarnings(as.numeric(df$latitude)) else rep(NA_real_, n_r)
  df$decimalLatitude <- lat_val; df$latitude <- lat_val
  valid_doy <- is.finite(df$doy) & df$doy >= 1.0 & df$doy <= 366.0
  valid_lat <- is.finite(df$decimalLatitude) & df$decimalLatitude >= -90.0 & df$decimalLatitude <= 90.0

  if ("reproductiveCondition" %in% names(df)) {
    rc <- tolower(trimws(as.character(df$reproductiveCondition))); rc[is.na(rc)] <- ""
    dwc_fruit_only <- grepl("fruit", rc) & !grepl("fl", rc); dwc_sterile_only <- grepl("sterile", rc) & !grepl("fl", rc)
  } else { dwc_fruit_only <- rep(FALSE, n_r); dwc_sterile_only <- rep(FALSE, n_r) }
  pass_dwc <- !(dwc_fruit_only | dwc_sterile_only)

  v_state <- if ("voucher_phenological_state" %in% names(df)) { s <- tolower(trimws(as.character(df$voucher_phenological_state))); s[is.na(s)] <- ""; s } else rep("", n_r)
  .as_fin <- function(cn) { if (!(cn %in% names(df))) return(rep(NA_real_, n_r)); v <- suppressWarnings(as.numeric(df[[cn]])); v[!is.finite(v)] <- NA_real_; v }
  .coal_cnt <- function(a, b) { o <- a; m <- (!is.finite(o) | o <= 0.0) & is.finite(b) & b > 0.0; o[m] <- b[m]; m2 <- !is.finite(o) & is.finite(b); o[m2] <- b[m2]; o[!is.finite(o)] <- 0.0; o }
  .as_flag01 <- function(cn) {
    if (!(cn %in% names(df))) return(rep(0L, n_r)); raw <- df[[cn]]
    if (is.logical(raw)) ifelse(!is.na(raw) & raw, 1L, 0L) else if (is.numeric(raw) || is.integer(raw)) ifelse(is.finite(raw) & raw != 0, 1L, 0L)
    else ifelse(tolower(trimws(as.character(raw))) %in% c("1", "true", "t", "yes"), 1L, 0L)
  }

  tot_caps <- .coal_cnt(.coal_cnt(.as_fin("total_capitula"), .as_fin("capitula_count")), .as_fin("corymb_head_count_est"))
  n_anth_vec <- .coal_cnt(.as_fin("n_anthesis"), .as_fin("anthesis_head_count"))
  n_fruit_vec <- .coal_cnt(.as_fin("n_fruiting"), .as_fin("fruiting_head_count"))
  n_bud_vec <- .as_fin("n_bud"); n_bud_vec[!is.finite(n_bud_vec)] <- 0.0
  has_head_phase <- (n_anth_vec + n_fruit_vec) > 0.0
  w_anth_raw <- ifelse(has_head_phase, n_anth_vec / (n_anth_vec + n_fruit_vec + 0.5 * pmax(n_bud_vec, 0.0)), NA_real_)
  is_late_fruit_sen <- (has_head_phase & (w_anth_raw < 0.20)) | ((.as_flag01("Seed_Fruit_MANY") == 1L) & (.as_flag01("Flower_MANY") == 0L) & (n_anth_vec <= 0.0))

  has_valid_cv <- ((v_state != "" & v_state != "unknown" & v_state != "na") & tot_caps > 0) | (n_anth_vec > 0) | (n_fruit_vec > 0)
  is_cv_anth_label <- (grepl("anthesis", v_state) | (n_anth_vec > 0)) & !is_late_fruit_sen
  raw_anthesis_cand <- valid_doy & valid_lat & pass_dwc & !is_late_fruit_sen & ((has_valid_cv & is_cv_anth_label) | (!has_valid_cv))

  cfg_b <- if (exists("load_step05_macroecology_config", mode = "function")) load_step05_macroecology_config(config_path) else read_config_phenology_bounds(config_path)
  min_d <- if (!is.null(cfg_b$min_flowering_doy) && is.finite(cfg_b$min_flowering_doy)) as.numeric(cfg_b$min_flowering_doy) else 30.0
  max_d <- if (!is.null(cfg_b$max_flowering_doy) && is.finite(cfg_b$max_flowering_doy)) as.numeric(cfg_b$max_flowering_doy) else 220.0
  raw_doy_pool <- if (any(raw_anthesis_cand)) df$doy[raw_anthesis_cand] else df$doy[valid_doy]
  circ_mean_raw <- if (length(raw_doy_pool) > 0L) compute_circular_mean_doy(raw_doy_pool) else compute_circular_mean_doy(c(min_d, max_d))
  circ_med_doy <- if (length(raw_doy_pool) > 0L && is.finite(circ_mean_raw)) wrap_doy_365(circ_mean_raw + stats::median(circular_signed_doy_offset(raw_doy_pool, circ_mean_raw), na.rm = TRUE)) else circ_mean_raw
  is_vernal_cohort <- is.finite(circ_med_doy) && (min_d <= max_d) && (circ_med_doy >= min_d) && (circ_med_doy < 170.0)
  in_vernal_win <- valid_doy & (df$doy >= min_d) & (df$doy <= min(max_d, 240.0))
  offseason_vernal_mask <- if (is_vernal_cohort) (valid_doy & !in_vernal_win) else rep(FALSE, n_r)

  lat_cand <- if (is_vernal_cohort && any(raw_anthesis_cand & in_vernal_win)) (raw_anthesis_cand & in_vernal_win) else raw_anthesis_cand
  lat_pool <- if (any(lat_cand)) df$decimalLatitude[lat_cand] else df$decimalLatitude[valid_lat]
  lat_min_gate <- if (length(lat_pool) > 0L) max(-90.0, min(lat_pool, na.rm = TRUE) - 3.0) else -90.0
  lat_max_gate <- if (length(lat_pool) > 0L) min(90.0, max(lat_pool, na.rm = TRUE) + 3.0) else 90.0

  if (is_vernal_cohort) {
    pass_doy_gate <- in_vernal_win
    vp <- if (any(raw_anthesis_cand & in_vernal_win)) df$doy[raw_anthesis_cand & in_vernal_win] else raw_doy_pool
    circ_mean_doy <- if (length(vp) > 0L) compute_circular_mean_doy(vp) else circ_mean_raw
  } else if (sum(raw_anthesis_cand) >= 5L) {
    circ_mean_doy <- circ_mean_raw; co <- circular_signed_doy_offset(df$doy[raw_anthesis_cand], circ_mean_doy)
    mad_d <- stats::mad(co, constant = 1.4826, na.rm = TRUE)
    h_arc <- min(120.0, max(45.0, 3.5 * (if (is.finite(mad_d) && mad_d > 0) mad_d else 21.0)))
    ao <- circular_signed_doy_offset(df$doy, circ_mean_doy)
    pass_doy_gate <- valid_doy & is.finite(ao) & (ao >= max(-h_arc, as.numeric(stats::quantile(co, 0.005, na.rm = TRUE)) - 21.0)) & (ao <= min(h_arc, as.numeric(stats::quantile(co, 0.995, na.rm = TRUE)) + 21.0))
  } else {
    pass_doy_gate <- if (min_d <= max_d) (valid_doy & df$doy >= min_d & df$doy <= max_d) else (valid_doy & (df$doy >= min_d | df$doy <= max_d))
    circ_mean_doy <- circ_mean_raw
  }
  pass_calendar <- pass_doy_gate & valid_lat & (df$decimalLatitude >= lat_min_gate) & (df$decimalLatitude <= lat_max_gate)
  cv_anthesis <- (has_valid_cv & is_cv_anth_label) | (!has_valid_cv & pass_calendar & pass_dwc & !is_late_fruit_sen)
  cv_fruit <- (has_valid_cv & !is_cv_anth_label & grepl("fruit", v_state)) | (!has_valid_cv & dwc_fruit_only)
  cv_sterile <- (has_valid_cv & !is_cv_anth_label & (v_state %in% c("sterile", "bud", "senescent") | grepl("sterile|bud", v_state))) | (!has_valid_cv & dwc_sterile_only)

  is_verified_anthesis <- pass_calendar & pass_dwc & cv_anthesis & !offseason_vernal_mask & !is_late_fruit_sen
  is_fruit_excluded <- (pass_calendar & (cv_fruit | dwc_fruit_only) & !is_verified_anthesis) | is_late_fruit_sen
  is_sterile_excluded <- pass_calendar & (cv_sterile | dwc_sterile_only) & !is_verified_anthesis & !is_fruit_excluded
  df$is_verified_anthesis <- is_verified_anthesis; df$is_active_anthesis <- is_verified_anthesis
  df$anthesis_dominance_weight <- ifelse(has_head_phase, round(w_anth_raw, 6), ifelse(is_verified_anthesis, 1.0, 0.0))
  phen_state <- ifelse(nzchar(v_state) & !(v_state %in% c("unknown", "na")), v_state, "unknown")
  phen_state[is_verified_anthesis & phen_state %in% c("unknown", "")] <- "anthesis"
  phen_state[is_sterile_excluded] <- "sterile"; phen_state[is_fruit_excluded] <- "fruiting_only"
  phen_state[!pass_calendar | offseason_vernal_mask] <- "offseason_or_fruiting_flush"
  phen_state[is_late_fruit_sen] <- "fruiting_only"; df$phenological_state <- phen_state
  afs <- rep("Excluded_Calendar_Gate", n_r); afs[pass_calendar] <- "Excluded_Other"; afs[is_sterile_excluded] <- "Excluded_Sterile"
  afs[is_fruit_excluded] <- "Excluded_Fruit_Only"; afs[is_verified_anthesis] <- "Verified_Anthesis"
  afs[is_late_fruit_sen] <- "Late_Fruiting_Senescent_Excluded"; df$anthesis_filter_status <- afs; df$phenophase_status <- afs

  if (!("authority_weight" %in% names(df)) && exists("resolve_ecological_authority_weights", mode = "function")) df$authority_weight <- resolve_ecological_authority_weights(df)
  if (!("spatial_balance_weight" %in% names(df)) && exists("compute_spatial_balance_weights", mode = "function")) df$spatial_balance_weight <- compute_spatial_balance_weights(df)
  is_active_ref <- rep(TRUE, n_r); is_t_ph <- if ("is_type_specimen" %in% names(df)) (.as_flag01("is_type_specimen") == 1L) else rep(FALSE, n_r)
  if ("determiner_tier" %in% names(df)) {
    dt_ph <- as.character(df$determiner_tier); dt_ph[is.na(dt_ph)] <- "Tier_3_Bronze"
    cand_act <- grepl("Tier_1|Gold|Tier_2|Silver", dt_ph, ignore.case = TRUE) | (.as_flag01("is_promoted_bronze") == 1L) | is_t_ph
    if (sum(cand_act & df$is_verified_anthesis) >= 5L || any(cand_act)) is_active_ref <- cand_act
  }
  if ("is_confident_label_noise" %in% names(df)) is_active_ref <- is_active_ref & ((.as_flag01("is_confident_label_noise") == 0L) | is_t_ph)
  quar_det_ph <- if (exists("detect_r_centroid_or_inst_quarantine", mode = "function")) as.logical(detect_r_centroid_or_inst_quarantine(df, is_type_vec = rep(FALSE, n_r))$is_quarantined) else rep(FALSE, n_r)
  quar_det_ph[is.na(quar_det_ph)] <- FALSE; quar_ph_mask <- (.as_flag01("is_centroid_or_inst_coord") == 1L) | quar_det_ph
  open_hyb_ph <- if (exists(".identify_open_or_hybrid_concept_mask", mode = "function")) as.logical(.identify_open_or_hybrid_concept_mask(df)) else (.as_flag01("is_open_or_hybrid_concept") == 1L)
  open_hyb_ph[is.na(open_hyb_ph)] <- FALSE
  ph_grp_v <- if ("species_standardized" %in% names(df)) as.character(df$species_standardized) else rep("ALL", n_r)
  ext_ph_mask <- if (exists(".is_step05_extralimital_non_type", mode = "function")) .is_step05_extralimital_non_type(df, group_vec = ph_grp_v, clean_base_mask = is_active_ref & !quar_ph_mask, is_type_vec = is_t_ph) else rep(FALSE, n_r)
  excl_base <- quar_ph_mask | open_hyb_ph; is_active_ref <- is_active_ref & !excl_base & !ext_ph_mask

  base_idx <- which(df$is_verified_anthesis & is_active_ref)
  if (length(base_idx) < 5L && sum(df$is_verified_anthesis & !excl_base) > length(base_idx)) base_idx <- which(df$is_verified_anthesis & !excl_base)
  if (length(base_idx) < 2L) base_idx <- which(pass_calendar & !excl_base & !is_late_fruit_sen)
  if (length(base_idx) == 0L) base_idx <- which(valid_doy & valid_lat & !excl_base & !is_late_fruit_sen)
  calc_row_wts <- function(idx_v) {
    if (length(idx_v) == 0L) return(numeric(0))
    wa <- if ("authority_weight" %in% names(df)) suppressWarnings(as.numeric(df$authority_weight[idx_v])) else rep(1.0, length(idx_v))
    ws <- if ("spatial_balance_weight" %in% names(df)) suppressWarnings(as.numeric(df$spatial_balance_weight[idx_v])) else rep(1.0, length(idx_v))
    wd <- suppressWarnings(as.numeric(df$anthesis_dominance_weight[idx_v]))
    wa[!is.finite(wa) | wa < 0.0] <- 0.0; ws[!is.finite(ws) | ws <= 0.0] <- 1.0; wd[!is.finite(wd) | wd < 0.0] <- 0.0
    if (all(wa <= 0.0)) wa <- rep(1.0, length(idx_v)); if (all(wd <= 0.0)) wd <- rep(1.0, length(idx_v))
    wt <- wa * ws * wd; if (sum(wt) <= 1e-12) rep(1.0, length(idx_v)) else wt
  }
  baseline <- model_latitudinal_spring_baseline(df[base_idx, , drop = FALSE], circ_mean_doy = circ_mean_doy, lat_bounds = c(lat_min_gate, lat_max_gate), obs_weights = calc_row_wts(base_idx))

  day_per_rad <- 365.25 / (2.0 * pi); theta_pred <- rep(NA_real_, n_r)
  df$expected_doy <- NA_real_; df$hopkins_expected <- NA_real_; df$is_pheno_uniform_null <- isTRUE(baseline$is_uniform_null)
  if (any(valid_lat) && is.function(baseline$predict_doy_rad)) {
    tp <- baseline$predict_doy_rad(df[valid_lat, , drop = FALSE]); theta_pred[valid_lat] <- tp
    ed <- ifelse(is.finite(tp), wrap_doy_365(tp * day_per_rad), NA_real_); df$expected_doy[valid_lat] <- ed; df$hopkins_expected[valid_lat] <- ed
  }
  df$delta_doy <- NA_real_; df$delta_doy_resid <- NA_real_; df$hopkins_adjusted_doy_resid <- NA_real_
  has_both <- valid_doy & df$is_verified_anthesis & !offseason_vernal_mask & is.finite(theta_pred) & is.finite(df$expected_doy)
  if (any(has_both)) { dd_b <- wrap_rad(doy_to_rad(df$doy[has_both]) - theta_pred[has_both]) * day_per_rad; df$delta_doy[has_both] <- dd_b; df$delta_doy_resid[has_both] <- dd_b; df$hopkins_adjusted_doy_resid[has_both] <- dd_b }
  med_lat <- if (any(valid_lat)) stats::median(df$decimalLatitude[valid_lat], na.rm = TRUE) else NA_real_
  ref_doy_base <- if (is.finite(med_lat) && is.function(baseline$predict_doy_rad)) wrap_doy_365(baseline$predict_doy_rad(data.frame(decimalLatitude = med_lat))[1L] * day_per_rad) else circ_mean_doy
  if (!is.finite(ref_doy_base)) ref_doy_base <- circ_mean_doy
  df$hopkins_adjusted_doy <- ifelse(has_both & is.finite(ref_doy_base), wrap_doy_365(ref_doy_base + df$delta_doy), ifelse(valid_doy & df$is_verified_anthesis, wrap_doy_365(df$doy), NA_real_))
  vh <- is.finite(df$hopkins_adjusted_doy) & df$is_verified_anthesis
  df$sin_doy <- ifelse(vh, round(sin(2.0 * pi * df$hopkins_adjusted_doy / 365.25), 6), NA_real_)
  df$cos_doy <- ifelse(vh, round(cos(2.0 * pi * df$hopkins_adjusted_doy / 365.25), 6), NA_real_)
  is_unif <- isTRUE(baseline$is_uniform_null); df$pheno_timing_category <- "Unknown"
  if (is_unif) df$pheno_timing_category[has_both] <- "Aseasonal_Uniform_Null" else {
    df$pheno_timing_category[has_both & df$delta_doy < -7.0] <- "Early_Flowering"; df$pheno_timing_category[has_both & df$delta_doy > 7.0] <- "Late_Flowering"
    df$pheno_timing_category[has_both & abs(df$delta_doy) <= 7.0] <- "Synchronous_with_Cline"
  }
  df$pheno_timing_category[!df$is_verified_anthesis] <- "Excluded_Non_Anthesis"

  active_tx <- if (exists("get_active_packera_taxa", mode = "function")) get_active_packera_taxa(df$species_standardized)
    else sort(unique(as.character(df$species_standardized[!is.na(df$species_standardized) & nzchar(df$species_standardized)])))
  tf_mask <- df$species_standardized %in% active_tx & is.finite(df$delta_doy) & df$is_verified_anthesis
  tf_idx <- if (sum(tf_mask & is_active_ref) >= 5L) which(tf_mask & is_active_ref) else which(tf_mask)
  target_flowering <- df[tf_idx, , drop = FALSE]
  bimod_cohort <- test_allochronic_bimodality_von_mises(df$hopkins_adjusted_doy_resid[tf_idx], weights = calc_row_wts(tf_idx))
  df$allochronic_bimodality_delta_bic <- rep(bimod_cohort$allochronic_bimodality_delta_bic, n_r)
  if (length(active_tx) >= 2L) for (tx_b in active_tx) {
    idx_tb <- tf_idx[as.character(df$species_standardized[tf_idx]) == tx_b]
    if (length(idx_tb) >= 8L) { bt <- test_allochronic_bimodality_von_mises(df$hopkins_adjusted_doy_resid[idx_tb], weights = calc_row_wts(idx_tb)); if (is.finite(bt$allochronic_bimodality_delta_bic)) df$allochronic_bimodality_delta_bic[as.character(df$species_standardized) == tx_b] <- bt$allochronic_bimodality_delta_bic }
  }

  .w_csv <- function(x_df, path) { if (exists("write_csv_atomic", mode = "function")) return(write_csv_atomic(x_df, path, demo = demo)); dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE); utils::write.csv(x_df, file = path, row.names = FALSE); TRUE }
  if (!is.null(output_csv) && nzchar(trimws(as.character(output_csv)))) {
    keep_cols <- c("catalogNumber", "plant_individual_id", "institutionCode", "scientificName", "species_raw", "species_standardized",
                   "decimalLatitude", "decimalLongitude", "eventDate", "year", "month", "day", "doy", "expected_doy", "hopkins_expected",
                   "hopkins_adjusted_doy", "hopkins_adjusted_doy_resid", "delta_doy", "delta_doy_resid", "anthesis_dominance_weight",
                   "allochronic_bimodality_delta_bic", "pheno_timing_category", "is_pheno_uniform_null", "is_verified_anthesis",
                   "is_active_anthesis", "phenological_state", "phenophase_status", "anthesis_filter_status", "voucher_phenological_state", "regional_group")
    .w_csv(df[, intersect(keep_cols, names(df)), drop = FALSE], output_csv)
  }
  circ_th <- if (nrow(target_flowering) > 0L && any(is.finite(target_flowering$hopkins_adjusted_doy))) doy_to_rad(target_flowering$hopkins_adjusted_doy)
    else if (nrow(target_flowering) > 0L) wrap_rad(doy_to_rad(target_flowering$delta_doy)) else numeric(0)
  circ_pt <- run_circular_multigroup_phase_test(circ_th, target_flowering$species_standardized, B = 999L, seed = 42L)
  mk_srow <- function(sec, cmp, n_s, md = NA_real_, med = NA_real_, s1 = NA_real_, s2 = NA_real_, eff = NA_real_, pv = "", intp = "") {
    data.frame(Analysis_Section = sec, Taxon_or_Comparison = cmp, Sample_Size = as.integer(n_s), Mean_Delta_DOY = md, SD_Delta_DOY = NA_real_,
               Median_Delta_DOY = med, IQR_Delta_DOY = NA_real_, CI_Lower_95 = NA_real_, CI_Upper_95 = NA_real_,
               Statistic_Primary = s1, Statistic_Secondary = s2, Effect_or_R2 = eff, P_Value = pv, Interpretation = intp, stringsAsFactors = FALSE)
  }
  base_row <- mk_srow("Baseline_Regression", "DOY ~ decimalLatitude (Verified Anthesis)", baseline$n, s1 = round(baseline$slope, 3), s2 = round(baseline$intercept, 3), eff = round(baseline$r_squared, 4),
                      pv = if (!is.na(baseline$p_value)) sprintf("%.4e", baseline$p_value) else "NA",
                      intp = if (is_unif) sprintf("Uniform circular phase null retained (Rayleigh p=%.3e, R_bar=%.3f; selected %s)", baseline$rayleigh_p, baseline$r_bar, baseline$selected_model)
                        else sprintf("Circular AICc selected [%s] slope = %+.2f d/deg N (Hopkins circular)", baseline$selected_model, baseline$slope))
  ab_int <- if (!is.finite(bimod_cohort$allochronic_bimodality_delta_bic)) "Insufficient active verified-anthesis sample size for K=1 vs K=2 von Mises mixture"
    else if (bimod_cohort$allochronic_bimodality_delta_bic >= 6.0 && is.finite(bimod_cohort$ashman_d) && bimod_cohort$ashman_d >= 2.0)
      sprintf("Strong bimodal allochronic isolation signal (Delta_BIC=%+.2f, Ashman_D=%.2f, modes %.1f vs %.1f d)", bimod_cohort$allochronic_bimodality_delta_bic, bimod_cohort$ashman_d, bimod_cohort$mu1_days, bimod_cohort$mu2_days)
    else sprintf("Unimodal / continuous phenological variation retained (Delta_BIC=%+.2f, Ashman_D=%.2f)", bimod_cohort$allochronic_bimodality_delta_bic, ifelse(is.finite(bimod_cohort$ashman_d), bimod_cohort$ashman_d, 0.0))
  bimod_row <- mk_srow("Allochronic_Bimodality_Test", "Cohort Verified Anthesis (K=1 vs K=2 von Mises)", bimod_cohort$n, md = bimod_cohort$mu1_days, med = bimod_cohort$mu2_days,
                       s1 = bimod_cohort$allochronic_bimodality_delta_bic, s2 = bimod_cohort$ashman_d,
                       eff = if (is.finite(bimod_cohort$mu1_days) && is.finite(bimod_cohort$mu2_days)) round(abs(bimod_cohort$mu2_days - bimod_cohort$mu1_days), 2) else NA_real_, intp = ab_int)
  c_rows <- summarize_weighted_circular_phenology(target_flowering, active_tx, circ_pt, is_unif = is_unif, B = 999L, seed = 42L)
  sr <- if (is.data.frame(c_rows) && nrow(c_rows) > 0L) rbind(base_row, c_rows, bimod_row) else rbind(base_row, bimod_row)
  sr$Selected_Hopkins_Model <- as.character(baseline$selected_model)
  sr$Hopkins_AICc_Weight <- if (!is.null(baseline$aicc_weight) && is.finite(baseline$aicc_weight)) round(as.numeric(baseline$aicc_weight), 6) else NA_real_
  sr$Rayleigh_R_Bar <- if (!is.null(baseline$r_bar) && is.finite(baseline$r_bar)) round(as.numeric(baseline$r_bar), 6) else NA_real_
  sr$Rayleigh_P <- if (!is.null(baseline$rayleigh_p) && is.finite(baseline$rayleigh_p)) as.numeric(baseline$rayleigh_p) else NA_real_
  sr$Circ_Watson_F <- if (!is.null(circ_pt$watson_f) && is.finite(circ_pt$watson_f)) round(as.numeric(circ_pt$watson_f), 4) else NA_real_
  sr$Circ_Perm_P <- if (!is.null(circ_pt$perm_p) && is.finite(circ_pt$perm_p)) as.numeric(circ_pt$perm_p) else NA_real_
  sr$Allochronic_Bimodality_Delta_BIC <- bimod_cohort$allochronic_bimodality_delta_bic; sr$Allochronic_Ashman_D <- bimod_cohort$ashman_d
  if (!is.null(summary_csv) && nzchar(trimws(as.character(summary_csv)))) {
    .w_csv(sr, summary_csv)
    if (!is.null(baseline$model_selection_table) && is.data.frame(baseline$model_selection_table) && nrow(baseline$model_selection_table) > 0L) .w_csv(baseline$model_selection_table, file.path(dirname(summary_csv), "hopkins_circular_model_selection.csv"))
  }
  list(df = df, baseline = baseline, anova = NULL, tukey = NULL, circular_phase_test = circ_pt, allochronic_bimodality = bimod_cohort, summary = sr)
}

export_phenological_latitudinal_figures <- function(df, baseline, out_path, pheno_summary = NULL, bimod_info = NULL) {
  if (!exists("render_step05b_hopkins_phenology_png", mode = "function")) {
    for (.f in c("scripts/figures/theme_packera_pub.R", "scripts/figures/step05_figures.R")) if (file.exists(.f)) source(.f)
  }
  out_png <- if (!is.null(out_path) && nzchar(out_path)) normalize_png_path(out_path) else "outputs/figures/fig05b_hopkins_circular_phenology.png"
  message(sprintf("Generating 300-DPI circular & Hopkins bioclimatic phenology figure -> %s", out_png))
  render_step05b_hopkins_phenology_png(df = df, baseline_info = baseline, out_pheno_path = out_png, pheno_summary = pheno_summary, bimod_info = bimod_info)
}
