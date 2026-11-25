if (grDevices::dev.cur() == 1L) grDevices::pdf(NULL)
# ==============================================================================
# Module: scripts/morphometrics/auth_osr_metrics.R
# Open-set recognition (OSR) congeneric Mahalanobis metrics and canonical GMM relabeling
# ==============================================================================

.compute_open_set_recognition_metrics <- function(df, X_mat, S_inv, d2_mat, active_taxa, taxa_vec, gmm_cl) {
  n <- nrow(X_mat)
  p <- max(1L, ncol(X_mat))
  min_congener_mahal_d2 <- rep(NA_real_, n)
  min_congener_taxon <- rep(NA_character_, n)
  open_set_novelty_p <- rep(NA_real_, n)
  open_set_df <- rep(as.integer(p), n)
  if (n == 0L) {
    return(list(min_congener_mahal_d2 = min_congener_mahal_d2, min_congener_taxon = min_congener_taxon, open_set_novelty_p = open_set_novelty_p, open_set_df = integer(0)))
  }

  sp_raw <- if (!is.null(df) && "species_standardized" %in% names(df)) as.character(df$species_standardized) else as.character(taxa_vec)
  sp_taxa <- active_taxa[grepl("^Packera[[:space:]]+", active_taxa)]
  if (length(sp_taxa) >= 2L && !is.null(d2_mat) && ncol(d2_mat) >= 2L) {
    sp_cols <- match(sp_taxa, active_taxa)
    d2_sp <- d2_mat[, sp_cols, drop = FALSE]
    min_all <- apply(d2_sp, 1L, min, na.rm = TRUE)
    min_all[!is.finite(min_all) | min_all < 0.0] <- 0.0
    open_set_novelty_p <- round(1.0 - stats::pchisq(min_all, df = p), 6)
    for (i in seq_len(n)) {
      yi_sp <- if (!is.na(sp_raw[i]) && grepl("^Packera[[:space:]]+", sp_raw[i])) sp_raw[i] else taxa_vec[i]
      oth_idx <- which(sp_taxa != yi_sp)
      if (length(oth_idx) >= 1L) {
        d_sub <- d2_sp[i, oth_idx]
        if (any(is.finite(d_sub))) {
          best_o <- oth_idx[which.min(d_sub)]
          min_congener_mahal_d2[i] <- round(as.numeric(d2_sp[i, best_o]), 4)
          min_congener_taxon[i] <- sp_taxa[best_o]
        }
      }
    }
    return(list(min_congener_mahal_d2 = min_congener_mahal_d2, min_congener_taxon = min_congener_taxon, open_set_novelty_p = open_set_novelty_p, open_set_df = open_set_df))
  }

  # Single-taxon intraspecific mode: query canonical SQLite for congeneric Packera centroids
  # in run-invariant standardized harmonic + vestiture + topology space,
  # falling back strictly to same-run_id PC1..3 when harmonics are absent.
  # Isolate early return(NULL) and on.exit cleanup inside an IIFE so failed DB reads fall back to local OSR.
  db_osr <- tryCatch((function() {
    if (!is.null(df) && exists("connect_packera_db", mode = "function") && requireNamespace("DBI", quietly = TRUE)) {
      con <- connect_packera_db(read_only = TRUE)
      if (is.null(con)) return(NULL)
      on.exit(tryCatch(DBI::dbDisconnect(con), error = function(e) NULL), add = TRUE)
      has_sheet_delta_inv <- (("blade_cie_l_q85_sheet_delta" %in% names(df)) &&
        (mean(is.finite(suppressWarnings(as.numeric(df$blade_cie_l_q85_sheet_delta)))) >= 0.50)) ||
        (("blade_cie_l_sheet_delta" %in% names(df)) &&
        (mean(is.finite(suppressWarnings(as.numeric(df$blade_cie_l_sheet_delta)))) >= 0.50))
      phot_inv_cands <- c(
        if (!has_sheet_delta_inv) "blade_cie_l_norm",
        "blade_cie_l_q85_sheet_delta", "basal_zone_cie_l_delta",
        "blade_cie_l_sheet_delta",
        "blade_cie_a_browning_delta", "blade_cie_b_yellowing_delta"
      )
      if (!("blade_cie_a_browning_delta" %in% names(df))) phot_inv_cands <- c(phot_inv_cands, "blade_cie_a_delta")
      if (!("blade_cie_b_yellowing_delta" %in% names(df))) phot_inv_cands <- c(phot_inv_cands, "blade_cie_b_delta")
      inv_cands <- intersect(
        c(
          paste0("A", 1:6), paste0("D", 1:6),
          "solidity", "dissection_index",
          "basal_cordate_index", "terminal_lobe_area_frac", "margin_serration_index",
          phot_inv_cands,
          "ECT_PC1", "ECT_PC2", "ect_lobe_complexity",
          "AR_inv", "radiate_head_fraction"
        ),
        names(df)
      )
      inv_cands <- inv_cands[vapply(inv_cands, function(cl) is.numeric(df[[cl]]) && mean(is.finite(df[[cl]])) >= 0.50, logical(1))]
      raw_pc_cands <- intersect(c("PC1", "PC2", "PC3"), names(df))
      raw_pc_cands <- raw_pc_cands[vapply(raw_pc_cands, function(cl) is.numeric(df[[cl]]) && mean(is.finite(df[[cl]])) >= 0.80, logical(1))]
      active_rid <- if (exists("resolve_packera_run_id", mode = "function")) tryCatch(as.character(resolve_packera_run_id(NULL)), error = function(e) "") else ""

      eval_db_space <- function(cand_cols, require_same_run = FALSE) {
        if (length(cand_cols) < 2L) return(NULL)
        db_df <- NULL
        v_geo_tbl <- if (DBI::dbExistsTable(con, "vouchers")) "vouchers" else if (DBI::dbExistsTable(con, "specimens")) "specimens" else ""
        has_geo_v <- nzchar(v_geo_tbl) && all(c("catalogNumber", "decimalLatitude", "decimalLongitude") %in% DBI::dbListFields(con, v_geo_tbl))
        if (DBI::dbExistsTable(con, "otu_morphometrics")) {
          flds <- DBI::dbListFields(con, "otu_morphometrics")
          use_c <- intersect(cand_cols, flds)
          if (length(use_c) >= 2L && "species_standardized" %in% flds) {
            geo_sel <- if (has_geo_v && "catalogNumber" %in% flds) ", v_geo.decimalLatitude AS _osr_lat, v_geo.decimalLongitude AS _osr_lon" else ""
            geo_join <- if (has_geo_v && "catalogNumber" %in% flds) sprintf(" LEFT JOIN %s v_geo ON o.catalogNumber = v_geo.catalogNumber", v_geo_tbl) else ""
            q_cols <- paste0(paste(c("o.species_standardized", paste0("o.", use_c)), collapse = ", "), geo_sel)
            where_parts <- c("o.species_standardized LIKE 'Packera %'")
            if ("authority_weight" %in% flds) where_parts <- c(where_parts, "COALESCE(o.authority_weight, 1.0) > 0")
            if ("morph_tier_source" %in% flds) where_parts <- c(where_parts, "COALESCE(o.morph_tier_source, 'Basal_Rosette') NOT IN ('Cauline_Only_Fallback', 'Transitional_Passive')")
            if (isTRUE(require_same_run) && "run_id" %in% flds) {
              if (!nzchar(active_rid)) return(NULL)
              sql <- sprintf("SELECT %s FROM otu_morphometrics o%s WHERE %s AND o.run_id = $1", q_cols, geo_join, paste(where_parts, collapse = " AND "))
              db_df <- DBI::dbGetQuery(con, sql, params = list(active_rid))
            } else if (!isTRUE(require_same_run)) {
              sql <- sprintf("SELECT %s FROM otu_morphometrics o%s WHERE %s", q_cols, geo_join, paste(where_parts, collapse = " AND "))
              db_df <- DBI::dbGetQuery(con, sql)
            }
          }
        }
        if ((is.null(db_df) || length(unique(db_df$species_standardized)) < 2L) && !isTRUE(require_same_run) &&
            DBI::dbExistsTable(con, "leaf_efa_harmonics") && (DBI::dbExistsTable(con, "specimens") || DBI::dbExistsTable(con, "vouchers"))) {
          v_tbl <- if (DBI::dbExistsTable(con, "specimens")) "specimens" else "vouchers"
          flds_h <- DBI::dbListFields(con, "leaf_efa_harmonics")
          flds_v <- DBI::dbListFields(con, v_tbl)
          use_c <- intersect(cand_cols, flds_h)
          if (length(use_c) >= 2L && "species_standardized" %in% flds_v) {
            geo_sel_h <- if (all(c("decimalLatitude", "decimalLongitude") %in% flds_v)) ", v.decimalLatitude AS _osr_lat, v.decimalLongitude AS _osr_lon" else ""
            q_cols <- paste0(paste(paste0("h.", use_c), collapse = ", "), geo_sel_h)
            lz_clause <- if ("leaf_zone" %in% flds_h) " AND COALESCE(LOWER(h.leaf_zone), 'basal') != 'cauline'" else ""
            db_df <- DBI::dbGetQuery(con, sprintf(
              "SELECT v.species_standardized, %s FROM leaf_efa_harmonics h JOIN %s v ON h.catalogNumber = v.catalogNumber WHERE v.species_standardized LIKE 'Packera %%'%s",
              q_cols, v_tbl, lz_clause
            ))
          }
        }
        if (is.null(db_df) || nrow(db_df) < 6L) return(NULL)
        use_c <- intersect(cand_cols, names(db_df))
        if (length(use_c) < 2L) return(NULL)
        for (uc in use_c) db_df[[uc]] <- suppressWarnings(as.numeric(db_df[[uc]]))
        cc <- stats::complete.cases(db_df[, use_c, drop = FALSE])
        db_df <- db_df[cc & !is.na(db_df$species_standardized) & grepl("^Packera[[:space:]]+", db_df$species_standardized), , drop = FALSE]
        focal_sp <- unique(sp_raw[!is.na(sp_raw) & grepl("^Packera[[:space:]]+", sp_raw)])
        if (all(c("_osr_lat", "_osr_lon") %in% names(db_df)) && all(c("decimalLatitude", "decimalLongitude") %in% names(df))) {
          d_lat <- suppressWarnings(as.numeric(df$decimalLatitude)); d_lon <- suppressWarnings(as.numeric(df$decimalLongitude))
          act_g <- is.finite(d_lat) & is.finite(d_lon) & (abs(d_lat) <= 90) & (abs(d_lon) <= 180)
          if (sum(act_g) >= 3L) {
            lat_b <- range(d_lat[act_g]) + c(-12.0, 12.0); lon_b <- range(d_lon[act_g]) + c(-15.0, 15.0)
            b_lat <- suppressWarnings(as.numeric(db_df[["_osr_lat"]])); b_lon <- suppressWarnings(as.numeric(db_df[["_osr_lon"]]))
            in_box <- is.finite(b_lat) & is.finite(b_lon) & b_lat >= lat_b[1] & b_lat <= lat_b[2] & b_lon >= lon_b[1] & b_lon <= lon_b[2]
            if (length(setdiff(unique(db_df$species_standardized[in_box]), focal_sp)) >= 2L) db_df <- db_df[in_box, , drop = FALSE]
          }
        }
        oth_sp <- setdiff(unique(db_df$species_standardized), focal_sp)
        if (length(oth_sp) < 2L) return(NULL)
        X_loc_raw <- as.matrix(df[, use_c, drop = FALSE]); storage.mode(X_loc_raw) <- "double"
        X_db_raw <- as.matrix(db_df[, use_c, drop = FALSE]); storage.mode(X_db_raw) <- "double"
        comb_sd <- vapply(seq_len(ncol(X_db_raw)), function(jc) stats::sd(c(X_loc_raw[, jc], X_db_raw[, jc]), na.rm = TRUE), numeric(1))
        keep_v <- which(is.finite(comb_sd) & comb_sd > 1e-6)
        if (length(keep_v) < 2L) return(NULL)
        use_c <- use_c[keep_v]
        sc_vec <- comb_sd[keep_v]
        X_loc_s <- sweep(X_loc_raw[, keep_v, drop = FALSE], 2L, sc_vec, "/")
        X_db_s <- sweep(X_db_raw[, keep_v, drop = FALSE], 2L, sc_vec, "/")
        X_loc_imp <- X_loc_s
        for (jc in seq_len(ncol(X_loc_imp))) {
          bad_c <- !is.finite(X_loc_imp[, jc])
          if (any(bad_c)) X_loc_imp[bad_c, jc] <- stats::median(c(X_loc_imp[!bad_c, jc], X_db_s[, jc]), na.rm = TRUE)
        }
        all_sp <- unique(c(focal_sp, oth_sp))
        cent_db <- matrix(0.0, nrow = length(all_sp), ncol = length(use_c), dimnames = list(all_sp, use_c))
        resid_list <- vector("list", length(all_sp))
        for (si in seq_along(all_sp)) {
          sp_nm <- all_sp[si]
          sub_m <- if (sp_nm %in% db_df$species_standardized) X_db_s[db_df$species_standardized == sp_nm, , drop = FALSE] else X_loc_imp
          med_s <- apply(sub_m, 2L, stats::median, na.rm = TRUE)
          cent_db[si, ] <- med_s
          resid_list[[si]] <- sweep(sub_m, 2L, med_s, "-")
        }
        R_pool <- do.call(rbind, resid_list)
        S_raw_db <- if (nrow(R_pool) >= 3L) stats::cov(R_pool) else diag(1.0, ncol(X_loc_imp))
        S_raw_db[!is.finite(S_raw_db)] <- 0.0
        S_pool_db <- compute_ledoit_wolf_cov(S_raw_db + diag(1e-5, ncol(X_loc_imp)), X_centered = R_pool, n_obs = nrow(R_pool))
        S_db_inv <- tryCatch(solve(S_pool_db), error = function(e) diag(1.0 / pmax(diag(S_pool_db), 1e-4), ncol(X_loc_imp)))
        d2_all <- matrix(NA_real_, nrow = n, ncol = length(all_sp), dimnames = list(NULL, all_sp))
        for (si in seq_along(all_sp)) {
          df_s <- sweep(X_loc_imp, 2L, cent_db[si, ], "-")
          d2_all[, si] <- pmax(0.0, rowSums((df_s %*% S_db_inv) * df_s))
        }
        df_vec <- rep(ncol(X_loc_imp), n)
        # Re-evaluate any multi-island rows lacking harmonics strictly on their observed topological subset
        topo_pos <- which(!(use_c %in% c(paste0("A", 1:24), paste0("D", 1:24))))
        incomplete_rows <- which(!stats::complete.cases(X_loc_s))
        if (length(incomplete_rows) > 0L && length(topo_pos) >= 2L) {
          S_sub <- S_pool_db[topo_pos, topo_pos, drop = FALSE]
          S_sub_inv <- tryCatch(solve(S_sub), error = function(e) diag(1.0 / pmax(diag(S_sub), 1e-4), length(topo_pos)))
          for (ri in incomplete_rows) {
            obs_r <- topo_pos[is.finite(X_loc_s[ri, topo_pos])]
            if (length(obs_r) >= 2L) {
              S_r <- S_pool_db[obs_r, obs_r, drop = FALSE]
              S_r_inv <- tryCatch(solve(S_r), error = function(e) diag(1.0 / pmax(diag(S_r), 1e-4), length(obs_r)))
              for (si in seq_along(all_sp)) {
                dv_r <- X_loc_s[ri, obs_r] - cent_db[si, obs_r]
                d2_all[ri, si] <- max(0.0, as.numeric(t(dv_r) %*% S_r_inv %*% dv_r))
              }
              df_vec[ri] <- length(obs_r)
            }
          }
        }
        oth_cols <- match(oth_sp, all_sp)
        d2_oth <- d2_all[, oth_cols, drop = FALSE]
        best_c <- max.col(-d2_oth, ties.method = "first")
        list(
          min_congener_mahal_d2 = round(d2_oth[cbind(seq_len(n), best_c)], 4),
          min_congener_taxon = oth_sp[best_c],
          open_set_novelty_p = {
            min_d2_r <- vapply(seq_len(n), function(ri) { vr <- d2_all[ri, ]; vr <- vr[is.finite(vr)]; if (length(vr) > 0L) min(vr) else NA_real_ }, numeric(1))
            round(1.0 - stats::pchisq(min_d2_r, df = df_vec), 6)
          },
          open_set_df = as.integer(df_vec)
        )
      }

      res_inv <- if (length(inv_cands) >= 3L) eval_db_space(inv_cands, require_same_run = FALSE) else NULL
      if (!is.null(res_inv)) {
        res_inv
      } else if (length(raw_pc_cands) >= 2L) {
        eval_db_space(raw_pc_cands, require_same_run = TRUE)
      } else {
        NULL
      }
    } else NULL
  })(), error = function(e) NULL)

  if (!is.null(db_osr)) return(db_osr)

  core_rows <- which(!is.na(gmm_cl) & gmm_cl == 1L)
  if (length(core_rows) < 2L) core_rows <- seq_len(n)
  mu_core <- colMeans(X_mat[core_rows, , drop = FALSE])
  diff_c <- sweep(X_mat, 2L, mu_core, "-")
  d2_core <- pmax(0.0, rowSums((diff_c %*% S_inv) * diff_c))
  open_set_novelty_p <- round(1.0 - stats::pchisq(d2_core, df = p), 6)
  list(min_congener_mahal_d2 = min_congener_mahal_d2, min_congener_taxon = min_congener_taxon, open_set_novelty_p = open_set_novelty_p, open_set_df = open_set_df)
}


.canonical_relabel_gmm_clusters <- function(sub_v, cls_sub, best_k, ref_cand, core_idx = NULL) {
  n_eff_ratio <- 1.0
  w_all <- if ("authority_weight" %in% names(sub_v) && any(is.finite(suppressWarnings(as.numeric(sub_v$authority_weight))))) {
    suppressWarnings(as.numeric(sub_v$authority_weight))
  } else {
    rep(1.0, nrow(sub_v))
  }
  w_all[!is.finite(w_all) | !ref_cand] <- 0.0
  if (all(w_all <= 0.0)) w_all[ref_cand] <- 1.0
  if (!is.null(core_idx) && length(core_idx) >= 2L) {
    w_c <- w_all[core_idx]
    w_c <- w_c[is.finite(w_c) & w_c > 0.0]
    if (length(w_c) >= 2L && sum(w_c^2) > 1e-12) {
      n_eff_ratio <- min(1.0, max(0.1, ((sum(w_c)^2) / sum(w_c^2)) / length(w_c)))
    }
  }
  if (!is.null(cls_sub) && best_k > 1L) {
    u_cl <- sort(unique(stats::na.omit(as.integer(cls_sub))))
    if (length(u_cl) > 1L) {
      cl_w <- tapply(w_all, factor(cls_sub, levels = u_cl), sum, na.rm = TRUE)
      cl_w[is.na(cl_w)] <- 0.0
      t_info <- .identify_type_specimen_vec(sub_v)
      t_idx <- which(t_info$is_type_specimen & !is.na(cls_sub))
      if (length(t_idx) > 0L) {
        t_wt <- pmax(1.0, w_all[t_idx]) * ifelse(t_info$type_rank[t_idx] == 2L, 2.0, 1.0)
        t_cl_w <- tapply(t_wt, factor(cls_sub[t_idx], levels = u_cl), sum, na.rm = TRUE)
        t_cl_w[is.na(t_cl_w)] <- 0.0
        ord_cl <- order(-t_cl_w, -cl_w, u_cl)
      } else {
        ord_cl <- order(-cl_w, u_cl)
      }
      cls_sub <- as.integer(match(as.integer(cls_sub), u_cl[ord_cl]))
    }
  }
  list(cls_sub = cls_sub, n_eff_ratio = n_eff_ratio)
}
