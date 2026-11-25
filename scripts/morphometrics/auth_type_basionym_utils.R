if (grDevices::dev.cur() == 1L) grDevices::pdf(NULL)
# ==============================================================================
# Module: scripts/morphometrics/auth_type_basionym_utils.R
# Basionym canonicalization, homotypic checks, and type/basionym centroid builder
# ==============================================================================

.canonicalize_type_basionym_vec <- function(df, fallback_taxa = NULL) {
  n <- if (is.null(df)) 0L else nrow(df)
  if (n == 0L) return(character(0))
  fb <- if (!is.null(fallback_taxa) && length(fallback_taxa) == n) {
    as.character(fallback_taxa)
  } else if ("species_standardized" %in% names(df)) {
    as.character(df$species_standardized)
  } else {
    rep("Unknown", n)
  }
  raw_v <- rep(NA_character_, n)
  if ("typified_basionym" %in% names(df)) {
    tb <- trimws(as.character(df$typified_basionym))
    hit_tb <- !is.na(tb) & nzchar(tb) & !(tolower(tb) %in% c("na", "nan", "null", "none", "unknown"))
    raw_v[hit_tb] <- tb[hit_tb]
  }
  for (tf in intersect(c("type_status", "typeStatus"), names(df))) {
    ts_s <- trimws(as.character(df[[tf]]))
    has_of <- is.na(raw_v) & !is.na(ts_s) & grepl("type[[:space:]]+of[[:space:]]+[A-Z][a-z]+[[:space:]]+[a-z][a-z-]+", ts_s, ignore.case = TRUE)
    if (any(has_of)) {
      raw_v[has_of] <- sub(".*type[[:space:]]+of[[:space:]]+", "", ts_s[has_of], ignore.case = TRUE)
    }
  }
  for (sc in intersect(c("species_raw", "verbatim_taxon_concept", "scientificName", "taxon", "species_standardized"), names(df))) {
    cand <- trimws(as.character(df[[sc]]))
    hit <- is.na(raw_v) & !is.na(cand) & nzchar(cand) & !(tolower(cand) %in% c("na", "nan", "null", "none", "unknown"))
    raw_v[hit] <- cand[hit]
  }
  raw_v[is.na(raw_v) | !nzchar(raw_v)] <- fb[is.na(raw_v) | !nzchar(raw_v)]
  parse_one <- function(s, fb_one) {
    if (is.na(s) || !nzchar(s)) return(fb_one)
    s_clean <- gsub("[[:space:]]+", " ", trimws(s))
    m_bin <- regmatches(s_clean, regexpr("^[A-Z][a-z]+[[:space:]]+[a-z][a-z-]+", s_clean))
    m_infra <- regmatches(s_clean, regexpr("^[A-Z][a-z]+[[:space:]]+[a-z][a-z-]+[[:space:]]+(var\\.|subsp\\.|ssp\\.)[[:space:]]+[a-z][a-z-]+", s_clean))
    if (length(m_infra) > 0L && nzchar(m_infra[1L]) && length(m_bin) > 0L) {
      # Strip botanical autonyms where infraspecific epithet duplicates specific epithet.
      toks_i <- strsplit(m_infra[1L], "[[:space:]]+")[[1L]]
      if (length(toks_i) >= 4L && tolower(toks_i[4L]) != tolower(toks_i[2L])) return(m_infra[1L])
    }
    m_hyb <- regmatches(s_clean, regexpr("^[A-Z][a-z]+[[:space:]]+[x×][[:space:]]*[a-z][a-z-]+", s_clean))
    if (length(m_hyb) > 0L && nzchar(m_hyb[1L])) return(gsub("[[:space:]]+", " ", m_hyb[1L]))
    if (length(m_bin) > 0L && nzchar(m_bin[1L])) return(m_bin[1L])
    fb_one
  }
  vapply(seq_len(n), function(i) parse_one(raw_v[i], fb[i]), character(1))
}

.extract_terminal_epithet_stem <- function(taxon_vec) {
  vapply(as.character(taxon_vec), function(s) {
    if (is.na(s) || !nzchar(trimws(s))) return("")
    toks <- strsplit(tolower(trimws(s)), "[[:space:]]+")[[1L]]
    toks <- toks[!grepl("^(var\\.?|subsp\\.?|ssp\\.?|f\\.?|x|×)$", toks)]
    if (length(toks) < 2L) return(toks[length(toks)])
    ep <- gsub("[^a-z]", "", toks[length(toks)])
    sub("(us|a|um|is|e)$", "", ep)
  }, character(1), USE.NAMES = FALSE)
}

# Check homotypic basionym identity honoring concept_synonymy_status, backbone_basionym_name, and replacement pairs.
.is_homotypic_sensu_stricto_basionym <- function(basionym_vec, accepted_vec, df = NULL) {
  n_b <- length(basionym_vec)
  if (n_b == 0L) return(logical(0))
  b_low <- gsub("[[:space:]]+", " ", tolower(trimws(as.character(basionym_vec))))
  a_low <- gsub("[[:space:]]+", " ", tolower(trimws(as.character(accepted_vec))))
  b_stem <- .extract_terminal_epithet_stem(basionym_vec)
  a_stem <- .extract_terminal_epithet_stem(accepted_vec)
  out <- nzchar(b_stem) & nzchar(a_stem) & (b_stem == a_stem)
  known_dubia <- (a_low == "packera dubia") & (b_low %in% c("cineraria dubia", "senecio dubius", "senecio tomentosus", "packera tomentosa"))
  out <- out | known_dubia
  if (!is.null(df) && nrow(df) == n_b) {
    if ("backbone_basionym_name" %in% names(df)) {
      bb <- gsub("[[:space:]]+", " ", tolower(trimws(as.character(df$backbone_basionym_name))))
      for (ii in seq_len(n_b)) {
        if (!is.na(bb[ii]) && nzchar(bb[ii]) && nzchar(b_low[ii]) && (bb[ii] == b_low[ii] || startsWith(bb[ii], paste0(b_low[ii], " ")))) {
          out[ii] <- TRUE
        }
      }
    }
    if ("concept_synonymy_status" %in% names(df)) {
      cs <- trimws(as.character(df$concept_synonymy_status))
      out[!is.na(cs) & (cs == "Heterotypic_Synonym")] <- FALSE
      out[!is.na(cs) & (cs %in% c("Homotypic_Synonym", "Accepted_Sensu_Stricto"))] <- TRUE
    }
  }
  out
}

.filter_homotypic_type_rows <- function(df, type_rows, accepted_taxon) {
  if (length(type_rows) <= 1L || is.null(df) || nrow(df) == 0L) return(type_rows)
  df_sub <- df[type_rows, , drop = FALSE]
  bas_sub <- .canonicalize_type_basionym_vec(df_sub, rep(accepted_taxon, length(type_rows)))
  is_homo <- .is_homotypic_sensu_stricto_basionym(bas_sub, rep(accepted_taxon, length(type_rows)), df = df_sub)
  if (length(unique(bas_sub)) <= 1L && all(is_homo == is_homo[1L])) return(type_rows)
  if (any(is_homo)) type_rows[is_homo] else type_rows
}

.compute_type_pool_centroid <- function(rows_sub, X_src, S_inv_src, df, has_efa_row, type_rank, prefer_efa = TRUE) {
  if (isTRUE(prefer_efa) && any(has_efa_row[rows_sub])) {
    rows_sub <- rows_sub[has_efa_row[rows_sub]]
  }
  prim_s <- rows_sub[type_rank[rows_sub] == 2L]
  pool_s <- if (length(rows_sub) >= 3L && length(prim_s) > 0L) prim_s else rows_sub
  mu_s <- if (length(rows_sub) == 1L) {
    as.numeric(X_src[rows_sub, ])
  } else if (length(rows_sub) == 2L) {
    w_t <- ifelse(type_rank[rows_sub] == 2L, 2.0, 1.0)
    as.numeric(colSums(X_src[rows_sub, , drop = FALSE] * w_t) / sum(w_t))
  } else {
    .compute_robust_taxon_medoid(X_src[pool_s, , drop = FALSE], S_inv = S_inv_src)
  }
  list(mu = mu_s, pool = pool_s, n_eff = .count_effective_type_gatherings(df, pool_s))
}

.build_taxon_type_basionym_centroids <- function(
  df, X_mat, S_inv, active_taxa, taxa_vec, is_type_specimen, is_t1_pure,
  has_efa_row, type_rank, det_tier, is_noise_vec, auth_w,
  basionym_vec_all, gmm_cl, u_clusters, cl_centroids
) {
  G <- length(active_taxa)
  p <- ncol(X_mat)
  compute_pool_centroid <- function(rows_sub, X_src = X_mat, S_inv_src = S_inv, prefer_efa = TRUE) {
    .compute_type_pool_centroid(rows_sub, X_src, S_inv_src, df, has_efa_row, type_rank, prefer_efa = prefer_efa)
  }
  type_centroids <- matrix(0.0, nrow = G, ncol = p, dimnames = list(active_taxa, colnames(X_mat)))
  has_emp_map <- stats::setNames(rep(FALSE, G), active_taxa)
  has_emp_homo_map <- stats::setNames(rep(FALSE, G), active_taxa)
  ss_type_rows_map <- stats::setNames(vector("list", G), active_taxa)
  typical_cl_map <- stats::setNames(rep(NA_integer_, G), active_taxa)
  typical_cl_list <- stats::setNames(vector("list", G), active_taxa)
  n_eff_type_vec <- stats::setNames(numeric(G), active_taxa)
  type_anchor_rows_list <- stats::setNames(vector("list", G), active_taxa)
  basionym_entries <- list()
  for (k in seq_len(G)) {
    tx_k <- active_taxa[k]
    tx_rows <- which(taxa_vec == tx_k)
    type_rows_k <- tx_rows[is_type_specimen[tx_rows] & is_t1_pure[tx_rows]]
    if (length(type_rows_k) == 0L) {
      type_rows_k <- tx_rows[is_type_specimen[tx_rows]]
    }
    ss_type_rows_k <- type_rows_k
    if (length(type_rows_k) > 0L) {
      has_emp_map[tx_k] <- TRUE
      bas_k <- basionym_vec_all[type_rows_k]
      u_bas_k <- unique(bas_k[!is.na(bas_k) & nzchar(bas_k)])
      if (length(u_bas_k) == 0L) u_bas_k <- tx_k
      bas_k_list <- lapply(u_bas_k, function(bn) {
        b_rows <- type_rows_k[bas_k == bn]
        if (length(b_rows) == 0L) b_rows <- type_rows_k
        cp <- compute_pool_centroid(b_rows)
        is_homo <- all(.is_homotypic_sensu_stricto_basionym(rep(bn, length(b_rows)), rep(tx_k, length(b_rows)), df = df[b_rows, , drop = FALSE]))
        list(
          taxon = tx_k,
          basionym = bn,
          is_homotypic = is_homo,
          rows = b_rows,
          pool = cp$pool,
          n_eff = cp$n_eff,
          max_rank = max(type_rank[b_rows], na.rm = TRUE),
          mu = cp$mu
        )
      })
      for (be in bas_k_list) {
        basionym_entries[[length(basionym_entries) + 1L]] <- be
      }
      homo_idx <- which(vapply(bas_k_list, function(x) isTRUE(x$is_homotypic), logical(1L)))
      if (length(homo_idx) >= 1L) has_emp_homo_map[tx_k] <- TRUE
      if (length(bas_k_list) == 1L) {
        mu_k <- bas_k_list[[1L]]$mu
        type_anchor_rows_list[[tx_k]] <- bas_k_list[[1L]]$pool
        n_eff_type_vec[tx_k] <- bas_k_list[[1L]]$n_eff
        ss_type_rows_k <- bas_k_list[[1L]]$rows
      } else {
        if (length(homo_idx) >= 1L) {
          ss_rows <- unlist(lapply(bas_k_list[homo_idx], `[[`, "rows"), use.names = FALSE)
          cp_ss <- compute_pool_centroid(ss_rows)
          mu_k <- cp_ss$mu
          type_anchor_rows_list[[tx_k]] <- cp_ss$pool
          n_eff_type_vec[tx_k] <- cp_ss$n_eff
          ss_type_rows_k <- ss_rows
        } else {
          mu_mat_b <- do.call(rbind, lapply(bas_k_list, `[[`, "mu"))
          max_pair_d2 <- 0.0
          for (bi in seq_len(nrow(mu_mat_b) - 1L)) {
            for (bj in (bi + 1L):nrow(mu_mat_b)) {
              dv_b <- mu_mat_b[bi, ] - mu_mat_b[bj, ]
              max_pair_d2 <- max(max_pair_d2, as.numeric(t(dv_b) %*% S_inv %*% dv_b))
            }
          }
          if (max_pair_d2 > stats::qchisq(0.95, df = max(1L, p))) {
            ord_b <- order(
              vapply(bas_k_list, `[[`, integer(1L), "max_rank"),
              vapply(bas_k_list, `[[`, numeric(1L), "n_eff"),
              decreasing = TRUE
            )
            best_be <- bas_k_list[[ord_b[1L]]]
            mu_k <- best_be$mu
            type_anchor_rows_list[[tx_k]] <- best_be$pool
            n_eff_type_vec[tx_k] <- best_be$n_eff
            ss_type_rows_k <- best_be$rows
          } else {
            cp_all <- compute_pool_centroid(type_rows_k)
            mu_k <- cp_all$mu
            type_anchor_rows_list[[tx_k]] <- cp_all$pool
            n_eff_type_vec[tx_k] <- cp_all$n_eff
            ss_type_rows_k <- type_rows_k
          }
        }
      }
      ss_type_rows_map[[tx_k]] <- ss_type_rows_k
    } else {
      has_emp_map[tx_k] <- FALSE
      type_anchor_rows_list[[tx_k]] <- integer(0)
      ss_type_rows_map[[tx_k]] <- integer(0)
      gold_rows_k <- tx_rows[(det_tier[tx_rows] == "Tier_1_Gold") & is_t1_pure[tx_rows] & !is_noise_vec[tx_rows]]
      if (any(has_efa_row[gold_rows_k])) gold_rows_k <- gold_rows_k[has_efa_row[gold_rows_k]]
      if (length(gold_rows_k) == 0L) {
        gold_rows_k <- tx_rows[(det_tier[tx_rows] == "Tier_1_Gold") & !is_noise_vec[tx_rows]]
        if (any(has_efa_row[gold_rows_k])) gold_rows_k <- gold_rows_k[has_efa_row[gold_rows_k]]
      }
      if (length(gold_rows_k) > 0L) {
        max_w_k <- max(auth_w[gold_rows_k], na.rm = TRUE)
        top_gold_k <- gold_rows_k[auth_w[gold_rows_k] >= 0.85 * max_w_k]
        n_eff_type_vec[tx_k] <- as.numeric(max(1L, length(top_gold_k)))
        mu_k <- .compute_robust_taxon_medoid(X_mat[top_gold_k, , drop = FALSE], S_inv = S_inv)
      } else {
        tx_rows_efa <- tx_rows[has_efa_row[tx_rows]]
        tx_pool_k <- if (length(tx_rows_efa) > 0L) tx_rows_efa else tx_rows
        n_eff_type_vec[tx_k] <- as.numeric(max(1L, length(tx_pool_k)))
        mu_k <- .compute_robust_taxon_medoid(X_mat[tx_pool_k, , drop = FALSE], S_inv = S_inv)
      }
      basionym_entries[[length(basionym_entries) + 1L]] <- list(
        taxon = tx_k,
        basionym = tx_k,
        is_homotypic = TRUE,
        rows = integer(0),
        pool = integer(0),
        n_eff = n_eff_type_vec[tx_k],
        max_rank = 0L,
        mu = mu_k
      )
    }
    type_centroids[k, ] <- mu_k

    # Identify per-species typical GMM morphotype cluster(s) for taxon k
    if (length(u_clusters) > 0L) {
      valid_cl_rows <- tx_rows[!is.na(gmm_cl[tx_rows]) & !is_noise_vec[tx_rows]]
      if (length(valid_cl_rows) == 0L) {
        valid_cl_rows <- tx_rows[!is.na(gmm_cl[tx_rows])]
      }
      gold_type_cl_rows <- valid_cl_rows[is_type_specimen[valid_cl_rows] | (det_tier[valid_cl_rows] == "Tier_1_Gold")]
      if (length(gold_type_cl_rows) >= 3L) {
        d2_gt <- rowSums((sweep(X_mat[gold_type_cl_rows, , drop = FALSE], 2L, mu_k, "-") %*% S_inv) *
                         sweep(X_mat[gold_type_cl_rows, , drop = FALSE], 2L, mu_k, "-"))
        inlier_gt <- gold_type_cl_rows[is_type_specimen[gold_type_cl_rows] | (is.finite(d2_gt) & d2_gt <= stats::qchisq(0.975, df = max(1L, p)))]
        if (length(inlier_gt) >= 2L) gold_type_cl_rows <- inlier_gt
      }
      ref_k <- if (length(gold_type_cl_rows) >= 2L) gold_type_cl_rows else valid_cl_rows
      typ_k <- integer(0)
      if (length(ref_k) > 0L) {
        cl_mass <- tapply(auth_w[ref_k], gmm_cl[ref_k], sum, na.rm = TRUE)
        cl_mass[is.na(cl_mass)] <- 0.0
        tot_mass <- sum(cl_mass, na.rm = TRUE)
        if (tot_mass > 0.0) {
          max_mass <- max(cl_mass, na.rm = TRUE)
          top_cls <- as.integer(names(cl_mass[cl_mass >= max_mass - 1e-9]))
          if (length(top_cls) == 1L) {
            modal_cl <- top_cls[1L]
          } else {
            d2_top <- vapply(top_cls, function(cv) {
              dv <- cl_centroids[as.character(cv), ] - mu_k
              as.numeric(t(dv) %*% S_inv %*% dv)
            }, numeric(1L))
            modal_cl <- top_cls[which.min(d2_top)]
          }
          sec_cls <- as.integer(names(cl_mass[cl_mass >= 0.30 * tot_mass]))
          typ_k <- unique(c(modal_cl, sec_cls))
        }
      }
      type_cl_vals <- as.integer(stats::na.omit(gmm_cl[ss_type_rows_k]))
      if (length(type_cl_vals) > 0L) {
        if (length(ss_type_rows_k) < length(type_rows_k)) {
          typ_k <- unique(type_cl_vals)
        } else {
          typ_k <- unique(c(type_cl_vals, typ_k))
        }
      }
      if (length(typ_k) == 0L) {
        tx_cls <- sort(unique(stats::na.omit(gmm_cl[tx_rows])))
        cand_cls <- if (length(tx_cls) > 0L) tx_cls else u_clusters
        d2_cls <- vapply(cand_cls, function(cv) {
          dv <- cl_centroids[as.character(cv), ] - mu_k
          as.numeric(t(dv) %*% S_inv %*% dv)
        }, numeric(1L))
        typ_k <- cand_cls[which.min(d2_cls)]
      }
      typical_cl_map[tx_k] <- typ_k[1L]
      typical_cl_list[[tx_k]] <- typ_k
    }
  }
  list(
    type_centroids = type_centroids,
    has_emp_map = has_emp_map,
    has_emp_homo_map = has_emp_homo_map,
    ss_type_rows_map = ss_type_rows_map,
    typical_cl_map = typical_cl_map,
    typical_cl_list = typical_cl_list,
    n_eff_type_vec = n_eff_type_vec,
    type_anchor_rows_list = type_anchor_rows_list,
    basionym_entries = basionym_entries
  )
}
