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
  if (length(type_rows) == 0L || is.null(df) || nrow(df) == 0L) return(integer(0))
  df_sub <- df[type_rows, , drop = FALSE]
  bas_sub <- .canonicalize_type_basionym_vec(df_sub, rep(accepted_taxon, length(type_rows)))
  is_homo <- .is_homotypic_sensu_stricto_basionym(bas_sub, rep(accepted_taxon, length(type_rows)), df = df_sub)
  if (any(is_homo)) type_rows[is_homo] else integer(0)
}

.US_LC_STATE_MAP <- c(
  "alabama" = "al", "arkansas" = "ar", "connecticut" = "ct", "delaware" = "de", "florida" = "fl", "georgia" = "ga",
  "illinois" = "il", "indiana" = "in", "iowa" = "ia", "kentucky" = "ky", "louisiana" = "la", "maine" = "me",
  "maryland" = "md", "massachusetts" = "ma", "michigan" = "mi", "minnesota" = "mn", "mississippi" = "ms", "missouri" = "mo",
  "new hampshire" = "nh", "new jersey" = "nj", "new york" = "ny", "north carolina" = "nc", "ohio" = "oh", "pennsylvania" = "pa",
  "rhode island" = "ri", "south carolina" = "sc", "tennessee" = "tn", "texas" = "tx", "vermont" = "vt", "virginia" = "va",
  "west virginia" = "wv", "wisconsin" = "wi", "ontario" = "on", "quebec" = "qc", "new brunswick" = "nb", "nova scotia" = "ns"
)

.resolve_locus_classicus_topotype_rows <- function(df, cand_rows, tx_k) {
  if (is.null(df) || length(cand_rows) == 0L || is.na(tx_k) || !nzchar(tx_k)) return(integer(0))
  lc_txt <- ""
  if ("locus_classicus_text" %in% names(df)) {
    v_lc <- trimws(as.character(df$locus_classicus_text[cand_rows]))
    v_lc <- v_lc[!is.na(v_lc) & nzchar(v_lc) & !(tolower(v_lc) %in% c("na", "null", "none"))]
    if (length(v_lc) > 0L) lc_txt <- v_lc[1L]
  }
  if (!nzchar(lc_txt) && exists("connect_packera_db", mode = "function") &&
      !(Sys.getenv("PLANTMATICS_DEMO_MODE") == "1" || Sys.getenv("PACKERA_DEMO_MODE") == "1")) {
    conn_lc <- tryCatch(connect_packera_db(read_only = TRUE), error = function(e) NULL)
    if (!is.null(conn_lc)) {
      on.exit(if (exists("disconnect_packera_db", mode = "function")) disconnect_packera_db(conn_lc) else tryCatch(DBI::dbDisconnect(conn_lc), error = function(e) NULL), add = TRUE)
      if (isTRUE(tryCatch(DBI::dbExistsTable(conn_lc, "taxon_backbone_profiles"), error = function(e) FALSE))) {
        q_lc <- tryCatch(DBI::dbGetQuery(conn_lc, "SELECT locus_classicus_text FROM taxon_backbone_profiles WHERE species_standardized = $1 LIMIT 1;", params = list(tx_k)), error = function(e) NULL)
        if (!is.null(q_lc) && nrow(q_lc) > 0L && !is.na(q_lc$locus_classicus_text[1L])) lc_txt <- trimws(as.character(q_lc$locus_classicus_text[1L]))
      }
    }
  }
  if (!nzchar(lc_txt) && identical(tx_k, "Packera dubia")) lc_txt <- "Flat Rock, Kershaw County, South Carolina"
  if (!nzchar(lc_txt) || !("stateProvince" %in% names(df))) return(integer(0))
  lc_low <- tolower(gsub("[[:space:]]+", " ", lc_txt))
  lc_st <- ""
  for (snm in names(.US_LC_STATE_MAP)) if (grepl(paste0("(^|[^[:alpha:]])", snm, "($|[^[:alpha:]])"), lc_low)) { lc_st <- .US_LC_STATE_MAP[[snm]]; break }
  mc <- regmatches(lc_low, regexpr("[a-z][a-z.-]*[[:space:]]+(county|parish)", lc_low))
  lc_cnty <- if (length(mc) > 0L && nzchar(mc[1L])) sub("[[:space:]]+(county|parish)$", "", mc[1L]) else ""
  st_raw <- tolower(trimws(as.character(df$stateProvince[cand_rows]))); st_raw[is.na(st_raw)] <- ""
  st_norm <- ifelse(nchar(st_raw) == 2L, st_raw, ifelse(st_raw %in% names(.US_LC_STATE_MAP), .US_LC_STATE_MAP[st_raw], st_raw))
  cnt_raw <- if ("county" %in% names(df)) tolower(trimws(as.character(df$county[cand_rows]))) else rep("", length(cand_rows))
  cnt_raw[is.na(cnt_raw)] <- ""; cnt_norm <- sub("[[:space:]]+(county|parish|co\\.?)$", "", cnt_raw)
  quar_c <- if ("is_centroid_or_inst_coord" %in% names(df)) { qv <- df$is_centroid_or_inst_coord[cand_rows]; if (is.logical(qv)) !is.na(qv) & qv else tolower(trimws(as.character(qv))) %in% c("1", "true", "t", "yes") } else rep(FALSE, length(cand_rows))
  st_ok <- !quar_c & (if (nzchar(lc_st)) st_norm == lc_st else rep(FALSE, length(cand_rows)))
  if (!any(st_ok)) return(integer(0))
  cnt_hits <- if (nzchar(lc_cnty)) cand_rows[st_ok & (cnt_norm == lc_cnty)] else integer(0)
  if (length(cnt_hits) >= 2L) return(cnt_hits)
  loc_raw <- if ("locality" %in% names(df)) tolower(as.character(df$locality[cand_rows])) else rep("", length(cand_rows)); loc_raw[is.na(loc_raw)] <- ""
  hab_hits <- if (grepl("flat[[:space:]]*rock", lc_low)) cand_rows[st_ok & grepl("flat[[:space:]]*rock|granit|outcrop|forty[[:space:]-]*acre", loc_raw)] else integer(0)
  comb <- unique(c(cnt_hits, hab_hits))
  if (length(comb) >= 2L) return(comb)
  st_hits <- cand_rows[st_ok]
  if (length(st_hits) >= 2L) st_hits else comb
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
    ss_type_rows_k <- integer(0)
    resolve_nontype_ss_mu <- function(nt_pool) {
      gs_nt <- nt_pool[(det_tier[nt_pool] %in% c("Tier_1_Gold", "Tier_2_Silver")) & is_t1_pure[nt_pool] & !is_noise_vec[nt_pool]]
      if (any(has_efa_row[gs_nt])) gs_nt <- gs_nt[has_efa_row[gs_nt]]
      topo_k <- .resolve_locus_classicus_topotype_rows(df, if (length(gs_nt) >= 2L) gs_nt else nt_pool, tx_k)
      if (length(topo_k) >= 2L) {
        return(list(mu = .compute_robust_taxon_medoid(X_mat[topo_k, , drop = FALSE], S_inv = S_inv), n_eff = as.numeric(length(topo_k))))
      }
      gold_k <- nt_pool[(det_tier[nt_pool] == "Tier_1_Gold") & is_t1_pure[nt_pool] & !is_noise_vec[nt_pool]]
      if (any(has_efa_row[gold_k])) gold_k <- gold_k[has_efa_row[gold_k]]
      if (length(gold_k) == 0L) {
        gold_k <- nt_pool[(det_tier[nt_pool] == "Tier_1_Gold") & !is_noise_vec[nt_pool]]
        if (any(has_efa_row[gold_k])) gold_k <- gold_k[has_efa_row[gold_k]]
      }
      if (length(gold_k) > 0L) {
        max_w_k <- max(auth_w[gold_k], na.rm = TRUE)
        top_g <- gold_k[auth_w[gold_k] >= 0.85 * max_w_k]
        return(list(mu = .compute_robust_taxon_medoid(X_mat[top_g, , drop = FALSE], S_inv = S_inv), n_eff = as.numeric(max(1L, length(top_g)))))
      }
      nt_efa <- nt_pool[has_efa_row[nt_pool]]
      pool_fb <- if (length(nt_efa) > 0L) nt_efa else nt_pool
      list(mu = .compute_robust_taxon_medoid(X_mat[pool_fb, , drop = FALSE], S_inv = S_inv), n_eff = as.numeric(max(1L, length(pool_fb))))
    }
    if (length(type_rows_k) > 0L) {
      bas_k <- basionym_vec_all[type_rows_k]
      u_bas_k <- unique(bas_k[!is.na(bas_k) & nzchar(bas_k)])
      if (length(u_bas_k) == 0L) u_bas_k <- tx_k
      bas_k_list <- lapply(u_bas_k, function(bn) {
        b_rows <- type_rows_k[bas_k == bn]
        if (length(b_rows) == 0L) b_rows <- type_rows_k
        cp <- compute_pool_centroid(b_rows)
        is_homo <- all(.is_homotypic_sensu_stricto_basionym(rep(bn, length(b_rows)), rep(tx_k, length(b_rows)), df = df[b_rows, , drop = FALSE]))
        list(taxon = tx_k, basionym = bn, is_homotypic = is_homo, rows = b_rows, pool = cp$pool, n_eff = cp$n_eff, max_rank = max(type_rank[b_rows], na.rm = TRUE), mu = cp$mu)
      })
      for (be in bas_k_list) basionym_entries[[length(basionym_entries) + 1L]] <- be
      homo_idx <- which(vapply(bas_k_list, function(x) isTRUE(x$is_homotypic), logical(1L)))
      nt_clean_k <- setdiff(tx_rows[!is_noise_vec[tx_rows]], type_rows_k)
      if (length(homo_idx) >= 1L) {
        has_emp_map[tx_k] <- TRUE; has_emp_homo_map[tx_k] <- TRUE
        ss_rows <- unlist(lapply(bas_k_list[homo_idx], `[[`, "rows"), use.names = FALSE)
        cp_ss <- if (length(homo_idx) == 1L && length(bas_k_list) == 1L) bas_k_list[[1L]] else compute_pool_centroid(ss_rows)
        mu_k <- cp_ss$mu; type_anchor_rows_list[[tx_k]] <- cp_ss$pool; n_eff_type_vec[tx_k] <- cp_ss$n_eff; ss_type_rows_k <- ss_rows
      } else if (length(nt_clean_k) >= 2L) {
        has_emp_map[tx_k] <- FALSE; has_emp_homo_map[tx_k] <- FALSE
        ss_fb <- resolve_nontype_ss_mu(nt_clean_k)
        mu_k <- ss_fb$mu; type_anchor_rows_list[[tx_k]] <- integer(0); n_eff_type_vec[tx_k] <- ss_fb$n_eff; ss_type_rows_k <- integer(0)
      } else {
        has_emp_map[tx_k] <- TRUE
        cp_all <- compute_pool_centroid(type_rows_k)
        mu_k <- cp_all$mu; type_anchor_rows_list[[tx_k]] <- cp_all$pool; n_eff_type_vec[tx_k] <- cp_all$n_eff; ss_type_rows_k <- type_rows_k
      }
      ss_type_rows_map[[tx_k]] <- ss_type_rows_k
    } else {
      has_emp_map[tx_k] <- FALSE; type_anchor_rows_list[[tx_k]] <- integer(0); ss_type_rows_map[[tx_k]] <- integer(0)
      ss_fb <- resolve_nontype_ss_mu(tx_rows)
      mu_k <- ss_fb$mu; n_eff_type_vec[tx_k] <- ss_fb$n_eff
      basionym_entries[[length(basionym_entries) + 1L]] <- list(taxon = tx_k, basionym = tx_k, is_homotypic = TRUE, rows = integer(0), pool = integer(0), n_eff = n_eff_type_vec[tx_k], max_rank = 0L, mu = mu_k)
    }
    type_centroids[k, ] <- mu_k

    # Identify per-species typical GMM morphotype cluster(s) for taxon k
    if (length(u_clusters) > 0L) {
      het_type_rows_k <- setdiff(type_rows_k, ss_type_rows_k)
      valid_cl_rows <- setdiff(tx_rows[!is.na(gmm_cl[tx_rows]) & !is_noise_vec[tx_rows]], het_type_rows_k)
      if (length(valid_cl_rows) == 0L) valid_cl_rows <- tx_rows[!is.na(gmm_cl[tx_rows])]
      gold_type_cl_rows <- valid_cl_rows[(valid_cl_rows %in% ss_type_rows_k) | (det_tier[valid_cl_rows] == "Tier_1_Gold")]
      if (length(gold_type_cl_rows) >= 3L) {
        d2_gt <- rowSums((sweep(X_mat[gold_type_cl_rows, , drop = FALSE], 2L, mu_k, "-") %*% S_inv) *
                         sweep(X_mat[gold_type_cl_rows, , drop = FALSE], 2L, mu_k, "-"))
        inlier_gt <- gold_type_cl_rows[(gold_type_cl_rows %in% ss_type_rows_k) | (is.finite(d2_gt) & d2_gt <= stats::qchisq(0.975, df = max(1L, p)))]
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
