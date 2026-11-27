# Step 05 Partitioned PCA-env Niche Tests, Weighted Edaphic Z-Scores & 3-Way Cross-Modal Consensus
if (grDevices::dev.cur() == 1L) grDevices::pdf(NULL)

.s05_bcol <- function(df, nm, def = FALSE) {
  n <- nrow(df); if (!(nm %in% names(df))) return(rep(def, n))
  b <- as.logical(df[[nm]]); b[is.na(b)] <- def; b
}
.s05_screen_point_mass_mat <- function(X_m, max_tie_frac = 0.25, min_keep = 2L) {
  if (is.null(X_m) || !is.matrix(X_m) || ncol(X_m) < 3L || !exists("screen_edaphic_point_mass_columns", mode = "function")) return(X_m)
  if (is.null(colnames(X_m))) colnames(X_m) <- paste0("V", seq_len(ncol(X_m)))
  res <- tryCatch(screen_edaphic_point_mass_columns(X_m, max_tie_frac = max_tie_frac), error = function(e) NULL); if (is.null(res)) return(X_m)
  sel <- if (is.matrix(res) || is.data.frame(res)) as.matrix(res) else if (is.list(res)) { kh <- intersect(c("kept_cols", "kept_vars", "selected_variables", "selected_vars", "keep_cols", "keep"), names(res)); if (length(kh) > 0L) res[[kh[1L]]] else if (is.matrix(res$X) || is.data.frame(res$X)) as.matrix(res$X) else NULL } else res
  out <- if (is.matrix(sel)) sel else if (is.character(sel)) X_m[, intersect(sel, colnames(X_m)), drop = FALSE] else if (is.logical(sel) && length(sel) == ncol(X_m)) X_m[, !is.na(sel) & sel, drop = FALSE] else if (is.numeric(sel)) X_m[, as.integer(sel[is.finite(sel) & sel >= 1L & sel <= ncol(X_m)]), drop = FALSE] else NULL
  if (!is.null(out) && ncol(out) >= min(as.integer(min_keep), ncol(X_m)) && nrow(out) == nrow(X_m)) out else X_m
}
.s05_resolve_anchor_flags <- function(df, grp_vec, active_taxa, is_intra = FALSE) {
  n <- nrow(df)
  fe <- if (exists("is_fine_edaphic_type_or_anchor_eligible", mode = "function")) is_fine_edaphic_type_or_anchor_eligible(df, 10000.0, 50000.0) else list(edaphic_ok = rep(TRUE, n), clim_ok = rep(TRUE, n), is_open_hyb = rep(FALSE, n), is_quarantined_type = rep(FALSE, n))
  elig <- .s05_bcol(df, "is_training_eligible_coord", TRUE) & fe$edaphic_ok
  is_t <- if (exists(".identify_type_specimen_vec", mode = "function")) as.logical(.identify_type_specimen_vec(df)$is_type_specimen) else .s05_bcol(df, "is_type_specimen", FALSE)
  is_p <- .s05_bcol(df, "is_promoted_bronze", FALSE); cl_n <- .s05_bcol(df, "is_confident_label_noise", FALSE)
  det_t <- if ("determiner_tier" %in% names(df)) ifelse(is.na(df$determiner_tier) | !nzchar(trimws(as.character(df$determiner_tier))), "Tier_3_Bronze", as.character(df$determiner_tier)) else rep("Tier_1_Gold", n)
  is_ht <- .s05_bcol(df, "is_heterotypic_type_audit", FALSE)
  if ("anchor_tier_source" %in% names(df)) is_ht <- is_ht | (!is.na(df$anchor_tier_source) & as.character(df$anchor_tier_source) == "Passive_Heterotypic_Type_Audit")
  is_tc <- is_t & !fe$is_quarantined_type & !is_ht
  tok <- (grepl("Tier_1|Gold|Tier_2|Silver", det_t, ignore.case = TRUE) | is_p | is_tc) & (!cl_n | is_tc) & !fe$is_open_hyb & !is_ht & !fe$is_quarantined_type
  if (exists(".is_step05_extralimital_non_type", mode = "function")) tok <- tok & !.is_step05_extralimital_non_type(df, grp_vec, tok & elig, is_tc)
  mg <- if (isTRUE(is_intra)) .s05_bcol(df, "is_missing_gmm_cluster", FALSE) else rep(FALSE, n)
  anc <- tok & elig & (grp_vec %in% active_taxa) & !mg
  if (sum(anc) < min(4L, n)) anc <- tok & elig & (grp_vec %in% active_taxa)
  if (!any(tok)) { anc <- !cl_n & !fe$is_open_hyb & !fe$is_quarantined_type & elig & (grp_vec %in% active_taxa); if (sum(anc) < min(2L, n)) anc <- !cl_n & elig & (grp_vec %in% active_taxa); if (sum(anc) == 0L) anc <- grp_vec %in% active_taxa }
  list(is_anchor = anc, tier_ok = tok, elig_coord = elig, is_type = is_t, is_type_clean = is_tc, cl_noise = cl_n, fine_elig = fe)
}
.s05_w_mean <- function(x, w) { ok <- is.finite(x) & is.finite(w) & (w > 0.0); if (!any(ok)) mean(as.numeric(x), na.rm = TRUE) else sum(x[ok] * w[ok]) / sum(w[ok]) }
.s05_w_sd <- function(x, w, ms = 0.25) {
  ok <- is.finite(x) & is.finite(w) & (w > 0.0); if (sum(ok) < 2L) return(max(stats::sd(as.numeric(x), na.rm = TRUE), ms, na.rm = TRUE))
  wm <- sum(x[ok] * w[ok]) / sum(w[ok]); dn <- sum(w[ok]) - (sum(w[ok]^2) / sum(w[ok]))
  max(sqrt(max(if (dn > 1e-6) sum(w[ok] * (x[ok] - wm)^2) / dn else stats::var(as.numeric(x), na.rm = TRUE), 0.0, na.rm = TRUE)), ms, na.rm = TRUE)
}
.s05_spatial_wts <- function(df, active_mask = NULL) {
  if (!exists("compute_spatial_balance_weights", mode = "function")) return(rep(1.0, nrow(df)))
  d_sb <- df; d_sb$authority_weight <- NULL; w <- if (!is.null(active_mask) && all(c("active_mask", "recompute") %in% names(formals(compute_spatial_balance_weights)))) compute_spatial_balance_weights(d_sb, active_mask = active_mask, recompute = TRUE) else compute_spatial_balance_weights(d_sb)
  w[!is.finite(w) | w <= 0.0] <- 1.0; w
}
.s05_scale_floors <- function(df, cols, anc_mask) stats::setNames(vapply(cols, function(cl) {
  va <- suppressWarnings(as.numeric(df[[cl]])); vf <- va[anc_mask]; vf <- vf[is.finite(vf)]; if (length(vf) < 2L) vf <- va[is.finite(va)]
  if (length(vf) < 2L) return(1e-3); dg <- max(c(stats::sd(vf), stats::IQR(vf) / 1.349), na.rm = TRUE); if (!is.finite(dg) || dg <= 1e-8) 1e-3 else max(0.15 * dg, 1e-4)
}, numeric(1)), cols)

.prepare_partition_scores_glob <- function(bg_mat, rf_var_weights = NULL, coords_occ = NULL) {
  if (is.null(bg_mat) || ncol(bg_mat) == 0L || nrow(bg_mat) < 6L) return(NULL)
  n_occ <- nrow(bg_mat)
  for (j in seq_len(ncol(bg_mat))) { cv <- bg_mat[, j]; if (any(!is.finite(cv))) { mj <- stats::median(cv[is.finite(cv)], na.rm = TRUE); bg_mat[!is.finite(cv), j] <- if (is.finite(mj)) mj else 0.0 } }
  act <- which(apply(bg_mat, 2L, stats::sd, na.rm = TRUE) > 1e-8); if (length(act) == 0L) return(NULL)
  bg_sub <- bg_mat[, act, drop = FALSE]
  cont <- which(apply(bg_sub, 2L, function(v) length(unique(round(v, 5L))) >= 4L)); if (length(cont) >= 1L) bg_sub <- bg_sub[, cont, drop = FALSE]
  bg_sub <- .s05_screen_point_mass_mat(bg_sub, max_tie_frac = 0.25, min_keep = 2L)
  has_c <- !is.null(coords_occ) && nrow(coords_occ) == n_occ && all(is.finite(as.numeric(coords_occ[, 1L]))) && all(is.finite(as.numeric(coords_occ[, 2L])))
  bg_ex <- NULL; bg_ec <- NULL
  if (has_c && n_occ >= 8L && exists(".generate_accessible_background_env", mode = "function")) {
    bg_g <- tryCatch(.generate_accessible_background_env(coords_occ, bg_sub, n_bg = min(500L, max(150L, 3L * n_occ)), seed = 108L), error = function(e) NULL)
    if (!is.null(bg_g) && nrow(bg_g$X_bg) >= 20L) { bg_ex <- bg_g$X_bg; bg_ec <- bg_g$coords }
  }
  comb_X <- if (!is.null(bg_ex)) rbind(bg_sub, bg_ex) else bg_sub
  sc <- tryCatch(if (ncol(comb_X) >= 2L) as.data.frame(stats::prcomp(comb_X, center = TRUE, scale. = TRUE)$x[, 1:2, drop = FALSE]) else { z1 <- as.numeric(scale(comb_X[, 1L])); data.frame(PC1 = z1, PC2 = 0.05 * sin(seq_along(z1) * 0.37)) }, error = function(e) NULL)
  if (is.null(sc)) return(NULL)
  list(bg_sub = bg_sub, scores_glob = sc[seq_len(n_occ), 1:2, drop = FALSE], scores_M_all = sc, scores_M_bg = if (!is.null(bg_ex)) sc[(n_occ + 1L):nrow(sc), 1:2, drop = FALSE] else NULL, coords_occ = if (has_c) as.matrix(coords_occ) else NULL, coords_M_bg = bg_ec)
}
.resolve_lineage_accessible_glob1 <- function(prep, idx_grp) {
  sc <- prep$scores_glob; if (is.null(prep$scores_M_bg) || is.null(prep$coords_occ) || is.null(prep$coords_M_bg)) return(sc)
  c_sp <- prep$coords_occ[idx_grp, , drop = FALSE]; c_bg <- prep$coords_M_bg; if (nrow(c_sp) == 0L || nrow(c_bg) == 0L) return(sc)
  md <- apply(compute_haversine_dist_km_matrix(c_bg, c_sp), 1L, min, na.rm = TRUE)
  hit <- which(md <= 125.0); if (length(hit) < 15L) hit <- order(md)[seq_len(min(30L, length(md)))]
  rbind(sc[idx_grp, 1:2, drop = FALSE], prep$scores_M_bg[hit, 1:2, drop = FALSE])
}
.eval_bivnorm_kde_grid <- function(xm, ext, dx, dy, R, gs, gr, em = NULL) {
  if (is.null(xm) || nrow(xm) < 2L) return(NULL)
  f <- is.finite(xm[, 1L]) & is.finite(xm[, 2L]); if (sum(f) < 2L) return(NULL); if (!all(f)) xm <- xm[f, , drop = FALSE]
  n <- nrow(xm); x1 <- (xm[, 1L] - ext[1L]) / dx; x2 <- (xm[, 2L] - ext[3L]) / dy
  vx <- stats::var(x1); vy <- stats::var(x2); h <- sqrt(0.5 * ((if (is.finite(vx)) vx else 0.0) + (if (is.finite(vy)) vy else 0.0))) * (n^(-1.0 / 6.0))
  if (!is.finite(h) || h <= 1e-6) return(NULL)
  ih2 <- -0.5 / (h * h); d <- tcrossprod(exp(ih2 * outer(gr, x2, "-")^2), exp(ih2 * outer(gs, x1, "-")^2)) / (n * 2.0 * pi * h * h)
  ci <- pmin(R, pmax(1L, as.integer(floor(x1 * R)) + 1L)); ri <- pmin(R, pmax(1L, as.integer(floor((1.0 - x2) * R)) + 1L))
  ci[xm[, 1L] >= ext[2L]] <- R; ri[xm[, 2L] <= ext[3L]] <- R; th <- suppressWarnings(min(d[cbind(ri, ci)], na.rm = TRUE))
  if (is.finite(th)) d[d < th] <- 0.0; if (!is.null(em)) d[!em] <- 0.0; d
}
.build_fast_pca_env_taxon_grid <- function(prep, idx_sp, R = 100L) {
  if (is.null(prep) || is.null(prep$scores_glob) || length(idx_sp) < 5L) return(NULL)
  sco <- as.matrix(prep$scores_glob[, 1:2, drop = FALSE]); scg <- as.matrix((if (!is.null(prep$scores_M_all)) prep$scores_M_all else prep$scores_glob)[, 1:2, drop = FALSE])
  ext <- c(min(scg[, 1L], na.rm = TRUE), max(scg[, 1L], na.rm = TRUE), min(scg[, 2L], na.rm = TRUE), max(scg[, 2L], na.rm = TRUE))
  dx <- abs(ext[2L] - ext[1L]); dy <- abs(ext[4L] - ext[3L]); if (!is.finite(dx) || dx <= 1e-12 || !is.finite(dy) || dy <= 1e-12) return(NULL)
  R <- as.integer(R); gs <- (0:(R - 1L)) / (R - 2L); gr <- gs[R:1L]
  g1m <- as.matrix(.resolve_lineage_accessible_glob1(prep, idx_sp)[, 1:2, drop = FALSE]); spm <- sco[idx_sp, 1:2, drop = FALSE]
  g1d <- .eval_bivnorm_kde_grid(g1m, ext, dx, dy, R, gs, gr, NULL); if (is.null(g1d) || !is.finite(sum(g1d)) || sum(g1d) <= 0.0) return(NULL)
  em <- g1d > 0.0; spd <- .eval_bivnorm_kde_grid(spm, ext, dx, dy, R, gs, gr, em); if (is.null(spd) || !is.finite(sum(spd)) || sum(spd) <= 0.0) return(NULL)
  Z <- g1d * (nrow(g1m) / sum(g1d)); z <- spd * (nrow(spm) / sum(spd)); if (!is.finite(max(z)) || max(z) <= 0.0) return(NULL)
  zu <- z / max(z); zc <- z / Z; zc[!is.finite(zc)] <- 0.0; if (!is.finite(max(zc)) || max(zc) <= 0.0) return(NULL)
  zc <- zc / max(zc); pc <- zc / sum(zc); crc <- which(Z > 0.0, arr.ind = TRUE)
  list(R = R, ext = ext, dx = dx, dy = dy, g_seq = gs, g_rev = gr, sp = spm, z = z, z.uncor = zu, z.cor = zc, p_cor = pc, Z = Z, env_mask = em,
       cent_rc = which(zu == max(zu), arr.ind = TRUE)[1L, , drop = TRUE], cand_rc = crc, cand_w = if (nrow(crc) > 0L) Z[crc] / max(Z) else numeric(0))
}
.run_fast_pca_env_pair_overlap <- function(g1, g2, n_perm) {
  p1 <- g1$p_cor; p2 <- g2$p_cor; sq1 <- sqrt(p1); De <- 1.0 - 0.5 * sum(abs(p1 - p2)); Ie <- 1.0 - 0.5 * sum((sq1 - sqrt(p2))^2)
  if (!is.finite(De) || !is.finite(Ie)) return(list(ok = FALSE))
  if (n_perm <= 0L) return(list(D_emp = De, I_emp = Ie, p_equiv_D = NA_real_, p_equiv_I = NA_real_, p_bg_D = NA_real_, p_bg_I = NA_real_, D_null = numeric(0), I_null = numeric(0), ok = TRUE))
  op <- rbind(g1$sp, g2$sp); n1 <- nrow(g1$sp); nt <- nrow(op)
  Z1i <- ifelse(g1$env_mask & (g1$Z > 0.0), 1.0 / g1$Z, 0.0); Z2i <- ifelse(g2$env_mask & (g2$Z > 0.0), 1.0 / g2$Z, 0.0)
  Dn <- rep(NA_real_, n_perm); In <- rep(NA_real_, n_perm)
  for (b in seq_len(n_perm)) {
    ip <- sample.int(nt, n1); d1 <- .eval_bivnorm_kde_grid(op[ip, , drop = FALSE], g1$ext, g1$dx, g1$dy, g1$R, g1$g_seq, g1$g_rev, g1$env_mask)
    d2 <- .eval_bivnorm_kde_grid(op[-ip, , drop = FALSE], g2$ext, g2$dx, g2$dy, g2$R, g2$g_seq, g2$g_rev, g2$env_mask); if (is.null(d1) || is.null(d2)) next
    z1 <- d1 * Z1i; z2 <- d2 * Z2i; s1 <- sum(z1); s2 <- sum(z2); if (!is.finite(s1) || s1 <= 0.0 || !is.finite(s2) || s2 <= 0.0) next
    p1s <- z1 / s1; p2s <- z2 / s2; Dn[b] <- 1.0 - 0.5 * sum(abs(p1s - p2s)); In[b] <- 1.0 - 0.5 * sum((sqrt(p1s) - sqrt(p2s))^2)
  }
  dfe <- Dn[is.finite(Dn)]; ife <- In[is.finite(In)]; pbD <- NA_real_; pbI <- NA_real_
  if (nrow(g2$cand_rc) > 0L && length(g2$cent_rc) == 2L) {
    R <- g2$R; cr <- as.integer(g2$cent_rc[1L]); cc <- as.integer(g2$cent_rc[2L])
    pks <- sample.int(nrow(g2$cand_rc), size = n_perm, replace = TRUE, prob = g2$cand_w); Db <- rep(NA_real_, n_perm); Ib <- rep(NA_real_, n_perm)
    for (b in seq_len(n_perm)) {
      dr <- as.integer(g2$cand_rc[pks[b], 1L]) - cr; dc <- as.integer(g2$cand_rc[pks[b], 2L]) - cc
      rs1 <- max(1L, 1L - dr); rs2 <- min(R, R - dr); cs1 <- max(1L, 1L - dc); cs2 <- min(R, R - dc); if (rs1 > rs2 || cs1 > cs2) next
      zs <- matrix(0.0, R, R); zs[(rs1 + dr):(rs2 + dr), (cs1 + dc):(cs2 + dc)] <- g2$z[rs1:rs2, cs1:cs2, drop = FALSE]
      zcs <- zs * Z2i; sb <- sum(zcs); if (!is.finite(sb) || sb <= 0.0) next
      p2b <- zcs / sb; Db[b] <- 1.0 - 0.5 * sum(abs(p1 - p2b)); Ib[b] <- 1.0 - 0.5 * sum((sq1 - sqrt(p2b))^2)
    }
    dfb <- Db[is.finite(Db)]; ifb <- Ib[is.finite(Ib)]
    pbD <- if (length(dfb) > 0L) (sum(dfb <= De) + 1.0) / (length(dfb) + 1.0) else NA_real_; pbI <- if (length(ifb) > 0L) (sum(ifb <= Ie) + 1.0) / (length(ifb) + 1.0) else NA_real_
  }
  list(D_emp = De, I_emp = Ie, p_equiv_D = if (length(dfe) > 0L) (sum(dfe <= De) + 1.0) / (length(dfe) + 1.0) else NA_real_,
       p_equiv_I = if (length(ife) > 0L) (sum(ife <= Ie) + 1.0) / (length(ife) + 1.0) else NA_real_, p_bg_D = pbD, p_bg_I = pbI, D_null = Dn, I_null = In, ok = TRUE)
}
.evaluate_partition_jobs_parallel <- function(job_specs, grp_col, n_perm) {
  nj <- length(job_specs); if (nj == 0L) return(list())
  np <- if (is.null(n_perm) || length(n_perm) == 0L || is.na(n_perm[1L])) 0L else as.integer(n_perm[1L])
  jseeds <- if (np > 0L) sample.int(2147483647L, size = nj, replace = TRUE) else rep(NA_integer_, nj)
  pmap <- list(); sprep <- list(); jpks <- character(nj); gspecs <- list()
  for (ji in seq_len(nj)) {
    jb <- job_specs[[ji]]; pk <- if (!is.null(jb$part_id) && nzchar(as.character(jb$part_id[1L]))) as.character(jb$part_id[1L]) else NULL
    pj <- if (!is.null(jb$prep)) jb$prep else .prepare_partition_scores_glob(jb$bg_mat, jb$rf_var_weights)
    if (is.null(pk)) { hk <- NULL; for (sk in names(sprep)) if (identical(sprep[[sk]], pj)) { hk <- sk; break }; if (is.null(hk)) { hk <- sprintf("part_%d", length(sprep) + 1L); sprep[[hk]] <- pj }; pk <- hk }
    if (is.null(pmap[[pk]])) pmap[[pk]] <- pj; jpks[ji] <- pk
    for (tx in c(jb$sp1, jb$sp2)) { gk <- paste0(pk, "::", tx); if (!(gk %in% names(gspecs))) gspecs[[gk]] <- list(part_key = pk, taxon = as.character(tx)) }
  }
  mw <- if (exists("resolve_packera_workers", mode = "function")) resolve_packera_workers(64L) else 1L; up <- (mw > 1L) && exists("packera_mclapply", mode = "function")
  bgf <- function(gs) { ix <- which(grp_col == gs$taxon); if (length(ix) < 5L) NULL else .build_fast_pca_env_taxon_grid(pmap[[gs$part_key]], ix, 100L) }
  gl <- if (up && length(gspecs) >= 4L) packera_mclapply(gspecs, bgf, mc.cores = resolve_packera_workers(64L, length(gspecs))) else lapply(gspecs, bgf); names(gl) <- names(gspecs)
  ej <- function(ji) { jb <- job_specs[[ji]]; pk <- jpks[ji]; compute_pca_env_partition_overlap(jb$bg_mat, grp_col, jb$sp1, jb$sp2, np, jb$rf_var_weights, pmap[[pk]], gl[[paste0(pk, "::", jb$sp1)]], gl[[paste0(pk, "::", jb$sp2)]], jseeds[ji]) }
  if (up && nj >= 2L) packera_mclapply(seq_len(nj), ej, mc.cores = resolve_packera_workers(64L, nj), mc.preschedule = FALSE) else lapply(seq_len(nj), ej)
}

compute_pca_env_partition_overlap <- function(bg_mat, grp_col, sp1, sp2, n_perm = 100, rf_var_weights = NULL, precomputed_prep = NULL, precomputed_g1 = NULL, precomputed_g2 = NULL, job_seed = NULL) {
  if (!is.null(job_seed) && length(job_seed) >= 1L && is.finite(job_seed[1L])) {
    hs <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE); os <- if (hs) get(".Random.seed", envir = .GlobalEnv, inherits = FALSE) else NULL
    on.exit({ if (hs && !is.null(os)) tryCatch(assign(".Random.seed", os, envir = .GlobalEnv), error = function(e) NULL) else if (!hs && exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) tryCatch(rm(".Random.seed", envir = .GlobalEnv), error = function(e) NULL) }, add = TRUE)
    set.seed(as.integer(job_seed[1L]))
  }
  np <- if (is.null(n_perm) || length(n_perm) == 0L || is.na(n_perm[1L])) 0L else as.integer(n_perm[1L])
  prep <- if (!is.null(precomputed_prep)) precomputed_prep else .prepare_partition_scores_glob(bg_mat, rf_var_weights)
  if (is.null(prep)) return(list(D_emp = NA_real_, I_emp = NA_real_, p_equiv_D = NA_real_, p_equiv_I = NA_real_, p_bg_D = NA_real_, p_bg_I = NA_real_, D_null = rep(NA_real_, max(0L, np)), I_null = rep(NA_real_, max(0L, np))))
  i1 <- which(grp_col == sp1); i2 <- which(grp_col == sp2); n1 <- length(i1); n2 <- length(i2)
  if (!is.null(prep$scores_glob) && n1 >= 5L && n2 >= 5L) {
    g1 <- if (!is.null(precomputed_g1)) precomputed_g1 else .build_fast_pca_env_taxon_grid(prep, i1, 100L)
    g2 <- if (!is.null(precomputed_g2)) precomputed_g2 else .build_fast_pca_env_taxon_grid(prep, i2, 100L)
    if (!is.null(g1) && !is.null(g2)) { ff <- tryCatch(.run_fast_pca_env_pair_overlap(g1, g2, np), error = function(e) list(ok = FALSE)); if (isTRUE(ff$ok) && (np <= 0L || (is.finite(ff$p_equiv_D) && is.finite(ff$p_bg_D)))) return(ff[c("D_emp", "I_emp", "p_equiv_D", "p_equiv_I", "p_bg_D", "p_bg_I", "D_null", "I_null")]) }
  }
  obg <- if (!is.null(prep$scores_glob)) as.matrix(prep$scores_glob) else as.matrix(prep$bg_sub)
  gp <- function(sm) { m <- colMeans(sm); s <- apply(sm, 2L, stats::sd) + 1e-4; p <- exp(-0.5 * rowSums(t((t(obg) - m) / s)^2)); p / max(sum(p), 1e-12) }
  p1 <- gp(obg[i1, , drop = FALSE]); p2 <- gp(obg[i2, , drop = FALSE]); De <- 1.0 - 0.5 * sum(abs(p1 - p2)); Ie <- 1.0 - 0.5 * sqrt(sum((sqrt(p1) - sqrt(p2))^2))
  if (np <= 0L) return(list(D_emp = De, I_emp = Ie, p_equiv_D = NA_real_, p_equiv_I = NA_real_, p_bg_D = NA_real_, p_bg_I = NA_real_, D_null = numeric(0), I_null = numeric(0)))
  pl <- rbind(obg[i1, , drop = FALSE], obg[i2, , drop = FALSE]); Dn <- numeric(np); In <- numeric(np); Db <- numeric(np); Ib <- numeric(np)
  for (k in seq_len(np)) {
    ip <- sample.int(n1 + n2, n1); pp1 <- gp(pl[ip, , drop = FALSE]); pp2 <- gp(pl[-ip, , drop = FALSE])
    Dn[k] <- 1.0 - 0.5 * sum(abs(pp1 - pp2)); In[k] <- 1.0 - 0.5 * sqrt(sum((sqrt(pp1) - sqrt(pp2))^2))
    bp2 <- gp(obg[sample.int(nrow(obg), min(n2, nrow(obg)), replace = TRUE), , drop = FALSE]); Db[k] <- 1.0 - 0.5 * sum(abs(p1 - bp2)); Ib[k] <- 1.0 - 0.5 * sqrt(sum((sqrt(p1) - sqrt(bp2))^2))
  }
  list(D_emp = De, I_emp = Ie, p_equiv_D = (sum(Dn <= De) + 1.0) / (np + 1.0), p_equiv_I = (sum(In <= Ie) + 1.0) / (np + 1.0), p_bg_D = (sum(Db <= De) + 1.0) / (np + 1.0), p_bg_I = (sum(Ib <= Ie) + 1.0) / (np + 1.0), D_null = Dn, I_null = In)
}
.prepare_niche_calib_df <- function(df, intraspecific_mode, primary_vars = NULL) {
  if (!("analysis_group" %in% names(df)) || all(is.na(df$analysis_group) | !nzchar(trimws(as.character(df$analysis_group))))) df <- resolve_analysis_group(df, intraspecific_mode = intraspecific_mode, primary_env_vars = primary_vars, rf_var_weights = NULL)
  fe <- if (exists("is_fine_edaphic_type_or_anchor_eligible", mode = "function")) is_fine_edaphic_type_or_anchor_eligible(df, 10000.0, 50000.0) else list(edaphic_ok = rep(TRUE, nrow(df)), is_open_hyb = rep(FALSE, nrow(df)), is_quarantined_type = rep(FALSE, nrow(df)))
  el <- .s05_bcol(df, "is_training_eligible_coord", TRUE) & fe$edaphic_ok & !fe$is_open_hyb & !fe$is_quarantined_type
  am <- el & !is.na(resolve_active_anchor_groups(df, intraspecific_mode = intraspecific_mode))
  dc <- if (sum(am) >= 10L) df[am, , drop = FALSE] else if (sum(el) >= 10L) df[el, , drop = FALSE] else df
  ii <- isTRUE(intraspecific_mode) || (length(unique(stats::na.omit(dc$species_standardized))) < 2L); gc <- if (ii) dc$analysis_group else dc$species_standardized
  cg <- if (ii) sort(unique(stats::na.omit(gc))) else get_active_packera_taxa(gc); if (length(cg) < 2L) cg <- sort(unique(stats::na.omit(dc$analysis_group)))
  list(df = df, df_calib = dc, is_intra = ii, grp_col = gc, cand_groups = cg, pairs = if (length(cg) >= 2L) utils::combn(cg, 2L, simplify = FALSE) else list())
}

compute_edaphic_schoener_d <- function(df, primary_edaph_vars = NULL, rf_var_weights = NULL, n_perm = 100, intraspecific_mode = FALSE) {
  cal <- .prepare_niche_calib_df(df, intraspecific_mode, primary_edaph_vars)
  ec <- setdiff(intersect(as.character(primary_edaph_vars), names(cal$df_calib)), NON_ENM_RAW_GEOCHEM_HYDRO_COLS)
  if (length(ec) == 0L) ec <- grep("^(fine_|ilr_texture_|ssurgo_|solus_|soil_|sgmc_|is_azonal|lithic_outcrop_|wetland_seepage_)", discover_candidate_eco_predictors(cal$df_calib, FALSE, 3L), value = TRUE, ignore.case = TRUE)
  emp <- data.frame(Pair = character(0), Taxon1 = character(0), Taxon2 = character(0), D_edaph = numeric(0), I_edaph = numeric(0), p_equiv_D_edaph = numeric(0), p_bg_D_edaph = numeric(0), stringsAsFactors = FALSE)
  if (length(ec) == 0L || length(cal$pairs) == 0L) return(emp)
  Xe <- as.matrix(cal$df_calib[, ec, drop = FALSE]); pe <- .prepare_partition_scores_glob(Xe, rf_var_weights); jobs <- list(); meta <- list()
  for (pr in cal$pairs) if (sum(cal$grp_col == pr[1L], na.rm = TRUE) >= 5L && sum(cal$grp_col == pr[2L], na.rm = TRUE) >= 5L) {
    jobs[[length(jobs) + 1L]] <- list(part_id = "edaph", bg_mat = Xe, sp1 = pr[1L], sp2 = pr[2L], rf_var_weights = rf_var_weights, prep = pe)
    meta[[length(meta) + 1L]] <- list(pn = paste0(gsub("Packera ", "P.", pr[1L]), " vs ", gsub("Packera ", "P.", pr[2L])), sp1 = pr[1L], sp2 = pr[2L])
  }
  if (length(jobs) == 0L) return(emp)
  outs <- .evaluate_partition_jobs_parallel(jobs, cal$grp_col, n_perm)
  do.call(rbind, lapply(seq_along(jobs), function(k) data.frame(Pair = meta[[k]]$pn, Taxon1 = meta[[k]]$sp1, Taxon2 = meta[[k]]$sp2, D_edaph = round(outs[[k]]$D_emp, 4), I_edaph = round(outs[[k]]$I_emp, 4), p_equiv_D_edaph = round(outs[[k]]$p_equiv_D, 4), p_bg_D_edaph = round(outs[[k]]$p_bg_D, 4), stringsAsFactors = FALSE)))
}

run_warren_niche_identity_tests <- function(df, env_vars = NULL, n_perm = 100, intraspecific_mode = FALSE, primary_env_vars = NULL, primary_clim_vars = NULL, primary_edaph_vars = NULL, rf_var_weights = NULL) {
  message(sprintf("Executing Partitioned PCA-env Niche Equivalency & Background Similarity Tests (%d permutations per pair)...", n_perm))
  cal <- .prepare_niche_calib_df(df, intraspecific_mode, primary_env_vars); dc <- cal$df_calib; gc <- cal$grp_col; cg <- cal$cand_groups; pairs <- cal$pairs
  psw <- is_post_selection_ecotype_groups(cg, cal$is_intra, cal$df)
  if (is.null(env_vars) || length(env_vars) == 0L) env_vars <- discover_candidate_eco_predictors(dc, FALSE, 3L)
  env_vars <- intersect(setdiff(unique(c(env_vars, primary_env_vars, primary_clim_vars, primary_edaph_vars)), NON_ENM_RAW_GEOCHEM_HYDRO_COLS), names(dc))
  Xr <- as.matrix(dc[, env_vars, drop = FALSE])
  hov <- any(vapply(list(primary_env_vars, primary_clim_vars, primary_edaph_vars), function(v) !is.null(v) && length(intersect(v, colnames(Xr))) > 0L, logical(1)))
  if (hov) {
    cv <- intersect(as.character(primary_clim_vars), colnames(Xr)); ev <- intersect(as.character(primary_edaph_vars), colnames(Xr)); jv <- intersect(as.character(primary_env_vars), colnames(Xr))
    if (length(cv) == 0L) cv <- grep("^(bio[0-9]+|delta_doy|doy)", jv, value = TRUE, ignore.case = TRUE)
    if (length(ev) == 0L) ev <- setdiff(jv, cv); if (length(jv) == 0L) jv <- unique(c(cv, ev)); bscr <- list(climatic_variables = cv, edaphic_variables = ev, selected_variables = jv, kept_vars = jv)
  } else {
    cfg <- if (exists("load_step05_macroecology_config", mode = "function")) load_step05_macroecology_config() else list(cor_threshold = 0.70, vif_threshold = 5.0)
    bscr <- if (exists("screen_collinear_predictors_by_block", mode = "function")) screen_collinear_predictors_by_block(Xr, cor.threshold = cfg$cor_threshold, vif.threshold = cfg$vif_threshold, anchor_groups = resolve_active_anchor_groups(dc, intraspecific_mode = cal$is_intra)) else { cc <- grep("^(bio|delta_doy|doy)", env_vars, value = TRUE, ignore.case = TRUE); list(climatic_variables = cc, edaphic_variables = setdiff(env_vars, cc), selected_variables = env_vars, kept_vars = env_vars) }
    cv <- intersect(bscr$climatic_variables, colnames(Xr)); if (length(cv) == 0L) cv <- intersect(grep("^(bio|delta_doy|doy)", colnames(Xr), value = TRUE, ignore.case = TRUE), colnames(Xr))
    ev <- intersect(unique(c(bscr$edaphic_variables, bscr$topo_variables)), colnames(Xr)); if (length(ev) == 0L) ev <- setdiff(colnames(Xr), cv); jv <- unique(c(cv, ev))
  }
  ccal <- if (all(c("longitude", "latitude") %in% names(dc))) as.matrix(dc[, c("longitude", "latitude"), drop = FALSE]) else NULL
  Xcl <- if (length(cv) > 0L) Xr[, cv, drop = FALSE] else Xr; Xed <- if (length(ev) > 0L) Xr[, ev, drop = FALSE] else Xr; Xjt <- Xr[, jv, drop = FALSE]
  pcl <- .prepare_partition_scores_glob(Xcl, NULL, ccal); ped <- .prepare_partition_scores_glob(Xed, NULL, ccal); pjt <- .prepare_partition_scores_glob(Xjt, NULL, ccal)
  mk_row <- function(pn, s1, s2, n1, n2, dc_, ic_, de_, ie_, dj_, ij_, ps_, dax_, outc_) data.frame(Pair = pn, Taxon1 = s1, Taxon2 = s2, N1 = as.integer(n1), N2 = as.integer(n2), D_clim = dc_, I_clim = ic_, p_equiv_D_clim = ps_[1L], p_bg_D_clim = ps_[2L], D_edaph = de_, I_edaph = ie_, p_equiv_D_edaph = ps_[3L], p_bg_D_edaph = ps_[4L], D_joint = dj_, I_joint = ij_, Schoener_D = dj_, Warren_D = dj_, Warren_I = ij_, p_value_D = ps_[5L], p_value_I = ps_[6L], p_bg_D_joint = ps_[7L], p_bg_I_joint = ps_[8L], Niche_Divergence_Axis = dax_, Taxonomic_Inference = outc_, stringsAsFactors = FALSE)
  if (length(pairs) == 0L && length(cg) == 1L) {
    tx1 <- cg[1L]; n1 <- sum(gc == tx1, na.rm = TRUE)
    return(list(summary = mk_row(sprintf("%s (Single Lineage K=1)", gsub("Packera ", "P.", tx1)), tx1, tx1, n1, n1, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, rep(1.0, 8L), "Unimodal_Niche_Breadth", "Single Unimodal Species (Facultative Niche Breadth; K=1)"), permutations = data.frame(), collinearity_blocks = bscr))
  }
  vm <- list(); pj <- list()
  for (pr in pairs) {
    n1 <- sum(gc == pr[1L], na.rm = TRUE); n2 <- sum(gc == pr[2L], na.rm = TRUE); if (n1 < 5L || n2 < 5L) next
    vm[[length(vm) + 1L]] <- list(sp1 = pr[1L], sp2 = pr[2L], pn = paste0(gsub("Packera ", "P.", pr[1L]), " vs ", gsub("Packera ", "P.", pr[2L])), n1 = n1, n2 = n2)
    pj <- c(pj, list(list(part_id = "clim", bg_mat = Xcl, sp1 = pr[1L], sp2 = pr[2L], rf_var_weights = rf_var_weights, prep = pcl), list(part_id = "edaph", bg_mat = Xed, sp1 = pr[1L], sp2 = pr[2L], rf_var_weights = rf_var_weights, prep = ped), list(part_id = "joint", bg_mat = Xjt, sp1 = pr[1L], sp2 = pr[2L], rf_var_weights = rf_var_weights, prep = pjt)))
  }
  po <- .evaluate_partition_jobs_parallel(pj, gc, if (isTRUE(psw)) 0L else n_perm); rl <- list(); pl <- list()
  for (pi in seq_along(vm)) {
    md <- vm[[pi]]; cr <- po[[(pi - 1L) * 3L + 1L]]; er <- po[[(pi - 1L) * 3L + 2L]]; jr <- po[[(pi - 1L) * 3L + 3L]]
    dax <- ifelse(!is.na(er$D_emp) && !is.na(cr$D_emp) && (er$D_emp + 0.08 < cr$D_emp), "Edaphic_Lithological_Divergence", ifelse(!is.na(cr$D_emp) && !is.na(er$D_emp) && (cr$D_emp + 0.08 < er$D_emp), "Macroclimatic_Divergence", "Joint_Climato_Edaphic"))
    if (isTRUE(psw)) { outc <- "Post-Clustering Ecological Regime Description (p-values suppressed)"; ps <- rep(NA_real_, 8L) } else {
      mbg <- suppressWarnings(min(c(jr$p_bg_D, er$p_bg_D, cr$p_bg_D), na.rm = TRUE)); meq <- suppressWarnings(min(c(jr$p_equiv_D, er$p_equiv_D, cr$p_equiv_D), na.rm = TRUE)); if (!is.finite(mbg)) mbg <- meq
      mde <- suppressWarnings(min(c(jr$D_emp, er$D_emp, cr$D_emp), na.rm = TRUE))
      outc <- ifelse(is.finite(mbg) && mbg < 0.01 && is.finite(mde) && mde < 0.30, "Distinct Ecological Niche (Background Similarity Divergence Confirmed)", ifelse((is.finite(mbg) && mbg < 0.05 && is.finite(mde) && mde < 0.50) || (is.finite(jr$p_equiv_D) && jr$p_equiv_D < 0.05 && is.finite(jr$D_emp) && jr$D_emp < 0.55) || (is.finite(meq) && meq < 0.05 && is.finite(mde) && mde < 0.45), "Significant Divergence (Subspecific Ecological Race)", "Conserved / Overlapping Niche (Shared Regional Background)"))
      ps <- round(c(cr$p_equiv_D, cr$p_bg_D, er$p_equiv_D, er$p_bg_D, jr$p_equiv_D, jr$p_equiv_I, jr$p_bg_D, jr$p_bg_I), 4)
    }
    rl[[md$pn]] <- mk_row(md$pn, md$sp1, md$sp2, md$n1, md$n2, round(cr$D_emp, 4), round(cr$I_emp, 4), round(er$D_emp, 4), round(er$I_emp, 4), round(jr$D_emp, 4), round(jr$I_emp, 4), ps, dax, outc)
    if (length(jr$D_null) > 0L) pl[[md$pn]] <- data.frame(Pair = md$pn, D_null = jr$D_null, I_null = jr$I_null, D_clim_null = cr$D_null, D_edaph_null = er$D_null)
  }
  sdf <- if (length(rl) > 0L) do.call(rbind, rl) else mk_row("", "", "", 0L, 0L, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, rep(NA_real_, 8L), "", "")[0L, , drop = FALSE]
  list(summary = sdf, permutations = if (length(pl) > 0L) do.call(rbind, pl) else data.frame(), collinearity_blocks = bscr)
}

compute_specimen_edaphic_z <- function(df, primary_edaph_vars = NULL, primary_clim_vars = NULL, primary_env_vars = NULL, rf_var_weights = NULL, intraspecific_mode = FALSE) {
  n <- if (is.data.frame(df)) nrow(df) else 0L
  if (n == 0L) return(list(edaphic_z_score = numeric(0), climatic_z_score = numeric(0), env_z_score = numeric(0), primary_edaph_vars = character(0), primary_clim_vars = character(0), primary_env_vars = character(0), eval_specimen_z = function(idx, target_tx, mode = "edaphic") 0.0))
  ii <- isTRUE(intraspecific_mode) || (("species_standardized" %in% names(df)) && length(unique(stats::na.omit(df$species_standardized))) < 2L)
  if (ii && !("analysis_group" %in% names(df))) df <- resolve_analysis_group(df, TRUE, primary_env_vars = primary_env_vars, rf_var_weights = rf_var_weights)
  gv <- if (ii && "analysis_group" %in% names(df)) as.character(df$analysis_group) else as.character(df$species_standardized)
  atx <- sort(unique(gv[!is.na(gv) & nzchar(gv) & gv != "Unknown"])); cp <- "^(bio[0-9]*_|summer_moisture_|maritime_continentality_)"; ca <- discover_candidate_eco_predictors(df, FALSE, 3L)
  if (length(primary_env_vars) == 0L) { primary_env_vars <- unique(c(as.character(primary_edaph_vars), as.character(primary_clim_vars))); if (length(primary_env_vars) == 0L) primary_env_vars <- ca }
  pev <- intersect(unique(as.character(primary_env_vars)), names(df)); pev <- pev[vapply(pev, function(cl) sum(is.finite(suppressWarnings(as.numeric(df[[cl]])))) >= min(3L, max(2L, n)), logical(1))]
  pcv <- if (length(primary_clim_vars) == 0L) { c1 <- grep(cp, pev, value = TRUE, ignore.case = TRUE); if (length(c1) > 0L) c1 else grep(cp, ca, value = TRUE, ignore.case = TRUE) } else intersect(unique(as.character(primary_clim_vars)), names(df))
  ped <- if (length(primary_edaph_vars) == 0L) { e1 <- setdiff(pev, grep(cp, pev, value = TRUE, ignore.case = TRUE)); if (length(e1) > 0L) e1 else setdiff(ca, grep(cp, ca, value = TRUE, ignore.case = TRUE)) } else intersect(unique(as.character(primary_edaph_vars)), names(df))
  azv <- unique(c(pev, ped, pcv)); azv <- azv[vapply(azv, function(cl) (cl %in% names(df)) && sum(is.finite(suppressWarnings(as.numeric(df[[cl]])))) >= min(3L, max(2L, n)), logical(1))]
  af <- .s05_resolve_anchor_flags(df, gv, atx, ii); aw <- resolve_ecological_authority_weights(df); aw[!is.finite(aw) | aw < 0.0] <- 0.0
  sw <- .s05_spatial_wts(df, af$is_anchor); rw <- ifelse(af$is_anchor, pmax(aw, 0.25), 0.0) * sw; sf <- .s05_scale_floors(df, azv, af$is_anchor); zp <- list()
  if (length(atx) > 0L && length(azv) > 0L) {
    for (tx in atx) {
      ix <- which(af$is_anchor & gv == tx); if (length(ix) == 0L) ix <- which(af$tier_ok & gv == tx); if (length(ix) == 0L) ix <- which(!af$cl_noise & af$elig_coord & gv == tx); if (length(ix) == 0L) ix <- which(gv == tx)
      wt <- rw[ix]; if (!any(wt > 0.0)) wt <- sw[ix]
      zp[[tx]] <- list(mu = stats::setNames(vapply(azv, function(cl) .s05_w_mean(suppressWarnings(as.numeric(df[[cl]][ix])), wt), numeric(1)), azv), sd = stats::setNames(vapply(azv, function(cl) .s05_w_sd(suppressWarnings(as.numeric(df[[cl]][ix])), wt, sf[[cl]]), numeric(1)), azv))
    }
    for (cl in azv) { bsd <- mean(vapply(zp, function(pr) pr$sd[[cl]], numeric(1)), na.rm = TRUE); if (is.finite(bsd) && bsd > 0.0) for (tx in atx) zp[[tx]]$sd[[cl]] <- max(0.5 * zp[[tx]]$sd[[cl]] + 0.5 * bsd, sf[[cl]]) }
  }
  rsw <- function(vs) {
    if (length(vs) == 0L) return(numeric(0)); wo <- stats::setNames(rep(1.0 / length(vs), length(vs)), vs)
    if (!is.null(rf_var_weights) && length(rf_var_weights) > 0L && !is.null(names(rf_var_weights)) && length(intersect(vs, names(rf_var_weights))) > 0L) {
      ht <- intersect(vs, names(rf_var_weights)); mp <- stats::median(as.numeric(rf_var_weights[rf_var_weights > 0.0]), na.rm = TRUE); if (!is.finite(mp) || mp <= 0.0) mp <- 1.0 / length(vs)
      wr <- stats::setNames(rep(0.25 * mp, length(vs)), vs); wh <- suppressWarnings(as.numeric(rf_var_weights[ht])); wh[!is.finite(wh) | wh <= 0.0] <- 0.25 * mp; wr[ht] <- wh; wo <- wr / sum(wr)
    }; wo
  }
  edu <- if (length(ped) > 0L) intersect(ped, azv) else intersect(pev, azv); clu <- intersect(pcv, azv); enu <- if (length(pev) > 0L) intersect(pev, azv) else azv
  wed <- rsw(edu); wcl <- rsw(clu); wen <- rsw(enu)
  ezf <- function(idx, ttx, mode = "edaphic") {
    if (is.na(ttx) || !(as.character(ttx) %in% names(zp))) return(if (identical(mode, "climatic")) NA_real_ else 0.0)
    vm <- if (identical(mode, "climatic")) clu else if (identical(mode, "env")) enu else edu; wm <- if (identical(mode, "climatic")) wcl else if (identical(mode, "env")) wen else wed
    if (length(vm) == 0L) return(if (identical(mode, "climatic")) NA_real_ else 0.0)
    pr <- zp[[as.character(ttx)]]; z2 <- ((vapply(vm, function(cl) suppressWarnings(as.numeric(df[[cl]][idx])), numeric(1)) - pr$mu[vm]) / pmax(pr$sd[vm], sf[vm], 1e-4))^2
    ok <- is.finite(z2) & is.finite(wm) & (wm > 0.0); if (!any(ok)) return(if (identical(mode, "climatic")) NA_real_ else 0.0); sqrt(sum(wm[ok] * z2[ok]) / sum(wm[ok]))
  }
  list(edaphic_z_score = round(vapply(seq_len(n), function(i) ezf(i, gv[i], "edaphic"), numeric(1)), 4), climatic_z_score = round(vapply(seq_len(n), function(i) ezf(i, gv[i], "climatic"), numeric(1)), 4), env_z_score = round(vapply(seq_len(n), function(i) ezf(i, gv[i], "env"), numeric(1)), 4), primary_edaph_vars = edu, primary_clim_vars = clu, primary_env_vars = enu, rf_var_weights = wen, eval_specimen_z = ezf)
}

.reconstruct_rf_prob_mat <- function(n, active_taxa, grp_vec, p_tx, p_pr, g_pr) {
  K <- length(active_taxa); if (is.null(p_tx) || is.null(p_pr) || length(p_tx) != n || length(p_pr) != n || K < 2L) return(NULL)
  mat <- matrix(NA_real_, nrow = n, ncol = K, dimnames = list(NULL, active_taxa))
  for (ri in seq_len(n)) {
    ptx <- as.character(p_tx[ri]); ppr <- suppressWarnings(as.numeric(p_pr[ri])); if (is.na(ptx) || !(ptx %in% active_taxa) || !is.finite(ppr)) next
    ppr <- min(0.999, max(1e-4, ppr)); gtx <- as.character(grp_vec[ri]); gpr <- if (!is.null(g_pr) && length(g_pr) == n) suppressWarnings(as.numeric(g_pr[ri])) else NA_real_
    if (!is.na(gtx) && (gtx %in% active_taxa) && gtx != ptx && is.finite(gpr)) {
      gpr <- min(1.0 - ppr, max(0.0, gpr)); rem <- setdiff(active_taxa, c(ptx, gtx))
      mat[ri, ] <- if (length(rem) > 0L) max(0.0, 1.0 - ppr - gpr) / length(rem) else 0.0; mat[ri, ptx] <- ppr; mat[ri, gtx] <- gpr
    } else { mat[ri, ] <- (1.0 - ppr) / max(1L, K - 1L); mat[ri, ptx] <- ppr }
  }
  mat
}
.compute_mahal_t_posteriors <- function(df, adf, ec, atx, gv, agrp, msd) {
  n <- nrow(df); K <- length(atx); out <- matrix(1.0 / max(1L, K), n, max(1L, K), dimnames = list(NULL, if (K > 0L) atx else "Unknown"))
  if (n == 0L || K < 2L || length(ec) < 1L) return(out)
  Ear <- do.call(cbind, lapply(ec, function(cl) suppressWarnings(as.numeric(df[[cl]])))); Ecr <- do.call(cbind, lapply(ec, function(cl) suppressWarnings(as.numeric(adf[[cl]]))))
  colnames(Ear) <- ec; colnames(Ecr) <- ec; mc <- colnames(.s05_screen_point_mass_mat(Ecr, 0.25, 2L)); pe <- length(mc)
  Ea <- Ear[, mc, drop = FALSE]; Ec <- Ecr[, mc, drop = FALSE]; wa <- ifelse(adf$case_wts > 0.0, adf$case_wts, adf$spatial_balance_weight); ms_v <- numeric(pe); sd_v <- numeric(pe)
  for (j in seq_len(pe)) {
    mj <- stats::median(Ec[is.finite(Ec[, j]), j], na.rm = TRUE); if (!is.finite(mj)) mj <- stats::median(Ea[is.finite(Ea[, j]), j], na.rm = TRUE); if (!is.finite(mj)) mj <- 0.0
    Ea[!is.finite(Ea[, j]), j] <- mj; Ec[!is.finite(Ec[, j]), j] <- mj; fl <- if (mc[j] %in% names(msd)) msd[[mc[j]]] else 0.25
    muj <- .s05_w_mean(Ec[, j], wa); ms_v[j] <- if (is.finite(muj)) muj else mean(Ea[, j]); sd_v[j] <- max(.s05_w_sd(Ec[, j], wa, fl), fl, 1e-4)
  }
  Za <- sweep(sweep(Ea, 2L, ms_v, "-"), 2L, sd_v, "/"); Zc <- sweep(sweep(Ec, 2L, ms_v, "-"), 2L, sd_v, "/"); Za[!is.finite(Za)] <- 0.0; Zc[!is.finite(Zc)] <- 0.0
  pca <- if (pe >= 2L && nrow(Zc) >= 3L) tryCatch(stats::prcomp(Zc, center = TRUE, scale. = FALSE), error = function(e) NULL) else NULL
  if (!is.null(pca) && ncol(pca$rotation) >= 1L) {
    ev <- pca$sdev^2; k9 <- which((if (sum(ev) > 1e-12) cumsum(ev) / sum(ev) else rep(1.0, length(ev))) >= 0.90)[1L]
    kp <- max(1L, min(if (is.na(k9)) ncol(pca$rotation) else k9, min(6L, ncol(pca$rotation), max(1L, nrow(Zc) - K))))
    rk <- pca$rotation[, seq_len(kp), drop = FALSE]; Zall <- sweep(Za, 2L, pca$center, "-") %*% rk; Zanc <- sweep(Zc, 2L, pca$center, "-") %*% rk
  } else { kp <- ncol(Za); Zall <- Za; Zanc <- Zc }
  mup <- matrix(0.0, K, kp, dimnames = list(atx, NULL)); rch <- list(); wtp <- 0.0; nop <- 0L
  for (tx in atx) {
    ia <- which(!is.na(agrp) & agrp == tx)
    if (length(ia) >= 2L) { Zk <- Zanc[ia, , drop = FALSE]; wk <- pmax(1e-3, wa[ia]) } else { ifu <- which(!is.na(gv) & gv == tx); Zk <- Zall[ifu, , drop = FALSE]; wk <- pmax(1e-3, df$spatial_balance_weight[ifu]) }
    if (nrow(Zk) >= 1L) { sw <- sum(wk); mk <- colSums(Zk * wk) / max(sw, 1e-8); mup[tx, ] <- mk; rch[[length(rch) + 1L]] <- sweep(Zk, 2L, mk, "-") * sqrt(wk); wtp <- wtp + sw; nop <- nop + nrow(Zk) }
  }
  Sw <- crossprod(if (length(rch) > 0L) do.call(rbind, rch) else Zanc) / max(wtp - K, 1.0); Sw[!is.finite(Sw)] <- 0.0; trw <- sum(diag(Sw)); tv <- if (is.finite(trw) && trw > 1e-8) trw / kp else 1.0
  lw <- if (kp == 1L) 0.10 else { tr2 <- sum(Sw * Sw); dn <- (max(nop, 3L) + 1.0 - 2.0 / kp) * (tr2 - (trw * trw) / kp); min(0.50, max(0.05, if (is.finite(dn) && dn > 1e-12) ((1.0 - 2.0 / kp) * tr2 + trw * trw) / dn else 0.15)) }
  Sr <- (1.0 - lw) * Sw + lw * diag(tv, kp, kp); Si <- tryCatch(solve(Sr), error = function(e) diag(1.0 / pmax(diag(Sr), 1e-4), kp, kp))
  d2m <- matrix(0.0, n, K, dimnames = list(NULL, atx)); for (tx in atx) { dfx <- sweep(Zall, 2L, mup[tx, ], "-"); d2m[, tx] <- pmax(0.0, rowSums((dfx %*% Si) * dfx)) }
  nu <- max(3.0, kp + 2.0); llm <- -0.5 * (nu + kp) * log1p(d2m / nu); pm <- exp(sweep(llm, 1L, apply(llm, 1L, max), "-")); sweep(pm, 1L, pmax(rowSums(pm), 1e-12), "/")
}

execute_crossmodal_consensus <- function(df, intraspecific_mode = FALSE, rf_res = NULL, primary_env_vars = NULL, primary_clim_vars = NULL, primary_edaph_vars = NULL, rf_var_weights = NULL) {
  message("Executing 3-way cross-modal consensus checks (Morpho + Fine-Grain Edaphic + Phenology)...")
  if (exists("flag_coordinate_uncertainty", mode = "function")) df <- flag_coordinate_uncertainty(df, 10000.0)
  if (!("analysis_group" %in% names(df)) || all(is.na(df$analysis_group) | !nzchar(trimws(as.character(df$analysis_group))))) df <- resolve_analysis_group(df, intraspecific_mode, primary_env_vars = primary_env_vars, rf_var_weights = NULL)
  ii <- isTRUE(intraspecific_mode) || (length(unique(stats::na.omit(df$species_standardized))) < 2L)
  atx <- if (ii) sort(unique(stats::na.omit(df$analysis_group))) else get_active_packera_taxa(df$species_standardized)
  if (length(atx) == 0L) atx <- sort(unique(as.character(stats::na.omit(df$analysis_group))))
  gv <- if (ii) df$analysis_group else df$species_standardized; n <- nrow(df); K <- length(atx)
  for (pc in list(c("fine_soil_ph", "soil_ph"), c("fine_soil_sand", "soil_sand"), c("fine_soil_clay", "soil_clay"), c("fine_bedrock_depth_cm", "ssurgo_bedrock_depth_cm"), c("fine_soil_cec", "soil_cec"))) if (!(pc[1L] %in% names(df))) df[[pc[1L]]] <- if (pc[2L] %in% names(df)) suppressWarnings(as.numeric(df[[pc[2L]]])) else rep(NA_real_, n)
  if (!("ssurgo_rock_outcrop_pct" %in% names(df))) df$ssurgo_rock_outcrop_pct <- rep(NA_real_, n)
  rdoy <- if ("doy" %in% names(df)) suppressWarnings(as.numeric(df$doy)) else rep(NA_real_, n); vdoy <- is.finite(rdoy) & rdoy >= 1.0 & rdoy <= 366.0
  vanth <- if ("is_verified_anthesis" %in% names(df)) .s05_bcol(df, "is_verified_anthesis", FALSE) & vdoy else { rc <- if ("reproductiveCondition" %in% names(df)) tolower(as.character(df$reproductiveCondition)) else rep("", n); rc[is.na(rc)] <- ""; vdoy & !(grepl("sterile|fruit", rc) & !grepl("fl", rc)) }
  hdoy <- if ("hopkins_adjusted_doy" %in% names(df)) suppressWarnings(as.numeric(df$hopkins_adjusted_doy)) else rep(NA_real_, n)
  dth <- ifelse(vanth & is.finite(hdoy), hdoy, ifelse(vanth, rdoy, NA_real_)); df$pheno_theta <- ifelse(is.finite(dth), 2.0 * pi * dth / 365.25, NA_real_)
  af <- .s05_resolve_anchor_flags(df, gv, atx, ii)
  df$authority_weight <- resolve_ecological_authority_weights(df); df$spatial_balance_weight <- .s05_spatial_wts(df, af$is_anchor)
  df$case_wts <- ifelse(af$is_anchor, pmax(df$authority_weight, 0.25), 0.0) * df$spatial_balance_weight
  if (length(primary_env_vars) == 0L && is.list(rf_res) && !is.null(rf_res$importance)) {
    sa <- select_primary_predictors_from_rf(rf_res, rf_res$collinearity_screening, df, intraspecific_mode = ii)
    if (length(sa$primary_env_vars) > 0L) { primary_env_vars <- sa$primary_env_vars; if (length(primary_clim_vars) == 0L) primary_clim_vars <- sa$primary_clim_vars; if (length(primary_edaph_vars) == 0L) primary_edaph_vars <- sa$primary_edaph_vars; if (length(rf_var_weights) == 0L) rf_var_weights <- sa$rf_var_weights }
  }
  adf <- df[af$is_anchor, , drop = FALSE]; if (nrow(adf) == 0L) adf <- df[gv %in% atx, , drop = FALSE]
  ce <- discover_candidate_eco_predictors(if (nrow(adf) >= 3L) adf else df, FALSE, 3L)
  ec <- if (length(intersect(primary_env_vars, ce)) > 0L) intersect(as.character(primary_env_vars), ce) else ce
  if (length(ec) == 0L) { fb <- intersect(c("fine_soil_ph", "fine_soil_sand", "fine_soil_cec", "soil_ph", "soil_sand", "soil_cec", unname(WC_BIO_LOOKUP)), names(df)); ec <- fb[vapply(fb, function(cl) any(is.finite(suppressWarnings(as.numeric(df[[cl]])))), logical(1))] }
  msd <- .s05_scale_floors(df, ec, af$is_anchor); agrp <- if (ii) as.character(adf$analysis_group) else as.character(adf$species_standardized)
  tprof <- stats::setNames(lapply(atx, function(tx) {
    ix <- which(af$is_anchor & gv == tx); if (length(ix) == 0L) ix <- which(gv == tx & (!af$cl_noise | af$is_type_clean)); if (length(ix) == 0L) ix <- which(gv == tx)
    ws <- df$case_wts[ix]; if (!any(ws > 0.0)) ws <- df$spatial_balance_weight[ix]
    sm <- .s05_w_mean(sin(df$pheno_theta[ix]), ws); cm <- .s05_w_mean(cos(df$pheno_theta[ix]), ws); list(mean_theta = atan2(sm, cm), R_bar = sqrt(sm^2 + cm^2))
  }), atx)
  zeng <- compute_specimen_edaphic_z(df, primary_edaph_vars, primary_clim_vars, ec, rf_var_weights, ii)
  roa <- function(d) {
    ns <- nrow(d); if (ns == 0L) return(logical(0))
    rp <- if ("ssurgo_rock_outcrop_pct" %in% names(d)) suppressWarnings(as.numeric(d$ssurgo_rock_outcrop_pct)) else rep(0.0, ns); rp[is.na(rp)] <- 0.0
    bd <- if ("fine_bedrock_depth_cm" %in% names(d)) suppressWarnings(as.numeric(d$fine_bedrock_depth_cm)) else rep(150.0, ns); bd[is.na(bd)] <- 150.0
    mb <- if ("min_bedrock_depth_buffer_cm" %in% names(d)) suppressWarnings(as.numeric(d$min_bedrock_depth_buffer_cm)) else bd; mb[is.na(mb)] <- bd[is.na(mb)]
    res_o <- (rp >= 15.0) | (bd <= 65.0) | (mb <= 50.0) | .s05_bcol(d, "is_rock_outcrop", FALSE) | .s05_bcol(d, "is_azonal_lithology", FALSE)
    if ("azonal_lithology_presence_buffer" %in% names(d)) res_o <- res_o | (!is.na(d$azonal_lithology_presence_buffer) & suppressWarnings(as.integer(d$azonal_lithology_presence_buffer)) > 0L); res_o
  }
  aoa <- roa(adf); foa <- roa(df)
  is_spec <- stats::setNames(vapply(atx, function(tx) { tm <- !is.na(agrp) & agrp == tx; if (sum(tm) >= 3L) { a <- sum(tm & aoa); b <- sum(tm & !aoa); c0 <- sum(!tm & aoa); d0 <- sum(!tm & !aoa) } else { fm <- !is.na(gv) & gv == tx; a <- sum(fm & foa); b <- sum(fm & !foa); c0 <- sum(!fm & foa); d0 <- sum(!fm & !foa) }; isTRUE((a / max(1L, a + b)) > 0.30 || (if ((c0 + d0) > 0L) ((a + 0.5) * (d0 + 0.5)) / ((b + 0.5) * (c0 + 0.5)) else ((a + 0.5) * 92.5) / ((b + 0.5) * 8.5)) > 2.0) }, logical(1)), atx)
  mpm <- .compute_mahal_t_posteriors(df, adf, ec, atx, gv, agrp, msd); rfm <- NULL; rfa <- NA_real_
  if (is.list(rf_res)) {
    for (an in c("Spatial_Block_CV_Accuracy", "Spatial_Block_CV_Acc", "spatial_block_cv_acc", "block_cv_acc", "OOB_Accuracy")) if (!is.null(rf_res[[an]]) && is.finite(suppressWarnings(as.numeric(rf_res[[an]][1L])))) { rfa <- as.numeric(rf_res[[an]][1L]); break }
    cm <- if (!is.null(rf_res$oof_prob_matrix)) rf_res$oof_prob_matrix else rf_res$full_prob_mat
    if (is.matrix(cm) && nrow(cm) == n && K >= 2L && length(intersect(colnames(cm), atx)) >= 1L) { rfm <- matrix(NA_real_, n, K, dimnames = list(NULL, atx)); for (ctx in intersect(colnames(cm), atx)) rfm[, ctx] <- suppressWarnings(as.numeric(cm[, ctx])) } else rfm <- .reconstruct_rf_prob_mat(n, atx, gv, if (!is.null(rf_res$rf_predicted_taxon)) rf_res$rf_predicted_taxon else rf_res$predicted_taxon, if (!is.null(rf_res$rf_predicted_prob)) rf_res$rf_predicted_prob else rf_res$predicted_prob, rf_res$rf_given_taxon_prob)
  } else if (K >= 2L && all(c("rf_predicted_taxon", "rf_predicted_prob") %in% names(df))) {
    if ("spatial_block_cv_acc" %in% names(df)) rfa <- stats::median(suppressWarnings(as.numeric(df$spatial_block_cv_acc)), na.rm = TRUE)
    rfm <- .reconstruct_rf_prob_mat(n, atx, gv, df$rf_predicted_taxon, df$rf_predicted_prob, if ("rf_given_taxon_prob" %in% names(df)) df$rf_given_taxon_prob else NULL)
  }
  chp <- 1.0 / max(2.0, K); if (!is.finite(rfa)) rfa <- 0.75
  wrf <- if (rfa > chp) min(0.85, max(0.35, (rfa - chp) / max(1e-6, 1.0 - chp))) else max(0.10, 0.35 * max(0.0, rfa / max(chp, 1e-6)))
  cf <- logical(n); tc <- character(n); cv_o <- numeric(n); nt <- character(n); eb <- character(n); ezv <- rep(NA_real_, n); czv <- rep(NA_real_, n); sf <- logical(n); nob <- integer(n)
  hm_v <- if ("cda_predicted_taxon" %in% names(df)) !is.na(df$cda_predicted_taxon) & nzchar(trimws(as.character(df$cda_predicted_taxon))) else rep(FALSE, n)
  he_v <- if (length(ec) > 0L) vapply(seq_len(n), function(ri) any(vapply(ec, function(cl) is.finite(suppressWarnings(as.numeric(df[[cl]][ri]))), logical(1))), logical(1)) else rep(FALSE, n)
  chyb <- !ii && (K >= 2L) && (sum(grepl("^Packera[[:space:]]+[a-z][a-z-]+$", atx)) >= 2L)
  ls <- if ("decimalLatitude" %in% names(df)) suppressWarnings(as.numeric(df$decimalLatitude)) else rep(NA_real_, n)
  lo <- if ("decimalLongitude" %in% names(df)) suppressWarnings(as.numeric(df$decimalLongitude)) else rep(NA_real_, n)
  cs <- if ("catalogNumber" %in% names(df)) as.character(df$catalogNumber) else as.character(seq_len(n))
  cbt <- if (chyb) stats::setNames(lapply(atx, function(tx) { ia <- which(gv == tx & af$is_anchor & is.finite(ls) & is.finite(lo)); if (length(ia) > 0L) ia else which(gv == tx & is.finite(ls) & is.finite(lo)) }), atx) else list()
  mpk <- function(ri, tx) { if (!is.finite(ls[ri]) || !is.finite(lo[ri])) return(Inf); ia <- cbt[[tx]]; ia <- ia[cs[ia] != cs[ri]]; if (length(ia) == 0L) return(Inf); dl <- (ls[ia] - ls[ri]) * pi / 180.0; do <- (lo[ia] - lo[ri]) * pi / 180.0; ah <- sin(dl * 0.5)^2 + cos(ls[ri] * pi / 180.0) * cos(ls[ia] * pi / 180.0) * sin(do * 0.5)^2; min(6371.0 * 2.0 * atan2(sqrt(pmax(0.0, ah)), sqrt(pmax(0.0, 1.0 - ah)))) }
  cdm <- if ("cda_posterior_margin" %in% names(df)) suppressWarnings(as.numeric(df$cda_posterior_margin)) else rep(NA_real_, n)
  for (i in seq_len(n)) {
    gsp <- gv[i]; hm <- isTRUE(hm_v[i]); he <- isTRUE(he_v[i]); hp <- isTRUE(vanth[i]) && is.finite(df$pheno_theta[i]) && !(("is_pheno_uniform_null" %in% names(df)) && isTRUE(as.logical(df$is_pheno_uniform_null[i])))
    msp <- if (hm) as.character(df$cda_predicted_taxon[i]) else gsp; if (ii && !(msp %in% atx)) msp <- gsp
    cp <- if ("cda_posterior_prob" %in% names(df) && !is.na(df$cda_posterior_prob[i])) as.numeric(df$cda_posterior_prob[i]) else NA_real_
    mo <- "morphometric_outlier_flag" %in% names(df) && isTRUE(as.logical(df$morphometric_outlier_flag[i]))
    if (!he) { esp <- NA_character_; ez <- NA_real_; cz <- NA_real_ } else {
      if (K <= 1L) esp <- if (K == 1L) atx[1L] else gsp else {
        pmi <- pmax(as.numeric(mpm[i, ]), 1e-6); pmi <- pmi / sum(pmi)
        if (!is.null(rfm) && any(is.finite(rfm[i, ]))) { pri <- as.numeric(rfm[i, ]); nar <- !is.finite(pri); if (any(nar)) pri[nar] <- max(0.0, 1.0 - sum(pri[!nar])) / sum(nar); pri <- pmax(pri, 1e-6); pri <- pri / sum(pri); sce <- wrf * log(pri) + (1.0 - wrf) * log(pmi) } else sce <- log(pmi)
        esp <- if (all(!is.finite(sce))) gsp else atx[which.max(sce)][1L]
      }
      ez <- zeng$eval_specimen_z(i, gsp, "edaphic"); cz <- zeng$eval_specimen_z(i, gsp, "climatic")
    }
    eb[i] <- esp; ezv[i] <- if (is.finite(ez)) round(ez, 3) else NA_real_; czv[i] <- if (is.finite(cz)) round(cz, 3) else NA_real_
    ap <- if (hm && msp != gsp) msp else if (he && !is.na(esp)) esp else gsp; eza <- if (he) zeng$eval_specimen_z(i, ap, "edaphic") else NA_real_
    sf[i] <- isTRUE(foa[i] && !is.na(gsp) && (as.character(gsp) %in% names(is_spec)) && isTRUE(is_spec[[as.character(gsp)]]))
    pg <- tprof[[gsp]]; ad <- if (hp && !is.null(pg) && is.finite(pg$mean_theta)) abs(atan2(sin(df$pheno_theta[i] - pg$mean_theta), cos(df$pheno_theta[i] - pg$mean_theta))) else NA_real_
    rgp <- if (!is.null(rfm) && !is.na(gsp) && (gsp %in% colnames(rfm))) suppressWarnings(as.numeric(rfm[i, gsp])) else if ("rf_given_taxon_prob" %in% names(df)) suppressWarnings(as.numeric(df$rf_given_taxon_prob[i])) else NA_real_
    pag <- if (is.finite(ad)) isTRUE(ad < (pi / 3.0)) else NA; if (is.finite(cz) && cz > 3.5 && (!is.finite(rgp) || rgp < 0.25)) pag <- NA
    om <- c(if (hm) (msp == gsp) else NA, if (he && !is.na(esp)) (esp == gsp) else NA, pag); noi <- sum(!is.na(om)); nob[i] <- as.integer(noi); cv_o[i] <- if (noi > 0L) round(sum(om, na.rm = TRUE) / noi, 2) else NA_real_
    edd <- he && !is.na(esp) && (esp != gsp); mod <- hm && (msp != gsp); phd <- isTRUE(!is.na(pag) && !pag)
    if (mod && edd && (esp == msp || phd || mo || (!is.na(cp) && cp >= 0.70))) { cf[i] <- TRUE; tc[i] <- "Severe_Triple_Stream_Conflict"; nt[i] <- sprintf("Given %s contradicted by Morphology (%s) and Fine-Grain Edaphic (%s)", gsp, msp, esp)
    } else if (he && is.finite(ez) && ez > 3.0 && !mod) { cf[i] <- TRUE; tc[i] <- "Edaphic_Envelope_Violation"; nt[i] <- sprintf("Specimen outside physiological primary-edaphic tolerance (%s; edaphic_z=%.2f)", paste(utils::head(zeng$primary_edaph_vars, 4L), collapse = ","), ez)
    } else if (phd && !mod && !edd) { cf[i] <- TRUE; tc[i] <- "Phenological_Anomaly_Outlier"; nt[i] <- sprintf("Verified-anthesis DOY %.0f diverges from seasonal peak (ang_diff=%.2f rad)", rdoy[i], ad)
    } else if (mod || edd) {
      mi <- hm && !is.na(cp) && (cp >= 0.32 && cp <= 0.65) && (is.na(cdm[i]) || cdm[i] <= 0.22) && !(msp == gsp && cp >= 0.62)
      syn <- chyb && (ap != gsp) && (mpk(i, gsp) <= 25.0) && (mpk(i, ap) <= 25.0)
      if (chyb && he && mi && syn && is.finite(ez) && is.finite(eza) && ez <= 2.2 && eza <= 2.2) { cf[i] <- TRUE; tc[i] <- "Putative_Hybrid_Zone_Intergrade"; nt[i] <- sprintf("Syntopic contact-zone intermediate morphometrics & edaphic overlap between %s and %s", gsp, ap)
      } else if (mod && ((he && esp == msp) || mo || (!is.na(cp) && cp >= 0.70))) { cf[i] <- TRUE; tc[i] <- "Severe_Triple_Stream_Conflict"; nt[i] <- sprintf("Cross-modal misidentification signal favoring %s over given %s", ap, gsp)
      } else { cf[i] <- FALSE; tc[i] <- "Clean_MultiModal_Consensus"; nt[i] <- "Concordant across observed morphological, phenological, and ecological streams" }
    } else { cf[i] <- FALSE; tc[i] <- "Clean_MultiModal_Consensus"; nt[i] <- if (noi == 0L) "No active traits observed on ungeoreferenced/sterile sheet" else "Concordant across observed morphological, phenological, and ecological streams" }
  }
  df$n_observed_modalities <- nob; df$multimodal_concordance <- cv_o; df$multimodal_conflict_flag <- cf; df$triage_category <- tc; df$conflict_discordance_notes <- nt
  df$edaphic_best_fit_taxon <- eb; df$edaphic_z_score <- ezv; df$climatic_z_score <- czv; df$substrate_specialist_flag <- sf
  df$primary_env_predictors <- rep(paste(zeng$primary_env_vars, collapse = ";"), n); df$primary_clim_predictors <- rep(paste(zeng$primary_clim_vars, collapse = ";"), n); df$primary_edaph_predictors <- rep(paste(zeng$primary_edaph_vars, collapse = ";"), n)
  for (nc in c("hopkins_adjusted_doy", "anthesis_dominance_weight", "allochronic_bimodality_delta_bic")) if (!(nc %in% names(df))) df[[nc]] <- rep(NA_real_, n)
  if (!("hopkins_adjusted_doy_resid" %in% names(df))) df$hopkins_adjusted_doy_resid <- if ("delta_doy" %in% names(df)) suppressWarnings(as.numeric(df$delta_doy)) else rep(NA_real_, n)
  tec <- compute_type_locality_ecological_control(df, primary_env_vars = zeng$primary_env_vars, primary_edaph_vars = zeng$primary_edaph_vars, rf_var_weights = zeng$rf_var_weights, intraspecific_mode = ii)
  df$type_edaphic_z_score <- tec$type_edaphic_z_score; df$nearest_type_niche_taxon <- tec$nearest_type_niche_taxon; df$matches_type_lithology <- tec$matches_type_lithology
  df$sensu_stricto_ecotype_regime <- rep(if (!is.null(tec$sensu_stricto_ecotype_regime) && length(tec$sensu_stricto_ecotype_regime) >= 1L) as.character(tec$sensu_stricto_ecotype_regime[1L]) else NA_character_, n)
  df$sensu_stricto_edaphic_z <- if (!is.null(tec$sensu_stricto_edaphic_z) && length(tec$sensu_stricto_edaphic_z) == n) tec$sensu_stricto_edaphic_z else rep(NA_real_, n)
  df$nearest_type_niche_basionym <- if (!is.null(tec$nearest_type_niche_basionym) && length(tec$nearest_type_niche_basionym) == n) tec$nearest_type_niche_basionym else if ("nearest_type_niche_basionym" %in% names(df)) df$nearest_type_niche_basionym else rep(NA_character_, n)
  df$nearest_type_basionym_edaphic_z <- if (!is.null(tec$nearest_type_basionym_edaphic_z) && length(tec$nearest_type_basionym_edaphic_z) == n) tec$nearest_type_basionym_edaphic_z else if ("nearest_type_basionym_edaphic_z" %in% names(df)) df$nearest_type_basionym_edaphic_z else rep(NA_real_, n)
  het_div <- (!is.null(tec$is_heterotypic_type_audit) && length(tec$is_heterotypic_type_audit) == n) && any(!is.na(tec$is_heterotypic_type_audit) & as.logical(tec$is_heterotypic_type_audit))
  if (isTRUE(het_div)) { hd <- !is.na(tec$is_heterotypic_type_audit) & as.logical(tec$is_heterotypic_type_audit); df$is_heterotypic_type_audit <- hd; df$multimodal_conflict_flag[hd] <- TRUE; df$triage_category[hd & df$triage_category == "Clean_MultiModal_Consensus"] <- "Heterotypic_Type_Niche_Divergence" }
  df
}
