if (grDevices::dev.cur() == 1L) grDevices::pdf(NULL)
# ==============================================================================
# Module: scripts/core/db_r_step05_env_cache.R
# Step 05 environmental_cache R*Tree queries, upserts, and spatial_rf_predictions persistence.
# ==============================================================================

query_environmental_cache_sqlite <- function(
  catalog_numbers = NULL,
  lat_lon_keys = NULL,
  bbox = NULL,
  cache_path = NULL,
  db_path = NULL,
  demo = FALSE
) {
  has_explicit_db <- !is.null(db_path) && nzchar(trimws(as.character(db_path)))
  if (isTRUE(demo) || (!has_explicit_db && !is_sqlite_active_for_inputs(cache_path, demo = demo))) {
    return(NULL)
  }
  conn <- connect_packera_db(db_path = if (has_explicit_db) as.character(db_path) else "", read_only = TRUE, demo = FALSE)
  if (is.null(conn)) {
    return(NULL)
  }
  on.exit(disconnect_packera_db(conn), add = TRUE)
  if (!db_table_exists_safe(conn, "environmental_cache")) {
    return(NULL)
  }

  if (!is.null(bbox) && length(bbox) == 4L && all(is.finite(as.numeric(bbox))) &&
    db_table_exists_safe(conn, "environmental_rtree")) {
    bb <- as.numeric(bbox)
    res <- tryCatch(
      DBI::dbGetQuery(
        conn,
        paste0(
          "SELECT ec.* FROM environmental_cache ec ",
          "INNER JOIN environmental_rtree rt ON ec.rowid = rt.id ",
          "WHERE rt.min_lon >= $1 AND rt.max_lon <= $2 AND rt.min_lat >= $3 AND rt.max_lat <= $4;"
        ),
        params = list(bb[1], bb[2], bb[3], bb[4])
      ),
      error = function(e) NULL
    )
    if (is.null(res) || nrow(res) == 0L) {
      return(NULL)
    }
    return(coerce_sqlite_logical_cols(res))
  }

  clauses <- character(0)
  if (!is.null(catalog_numbers) && length(catalog_numbers) > 0L) {
    cats <- unique(trimws(as.character(catalog_numbers)))
    cats <- cats[!is.na(cats) & nzchar(cats)]
    if (length(cats) > 0L) {
      qc <- paste(vapply(cats, function(x) as.character(DBI::dbQuoteString(conn, x)), character(1)), collapse = ", ")
      clauses <- c(clauses, sprintf("catalogNumber IN (%s)", qc))
    }
  }
  if (!is.null(lat_lon_keys) && length(lat_lon_keys) > 0L) {
    keys <- unique(trimws(as.character(lat_lon_keys)))
    keys <- keys[!is.na(keys) & nzchar(keys)]
    if (length(keys) > 0L) {
      qk <- paste(vapply(keys, function(x) as.character(DBI::dbQuoteString(conn, x)), character(1)), collapse = ", ")
      clauses <- c(clauses, sprintf("lat_lon_key IN (%s)", qk))
    }
  }

  sql <- if (length(clauses) > 0L) {
    sprintf("SELECT * FROM environmental_cache WHERE %s;", paste(clauses, collapse = " OR "))
  } else {
    "SELECT * FROM environmental_cache;"
  }
  res <- tryCatch(DBI::dbGetQuery(conn, sql), error = function(e) NULL)
  if (is.null(res) || nrow(res) == 0L) {
    return(NULL)
  }
  coerce_sqlite_logical_cols(res)
}

upsert_environmental_cache_sqlite <- function(
  df,
  cache_path = NULL,
  demo = FALSE
) {
  if (isTRUE(demo) || is.null(df) || nrow(df) == 0L) {
    return(invisible(0L))
  }
  if (!is_sqlite_active_for_inputs(cache_path, demo = demo)) {
    return(invisible(0L))
  }

  conn <- connect_packera_db(read_only = FALSE, demo = FALSE)
  if (is.null(conn)) {
    return(invisible(0L))
  }
  on.exit(disconnect_packera_db(conn), add = TRUE)
  if (!db_table_exists_safe(conn, "environmental_cache")) {
    return(invisible(0L))
  }

  env_df <- as.data.frame(df, stringsAsFactors = FALSE)
  for (bc in intersect(c("ssurgo_lithic_series_pct", PACKERA_DB_WC_BIO_COLS), names(env_df))) {
    v_bc <- env_df[[bc]]
    env_df[[bc]] <- if (is.factor(v_bc)) suppressWarnings(as.numeric(as.character(v_bc))) else suppressWarnings(as.numeric(v_bc))
  }
  if (!("latitude" %in% names(env_df)) && "decimalLatitude" %in% names(env_df)) {
    env_df$latitude <- suppressWarnings(as.numeric(env_df$decimalLatitude))
  }
  if (!("longitude" %in% names(env_df)) && "decimalLongitude" %in% names(env_df)) {
    env_df$longitude <- suppressWarnings(as.numeric(env_df$decimalLongitude))
  }
  lat_num <- suppressWarnings(as.numeric(env_df$latitude))
  lon_num <- suppressWarnings(as.numeric(env_df$longitude))
  valid_xy <- !is.na(lat_num) & !is.na(lon_num) & is.finite(lat_num) & is.finite(lon_num)

  if (!("lat_lon_key" %in% names(env_df))) {
    env_df$lat_lon_key <- NA_character_
  }
  missing_key <- is.na(env_df$lat_lon_key) | !nzchar(trimws(as.character(env_df$lat_lon_key)))
  fill_idx <- which(missing_key & valid_xy)
  if (length(fill_idx) > 0L) {
    env_df$lat_lon_key[fill_idx] <- sprintf("%.5f_%.5f", lat_num[fill_idx], lon_num[fill_idx])
  }
  env_df <- env_df[!is.na(env_df$lat_lon_key) & nzchar(trimws(as.character(env_df$lat_lon_key))), , drop = FALSE]
  if (nrow(env_df) == 0L) {
    return(invisible(0L))
  }

  n_up <- db_upsert_dataframe(
    conn,
    table_name = "environmental_cache",
    df = env_df,
    pk_cols = c("lat_lon_key"),
    preserve_existing_on_null = TRUE
  )
  invisible(n_up)
}

persist_step05_to_sqlite <- function(
  vouchers_df,
  vouchers_path = NULL,
  run_id = NULL,
  demo = FALSE
) {
  if (isTRUE(demo) || is.null(vouchers_df) || nrow(vouchers_df) == 0L) {
    return(invisible(FALSE))
  }
  if (!is_sqlite_active_for_inputs(vouchers_path, demo = demo)) {
    return(invisible(FALSE))
  }

  srf_df <- as.data.frame(vouchers_df, stringsAsFactors = FALSE)
  upsert_environmental_cache_sqlite(srf_df, cache_path = vouchers_path, demo = demo)

  conn <- connect_packera_db(read_only = FALSE, demo = FALSE)
  if (is.null(conn)) {
    return(invisible(FALSE))
  }
  on.exit(disconnect_packera_db(conn), add = TRUE)

  rid <- resolve_packera_run_id(run_id)
  ensure_pipeline_run_r(conn, run_id = rid)

  srf_df$run_id <- rid
  if (!("plant_individual_id" %in% names(srf_df))) srf_df$plant_individual_id <- 0L
  srf_df$plant_individual_id <- suppressWarnings(as.integer(srf_df$plant_individual_id))
  srf_df$plant_individual_id[is.na(srf_df$plant_individual_id)] <- 0L

  for (pc in intersect(c("primary_env_predictors", "primary_clim_predictors", "primary_edaph_predictors"), names(srf_df))) {
    srf_df[[pc]] <- .normalize_semicolon_col_r(srf_df[[pc]])
  }
  for (sc in intersect(c("ecotype_regime", "nearest_type_niche_basionym"), names(srf_df))) {
    srf_df[[sc]] <- as.character(srf_df[[sc]])
  }
  for (nc in intersect(c(
    "edaphic_z_score", "climatic_z_score", "type_edaphic_z_score", "env_mahal_d2",
    "lithic_outcrop_index", "summer_moisture_deficit_idx", "maritime_continentality_idx",
    "hopkins_adjusted_doy", "hopkins_adjusted_doy_resid", "anthesis_dominance_weight",
    "allochronic_bimodality_delta_bic", "nearest_type_basionym_edaphic_z",
    PACKERA_DB_WC_BIO_COLS
  ), names(srf_df))) {
    v_nc <- srf_df[[nc]]
    srf_df[[nc]] <- if (is.factor(v_nc)) suppressWarnings(as.numeric(as.character(v_nc))) else suppressWarnings(as.numeric(v_nc))
  }
  if ("env_mcd_outlier_flag" %in% names(srf_df)) {
    v_mcd <- srf_df$env_mcd_outlier_flag
    if (is.logical(v_mcd)) {
      srf_df$env_mcd_outlier_flag <- as.integer(v_mcd)
    } else if (is.character(v_mcd) || is.factor(v_mcd)) {
      s_mcd <- trimws(as.character(v_mcd))
      out_mcd <- suppressWarnings(as.integer(s_mcd))
      u_mcd <- toupper(s_mcd)
      out_mcd[is.na(out_mcd) & u_mcd %in% c("TRUE", "T", "YES", "Y")] <- 1L
      out_mcd[is.na(out_mcd) & u_mcd %in% c("FALSE", "F", "NO", "N")] <- 0L
      srf_df$env_mcd_outlier_flag <- out_mcd
    } else {
      v_num <- suppressWarnings(as.numeric(v_mcd))
      out_mcd <- rep(NA_integer_, length(v_num))
      ok_mcd <- !is.na(v_num) & is.finite(v_num)
      out_mcd[ok_mcd] <- as.integer(v_num[ok_mcd] != 0)
      srf_df$env_mcd_outlier_flag <- out_mcd
    }
  }

  if (!("lat_lon_key" %in% names(srf_df))) {
    lat_src <- if ("latitude" %in% names(srf_df)) srf_df$latitude else if ("decimalLatitude" %in% names(srf_df)) srf_df$decimalLatitude else NULL
    lon_src <- if ("longitude" %in% names(srf_df)) srf_df$longitude else if ("decimalLongitude" %in% names(srf_df)) srf_df$decimalLongitude else NULL
    if (!is.null(lat_src) && !is.null(lon_src)) {
      lat_n <- suppressWarnings(as.numeric(lat_src))
      lon_n <- suppressWarnings(as.numeric(lon_src))
      ok_ll <- !is.na(lat_n) & !is.na(lon_n) & is.finite(lat_n) & is.finite(lon_n)
      srf_df$lat_lon_key <- NA_character_
      srf_df$lat_lon_key[ok_ll] <- sprintf("%.5f_%.5f", lat_n[ok_ll], lon_n[ok_ll])
    }
  }

  db_upsert_dataframe(
    conn,
    table_name = "spatial_rf_predictions",
    df = srf_df,
    pk_cols = c("run_id", "catalogNumber", "plant_individual_id"),
    preserve_existing_on_null = FALSE
  )
  .prune_stale_run_rows_r(conn, "spatial_rf_predictions", rid, srf_df)
  invisible(TRUE)
}
