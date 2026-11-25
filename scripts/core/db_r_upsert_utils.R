if (grDevices::dev.cur() == 1L) grDevices::pdf(NULL)
# ==============================================================================
# Module: scripts/core/db_r_upsert_utils.R
# SQLite table checks, pipeline_runs registration, logical coercion, and batched upserts.
# ==============================================================================

db_table_exists_safe <- function(conn, table_name) {
  if (is.null(conn)) {
    return(FALSE)
  }
  isTRUE(tryCatch(DBI::dbExistsTable(conn, table_name), error = function(e) FALSE))
}

ensure_pipeline_run_r <- function(conn, run_id = resolve_packera_run_id(), focal_taxon = NULL, intraspecific_mode = NULL) {
  if (is.null(conn) || !db_table_exists_safe(conn, "pipeline_runs")) {
    return(invisible(FALSE))
  }
  rid <- resolve_packera_run_id(run_id)
  now_iso <- format(as.POSIXct(Sys.time(), tz = "UTC"), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")
  rdir <- .get_env_run_dir_r()
  if (!nzchar(rdir)) rdir <- NA_character_
  ft <- if (!is.null(focal_taxon) && nzchar(trimws(focal_taxon))) trimws(focal_taxon) else NA_character_
  im <- if (is.null(intraspecific_mode) || is.na(intraspecific_mode)) NA_integer_ else if (isTRUE(intraspecific_mode)) 1L else 0L

  tryCatch(
    DBI::dbExecute(
      conn,
      paste0(
        "INSERT INTO pipeline_runs (run_id, run_dir, start_timestamp, focal_taxon, intraspecific_mode) ",
        "VALUES ($1, $2, $3, $4, COALESCE($5, 0)) ",
        "ON CONFLICT(run_id) DO UPDATE SET ",
        "run_dir = COALESCE(excluded.run_dir, pipeline_runs.run_dir), ",
        "focal_taxon = COALESCE(excluded.focal_taxon, pipeline_runs.focal_taxon), ",
        "intraspecific_mode = COALESCE($5, pipeline_runs.intraspecific_mode);"
      ),
      params = list(rid, rdir, now_iso, ft, im)
    ),
    error = function(e) NULL
  )
  invisible(TRUE)
}

coerce_sqlite_logical_cols <- function(df, cols = PACKERA_DB_LOGICAL_COLS) {
  if (is.null(df) || nrow(df) == 0L) {
    return(df)
  }
  hit_cols <- intersect(cols, names(df))
  for (cn in hit_cols) {
    v <- df[[cn]]
    if (is.logical(v)) next
    if (is.numeric(v) || is.integer(v)) {
      df[[cn]] <- ifelse(is.na(v), NA, v != 0)
    } else {
      s <- toupper(trimws(as.character(v)))
      out <- rep(NA, length(s))
      out[s %in% c("TRUE", "T", "1", "1.0", "YES")] <- TRUE
      out[s %in% c("FALSE", "F", "0", "0.0", "NO")] <- FALSE
      df[[cn]] <- out
    }
  }
  df
}

db_upsert_dataframe <- function(
  conn,
  table_name,
  df,
  pk_cols,
  preserve_existing_on_null = TRUE
) {
  if (is.null(conn) || is.null(df) || nrow(df) == 0L) {
    return(0L)
  }
  if (!db_table_exists_safe(conn, table_name)) {
    return(0L)
  }

  .migrate_dynamic_sqlite_table_r(conn, table_name)
  tbl_fields <- DBI::dbListFields(conn, table_name)
  keep_cols <- intersect(names(df), tbl_fields)
  if (!all(pk_cols %in% keep_cols)) {
    stop(sprintf(
      "db_upsert_dataframe: missing primary key columns (%s) for table '%s'",
      paste(setdiff(pk_cols, keep_cols), collapse = ", "), table_name
    ))
  }

  sub_df <- as.data.frame(df[, keep_cols, drop = FALSE], stringsAsFactors = FALSE)
  for (cn in keep_cols) {
    col_v <- sub_df[[cn]]
    if (is.list(col_v)) {
      col_v <- .normalize_semicolon_col_r(col_v)
      sub_df[[cn]] <- col_v
    }
    if (is.factor(col_v)) {
      sub_df[[cn]] <- as.character(col_v)
    } else if (is.logical(col_v)) {
      sub_df[[cn]] <- ifelse(is.na(col_v), NA_integer_, ifelse(col_v, 1L, 0L))
    } else if (is.numeric(col_v)) {
      bad_num <- is.nan(col_v) | is.infinite(col_v)
      if (any(bad_num)) col_v[bad_num] <- NA_real_
      sub_df[[cn]] <- col_v
    } else if (is.character(col_v)) {
      trimmed <- trimws(col_v)
      trimmed[!nzchar(trimmed) | tolower(trimmed) %in% c("na", "nan", "null")] <- NA_character_
      sub_df[[cn]] <- trimmed
    }
  }

  # Deduplicate on pk_cols keeping the last occurrence so transactional batch binds never collide
  pk_sig <- do.call(paste, c(lapply(pk_cols, function(k) as.character(sub_df[[k]])), sep = "\r"))
  sub_df <- sub_df[!duplicated(pk_sig, fromLast = TRUE), , drop = FALSE]

  non_pk <- setdiff(keep_cols, pk_cols)
  q_table <- DBI::dbQuoteIdentifier(conn, table_name)
  q_cols <- vapply(keep_cols, function(x) as.character(DBI::dbQuoteIdentifier(conn, x)), character(1))
  q_pks <- vapply(pk_cols, function(x) as.character(DBI::dbQuoteIdentifier(conn, x)), character(1))
  placeholders <- paste0("$", seq_along(keep_cols), collapse = ", ")

  if (length(non_pk) > 0L) {
    q_non_pk <- vapply(non_pk, function(x) as.character(DBI::dbQuoteIdentifier(conn, x)), character(1))
    if (isTRUE(preserve_existing_on_null)) {
      upd_clause <- paste(
        sprintf("%s = COALESCE(excluded.%s, %s.%s)", q_non_pk, q_non_pk, q_table, q_non_pk),
        collapse = ", "
      )
    } else {
      upd_clause <- paste(sprintf("%s = excluded.%s", q_non_pk, q_non_pk), collapse = ", ")
    }
    sql <- sprintf(
      "INSERT INTO %s (%s) VALUES (%s) ON CONFLICT (%s) DO UPDATE SET %s;",
      q_table,
      paste(q_cols, collapse = ", "),
      placeholders,
      paste(q_pks, collapse = ", "),
      upd_clause
    )
  } else {
    sql <- sprintf(
      "INSERT OR IGNORE INTO %s (%s) VALUES (%s);",
      q_table,
      paste(q_cols, collapse = ", "),
      placeholders
    )
  }

  n_sub <- nrow(sub_df)
  chunk_sz <- 5000L
  DBI::dbWithTransaction(conn, {
    if (n_sub <= chunk_sz) {
      DBI::dbExecute(conn, sql, params = unname(as.list(sub_df)))
    } else {
      starts <- seq.int(1L, n_sub, by = chunk_sz)
      for (s_idx in starts) {
        e_idx <- min(n_sub, s_idx + chunk_sz - 1L)
        DBI::dbExecute(conn, sql, params = unname(as.list(sub_df[s_idx:e_idx, , drop = FALSE])))
      }
    }
  })
  n_sub
}
