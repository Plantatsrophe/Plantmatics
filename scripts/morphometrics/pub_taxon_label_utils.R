if (grDevices::dev.cur() == 1L) grDevices::pdf(NULL)
# ==============================================================================
# Module: scripts/morphometrics/pub_taxon_label_utils.R
# Predictor labels, Packera synonym standardization, palettes, and plot-group resolution
# ==============================================================================

EXTENDED_COLORBLIND_HEX <- c(
  "#E69F00", "#56B4E9", "#F0E442", "#333333", "#882255",
  "#44AA99", "#999933", "#332288", "#117733", "#AA4499",
  "#DDCC77", "#661100", "#6699CC", "#D95F02", "#7570B3",
  "#E7298A", "#66A61E", "#E6AB02", "#A6761D", "#008080"
)

format_predictor_labels <- function(var_vec) {
  dict <- c(
    "fine_soil_ph"                     = "Rhizosphere Soil pH in H2O (0-15cm)",
    "ilr_texture_1"                    = "Soil Texture ILR-1: Sand vs Silt-Clay (0-15cm)",
    "ilr_texture_2"                    = "Soil Texture ILR-2: Silt vs Clay (0-15cm)",
    "fine_soil_sand"                   = "Soil Sand Fraction (%, 0-15cm)",
    "fine_soil_clay"                   = "Soil Clay Fraction (%, 0-15cm)",
    "fine_soil_cec"                    = "Rhizosphere Soil CEC (cmol+/kg, 0-15cm)",
    "fine_cfvo_pct"                    = "Coarse Fragment Volume (%, 0-15cm)",
    "fine_bedrock_depth_cm"            = "Depth to Lithic Bedrock (cm)",
    "is_azonal_lithology"              = "USGS SGMC Azonal Lithology Indicator",
    "ssurgo_rock_outcrop_pct"          = "SSURGO Rock Outcrop Component (%)",
    "ssurgo_bedrock_depth_cm"          = "SSURGO Lithic Bedrock Depth (cm)",
    "ssurgo_aws_0_50cm"                = "Available Water Storage 0-50cm (cm)",
    "topo_slope_deg"                   = "30m Terrain Slope (deg)",
    "topo_tpi_30m"                     = "30m Topographic Position Index (TPI, m)",
    "topo_hli_30m"                     = "30m Heat Load Index (HLI, McCune & Keon 2002)",
    "topo_twi_30m"                     = "30m D8 Topographic Wetness Index (TWI)",
    "soil_ph"                          = "SoilGrids 250m Soil pH (H2O)",
    "soil_sand"                        = "SoilGrids 250m Sand Fraction (%)",
    "soil_clay"                        = "SoilGrids 250m Clay Fraction (%)",
    "soil_silt"                        = "SoilGrids 250m Silt Fraction (%)",
    "soil_cec"                         = "SoilGrids 250m Cation Exchange (cmol+/kg)",
    "soil_cfvo"                        = "SoilGrids 250m Coarse Fragments (%)",
    "soil_bulk_density"                = "SoilGrids 250m Bulk Density (g/cm3)",
    "bio1_temp_mean"                   = "Mean Annual Temp (BIO1, deg C)",
    "bio2_mean_diurnal_range"          = "Mean Diurnal Temp Range (BIO2, deg C)",
    "bio3_isothermality"               = "Isothermality (BIO3, BIO2/BIO7 x 100)",
    "bio4_temp_seasonality"            = "Temp Seasonality (BIO4, SD x 100)",
    "bio5_max_temp_warmest_month"      = "Max Temp Warmest Month (BIO5, deg C)",
    "bio6_min_temp_coldest_month"      = "Min Temp Coldest Month (BIO6, deg C)",
    "bio7_temp_annual_range"           = "Temp Annual Range (BIO7, deg C)",
    "bio8_mean_temp_wettest_qtr"       = "Mean Temp Wettest Quarter (BIO8, deg C)",
    "bio9_mean_temp_driest_qtr"        = "Mean Temp Driest Quarter (BIO9, deg C)",
    "bio10_mean_temp_warmest_qtr"      = "Mean Temp Warmest Quarter (BIO10, deg C)",
    "bio11_mean_temp_coldest_qtr"      = "Mean Temp Coldest Quarter (BIO11, deg C)",
    "bio12_precip_annual"              = "Annual Precipitation (BIO12, mm)",
    "bio13_precip_wettest_month"       = "Wettest-Month Precip (BIO13, mm)",
    "bio14_precip_driest_month"        = "Driest-Month Precip (BIO14, mm)",
    "bio15_precip_seasonality"         = "Precipitation Seasonality (BIO15, CV)",
    "bio16_precip_wettest_qtr"         = "Wettest-Quarter Precip (BIO16, mm)",
    "bio17_precip_driest_qtr"          = "Driest-Quarter Precip (BIO17, mm)",
    "bio18_precip_warmest_qtr"         = "Warmest-Quarter Precip (BIO18, mm)",
    "bio19_precip_coldest_qtr"         = "Coldest-Quarter Precip (BIO19, mm)",
    "lithic_outcrop_index"             = "Lithic Outcrop Exposure Index (0-1)",
    "summer_moisture_deficit_idx"      = "Summer Moisture Deficit Index",
    "maritime_continentality_idx"      = "Maritime Continentality Index",
    "min_bedrock_depth_buffer_cm"      = "Min Buffered Bedrock Depth (cm)",
    "ph_q10_buffer"                    = "Buffered Soil pH 10th Pct",
    "ph_q90_buffer"                    = "Buffered Soil pH 90th Pct",
    "max_cfvo_buffer_pct"              = "Max Buffered Coarse Fragments (%)",
    "azonal_lithology_presence_buffer" = "Buffered Azonal Lithology Presence",
    "fine_soil_soc"                    = "Soil Organic Carbon (g/kg, 0-15cm)",
    "topo_elev_m"                      = "30m Elevation (m a.s.l.)",
    "delta_doy"                        = "Hopkins Phenological Anomaly (Delta DOY)",
    "doy"                              = "Flowering Day of Year (DOY)",
    "soil_ca_mg_molar_ratio"           = "Soil Ca:Mg Molar Ratio",
    "soil_ultramafic_ni_cr_index"      = "Soil Ultramafic Ni-Cr Index",
    "soil_carbonate_c_pct"             = "Soil Carbonate Carbon (%)",
    "nwi_wetland_dist_m"               = "NWI Wetland Distance (m)",
    "ultramafic_calcareous_gradient"   = "Ultramafic-Calcareous Geochem Gradient",
    "wetland_seepage_index"            = "Wetland & Hydrologic Seepage Index"
  )
  vapply(as.character(var_vec), function(v) {
    if (is.na(v) || !nzchar(v)) return("Unknown Predictor")
    if (v %in% names(dict)) return(unname(dict[[v]]))
    if (grepl("^MEM[0-9]+$", v)) return(sub("^MEM", "Spatial Eigenvector MEM-", v))
    gsub("_", " ", v)
  }, character(1), USE.NAMES = FALSE)
}

.pub_figures_script_dir <- tryCatch({
  ofiles <- unlist(lapply(sys.frames(), function(fr) fr$ofile))
  cand <- if (length(ofiles) > 0L) tail(ofiles, 1L) else sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))[1L]
  if (!is.na(cand) && nzchar(cand)) dirname(normalizePath(cand, mustWork = FALSE)) else ""
}, error = function(e) "")

.packera_taxon_syn_env <- new.env(parent = emptyenv())
.packera_taxon_syn_env$db_sig <- ""
.packera_taxon_syn_env$syn_map <- character(0)
.packera_taxon_syn_env$py_memo <- character(0)

.expand_synonym_binomial_keys <- function(syn_vec) {
  if (length(syn_vec) == 0L) return(character(0))
  nms <- tolower(gsub("[[:space:]]+", " ", trimws(names(syn_vec))))
  vals <- gsub("^Packera[[:space:]]+(x|×)[[:space:]]+", "Packera ", trimws(as.character(syn_vec)), ignore.case = TRUE)
  keep <- !is.na(nms) & nzchar(nms) & !is.na(vals) & nzchar(vals)
  if (!any(keep)) return(character(0))
  out <- stats::setNames(vals[keep], nms[keep])
  bin_pat <- "^([a-z]+[[:space:]]+((x|×)[[:space:]]*)?[a-z][a-z-]*)([^A-Za-z0-9_].*|$)"
  is_bin <- grepl(bin_pat, names(out))
  if (any(is_bin)) {
    b_keys <- sub(bin_pat, "\\1", names(out)[is_bin])
    b_vals <- unname(out[is_bin])
    new_m <- !(b_keys %in% names(out)) & !duplicated(b_keys)
    if (any(new_m)) out <- c(out, stats::setNames(b_vals[new_m], b_keys[new_m]))
  }
  out
}

.get_packera_synonym_map <- function() {
  cand_dbs <- character(0)
  env_db <- trimws(Sys.getenv("PLANTMATICS_DB_PATH", unset = ""))
  if (!nzchar(env_db)) env_db <- trimws(Sys.getenv("PACKERA_DB_PATH", unset = ""))
  if (nzchar(env_db) && env_db != ":memory:") cand_dbs <- c(cand_dbs, env_db)
  if (exists("resolve_plantmatics_db_path", mode = "function")) {
    rdb <- tryCatch(as.character(resolve_plantmatics_db_path()), error = function(e) "")
    if (!is.na(rdb) && nzchar(rdb) && rdb != ":memory:") cand_dbs <- c(cand_dbs, rdb)
  } else if (exists("resolve_packera_db_path", mode = "function")) {
    rdb <- tryCatch(as.character(resolve_packera_db_path()), error = function(e) "")
    if (!is.na(rdb) && nzchar(rdb) && rdb != ":memory:") cand_dbs <- c(cand_dbs, rdb)
  }
  cand_dbs <- c(cand_dbs, "data/plantmatics.sqlite")
  if (!file.exists("data/plantmatics.sqlite") && file.exists("data/packera_pipeline.sqlite")) {
    cand_dbs <- c(cand_dbs, "data/packera_pipeline.sqlite")
  }
  if (nzchar(.pub_figures_script_dir)) {
    rel_new <- file.path(.pub_figures_script_dir, "..", "..", "data", "plantmatics.sqlite")
    rel_old <- file.path(.pub_figures_script_dir, "..", "..", "data", "packera_pipeline.sqlite")
    cand_dbs <- c(cand_dbs, rel_new)
    if (!file.exists(rel_new) && file.exists(rel_old)) {
      cand_dbs <- c(cand_dbs, rel_old)
    }
  }
  cand_dbs <- unique(normalizePath(cand_dbs[file.exists(cand_dbs)], mustWork = FALSE))
  sig <- if (length(cand_dbs) == 0L) "none" else paste(cand_dbs, as.numeric(file.info(cand_dbs)$mtime), sep = "@", collapse = "|")
  if (identical(sig, .packera_taxon_syn_env$db_sig)) return(.packera_taxon_syn_env$syn_map)

  merged <- character(0)
  if (length(cand_dbs) > 0L && requireNamespace("DBI", quietly = TRUE) && requireNamespace("RSQLite", quietly = TRUE)) {
    for (dbp in cand_dbs) {
      res <- tryCatch({
        con <- DBI::dbConnect(RSQLite::SQLite(), dbname = dbp, flags = RSQLite::SQLITE_RO)
        on.exit(tryCatch(DBI::dbDisconnect(con), error = function(e) NULL), add = TRUE)
        if ("taxon_synonyms" %in% DBI::dbListTables(con)) {
          DBI::dbGetQuery(con, "SELECT raw_name_lower, canonical_packera FROM taxon_synonyms;")
        } else NULL
      }, error = function(e) NULL)
      if (!is.null(res) && nrow(res) > 0L) {
        cur <- .expand_synonym_binomial_keys(stats::setNames(as.character(res$canonical_packera), as.character(res$raw_name_lower)))
        add_nms <- setdiff(names(cur), names(merged))
        if (length(add_nms) > 0L) merged <- c(merged, cur[add_nms])
      }
    }
  }
  .packera_taxon_syn_env$db_sig <- sig
  .packera_taxon_syn_env$syn_map <- merged
  merged
}

.resolve_taxa_via_harvester_py <- function(query_vec) {
  if (length(query_vec) == 0L) return(character(0))
  u_q <- unique(as.character(query_vec))
  memo <- .packera_taxon_syn_env$py_memo
  miss <- setdiff(u_q, names(memo))
  if (length(miss) > 0L) {
    ws_root <- if (nzchar(.pub_figures_script_dir)) normalizePath(file.path(.pub_figures_script_dir, "..", ".."), mustWork = FALSE) else normalizePath(".", mustWork = FALSE)
    py_cands <- c(file.path(ws_root, ".venv", "bin", "python"), file.path(".venv", "bin", "python"), Sys.which("python3"))
    py_bin <- py_cands[nzchar(py_cands) & file.exists(py_cands)][1L]
    resolved_miss <- stats::setNames(rep(NA_character_, length(miss)), miss)
    if (!is.na(py_bin) && nzchar(py_bin) && file.exists(file.path(ws_root, "scripts", "core", "harvester.py"))) {
      tf_in <- tempfile("packera_tax_in_"); on.exit(unlink(tf_in), add = TRUE)
      writeLines(miss, con = tf_in, useBytes = TRUE)
      py_code <- sprintf(
        "import sys; sys.path.insert(0, %s); from scripts.core.harvester import standardize_packera_taxon as f\nwith open(%s, 'r', encoding='utf-8') as fh:\n for ln in fh:\n  r = f(ln.rstrip('\\r\\n'))\n  print(r if r else '')",
        shQuote(ws_root), shQuote(tf_in)
      )
      py_out <- tryCatch(suppressWarnings(system2(py_bin, args = c("-c", shQuote(py_code)), stdout = TRUE, stderr = FALSE)), error = function(e) character(0))
      if (length(py_out) == length(miss)) {
        clean_out <- trimws(as.character(py_out))
        resolved_miss <- stats::setNames(ifelse(nzchar(clean_out), clean_out, NA_character_), miss)
      }
    }
    memo <- c(memo, resolved_miss)
    .packera_taxon_syn_env$py_memo <- memo
  }
  unname(memo[as.character(query_vec)])
}

standardize_taxon_name <- function(raw_vec) {
  if (is.null(raw_vec) || length(raw_vec) == 0L) return(character(0))
  chr_vec <- as.character(raw_vec)
  u_vals <- unique(chr_vec)
  syn_map <- .get_packera_synonym_map()
  has_syn <- length(syn_map) > 0L
  need_py_idx <- integer(0)
  fb_packera <- rep(NA_character_, length(u_vals))
  open_ep <- c("sp", "spp", "cf", "aff", "indet", "unknown", "na")
  syn_gen_pat <- "^(packera|p\\.|senecio|s\\.|cineraria|jacobaea|aureolaria)[[:space:]]*(cf\\.?[[:space:]]+|aff\\.?[[:space:]]+)?(x[[:space:]]+|×[[:space:]]*)?([a-z][a-z-]*).*$"

  std_u <- character(length(u_vals))
  for (i in seq_along(u_vals)) {
    val <- u_vals[i]
    if (is.na(val)) { std_u[i] <- NA_character_; next }
    s_orig <- trimws(val)
    s <- tolower(gsub("[[:space:]]+", " ", s_orig))
    if (!nzchar(s_orig) || s %in% c("none", "nan", "null", "na", "unknown", "unassigned", "unassigned_packera", "packera")) {
      std_u[i] <- NA_character_; next
    }
    if (grepl("^Morphotype_K[0-9]+$|^Ecotype_|^LatTertile_|_Coastal$|_Outcrop$|\\[K[0-9]+\\]$", s_orig, ignore.case = TRUE)) {
      std_u[i] <- s_orig; next
    }
    if (has_syn && (s %in% names(syn_map))) {
      std_u[i] <- unname(syn_map[[s]]); next
    }
    if (has_syn && grepl("^[a-z]+[[:space:]]+[a-z-]+[[:space:]]+(var\\.?|subsp\\.?|ssp\\.?|f\\.?)[[:space:]]+[a-z-]+", s)) {
      v_raw <- sub("^([a-z]+[[:space:]]+[a-z-]+[[:space:]]+(var\\.?|subsp\\.?|ssp\\.?|f\\.?)[[:space:]]+[a-z-]+).*$", "\\1", s)
      v_dot <- gsub("(^|[^A-Za-z0-9_])(var|subsp|ssp|f)\\.?[[:space:]]+", "\\1\\2. ", v_raw)
      if (v_raw %in% names(syn_map)) { std_u[i] <- unname(syn_map[[v_raw]]); next }
      if (v_dot %in% names(syn_map)) { std_u[i] <- unname(syn_map[[v_dot]]); next }
    }
    if (grepl(syn_gen_pat, s)) {
      gen_tok <- sub(syn_gen_pat, "\\1", s)
      hyb_tok <- trimws(sub(syn_gen_pat, "\\3", s))
      ep <- sub("^[×-]+|[-]+$", "", sub(syn_gen_pat, "\\4", s))
      if (!nzchar(ep) || nchar(ep) < 2L || (ep %in% open_ep)) { std_u[i] <- NA_character_; next }
      gen_full <- if (gen_tok == "p.") "packera" else if (gen_tok == "s.") "senecio" else gen_tok
      is_hyb <- nzchar(hyb_tok)
      b_plain <- paste(gen_full, ep)
      b_hyb <- if (is_hyb) paste(gen_full, "×", ep) else b_plain
      if (has_syn && (b_hyb %in% names(syn_map))) { std_u[i] <- unname(syn_map[[b_hyb]]); next }
      if (has_syn && (b_plain %in% names(syn_map))) { std_u[i] <- unname(syn_map[[b_plain]]); next }
      if (gen_full == "packera") {
        p_can <- paste0("Packera ", ep)
        if (has_syn) { std_u[i] <- p_can; next }
        fb_packera[i] <- p_can
      }
      need_py_idx <- c(need_py_idx, i)
      std_u[i] <- NA_character_
      next
    }
    std_u[i] <- NA_character_
  }

  if (length(need_py_idx) > 0L) {
    py_res <- .resolve_taxa_via_harvester_py(u_vals[need_py_idx])
    std_u[need_py_idx] <- ifelse(!is.na(py_res) & nzchar(py_res), py_res, fb_packera[need_py_idx])
  }
  std_u[match(chr_vec, u_vals)]
}
standardize_packera_taxon <- function(species_vec, na_value = "Unknown") {
  res <- standardize_taxon_name(species_vec)
  if (!is.null(na_value)) {
    res[is.na(res)] <- na_value
  }
  res
}

.order_unique_packera_taxa <- function(taxa_vec, preserve_order = FALSE) {
  clean_v <- trimws(as.character(taxa_vec[!is.na(taxa_vec) & nzchar(trimws(as.character(taxa_vec)))]))
  if (length(clean_v) == 0L) return(character(0))
  u_raw <- unique(clean_v)
  u_non_unk <- setdiff(u_raw, "Unknown")
  if (length(u_non_unk) == 0L) return(character(0))
  if (isTRUE(preserve_order)) {
    return(u_non_unk)
  }
  sort(u_non_unk)
}

get_dynamic_packera_palette <- function(taxa_vec = NULL, preserve_order = FALSE) {
  if (is.null(taxa_vec) || length(taxa_vec) == 0L) return(c("Unknown" = "#7F7F7F"))
  u_non_unk <- .order_unique_packera_taxa(taxa_vec, preserve_order = preserve_order)
  if (length(u_non_unk) == 0L) return(c("Unknown" = "#7F7F7F"))
  pal_pool <- c("#0072B2", "#D55E00", "#009E73", "#CC79A7", EXTENDED_COLORBLIND_HEX)
  if (length(u_non_unk) > length(pal_pool)) {
    more_n <- length(u_non_unk) - length(pal_pool)
    hcl_ext <- toupper(grDevices::hcl.colors(more_n + 8L, palette = "Dark 3"))
    hcl_ext <- setdiff(hcl_ext, c(toupper(pal_pool), "#7F7F7F"))
    pal_pool <- c(pal_pool, hcl_ext)
  }
  cols <- rep_len(pal_pool, length(u_non_unk))
  names(cols) <- u_non_unk
  c(cols, "Unknown" = "#7F7F7F")
}

resolve_plot_group_column <- function(df, cda_res = NULL, gmm_res = NULL) {
  if (is.null(df) || nrow(df) == 0L) return(character(0))
  n <- nrow(df)
  raw_sp <- if ("species_standardized" %in% names(df)) {
    as.character(df$species_standardized)
  } else if ("taxon" %in% names(df)) {
    as.character(df$taxon)
  } else if ("scientificName" %in% names(df)) {
    as.character(df$scientificName)
  } else if ("species_raw" %in% names(df)) {
    as.character(df$species_raw)
  } else {
    rep(NA_character_, n)
  }
  std_sp <- standardize_taxon_name(raw_sp)
  cand_sp <- ifelse(!is.na(std_sp), std_sp, ifelse(!is.na(raw_sp) & nzchar(trimws(raw_sp)), trimws(raw_sp), NA_character_))
  valid_sp <- cand_sp[!is.na(cand_sp) & nzchar(cand_sp) & !(tolower(cand_sp) %in% c("unknown", "na", "none"))]
  fallback_sp <- if (length(valid_sp) > 0L) valid_sp[1L] else "Unassigned_Packera"
  sp_vec <- ifelse(!is.na(cand_sp) & nzchar(cand_sp), cand_sp, fallback_sp)
  sp_valid <- unique(sp_vec[!is.na(sp_vec) & nzchar(sp_vec) & sp_vec != "Unknown"])

  clean_cand <- function(v) {
    if (is.null(v) || length(v) != n) return(NULL)
    vc <- trimws(as.character(v))
    vc[is.na(vc) | !nzchar(vc)] <- NA_character_
    vc
  }
  has_valid_multi <- function(vc) {
    !is.null(vc) && length(unique(vc[!is.na(vc) & vc != "Unknown"])) >= 2L
  }
  is_subcluster_vec <- function(vc) {
    if (!has_valid_multi(vc)) return(FALSE)
    vv <- unique(vc[!is.na(vc) & vc != "Unknown"])
    any(grepl("Morphotype_K|\\[K[0-9]+\\]|Ecotype|Subcluster|Region_", vv, ignore.case = TRUE)) ||
      !all(grepl("^Packera[[:space:]]+", vv))
  }

  cand_cda <- if (!is.null(cda_res) && !is.null(cda_res$grouping_col)) clean_cand(cda_res$grouping_col) else NULL
  cand_grp <- if ("grouping_col" %in% names(df)) clean_cand(df$grouping_col) else NULL
  cand_ana <- if ("analysis_group" %in% names(df)) clean_cand(df$analysis_group) else NULL
  cand_gmm <- NULL
  if ("gmm_cluster" %in% names(df) && length(unique(stats::na.omit(df$gmm_cluster))) >= 2L) {
    cand_gmm <- paste0("Morphotype_K", ifelse(is.na(df$gmm_cluster), 1L, as.integer(df$gmm_cluster)))
  } else if (!is.null(gmm_res) && !is.null(gmm_res$classification) && length(gmm_res$classification) == n &&
             length(unique(stats::na.omit(gmm_res$classification))) >= 2L) {
    cand_gmm <- paste0("Morphotype_K", ifelse(is.na(gmm_res$classification), 1L, as.integer(gmm_res$classification)))
  }

  use_subclusters <- (length(sp_valid) < 2L) ||
    is_subcluster_vec(cand_cda) ||
    is_subcluster_vec(cand_grp) ||
    is_subcluster_vec(cand_ana)

  if (use_subclusters) {
    for (cand in list(cand_cda, cand_grp, cand_ana, cand_gmm)) {
      if (has_valid_multi(cand)) {
        if (length(sp_valid) >= 2L && !is_subcluster_vec(cand)) next
        def_lbl <- stats::na.omit(cand)[1]
        return(ifelse(!is.na(cand), cand, def_lbl))
      }
    }
    for (col_fb in c("cda_predicted_taxon", "edaphic_best_fit_taxon", "recommended_determination")) {
      if (col_fb %in% names(df)) {
        cf <- clean_cand(df[[col_fb]])
        if (is_subcluster_vec(cf)) {
          def_lbl <- stats::na.omit(cf)[1]
          return(ifelse(!is.na(cf), cf, def_lbl))
        }
      }
    }
    if ("gmm_cluster" %in% names(df) && length(unique(stats::na.omit(df$gmm_cluster))) >= 2L) {
      return(paste0("Morphotype_K", ifelse(is.na(df$gmm_cluster), 1L, as.integer(df$gmm_cluster))))
    }
  }
  sp_vec
}
