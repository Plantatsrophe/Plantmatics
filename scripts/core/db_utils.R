if (grDevices::dev.cur() == 1L) grDevices::pdf(NULL)
# ==============================================================================
# Module: scripts/core/db_utils.R
# Project: Plantmatics / Packera Species Delimitation & Morphometrics Pipeline
# Affiliation: University of North Carolina at Chapel Hill Herbarium (NCU)
#
# Description:
#   Cross-language SQLite (WAL mode) DBI facade for R pipeline stages (03-07).
#   Auto-sources modular `scripts/core/db_r_*.R` companions into the caller
#   environment, exposing `RSQLite::SQLite()` connection management, dynamic
#   schema migration, transactional upserts, Step 03 contour I/O, Step 04
#   morphometric restoration, Step 05 environmental cache and
#   `type_edaphic_z_score` persistence, Step 06 genomics persistence, and
#   Step 07 `v_integrative_triage_master` synthesis queries.
# ==============================================================================

PLANTMATICS_SCHEMA_USER_VERSION <- 20261127L
PACKERA_SCHEMA_USER_VERSION <- PLANTMATICS_SCHEMA_USER_VERSION

.db_r_env <- environment()
.db_r_resolve_dir <- function() {
  frame_ofiles <- unlist(lapply(sys.frames(), function(fr) fr$ofile))
  frame_dirs <- if (length(frame_ofiles) > 0L) rev(dirname(as.character(frame_ofiles))) else character(0)
  cmd_files <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  cmd_dirs <- if (length(cmd_files) > 0L) dirname(cmd_files) else character(0)
  rel_frame_core <- if (length(frame_dirs) > 0L) file.path(frame_dirs, "..", "core") else character(0)
  rel_cmd_core <- if (length(cmd_dirs) > 0L) file.path(cmd_dirs, "..", "core") else character(0)
  cands <- unique(c(
    frame_dirs,
    cmd_dirs,
    rel_frame_core,
    rel_cmd_core,
    "scripts/core",
    file.path(getwd(), "scripts", "core")
  ))
  cands <- cands[!is.na(cands) & nzchar(cands)]
  for (d in cands) {
    if (file.exists(file.path(d, "db_r_env_parallel.R"))) {
      return(d)
    }
  }
  stop("scripts/core/db_utils.R: unable to locate sibling db_r_env_parallel.R module.")
}

.db_r_dir <- .db_r_resolve_dir()
.db_r_modules <- c(
  "db_r_env_parallel.R",
  "db_r_schema_specs_a.R",
  "db_r_schema_specs_b.R",
  "db_r_path_conn.R",
  "db_r_upsert_utils.R",
  "db_r_step03_qc_fetch.R",
  "db_r_step03_contour_io.R",
  "db_r_step03_persist.R",
  "db_r_step04_restore_a.R",
  "db_r_step04_restore_b.R",
  "db_r_step04_persist.R",
  "db_r_step05_env_cache.R",
  "db_r_step06_07_master.R"
)
for (.db_r_f in .db_r_modules) {
  .db_r_path <- file.path(.db_r_dir, .db_r_f)
  if (!file.exists(.db_r_path)) {
    stop(sprintf("scripts/core/db_utils.R: missing required companion module '%s'.", .db_r_path))
  }
  source(.db_r_path, local = .db_r_env)
}
