if (grDevices::dev.cur() == 1L) grDevices::pdf(NULL)
# [WORKER1_MODULAR_FACADE]
# ==============================================================================
# Script: scripts/morphometrics/pub_ecoregion_manova_utils.R
# Facade for EPA Level III Ecoregion spatial lookup, Track 2 supervised CDA,
# weighted MANOVA/Hotelling T^2 audits, and Fig 04c plotting delegation.
# ==============================================================================

.resolve_self_dir <- function(default_rel, probe_file) {
  f_frm <- as.character(unlist(lapply(sys.frames(), function(fr) fr$ofile)))
  f_cmd <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  cand_d <- unique(c(rev(dirname(f_frm)), dirname(f_cmd), default_rel))
  cand_d <- cand_d[!is.na(cand_d) & nzchar(cand_d)]
  for (d in cand_d) {
    if (file.exists(file.path(d, probe_file))) return(d)
  }
  default_rel
}

.load_env <- environment()
.eco_manova_dir <- .resolve_self_dir("scripts/morphometrics", "ecoregion_epa_spatial.R")
for (.em_mod in c(
  "ecoregion_epa_spatial.R",
  "ecoregion_supervised_cda.R",
  "ecoregion_manova_stats.R",
  "ecoregion_manova_plot.R"
)) {
  source(file.path(.eco_manova_dir, .em_mod), local = .load_env)
}
