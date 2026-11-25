if (grDevices::dev.cur() == 1L) grDevices::pdf(NULL)
# [WORKER1_MODULAR_FACADE]
# ==============================================================================
# Module: scripts/morphometrics/dissected_leaf_utils.R
# Facade for positional zone routing, harmonic calibration, MCD/ECT, cauline
# track aggregation, collinearity pruning, and O(N*P) trace Euclidean PERMANOVA.
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
.dissect_dir <- .resolve_self_dir("scripts/morphometrics", "dissect_zone_and_harmonics.R")
for (.d_mod in c(
  "dissect_zone_and_harmonics.R",
  "dissect_mcd_and_ect.R",
  "dissect_cauline_track.R",
  "dissect_scalars_permanova.R"
)) {
  source(file.path(.dissect_dir, .d_mod), local = .load_env)
}
