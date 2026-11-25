if (grDevices::dev.cur() == 1L) grDevices::pdf(NULL)
# [WORKER1_MODULAR_FACADE]
# ==============================================================================
# Script: scripts/morphometrics/authority_weighting_utils.R
# Unified Empirical Bayes / Dawid-Skene Authority Weighting & Adaptive Anchor
# Facade (auto-sources bounded sibling modules in deterministic load order).
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
.auth_dir <- .resolve_self_dir("scripts/morphometrics", "auth_concept_type_masks.R")

if (!exists("prune_collinear_dissection_scalars", mode = "function")) {
  for (.dl_p in unique(c("scripts/morphometrics/dissected_leaf_utils.R", file.path(.auth_dir, "dissected_leaf_utils.R")))) {
    if (!is.na(.dl_p) && nzchar(.dl_p) && file.exists(.dl_p)) {
      source(.dl_p, local = .load_env)
      break
    }
  }
}
if (!exists("packera_mclapply", mode = "function")) {
  for (.db_p in unique(c("scripts/core/db_utils.R", file.path(.auth_dir, "..", "core", "db_utils.R")))) {
    if (!is.na(.db_p) && nzchar(.db_p) && file.exists(.db_p)) {
      source(.db_p, local = .load_env)
      break
    }
  }
}

for (.auth_mod in c(
  "auth_concept_type_masks.R",
  "auth_subspace_mcd_cov.R",
  "auth_anchor_pool_helpers.R",
  "auth_adaptive_anchors.R",
  "auth_deme_ledoit_cov.R",
  "auth_shrinkage_cda.R",
  "auth_eb_likelihoods.R",
  "auth_eb_weights_helpers.R",
  "auth_eb_weights.R",
  "auth_gmm_mclust_fit.R",
  "auth_type_basionym_utils.R",
  "auth_type_shared_subspace.R",
  "auth_type_morph_control.R",
  "auth_osr_metrics.R",
  "auth_oof_cda_audit.R"
)) {
  source(file.path(.auth_dir, .auth_mod), local = .load_env)
}
