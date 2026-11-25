"""Preload worker R loop, managed env keys, and stage dispatch maps for `scripts/tests/_r_session.py`."""

from __future__ import annotations

from typing import Dict

_MANAGED_ENV_KEYS = (
    "PLANTMATICS_DB_PATH",
    "PLANTMATICS_RUN_DIR",
    "PLANTMATICS_GENOMICS_DIR",
    "PLANTMATICS_SKIP_FIGURES",
    "PLANTMATICS_DEMO_MODE",
    "PLANTMATICS_MAX_WORKERS",
    "PLANTMATICS_IN_FORK_WORKER",
    "PLANTMATICS_CYTOLOGY_DIR",
    "PACKERA_DB_PATH",
    "PACKERA_RUN_DIR",
    "PACKERA_GENOMICS_DIR",
    "PACKERA_SKIP_FIGURES",
    "PACKERA_DEMO_MODE",
    "PACKERA_MAX_WORKERS",
    "PACKERA_CYTOLOGY_DIR",
    "PACKERA_SQLITE_DB",
    "RUN_ID",
)

_WORKER_R_CODE = r"""
local({
  options(rgl.useNULL = TRUE)
  boot_wd <- getwd()
  preload_script <- Sys.getenv("PACKERA_PRELOAD_R_SCRIPT", unset = "")
  Sys.unsetenv("PACKERA_PRELOAD_R_SCRIPT")

  if (nzchar(preload_script)) {
    boot_err <- NULL
    tryCatch(
      suppressMessages(source(preload_script, local = .GlobalEnv)),
      error = function(e) { boot_err <<- conditionMessage(e) }
    )
    if (!is.null(boot_err)) {
      cat("<<<PACKERA_BOOT_ERROR>>>:", boot_err, "\n", sep = "")
      flush(stdout())
      quit(save = "no", status = 2L)
    }
  }

  while (grDevices::dev.cur() > 1L) {
    tryCatch(grDevices::dev.off(), error = function(e) break)
  }
  if (grDevices::dev.cur() == 1L) {
    grDevices::pdf(NULL)
  }

  boot_snapshot <- as.list.environment(.GlobalEnv, all.names = TRUE)
  boot_names <- names(boot_snapshot)

  restore_global_state <- function() {
    while (grDevices::dev.cur() > 1L) {
      tryCatch(grDevices::dev.off(), error = function(e) break)
    }
    if (grDevices::dev.cur() == 1L) {
      grDevices::pdf(NULL)
    }
    tryCatch(setwd(boot_wd), error = function(e) NULL)
    cur_names <- ls(envir = .GlobalEnv, all.names = TRUE)
    extra_names <- setdiff(cur_names, boot_names)
    if (length(extra_names) > 0L) {
      rm(list = extra_names, envir = .GlobalEnv)
    }
    for (nm in boot_names) {
      if (!exists(nm, envir = .GlobalEnv, inherits = FALSE) ||
          !identical(get(nm, envir = .GlobalEnv, inherits = FALSE), boot_snapshot[[nm]])) {
        assign(nm, boot_snapshot[[nm]], envir = .GlobalEnv)
      }
    }
    if (exists(".packera_taxon_syn_env", envir = .GlobalEnv, inherits = FALSE)) {
      syn_e <- get(".packera_taxon_syn_env", envir = .GlobalEnv, inherits = FALSE)
      if (is.environment(syn_e)) {
        syn_e$db_sig <- ""
        syn_e$syn_map <- character(0)
        syn_e$py_memo <- character(0)
      }
    }
    if (exists(".AUTH_WEIGHT_CFG_ENV", envir = .GlobalEnv, inherits = FALSE)) {
      aw_e <- get(".AUTH_WEIGHT_CFG_ENV", envir = .GlobalEnv, inherits = FALSE)
      if (is.environment(aw_e)) {
        rm(list = ls(envir = aw_e, all.names = TRUE), envir = aw_e)
      }
    }
  }

  cat("<<<PACKERA_BOOT_READY>>>\n")
  flush(stdout())

  con_in <- file("stdin", open = "r")
  on.exit(tryCatch(close(con_in), error = function(e) NULL), add = TRUE)

  repeat {
    spec_line <- readLines(con_in, n = 1L, warn = FALSE)
    if (length(spec_line) == 0L) break
    spec_line <- trimws(spec_line)
    if (!nzchar(spec_line)) next
    if (identical(spec_line, "__QUIT__")) break

    hdr <- readLines(spec_line, warn = FALSE, encoding = "UTF-8")
    mode <- hdr[1L]
    entry_fn <- hdr[2L]
    req_cwd <- hdr[3L]
    script_path <- hdr[4L]
    args_file <- hdr[5L]
    env_keys_file <- hdr[6L]
    env_vals_file <- hdr[7L]
    out_file <- hdr[8L]
    err_file <- hdr[9L]
    status_file <- hdr[10L]
    done_token <- hdr[11L]

    cli_args <- if (nzchar(args_file) && file.exists(args_file)) {
      readLines(args_file, warn = FALSE, encoding = "UTF-8")
    } else {
      character(0)
    }
    env_keys <- if (nzchar(env_keys_file) && file.exists(env_keys_file)) {
      readLines(env_keys_file, warn = FALSE, encoding = "UTF-8")
    } else {
      character(0)
    }
    env_vals <- if (nzchar(env_vals_file) && file.exists(env_vals_file)) {
      readLines(env_vals_file, warn = FALSE, encoding = "UTF-8")
    } else {
      character(0)
    }

    prev_env <- stats::setNames(vector("list", length(env_keys)), env_keys)
    for (idx in seq_along(env_keys)) {
      k <- env_keys[idx]
      prev_env[[k]] <- Sys.getenv(k, unset = NA_character_)
      v <- if (idx <= length(env_vals)) env_vals[idx] else "__UNSET__"
      if (identical(v, "__UNSET__")) {
        Sys.unsetenv(k)
      } else {
        do.call(Sys.setenv, stats::setNames(list(v), k))
      }
    }

    err_chunks <- character(0)
    job_ok <- TRUE

    tryCatch(
      withCallingHandlers(
        utils::capture.output({
          if (nzchar(req_cwd)) setwd(req_cwd)
          if (identical(mode, "cli")) {
            if (identical(entry_fn, "main")) {
              main(cli_args)
            } else {
              fn_obj <- get(entry_fn, envir = .GlobalEnv, mode = "function")
              opts_obj <- parse_args_robust(cli_args)
              fn_obj(opts_obj)
            }
          } else {
            assign("commandArgs", function(trailingOnly = FALSE) {
              if (isTRUE(trailingOnly)) cli_args else base::commandArgs(FALSE)
            }, envir = .GlobalEnv)
            exprs <- parse(file = script_path, encoding = "UTF-8")
            if (length(exprs) > 0L) {
              eval(as.call(c(quote(`{`), as.list(exprs))), envir = .GlobalEnv)
            }
          }
          invisible(NULL)
        }, file = out_file),
        message = function(m) {
          if (sink.number(type = "message") == 2L) {
            err_chunks <<- c(err_chunks, conditionMessage(m))
            invokeRestart("muffleMessage")
          }
        },
        warning = function(w) {
          if (sink.number(type = "message") == 2L) {
            err_chunks <<- c(err_chunks, paste0("Warning message: ", conditionMessage(w), "\n"))
            invokeRestart("muffleWarning")
          }
        }
      ),
      error = function(e) {
        job_ok <<- FALSE
        err_chunks <<- c(err_chunks, paste0("Error: ", conditionMessage(e), "\n"))
      }
    )

    for (k in names(prev_env)) {
      pv <- prev_env[[k]]
      if (is.na(pv)) {
        Sys.unsetenv(k)
      } else {
        do.call(Sys.setenv, stats::setNames(list(pv), k))
      }
    }
    restore_global_state()

    writeLines(paste(err_chunks, collapse = ""), con = err_file, sep = "")
    writeLines(if (job_ok) "0" else "1", con = status_file)
    cat(done_token, "\n", sep = "")
    flush(stdout())
  }
})
"""

_STAGE_CONFIG: Dict[str, tuple[str, str]] = {
    "step03": ("scripts/morphometrics/03_fourier_extractor.R", "run_fourier_extraction"),
    "step04": ("scripts/morphometrics/04_gmm_morphotools.R", "run_gmm_morphotools_pipeline"),
    "step05": ("scripts/analysis/05_multimodal_spatial_rf.R", "main"),
    "step06": ("scripts/analysis/06_genomic_phylogenomics.R", "main"),
    "step07": ("scripts/analysis/07_triage_dashboard_synthesis.R", "main"),
}

_SCRIPT_TO_STAGE: Dict[str, str] = {
    "03_fourier_extractor.R": "step03",
    "04_gmm_morphotools.R": "step04",
    "05_multimodal_spatial_rf.R": "step05",
    "06_genomic_phylogenomics.R": "step06",
    "07_triage_dashboard_synthesis.R": "step07",
}

__all__ = ["_MANAGED_ENV_KEYS", "_WORKER_R_CODE", "_STAGE_CONFIG", "_SCRIPT_TO_STAGE"]
