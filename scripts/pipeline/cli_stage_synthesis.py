"""CLI execution handlers for synthesis (Steps 05 to 07), `run-all`, and Tier A/B purge operations."""

from __future__ import annotations

import argparse
from datetime import datetime
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
from typing import Optional

from scripts.core.config import (
    PipelineConfig,
    generate_default_run_id,
    resolve_canonical_run_id,
    write_run_metadata,
)
from scripts.core.db import purge_runs_only, purge_segments_and_runs
from scripts.pipeline.cli_run_context import (
    PROJECT_ROOT,
    _apply_explicit_db_arg,
    _main_attr,
    _resolve_active_run_context,
    _resolve_focal_and_intraspecific,
    _sync_run_cohort_in_db,
    logger,
    verify_file_exists,
)
from scripts.pipeline.cli_stage_genomics_dl import run_genomics_stream
from scripts.pipeline.cli_stage_harvest_seg import run_harvest, run_segment
from scripts.pipeline.cli_stage_morph_layers import (
    _maybe_ensure_cytology_records,
    _maybe_ensure_environmental_layers,
    _maybe_ensure_taxonomic_backbone,
    run_morphometrics,
)


def run_synthesis(args: argparse.Namespace, cfg: PipelineConfig) -> None:
    """Executes Phases 5 and 6: Spatial Macroecology and Synthesis Triage."""
    _logger = _main_attr("logger", logger)
    _proj_root = _main_attr("PROJECT_ROOT", PROJECT_ROOT)
    _shutil = _main_attr("shutil", shutil)
    _subprocess = _main_attr("subprocess", subprocess)
    _verify_file_exists = _main_attr("verify_file_exists", verify_file_exists)
    _resolve_ctx = _main_attr("_resolve_active_run_context", _resolve_active_run_context)
    _resolve_focal = _main_attr("_resolve_focal_and_intraspecific", _resolve_focal_and_intraspecific)
    _ensure_env = _main_attr("_maybe_ensure_environmental_layers", _maybe_ensure_environmental_layers)
    _ensure_cyto = _main_attr("_maybe_ensure_cytology_records", _maybe_ensure_cytology_records)
    _ensure_bb = _main_attr("_maybe_ensure_taxonomic_backbone", _maybe_ensure_taxonomic_backbone)
    _run_geno_stream = _main_attr("run_genomics_stream", run_genomics_stream)
    _logger.info("=== Starting Synthesis Workflow (Phases 5, 6) ===")
    start_ts = datetime.now().astimezone().isoformat()
    run_dir, cfg = _resolve_ctx(args, cfg)

    rscript_bin = _shutil.which("Rscript")
    if not rscript_bin:
        _logger.error(
            "Rscript binary not found in system PATH.\n"
            "Spatial macroecology and synthesis dashboard require R (>= 4.3) with spatialRF, ENMTools, and terra."
        )
        sys.exit(1)

    tables_dir = (
        run_dir / "data" / "tables"
        if run_dir is not None
        else Path(cfg.get("paths", {}).get("tables_dir", "data/tables"))
    )
    figures_dir = (
        run_dir / "outputs" / "figures"
        if run_dir is not None
        else Path(cfg.get("paths", {}).get("figures_dir", "outputs/figures"))
    )
    tables_dir.mkdir(parents=True, exist_ok=True)
    figures_dir.mkdir(parents=True, exist_ok=True)
    focal_taxon, intraspecific_mode = _resolve_focal(args, cfg)

    default_vouchers = (
        tables_dir / "curated_vouchers.csv"
        if run_dir is not None
        else cfg.get("paths", {}).get("curated_vouchers_csv", str(tables_dir / "curated_vouchers.csv"))
    )
    vouchers_csv = Path(getattr(args, "vouchers", None) or default_vouchers)
    _verify_file_exists(vouchers_csv, "curated vouchers table", "python main.py harvest")

    morph_flags = Path(
        getattr(args, "morphometrics_csv", None)
        or getattr(args, "morphometrics", None)
        or (tables_dir / "morphometric_flags.csv")
    )

    # -------------------------------------------------------------------------
    # Stage 5: Multimodal Spatial Random Forests & Warren's Niche Identity Tests
    # -------------------------------------------------------------------------
    spatial_rf_script = _proj_root / "scripts" / "analysis" / "05_multimodal_spatial_rf.R"
    _verify_file_exists(spatial_rf_script, "Spatial RF R script")

    env_dir = Path(getattr(args, "env_dir", None) or "data/environmental")
    _ensure_env(args, env_dir=env_dir, vouchers_csv=vouchers_csv)
    _ensure_cyto(args)
    _ensure_bb(args)
    reports_dir_s05 = (
        run_dir / "outputs" / "reports"
        if run_dir is not None
        else Path(cfg.get("paths", {}).get("reports_dir", str(figures_dir.parent / "reports")))
    )
    reports_dir_s05.mkdir(parents=True, exist_ok=True)
    pheno_states = Path(getattr(args, "pheno_states", None) or (tables_dir / "voucher_phenological_states.csv"))
    ssurgo_table = Path(getattr(args, "ssurgo_table", None) or (tables_dir / "ssurgo_edaphic_validation.csv"))
    contingency_table = Path(getattr(args, "contingency_table", None) or (reports_dir_s05 / "micro_edaphic_outcrop_contingency.csv"))
    pheno_anomalies = Path(getattr(args, "pheno_anomalies", None) or (tables_dir / "phenological_anomalies.csv"))
    pheno_summary = Path(getattr(args, "pheno_summary", None) or (reports_dir_s05 / "phenological_anomaly_summary.csv"))
    pheno_plot = Path(getattr(args, "pheno_plot", None) or (figures_dir / "fig05b_hopkins_circular_phenology.png"))
    conflict_flags = Path(
        getattr(args, "conflict_flags", None)
        or (tables_dir / "multimodal_conflict_flags.csv")
    )
    spatial_plot = Path(
        getattr(args, "spatial_plot", None)
        or (figures_dir / "fig05a_ecological_niche_and_rf_importance.png")
    )
    sdm_plot = Path(
        getattr(args, "sdm_plot", None)
        or (figures_dir / "fig05c_sdm_suitability_and_distribution_map.png")
    )
    spatial_summary = Path(
        getattr(args, "spatial_summary", None)
        or (reports_dir_s05 / "multimodal_spatial_rf_summary.csv")
    )

    cmd_spatial = [
        rscript_bin, str(spatial_rf_script),
        "--vouchers", str(vouchers_csv),
        "--morphometrics", str(morph_flags),
        "--pheno-states", str(pheno_states),
        "--env-dir", str(env_dir),
        "--ssurgo-table", str(ssurgo_table),
        "--contingency-table", str(contingency_table),
        "--pheno-anomalies", str(pheno_anomalies),
        "--pheno-summary", str(pheno_summary),
        "--pheno-plot", str(pheno_plot),
        "--output-flags", str(conflict_flags),
        "--output-plot", str(spatial_plot),
        "--sdm-plot", str(sdm_plot),
        "--output-summary", str(spatial_summary),
    ]
    if getattr(args, "permutations", None) is not None:
        cmd_spatial.extend(["--permutations", str(args.permutations)])
    if getattr(args, "n_trees", None) is not None:
        cmd_spatial.extend(["--n-trees", str(args.n_trees)])
    if getattr(args, "ssurgo_batch_size", None) is not None:
        cmd_spatial.extend(["--ssurgo-batch-size", str(args.ssurgo_batch_size)])
    if focal_taxon:
        cmd_spatial.extend(["--focal-taxon", focal_taxon])
    if intraspecific_mode:
        cmd_spatial.append("--intraspecific-mode")
    if getattr(args, "demo", False):
        cmd_spatial.append("--demo")

    _logger.info(f"Running Phase 5 Multimodal Spatial RF: {' '.join(cmd_spatial)}")
    _subprocess.run(cmd_spatial, check=True)
    _verify_file_exists(conflict_flags, "multimodal conflict flags CSV")

    # -------------------------------------------------------------------------
    # Stage 6 (Automatic when data/genomics/ populated): Standalone Step 06 Genomic & Phylogenomic Stream
    # -------------------------------------------------------------------------
    geno_out_dir = _run_geno_stream(
        args,
        run_dir,
        tables_dir,
        figures_dir,
        rscript_bin,
        reports_dir=reports_dir_s05,
    )

    # -------------------------------------------------------------------------
    # Stage 7: Multi-Evidence Taxonomic Decision Matrix & Triage Dashboard
    # -------------------------------------------------------------------------
    triage_script = _proj_root / "scripts" / "analysis" / "07_triage_dashboard_synthesis.R"
    _verify_file_exists(triage_script, "Step 07 Triage Dashboard R script")

    default_gmm_summary = (
        run_dir / "outputs" / "reports" / "gmm_bayes_factors_summary.csv"
        if run_dir is not None
        else cfg.get("paths", {}).get("gmm_report_csv", "outputs/reports/gmm_bayes_factors_summary.csv")
    )
    gmm_summary = Path(getattr(args, "gmm_summary", None) or default_gmm_summary)
    output_queue = Path(
        getattr(args, "output_queue", None)
        or ((run_dir / "data" / "tables" / "triage_queue.csv") if run_dir is not None else str(tables_dir / "triage_queue.csv"))
    )
    focal_taxon_raw = focal_taxon or ""
    if focal_taxon_raw:
        focal_slug = re.sub(r"[^A-Za-z0-9]+", "_", focal_taxon_raw).strip("_") or "Packera"
        default_fig_name = f"Figure_Integrative_{focal_slug}_Revision.png"
        default_rep_name = f"{focal_slug}_Taxonomic_Revision_Summary.md"
        fallback_fig = str(figures_dir / default_fig_name)
        fallback_rep = f"outputs/reports/{default_rep_name}"
    else:
        default_fig_name = "fig07a_integrative_taxonomic_synthesis.png"
        default_rep_name = "Packera_Taxonomic_Revision_Summary.md"
        fallback_fig = cfg.get("paths", {}).get("cda_biplot_png", str(figures_dir / default_fig_name))
        fallback_rep = cfg.get("paths", {}).get("revision_report_md", f"outputs/reports/{default_rep_name}")

    synthesis_plot = Path(
        getattr(args, "synthesis_plot", None)
        or ((figures_dir / default_fig_name) if run_dir is not None else fallback_fig)
    )
    synthesis_report = Path(
        getattr(args, "synthesis_report", None)
        or ((run_dir / "outputs" / "reports" / default_rep_name) if run_dir is not None else fallback_rep)
    )
    synthesis_html = Path(
        getattr(args, "synthesis_html", None)
        or (reports_dir_s05 / "Triage_Dashboard_Interactive.html")
    )
    harmonics_csv = tables_dir / "leaf_efa_harmonics.csv"
    output_complexes = Path(getattr(args, "output_complexes", None) or (reports_dir_s05 / "empirical_taxon_complexes.csv"))
    output_affinity = Path(getattr(args, "output_affinity", None) or (reports_dir_s05 / "pairwise_taxon_affinity_matrix.csv"))
    output_neighbors = Path(getattr(args, "output_neighbors", None) or (reports_dir_s05 / "focal_auto_neighbors.txt"))
    output_complex_plot = Path(getattr(args, "output_complex_plot", None) or (figures_dir / "fig07b_taxon_complex_delimitation_network.png"))
    min_affinity_val = getattr(args, "min_affinity", None)
    if min_affinity_val is None:
        min_affinity_val = float((cfg.get("synthesis", {}) or {}).get("min_affinity", 0.03))

    cmd_triage = [
        rscript_bin, str(triage_script),
        "--vouchers", str(vouchers_csv),
        "--harmonics", str(harmonics_csv),
        "--morphometrics", str(morph_flags),
        "--multimodal-flags", str(conflict_flags),
        "--gmm-summary", str(gmm_summary),
        "--niche-summary", str(spatial_summary),
        "--output-queue", str(output_queue),
        "--output-plot", str(synthesis_plot),
        "--output-report", str(synthesis_report),
        "--output-html", str(synthesis_html),
        "--min-affinity", str(min_affinity_val),
        "--output-complexes", str(output_complexes),
        "--output-affinity", str(output_affinity),
        "--output-neighbors", str(output_neighbors),
        "--output-complex-plot", str(output_complex_plot),
    ]
    if focal_taxon:
        cmd_triage.extend(["--focal-taxon", focal_taxon])
    if geno_out_dir is not None:
        cmd_triage.extend(["--genomic-dir", str(geno_out_dir)])
    if intraspecific_mode:
        cmd_triage.append("--intraspecific-mode")
    if getattr(args, "demo", False):
        cmd_triage.append("--demo")
    _logger.info(f"Running Step 07 Triage Dashboard Synthesis: {' '.join(cmd_triage)}")
    _subprocess.run(cmd_triage, check=True)
    _verify_file_exists(output_queue, "triage queue CSV")
    _verify_file_exists(synthesis_report, "taxonomic revision summary markdown report")
    _verify_file_exists(synthesis_html, "interactive triage HTML dashboard")

    if run_dir is not None:
        write_run_metadata(run_dir, start_time=start_ts)
    _logger.info("=== Synthesis Workflow (Phases 5, 6) Completed Successfully ===")


def run_all(args: argparse.Namespace, cfg: PipelineConfig) -> None:
    """Executes the end-to-end production workflow sequentially inside a run-scoped directory."""
    _logger = _main_attr("logger", logger)
    _proj_root = _main_attr("PROJECT_ROOT", PROJECT_ROOT)
    _resolve_ctx = _main_attr("_resolve_active_run_context", _resolve_active_run_context)
    _sync_cohort = _main_attr("_sync_run_cohort_in_db", _sync_run_cohort_in_db)
    _run_harvest = _main_attr("run_harvest", run_harvest)
    _run_segment = _main_attr("run_segment", run_segment)
    _run_morph = _main_attr("run_morphometrics", run_morphometrics)
    _run_synth = _main_attr("run_synthesis", run_synthesis)
    start_ts = datetime.now().astimezone().isoformat()
    run_id = getattr(args, "run_id", None) or os.environ.get("RUN_ID") or generate_default_run_id()
    setattr(args, "run_id", run_id)
    run_dir, cfg = _resolve_ctx(args, cfg)
    eff_run_id = (
        resolve_canonical_run_id(run_dir=run_dir, workspace_root=_proj_root)
        if run_dir is not None
        else run_id
    )
    if run_dir is not None:
        write_run_metadata(run_dir, run_id=eff_run_id, start_time=start_ts)

    _logger.info("==================================================================")
    _logger.info(f"Executing End-to-End Production Pipeline (Run: {eff_run_id})")
    _logger.info("  (harvest -> segment -> morphometrics -> synthesis)")
    _logger.info("==================================================================")
    _run_harvest(args, cfg)
    _run_segment(args, cfg)
    _run_morph(args, cfg)
    _run_synth(args, cfg)
    if run_dir is not None:
        end_ts = datetime.now().astimezone().isoformat()
        _sync_cohort(run_dir, args=args, cfg=cfg, end_timestamp=end_ts)
        write_run_metadata(run_dir, run_id=eff_run_id, start_time=start_ts, end_time=end_ts)
    _logger.info("==================================================================")
    _logger.info("End-to-End Pipeline Execution Completed Successfully!")
    _logger.info("==================================================================")


def run_purge_segments(args: argparse.Namespace, cfg: Optional[PipelineConfig] = None) -> None:
    """Purges all Step 02 to 07 segments, masks, contours, run directories, and SQLite tables."""
    _logger = _main_attr("logger", logger)
    _proj_root = _main_attr("PROJECT_ROOT", PROJECT_ROOT)
    _apply_db = _main_attr("_apply_explicit_db_arg", _apply_explicit_db_arg)
    _purge_seg = _main_attr("purge_segments_and_runs", purge_segments_and_runs)
    _ = cfg
    _apply_db(args)
    summary = _purge_seg(
        workspace_root=_proj_root,
        db_path=getattr(args, "db", None),
        output_dir=getattr(args, "output_dir", None),
        purge_runs=True,
        vacuum=True,
    )
    _logger.info(
        "Purge complete: removed %d artifact files and %d SQLite rows (Step 01 vouchers & raw images preserved).",
        summary["removed_files"],
        summary["purged_db_rows"],
    )


def run_purge_runs_only(args: argparse.Namespace, cfg: Optional[PipelineConfig] = None) -> None:
    """Purges only Step 03 to 07 run outputs and Tier B SQLite tables while preserving Step 02 segmentation."""
    _logger = _main_attr("logger", logger)
    _proj_root = _main_attr("PROJECT_ROOT", PROJECT_ROOT)
    _apply_db = _main_attr("_apply_explicit_db_arg", _apply_explicit_db_arg)
    _purge_runs = _main_attr("purge_runs_only", purge_runs_only)
    _ = cfg
    _apply_db(args)
    summary = _purge_runs(
        workspace_root=_proj_root,
        db_path=getattr(args, "db", None),
        output_dir=getattr(args, "output_dir", None),
        refresh_env_buffers=True,
        vacuum=True,
    )
    _logger.info(
        "Run-output purge complete: removed %d run/output files and %d Tier B SQLite rows (Step 01 & Step 02 segmentation preserved).",
        summary["removed_files"],
        summary["purged_db_rows"],
    )
