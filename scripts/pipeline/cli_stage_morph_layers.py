"""CLI execution handlers for Phase 3 (`morphometrics`) and environmental/cytology/backbone provisioning."""

from __future__ import annotations

import argparse
from datetime import datetime
import os
from pathlib import Path
import shutil
import subprocess
import sys
from typing import Optional

from scripts.core.config import PipelineConfig, write_run_metadata
from scripts.pipeline.cli_run_context import (
    PROJECT_ROOT,
    _db_has_extracted_contours,
    _main_attr,
    _resolve_active_run_context,
    _resolve_focal_and_intraspecific,
    logger,
    verify_dir_has_files,
    verify_file_exists,
)


def run_morphometrics(args: argparse.Namespace, cfg: PipelineConfig) -> None:
    """Executes Phase 3: EFA, GMM, and passive sample CDA via R scripts."""
    _logger = _main_attr("logger", logger)
    _proj_root = _main_attr("PROJECT_ROOT", PROJECT_ROOT)
    _shutil = _main_attr("shutil", shutil)
    _subprocess = _main_attr("subprocess", subprocess)
    _verify_file_exists = _main_attr("verify_file_exists", verify_file_exists)
    _verify_dir_has_files = _main_attr("verify_dir_has_files", verify_dir_has_files)
    _resolve_ctx = _main_attr("_resolve_active_run_context", _resolve_active_run_context)
    _has_contours = _main_attr("_db_has_extracted_contours", _db_has_extracted_contours)
    _resolve_focal = _main_attr("_resolve_focal_and_intraspecific", _resolve_focal_and_intraspecific)
    _logger.info("=== Starting Phase 3: Morphometrics, EFA & Discriminant Analysis ===")
    start_ts = datetime.now().astimezone().isoformat()
    run_dir, cfg = _resolve_ctx(args, cfg)

    rscript_bin = _shutil.which("Rscript")
    if not rscript_bin:
        _logger.error(
            "Rscript binary not found in system PATH.\n"
            "Phase 3 requires R (>= 4.3) with Momocs, MorphoTools2, and mclust installed."
        )
        sys.exit(1)

    default_contours = (
        (run_dir / "data" / "contours")
        if run_dir is not None
        else cfg.get("paths", {}).get("contours_dir", "data/contours")
    )
    contours_dir = Path(getattr(args, "input", None) or args.contours_dir or default_contours)
    if not _has_contours(contours_dir):
        _verify_dir_has_files(contours_dir, "*.csv", "leaf contours", "python main.py segment")

    default_vouchers = (
        run_dir / "data" / "tables" / "curated_vouchers.csv"
        if run_dir is not None
        else cfg.get("paths", {}).get("curated_vouchers_csv", "data/tables/curated_vouchers.csv")
    )
    vouchers_csv = Path(args.vouchers or default_vouchers)
    _verify_file_exists(vouchers_csv, "curated vouchers table", "python main.py harvest")

    default_manifest = (
        run_dir / "data" / "tables" / "extracted_leaf_manifest.csv"
        if run_dir is not None
        else "data/tables/extracted_leaf_manifest.csv"
    )
    manifest_csv = Path(getattr(args, "manifest", None) or default_manifest)

    efa_script = _proj_root / "scripts" / "morphometrics" / "03_fourier_extractor.R"
    _verify_file_exists(efa_script, "EFA extractor R script")

    default_efa_out = (
        run_dir / "data" / "tables" / "leaf_efa_harmonics.csv"
        if run_dir is not None
        else cfg.get("paths", {}).get("leaf_efa_harmonics_csv", "data/tables/leaf_efa_harmonics.csv")
    )
    efa_out = Path(args.harmonics_out or default_efa_out)
    efa_out.parent.mkdir(parents=True, exist_ok=True)
    efa_ind_out = Path(
        getattr(args, "output_individual", None)
        or (efa_out.parent / "leaf_efa_harmonics_individual.csv")
    )
    default_masks = (
        (run_dir / "data" / "masks")
        if run_dir is not None
        else cfg.get("paths", {}).get("masks_dir", "data/masks")
    )
    masks_dir = Path(getattr(args, "masks_dir", None) or default_masks)
    efa_figures_dir = (
        run_dir / "outputs" / "figures"
        if run_dir is not None
        else Path(cfg.get("paths", {}).get("figures_dir", "outputs/figures"))
    )

    m_cfg = cfg.morphometrics
    harmonics = args.harmonics if getattr(args, "harmonics", None) is not None else m_cfg.harmonics
    num_pcs = args.num_pcs if getattr(args, "num_pcs", None) is not None else m_cfg.num_pcs
    sym_flag = "TRUE" if bool(getattr(m_cfg, "symmetrize_efa_harmonics", True)) else "FALSE"

    # Step 3A: Fourier EFA
    cmd_efa = [
        rscript_bin, str(efa_script),
        "--input", str(contours_dir),
        "--masks-dir", str(masks_dir),
        "--manifest", str(manifest_csv),
        "--vouchers", str(vouchers_csv),
        "--output", str(efa_out),
        "--output-individual", str(efa_ind_out),
        "--figures-dir", str(efa_figures_dir),
        "--harmonics", str(harmonics),
        "--num-pcs", str(num_pcs),
        "--symmetrize-efa-harmonics", sym_flag,
        "--seed", str(args.seed if getattr(args, "seed", None) is not None else getattr(m_cfg, "random_seed", 42)),
    ]
    _logger.info(f"Running EFA: {' '.join(cmd_efa)}")
    _subprocess.run(cmd_efa, check=True)
    _verify_file_exists(efa_out, "EFA harmonics table")

    # Step 3B: GMM & CDA
    gmm_script = _proj_root / "scripts" / "morphometrics" / "04_gmm_morphotools.R"
    _verify_file_exists(gmm_script, "GMM / MorphoTools2 R script")

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

    default_flags_out = (
        tables_dir / "morphometric_flags.csv"
        if run_dir is not None
        else cfg.get("paths", {}).get("morphometric_flags_csv", str(tables_dir / "morphometric_flags.csv"))
    )
    default_plot_out = figures_dir / "fig04b_gmm_cda_morphometrics.png"
    default_report_out = (
        run_dir / "outputs" / "reports" / "gmm_bayes_factors_summary.csv"
        if run_dir is not None
        else cfg.get("paths", {}).get("gmm_report_csv", "outputs/reports/gmm_bayes_factors_summary.csv")
    )
    flags_out = Path(args.flags_out or default_flags_out)
    plot_out = Path(args.plot_out or default_plot_out)
    report_out = Path(args.report_out or default_report_out)
    reproductive_csv = Path(getattr(args, "reproductive", None) or (tables_dir / "voucher_reproductive_metrics.csv"))
    allometry_report = Path(getattr(args, "allometry_report", None) or (tables_dir / "allometry_permanova_summary.csv"))
    allometry_pairwise = Path(getattr(args, "allometry_pairwise", None) or (tables_dir / "allometry_pairwise_slopes.csv"))
    allometry_plot = Path(getattr(args, "allometry_plot", None) or (figures_dir / "fig04a_foliar_allometry_audit.png"))
    ecoregion_plot = Path(getattr(args, "ecoregion_plot", None) or (figures_dir / "fig04c_ecoregion_and_batch_audit.png"))
    batch_audit_out = Path(getattr(args, "output_batch_audit", None) or (tables_dir / "batch_variance_audit.csv"))
    flags_out.parent.mkdir(parents=True, exist_ok=True)
    plot_out.parent.mkdir(parents=True, exist_ok=True)
    report_out.parent.mkdir(parents=True, exist_ok=True)

    m_cfg = cfg.get("morphometrics", {}) or {}
    max_k = args.max_k if getattr(args, "max_k", None) is not None else m_cfg.get("max_k", 8)
    correct_allometry = getattr(args, "correct_allometry", None)
    if correct_allometry is None:
        correct_allometry = m_cfg.get("correct_allometry", "auto")
    gmm_seed = getattr(args, "seed", None)
    if gmm_seed is None:
        gmm_seed = m_cfg.get("random_seed", 42)
    gmm_perms = getattr(args, "permutations", None)
    if gmm_perms is None:
        gmm_perms = m_cfg.get("permutations", None)

    cmd_gmm = [
        rscript_bin, str(gmm_script),
        "--vouchers", str(vouchers_csv),
        "--harmonics", str(efa_out),
        "--reproductive", str(reproductive_csv),
        "--correct-allometry", str(correct_allometry).lower() if isinstance(correct_allometry, bool) else str(correct_allometry),
        "--allometry-report", str(allometry_report),
        "--allometry-pairwise", str(allometry_pairwise),
        "--allometry-plot", str(allometry_plot),
        "--output-batch-audit", str(batch_audit_out),
        "--output-flags", str(flags_out),
        "--output-plot", str(plot_out),
        "--ecoregion-plot", str(ecoregion_plot),
        "--output-report", str(report_out),
        "--max-k", str(max_k),
        "--num-pcs", str(num_pcs),
        "--seed", str(gmm_seed),
    ]
    if gmm_perms is not None:
        cmd_gmm.extend(["--permutations", str(gmm_perms)])
    if getattr(args, "outline_only", False) or m_cfg.get("outline_only", False):
        cmd_gmm.append("--outline-only")
    if getattr(args, "use_clustvarsel", False) or m_cfg.get("use_clustvarsel", False):
        cmd_gmm.append("--use-clustvarsel")
    focal_taxon, intraspecific_mode = _resolve_focal(args, cfg)
    if focal_taxon:
        cmd_gmm.extend(["--focal-taxon", focal_taxon])
    if intraspecific_mode:
        cmd_gmm.append("--intraspecific-mode")
    _logger.info(f"Running GMM/CDA: {' '.join(cmd_gmm)}")
    _subprocess.run(cmd_gmm, check=True)
    _verify_file_exists(flags_out, "morphometric flags CSV")
    if run_dir is not None:
        write_run_metadata(run_dir, start_time=start_ts)
    _logger.info("Phase 3 completed successfully.")


def _maybe_ensure_environmental_layers(
    args: argparse.Namespace,
    env_dir: Path,
    vouchers_csv: Optional[Path] = None,
) -> None:
    """Ensure SOLUS100, USGS SGMC, 30 m DEM, DS-801, and NWI/NWPL layers exist before Step 05 unless in demo/test mode."""
    _logger = _main_attr("logger", logger)
    _proj_root = _main_attr("PROJECT_ROOT", PROJECT_ROOT)
    _subprocess = _main_attr("subprocess", subprocess)
    if (
        getattr(args, "demo", False)
        or getattr(args, "skip_env_download", False)
        or os.environ.get("PLANTMATICS_DEMO_MODE") == "1"
        or os.environ.get("PACKERA_DEMO_MODE") == "1"
        or os.environ.get("PLANTMATICS_SKIP_ENV_DOWNLOAD") == "1"
        or os.environ.get("PACKERA_SKIP_ENV_DOWNLOAD") == "1"
        or "PYTEST_CURRENT_TEST" in os.environ
        or getattr(_subprocess.run, "__module__", "") != "subprocess"
    ):
        return
    try:
        from scripts.data_prep.download_environmental_layers import ensure_environmental_layers
        ensure_environmental_layers(
            env_dir=env_dir,
            vouchers_csv=vouchers_csv,
            workspace_root=_proj_root,
            db_path=getattr(args, "db", None),
            force=bool(getattr(args, "force_env_download", False)),
        )
    except Exception as exc:
        _logger.warning(f"Environmental layer preflight check could not complete: {exc}")
    try:
        from scripts.data_prep.download_geochem_wetland_layers import download_geochem_wetland_layers
        download_geochem_wetland_layers(
            env_dir=env_dir,
            vouchers_csv=vouchers_csv,
            db_path=getattr(args, "db", None),
            genus=getattr(args, "genus", None),
            skip_nwi_network=bool(getattr(args, "skip_nwi_network", False) or os.environ.get("PLANTMATICS_SKIP_NWI_NETWORK") == "1"),
            force=bool(getattr(args, "force_env_download", False)),
        )
    except Exception as exc:
        _logger.warning(f"Geochemistry/wetland layer preflight check could not complete: {exc}")


def _maybe_ensure_cytology_records(args: argparse.Namespace) -> None:
    """Ensure flat data/cytology/ records exist and are synced to SQLite unless skipped or in demo/test mode."""
    _logger = _main_attr("logger", logger)
    _proj_root = _main_attr("PROJECT_ROOT", PROJECT_ROOT)
    _subprocess = _main_attr("subprocess", subprocess)
    if (
        getattr(args, "demo", False)
        or getattr(args, "skip_cytology_download", False)
        or os.environ.get("PLANTMATICS_DEMO_MODE") == "1"
        or os.environ.get("PACKERA_DEMO_MODE") == "1"
        or os.environ.get("PLANTMATICS_SKIP_CYTOLOGY_DOWNLOAD") == "1"
        or os.environ.get("PACKERA_SKIP_CYTOLOGY_DOWNLOAD") == "1"
        or "PYTEST_CURRENT_TEST" in os.environ
        or getattr(_subprocess.run, "__module__", "") != "subprocess"
    ):
        return
    dl_script = _proj_root / "scripts" / "data_prep" / "download_public_cytology.py"
    if not dl_script.is_file():
        return
    try:
        cyto_dir = Path(
            os.environ.get("PLANTMATICS_CYTOLOGY_DIR")
            or os.environ.get("PACKERA_CYTOLOGY_DIR")
            or (_proj_root / "data" / "cytology")
        )
        rec_csv = cyto_dir / "public_packera_cytology_records.csv"
        force_dl = bool(getattr(args, "force_cytology_download", False))
        cmd = [sys.executable, str(dl_script), "--cytology-dir", str(cyto_dir)]
        if getattr(args, "db", None):
            cmd.extend(["--db", str(args.db)])
        if force_dl:
            cmd.append("--force")
        elif rec_csv.is_file() and rec_csv.stat().st_size > 0:
            cmd.append("--sync-db-only")
        _subprocess.run(cmd, check=False)
    except Exception as exc:
        _logger.warning(f"Cytology preflight check could not complete: {exc}")


def _maybe_ensure_taxonomic_backbone(args: argparse.Namespace) -> None:
    """Ensure flat data/taxonomy/ taxonomic backbone exists and is synced unless skipped or in demo/test mode."""
    _logger = _main_attr("logger", logger)
    _proj_root = _main_attr("PROJECT_ROOT", PROJECT_ROOT)
    _subprocess = _main_attr("subprocess", subprocess)
    if (
        getattr(args, "demo", False)
        or getattr(args, "skip_backbone_download", False)
        or os.environ.get("PLANTMATICS_DEMO_MODE") == "1"
        or os.environ.get("PACKERA_DEMO_MODE") == "1"
        or os.environ.get("PLANTMATICS_SKIP_BACKBONE_DOWNLOAD") == "1"
        or os.environ.get("PACKERA_SKIP_BACKBONE_DOWNLOAD") == "1"
        or "PYTEST_CURRENT_TEST" in os.environ
        or getattr(_subprocess.run, "__module__", "") != "subprocess"
    ):
        return
    dl_script = _proj_root / "scripts" / "data_prep" / "download_taxonomic_backbone.py"
    if not dl_script.is_file():
        return
    try:
        tax_dir = Path(
            os.environ.get("PLANTMATICS_TAXONOMY_DIR")
            or os.environ.get("PACKERA_TAXONOMY_DIR")
            or (_proj_root / "data" / "taxonomy")
        )
        script_txt = dl_script.read_text(encoding="utf-8", errors="ignore")
        force_dl = bool(getattr(args, "force_backbone_download", False))
        cmd = [sys.executable, str(dl_script)]
        if getattr(args, "db", None) and "--db" in script_txt:
            cmd.extend(["--db", str(args.db)])
        if force_dl and "--force" in script_txt:
            cmd.append("--force")
        elif any(tax_dir.glob("*.csv")) or any(tax_dir.glob("*.json")):
            if "--sync-db-only" in script_txt:
                cmd.append("--sync-db-only")
            elif "--sync-only" in script_txt:
                cmd.append("--sync-only")
            else:
                return
        _subprocess.run(cmd, check=False)
    except Exception as exc:
        _logger.warning(f"Taxonomic backbone preflight check could not complete: {exc}")


def run_download_env(args: argparse.Namespace, cfg: Optional[PipelineConfig] = None) -> None:
    """Verify or download SOLUS100, USGS SGMC, 3DEP 30 m DEM, DS-801, and NWI/NWPL layers."""
    _logger = _main_attr("logger", logger)
    _proj_root = _main_attr("PROJECT_ROOT", PROJECT_ROOT)
    _ = cfg
    from scripts.data_prep.download_environmental_layers import ensure_environmental_layers
    env_dir = Path(getattr(args, "env_dir", None) or (_proj_root / "data" / "environmental"))
    vouchers_csv = getattr(args, "vouchers", None)
    ensure_environmental_layers(
        env_dir=env_dir,
        vouchers_csv=vouchers_csv,
        workspace_root=_proj_root,
        db_path=getattr(args, "db", None),
        force=bool(getattr(args, "force", False)),
        demo=bool(getattr(args, "demo", False)),
    )
    if not getattr(args, "demo", False) and os.environ.get("PLANTMATICS_DEMO_MODE") != "1" and os.environ.get("PACKERA_DEMO_MODE") != "1":
        try:
            from scripts.data_prep.download_geochem_wetland_layers import download_geochem_wetland_layers
            download_geochem_wetland_layers(
                env_dir=env_dir,
                vouchers_csv=vouchers_csv,
                db_path=getattr(args, "db", None),
                genus=getattr(args, "genus", None),
                skip_nwi_network=bool(getattr(args, "skip_nwi_network", False)),
                force=bool(getattr(args, "force", False)),
            )
        except Exception as exc:
            _logger.warning(f"Geochemistry/wetland layer download could not complete: {exc}")
