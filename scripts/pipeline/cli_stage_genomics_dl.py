"""CLI execution handlers for Step 06 genomics, organellar skimming, and public genomics/backbone/cytology downloaders."""

from __future__ import annotations

import argparse
import os
from pathlib import Path
import shutil
import subprocess
import sys
from typing import Optional

from scripts.core.config import PipelineConfig, resolve_canonical_run_id
from scripts.pipeline.cli_run_context import (
    PROJECT_ROOT,
    _main_attr,
    _resolve_active_run_context,
    logger,
    verify_file_exists,
)


def run_genomics_stream(
    args: argparse.Namespace,
    run_dir: Optional[Path],
    tables_dir: Path,
    figures_dir: Path,
    rscript_bin: str,
    reports_dir: Optional[Path] = None,
    *,
    respect_test_mock_guard: bool = True,
) -> Optional[Path]:
    """Execute standalone Step 06 genomic ingestion and phylogenomic integration when data/genomics/ has inputs."""
    _logger = _main_attr("logger", logger)
    _proj_root = _main_attr("PROJECT_ROOT", PROJECT_ROOT)
    _subprocess = _main_attr("subprocess", subprocess)
    _verify_file_exists = _main_attr("verify_file_exists", verify_file_exists)
    env_geno = (os.environ.get("PLANTMATICS_GENOMICS_DIR") or os.environ.get("PACKERA_GENOMICS_DIR") or "").strip()
    if respect_test_mock_guard and not env_geno and getattr(_subprocess.run, "__module__", "") != "subprocess":
        return None

    import importlib

    ingest_mod = importlib.import_module("scripts.data_prep.06_ingest_genetic_data")
    genomics_dir = Path(getattr(args, "genomics_dir", None) or env_geno or (_proj_root / "data" / "genomics"))
    assets = ingest_mod.discover_genomic_inputs(genomics_dir)
    has_pub = ingest_mod.has_valid_public_genomic_records(assets)
    focal_arg = getattr(args, "focal_taxon", None)
    intra_flag = bool(getattr(args, "intraspecific_mode", False))
    vouchers_csv = tables_dir / "curated_vouchers.csv"
    if not (assets.has_genetic_payload or (has_pub and (run_dir is not None or bool(focal_arg)))):
        _logger.info("[Step 06/07] No genomic files staged in data/genomics/ - skipping standalone Step 06")
        print("[Step 06/07] No genomic files staged in data/genomics/ - skipping standalone Step 06")
        return None

    eff_reports_dir = (
        Path(reports_dir)
        if reports_dir is not None
        else ((run_dir / "outputs" / "reports") if run_dir is not None else (_proj_root / "outputs" / "reports"))
    )
    eff_reports_dir.mkdir(parents=True, exist_ok=True)
    geno_out_dir = (
        (run_dir / "outputs" / "genomics")
        if run_dir is not None
        else (_proj_root / "outputs" / "genomics")
    )
    geno_out_dir.mkdir(parents=True, exist_ok=True)
    db_arg = (
        getattr(args, "db", None)
        or os.environ.get("PLANTMATICS_DB_PATH")
        or os.environ.get("PACKERA_SQLITE_DB")
        or os.environ.get("PACKERA_DB_PATH")
        or None
    )
    is_demo = bool(getattr(args, "demo", False))
    eff_run_id = (
        resolve_canonical_run_id(run_dir=run_dir, workspace_root=_proj_root)
        if run_dir is not None
        else (getattr(args, "run_id", None) or "canonical")
    )

    ingest_script = _proj_root / "scripts" / "data_prep" / "06_ingest_genetic_data.py"
    _verify_file_exists(ingest_script, "Step 06 Genomic Ingestion Python script")
    cmd_ingest = [
        sys.executable,
        str(ingest_script),
        "--genomics-dir",
        str(genomics_dir),
        "--out-dir",
        str(eff_reports_dir),
    ]
    if vouchers_csv.is_file():
        cmd_ingest.extend(["--vouchers", str(vouchers_csv)])
    if vouchers_csv.is_file() or focal_arg:
        cmd_ingest.extend(["--genomic-tier", "auto"])
    if focal_arg:
        cmd_ingest.extend(["--focal-taxon", str(focal_arg)])
    if intra_flag:
        cmd_ingest.append("--intraspecific-mode")
    if db_arg and str(db_arg) != ":memory:":
        cmd_ingest.extend(["--db-path", str(db_arg)])
    if is_demo:
        cmd_ingest.append("--demo")

    _logger.info("Detected staged genomic inputs in %s; running standalone Step 06.", genomics_dir)
    _subprocess.run(cmd_ingest, check=True)

    geno_script = _proj_root / "scripts" / "analysis" / "06_genomic_phylogenomics.R"
    _verify_file_exists(geno_script, "Step 06 Genomic Phylogenomics R script")
    cmd_geno = [
        rscript_bin,
        str(geno_script),
        "--genomics-dir",
        str(genomics_dir),
        "--parsed-dir",
        str(eff_reports_dir),
        "--out-dir",
        str(geno_out_dir),
        "--tables-dir",
        str(tables_dir),
        "--figure-a",
        str(figures_dir / "fig06a_genomic_pcoa_tree_concordance.png"),
        "--figure-b",
        str(figures_dir / "fig06b_reticulation_gdi_landscape.png"),
        "--figure-c",
        str(figures_dir / "fig06c_cytonuclear_paralog_cytotype_audit.png"),
        "--run-id",
        eff_run_id,
    ]
    if vouchers_csv.is_file():
        cmd_geno.extend(["--vouchers", str(vouchers_csv)])
    if vouchers_csv.is_file() or focal_arg:
        cmd_geno.extend(["--genomic-tier", "auto"])
    if focal_arg:
        cmd_geno.extend(["--focal-taxon", str(focal_arg)])
    if intra_flag:
        cmd_geno.append("--intraspecific-mode")
    if db_arg and str(db_arg) != ":memory:":
        cmd_geno.extend(["--db", str(db_arg)])
    if is_demo:
        cmd_geno.append("--demo")

    _logger.info(f"Running Step 06 Genomic Phylogenomics: {' '.join(cmd_geno)}")
    _subprocess.run(cmd_geno, check=True)
    return geno_out_dir


def run_genomics_cli(args: argparse.Namespace, cfg: PipelineConfig) -> None:
    """CLI entrypoint for standalone Step 06 genomic ingestion & phylogenomics."""
    _logger = _main_attr("logger", logger)
    _shutil = _main_attr("shutil", shutil)
    _resolve_ctx = _main_attr("_resolve_active_run_context", _resolve_active_run_context)
    _run_geno_stream = _main_attr("run_genomics_stream", run_genomics_stream)
    run_dir, cfg = _resolve_ctx(args, cfg)
    rscript_bin = _shutil.which("Rscript")
    if not rscript_bin:
        _logger.error("Rscript binary not found in system PATH.")
        sys.exit(1)
    tables_dir = (run_dir / "data" / "tables") if run_dir is not None else Path(cfg.get("paths", {}).get("tables_dir", "data/tables"))
    figures_dir = (run_dir / "outputs" / "figures") if run_dir is not None else Path(cfg.get("paths", {}).get("figures_dir", "outputs/figures"))
    reports_dir = (run_dir / "outputs" / "reports") if run_dir is not None else Path(cfg.get("paths", {}).get("reports_dir", str(figures_dir.parent / "reports")))
    tables_dir.mkdir(parents=True, exist_ok=True)
    figures_dir.mkdir(parents=True, exist_ok=True)
    reports_dir.mkdir(parents=True, exist_ok=True)
    _run_geno_stream(
        args,
        run_dir,
        tables_dir,
        figures_dir,
        rscript_bin,
        reports_dir=reports_dir,
        respect_test_mock_guard=False,
    )


def run_skim_organellar_cli(args: argparse.Namespace, cfg: PipelineConfig) -> None:
    """CLI entrypoint for offline off-target plastid/ITS/waxy skimming into data/genomics/public_organellar_and_loci/."""
    _proj_root = _main_attr("PROJECT_ROOT", PROJECT_ROOT)
    _subprocess = _main_attr("subprocess", subprocess)
    _verify_file_exists = _main_attr("verify_file_exists", verify_file_exists)
    _ = cfg
    skim_script = _proj_root / "scripts" / "data_prep" / "skim_offtarget_organellar.py"
    _verify_file_exists(skim_script, "Off-target Organellar Skimmer Python script")
    cmd = [sys.executable, str(skim_script)]
    if getattr(args, "fastq_source", None):
        cmd.extend(["--fastq-source", str(args.fastq_source)])
    if getattr(args, "genomics_dir", None):
        cmd.extend(["--genomics-dir", str(args.genomics_dir)])
    if getattr(args, "nuclear_tree", None):
        cmd.extend(["--nuclear-tree", str(args.nuclear_tree)])
    if getattr(args, "out_dir", None):
        cmd.extend(["--out-dir", str(args.out_dir)])
    if getattr(args, "threads", None) is not None:
        cmd.extend(["--threads", str(args.threads)])
    _subprocess.run(cmd, check=True)


def run_download_public_genomics_cli(args: argparse.Namespace, cfg: PipelineConfig) -> None:
    """CLI entrypoint for downloading public Packera sequences and trees into data/genomics/ subdirectories."""
    _proj_root = _main_attr("PROJECT_ROOT", PROJECT_ROOT)
    _subprocess = _main_attr("subprocess", subprocess)
    _verify_file_exists = _main_attr("verify_file_exists", verify_file_exists)
    _ = cfg
    dl_script = _proj_root / "scripts" / "data_prep" / "download_public_packera_sequences.py"
    _verify_file_exists(dl_script, "Public Packera Sequence Downloader Python script")
    cmd = [sys.executable, str(dl_script)]
    if getattr(args, "genomics_dir", None):
        cmd.extend(["--genomics-dir", str(args.genomics_dir)])
    if getattr(args, "sources", None):
        cmd.extend(["--sources", *[str(s) for s in args.sources]])
    if getattr(args, "dryad_zip_dir", None):
        cmd.extend(["--dryad-zip-dir", str(args.dryad_zip_dir)])
    if getattr(args, "genus", None):
        cmd.extend(["--genus", str(args.genus)])
    if getattr(args, "family", None):
        cmd.extend(["--family", str(args.family)])
    if getattr(args, "skim_sra_accessions", None):
        cmd.extend(["--skim-sra-accessions", *[str(a) for a in args.skim_sra_accessions]])
    if getattr(args, "force", False):
        cmd.append("--force")
    _subprocess.run(cmd, check=True)


def run_download_taxonomic_backbone_cli(args: argparse.Namespace, cfg: PipelineConfig) -> None:
    """CLI entrypoint for downloading dynamic taxonomic backbone profiles into data/taxonomy/."""
    _proj_root = _main_attr("PROJECT_ROOT", PROJECT_ROOT)
    _subprocess = _main_attr("subprocess", subprocess)
    _verify_file_exists = _main_attr("verify_file_exists", verify_file_exists)
    _ = cfg
    dl_script = _proj_root / "scripts" / "data_prep" / "download_taxonomic_backbone.py"
    _verify_file_exists(dl_script, "Taxonomic Backbone Downloader Python script")
    script_txt = dl_script.read_text(encoding="utf-8", errors="ignore")
    cmd = [sys.executable, str(dl_script)]
    if getattr(args, "taxonomy_dir", None):
        if "--taxonomy-dir" in script_txt:
            cmd.extend(["--taxonomy-dir", str(args.taxonomy_dir)])
        elif "--out-dir" in script_txt:
            cmd.extend(["--out-dir", str(args.taxonomy_dir)])
    if getattr(args, "genus", None) and "--genus" in script_txt:
        cmd.extend(["--genus", str(args.genus)])
    if getattr(args, "family", None) and "--family" in script_txt:
        cmd.extend(["--family", str(args.family)])
    if getattr(args, "sources", None) and "--sources" in script_txt:
        cmd.extend(["--sources", *[str(s) for s in args.sources]])
    if getattr(args, "sync_only", False):
        if "--sync-only" in script_txt:
            cmd.append("--sync-only")
        elif "--sync-db-only" in script_txt:
            cmd.append("--sync-db-only")
    if getattr(args, "db", None) and "--db" in script_txt:
        cmd.extend(["--db", str(args.db)])
    if getattr(args, "no_db", False) and "--no-db" in script_txt:
        cmd.append("--no-db")
    if getattr(args, "force", False) and "--force" in script_txt:
        cmd.append("--force")
    _subprocess.run(cmd, check=True)


def run_download_public_cytology_cli(args: argparse.Namespace, cfg: PipelineConfig) -> None:
    """CLI entrypoint for downloading and cross-linking public Packera cytological records into data/cytology/."""
    _proj_root = _main_attr("PROJECT_ROOT", PROJECT_ROOT)
    _subprocess = _main_attr("subprocess", subprocess)
    _verify_file_exists = _main_attr("verify_file_exists", verify_file_exists)
    _ = cfg
    dl_script = _proj_root / "scripts" / "data_prep" / "download_public_cytology.py"
    _verify_file_exists(dl_script, "Public Packera Cytology Downloader Python script")
    cmd = [sys.executable, str(dl_script)]
    if getattr(args, "cytology_dir", None):
        cmd.extend(["--cytology-dir", str(args.cytology_dir)])
    if getattr(args, "sources", None):
        cmd.extend(["--sources", *[str(s) for s in args.sources]])
    if getattr(args, "local_cytology_csv", None):
        cmd.extend(["--local-cytology-csv", str(args.local_cytology_csv)])
    if getattr(args, "sync_db_only", False):
        cmd.append("--sync-db-only")
    if getattr(args, "db", None):
        cmd.extend(["--db", str(args.db)])
    if getattr(args, "no_db", False):
        cmd.append("--no-db")
    if getattr(args, "force", False):
        cmd.append("--force")
    _subprocess.run(cmd, check=True)
