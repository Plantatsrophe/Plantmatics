"""Environment diagnostics (`check-env`) and model weights acquisition (`download-weights`) CLI handlers."""

from __future__ import annotations

import argparse
import importlib
from pathlib import Path
import shutil
import subprocess
import sys
from typing import Any, Dict, List, Optional

import yaml

from scripts.core.artifact_manager import ensure_model_weights
from scripts.core.config import PipelineConfig
from scripts.pipeline.cli_run_context import PROJECT_ROOT, _main_attr, logger


def check_environment(args: Optional[argparse.Namespace] = None) -> bool:
    """Performs comprehensive diagnostic checks across Python, CUDA, R,

    configuration files, filesystem directories, and model checkpoints.

    Returns:
        bool: True if all critical checks pass, False if any critical check fails.
    """
    _proj_root = _main_attr("PROJECT_ROOT", PROJECT_ROOT)
    _importlib = _main_attr("importlib", importlib)
    _shutil = _main_attr("shutil", shutil)
    _subprocess = _main_attr("subprocess", subprocess)
    _PipelineConfig = _main_attr("PipelineConfig", PipelineConfig)
    _ensure_weights = _main_attr("ensure_model_weights", ensure_model_weights)
    strict = getattr(args, "strict", False) if args else False
    custom_config = getattr(args, "config", None) if args else None
    critical_failure = False
    warning_count = 0

    print("=" * 79)
    print("Packera Morphometrics Pipeline: Environment & Dependency Diagnostic")
    print("=" * 79)

    # -------------------------------------------------------------------------
    # 1. Python & CUDA Environment
    # -------------------------------------------------------------------------
    print("\n[1/4] Python & CUDA Environment:")

    # Python version check (>= 3.10)
    py_ver_str = f"{sys.version_info.major}.{sys.version_info.minor}.{sys.version_info.micro}"
    if sys.version_info >= (3, 10):
        print(f"  [PASS] Python version: {py_ver_str} (>= 3.10 required)")
    else:
        print(f"  [FAIL] Python version: {py_ver_str} is unsupported (>= 3.10 required).")
        critical_failure = True

    # PyTorch & CUDA check
    try:
        import torch
        if torch.cuda.is_available():
            dev_name = torch.cuda.get_device_name(0)
            props = torch.cuda.get_device_properties(0)
            vram_gb = props.total_memory / (1024 ** 3)
            print(f"  [PASS] PyTorch CUDA acceleration: {dev_name} ({vram_gb:.2f} GB VRAM)")
        else:
            print("  [WARN] CUDA acceleration unavailable; PyTorch operations will run on CPU.")
            warning_count += 1
            if strict:
                critical_failure = True
    except ImportError:
        print("  [WARN] PyTorch is not installed in current environment; cannot probe CUDA.")
        warning_count += 1
        if strict:
            critical_failure = True

    # Core package imports
    core_packages = [
        ("cv2", "OpenCV"),
        ("torch", "PyTorch"),
        ("pygbif", "pygbif"),
        ("yaml", "PyYAML"),
    ]
    for pkg_module, pkg_display in core_packages:
        try:
            mod = _importlib.import_module(pkg_module)
            ver = getattr(mod, "__version__", "installed")
            print(f"  [PASS] Python package '{pkg_module}' ({pkg_display}): version {ver}")
        except ImportError:
            print(f"  [FAIL] Python package '{pkg_module}' ({pkg_display}) is NOT installed.")
            print(f"         Remediation: pip install {pkg_module}")
            critical_failure = True

    # -------------------------------------------------------------------------
    # 2. R Runtime & Statistical Packages
    # -------------------------------------------------------------------------
    print("\n[2/4] R Runtime & Statistical Packages:")
    rscript_bin = _shutil.which("Rscript")
    required_r_pkgs = ['Momocs', 'mclust', 'MorphoTools2', 'spatialRF', 'terra', 'tidyverse', 'optparse']
    optional_r_pkgs = ['RRPP', 'clustvarsel', 'geodata', 'ranger', 'adespatial', 'blockCV', 'ecospat']
    cran_r_pkgs = [
        "remotes", "tidyverse", "optparse", "mclust", "clustvarsel", "MorphoTools2",
        "RRPP", "vegan", "terra", "sf", "geodata", "soilDB", "spatialRF", "ranger",
        "adespatial", "blockCV", "ecospat", "ENMTools", "patchwork", "viridis",
    ]
    cran_pkg_str = ", ".join(f'"{p}"' for p in cran_r_pkgs)
    remediation_two_step = (
        f"Rscript -e 'install.packages(c({cran_pkg_str}), repos=\"https://cloud.r-project.org\"); "
        "if (!requireNamespace(\"Momocs\", quietly=TRUE)) remotes::install_github(\"MomX/Momocs\")'"
    )

    if not rscript_bin:
        print("  [FAIL] 'Rscript' binary not found in system $PATH.")
        print("         R (>= 4.3) is required for Momocs, MorphoTools2, and spatialRF.")
        print("         Remediation:")
        print("           1. Install R (>= 4.3) via system package manager (e.g., apt-get install r-base r-base-dev).")
        print("           2. Install required CRAN packages and Momocs (from GitHub):")
        print(f"              {remediation_two_step}")
        print("              Rscript -e 'remotes::install_github(\"MomX/Momocs\")'")
        critical_failure = True
    else:
        try:
            r_ver_proc = _subprocess.run(
                [rscript_bin, "-e", "cat(as.character(getRversion()))"],
                capture_output=True, text=True, timeout=5,
            )
            r_ver_text = r_ver_proc.stdout.strip()
            if r_ver_proc.returncode == 0 and r_ver_text and r_ver_text[0].isdigit():
                print(f"  [PASS] Rscript executable: {rscript_bin} (R {r_ver_text})")
            else:
                print(f"  [PASS] Rscript executable: {rscript_bin}")
        except Exception:
            print(f"  [PASS] Rscript executable: {rscript_bin}")

        r_probe_cmd = (
            "pkgs <- c('Momocs', 'mclust', 'MorphoTools2', 'spatialRF', 'terra', 'tidyverse', 'optparse'); "
            "opt_pkgs <- c('RRPP', 'clustvarsel', 'geodata', 'ranger', 'adespatial', 'blockCV', 'ecospat'); "
            "missing <- pkgs[!sapply(pkgs, requireNamespace, quietly=TRUE)]; "
            "opt_ok <- opt_pkgs[sapply(opt_pkgs, requireNamespace, quietly=TRUE)]; "
            "opt_miss <- opt_pkgs[!sapply(opt_pkgs, requireNamespace, quietly=TRUE)]; "
            "cat('OPT_INSTALLED:', paste(opt_ok, collapse=','), fill=TRUE); "
            "cat('OPT_MISSING:', paste(opt_miss, collapse=','), fill=TRUE); "
            "if(length(missing)>0) { cat('MISSING:', paste(missing, collapse=','), fill=TRUE); quit(status=1) } "
            "else { cat('ALL_INSTALLED', fill=TRUE); quit(status=0) }"
        )
        try:
            probe_proc = _subprocess.run(
                [rscript_bin, "-e", r_probe_cmd],
                capture_output=True, text=True, timeout=10,
            )
            combined_out = probe_proc.stdout + "\n" + probe_proc.stderr
            opt_installed: List[str] = list(optional_r_pkgs)
            opt_missing: List[str] = []
            for line in combined_out.splitlines():
                if "OPT_INSTALLED:" in line:
                    part = line.split("OPT_INSTALLED:", 1)[1].strip()
                    opt_installed = [p.strip() for p in part.split(",") if p.strip()]
                elif "OPT_MISSING:" in line:
                    part = line.split("OPT_MISSING:", 1)[1].strip()
                    opt_missing = [p.strip() for p in part.split(",") if p.strip()]

            if probe_proc.returncode == 0 and "ALL_INSTALLED" in probe_proc.stdout:
                print(f"  [PASS] Required R packages verified: {', '.join(required_r_pkgs)}")
                if not opt_missing:
                    print(f"  [PASS] Optional high-performance R packages verified: {', '.join(opt_installed)}")
                else:
                    print(
                        f"  [INFO] Optional high-performance R packages installed: "
                        f"{', '.join(opt_installed) if opt_installed else 'none'} "
                        f"(optional missing: {', '.join(opt_missing)})"
                    )
            else:
                missing_pkgs: List[str] = []
                for line in combined_out.splitlines():
                    if "MISSING:" in line:
                        pkgs_part = line.split("MISSING:", 1)[1].strip()
                        missing_pkgs = [p.strip() for p in pkgs_part.split(",") if p.strip()]
                        break
                if not missing_pkgs:
                    missing_pkgs = required_r_pkgs

                print(f"  [FAIL] Missing R package(s): {', '.join(missing_pkgs)}")
                print(
                    f"  [INFO] Optional high-performance R packages installed: "
                    f"{', '.join(opt_installed) if opt_installed else 'none'}"
                )
                print("         Remediation (install CRAN packages, then Momocs from GitHub):")
                print(f"           {remediation_two_step}")
                if "Momocs" in missing_pkgs:
                    print("         Momocs remediation (archived on CRAN): Rscript -e 'remotes::install_github(\"MomX/Momocs\")'")
                critical_failure = True
        except _subprocess.TimeoutExpired:
            print("  [FAIL] R package validation probe timed out (>10s).")
            critical_failure = True
        except Exception as exc:
            print(f"  [FAIL] Error executing R package probe: {exc}")
            critical_failure = True

    # -------------------------------------------------------------------------
    # 3. Directory & Configuration Verification
    # -------------------------------------------------------------------------
    print("\n[3/4] Directory & Configuration Verification:")
    cfg_file = Path(custom_config or (_proj_root / "config" / "config.yaml"))
    parsed_config: Optional[Dict[str, Any]] = None

    if not cfg_file.exists():
        print(f"  [FAIL] Pipeline configuration file not found at: {cfg_file}")
        critical_failure = True
    else:
        try:
            with open(cfg_file, "r", encoding="utf-8") as f:
                parsed_config = yaml.safe_load(f)
            if not isinstance(parsed_config, dict) or "paths" not in parsed_config:
                print(f"  [FAIL] Configuration at {cfg_file} is invalid or missing 'paths' section.")
                critical_failure = True
            else:
                rel_cfg = cfg_file.relative_to(_proj_root) if cfg_file.is_relative_to(_proj_root) else cfg_file
                print(f"  [PASS] Configuration file parsed & valid: {rel_cfg}")
        except Exception as exc:
            print(f"  [FAIL] Error parsing YAML configuration at {cfg_file}: {exc}")
            critical_failure = True

    target_dirs = [
        "data/raw_vouchers",
        "data/tables",
        "data/contours",
        "outputs/figures",
    ]
    for dir_rel in target_dirs:
        dir_path = _proj_root / dir_rel
        try:
            dir_path.mkdir(parents=True, exist_ok=True)
            probe_file = dir_path / ".perm_check_probe"
            probe_file.touch()
            probe_file.unlink()
            print(f"  [PASS] Directory exists and writable: {dir_rel}/")
        except Exception as exc:
            print(f"  [FAIL] Directory missing or write-permission denied: {dir_rel}/ ({exc})")
            critical_failure = True

    # -------------------------------------------------------------------------
    # 4. Model Checkpoint & Sub-Environments
    # -------------------------------------------------------------------------
    print("\n[4/4] Model Checkpoint & Sub-Environments:")
    venv_lm3_candidates = [
        _proj_root / "LeafMachine3" / ".venv",
        _proj_root / ".venv_LM3",
    ]
    venv_lm3 = next(
        (
            cand
            for cand in venv_lm3_candidates
            if cand.is_dir()
            and ((cand / "bin" / "python").exists() or (cand / "Scripts" / "python.exe").exists())
        ),
        None,
    )
    if venv_lm3 is not None:
        rel_lm3 = venv_lm3.relative_to(_proj_root) if venv_lm3.is_relative_to(_proj_root) else venv_lm3
        print(f"  [PASS] LeafMachine3 virtual environment found: {rel_lm3}/")
    else:
        try:
            ort_mod = _importlib.import_module("onnxruntime")
            ort_ver = getattr(ort_mod, "__version__", "installed")
            print(f"  [PASS] LeafMachine3 ONNX Runtime available in active environment (onnxruntime {ort_ver})")
        except ImportError:
            print(
                f"  [WARN] LeafMachine3 environment (LeafMachine3/.venv or .venv_LM3) not found at {venv_lm3_candidates[0]}. "
                "Run `bash setup_leafmachine3.sh` to provision LeafMachine3."
            )
            warning_count += 1
            if strict:
                critical_failure = True

    default_seg_weights = (
        "LeafMachine3/leafmachine3/models/"
        "lm3_leaf_segmenter__yolo26x_seg_1024.onnx"
    )
    model_rel = (
        (
            parsed_config.get("segmentation", {}).get("model_weights")
            or parsed_config.get("paths", {}).get("model_weights", default_seg_weights)
        )
        if parsed_config
        else default_seg_weights
    )
    model_path = Path(model_rel) if Path(model_rel).is_absolute() else (_proj_root / model_rel)
    if not (model_path.exists() and model_path.is_file() and model_path.stat().st_size > 0):
        try:
            from scripts.vision.configure_leafmachine3 import resolve_lm3_weights
            resolved_cand = resolve_lm3_weights(Path(model_rel).name, explicit_path=model_path)
            if resolved_cand.exists() and resolved_cand.is_file() and resolved_cand.stat().st_size > 0:
                model_path = resolved_cand
        except Exception:
            pass
    if not (model_path.exists() and model_path.is_file() and model_path.stat().st_size > 0):
        print(f"  [INFO] Model checkpoint missing at {model_path}. Invoking automated artifact manager...")
        try:
            cfg_obj = _PipelineConfig.from_yaml(cfg_file) if cfg_file.exists() else (parsed_config or {})
            ensured = _ensure_weights(cfg_obj)
            if isinstance(ensured, Path) and ensured.exists() and ensured.is_file() and ensured.stat().st_size > 0:
                model_path = ensured
        except Exception as exc:
            print(f"  [WARN] Automatic weights acquisition failed: {exc}")

    if model_path.exists() and model_path.is_file() and model_path.stat().st_size > 0:
        size_mb = model_path.stat().st_size / (1024 * 1024)
        print(f"  [PASS] Fine-tuned model checkpoint verified: {model_rel} ({size_mb:.2f} MB)")
    else:
        print(f"  [FAIL] Fine-tuned model checkpoint not found at: {model_path}")
        print("         Remediation: Run `python main.py download-weights` or place trained weights manually.")
        critical_failure = True

    # -------------------------------------------------------------------------
    # Readiness Verdict
    # -------------------------------------------------------------------------
    print("=" * 79)
    if critical_failure:
        print("Readiness Verdict: [FAIL] Environment is NOT ready for production.")
        print("Please resolve the [FAIL] dependencies listed above before starting runs.")
        print("=" * 79)
        return False
    else:
        if warning_count > 0:
            print(f"Readiness Verdict: [PASS] Ready with {warning_count} warning(s).")
        else:
            print("Readiness Verdict: [PASS] Complete pipeline environment verified & ready!")
        print("=" * 79)
        return True


def run_check_env(args: argparse.Namespace, cfg: Optional[PipelineConfig] = None) -> None:
    """CLI handler executing environment diagnostics and terminating with status code."""
    success = _main_attr("check_environment", check_environment)(args)
    sys.exit(0 if success else 1)


def run_download_weights(args: argparse.Namespace, cfg: PipelineConfig) -> None:
    """Downloads and verifies fine-tuned model checkpoint weights."""
    _logger = _main_attr("logger", logger)
    _ensure_weights = _main_attr("ensure_model_weights", ensure_model_weights)
    _logger.info("=== Checking & Downloading Model Weights Checkpoint ===")
    try:
        saved_path = _ensure_weights(cfg, force=getattr(args, "force", False))
        _logger.info(f"Model weights checkpoint ready at: {saved_path}")
    except Exception as exc:
        _logger.error(f"Failed to acquire model weights: {exc}")
        sys.exit(1)
