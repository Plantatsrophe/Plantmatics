"""Run-directory lifecycle, `.env` credential loading, and atomic run metadata persistence."""

from __future__ import annotations

import csv
from datetime import datetime
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
from typing import Any, Dict, List, Optional, Sequence, Set, Tuple, Union

PROJECT_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_CONFIG_PATH = PROJECT_ROOT / "config" / "config.yaml"
RUNS_ROOT = PROJECT_ROOT / "runs"


def load_dotenv_defaults(paths: Optional[Sequence[Path]] = None) -> None:
    """Load missing KEY=VALUE credentials from local `.env` files without overriding existing env vars."""
    if os.environ.get("PLANTMATICS_SKIP_DOTENV", "").strip() == "1":
        return
    candidates = tuple(paths) if paths is not None else (PROJECT_ROOT / ".env", Path.home() / ".env")
    for env_path in candidates:
        if not env_path.is_file():
            continue
        try:
            lines = env_path.read_text(encoding="utf-8").splitlines()
        except OSError:
            continue
        for raw_line in lines:
            line = raw_line.strip()
            if not line or line.startswith("#"):
                continue
            if line.startswith("export "):
                line = line[7:].strip()
            if "=" not in line:
                continue
            key, val = line.split("=", 1)
            key, val = key.strip(), val.strip()
            if len(val) >= 2 and val[0] == val[-1] and val[0] in ("'", '"'):
                val = val[1:-1]
            if key and key not in os.environ:
                os.environ[key] = val


load_dotenv_defaults()

# Path fields in PathsConfig that are scoped to the active run directory
# (shared read-only caches like raw_vouchers_dir, quarantine_dir, models_dir, and model_weights remain at PROJECT_ROOT)
RUN_SCOPED_PATH_FIELDS: Set[str] = {
    "curated_vouchers_csv",
    "failed_qc_vouchers_csv",
    "extraction_manifest_csv",
    "contours_dir",
    "masks_dir",
    "leaf_efa_harmonics_csv",
    "morphometric_flags_csv",
    "cda_biplot_png",
    "revision_report_md",
    "gmm_report_csv",
    "summary_log",
    "tables_dir",
    "figures_dir",
}


def generate_default_run_id() -> str:
    """Generate a timestamped run identifier: run_YYYYMMDD_HHMMSS."""
    return f"run_{datetime.now().strftime('%Y%m%d_%H%M%S')}"


def resolve_canonical_run_id(
    run_dir: Optional[Path] = None,
    explicit_run_id: Optional[str] = None,
    workspace_root: Path = PROJECT_ROOT,
) -> str:
    """
    Resolve a canonical POSIX run identifier, preserving nested subpaths under `runs/`
    (e.g., `eastern_na_20261114_a7f3c2/Packera_dubia/03_intraspecific_clean`).
    """
    ws_root = Path(workspace_root)
    runs_root = (ws_root / "runs").resolve()

    if explicit_run_id is not None:
        clean = str(explicit_run_id).strip().replace("\\", "/")
        while clean.startswith("./"):
            clean = clean[2:]
        if clean.startswith("runs/"):
            clean = clean[len("runs/"):]
        clean = clean.rstrip("/")
        if clean:
            if Path(clean).is_absolute():
                tp = Path(clean).resolve()
                if runs_root in tp.parents:
                    return tp.relative_to(runs_root).as_posix()
                if run_dir is None:
                    return os.environ.get("RUN_ID", "").strip() or tp.name
            else:
                return clean

    if run_dir is not None:
        rd = Path(run_dir).resolve()
        if runs_root in rd.parents:
            return rd.relative_to(runs_root).as_posix()
        return os.environ.get("RUN_ID", "").strip() or rd.name

    return os.environ.get("RUN_ID", "").strip() or "canonical"


def get_active_run_dir(
    run_id: Optional[str] = None,
    workspace_root: Path = PROJECT_ROOT,
) -> Optional[Path]:
    """
    Return the resolved run directory if `run_id` is supplied or `PLANTMATICS_RUN_DIR`
    (or `PACKERA_RUN_DIR`) is set in the environment; otherwise return None (default legacy/test mode).
    """
    if run_id is not None and str(run_id).strip():
        clean_id = str(run_id).strip()
        p = Path(clean_id)
        if p.is_absolute():
            return p.resolve()
        if p.parts and p.parts[0] == "runs":
            return (workspace_root / p).resolve()
        return (workspace_root / "runs" / clean_id).resolve()

    env_run_dir = (
        os.environ.get("PLANTMATICS_RUN_DIR", "").strip()
        or os.environ.get("PACKERA_RUN_DIR", "").strip()
    )
    if env_run_dir:
        p = Path(env_run_dir)
        return p.resolve() if p.is_absolute() else (workspace_root / p).resolve()

    return None


def initialize_run_directory(
    run_id: Optional[str] = None,
    workspace_root: Path = PROJECT_ROOT,
    update_latest_symlink: bool = True,
    export_env: bool = True,
) -> Path:
    """
    Create the standardized self-contained run directory structure under `runs/<run_id>/`:
      - data/masks
      - data/contours
      - data/cropped_patches/capitula
      - data/tables
      - outputs/{figures,reports}
      - logs
    Also maintains the `runs/latest` symlink pointing to the active run directory.
    """
    run_dir = get_active_run_dir(run_id=run_id, workspace_root=workspace_root)
    if run_dir is None:
        eff_id = os.environ.get("RUN_ID", "").strip() or generate_default_run_id()
        run_dir = (workspace_root / "runs" / eff_id).resolve()

    subdirs = [
        run_dir / "data" / "masks",
        run_dir / "data" / "contours",
        run_dir / "data" / "cropped_patches" / "capitula",
        run_dir / "data" / "tables",
        run_dir / "outputs" / "figures",
        run_dir / "outputs" / "reports",
        run_dir / "logs",
    ]
    for d in subdirs:
        d.mkdir(parents=True, exist_ok=True)

    if update_latest_symlink:
        runs_parent = workspace_root / "runs"
        runs_parent.mkdir(parents=True, exist_ok=True)
        latest_link = runs_parent / "latest"
        try:
            tmp_link = runs_parent / f".latest_{os.getpid()}.tmp"
            if tmp_link.exists() or tmp_link.is_symlink():
                tmp_link.unlink()
            rd_resolved = run_dir.resolve()
            runs_parent_resolved = runs_parent.resolve()
            target_rel = (
                rd_resolved.relative_to(runs_parent_resolved).as_posix()
                if runs_parent_resolved in rd_resolved.parents
                else str(run_dir)
            )
            tmp_link.symlink_to(target_rel)
            tmp_link.replace(latest_link)
        except OSError:
            pass

    if export_env:
        os.environ["PLANTMATICS_RUN_DIR"] = str(run_dir)
        os.environ["PACKERA_RUN_DIR"] = str(run_dir)
        os.environ["RUN_ID"] = resolve_canonical_run_id(
            run_dir=run_dir,
            explicit_run_id=run_id,
            workspace_root=workspace_root,
        )

    return run_dir


def ensure_run_curated_vouchers_snapshot(
    run_dir: Path,
    workspace_root: Path = PROJECT_ROOT,
) -> Path:
    """
    Ensure `${RUN_DIR}/data/tables/curated_vouchers.csv` exists before Step 02.
    If absent in the run directory, copy the shared `data/tables/curated_vouchers.csv`
    into the run directory atomically so the run folder is a self-contained snapshot.
    """
    run_vouchers = run_dir / "data" / "tables" / "curated_vouchers.csv"
    if run_vouchers.exists() and run_vouchers.stat().st_size > 0:
        return run_vouchers

    shared_vouchers = workspace_root / "data" / "tables" / "curated_vouchers.csv"
    if shared_vouchers.exists() and shared_vouchers.stat().st_size > 0:
        run_vouchers.parent.mkdir(parents=True, exist_ok=True)
        tmp_vouchers = run_vouchers.with_suffix(".csv.tmp")
        shutil.copy2(shared_vouchers, tmp_vouchers)
        tmp_vouchers.replace(run_vouchers)

    return run_vouchers


def _count_csv_rows(csv_path: Path) -> int:
    """Count data rows (excluding header) in a CSV file."""
    if not csv_path.exists() or csv_path.stat().st_size == 0:
        return 0
    try:
        with open(csv_path, "r", encoding="utf-8", errors="replace") as f:
            reader = csv.reader(f)
            total = sum(1 for _ in reader)
            return max(0, total - 1)
    except Exception:
        return 0


def write_run_metadata(
    run_dir: Path,
    run_id: Optional[str] = None,
    start_time: Optional[str] = None,
    end_time: Optional[str] = None,
    cli_args: Optional[List[str]] = None,
    workspace_root: Path = PROJECT_ROOT,
    status: Optional[str] = None,
    exit_code: Optional[int] = None,
) -> Path:
    """
    Atomically write `${RUN_DIR}/outputs/reports/run_metadata.json` recording:
      - run_id
      - run_dir
      - start_timestamp / end_timestamp
      - status / exit_code
      - git_commit_sha
      - cli_arguments
      - table_row_counts (for all generated CSV tables in data/tables and outputs/reports)
    """
    reports_dir = run_dir / "outputs" / "reports"
    reports_dir.mkdir(parents=True, exist_ok=True)
    metadata_path = reports_dir / "run_metadata.json"

    existing: Dict[str, Any] = {}
    if metadata_path.exists():
        try:
            with open(metadata_path, "r", encoding="utf-8") as f:
                loaded = json.load(f)
                if isinstance(loaded, dict):
                    existing = loaded
        except Exception:
            existing = {}

    git_sha = existing.get("git_commit_sha", "unknown")
    if git_sha == "unknown":
        try:
            res = subprocess.run(
                ["git", "rev-parse", "HEAD"],
                cwd=str(workspace_root),
                capture_output=True,
                text=True,
                timeout=3,
                check=False,
            )
            if res.returncode == 0 and res.stdout.strip():
                git_sha = res.stdout.strip()
        except Exception:
            pass

    now_iso = datetime.now().astimezone().isoformat()
    eff_start = start_time or existing.get("start_timestamp") or now_iso
    eff_end = end_time or now_iso
    eff_run_id = (
        resolve_canonical_run_id(
            run_dir=run_dir,
            explicit_run_id=run_id,
            workspace_root=workspace_root,
        )
        if (run_id is not None and str(run_id).strip())
        else (
            existing.get("run_id")
            or resolve_canonical_run_id(
                run_dir=run_dir,
                explicit_run_id=None,
                workspace_root=workspace_root,
            )
        )
    )
    eff_cli = cli_args if cli_args is not None else existing.get("cli_arguments", list(sys.argv))
    if status is not None:
        eff_status = str(status)
    elif exit_code is not None:
        eff_status = "completed" if int(exit_code) == 0 else "failed"
    elif end_time is not None:
        eff_status = "completed"
    else:
        eff_status = existing.get("status", "running")

    if exit_code is not None:
        eff_exit_code: Optional[int] = int(exit_code)
    elif status == "completed" or (status is None and end_time is not None):
        eff_exit_code = 0
    else:
        raw_existing_exit = existing.get("exit_code")
        eff_exit_code = int(raw_existing_exit) if raw_existing_exit is not None else None

    table_counts: Dict[str, int] = {}
    for search_dir in (run_dir / "data" / "tables", run_dir / "outputs" / "reports"):
        if search_dir.exists():
            for csv_file in sorted(search_dir.glob("*.csv")):
                rel_key = str(csv_file.relative_to(run_dir))
                table_counts[rel_key] = _count_csv_rows(csv_file)
                table_counts[csv_file.name] = table_counts[rel_key]

    payload = {
        "run_id": eff_run_id,
        "run_dir": str(run_dir),
        "start_timestamp": eff_start,
        "end_timestamp": eff_end,
        "status": eff_status,
        "exit_code": eff_exit_code,
        "git_commit_sha": git_sha,
        "cli_arguments": eff_cli,
        "table_row_counts": table_counts,
    }

    tmp_path = metadata_path.with_suffix(".json.tmp")
    with open(tmp_path, "w", encoding="utf-8") as f:
        json.dump(payload, f, indent=2, sort_keys=False)
        f.write("\n")
        f.flush()
        os.fsync(f.fileno())
    tmp_path.replace(metadata_path)
    return metadata_path
