#!/usr/bin/env python3
"""
scripts/pipeline/step07_focal_parent_run.py
===========================================
Unique parent-directory generation, resume symlink management, checkpoint
detection, and JSON run-manifest writers for the six-stage focal species-complex
orchestrator.
"""

from __future__ import annotations

import csv
import json
import os
import re
import secrets
import subprocess
import sys
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, Mapping, Optional, Sequence, Tuple

PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from scripts.core.harvester import standardize_packera_taxon


def slugify_focal_taxon(taxon_name: str) -> str:
    clean = re.sub(r"[^A-Za-z0-9]+", "_", taxon_name.strip()).strip("_")
    return clean or "Packera_focal"


def default_focal_prefix(taxon_name: str) -> str:
    parts = taxon_name.strip().split()
    if len(parts) >= 2 and parts[0].lower() == "packera":
        return re.sub(r"[^a-z0-9]+", "_", parts[1].lower()).strip("_")
    return slugify_focal_taxon(taxon_name).lower()


def _extract_uid_hex6(run_id: str) -> str:
    m = re.search(r"_([0-9a-f]{6})(?:$|/)", run_id.strip())
    return m.group(1) if m else ""


def generate_focal_parent_run_id(prefix: str = "dubia") -> str:
    """Return a unique timestamped parent run ID matching '<prefix>_YYYYMMDD_<hex6>'."""
    clean = (prefix or "dubia").strip().strip("/") or "dubia"
    clean = re.sub(r"[^A-Za-z0-9_-]+", "_", clean).strip("_") or "dubia"
    return f"{clean}_{datetime.now().strftime('%Y%m%d')}_{secrets.token_hex(3)}"


def _six_stage_run_ids(
    prefix: str, discovery_run_id: Optional[str] = None
) -> Tuple[str, str, str, str, str, str]:
    """Return ordered (r1..r6) run IDs for nested ('parent/') or flat ('dubia') prefixes."""
    clean = prefix.strip()
    sep = "" if clean.endswith("/") else "_"
    r1 = f"{clean}{sep}01_eastern_discovery"
    r2 = f"{clean}{sep}02_complex_triage"
    r3 = f"{clean}{sep}03_intraspecific_clean"
    r4 = f"{clean}{sep}04_spatial_thin_field_audit"
    r5 = f"{clean}{sep}05_monographer_gold_audit"
    r6 = f"{clean}{sep}06_lineage_diagnosability"
    if discovery_run_id is not None and str(discovery_run_id).strip():
        r1 = str(discovery_run_id).strip()
    return r1, r2, r3, r4, r5, r6


def resolve_focal_parent_run_id(
    workspace_root: Path,
    focal_taxon: str,
    prefix: str = "",
    parent_run_id: Optional[str] = None,
    resume: bool = False,
) -> Tuple[str, str, str, Path]:
    """
    Resolve `(parent_run_id, unique_id, stage_prefix, default_output_dir)` so all
    stage folders and synthesis reports live inside a single parent folder under `runs/`.
    """
    focal_std = standardize_packera_taxon(focal_taxon) or focal_taxon.strip()
    slug = slugify_focal_taxon(focal_std)
    epithet = default_focal_prefix(focal_std)
    base_prefix = prefix.strip().strip("/") or epithet
    runs_root = (workspace_root / "runs").resolve()

    if parent_run_id is not None and str(parent_run_id).strip():
        clean = str(parent_run_id).strip().strip("/")
        if clean.startswith("runs/"):
            clean = clean[len("runs/"):].strip("/")
        if not clean:
            raise ValueError("Provided --parent-run-id resolved to an empty string.")
        low_c = clean.lower()
        is_direct_focal = (
            clean.endswith(f"/{slug}")
            or clean == slug
            or low_c == epithet
            or low_c.startswith(f"{epithet}_")
            or (base_prefix and (low_c == base_prefix.lower() or low_c.startswith(f"{base_prefix.lower()}_")))
        )
        if is_direct_focal:
            return clean, _extract_uid_hex6(clean), f"{clean}/", workspace_root / "runs" / clean
        nested = f"{clean}/{slug}"
        return clean, _extract_uid_hex6(clean), f"{nested}/", workspace_root / "runs" / clean / slug

    if resume:
        cand_links = [runs_root / f"latest_{base_prefix}", runs_root / f"latest_{epithet}"]
        latest_link = next((lk for lk in cand_links if lk.is_symlink() or lk.exists()), None)
        if latest_link is None:
            raise ValueError(
                f"Cannot --resume without --parent-run-id: {cand_links[0]} is missing."
            )
        raw_t = Path(os.readlink(latest_link)) if latest_link.is_symlink() else latest_link
        resolved_dir = (
            (latest_link.parent / raw_t).resolve() if not raw_t.is_absolute() else raw_t.resolve()
        )
        if not resolved_dir.is_dir() or runs_root not in resolved_dir.parents:
            raise ValueError(f"Cannot --resume: {latest_link} target {resolved_dir} is invalid.")
        rel_id = resolved_dir.relative_to(runs_root).as_posix()
        return rel_id, _extract_uid_hex6(rel_id), f"{rel_id}/", workspace_root / "runs" / rel_id

    # Backward-compatible check for summarize_only when an exact prefix folder already exists
    if prefix.strip():
        cand_existing = workspace_root / "runs" / base_prefix
        if cand_existing.is_dir() and any(
            (cand_existing / sub).is_dir()
            for sub in ("01_eastern_discovery", "02_complex_triage", "03_intraspecific_clean")
        ):
            return base_prefix, _extract_uid_hex6(base_prefix), f"{base_prefix}/", cand_existing

    fresh_id = generate_focal_parent_run_id(base_prefix)
    return fresh_id, _extract_uid_hex6(fresh_id), f"{fresh_id}/", workspace_root / "runs" / fresh_id


def update_latest_focal_symlink(
    workspace_root: Path, focal_taxon: str, parent_run_id: str
) -> None:
    """Point `runs/latest_<epithet>` atomically at `runs/<parent_run_id>` after ensuring the directory exists."""
    if os.environ.get("PLANTMATICS_NO_LATEST_SYMLINK") == "1" or "/" in parent_run_id.strip("/"):
        return
    focal_std = standardize_packera_taxon(focal_taxon) or focal_taxon.strip()
    epithet = default_focal_prefix(focal_std)
    runs_dir = workspace_root / "runs"
    target_dir = runs_dir / parent_run_id
    target_dir.mkdir(parents=True, exist_ok=True)
    link_path = runs_dir / f"latest_{epithet}"
    tmp_link = runs_dir / f".latest_{epithet}_{os.getpid()}.tmp"
    if tmp_link.exists() or tmp_link.is_symlink():
        tmp_link.unlink()
    tmp_link.symlink_to(parent_run_id)
    tmp_link.replace(link_path)


def stage_is_complete_focal(workspace_root: Path, run_id: str) -> bool:
    triage_csv = workspace_root / "runs" / run_id / "data" / "tables" / "triage_queue.csv"
    if not triage_csv.is_file() or triage_csv.stat().st_size == 0:
        return False
    with open(triage_csv, "r", encoding="utf-8", newline="") as fh:
        reader = csv.reader(fh)
        return bool(next(reader, None) and next(reader, None))


def _git_commit_sha_focal(workspace_root: Path) -> str:
    try:
        res = subprocess.run(
            ["git", "rev-parse", "HEAD"],
            cwd=str(workspace_root),
            capture_output=True,
            text=True,
            timeout=3,
            check=False,
        )
        return res.stdout.strip() or "unknown" if res.returncode == 0 else "unknown"
    except Exception:
        return "unknown"


def write_focal_run_manifest(
    *,
    output_dir: Path,
    parent_run_id: str,
    unique_id: str,
    focal_taxon: str,
    cli_arguments: Sequence[str],
    stage_run_ids: Mapping[str, str],
    stage_statuses: Mapping[str, str],
    workspace_root: Path = PROJECT_ROOT,
) -> Path:
    output_dir.mkdir(parents=True, exist_ok=True)
    manifest_path = output_dir / "focal_complex_run_manifest.json"
    payload: Dict[str, Any] = {
        "parent_run_id": parent_run_id,
        "unique_id": unique_id,
        "focal_taxon": focal_taxon,
        "created_date": datetime.now().astimezone().isoformat(),
        "git_commit_sha": _git_commit_sha_focal(workspace_root),
        "cli_arguments": list(cli_arguments),
        "stage_run_ids": dict(stage_run_ids),
        "stage_completion_statuses": dict(stage_statuses),
    }
    tmp_p = manifest_path.with_suffix(".json.tmp")
    with open(tmp_p, "w", encoding="utf-8", newline="\n") as fh:
        fh.write(json.dumps(payload, indent=2) + "\n")
        fh.flush()
        os.fsync(fh.fileno())
    tmp_p.replace(manifest_path)
    return manifest_path

