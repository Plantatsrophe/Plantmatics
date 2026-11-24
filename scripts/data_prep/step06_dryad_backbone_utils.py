#!/usr/bin/env python3
"""Dryad API authentication, archive routing, SHA-256 verification, and safe extraction for Packera backbones."""
from __future__ import annotations

import gzip
import hashlib
import importlib.util
import json
import os
import re
import shutil
import sys
import urllib.error
import urllib.parse
import urllib.request
import zipfile
from pathlib import Path
from typing import Any, Callable, Dict, List, Optional, Sequence, Tuple

_UA = "Plantmatics/1.0 (https://github.com/plantmatics; mailto:plantmatics@users.noreply.github.com)"
_DRYAD_DOI = "doi%3A10.5061%2Fdryad.bg79cnppn"
_DRYAD_ZIPS = ("speciestree_files.zip", "PhyTop_results.zip", "RelTime_trees.zip", "alignment_files.zip")
_DRYAD_EXACT_MAP: Dict[str, str] = {
    "18complete_packera_plastomes.fasta": "public_packera_18complete_plastomes.fasta",
    "111mapping_packera_plastomes.fasta": "public_packera_111mapping_plastomes.fasta",
    "mo_astral_lpp.tre": "public_packera_MO_Astral_lpp.tre",
    "astralpro_astral_lpp.tre": "public_packera_ASTRALPro_Astral_lpp.tre",
    "wastral_astral_lpp.tre": "public_packera_wASTRAL_Astral_lpp.tre",
    "111mapping_packera_plastomes.tre": "public_packera_111mapping_plastomes.tre",
    "mo_concat_best.fasta": "public_packera_MO_concat_best.fasta",
    "1-to-1_concat_best.fasta": "public_packera_1to1_concat_best.fasta",
}


class _CrossHostAuthStripRedirectHandler(urllib.request.HTTPRedirectHandler):
    """Strip Authorization headers on redirects across domains (e.g. Dryad -> S3)."""

    def redirect_request(self, req: Any, fp: Any, code: int, msg: str, headers: Any, newurl: str) -> Any:
        new_req = super().redirect_request(req, fp, code, msg, headers, newurl)
        if new_req is not None and urllib.parse.urlparse(req.full_url).netloc.lower() != urllib.parse.urlparse(newurl).netloc.lower():
            for hdrs in (new_req.headers, new_req.unredirected_hdrs):
                for k in [k for k in hdrs if k.lower() == "authorization"]:
                    hdrs.pop(k, None)
        return new_req


def _default_http_request(
    url: str, *, headers: Optional[Dict[str, str]] = None, timeout: int = 60,
    opener: Optional[urllib.request.OpenerDirector] = None, data: Optional[bytes] = None, dest: Optional[Path] = None,
) -> bytes:
    mod = sys.modules.get("scripts.data_prep.download_public_packera_sequences")
    fn = getattr(mod, "_http_request", None) if mod is not None else None
    if callable(fn):
        return fn(url, headers=headers, timeout=timeout, opener=opener, data=data, dest=dest)
    req = urllib.request.Request(url, data=data, headers={"User-Agent": _UA, **(headers or {})}, method="POST" if data is not None else "GET")
    open_fn = opener.open if opener is not None else urllib.request.urlopen
    with open_fn(req, timeout=timeout) as resp:
        if dest is not None:
            with open(dest, "wb") as out:
                shutil.copyfileobj(resp, out, length=65536)
            return b""
        raw = resp.read()
    return gzip.decompress(raw) if raw[:2] == b"\x1f\x8b" else raw


def _load_dotenv_defaults(paths: Optional[Sequence[Path]] = None) -> None:
    """Load missing KEY=VALUE secrets from `.env` files without overriding active environment variables."""
    if os.environ.get("PLANTMATICS_SKIP_DOTENV", "").strip() == "1":
        return
    for env_path in (tuple(paths) if paths is not None else (Path(__file__).resolve().parents[2] / ".env", Path.home() / ".env")):
        if not env_path.is_file():
            continue
        for raw in env_path.read_text(encoding="utf-8", errors="ignore").splitlines():
            ln = raw.strip()
            if not ln or ln.startswith("#"):
                continue
            ln = ln[7:].strip() if ln.startswith("export ") else ln
            if "=" in ln:
                k, v = (p.strip() for p in ln.split("=", 1))
                v = v[1:-1] if len(v) >= 2 and v[0] == v[-1] and v[0] in ("'", '"') else v
                if k and k not in os.environ:
                    os.environ[k] = v


def _resolve_dryad_bearer_token(http_request_fn: Optional[Callable[..., bytes]] = None) -> str:
    """Resolve a Dryad API bearer token from environment/dotenv or OAuth2 client credentials."""
    _load_dotenv_defaults()
    if tok := os.environ.get("DRYAD_API_TOKEN", "").strip():
        return tok
    cid, csec = os.environ.get("DRYAD_CLIENT_ID", "").strip(), os.environ.get("DRYAD_CLIENT_SECRET", "").strip()
    if cid and csec:
        req_fn = http_request_fn or _default_http_request
        body = urllib.parse.urlencode({"client_id": cid, "client_secret": csec, "grant_type": "client_credentials"}).encode("utf-8")
        raw = req_fn("https://datadryad.org/oauth/token", headers={"Content-Type": "application/x-www-form-urlencoded"}, data=body)
        return str(json.loads(raw.decode("utf-8")).get("access_token") or "").strip()
    return ""


def _route_dryad_file(src: Path, zip_name: str) -> Optional[str]:
    """Map an extracted Dryad member file to its canonical `public_packera_*` destination basename."""
    low, zlow = src.name.lower(), zip_name.lower()
    if "__macosx" in str(src).lower() or src.name.startswith("._"):
        return None
    if low in _DRYAD_EXACT_MAP:
        return _DRYAD_EXACT_MAP[low]
    if ("phytop" in zlow or "phytop" in low) and low.endswith((".tsv", ".csv", ".tree", ".tre", ".nwk")):
        stem_clean = re.sub(r"_phytop(\.info)?$", "", src.stem, flags=re.IGNORECASE)
        method = (
            "ASTRALPro" if ("astralpro" in low or "astral-pro" in low)
            else "wASTRAL" if "wastral" in low
            else "HybPiper-df" if "hybpiper" in low
            else "MO" if re.search(r"(^|[^a-z0-9])mo([^a-z0-9]|$)", low)
            else (re.sub(r"[^A-Za-z0-9_-]+", "_", stem_clean).strip("_") or "table")
        )
        ext = ".tre" if low.endswith((".tree", ".tre", ".nwk")) else ".tsv"
        return f"public_packera_PhyTop_{method}{ext}"
    if ("reltime" in zlow or "reltime" in low) and low.endswith((".tre", ".nwk", ".tree")):
        return f"public_packera_RelTime_{re.sub(r'[^A-Za-z0-9_.-]+', '_', src.stem).strip('_') or 'tree'}.tre"
    return None


def _fetch_dryad_packera_files(
    genomics_dir: Path, tmpd: Path, dryad_zip_dir: Optional[Path], force: bool,
    *,
    http_request_fn: Optional[Callable[..., bytes]] = None,
    resolve_dest_fn: Optional[Callable[[Path, str], Path]] = None,
    atomic_save_fn: Optional[Callable[..., Tuple[Path, bool]]] = None,
    token_fn: Optional[Callable[[], str]] = None,
) -> Tuple[str, int]:
    """Retrieve, hash-verify, zip-slip-guard, and route Packera Dryad archives into `genomics_dir`."""
    req_fn = http_request_fn or _default_http_request
    if resolve_dest_fn is None or atomic_save_fn is None:
        from scripts.data_prep import download_public_packera_sequences as _dl_pub
        resolve_dest_fn = resolve_dest_fn or _dl_pub._resolve_public_dest
        atomic_save_fn = atomic_save_fn or _dl_pub._atomic_save
    manifest_files: List[Dict[str, Any]] = []
    default_m = resolve_dest_fn(genomics_dir, "public_packera_dryad_manifest.json")
    manifest_path = next((p for p in (default_m, genomics_dir / "public_packera_dryad_manifest.json") if p.is_file()), default_m)
    try:
        v_raw = req_fn(f"https://datadryad.org/api/v2/datasets/{_DRYAD_DOI}/versions", timeout=45)
        if v_list := json.loads(v_raw.decode("utf-8")).get("_embedded", {}).get("stash:versions", []):
            v_href = v_list[-1]["_links"]["self"]["href"]
            f_raw = req_fn(urllib.parse.urljoin("https://datadryad.org", v_href.rstrip("/") + "/files?per_page=100"), timeout=45)
            for item in json.loads(f_raw.decode("utf-8")).get("_embedded", {}).get("stash:files", []):
                dl_href = item.get("_links", {}).get("stash:download", {}).get("href") or (f"/api/v2/files/{item['id']}/download" if item.get("id") is not None else "")
                manifest_files.append({k: item.get(k, "") for k in ("id", "path", "size", "digestType", "digest")} | {"download_href": dl_href})
            atomic_save_fn(genomics_dir, "public_packera_dryad_manifest.json", json.dumps({"doi": "10.5061/dryad.bg79cnppn", "files": manifest_files}, indent=2).encode("utf-8"), force=force)
    except Exception as exc:
        print(f"[WARN] Dryad manifest fetch failed ({exc}).", file=sys.stderr)
        if manifest_path.is_file():
            try:
                manifest_files = json.loads(manifest_path.read_text(encoding="utf-8")).get("files", [])
            except Exception:
                manifest_files = []

    by_name = {Path(str(m.get("path", ""))).name: m for m in manifest_files if m.get("path")}
    staged_zips: List[Path] = []
    status = "manifest_only_auth_required"
    if dryad_zip_dir is not None and dryad_zip_dir.is_dir():
        for zname in _DRYAD_ZIPS:
            cand = next((p for p in dryad_zip_dir.iterdir() if p.is_file() and p.name.lower() == zname.lower()), None)
            if cand is not None:
                shutil.copy2(cand, tmpd / zname)
                staged_zips.append(tmpd / zname)
        if staged_zips:
            status = "local_zip_dir"

    if not staged_zips:
        try:
            tok = token_fn() if token_fn is not None else _resolve_dryad_bearer_token(http_request_fn=req_fn)
            if tok:
                opener = urllib.request.build_opener(_CrossHostAuthStripRedirectHandler())
                for zname in _DRYAD_ZIPS:
                    if (meta := by_name.get(zname)) and meta.get("download_href"):
                        req_fn(urllib.parse.urljoin("https://datadryad.org", str(meta["download_href"])), headers={"Authorization": f"Bearer {tok}"}, timeout=180, opener=opener, dest=tmpd / zname)
                        staged_zips.append(tmpd / zname)
                if staged_zips:
                    status = "downloaded"
        except urllib.error.HTTPError as exc:
            if exc.code in (401, 403):
                staged_zips.clear()
            else:
                raise

    if not staged_zips:
        pw_avail = importlib.util.find_spec("playwright") is not None
        print(f'[INFO] Dryad archives require authentication (playwright_installed={pw_avail}). Set DRYAD_API_TOKEN and run: export DRYAD_API_TOKEN="<token>" && .venv/bin/python main.py download-public-genomics --genomics-dir data/genomics or pass --dryad-zip-dir <DIR>.')
        return status, 0

    n_extracted = 0
    for zp in staged_zips:
        meta = by_name.get(zp.name, {})
        if str(meta.get("digestType", "")).lower() in ("sha-256", "sha256") and meta.get("digest"):
            h = hashlib.sha256()
            with open(zp, "rb") as fh:
                for chunk in iter(lambda: fh.read(65536), b""):
                    h.update(chunk)
            if h.hexdigest().lower() != str(meta["digest"]).strip().lower():
                raise ValueError(f"SHA-256 mismatch for {zp.name}: expected {meta['digest']}, got {h.hexdigest()}")
        unpack_dir = tmpd / f"unpack_{zp.stem}"
        unpack_dir.mkdir(parents=True, exist_ok=True)
        unpack_res = unpack_dir.resolve()
        with zipfile.ZipFile(zp, "r") as zf:
            for member in zf.infolist():
                target = (unpack_dir / member.filename).resolve()
                if target != unpack_res and unpack_res not in target.parents:
                    raise ValueError(f"Unsafe zip-slip entry in {zp.name}: {member.filename!r}")
            zf.extractall(unpack_dir)
        for fp in sorted(unpack_dir.rglob("*")):
            if fp.is_file() and (out_name := _route_dryad_file(fp, zp.name)):
                atomic_save_fn(genomics_dir, out_name, fp, force=force)
                n_extracted += 1
    return status, n_extracted


_obtain_dryad_token = _resolve_dryad_bearer_token
_fetch_dryad = _fetch_dryad_packera_files
