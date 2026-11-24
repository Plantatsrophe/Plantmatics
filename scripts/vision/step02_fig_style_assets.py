#!/usr/bin/env python3
"""
Step 02 figure styling, ICN botanical mathtext formatting, atomic 300-DPI PNG
export, path resolution, closed-mask validation, and LeafMachine3 cutout loading.
"""

from __future__ import annotations

import math
import os
import re
import sys
from pathlib import Path
from typing import Any, Dict, Optional, Sequence, Set, Tuple

import cv2
import numpy as np
import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

import scripts.vision.lm3_compat  # noqa: F401

# 12-hue Tol + Okabe-Ito colorblind-safe categorical palette
TOL_OKABE_12: Tuple[str, ...] = (
    "#0072B2", "#D55E00", "#009E73", "#CC79A7",  # Blue, Vermilion, Bluish green, Reddish purple
    "#E69F00", "#56B4E9", "#332288", "#882255",  # Orange, Sky blue, Indigo, Wine
    "#44AA99", "#999933", "#888888", "#117733",  # Teal, Olive, Slate, Forest green
)

# 5 unfilled/hollow marker geometries -> 12 x 5 = 60 unique taxon combinations
HOLLOW_MARKERS_5: Tuple[str, ...] = ("o", "^", "D", "s", "v")

ICN_ROMAN_TOKENS: Set[str] = {
    "subsp.", "subsp", "ssp.", "ssp", "var.", "var", "f.", "forma",
    "×", "x", "aff.", "aff", "cf.", "cf", "sp.", "sp", "s.s.", "s.l.", "ined.",
}

_ASSET_ANCHOR_DIRS = ("masks", "contours", "cropped_patches", "raw_vouchers")
_NULL_TAXON_STRINGS: Set[str] = {"none", "nan", "null", "na", "unknown", "unassigned_packera"}
_INFRA_RANK_RE = re.compile(r"\b(subsp\.?|ssp\.?|var\.?|f\.?|forma)\s+([a-z][a-z-]+)\b", flags=re.I)
_ICN_RANK_MAP: Dict[str, str] = {
    "x": "×", "×": "×", "subsp": "subsp.", "ssp": "subsp.", "subsp.": "subsp.",
    "ssp.": "subsp.", "var": "var.", "var.": "var.", "f": "f.", "f.": "f.", "forma": "f.",
}


def _resolve_facade_hook(name: str, local_obj: Any, orig_obj: Any) -> Any:
    """Honor monkeypatches applied either to the `step02_figures` facade or local submodule."""
    if local_obj is not orig_obj:
        return local_obj
    facade = sys.modules.get("scripts.vision.step02_figures")
    if facade is not None:
        patched = getattr(facade, name, orig_obj)
        if patched is not orig_obj:
            return patched
    return local_obj


def _atomic_save_matplotlib_fig(fig: Any, out_base_path: Path, dpi: int = 300) -> Path:
    """Atomically save a Matplotlib figure strictly as a 300-DPI PNG and remove stale PDF siblings."""
    out_base = Path(out_base_path)
    out_base.parent.mkdir(parents=True, exist_ok=True)
    png_path = out_base.with_suffix(".png")
    stale_pdf = out_base.with_suffix(".pdf")
    tmp_png = out_base.with_name(f".{out_base.stem}_{os.getpid()}_tmp.png")
    try:
        fig.savefig(str(tmp_png), format="png", dpi=int(dpi), bbox_inches="tight", facecolor="white")
        tmp_png.replace(png_path)
        if stale_pdf.exists():
            stale_pdf.unlink(missing_ok=True)
    finally:
        if tmp_png.exists():
            tmp_png.unlink(missing_ok=True)
    return png_path


def _extract_infra_suffix(text: str) -> str:
    """Parse an explicit botanical infraspecific rank + epithet suffix from text."""
    m_inf = _INFRA_RANK_RE.search(text)
    if not m_inf:
        return ""
    rk = m_inf.group(1).lower().rstrip(".")
    rk_std = "subsp." if rk in {"subsp", "ssp"} else ("var." if rk == "var" else "f.")
    return f" {rk_std} {m_inf.group(2).lower()}"


def _standardize_packera_species_py(raw_val: Any, infraspecific_val: Any = None) -> str:
    """Standardize Packera species, nothospecies, infraspecific ranks (subsp./var./f.), or morphotype IDs."""
    if raw_val is None or pd.isna(raw_val):
        return "Unassigned_Packera"
    s_raw = str(raw_val).strip()
    if not s_raw or s_raw.lower() in _NULL_TAXON_STRINGS:
        return "Unassigned_Packera"
    if re.match(r"^(?:Morphotype|Ecotype|Complex)_[A-Za-z0-9_.-]+$", s_raw):
        return s_raw

    # Extract explicit infraspecific rank from raw string or companion column before binomial standardization
    infra_suffix = ""
    if infraspecific_val is not None and not pd.isna(infraspecific_val):
        inf_str = str(infraspecific_val).strip()
        if inf_str and inf_str.lower() not in {"none", "nan", "null", "na"}:
            infra_suffix = _extract_infra_suffix(inf_str)
            if not infra_suffix and re.match(r"^[a-z][a-z-]+$", inf_str.lower()):
                infra_suffix = f" var. {inf_str.lower()}"
    if not infra_suffix:
        infra_suffix = _extract_infra_suffix(s_raw)

    base_std = None
    try:
        from scripts.core.harvester import standardize_packera_taxon

        std = standardize_packera_taxon(s_raw)
        if std and str(std).strip() and str(std).strip().lower() not in _NULL_TAXON_STRINGS:
            base_std = str(std).strip()
    except Exception:
        base_std = None

    if base_std is None:
        m_binom = re.match(r"^([A-Za-z][a-z]+)\s+(?:([x×])\s+)?([a-z][a-z-]+)\b", s_raw.split("(")[0].strip())
        if m_binom:
            hyb = "× " if m_binom.group(2) else ""
            base_std = f"{m_binom.group(1).capitalize()} {hyb}{m_binom.group(3).lower()}"

    if base_std and re.search(r"(?:^|\s)[×x]\s+[a-z]", s_raw) and "×" not in base_std:
        b_parts = base_std.split()
        if len(b_parts) >= 2:
            base_std = f"{b_parts[0]} × {' '.join(b_parts[1:])}"

    if base_std:
        if infra_suffix and infra_suffix.split()[-1] not in base_std.split()[1:2]:
            return f"{base_std}{infra_suffix}"
        # Retain autonym if explicitly present in raw input
        if infra_suffix and re.search(r"\b(subsp\.?|ssp\.?|var\.?|f\.?)\s+", s_raw, flags=re.I):
            return f"{base_std}{infra_suffix}"
        return base_std
    return "Unassigned_Packera"


def format_icn_taxon_mathtext(label_raw: Any, abbreviate_genus: bool = False) -> str:
    """Format a botanical name using ICN typographic conventions in Matplotlib mathtext."""
    s = str(label_raw or "").strip()
    if not s:
        return "Unknown"
    if re.match(r"^(?:Morphotype|Ecotype|Complex|Unassigned|Cluster|Unknown)\b", s, flags=re.I):
        return s.replace("_", " ")
    tokens = s.split()
    if len(tokens) < 2 or not re.match(r"^[A-Z][a-z]+$", tokens[0]):
        return s.replace("_", " ")

    out_parts = []
    for idx, tok in enumerate(tokens):
        low = tok.lower()
        if idx == 0:
            g_txt = f"{tok[0]}." if abbreviate_genus else tok
            out_parts.append(rf"$\it{{{g_txt}}}$")
        elif low in ICN_ROMAN_TOKENS:
            out_parts.append(_ICN_RANK_MAP.get(low, tok))
        elif re.match(r"^[a-z][a-z-]+$", tok):
            out_parts.append(rf"$\it{{{tok}}}$")
        else:
            out_parts.append(tok)
    return " ".join(out_parts)


def get_combinatorial_taxon_style(taxa_seq: Sequence[str]) -> Dict[str, Tuple[str, str]]:
    """Map an arbitrary number of taxa (1..60+) to distinct (color, hollow_marker) pairs."""
    unique_taxa = list(dict.fromkeys(str(t) for t in taxa_seq if str(t).strip()))
    n_c, n_m = len(TOL_OKABE_12), len(HOLLOW_MARKERS_5)
    return {
        tx: (
            TOL_OKABE_12[idx % n_c],
            HOLLOW_MARKERS_5[(idx // n_c) % n_m] if len(unique_taxa) > 8 else HOLLOW_MARKERS_5[idx % n_m],
        )
        for idx, tx in enumerate(unique_taxa)
    }


def _parse_bool_val(val: Any, default: bool = False) -> bool:
    """Parse boolean flags exported by Python or R ('TRUE'/'FALSE', True/False, 1/0)."""
    if val is None or (isinstance(val, float) and np.isnan(val)):
        return default
    if isinstance(val, (bool, np.bool_)):
        return bool(val)
    s = str(val).strip().upper()
    if s in {"TRUE", "T", "1", "1.0", "YES"}:
        return True
    return False if s in {"FALSE", "F", "0", "0.0", "NO"} else default


def _resolve_asset_path(p_str: Any, out_dir: Path) -> Optional[Path]:
    """Resolve mask, contour, crop, or raw voucher scan path across run-scoped and canonical stores."""
    raw = str(p_str or "").strip()
    if not raw or raw.lower() in {"nan", "none", "na"}:
        return None
    p = Path(raw)
    if p.exists():
        return p
    for base in (out_dir.parent.parent, out_dir.parent, PROJECT_ROOT):
        cand = base / raw
        if cand.exists():
            return cand
    parts = p.parts
    anchor_idx = next((i for i in range(len(parts) - 1, -1, -1) if parts[i] in _ASSET_ANCHOR_DIRS), None)
    if anchor_idx is not None and anchor_idx < len(parts) - 1:
        rel_tail = Path(*parts[anchor_idx:])
        for data_root in (out_dir.parent.parent / "data", out_dir.parent / "data", PROJECT_ROOT / "data"):
            cand = data_root / rel_tail
            if cand.exists():
                return cand
    return None


_ORIG_RESOLVE_ASSET_PATH = _resolve_asset_path


def _is_valid_smooth_closed_mask(m_path_str: Any, out_dir: Path) -> bool:
    """Verify binary mask has unclipped borders, coherent lamina islands, and clean geometry."""
    resolve_fn = _resolve_facade_hook("_resolve_asset_path", _resolve_asset_path, _ORIG_RESOLVE_ASSET_PATH)
    resolved = resolve_fn(m_path_str, out_dir)
    if resolved is None:
        return not bool(str(m_path_str or "").strip())
    mask_img = cv2.imread(str(resolved), cv2.IMREAD_GRAYSCALE)
    if mask_img is None:
        return False
    if mask_img.ndim == 3:
        mask_img = mask_img[:, :, 0]
    fg = (mask_img > 0).astype(np.uint8)
    total_fg = int(np.count_nonzero(fg))
    if total_fg < 50:
        return False

    h_m, w_m = fg.shape[:2]
    max_edge_touch = max(3, int(0.02 * max(h_m, w_m)))
    if any(int(np.count_nonzero(edge)) > max_edge_touch for edge in (fg[0, :], fg[-1, :], fg[:, 0], fg[:, -1])):
        return False

    num_l, lbls_cc, stats_cc, _ = cv2.connectedComponentsWithStats(fg, connectivity=8)
    if num_l <= 1:
        return False
    min_sig_area = max(25, int(0.01 * total_fg))
    cc_areas = stats_cc[1:, cv2.CC_STAT_AREA]
    sig_comps = int(np.sum(cc_areas >= min_sig_area))
    if sig_comps != 1:
        best_cc = 1 + int(np.argmax(cc_areas))
        best_cc_area = int(stats_cc[best_cc, cv2.CC_STAT_AREA])
        if float(best_cc_area) / float(max(total_fg, 1)) >= 0.85:
            fg = (lbls_cc == best_cc).astype(np.uint8)
            total_fg = best_cc_area
        else:
            # Allow multi-island deeply pinnatifid Packera leaves when clean and unclipped
            if sig_comps > 12:
                return False
            fg = (lbls_cc == best_cc).astype(np.uint8)
            total_fg = int(np.count_nonzero(fg))

    cnts, hier = cv2.findContours(fg, cv2.RETR_CCOMP, cv2.CHAIN_APPROX_NONE)
    if not cnts or hier is None:
        return False
    outer_cnts = [
        cnts[i] for i in range(len(cnts))
        if int(hier[0][i][3]) == -1 and float(cv2.contourArea(cnts[i])) >= min(25, total_fg // 2)
    ]
    if not outer_cnts:
        return False
    main_cnt = max(outer_cnts, key=cv2.contourArea)
    if len(main_cnt) < 12:
        return False
    area_c = float(cv2.contourArea(main_cnt))
    perim_c = float(cv2.arcLength(main_cnt, True))
    if area_c <= 0.0 or perim_c <= 0.0:
        return False
    sol_c = area_c / max(float(cv2.contourArea(cv2.convexHull(main_cnt))), 1.0)
    min_hole_area = max(12.0, (0.025 if sol_c < 0.80 else 0.005) * float(total_fg))
    for i, c in enumerate(cnts):
        if int(hier[0][i][3]) != -1 and float(cv2.contourArea(c)) >= min_hole_area:
            return False
    circ = (4.0 * math.pi * area_c) / (perim_c * perim_c)
    return not (circ < 0.12 and not (circ >= 0.020 and sol_c < 0.78))


def _discover_downstream_table(filename: str, out_dir: Path, active_cats: Set[str]) -> Optional[Path]:
    """Locate downstream QC/morphometric CSV tables in run-scoped or workspace-root table directories."""
    for d in (out_dir.parent.parent / "data" / "tables", out_dir.parent / "data" / "tables", out_dir.parent / "tables", out_dir / "tables"):
        cand = d / filename
        if cand.exists() and cand.stat().st_size > 0:
            return cand

    root_cand = PROJECT_ROOT / "data" / "tables" / filename
    if root_cand.exists() and root_cand.stat().st_size > 0:
        try:
            out_res, root_res = out_dir.resolve(), PROJECT_ROOT.resolve()
            if out_res == root_res or root_res in out_res.parents:
                return root_cand
        except Exception:
            pass
        if active_cats:
            try:
                tdf = pd.read_csv(root_cand, usecols=lambda c: c == "catalogNumber", dtype=str)
                if "catalogNumber" in tdf.columns and (set(tdf["catalogNumber"].astype(str).str.strip()) & active_cats):
                    return root_cand
            except Exception:
                pass
    return None


def _load_mask_or_contour_u8(row: pd.Series, out_dir: Path) -> Tuple[Optional[np.ndarray], Optional[np.ndarray]]:
    """Load upright binary mask (`uint8` 0/255) and optional RGB patch via LM3 `render_leaf_products` rules."""
    resolve_fn = _resolve_facade_hook("_resolve_asset_path", _resolve_asset_path, _ORIG_RESOLVE_ASSET_PATH)
    mask_u8: Optional[np.ndarray] = None
    rgb_patch: Optional[np.ndarray] = None

    m_path = resolve_fn(row.get("mask_path", ""), out_dir)
    if m_path is not None and m_path.exists():
        m_in = cv2.imread(str(m_path), cv2.IMREAD_GRAYSCALE)
        if m_in is not None and np.count_nonzero(m_in) > 0:
            mask_u8 = ((m_in > 0).astype(np.uint8)) * 255

    if mask_u8 is None:
        c_path = resolve_fn(row.get("contour_path", ""), out_dir)
        if c_path is not None and c_path.exists():
            try:
                c_df = pd.read_csv(c_path)
                if {"x", "y"}.issubset(c_df.columns) and len(c_df) >= 6:
                    pts = c_df[["x", "y"]].to_numpy(dtype=np.float32)
                    min_xy = pts.min(axis=0)
                    scale_f = 320.0 / float(np.max(np.maximum(pts.max(axis=0) - min_xy, 1.0)))
                    norm_pts = ((pts - min_xy) * scale_f + 32.0).astype(np.int32)
                    canvas = np.zeros((384, 384), dtype=np.uint8)
                    cv2.fillPoly(canvas, [norm_pts], 255)
                    mask_u8 = canvas
            except Exception:
                pass

    if mask_u8 is None or int(np.count_nonzero(mask_u8)) == 0:
        return None, None

    is_upright = _parse_bool_val(row.get("mask_is_upright", True)) if ("mask_is_upright" in row and pd.notna(row.get("mask_is_upright"))) else True
    rot_cw: Optional[float] = None
    m_ang = row.get("midrib_angle_deg", np.nan)
    if pd.notna(m_ang) and str(m_ang).strip() != "":
        try:
            scw = (-(float(m_ang) + 90.0)) % 360.0
            if 0.25 < scw < 359.75:
                rot_cw = float(scw)
        except (TypeError, ValueError):
            rot_cw = None
    angle_cw: Optional[float] = None if is_upright else rot_cw

    # Attempt to extract authentic RGB crop from voucher sheet if bbox and scan exist
    img_cand = resolve_fn(row.get("image_path", ""), out_dir)
    if img_cand is None:
        cat_s = str(row.get("catalogNumber", "")).strip()
        if cat_s:
            for ext in (".webp", ".jpg", ".jpeg", ".png", ".tif"):
                cand_raw = resolve_fn(f"raw_vouchers/{cat_s}{ext}", out_dir)
                if cand_raw is not None:
                    img_cand = cand_raw
                    break

    b_cols = ("bbox_x1", "bbox_y1", "bbox_x2", "bbox_y2")
    if img_cand is not None and img_cand.exists() and all(c in row and pd.notna(row.get(c)) for c in b_cols):
        try:
            bgr_full = cv2.imread(str(img_cand), cv2.IMREAD_COLOR)
            if bgr_full is not None:
                H_f, W_f = bgr_full.shape[:2]
                x1, y1 = max(0, min(W_f - 1, int(round(float(row["bbox_x1"]))))), max(0, min(H_f - 1, int(round(float(row["bbox_y1"])))))
                x2, y2 = max(x1 + 1, min(W_f, int(round(float(row["bbox_x2"]))))), max(y1 + 1, min(H_f, int(round(float(row["bbox_y2"])))))
                crop_rgb = cv2.cvtColor(bgr_full[y1:y2, x1:x2], cv2.COLOR_BGR2RGB)
                if is_upright and rot_cw is not None:
                    from leafmachine3.core.imaging import rotate_image
                    crop_rgb = rotate_image(crop_rgb, rot_cw, bg=255, nearest=False)
                mh, mw = mask_u8.shape[:2]
                ch_r, cw_r = crop_rgb.shape[:2]
                if (ch_r, cw_r) != (mh, mw):
                    if is_upright and mh <= ch_r + 20 and mw <= cw_r + 20 and ch_r <= max(int(2.2 * mh), mh + 48) and cw_r <= max(int(2.2 * mw), mw + 48):
                        aligned = np.full((mh, mw, 3), 255, dtype=np.uint8)
                        sy0, sx0 = max(0, (ch_r - mh) // 2), max(0, (cw_r - mw) // 2)
                        dy0, dx0 = max(0, (mh - ch_r) // 2), max(0, (mw - cw_r) // 2)
                        cp_h, cp_w = min(ch_r - sy0, mh - dy0), min(cw_r - sx0, mw - dx0)
                        if cp_h > 0 and cp_w > 0:
                            aligned[dy0 : dy0 + cp_h, dx0 : dx0 + cp_w] = crop_rgb[sy0 : sy0 + cp_h, sx0 : sx0 + cp_w]
                        crop_rgb = aligned
                    else:
                        crop_rgb = cv2.resize(crop_rgb, (mw, mh), interpolation=cv2.INTER_AREA)
                rgb_patch = crop_rgb
        except Exception:
            rgb_patch = None

    # Delegate orientation & cutout generation directly to LeafMachine3 render_leaf_products
    from leafmachine3.reporting.leaf_products import render_leaf_products

    base_rgb = rgb_patch if rgb_patch is not None else np.full((*mask_u8.shape[:2], 3), 42, dtype=np.uint8)
    prods = render_leaf_products(
        base_rgb, {"lamina": mask_u8 > 0, "silhouette": mask_u8 > 0, "silhouettepet": mask_u8 > 0, "hole": None},
        angle_cw=angle_cw, want={"lamina", "rgb_cropped"} if rgb_patch is not None else {"lamina"},
    )
    sil_prod = prods.get("lamina") if prods.get("lamina") is not None else prods.get("silhouette")
    if sil_prod is not None and np.count_nonzero(sil_prod) > 0:
        mask_u8 = (np.asarray(sil_prod) > 0).astype(np.uint8) * 255
    if rgb_patch is not None:
        rgb_prod = prods.get("rgb_cropped")
        if isinstance(rgb_prod, np.ndarray) and rgb_prod.ndim == 3 and rgb_prod.shape[:2] == mask_u8.shape[:2]:
            rgb_patch = rgb_prod
        elif angle_cw is not None:
            from leafmachine3.core.imaging import rotate_image
            rgb_patch = rotate_image(rgb_patch, angle_cw, bg=255, nearest=False)

    # Tightly crop around foreground with clean margin
    ys, xs = np.where(mask_u8 > 0)
    if len(xs) > 0 and len(ys) > 0:
        pad = max(10, int(0.08 * max(ys.max() - ys.min() + 1, xs.max() - xs.min() + 1)))
        y0, y1 = max(0, int(ys.min()) - pad), min(mask_u8.shape[0], int(ys.max()) + pad + 1)
        x0, x1 = max(0, int(xs.min()) - pad), min(mask_u8.shape[1], int(xs.max()) + pad + 1)
        mask_u8 = mask_u8[y0:y1, x0:x1]
        if rgb_patch is not None and rgb_patch.shape[:2] != mask_u8.shape[:2]:
            rgb_patch = rgb_patch[y0:y1, x0:x1]
    return mask_u8, rgb_patch


__all__ = [
    "PROJECT_ROOT", "TOL_OKABE_12", "HOLLOW_MARKERS_5", "ICN_ROMAN_TOKENS",
    "_ASSET_ANCHOR_DIRS", "_resolve_facade_hook", "_atomic_save_matplotlib_fig",
    "_standardize_packera_species_py", "format_icn_taxon_mathtext", "get_combinatorial_taxon_style",
    "_parse_bool_val", "_resolve_asset_path", "_is_valid_smooth_closed_mask",
    "_discover_downstream_table", "_load_mask_or_contour_u8",
]
