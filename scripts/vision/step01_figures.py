#!/usr/bin/env python3
"""
Step 01 publication-grade cohort biogeography, determiner authority, EPA Level III
ecoregion, and temporal/optical QC audit figure (`fig01_cohort_biogeography_and_authority.png`).
"""

from __future__ import annotations

import math, os, sqlite3, sys
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple, Union

import numpy as np
import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from scripts.vision.step02_fig_style_assets import (
    _atomic_save_matplotlib_fig, _parse_bool_val, _standardize_packera_species_py,
    format_icn_taxon_mathtext, get_combinatorial_taxon_style,
)

DEFAULT_FIG01_PNG = "outputs/figures/fig01_cohort_biogeography_and_authority.png"

_NA_LAND_LONLAT: Tuple[Tuple[float, float], ...] = (
    (-130.0, 54.0), (-125.0, 49.0), (-124.5, 42.0), (-122.5, 37.8), (-117.2, 32.5), (-114.5, 31.8), (-110.5, 24.0), (-105.0, 20.5),
    (-97.5, 22.0), (-97.2, 26.0), (-94.0, 29.7), (-89.2, 29.1), (-84.0, 30.0), (-82.5, 27.5), (-80.2, 25.2), (-80.0, 27.0),
    (-81.4, 30.8), (-78.0, 33.9), (-75.5, 35.2), (-75.0, 38.8), (-74.0, 40.5), (-70.0, 41.8), (-67.0, 44.8), (-64.0, 44.6),
    (-60.0, 46.5), (-64.5, 48.8), (-67.5, 49.2), (-61.0, 50.2), (-55.5, 52.0), (-60.5, 56.5), (-65.0, 60.0), (-78.0, 62.0),
    (-79.5, 52.0), (-82.5, 55.2), (-94.5, 59.0), (-95.0, 64.0), (-115.0, 66.0), (-132.0, 62.0), (-130.0, 54.0),
)

# Ellipsoidal GRS80/NAD83 EPSG:5070 Albers Equal-Area Conic constants (km)
_A_KM, _F = 6378.137, 1.0 / 298.257222101
_E2 = 2.0 * _F - _F * _F
_E = math.sqrt(_E2)
_PHI1, _PHI2, _PHI0, _LAM0 = map(math.radians, (29.5, 45.5, 23.0, -96.0))


def _albers_q(phi: np.ndarray) -> np.ndarray:
    s = np.sin(phi)
    return (1.0 - _E2) * (s / (1.0 - _E2 * s * s) - (0.5 / _E) * np.log((1.0 - _E * s) / (1.0 + _E * s)))


_Q1, _Q2, _Q0 = (float(v) for v in _albers_q(np.array([_PHI1, _PHI2, _PHI0])))
_M1, _M2 = (math.cos(p) / math.sqrt(1.0 - _E2 * math.sin(p) ** 2) for p in (_PHI1, _PHI2))
_ALB_N = (_M1 * _M1 - _M2 * _M2) / (_Q2 - _Q1)
_ALB_C = _M1 * _M1 + _ALB_N * _Q1
_ALB_RHO0 = (_A_KM / _ALB_N) * math.sqrt(max(_ALB_C - _ALB_N * _Q0, 0.0))


def _project_epsg5070_km(lon_deg: Any, lat_deg: Any) -> Tuple[np.ndarray, np.ndarray]:
    """Project (lon, lat) degrees to EPSG:5070 Albers Equal-Area Conic coordinates (km)."""
    lon, lat = np.radians(np.asarray(lon_deg, dtype=float)), np.radians(np.clip(np.asarray(lat_deg, dtype=float), -85.0, 85.0))
    rho = (_A_KM / _ALB_N) * np.sqrt(np.maximum(_ALB_C - _ALB_N * _albers_q(lat), 0.0))
    theta = _ALB_N * (lon - _LAM0)
    return rho * np.sin(theta), _ALB_RHO0 - rho * np.cos(theta)


def _load_vouchers_input(vouchers_input: Union[pd.DataFrame, str, Path]) -> Optional[pd.DataFrame]:
    """Load Step 01 vouchers from a DataFrame, CSV file, or SQLite database."""
    if isinstance(vouchers_input, pd.DataFrame):
        return vouchers_input.copy()
    p = Path(vouchers_input)
    if not p.exists() or p.stat().st_size == 0:
        return None
    if p.suffix.lower() in {".sqlite", ".db", ".sqlite3"}:
        conn = sqlite3.connect(str(p))
        try:
            return pd.read_sql_query("SELECT * FROM vouchers;", conn)
        except Exception:
            return None
        finally:
            conn.close()
    return pd.read_csv(p, low_memory=False)


def _resolve_fig01_output_path(vouchers_input: Union[pd.DataFrame, str, Path], output_png: Union[str, Path]) -> Path:
    """Auto-resolve run-scoped fig01 path when vouchers live under <run_dir>/data/tables/."""
    if str(output_png) == DEFAULT_FIG01_PNG and not isinstance(vouchers_input, pd.DataFrame):
        vp = Path(vouchers_input)
        if len(vp.parents) >= 3 and vp.parent.name == "tables" and vp.parents[1].name == "data":
            return vp.parents[2] / "outputs" / "figures" / Path(DEFAULT_FIG01_PNG).name
    return Path(output_png)


def _prepare_cohort_flags(df: pd.DataFrame) -> pd.DataFrame:
    """Attach standardized species, centroid quarantine, type, and 5-tier authority columns."""
    sp_col = next((c for c in ("species_standardized", "scientificName", "species_raw", "taxon") if c in df.columns), None)
    inf_col = "infraspecific_taxon" if "infraspecific_taxon" in df.columns else None
    df["_species"] = (
        [_standardize_packera_species_py(s, inf) for s, inf in zip(df[sp_col], df[inf_col])]
        if (sp_col is not None and inf_col is not None)
        else (df[sp_col].apply(_standardize_packera_species_py) if sp_col is not None else "Unassigned_Packera")
    )
    c_mask = df["is_centroid_or_inst_coord"].apply(_parse_bool_val) if "is_centroid_or_inst_coord" in df.columns else pd.Series(False, index=df.index, dtype=bool)
    if "coord_quarantine_reason" in df.columns:
        c_mask |= ~df["coord_quarantine_reason"].astype(str).str.strip().str.lower().isin({"", "none", "nan", "na", "false", "0"})
    if "georef_qc_status" in df.columns:
        c_mask |= df["georef_qc_status"].astype(str).str.strip().str.lower().str.contains("centroid|inst|uncertain_coord|quarantine", regex=True)
    df["_is_centroid"] = c_mask

    t_mask = df["is_type_specimen"].apply(_parse_bool_val) if "is_type_specimen" in df.columns else pd.Series(False, index=df.index, dtype=bool)
    if "type_status" in df.columns:
        t_mask |= ~df["type_status"].astype(str).str.strip().str.lower().isin({"", "none", "nan", "na", "not a type", "false", "0"})
    df["_is_type"] = t_mask

    def _tier_of(row: pd.Series, is_t: bool) -> str:
        t_raw = str(row.get("determiner_tier", "")).strip()
        if is_t:
            return "Tier_0_Type_Anchor"
        if _parse_bool_val(row.get("is_open_or_hybrid_concept", False)) or "passive_open" in t_raw.lower():
            return "Passive_Open_or_Hybrid_Concept"
        return "Tier_1_Gold" if ("Tier_1" in t_raw or "Gold" in t_raw) else ("Tier_2_Silver" if ("Tier_2" in t_raw or "Silver" in t_raw) else "Tier_3_Bronze")

    df["_auth_tier"] = [_tier_of(row, bool(is_t)) for (_, row), is_t in zip(df.iterrows(), t_mask)]
    return df


def _draw_panel_a_map(ax: Any, df: pd.DataFrame, active_taxa: List[str], style_map: Dict[str, Tuple[str, str]]) -> None:
    """Panel A: EPSG:5070 Albers Equal-Area map with centroid quarantine and type anchors."""
    from matplotlib.lines import Line2D

    ax.set_facecolor("#EEF4F8")
    lat_col = next((c for c in ("decimalLatitude", "latitude", "lat") if c in df.columns), None)
    lon_col = next((c for c in ("decimalLongitude", "longitude", "lon") if c in df.columns), None)
    n_quar = 0
    plot_df = df.iloc[0:0].copy()
    if lat_col is not None and lon_col is not None:
        lat, lon = pd.to_numeric(df[lat_col], errors="coerce"), pd.to_numeric(df[lon_col], errors="coerce")
        valid = lat.notna() & lon.notna() & np.isfinite(lat) & np.isfinite(lon) & (lat.abs() <= 85.0) & (lon.abs() <= 180.0)
        sub = df.loc[valid].assign(_lat=lat.loc[valid].to_numpy(dtype=float), _lon=lon.loc[valid].to_numpy(dtype=float))
        n_quar = int(((~sub["_is_type"]) & sub["_is_centroid"]).sum())
        plot_df = sub[~((~sub["_is_type"]) & sub["_is_centroid"])].copy()
    if len(plot_df) == 0:
        msg = "No geographic coordinates\nin cohort" if (lat_col is None or lon_col is None) else f"No clean georeferenced records\n({n_quar} centroid/inst. quarantined)"
        ax.text(0.5, 0.5, msg, ha="center", va="center", transform=ax.transAxes, fontsize=8.5, color="#64748B")
        ax.set_title("A. Cohort Biogeography (EPSG:5070)", fontsize=10, fontweight="bold", loc="left")
        return

    px, py = _project_epsg5070_km(plot_df["_lon"], plot_df["_lat"])
    plot_df["_x"], plot_df["_y"] = px, py
    ax.fill(*_project_epsg5070_km(*zip(*_NA_LAND_LONLAT)), facecolor="#FAFAF8", edgecolor="#94A3B8", lw=0.8, zorder=1)
    glons, glats = np.linspace(-135.0, -52.0, 80), np.linspace(18.0, 68.0, 80)
    for lat_g in range(20, 70, 5):
        ax.plot(*_project_epsg5070_km(glons, np.full_like(glons, float(lat_g))), color="#CBD5E1", lw=0.45, ls=":", zorder=2)
    for lon_g in range(-130, -50, 10):
        ax.plot(*_project_epsg5070_km(np.full_like(glats, float(lon_g)), glats), color="#CBD5E1", lw=0.45, ls=":", zorder=2)

    non_t = plot_df[~plot_df["_is_type"]]
    alpha_pt = 0.85 if len(non_t) <= 200 else (0.68 if len(non_t) <= 1000 else 0.52)
    for tx in active_taxa:
        tr = non_t[non_t["_species"] == tx]
        if len(tr) > 0:
            col, mk = style_map[tx]
            ax.scatter(tr["_x"], tr["_y"], s=28, marker=mk, facecolors="none", edgecolors=col, linewidths=0.95, alpha=alpha_pt, zorder=4)

    flags = {"conf": False, "syn": False, "unc": False}
    for _, tr in plot_df[plot_df["_is_type"]].iterrows():
        s_blob = " ".join(str(tr.get(c, "")) for c in ("type_control_status", "concept_synonymy_status", "taxonomic_treatment_subtype", "determiner_tier")).lower()
        if bool(tr["_is_centroid"]) or "passive_type_uncertain_coord" in s_blob:
            flags["unc"] = True
            ax.scatter([tr["_x"]], [tr["_y"]], s=76, marker="o", facecolors="none", edgecolors="#64748B", linewidths=1.3, zorder=6)
            ax.scatter([tr["_x"]], [tr["_y"]], s=36, marker="D", facecolors="none", edgecolors="#64748B", linewidths=1.1, zorder=7)
        elif "synonymy" in s_blob or "heterotypic" in s_blob:
            flags["syn"] = True
            ax.scatter([tr["_x"]], [tr["_y"]], s=76, marker="^", facecolors="none", edgecolors="#882255", linewidths=1.5, zorder=7)
        else:
            flags["conf"] = True
            ax.scatter([tr["_x"]], [tr["_y"]], s=86, marker="D", facecolors="none", edgecolors="#332288", linewidths=1.5, zorder=6)
            ax.scatter([tr["_x"]], [tr["_y"]], s=30, marker="D", facecolors="none", edgecolors="#111827", linewidths=1.0, zorder=7)

    x_lo, x_hi, y_lo, y_hi = float(np.min(px)), float(np.max(px)), float(np.min(py)), float(np.max(py))
    x_pad, y_pad = max(180.0, 0.14 * max(x_hi - x_lo, 200.0)), max(160.0, 0.14 * max(y_hi - y_lo, 200.0))
    ax.set(xlim=(x_lo - x_pad, x_hi + x_pad), ylim=(y_lo - y_pad, y_hi + y_pad), xlabel="Albers Easting (km, EPSG:5070)", ylabel="Albers Northing (km, EPSG:5070)")
    ax.set_aspect("equal", adjustable="box")

    (cx0, cx1), (cy0, cy1) = ax.get_xlim(), ax.get_ylim()
    sx, sy = cx1 - cx0, cy1 - cy0
    bar_km = min((50, 100, 200, 250, 500, 1000, 2000), key=lambda v: abs(v - 0.22 * sx))
    bx0, by0 = cx0 + 0.05 * sx, cy0 + 0.07 * sy
    ax.plot([bx0, bx0 + bar_km], [by0, by0], color="#111827", lw=2.0, zorder=9)
    ax.text(bx0 + 0.5 * bar_km, by0 + 0.025 * sy, f"{bar_km} km", ha="center", va="bottom", fontsize=7.2, color="#111827", fontweight="bold", zorder=9)
    ax.annotate("N", xy=(0.94, 0.94), xytext=(0.94, 0.83), xycoords="axes fraction", ha="center", va="center", fontsize=8.5, fontweight="bold", color="#111827", arrowprops=dict(arrowstyle="-|>", color="#111827", lw=1.2), zorder=9)

    leg_specs = [("conf", "D", "#332288", "Type Anchor (Confirmed/Homotypic)"), ("syn", "^", "#882255", "Type Synonymy Candidate"), ("unc", "o", "#64748B", "Type (Uncertain/Centroid Coord)")]
    handles = [Line2D([0], [0], marker=m, color="none", markerfacecolor="none", markeredgecolor=c, markeredgewidth=1.4, markersize=6.5, label=lb) for k, m, c, lb in leg_specs if flags[k]]
    if handles:
        ax.legend(handles=handles, loc="lower right", fontsize=6.7, frameon=True, facecolor="white", framealpha=0.90, edgecolor="#CBD5E1")
    ax.set_title(f"A. Cohort Biogeography (EPSG:5070)\n(n={len(plot_df)} mapped | {n_quar} non-type centroid/inst. quarantined)", fontsize=9.3, fontweight="bold", loc="left")


def _draw_panel_b_authority(ax: Any, df: pd.DataFrame) -> None:
    """Panel B: Taxon-stratified stacked horizontal bars of authority and concept tiers."""
    sp_counts = df["_species"].value_counts()
    top_taxa = list(sp_counts.head(18).index)
    tiers = [("Tier_0_Type_Anchor", "Type Anchor", "#332288"), ("Tier_1_Gold", "Tier 1 Gold", "#E69F00"), ("Tier_2_Silver", "Tier 2 Silver", "#56B4E9"), ("Tier_3_Bronze", "Tier 3 Bronze", "#94A3B8"), ("Passive_Open_or_Hybrid_Concept", "Passive Open/Hybrid", "#CC79A7")]
    y_pos, left_acc = np.arange(len(top_taxa)), np.zeros(len(top_taxa), dtype=float)
    for t_key, t_lbl, t_col in tiers:
        vals = np.array([int(((df["_species"] == tx) & (df["_auth_tier"] == t_key)).sum()) for tx in top_taxa], dtype=float)
        if vals.sum() > 0:
            ax.barh(y_pos, vals, left=left_acc, height=0.64, color=t_col, alpha=0.90, label=t_lbl)
            left_acc += vals
    ax.set(yticks=y_pos, xlabel="Curated Voucher Count (n)")
    ax.set_yticklabels([f"{format_icn_taxon_mathtext(tx, abbreviate_genus=(len(top_taxa) > 5))} (n={int(sp_counts[tx])})" for tx in top_taxa], fontsize=7.6)
    ax.invert_yaxis()
    ax.legend(loc="lower center", bbox_to_anchor=(0.5, 1.01), ncol=3, fontsize=6.7, frameon=False)
    ax.grid(True, axis="x", linestyle=":", lw=0.5, alpha=0.55)
    ax.set_title(f"B. Taxon Authority & Concept Tiers (N={len(df)})", fontsize=9.6, fontweight="bold", loc="left", pad=32)


def _draw_panel_c_ecoregions(ax: Any, df: pd.DataFrame) -> None:
    """Panel C: Top 12 US EPA Level III Ecoregion (regional_group) voucher counts."""
    col = "regional_group" if "regional_group" in df.columns else ("stateProvince" if "stateProvince" in df.columns else None)
    vals = df[col].astype(str).str.strip() if col is not None else pd.Series(dtype=str)
    vals = vals[~vals.str.lower().isin({"", "nan", "none", "na", "unknown"})]
    if col is None or len(vals) == 0:
        ax.text(0.5, 0.5, "No populated EPA Level III\necoregions in cohort", ha="center", va="center", transform=ax.transAxes, fontsize=8.5, color="#64748B")
        ax.set_title("C. Top US EPA Level III Ecoregions", fontsize=9.6, fontweight="bold", loc="left")
        return
    top_eco = vals.value_counts().head(12)
    y_pos, max_v = np.arange(len(top_eco)), max(1.0, float(top_eco.max()))
    ax.barh(y_pos, top_eco.to_numpy(dtype=float), height=0.62, color="#44AA99", alpha=0.88)
    for idx_e, cnt in enumerate(top_eco.to_numpy()):
        ax.text(float(cnt) + 0.015 * max_v, idx_e, str(int(cnt)), va="center", ha="left", fontsize=7.4, color="#1F2937")
    ax.set(yticks=y_pos, xlim=(0, max_v * 1.14), xlabel="Voucher Count (n)")
    ax.set_yticklabels([lb[:34] + ("..." if len(lb) > 34 else "") for lb in top_eco.index], fontsize=7.6)
    ax.invert_yaxis()
    ax.grid(True, axis="x", linestyle=":", lw=0.5, alpha=0.55)
    ax.set_title(f"C. Top {len(top_eco)} {'US EPA Level III Ecoregions' if col == 'regional_group' else 'States / Provinces'}", fontsize=9.6, fontweight="bold", loc="left")


def _draw_panel_d_temporal_qc(ax: Any, df: pd.DataFrame) -> None:
    """Panel D: Collection-year and/or P90 Laplacian sharpness optical QC distribution."""
    yr_s = pd.to_numeric(df["year"], errors="coerce") if "year" in df.columns else pd.Series(np.nan, index=df.index)
    if yr_s.isna().all() and "eventDate" in df.columns:
        yr_s = pd.to_numeric(df["eventDate"].astype(str).str.extract(r"(\d{4})", expand=False), errors="coerce")
    yr_s = yr_s.where((yr_s >= 1750) & (yr_s <= 2035))
    sh_col = next((c for c in ("sharpness_score", "laplacian_sharpness", "p90_sharpness") if c in df.columns), None)
    sh_s = pd.to_numeric(df[sh_col], errors="coerce") if sh_col is not None else pd.Series(np.nan, index=df.index)
    both = yr_s.notna() & sh_s.notna() & np.isfinite(yr_s) & np.isfinite(sh_s)

    if int(both.sum()) >= 3:
        t_pal = {"Tier_0_Type_Anchor": "#332288", "Tier_1_Gold": "#E69F00", "Tier_2_Silver": "#56B4E9", "Tier_3_Bronze": "#64748B", "Passive_Open_or_Hybrid_Concept": "#CC79A7"}
        for tk, col in t_pal.items():
            m = both & (df["_auth_tier"] == tk)
            if m.any():
                ax.scatter(yr_s[m], sh_s[m], s=28, marker="o", facecolors="none", edgecolors=col, linewidths=0.95, alpha=0.80)
        ax.axvline(float(yr_s[both].median()), color="#475569", ls="--", lw=0.8)
        ax.axhline(float(sh_s[both].median()), color="#475569", ls=":", lw=0.8)
        ax.set(xlabel="Collection Year", ylabel="P90 Laplacian Sharpness")
        ax.set_title(f"D. Temporal & Optical Sharpness QC (med yr={int(round(float(yr_s[both].median())))}, P90={float(sh_s[both].median()):.1f})", fontsize=9.1, fontweight="bold", loc="left")
    elif int(sh_s.notna().sum()) >= 2:
        vals = sh_s.dropna().to_numpy(dtype=float)
        bp = ax.boxplot([vals], orientation="horizontal", widths=0.45, patch_artist=True, showfliers=False, medianprops=dict(color="#111827", lw=1.4))
        for b in bp["boxes"]:
            b.set(facecolor="none", edgecolor="#475569")
        ax.scatter(vals, np.ones(len(vals)) + np.random.default_rng(42).uniform(-0.15, 0.15, len(vals)), s=28, marker="o", facecolors="none", edgecolors="#0072B2", linewidths=0.95, alpha=0.8)
        ax.set(yticks=[1], yticklabels=[f"Cohort\n(n={len(vals)})"], xlabel="P90 Laplacian Sharpness")
        ax.set_title(f"D. Optical Sharpness QC Distribution (med={float(np.median(vals)):.1f})", fontsize=9.4, fontweight="bold", loc="left")
    elif int(yr_s.notna().sum()) >= 2:
        vals = yr_s.dropna().to_numpy(dtype=float)
        ax.hist(vals, bins=min(24, max(6, int(np.sqrt(len(vals))))), color="#0072B2", alpha=0.78, edgecolor="white")
        ax.axvline(float(np.median(vals)), color="#D55E00", ls="--", lw=1.2, label=f"Median ({int(round(float(np.median(vals))))})")
        ax.legend(loc="upper left", fontsize=7.4, frameon=True, facecolor="white")
        ax.set(xlabel="Collection Year", ylabel="Voucher Count (n)")
        ax.set_title(f"D. Historical Collection Year Distribution (n={len(vals)})", fontsize=9.4, fontweight="bold", loc="left")
    else:
        ax.text(0.5, 0.5, "No collection-year or P90 sharpness\nrecords available in cohort", ha="center", va="center", transform=ax.transAxes, fontsize=8.5, color="#64748B")
        ax.set(xticks=[], yticks=[])
        ax.set_title("D. Collection Year & Optical Sharpness QC", fontsize=9.4, fontweight="bold", loc="left")
        return
    ax.grid(True, linestyle=":", lw=0.5, alpha=0.55)


def generate_step01_cohort_audit_figure(
    vouchers_input: Union[pd.DataFrame, str, Path], output_png: Union[str, Path] = DEFAULT_FIG01_PNG, dpi: int = 300,
) -> Optional[Path]:
    """Render 300-DPI 2x2 Step 01 cohort biogeography, authority, ecoregion, and optical QC figure."""
    if os.environ.get("PLANTMATICS_SKIP_FIGURES") == "1" or os.environ.get("PACKERA_SKIP_FIGURES") == "1":
        return None
    df = _load_vouchers_input(vouchers_input)
    if df is None or len(df) == 0:
        return None
    import matplotlib
    if "matplotlib.pyplot" not in sys.modules:
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    df = _prepare_cohort_flags(df)
    active_taxa = sorted(t for t in df["_species"].unique() if str(t).strip())
    fig = plt.figure(figsize=(14.2, 10.0), dpi=int(dpi), layout="constrained")
    gs = fig.add_gridspec(2, 2)
    _draw_panel_a_map(fig.add_subplot(gs[0, 0]), df, active_taxa, get_combinatorial_taxon_style(active_taxa))
    _draw_panel_b_authority(fig.add_subplot(gs[0, 1]), df)
    _draw_panel_c_ecoregions(fig.add_subplot(gs[1, 0]), df)
    _draw_panel_d_temporal_qc(fig.add_subplot(gs[1, 1]), df)
    fig.suptitle(f"Step 01: Curated Herbarium Cohort Biogeography, Determiner Authority & QC Audit (N={len(df)} vouchers, {len(active_taxa)} taxa)", fontsize=11.2, fontweight="bold")
    saved = _atomic_save_matplotlib_fig(fig, _resolve_fig01_output_path(vouchers_input, output_png), dpi=int(dpi))
    plt.close(fig)
    return saved


__all__ = ["DEFAULT_FIG01_PNG", "generate_step01_cohort_audit_figure", "_project_epsg5070_km"]
