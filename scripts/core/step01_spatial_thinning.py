"""step01_spatial_thinning.py - Spatial balancing, quota thinning, and inverse geodesic density weights."""

from __future__ import annotations

import math
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Set, Tuple

try:
    import numpy as np
except ImportError:
    np = None  # type: ignore[assignment]
try:
    import pandas as pd
except ImportError:
    pd = None  # type: ignore[assignment]

from scripts.core.step01_determiner_surnames import classify_nomenclatural_type_rank, classify_nomenclatural_type_subrank
from scripts.core.step01_duplicate_stratifier import stratify_duplicates
from scripts.core.step01_exsiccatae_fingerprint import (
    _as_bool_flag, get_image_resolution_score, haversine_distance_meters, resolve_voucher_image_path,
)


def spatially_balance_and_thin_vouchers(
    df: pd.DataFrame,
    spatial_thin_km: float = 0.0,
    max_per_state: Optional[int] = None,
    max_per_region: Optional[int] = None,
    primary_duplicates_only: bool = True,
    preserve_all_gold: bool = True,
    manifest_df: Optional[pd.DataFrame] = None,
    workspace_dir: Optional[Path] = None,
    preserve_catalog_numbers: Optional[Iterable[str]] = None,
) -> pd.DataFrame:
    """
    Stratifies, deduplicates, spatially thins, and computes continuous inverse
    spatial-density weights (`spatial_balance_weight`) for botanical vouchers.
    """
    if df is None:
        from scripts.core.taxonomy_authority import EXPORT_COLUMNS
        return pd.DataFrame(columns=EXPORT_COLUMNS + ["spatial_balance_weight"])
    if df.empty:
        out_empty = df.copy()
        out_empty["spatial_balance_weight"] = pd.Series(dtype=float)
        return out_empty

    from scripts.core.harvester import standardize_packera_taxon

    out_df = df.copy()

    # Ensure standardized taxon column exists
    if "species_standardized" not in out_df.columns:
        src_col = "scientificName" if "scientificName" in out_df.columns else (
            "species_raw" if "species_raw" in out_df.columns else None
        )
        if src_col is not None:
            out_df["species_standardized"] = out_df[src_col].map(
                lambda x: standardize_packera_taxon(str(x)) or str(x).strip() or "Unknown"
                if not pd.isna(x) else "Unknown"
            )
        else:
            out_df["species_standardized"] = "Unknown"
    else:
        src_col = "scientificName" if "scientificName" in out_df.columns else (
            "species_raw" if "species_raw" in out_df.columns else None
        )
        def _fill_sp(row: pd.Series) -> str:
            val = row.get("species_standardized")
            if val is not None and not pd.isna(val) and str(val).strip():
                return str(val).strip()
            if src_col is not None:
                raw_v = row.get(src_col)
                if raw_v is not None and not pd.isna(raw_v) and str(raw_v).strip():
                    return standardize_packera_taxon(str(raw_v)) or str(raw_v).strip()
            return "Unknown"
        out_df["species_standardized"] = out_df.apply(_fill_sp, axis=1)

    preserve_cats_set: Set[str] = {
        str(c).strip() for c in (preserve_catalog_numbers or [])
        if c is not None and str(c).strip()
    }

    # Step A: Exsiccatae Deduplication (unconditionally retaining all nomenclatural type specimens and preserved vouchers)
    if primary_duplicates_only:
        if "is_primary_duplicate" not in out_df.columns:
            if any(c in out_df.columns for c in ("recordedBy", "recordNumber", "exsiccatae_key")):
                out_df = stratify_duplicates(out_df, workspace_dir=workspace_dir)
            else:
                out_df["is_primary_duplicate"] = True
        prim_mask = out_df["is_primary_duplicate"].map(
            lambda v: v is True or str(v).strip().lower() in {"true", "1", "t", "yes"}
        )
        pres_mask_a = (
            out_df["catalogNumber"].astype(str).str.strip().isin(preserve_cats_set)
            if (preserve_cats_set and "catalogNumber" in out_df.columns)
            else pd.Series(False, index=out_df.index)
        )
        ts_col_a = "type_status" if "type_status" in out_df.columns else (
            "typeStatus" if "typeStatus" in out_df.columns else None
        )
        if ts_col_a is not None:
            type_mask_a = out_df[ts_col_a].map(lambda t: classify_nomenclatural_type_rank(t)[0] >= 1)
            keep_mask = prim_mask | type_mask_a | pres_mask_a
        else:
            keep_mask = prim_mask | pres_mask_a
        out_df = out_df.loc[keep_mask].reset_index(drop=True)
        if out_df.empty:
            out_df["spatial_balance_weight"] = pd.Series(dtype=float)
            return out_df

    # Step B: Candidate Priority Ranking
    t1_leaf_counts: Dict[str, int] = {}
    if manifest_df is not None and not manifest_df.empty and "catalogNumber" in manifest_df.columns:
        for cat_k, cnt in manifest_df["catalogNumber"].astype(str).value_counts().items():
            t1_leaf_counts[str(cat_k)] = int(cnt)

    lat_col = "decimalLatitude" if "decimalLatitude" in out_df.columns else None
    lon_col = "decimalLongitude" if "decimalLongitude" in out_df.columns else None

    n_rows = len(out_df)
    lats: List[Optional[float]] = [None] * n_rows
    lons: List[Optional[float]] = [None] * n_rows
    priority_tuples: List[Tuple[int, float, int, int, float, float, str]] = []
    is_gold_flags: List[bool] = []
    is_type_flags: List[bool] = []
    is_preserved_flags: List[bool] = []

    for i in range(n_rows):
        row = out_df.iloc[i]
        cat_str = str(row.get("catalogNumber", i))
        is_preserved_flags.append(cat_str.strip() in preserve_cats_set)

        # 1. Nomenclatural type rank (2 = Primary Type, 1 = Secondary Type, 0 = Non-type) + ICN Art. 9.6 sub-rank
        raw_ts = row.get("type_status")
        if raw_ts is None or (pd is not None and pd.isna(raw_ts)) or not str(raw_ts).strip():
            raw_ts = row.get("typeStatus")
        has_type, _ = classify_nomenclatural_type_rank(raw_ts)
        t_subrank = classify_nomenclatural_type_subrank(raw_ts)
        is_type_flags.append(has_type >= 1)

        # 2. determiner_tier
        tier_str = str(row.get("determiner_tier", "") if not pd.isna(row.get("determiner_tier", "")) else "").strip()
        tier_score = 3 if (tier_str == "Tier_1_Gold" or has_type >= 1) else (2 if tier_str == "Tier_2_Silver" else 1)
        is_gold_flags.append(tier_str == "Tier_1_Gold" or has_type >= 1)

        # 3. Extracted leaf count
        if cat_str in t1_leaf_counts:
            n_t1 = t1_leaf_counts[cat_str]
        elif "leaf_count" in out_df.columns and not pd.isna(row.get("leaf_count")):
            try:
                n_t1 = int(float(row.get("leaf_count", 0)))
            except (ValueError, TypeError):
                n_t1 = 0
        else:
            n_t1 = 0

        # Coordinates & 4. coordinate_precision
        lat_v: Optional[float] = None
        lon_v: Optional[float] = None
        if lat_col and lon_col:
            try:
                lv = float(row.get(lat_col))
                lnv = float(row.get(lon_col))
                if math.isfinite(lv) and math.isfinite(lnv):
                    lat_v, lon_v = lv, lnv
            except (ValueError, TypeError):
                pass
        lats[i] = lat_v
        lons[i] = lon_v

        coord_prec = 0.0
        if lat_v is not None and lon_v is not None:
            unc_raw = row.get("coordinateUncertainty")
            try:
                unc_f = float(unc_raw)
                if math.isfinite(unc_f) and unc_f >= 0.0:
                    coord_prec = 1.0 / (1.0 + unc_f)
                else:
                    coord_prec = 1.0 / (1.0 + 1000.0)
            except (ValueError, TypeError):
                coord_prec = 1.0 / (1.0 + 1000.0)

        # 5. image_resolution_score
        img_p = str(row.get("image_path", "") if not pd.isna(row.get("image_path", "")) else "")
        resolved_img_p = resolve_voucher_image_path(img_p, workspace_dir=workspace_dir) if img_p else Path("")
        res_score = get_image_resolution_score(str(resolved_img_p) if img_p else "", workspace_dir=workspace_dir)

        priority_tuples.append((has_type, t_subrank, tier_score, n_t1, round(coord_prec, 8), round(res_score, 4), cat_str))

    active_indices = list(range(n_rows))

    # Step C: Greedy Priority-Ordered Pairwise Haversine Thinning (spatial_thin_km > 0)
    if spatial_thin_km is not None and float(spatial_thin_km) > 0.0:
        thin_m = float(spatial_thin_km) * 1000.0
        sp_geo_groups: Dict[str, List[int]] = {}
        non_geo_indices: List[int] = []
        has_hab_cols = any(
            c in out_df.columns
            for c in ("_thin_substrate_stratum", "is_rock_outcrop", "lithic_outcrop_index", "ssurgo_rock_outcrop_pct", "ecotype_regime")
        )

        def _row_hab_stratum(r_idx: int) -> str:
            if not has_hab_cols:
                return ""
            r_sub = out_df.iloc[r_idx]
            raw_s = str(r_sub.get("_thin_substrate_stratum", "") or "").strip()
            if raw_s and raw_s.lower() not in ("nan", "none", "na"):
                return raw_s
            er_s = str(r_sub.get("ecotype_regime", "") or "").strip()
            if er_s and er_s.startswith("Ecotype_"):
                return er_s
            ro_v = r_sub.get("is_rock_outcrop")
            if ro_v is True or str(ro_v).strip().lower() in ("true", "1", "t", "yes"):
                return "azonal_outcrop"
            for nc_h, thr_h in (("lithic_outcrop_index", 0.45), ("ssurgo_rock_outcrop_pct", 15.0)):
                if nc_h in out_df.columns:
                    try:
                        fv = float(r_sub.get(nc_h))
                        if math.isfinite(fv):
                            return "azonal_outcrop" if fv >= thr_h else "zonal_matrix"
                    except (TypeError, ValueError):
                        pass
            return "zonal_matrix" if "is_rock_outcrop" in out_df.columns and str(ro_v).strip().lower() in ("false", "0", "f", "no") else ""

        for idx in active_indices:
            lat_v = lats[idx]
            lon_v = lons[idx]
            sp_val = str(out_df.iloc[idx]["species_standardized"])
            if lat_v is None or lon_v is None:
                non_geo_indices.append(idx)
            else:
                hab_k = _row_hab_stratum(idx)
                grp_k = f"{sp_val}||{hab_k}" if hab_k else sp_val
                sp_geo_groups.setdefault(grp_k, []).append(idx)

        thinned_indices: List[int] = []
        for _sp_k, sp_idxs in sp_geo_groups.items():
            sp_sorted = sorted(sp_idxs, key=lambda i: priority_tuples[i], reverse=True)
            kept_sp: List[int] = []
            for cand_i in sp_sorted:
                if is_type_flags[cand_i] or is_preserved_flags[cand_i] or (preserve_all_gold and is_gold_flags[cand_i]):
                    kept_sp.append(cand_i)
                    continue
                c_lat = lats[cand_i]
                c_lon = lons[cand_i]
                too_close = any(
                    haversine_distance_meters(c_lat, c_lon, lats[k_i], lons[k_i]) < thin_m
                    for k_i in kept_sp
                )
                if not too_close:
                    kept_sp.append(cand_i)
            thinned_indices.extend(kept_sp)

        thinned_indices.extend(non_geo_indices)
        active_indices = sorted(set(thinned_indices))

    def _select_spatially_spread_quota(
        cand_pool: List[int],
        quota: int,
        stratum_values: Dict[int, str],
    ) -> List[int]:
        """
        Selects up to `quota` indices from `cand_pool`, unconditionally preserving
        all nomenclatural type specimens (and all Gold anchors when `preserve_all_gold=True`)
        and spreading remaining selections across `stratum_values` (counties or states)
        and farthest Haversine distance.
        """
        if len(cand_pool) <= quota:
            return sorted(cand_pool)

        if preserve_all_gold:
            gold_idxs = [i for i in cand_pool if is_gold_flags[i] or is_preserved_flags[i]]
            rem_pool = [i for i in cand_pool if not (is_gold_flags[i] or is_preserved_flags[i])]
        else:
            gold_idxs = [i for i in cand_pool if is_type_flags[i] or is_preserved_flags[i]]
            rem_pool = [i for i in cand_pool if not (is_type_flags[i] or is_preserved_flags[i])]

        rem_pool.sort(key=lambda i: priority_tuples[i], reverse=True)
        selected: List[int] = list(gold_idxs)
        if len(selected) >= quota:
            return (
                sorted(selected)
                if (preserve_all_gold or any(is_type_flags[i] or is_preserved_flags[i] for i in selected))
                else sorted(selected[:quota])
            )

        stratum_counts: Dict[str, int] = {}
        for i in selected:
            st_key = stratum_values[i]
            stratum_counts[st_key] = stratum_counts.get(st_key, 0) + 1

        while len(selected) < quota and rem_pool:
            min_cnt = min(stratum_counts.get(stratum_values[i], 0) for i in rem_pool)
            eligible = [i for i in rem_pool if stratum_counts.get(stratum_values[i], 0) == min_cnt]

            sel_geo = [i for i in selected if lats[i] is not None and lons[i] is not None]
            if not sel_geo:
                best_i = eligible[0]
            else:
                best_i = eligible[0]
                best_score: Tuple[float, Tuple[int, int, int, float, float, str]] = (-1.0, priority_tuples[best_i])
                for cand_i in eligible:
                    if lats[cand_i] is not None and lons[cand_i] is not None:
                        min_d_km = min(
                            haversine_distance_meters(
                                lats[cand_i], lons[cand_i], lats[s_i], lons[s_i]
                            ) / 1000.0
                            for s_i in sel_geo
                        )
                    else:
                        min_d_km = 0.0
                    cand_score = (round(min_d_km, 2), priority_tuples[cand_i])
                    if cand_score > best_score:
                        best_score = cand_score
                        best_i = cand_i

            selected.append(best_i)
            st_key = stratum_values[best_i]
            stratum_counts[st_key] = stratum_counts.get(st_key, 0) + 1
            rem_pool.remove(best_i)

        return sorted(selected)

    # Step D.1: Per-State Quota Balancing (max_per_state > 0)
    if max_per_state is not None and int(max_per_state) > 0:
        state_quota = int(max_per_state)
        state_groups: Dict[Tuple[str, str], List[int]] = {}
        county_labels: Dict[int, str] = {}
        for idx in active_indices:
            row = out_df.iloc[idx]
            sp_val = str(row["species_standardized"])
            st_raw = row.get("stateProvince", "")
            st_val = str(st_raw if not pd.isna(st_raw) else "").strip().upper() or "UNKNOWN_STATE"
            cnty_raw = row.get("county", "")
            cnty_val = str(cnty_raw if not pd.isna(cnty_raw) else "").strip().upper()
            if not cnty_val:
                if lats[idx] is not None and lons[idx] is not None:
                    cnty_val = f"GRID_{int(round(lats[idx] * 4.0))}_{int(round(lons[idx] * 4.0))}"
                else:
                    cnty_val = f"UNSPECIFIED_{idx}"
            county_labels[idx] = cnty_val
            state_groups.setdefault((sp_val, st_val), []).append(idx)

        state_balanced: List[int] = []
        for _grp_key, grp_idxs in state_groups.items():
            state_balanced.extend(_select_spatially_spread_quota(grp_idxs, state_quota, county_labels))
        active_indices = sorted(set(state_balanced))

    # Step D.2: Per-Biogeographic-Region Quota Balancing (max_per_region > 0)
    if max_per_region is not None and int(max_per_region) > 0:
        region_quota = int(max_per_region)
        region_groups: Dict[Tuple[str, str], List[int]] = {}
        state_labels: Dict[int, str] = {}
        for idx in active_indices:
            row = out_df.iloc[idx]
            sp_val = str(row["species_standardized"])
            reg_raw = row.get("regional_group", "")
            reg_val = str(reg_raw if not pd.isna(reg_raw) else "").strip() or "Unassigned_NA"
            st_raw = row.get("stateProvince", "")
            st_val = str(st_raw if not pd.isna(st_raw) else "").strip().upper() or "UNKNOWN_STATE"
            state_labels[idx] = st_val
            region_groups.setdefault((sp_val, reg_val), []).append(idx)

        region_balanced: List[int] = []
        for _grp_key, grp_idxs in region_groups.items():
            region_balanced.extend(_select_spatially_spread_quota(grp_idxs, region_quota, state_labels))
        active_indices = sorted(set(region_balanced))

    out_df = out_df.iloc[active_indices].reset_index(drop=True)
    if out_df.empty:
        out_df["spatial_balance_weight"] = pd.Series(dtype=float)
        return out_df

    # Step E: Compute Continuous `spatial_balance_weight` strictly over non-quarantined active anchors
    ret_lats = np.array(
        [lats[i] if lats[i] is not None else np.nan for i in active_indices],
        dtype=float,
    )
    ret_lons = np.array(
        [lons[i] if lons[i] is not None else np.nan for i in active_indices],
        dtype=float,
    )
    ret_taxa = out_df["species_standardized"].astype(str).tolist()
    has_auth_wt = "authority_weight" in out_df.columns
    eligible_anchor_mask = np.ones(len(out_df), dtype=bool)
    for k in range(len(out_df)):
        row_k = out_df.iloc[k]
        if _as_bool_flag(row_k.get("is_centroid_or_inst_coord")):
            eligible_anchor_mask[k] = False
            continue
        raw_ts_k = row_k.get("type_status") if not pd.isna(row_k.get("type_status", "")) else row_k.get("typeStatus")
        if classify_nomenclatural_type_rank(raw_ts_k)[0] >= 1:
            continue
        if has_auth_wt and not pd.isna(row_k.get("authority_weight")):
            try:
                if float(row_k.get("authority_weight")) <= 0.0:
                    eligible_anchor_mask[k] = False
                    continue
            except (ValueError, TypeError):
                pass
        tier_k = str(row_k.get("determiner_tier", "") if not pd.isna(row_k.get("determiner_tier", "")) else "").strip()
        if tier_k == "Tier_3_Bronze" and not _as_bool_flag(row_k.get("is_promoted_bronze")):
            eligible_anchor_mask[k] = False

    spatial_weights = np.ones(len(out_df), dtype=float)
    bandwidth_km = 25.0

    taxon_to_rows: Dict[str, List[int]] = {}
    for k, sp_name in enumerate(ret_taxa):
        taxon_to_rows.setdefault(sp_name, []).append(k)

    for _sp_name, t_idxs in taxon_to_rows.items():
        t_arr = np.array(t_idxs, dtype=int)
        valid_geo_mask = np.isfinite(ret_lats[t_arr]) & np.isfinite(ret_lons[t_arr]) & eligible_anchor_mask[t_arr]
        geo_idxs = t_arr[valid_geo_mask]
        if len(geo_idxs) == 0:
            continue

        # Continuous geodesic 25 km + 75 km multi-scale Gaussian kernel specimen density d_local
        # (free of arbitrary US state political-area penalties)
        rad_lat = np.radians(ret_lats[geo_idxs])
        rad_lon = np.radians(ret_lons[geo_idxs])
        dlat_mat = rad_lat[:, None] - rad_lat[None, :]
        dlon_mat = rad_lon[:, None] - rad_lon[None, :]
        a_mat = (
            np.sin(dlat_mat / 2.0) ** 2
            + np.cos(rad_lat[:, None]) * np.cos(rad_lat[None, :]) * (np.sin(dlon_mat / 2.0) ** 2)
        )
        dist_km_mat = 2.0 * 6371.0 * np.arcsin(np.clip(np.sqrt(a_mat), 0.0, 1.0))
        kern_fine = np.exp(-(dist_km_mat ** 2) / (2.0 * (bandwidth_km ** 2)))
        kern_regional = 0.35 * np.exp(-(dist_km_mat ** 2) / (2.0 * ((3.0 * bandwidth_km) ** 2)))
        kern_mat = kern_fine + kern_regional
        np.fill_diagonal(kern_mat, 0.0)
        d_local = 1.0 + np.sum(kern_mat, axis=1)

        w_raw = 1.0 / np.sqrt(np.maximum(d_local, 1e-8))
        mean_raw = float(np.mean(w_raw))
        if mean_raw > 0.0:
            w_norm = w_raw / mean_raw
        else:
            w_norm = np.ones_like(w_raw)

        w_spatial = np.clip(w_norm, 0.20, 3.00)
        for _ in range(3):
            m_cur = float(np.mean(w_spatial))
            if m_cur > 0.0:
                w_spatial = np.clip(w_spatial / m_cur, 0.20, 3.00)

        spatial_weights[geo_idxs] = np.round(w_spatial, 6)

    out_df["spatial_balance_weight"] = spatial_weights
    return out_df
