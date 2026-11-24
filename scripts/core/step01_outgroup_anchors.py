"""step01_outgroup_anchors.py - Cascading type-, authority-, and geography-stratified outgroup anchor selection."""

from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, List, Optional, Set, Tuple

try:
    import pandas as pd
except ImportError:
    pd = None  # type: ignore[assignment]

from scripts.core.config import PROJECT_ROOT
from scripts.core.step01_determiner_surnames import (
    classify_nomenclatural_type_rank,
    classify_nomenclatural_type_subrank,
)
from scripts.core.step01_exsiccatae_fingerprint import (
    _as_bool_flag,
    _is_extralimital_status,
    get_image_resolution_score,
    resolve_voucher_image_path,
)


def select_cascading_outgroup_anchors(
    df: pd.DataFrame,
    outgroup_taxa: List[str],
    max_anchors: int = 30,
    workspace_dir: Optional[Path] = None,
    focal_df: Optional[pd.DataFrame] = None,
) -> Tuple[List[pd.DataFrame], Dict[str, int], Dict[str, Dict[str, int]]]:
    """
    Selects up to `max_anchors` reference sheets per outgroup taxon by cascading across
    nomenclatural types, range-clean standing, and determiner authority tiers
    (`Tier_1_Gold` -> `Tier_2_Silver` -> `Tier_3_Bronze`), with geographic round-robin
    stratification across `(stateProvince, regional_group)` within each stratum
    prioritizing regions overlapping `focal_df`.
    """
    outgroup_dfs: List[pd.DataFrame] = []
    outgroup_counts: Dict[str, int] = {}
    outgroup_tier_breakdown: Dict[str, Dict[str, int]] = {}

    if df is None or df.empty or not outgroup_taxa or max_anchors <= 0:
        for tx in outgroup_taxa or []:
            outgroup_counts[tx] = 0
            outgroup_tier_breakdown[tx] = {}
        return outgroup_dfs, outgroup_counts, outgroup_tier_breakdown

    tier_rank = {"Tier_1_Gold": 3, "Tier_2_Silver": 2, "Tier_3_Bronze": 1}
    ws_root = Path(workspace_dir).resolve() if workspace_dir is not None else PROJECT_ROOT

    def _norm_tok(v: Any) -> str:
        if v is None or (pd is not None and pd.isna(v)):
            return ""
        return str(v).strip().lower()

    focal_states: Set[str] = set()
    focal_regions: Set[str] = set()
    if focal_df is not None and not focal_df.empty:
        if "stateProvince" in focal_df.columns:
            focal_states = {_norm_tok(v) for v in focal_df["stateProvince"].tolist()} - {""}
        if "regional_group" in focal_df.columns:
            focal_regions = {_norm_tok(v) for v in focal_df["regional_group"].tolist()} - {""}

    for tx in outgroup_taxa:
        tx_df = df[df["species_standardized"] == tx].copy()
        if tx_df.empty:
            outgroup_counts[tx] = 0
            outgroup_tier_breakdown[tx] = {}
            continue

        def _has_local(row: pd.Series) -> int:
            img_p = str(row.get("image_path", "") if not pd.isna(row.get("image_path", "")) else "").strip()
            cat_str = str(row.get("catalogNumber", "") if not pd.isna(row.get("catalogNumber", "")) else "").strip()
            cat_rel = f"data/raw_vouchers/{cat_str}.webp" if cat_str else ""
            if img_p and resolve_voucher_image_path(img_p, workspace_dir=ws_root).is_file():
                return 1
            if cat_rel and resolve_voucher_image_path(cat_rel, workspace_dir=ws_root).is_file():
                return 1
            return 0

        def _has_url_or_path(row: pd.Series) -> bool:
            for c in ("image_path", "_image_url", "image_url", "identifier", "good_quality_AccessURI"):
                v = row.get(c, None)
                if v is not None and not pd.isna(v) and str(v).strip():
                    return True
            return False

        def _is_prim(row: pd.Series) -> int:
            if "is_primary_duplicate" not in row.index or pd.isna(row.get("is_primary_duplicate")):
                return 1
            v = row.get("is_primary_duplicate")
            return 1 if (v is True or str(v).strip().lower() in {"true", "1", "t", "yes"}) else 0

        tx_df["_has_local"] = tx_df.apply(_has_local, axis=1)
        tx_df["_has_valid"] = (tx_df["_has_local"] > 0) | tx_df.apply(_has_url_or_path, axis=1)
        tx_df = tx_df[tx_df["_has_valid"]].copy()
        if tx_df.empty:
            outgroup_counts[tx] = 0
            outgroup_tier_breakdown[tx] = {}
            continue

        ts_col = "type_status" if "type_status" in tx_df.columns else ("typeStatus" if "typeStatus" in tx_df.columns else None)
        tx_df["_type_rank"] = (
            tx_df[ts_col].map(lambda t: classify_nomenclatural_type_rank(t)[0])
            if ts_col
            else 0
        )
        tx_df["_type_subrank"] = (
            tx_df[ts_col].map(classify_nomenclatural_type_subrank)
            if ts_col
            else 0.0
        )
        tx_df["_tier_rank"] = (
            tx_df["determiner_tier"].astype(str).str.strip().map(lambda t: tier_rank.get(t, 0))
            if "determiner_tier" in tx_df.columns
            else 0
        )
        tx_df["_is_prim"] = tx_df.apply(_is_prim, axis=1)
        tx_df["_range_clean_rank"] = [
            0 if (_is_extralimital_status(r) or _as_bool_flag(r.get("exsiccatae_taxonomic_conflict"))) else 1
            for _, r in tx_df.iterrows()
        ]
        tx_df["_img_res"] = tx_df["image_path"].map(
            lambda p: get_image_resolution_score(str(p) if not pd.isna(p) else "", workspace_dir=ws_root)
        ) if "image_path" in tx_df.columns else 0.0

        tx_df = tx_df.sort_values(
            by=["_type_rank", "_type_subrank", "_range_clean_rank", "_tier_rank", "_is_prim", "_has_local", "_img_res"],
            ascending=[False, False, False, False, False, False, False],
            kind="mergesort",
        )

        # Within each (type_rank, range_clean_rank, tier_rank) cascade stratum, round-robin
        # across (stateProvince, regional_group) buckets with focal-overlapping buckets first
        selected_positions: List[int] = []
        cap_n = int(max_anchors)
        strata_groups: Dict[Tuple[int, int, int], List[int]] = {}
        strata_order: List[Tuple[int, int, int]] = []
        for pos_i, (_, row_i) in enumerate(tx_df.iterrows()):
            s_key = (int(row_i["_type_rank"]), int(row_i["_range_clean_rank"]), int(row_i["_tier_rank"]))
            if s_key not in strata_groups:
                strata_groups[s_key] = []
                strata_order.append(s_key)
            strata_groups[s_key].append(pos_i)

        for s_key in strata_order:
            if len(selected_positions) >= cap_n:
                break
            stratum_positions = strata_groups[s_key]
            buckets: Dict[Tuple[str, str], List[int]] = {}
            bucket_order: List[Tuple[str, str]] = []
            for p_i in stratum_positions:
                r_p = tx_df.iloc[p_i]
                b_key = (_norm_tok(r_p.get("stateProvince")), _norm_tok(r_p.get("regional_group")))
                if b_key not in buckets:
                    buckets[b_key] = []
                    bucket_order.append(b_key)
                buckets[b_key].append(p_i)
            if focal_states or focal_regions:
                focal_buckets = [
                    bk for bk in bucket_order
                    if (bk[0] and bk[0] in focal_states) or (bk[1] and bk[1] in focal_regions)
                ]
                other_buckets = [bk for bk in bucket_order if bk not in focal_buckets]
                bucket_order = focal_buckets + other_buckets
            wave_idx = 0
            while len(selected_positions) < cap_n:
                added_wave = False
                for bk in bucket_order:
                    b_list = buckets[bk]
                    if wave_idx < len(b_list):
                        selected_positions.append(b_list[wave_idx])
                        added_wave = True
                        if len(selected_positions) >= cap_n:
                            break
                if not added_wave:
                    break
                wave_idx += 1

        tx_head = tx_df.iloc[selected_positions]
        n_type_sel = int((tx_head["_type_rank"] >= 1).sum())
        n_prim_type_sel = int((tx_head["_type_rank"] == 2).sum())
        n_sec_type_sel = int((tx_head["_type_rank"] == 1).sum())
        tx_sel = tx_head.drop(
            columns=["_has_local", "_has_valid", "_type_rank", "_type_subrank", "_range_clean_rank", "_tier_rank", "_is_prim", "_img_res"]
        )
        outgroup_dfs.append(tx_sel)
        outgroup_counts[tx] = len(tx_sel)
        if "determiner_tier" in tx_sel.columns:
            breakdown = {
                str(k): int(v) for k, v in tx_sel["determiner_tier"].astype(str).value_counts().items()
            }
        else:
            breakdown = {"Unspecified": len(tx_sel)}
        breakdown["Type_Specimens"] = n_type_sel
        if n_prim_type_sel > 0:
            breakdown["Tier_0_Primary_Type"] = n_prim_type_sel
        if n_sec_type_sel > 0:
            breakdown["Tier_0_Secondary_Type"] = n_sec_type_sel
        outgroup_tier_breakdown[tx] = breakdown

    return outgroup_dfs, outgroup_counts, outgroup_tier_breakdown
