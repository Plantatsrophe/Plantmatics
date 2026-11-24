"""
Module: step01_usda_voucher_filter.py
Description:
    USDA PLANTS county-range voucher filtering with 1-ring adjacency, geodesic buffer,
    type/Gold exemptions, and NatureServe subnational fallback evaluation.
"""

from __future__ import annotations

import logging
import math
from pathlib import Path
import re
import sqlite3
from typing import Any, Callable, Dict, List, Mapping, Optional, Set, Tuple, Union

try:
    import pandas as pd
except ImportError:
    pd = None  # type: ignore[assignment]
try:
    from shapely.ops import unary_union
except ImportError:
    unary_union = None  # type: ignore[assignment]

from scripts.core.step01_usda_natureserve_utils import (  # noqa: F401
    _evaluate_natureserve_voucher,
    _load_natureserve_ranges_for_taxa,
    _resolve_voucher_country_and_subnation,
)

PROJECT_ROOT = Path(__file__).resolve().parents[2]


def filter_vouchers_by_usda_county_range(
    df: "pd.DataFrame",
    standardize_taxon_fn: Callable[[Any], Optional[str]],
    synonym_to_canonical: Optional[Mapping[str, str]] = None,
    county_buffer_km: float = 50.0,
    exempt_gold: bool = True,
    cache_dir: Optional[Path] = None,
    county_geojson_path: Optional[Path] = None,
    db_path: Optional[Union[str, Path]] = None,
    workspace_dir: Path = PROJECT_ROOT,
    allow_network: bool = True,
    force_refresh: bool = False,
    drop_out_of_range: bool = True,
    logger: Optional[logging.Logger] = None,
    target_genus: Optional[str] = None,
    conn: Optional[sqlite3.Connection] = None,
) -> Tuple["pd.DataFrame", Dict[str, int]]:
    """
    Evaluate harvested vouchers against USDA PLANTS county distributions (with 1-ring
    county adjacency and `county_buffer_km` geodesic buffer) and SQLite
    `natureserve_subnational_ranges` fallback for unmapped taxa or Canadian vouchers.
    Records `usda_symbol` and `usda_county_status` (`in_range`, `buffered_adjacent`,
    `natureserve_state_in_range`, `natureserve_exotic_or_waif`, `out_of_range`,
    `unmapped_taxon`, `non_us_or_unresolved`) for SQLite persistence. When
    `drop_out_of_range=False`, all vouchers are retained in soft audit mode; when
    `drop_out_of_range=True`, non-exempt `out_of_range` records are removed prior to
    Step 01 image downloading.

    Exemptions (when `drop_out_of_range=True`):
      - Verified nomenclatural type specimens (`type_rank >= 1`) are unconditionally retained.
      - `Tier_1_Gold` specialist determinations are retained when `exempt_gold=True`.
      - Taxa with zero mapped US counties and zero NatureServe subnational records
        (`unmapped_taxon`) and non-US/CA or unresolved occurrences (`non_us_or_unresolved`)
        fail open so valid records are never dropped due to missing checklist coverage.
    """
    import scripts.core.usda_plants_ranges as _usda

    stats: Dict[str, int] = {
        "evaluated": 0,
        "in_range": 0,
        "buffered_adjacent": 0,
        "natureserve_state_in_range": 0,
        "natureserve_exotic_or_waif": 0,
        "exempt_retained": 0,
        "unmapped_or_non_us": 0,
        "dropped_out_of_range": 0,
    }
    if df is None or df.empty:
        return df, stats

    out_df = df.copy()
    stats["evaluated"] = len(out_df)

    # Build reverse synonym map: canonical_taxon -> sorted synonym strings
    canon_to_syns: Dict[str, List[str]] = {}
    for raw_low, canon_val in (synonym_to_canonical or {}).items():
        if canon_val:
            canon_to_syns.setdefault(canon_val, []).append(raw_low)

    # Resolve canonical target-genus taxon per row
    eff_genus = _usda.resolve_target_genus(explicit_genus=target_genus)
    genus_prefix = f"{eff_genus.strip()} "
    resolved_taxa: List[str] = []
    for _, row in out_df.iterrows():
        std_tx: Optional[str] = None
        for col in ("species_standardized", "scientificName", "species_raw"):
            if col in row and row[col] is not None and not pd.isna(row[col]) and str(row[col]).strip():
                std_tx = standardize_taxon_fn(row[col])
                if std_tx:
                    break
        resolved_taxa.append(std_tx or "")

    unique_taxa = sorted({tx for tx in resolved_taxa if tx.startswith(genus_prefix)})
    if not unique_taxa:
        out_df["usda_symbol"] = ""
        out_df["usda_county_status"] = "unmapped_taxon"
        stats["unmapped_or_non_us"] = len(out_df)
        return out_df, stats

    taxon_ranges: Dict[str, Dict[str, Any]] = {}
    for tx in unique_taxa:
        taxon_ranges[tx] = _usda.fetch_or_load_usda_taxon_range(
            canonical_packera=tx,
            synonym_names=canon_to_syns.get(tx, ()),
            cache_dir=cache_dir,
            db_path=db_path,
            workspace_dir=workspace_dir,
            allow_network=allow_network,
            force_refresh=force_refresh,
            logger=logger,
        )

    ns_ranges_by_taxon = _usda._load_natureserve_ranges_for_taxa(
        unique_taxa,
        conn=conn,
        db_path=db_path,
        workspace_dir=workspace_dir,
        cache_dir=cache_dir,
    )

    has_any_mapped = any(
        bool(r["counties"] or r["state_only_states"])
        for r in taxon_ranges.values()
    ) or any(bool(v) for v in ns_ranges_by_taxon.values())
    if not has_any_mapped:
        out_df["usda_symbol"] = ""
        out_df["usda_county_status"] = "unmapped_taxon"
        stats["unmapped_or_non_us"] = len(out_df)
        return out_df, stats

    spatial_graph = _usda.load_us_county_spatial_graph(
        county_path=county_geojson_path,
        allow_download=allow_network,
    )
    counties_dict = spatial_graph["counties"]
    neighbors_map: Mapping[str, Set[str]] = spatial_graph["neighbors"]

    # Precompute 1-ring neighbor sets and unioned Shapely geometries per mapped taxon
    taxon_expanded_neighbors: Dict[str, Set[str]] = {}
    taxon_doc_union_geom: Dict[str, Any] = {}
    for tx, r_info in taxon_ranges.items():
        doc_cks = _usda._expand_virginia_city_aliases(r_info["counties"])
        adj_set: Set[str] = set()
        doc_geoms_tx: List[Any] = []
        for ck in doc_cks:
            adj_set.update(neighbors_map.get(ck, ()))
            c_rec = counties_dict.get(ck)
            if c_rec and c_rec.get("geom") is not None:
                doc_geoms_tx.append(c_rec["geom"])
        adj_set -= doc_cks
        taxon_expanded_neighbors[tx] = adj_set
        taxon_doc_union_geom[tx] = (
            unary_union(doc_geoms_tx) if (unary_union is not None and doc_geoms_tx) else None
        )

    keep_mask: List[bool] = []
    usda_symbols_col: List[str] = []
    usda_status_col: List[str] = []
    dropped_by_taxon: Dict[str, int] = {}

    for idx_pos, (_, row) in enumerate(out_df.iterrows()):
        tx = resolved_taxa[idx_pos]
        r_info = taxon_ranges.get(tx)
        ns_map = ns_ranges_by_taxon.get(tx, {})
        has_usda = bool(r_info and (r_info["counties"] or r_info["state_only_states"]))
        sym_str = ";".join(r_info["symbols"]) if (r_info and r_info["symbols"]) else ""
        usda_symbols_col.append(sym_str)

        if not has_usda and not ns_map:
            usda_status_col.append("unmapped_taxon")
            stats["unmapped_or_non_us"] += 1
            keep_mask.append(True)
            continue

        raw_ts = row.get("type_status")
        if raw_ts is None or pd.isna(raw_ts) or not str(raw_ts).strip():
            raw_ts = row.get("typeStatus")
        type_rank, _ = _usda.classify_nomenclatural_type_rank(raw_ts)
        is_type = type_rank >= 1
        det_tier = str(row.get("determiner_tier") or "").strip()
        is_exempt = is_type or (exempt_gold and det_tier == "Tier_1_Gold")

        lat_v: Optional[float] = None
        lon_v: Optional[float] = None
        for c_lat in ("decimalLatitude",):
            if c_lat in row and row[c_lat] is not None and not pd.isna(row[c_lat]):
                try:
                    fv = float(row[c_lat])
                    if math.isfinite(fv):
                        lat_v = fv
                        break
                except (TypeError, ValueError):
                    pass
        for c_lon in ("decimalLongitude",):
            if c_lon in row and row[c_lon] is not None and not pd.isna(row[c_lon]):
                try:
                    fv = float(row[c_lon])
                    if math.isfinite(fv):
                        lon_v = fv
                        break
                except (TypeError, ValueError):
                    pass

        st_raw = None if pd.isna(row.get("stateProvince")) else str(row.get("stateProvince") or "")
        co_raw = None if pd.isna(row.get("county")) else str(row.get("county") or "")
        st_norm, co_norm, is_non_us = _usda.resolve_voucher_us_county(
            lat=lat_v,
            lon=lon_v,
            state_province=st_raw,
            county=co_raw,
            spatial_graph=spatial_graph,
        )
        c_code, sub_code = _usda._resolve_voucher_country_and_subnation(row, st_raw=st_raw, st_norm=st_norm)

        # Fallback to NatureServe subnational ranges when USDA lacks taxon counties OR voucher is Canadian
        if (not has_usda) or (c_code == "CA" and ns_map):
            ns_status = _usda._evaluate_natureserve_voucher(
                c_code=c_code,
                sub_code=sub_code,
                ns_map=ns_map,
                extra_us_states=(r_info or {}).get("documented_states") if has_usda else None,
            )
            if ns_status == "natureserve_state_in_range":
                usda_status_col.append(ns_status)
                stats["natureserve_state_in_range"] += 1
                keep_mask.append(True)
            elif ns_status == "natureserve_exotic_or_waif":
                usda_status_col.append(ns_status)
                stats["natureserve_exotic_or_waif"] += 1
                keep_mask.append(True)
            elif ns_status == "buffered_adjacent":
                usda_status_col.append("buffered_adjacent")
                stats["buffered_adjacent"] += 1
                keep_mask.append(True)
            elif ns_status == "non_us_or_unresolved":
                usda_status_col.append("non_us_or_unresolved")
                stats["unmapped_or_non_us"] += 1
                keep_mask.append(True)
            else:
                usda_status_col.append("out_of_range")
                if is_exempt:
                    stats["exempt_retained"] += 1
                    keep_mask.append(True)
                elif not drop_out_of_range:
                    keep_mask.append(True)
                else:
                    stats["dropped_out_of_range"] += 1
                    dropped_by_taxon[tx] = dropped_by_taxon.get(tx, 0) + 1
                    keep_mask.append(False)
            continue

        if is_non_us or (st_norm is None and lat_v is None and lon_v is None):
            usda_status_col.append("non_us_or_unresolved")
            stats["unmapped_or_non_us"] += 1
            keep_mask.append(True)
            continue

        doc_cks = _usda._expand_virginia_city_aliases(r_info["counties"])
        st_only = r_info["state_only_states"]
        doc_states = r_info["documented_states"]
        voucher_ck = f"{st_norm}|{co_norm}" if (st_norm and co_norm) else None
        va_alias_hit = False
        if st_norm == "virginia" and co_norm:
            co_va_base = re.sub(r"\s+city$", "", co_norm.strip().lower()).strip()
            va_alias_hit = (
                f"virginia|{co_va_base}" in doc_cks
                or f"virginia|{co_va_base} city" in doc_cks
            )

        if (st_norm and st_norm in st_only) or (voucher_ck and voucher_ck in doc_cks) or va_alias_hit:
            usda_status_col.append("in_range")
            stats["in_range"] += 1
            keep_mask.append(True)
            continue

        # Ungeoreferenced state-only label without county name: retain if state has documented presence
        if st_norm and co_norm is None and lat_v is None and st_norm in doc_states:
            usda_status_col.append("in_range")
            stats["in_range"] += 1
            keep_mask.append(True)
            continue

        is_adj = bool(voucher_ck and voucher_ck in taxon_expanded_neighbors.get(tx, set()))
        if not is_adj and st_norm == "virginia" and co_norm:
            co_va_base = re.sub(r"\s+city$", "", co_norm.strip().lower()).strip()
            is_adj = (
                f"virginia|{co_va_base}" in taxon_expanded_neighbors.get(tx, set())
                or f"virginia|{co_va_base} city" in taxon_expanded_neighbors.get(tx, set())
            )
        if not is_adj and county_buffer_km > 0.0:
            is_adj = _usda._is_within_geodesic_buffer(
                lat_f=lat_v,
                lon_f=lon_v,
                voucher_ck=voucher_ck,
                doc_counties=doc_cks,
                doc_union_geom=taxon_doc_union_geom.get(tx),
                spatial_graph=spatial_graph,
                buffer_km=county_buffer_km,
            )

        if is_adj:
            usda_status_col.append("buffered_adjacent")
            stats["buffered_adjacent"] += 1
            keep_mask.append(True)
            continue

        if c_code == "US" and sub_code and ns_map:
            ns_rec = ns_map.get(("US", sub_code))
            if ns_rec is not None and not ns_rec.get("is_exotic"):
                usda_status_col.append("natureserve_state_in_range")
                stats["natureserve_state_in_range"] += 1
                keep_mask.append(True)
                continue

        if st_norm == "texas" and "texas" in doc_states and (lon_v is None or lon_v >= -96.2):
            usda_status_col.append("in_range")
            stats["in_range"] += 1
            keep_mask.append(True)
            continue

        usda_status_col.append("out_of_range")
        if is_exempt:
            stats["exempt_retained"] += 1
            keep_mask.append(True)
        elif not drop_out_of_range:
            keep_mask.append(True)
        else:
            stats["dropped_out_of_range"] += 1
            dropped_by_taxon[tx] = dropped_by_taxon.get(tx, 0) + 1
            keep_mask.append(False)

    out_df["usda_symbol"] = usda_symbols_col
    out_df["usda_county_status"] = usda_status_col
    filtered_df = out_df.loc[keep_mask].reset_index(drop=True)

    if logger:
        logger.info(
            f"USDA PLANTS County Range Filter (buffer={county_buffer_km:.1f} km, exempt_gold={exempt_gold}): "
            f"Evaluated {stats['evaluated']} -> Retained {len(filtered_df)} "
            f"[in_range={stats['in_range']}, buffered_adjacent={stats['buffered_adjacent']}, "
            f"exempt={stats['exempt_retained']}, unmapped/non_us={stats['unmapped_or_non_us']}] | "
            f"Dropped out-of-range={stats['dropped_out_of_range']}."
        )
        if dropped_by_taxon:
            top_drops = ", ".join(f"{k}: {v}" for k, v in sorted(dropped_by_taxon.items()))
            logger.info(f"  Out-of-range vouchers dropped by taxon: {top_drops}")

    return filtered_df, stats
