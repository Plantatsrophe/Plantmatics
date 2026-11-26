"""step01_duplicate_stratifier.py - Exsiccatae duplicate stratification and authority-first ranking."""

from __future__ import annotations

import logging
from pathlib import Path
import re
from typing import Any, Dict, List, Optional, Tuple

try:
    import pandas as pd
except ImportError:
    pd = None  # type: ignore[assignment]

from scripts.core.config import DEFAULT_HOME_HERBARIUM
from scripts.core.step01_determiner_surnames import classify_nomenclatural_type_rank, classify_nomenclatural_type_subrank
from scripts.core.step01_exsiccatae_fingerprint import (
    _as_bool_flag, _is_extralimital_status, _resolve_row_standardized_packera,
    get_image_resolution_score, haversine_distance_meters, normalize_collector,
    normalize_event_date, normalize_record_number, resolve_voucher_image_path,
)


def stratify_duplicates(
    df: pd.DataFrame,
    max_distance_meters: float = 50.0,
    workspace_dir: Optional[Path] = None,
    logger: Optional[logging.Logger] = None,
    home_herbarium: str = DEFAULT_HOME_HERBARIUM,
) -> pd.DataFrame:
    """
    Identifies duplicate voucher specimens (exsiccatae) representing the same
    field collection event distributed across multiple herbaria.

    Duplicates are resolved by:
      1. Standardized exsiccatae fingerprint (norm_collector + norm_number + norm_date).
      2. If collection number is missing, falls back to geographic coordinate proximity (< 50 m)
         within identical collection dates.

    Designates exactly one primary voucher per duplicate group using authority-first hierarchy:
      Priority 1: Verified Type Status (and ICN Art. 9.6 sub-rank)
      Priority 2: Documented in-range standing
      Priority 3: Determiner tier (`Tier_1_Gold` > `Tier_2_Silver` > `Tier_3_Bronze`)
      Priority 4: Home Herbarium (default: NCU; configurable via home_herbarium)
      Priority 5: Highest specimen image resolution

    Adds/updates columns:
      - exsiccatae_key: Standardized signature string
      - is_primary_duplicate: Boolean flag (True for primary specimen, False for duplicates)
      - duplicate_count: Total sheets in this exsiccatae group
      - exsiccatae_taxonomic_conflict: Boolean flag (True when duplicates disagree in Packera taxon)
      - exsiccatae_gold_taxon: Standardized Packera determination held by Tier_1_Gold duplicate(s)
    """
    if df.empty:
        df["exsiccatae_key"] = pd.Series(dtype=str)
        df["is_primary_duplicate"] = pd.Series(dtype=bool)
        df["duplicate_count"] = pd.Series(dtype=int)
        df["exsiccatae_taxonomic_conflict"] = pd.Series(dtype=bool)
        df["exsiccatae_gold_taxon"] = pd.Series(dtype=str)
        return df

    from scripts.core.taxonomy_authority import sanitize_filename

    orig_index = df.index
    df = df.reset_index(drop=True).copy()

    # Ensure metadata columns exist
    if "recordedBy" not in df.columns:
        df["recordedBy"] = df.get("identifiedBy", "")
    if "recordNumber" not in df.columns:
        df["recordNumber"] = ""
    if "image_path" not in df.columns:
        df["image_path"] = ""

    norm_collectors: List[str] = []
    norm_numbers: List[str] = []
    norm_dates: List[str] = []
    has_coords: List[bool] = []
    lats: List[Optional[float]] = []
    lons: List[Optional[float]] = []

    for _, row in df.iterrows():
        rec_by = str(row.get("recordedBy", "") if not pd.isna(row.get("recordedBy", "")) else "")
        rec_num = row.get("recordNumber", "") if not pd.isna(row.get("recordNumber", "")) else ""
        ev_date = str(row.get("eventDate", "") if not pd.isna(row.get("eventDate", "")) else "")
        yr = row.get("year")
        mo = row.get("month")
        dy = row.get("day")

        c_norm = normalize_collector(rec_by)
        n_norm = normalize_record_number(rec_num)
        d_norm = normalize_event_date(ev_date, year=yr, month=mo, day=dy)

        lat_val = None
        lon_val = None
        for col_lat in ["decimalLatitude"]:
            if col_lat in row and not pd.isna(row[col_lat]):
                try:
                    lat_val = float(row[col_lat])
                    break
                except (ValueError, TypeError):
                    pass
        for col_lon in ["decimalLongitude"]:
            if col_lon in row and not pd.isna(row[col_lon]):
                try:
                    lon_val = float(row[col_lon])
                    break
                except (ValueError, TypeError):
                    pass

        norm_collectors.append(c_norm)
        norm_numbers.append(n_norm)
        norm_dates.append(d_norm)
        has_coords.append(lat_val is not None and lon_val is not None)
        lats.append(lat_val)
        lons.append(lon_val)

    n_rows = len(df)
    group_labels: List[Optional[str]] = [None] * n_rows
    row_std_taxa: List[str] = [_resolve_row_standardized_packera(df.iloc[i]) for i in range(n_rows)]

    # Pass 1: Assign keys for records with explicit collection numbers
    for idx in range(n_rows):
        if norm_numbers[idx]:
            group_labels[idx] = f"{norm_collectors[idx]}_{norm_numbers[idx]}_{norm_dates[idx]}"

    # Pass 2: Missing collection number -> fallback to coordinate proximity (< 50m) on identical dates
    missing_num_indices = [idx for idx in range(n_rows) if not norm_numbers[idx]]
    date_to_indices: Dict[str, List[int]] = {}
    for idx in missing_num_indices:
        d = norm_dates[idx]
        date_to_indices.setdefault(d, []).append(idx)

    for d_val, idx_list in date_to_indices.items():
        if d_val == "unknown_date" or d_val.endswith("-00-00"):
            for i in idx_list:
                cat = sanitize_filename(str(df.iloc[i].get("catalogNumber", i)))
                group_labels[i] = f"{norm_collectors[i]}_sn_{cat}_{d_val}"
            continue

        coord_indices: List[int] = []
        no_coord_indices: List[int] = []
        for i in idx_list:
            is_cent = _as_bool_flag(df.iloc[i].get("is_centroid_or_inst_coord"))
            if not has_coords[i] or norm_collectors[i] == "unknown" or is_cent:
                no_coord_indices.append(i)
            else:
                coord_indices.append(i)

        # Records lacking exact coords, named collector, or centroid hygiene stay unique singletons
        for i in no_coord_indices:
            cat = sanitize_filename(str(df.iloc[i].get("catalogNumber", i)))
            group_labels[i] = f"{norm_collectors[i]}_sn_{cat}_{d_val}"

        # Match unnumbered records to any numbered collection on the same date within max_distance_meters
        numbered_date_indices = [
            i for i in range(n_rows)
            if norm_numbers[i] and norm_dates[i] == d_val and has_coords[i]
            and norm_collectors[i] != "unknown" and not _as_bool_flag(df.iloc[i].get("is_centroid_or_inst_coord"))
        ]

        matched_to_numbered = set()
        for i in coord_indices:
            for num_i in numbered_date_indices:
                same_tx = not row_std_taxa[i] or not row_std_taxa[num_i] or row_std_taxa[i] == row_std_taxa[num_i]
                if same_tx and norm_collectors[i] == norm_collectors[num_i]:
                    dist = haversine_distance_meters(lats[i], lons[i], lats[num_i], lons[num_i])
                    if dist <= max_distance_meters:
                        group_labels[i] = group_labels[num_i]
                        matched_to_numbered.add(i)
                        break

        remaining_coord_indices = [i for i in coord_indices if i not in matched_to_numbered]
        num_rem = len(remaining_coord_indices)
        if num_rem == 0:
            continue

        visited = [False] * num_rem
        for a in range(num_rem):
            if visited[a]:
                continue
            cluster = [a]
            visited[a] = True
            q = [a]
            while q:
                curr = q.pop(0)
                curr_idx = remaining_coord_indices[curr]
                for b in range(num_rem):
                    if not visited[b]:
                        b_idx = remaining_coord_indices[b]
                        same_coll = norm_collectors[curr_idx] == norm_collectors[b_idx] and norm_collectors[curr_idx] != "unknown"
                        same_tx = not row_std_taxa[curr_idx] or not row_std_taxa[b_idx] or row_std_taxa[curr_idx] == row_std_taxa[b_idx]
                        if same_coll and same_tx:
                            dist = haversine_distance_meters(
                                lats[curr_idx], lons[curr_idx], lats[b_idx], lons[b_idx]
                            )
                            if dist <= max_distance_meters:
                                visited[b] = True
                                cluster.append(b)
                                q.append(b)

            actual_indices = [remaining_coord_indices[c_idx] for c_idx in cluster]
            best_coll = "unknown"
            for act_i in actual_indices:
                if norm_collectors[act_i] != "unknown":
                    best_coll = norm_collectors[act_i]
                    break

            if len(actual_indices) > 1:
                anchor_cat = sanitize_filename(str(df.iloc[actual_indices[0]].get("catalogNumber", actual_indices[0])))
                cluster_key = f"{best_coll}_geo_{anchor_cat}_{d_val}"
            else:
                act_i = actual_indices[0]
                cat = sanitize_filename(str(df.iloc[act_i].get("catalogNumber", act_i)))
                cluster_key = f"{best_coll}_sn_{cat}_{d_val}"

            for act_i in actual_indices:
                group_labels[act_i] = cluster_key

    # Pass 2.5: Merge same-institution portal barcode aliases (e.g., 19398 vs WILLI_19398)
    # when stripped catalog suffixes match and collection metadata is compatible.
    alias_groups: Dict[str, List[int]] = {}
    row_cat_prefixes: List[str] = [""] * n_rows
    row_inst_codes: List[str] = [""] * n_rows
    for idx in range(n_rows):
        row = df.iloc[idx]
        raw_inst = row.get("institutionCode", "")
        inst_clean = "" if (raw_inst is None or pd.isna(raw_inst)) else str(raw_inst).strip().upper()
        row_inst_codes[idx] = inst_clean
        raw_cat = row.get("catalogNumber", "")
        cat_clean = "" if (raw_cat is None or pd.isna(raw_cat)) else str(raw_cat).strip()
        m_pref = re.match(r"^([A-Za-z]{2,12})[-_ ]+(.+)$", cat_clean)
        if m_pref:
            row_cat_prefixes[idx] = m_pref.group(1).strip().upper()
            stripped_cat = m_pref.group(2).strip().upper()
        else:
            stripped_cat = cat_clean.upper()
        if len(stripped_cat) >= 3:
            alias_groups.setdefault(stripped_cat, []).append(idx)

    for _stripped_cat, cand_idxs in alias_groups.items():
        if len(cand_idxs) < 2:
            continue
        for pos_a in range(len(cand_idxs)):
            i = cand_idxs[pos_a]
            for pos_b in range(pos_a + 1, len(cand_idxs)):
                j = cand_idxs[pos_b]
                if group_labels[i] == group_labels[j]:
                    continue
                inst_i = row_inst_codes[i]
                inst_j = row_inst_codes[j]
                pref_i = row_cat_prefixes[i]
                pref_j = row_cat_prefixes[j]
                same_portal = (
                    (inst_i != "" and inst_i == inst_j)
                    or (pref_i != "" and pref_i == inst_j)
                    or (pref_j != "" and pref_j == inst_i)
                    or (pref_i != "" and pref_i == pref_j)
                )
                if not same_portal:
                    continue
                coll_ok = (
                    norm_collectors[i] == norm_collectors[j]
                    or norm_collectors[i] == "unknown"
                    or norm_collectors[j] == "unknown"
                )
                num_ok = (
                    norm_numbers[i] == norm_numbers[j]
                    or not norm_numbers[i]
                    or not norm_numbers[j]
                )
                date_ok = (
                    norm_dates[i] == norm_dates[j]
                    or norm_dates[i] == "unknown_date"
                    or norm_dates[j] == "unknown_date"
                )
                coords_ok = True
                if has_coords[i] and has_coords[j]:
                    dist_ij = haversine_distance_meters(lats[i], lons[i], lats[j], lons[j])
                    coords_ok = dist_ij <= max(float(max_distance_meters), 5000.0)
                if coll_ok and num_ok and date_ok and coords_ok:
                    lbl_i = group_labels[i] or ""
                    lbl_j = group_labels[j] or ""
                    if ("_sn_" in lbl_i and "_sn_" not in lbl_j) or (
                        ("_sn_" in lbl_i) == ("_sn_" in lbl_j) and len(lbl_j) < len(lbl_i)
                    ):
                        canon_lbl, old_lbl = lbl_j, lbl_i
                    else:
                        canon_lbl, old_lbl = lbl_i, lbl_j
                    if canon_lbl and old_lbl and canon_lbl != old_lbl:
                        for k_idx in range(n_rows):
                            if group_labels[k_idx] == old_lbl:
                                group_labels[k_idx] = canon_lbl

    for idx in range(n_rows):
        if group_labels[idx] is None:
            cat = sanitize_filename(str(df.iloc[idx].get("catalogNumber", idx)))
            group_labels[idx] = f"{norm_collectors[idx]}_sn_{cat}_{norm_dates[idx]}"

    df["exsiccatae_key"] = group_labels

    # Pass 3: Cross-duplicate taxonomic conflict auditing and authority-first duplicate ranking
    grouped = df.groupby("exsiccatae_key", sort=False)
    df["duplicate_count"] = grouped["catalogNumber"].transform("count")

    idx_to_pos: Dict[Any, int] = {idx_lbl: pos for pos, idx_lbl in enumerate(df.index)}
    is_primary_list = [False] * n_rows
    pre_conf_col = df["exsiccatae_taxonomic_conflict"] if "exsiccatae_taxonomic_conflict" in df.columns else [False] * n_rows
    pre_gold_col = df["exsiccatae_gold_taxon"] if "exsiccatae_gold_taxon" in df.columns else [""] * n_rows
    conflict_list = [_as_bool_flag(v) for v in pre_conf_col]
    gold_taxon_list = [
        "" if (v is None or (pd is not None and pd.isna(v))) else str(v).strip()
        for v in pre_gold_col
    ]

    for key, group_df in grouped:
        indices = list(group_df.index)
        if len(indices) == 1:
            is_primary_list[idx_to_pos[indices[0]]] = True
            continue

        scored_candidates = []
        idx_to_std_tx: Dict[Any, str] = {}
        idx_to_tier: Dict[Any, str] = {}
        idx_to_score: Dict[Any, Tuple[Any, ...]] = {}
        idx_to_extra: Dict[Any, bool] = {}

        for idx in indices:
            row = group_df.loc[idx]

            # Priority 1: Verified Nomenclatural Type Rank (2 = Primary Type, 1 = Secondary Type, 0 = Non-type)
            # and ICN Art. 9.6 fine-grained float sub-rank (2.2 = Holotype-tier, 2.0 = Isotype-tier,
            # 1.5 = Syntype-tier, 1.0 = Paratype-tier, 0.0 = Non-type)
            raw_ts = row.get("type_status")
            if raw_ts is None or (pd is not None and pd.isna(raw_ts)) or not str(raw_ts).strip():
                raw_ts = row.get("typeStatus")
            type_rank, _ = classify_nomenclatural_type_rank(raw_ts)
            type_subrank = classify_nomenclatural_type_subrank(raw_ts)
            has_type = int(type_rank)

            # Priority 1.5: Documented in-range status ahead of extralimital sheets
            is_extra_i = _is_extralimital_status(row)
            in_range_rank = 0 if is_extra_i else 1

            # Priority 2: Home Herbarium (default: NCU, configurable via home_herbarium)
            target_home = str(home_herbarium).strip().upper() if home_herbarium else ""
            inst = str(row.get("institutionCode", "") if not pd.isna(row.get("institutionCode", "")) else "").strip().upper()
            is_home_herb = 1 if (target_home and inst == target_home) else 0

            # Priority 3: Tier 1 Specialist determinations
            tier = str(row.get("determiner_tier", "") if not pd.isna(row.get("determiner_tier", "")) else "").strip()
            tier_score = 3 if tier == "Tier_1_Gold" else (2 if tier == "Tier_2_Silver" else 1)

            # Priority 4: Highest specimen image resolution
            img_path = str(row.get("image_path", "") if not pd.isna(row.get("image_path", "")) else "")
            resolved_img = resolve_voucher_image_path(img_path, workspace_dir=workspace_dir) if img_path else Path("")
            res_score = get_image_resolution_score(str(resolved_img) if img_path else "", workspace_dir=workspace_dir)

            cat_tiebreak = str(row.get("catalogNumber", idx))

            score_tuple = (has_type, type_subrank, tier_score, in_range_rank, is_home_herb, res_score, cat_tiebreak)
            scored_candidates.append((score_tuple, idx))
            idx_to_std_tx[idx] = _resolve_row_standardized_packera(row)
            idx_to_tier[idx] = tier
            idx_to_score[idx] = score_tuple
            idx_to_extra[idx] = is_extra_i

        distinct_taxa = {tx for tx in idx_to_std_tx.values() if tx}
        inrange_taxa = {idx_to_std_tx[idx] for idx in indices if idx_to_std_tx.get(idx) and not idx_to_extra.get(idx, False)}
        has_conflict = len(distinct_taxa) >= 2

        inrange_gold_counts: Dict[str, int] = {}
        inrange_gold_best_idx: Dict[str, Any] = {}
        all_gold_counts: Dict[str, int] = {}
        all_gold_best_idx: Dict[str, Any] = {}
        for idx in indices:
            tx_i = idx_to_std_tx.get(idx, "")
            if idx_to_tier.get(idx) == "Tier_1_Gold" and tx_i:
                all_gold_counts[tx_i] = all_gold_counts.get(tx_i, 0) + 1
                if tx_i not in all_gold_best_idx or idx_to_score[idx] > idx_to_score[all_gold_best_idx[tx_i]]:
                    all_gold_best_idx[tx_i] = idx
                if not idx_to_extra.get(idx, False) and not conflict_list[idx_to_pos[idx]]:
                    inrange_gold_counts[tx_i] = inrange_gold_counts.get(tx_i, 0) + 1
                    if tx_i not in inrange_gold_best_idx or idx_to_score[idx] > idx_to_score[inrange_gold_best_idx[tx_i]]:
                        inrange_gold_best_idx[tx_i] = idx

        best_gold_taxon = ""
        best_gold_source_idx: Any = None
        best_gold_from_inrange = False
        active_gold_counts = inrange_gold_counts if inrange_gold_counts else all_gold_counts
        active_gold_best = inrange_gold_best_idx if inrange_gold_counts else all_gold_best_idx
        if active_gold_counts:
            best_gold_from_inrange = bool(inrange_gold_counts)
            ranked_gold = sorted(
                active_gold_counts.keys(),
                key=lambda tx: (active_gold_counts[tx], idx_to_score[active_gold_best[tx]]),
                reverse=True,
            )
            best_gold_taxon = ranked_gold[0]
            best_gold_source_idx = active_gold_best[best_gold_taxon]

        multi_gold_disagreement = len(active_gold_counts) >= 2
        unanimous_inrange_tx = next(iter(inrange_taxa)) if len(inrange_taxa) == 1 else ""
        for idx in indices:
            pos_i = idx_to_pos[idx]
            row_tx_i = idx_to_std_tx.get(idx, "")
            if not has_conflict:
                row_conflict = False
            elif best_gold_taxon and not multi_gold_disagreement:
                # Unanimous Tier_1_Gold consensus exists: flag only sheets whose label disagreed with Gold consensus
                row_conflict = bool(row_tx_i and row_tx_i != best_gold_taxon)
            elif unanimous_inrange_tx:
                # Unanimous in-range determination across duplicates: flag only disagreeing extralimital duplicates
                row_conflict = bool(idx_to_extra.get(idx, False) and row_tx_i and row_tx_i != unanimous_inrange_tx)
            else:
                row_conflict = True
            conflict_list[pos_i] = bool(conflict_list[pos_i] or row_conflict)
            if best_gold_taxon and (not multi_gold_disagreement or best_gold_from_inrange or idx_to_extra.get(idx, False)):
                gold_taxon_list[pos_i] = best_gold_taxon

        scored_candidates.sort(key=lambda x: x[0], reverse=True)
        primary_idx = scored_candidates[0][1]
        is_primary_list[idx_to_pos[primary_idx]] = True
        if not has_conflict and "_sn_" not in str(key) and "decimalLatitude" in df.columns and "decimalLongitude" in df.columns:
            donors = [
                j for _, j in scored_candidates
                if has_coords[idx_to_pos[j]] and not _as_bool_flag(df.at[j, "is_centroid_or_inst_coord"] if "is_centroid_or_inst_coord" in df.columns else False)
            ]
            if donors:
                d0 = min(donors, key=lambda j: float(df.at[j, "coordinateUncertainty"]) if ("coordinateUncertainty" in df.columns and pd.notna(df.at[j, "coordinateUncertainty"])) else 99999.0)
                d_st = str(df.at[d0, "stateProvince"] if "stateProvince" in df.columns and pd.notna(df.at[d0, "stateProvince"]) else "").strip().lower()
                for idx in indices:
                    if not has_coords[idx_to_pos[idx]]:
                        r_st = str(df.at[idx, "stateProvince"] if "stateProvince" in df.columns and pd.notna(df.at[idx, "stateProvince"]) else "").strip().lower()
                        if not r_st or not d_st or r_st == d_st:
                            df.at[idx, "decimalLatitude"], df.at[idx, "decimalLongitude"] = lats[idx_to_pos[d0]], lons[idx_to_pos[d0]]
                            if "coordinateUncertainty" in df.columns and pd.isna(df.at[idx, "coordinateUncertainty"]) and pd.notna(df.at[d0, "coordinateUncertainty"]):
                                df.at[idx, "coordinateUncertainty"] = df.at[d0, "coordinateUncertainty"]
                            if "regional_group" in df.columns and (pd.isna(df.at[idx, "regional_group"]) or str(df.at[idx, "regional_group"]).strip() in {"", "NA", "Unknown"}) and pd.notna(df.at[d0, "regional_group"]):
                                df.at[idx, "regional_group"] = df.at[d0, "regional_group"]

    df["is_primary_duplicate"] = is_primary_list
    df["exsiccatae_taxonomic_conflict"] = conflict_list
    df["exsiccatae_gold_taxon"] = gold_taxon_list
    df.index = orig_index

    if logger:
        n_prim = sum(is_primary_list)
        n_dup = n_rows - n_prim
        n_conf = sum(1 for v in conflict_list if v)
        logger.info(
            f"Duplicate Stratification: Processed {n_rows} vouchers -> "
            f"{n_prim} primary duplicate vouchers designated ({n_dup} duplicate sheets flagged, "
            f"{n_conf} sheets in taxonomically conflicted duplicate groups)."
        )

    return df
