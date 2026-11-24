"""step01_exsiccatae_reconcile.py - Cross-taxon SQLite exsiccatae conflict and Gold-consensus reconciliation."""

from __future__ import annotations

from typing import Any, Dict, List, Optional

try:
    import pandas as pd
except ImportError:
    pd = None  # type: ignore[assignment]

from scripts.core.step01_exsiccatae_fingerprint import (
    _as_bool_flag,
    _is_extralimital_status,
    _resolve_row_standardized_packera,
)


def reconcile_cross_taxon_exsiccatae_with_db(
    df: pd.DataFrame,
    conn: Any,
    standardize_fn: Optional[Any] = None,
) -> pd.DataFrame:
    """
    Cross-reference `exsiccatae_key` across the full SQLite `vouchers` table so duplicate
    sheets filed under separate online species still detect taxonomic disagreement and
    in-range specialist annotations.
    """
    if pd is None or df is None or df.empty or conn is None or "exsiccatae_key" not in df.columns:
        return df

    try:
        has_tbl = conn.execute(
            "SELECT 1 FROM sqlite_master WHERE type='table' AND name='vouchers' LIMIT 1"
        ).fetchone()
    except Exception:
        return df
    if not has_tbl:
        return df

    raw_keys = [
        str(v).strip()
        for v in df["exsiccatae_key"].dropna().unique().tolist()
        if str(v).strip() and str(v).strip().lower() not in {"none", "nan"}
    ]
    if not raw_keys:
        return df

    try:
        v_cols = {r[1] for r in conn.execute("PRAGMA table_info(vouchers)").fetchall()}
    except Exception:
        return df

    want_cols = [
        c
        for c in (
            "catalogNumber",
            "exsiccatae_key",
            "species_standardized",
            "scientificName",
            "species_raw",
            "determiner_tier",
            "usda_county_status",
        )
        if c in v_cols
    ]
    if "exsiccatae_key" not in want_cols:
        return df

    db_rows: List[Dict[str, Any]] = []
    sel_sql = ", ".join(want_cols)
    for start in range(0, len(raw_keys), 500):
        chunk = raw_keys[start : start + 500]
        ph = ", ".join(["?"] * len(chunk))
        try:
            cur = conn.execute(
                f"SELECT {sel_sql} FROM vouchers WHERE exsiccatae_key IN ({ph})",
                chunk,
            )
            for row_t in cur.fetchall():
                db_rows.append(dict(zip(want_cols, row_t)))
        except Exception:
            continue

    if not db_rows:
        return df

    def _std_tx(rec: Any) -> str:
        raw_val = ""
        for fld in ("species_standardized", "scientificName", "species_raw"):
            v = rec.get(fld) if isinstance(rec, dict) else rec.get(fld, None)
            if v is not None and not (pd is not None and pd.isna(v)) and str(v).strip():
                raw_val = str(v).strip()
                break
        if not raw_val:
            return ""
        if standardize_fn is not None:
            try:
                res = standardize_fn(raw_val)
                if res:
                    return str(res).strip()
            except Exception:
                pass
        return _resolve_row_standardized_packera(rec if isinstance(rec, dict) else rec)

    # Index DB records grouped by exsiccatae_key
    db_by_key: Dict[str, List[Dict[str, Any]]] = {}
    for drec in db_rows:
        ek = str(drec.get("exsiccatae_key") or "").strip()
        if ek:
            db_by_key.setdefault(ek, []).append(drec)

    out_df = df.copy()
    if "exsiccatae_taxonomic_conflict" not in out_df.columns:
        out_df["exsiccatae_taxonomic_conflict"] = False
    if "exsiccatae_gold_taxon" not in out_df.columns:
        out_df["exsiccatae_gold_taxon"] = ""

    conf_vals = [_as_bool_flag(v) for v in out_df["exsiccatae_taxonomic_conflict"].tolist()]
    gold_vals = [
        "" if (v is None or pd.isna(v)) else str(v).strip()
        for v in out_df["exsiccatae_gold_taxon"].tolist()
    ]

    for key, idx_lbls in out_df.groupby("exsiccatae_key", sort=False).groups.items():
        ek = str(key or "").strip()
        if not ek or ek not in db_by_key:
            continue
        batch_indices = list(idx_lbls)
        batch_cats = {
            str(out_df.at[i_l, "catalogNumber"]).strip()
            for i_l in batch_indices
            if "catalogNumber" in out_df.columns and not pd.isna(out_df.at[i_l, "catalogNumber"])
        }
        siblings: List[Dict[str, Any]] = []
        for i_l in batch_indices:
            r_s = out_df.loc[i_l]
            siblings.append({
                "catalogNumber": str(r_s.get("catalogNumber") or ""),
                "tx": _std_tx(r_s),
                "tier": str(r_s.get("determiner_tier") or "").strip(),
                "extra": _is_extralimital_status(r_s),
                "pre_conf": _as_bool_flag(r_s.get("exsiccatae_taxonomic_conflict")),
            })
        for drec in db_by_key[ek]:
            d_cat = str(drec.get("catalogNumber") or "").strip()
            d_tx = _std_tx(drec)
            # Include DB siblings under other catalogNumbers (or differing taxon snapshots)
            if d_cat and d_cat in batch_cats:
                same_cat_tx = [s["tx"] for s in siblings if s["catalogNumber"].strip() == d_cat]
                if any(t == d_tx for t in same_cat_tx):
                    continue
            siblings.append({
                "catalogNumber": d_cat,
                "tx": d_tx,
                "tier": str(drec.get("determiner_tier") or "").strip(),
                "extra": _is_extralimital_status(drec),
                "pre_conf": False,
            })

        distinct_taxa = {s["tx"] for s in siblings if s["tx"]}
        if len(distinct_taxa) < 2:
            continue

        inrange_taxa = {s["tx"] for s in siblings if s["tx"] and not s["extra"]}
        inrange_gold_counts: Dict[str, int] = {}
        all_gold_counts: Dict[str, int] = {}
        for s in siblings:
            if s["tier"] == "Tier_1_Gold" and s["tx"]:
                all_gold_counts[s["tx"]] = all_gold_counts.get(s["tx"], 0) + 1
                if not s["extra"] and not s["pre_conf"]:
                    inrange_gold_counts[s["tx"]] = inrange_gold_counts.get(s["tx"], 0) + 1

        active_gold = inrange_gold_counts if inrange_gold_counts else all_gold_counts
        best_gold_tx = ""
        best_gold_from_inrange = False
        multi_gold = len(active_gold) >= 2
        if active_gold:
            best_gold_from_inrange = bool(inrange_gold_counts)
            best_gold_tx = sorted(active_gold.items(), key=lambda kv: (kv[1], kv[0]), reverse=True)[0][0]
        unanimous_inrange_tx = next(iter(inrange_taxa)) if len(inrange_taxa) == 1 else ""

        pos_map = {idx_lbl: pos for pos, idx_lbl in enumerate(out_df.index)}
        for i_l in batch_indices:
            pos_i = pos_map[i_l]
            r_s = out_df.loc[i_l]
            r_tx = _std_tx(r_s)
            r_extra = _is_extralimital_status(r_s)
            if best_gold_tx and not multi_gold:
                row_conf = bool(r_tx and r_tx != best_gold_tx)
            elif unanimous_inrange_tx:
                row_conf = bool(r_extra and r_tx and r_tx != unanimous_inrange_tx)
            else:
                row_conf = True
            conf_vals[pos_i] = bool(conf_vals[pos_i] or row_conf)
            if best_gold_tx and not gold_vals[pos_i] and (not multi_gold or best_gold_from_inrange or r_extra):
                gold_vals[pos_i] = best_gold_tx

    out_df["exsiccatae_taxonomic_conflict"] = conf_vals
    out_df["exsiccatae_gold_taxon"] = gold_vals
    return out_df
