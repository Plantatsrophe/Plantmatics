"""Step 01 curated table export, summary logging, and verbatim conflict detection."""

from __future__ import annotations

import logging
import math
import os
import re
import tempfile
from pathlib import Path
from typing import TYPE_CHECKING, Any, Dict, List, Optional, Set, Tuple, Union

from scripts.core.config import DEFAULT_HOME_HERBARIUM, resolve_target_genus
from scripts.core.step01_harvester_filters import _SYNONYM_GENERA, _extract_canonical_binomial

if TYPE_CHECKING:
    import pandas as pd


def export_curated_table(
    df: pd.DataFrame,
    output_path: Path,
    logger: Optional[logging.Logger] = None,
    home_herbarium: str = DEFAULT_HOME_HERBARIUM,
    db_path: Optional[Union[str, Path]] = None,
) -> Path:
    """
    Atomically exports the 23-column Darwin Core `EXPORT_COLUMNS` voucher table to `output_path`
    (merging by `catalogNumber` on resume) while synchronizing extended diagnostic metadata
    (`infraspecific_taxon`, `georef_qc_status`, `usda_symbol`, `usda_county_status`,
    `exsiccatae_taxonomic_conflict`, `exsiccatae_gold_taxon`, `spatial_balance_weight`)
    and dynamic synonym mappings into SQLite (`vouchers`, `taxon_synonyms`).
    """
    import scripts.core.harvester as _harv

    pd_mod = _harv.pd
    np_mod = _harv.np
    export_cols = _harv.EXPORT_COLUMNS
    output_path = Path(output_path).resolve()
    output_path.parent.mkdir(parents=True, exist_ok=True)

    # Merge with existing records if table exists to prevent data loss on resumption
    if output_path.exists() and output_path.stat().st_size > 0:
        try:
            df_existing = pd_mod.read_csv(output_path)
            if not df_existing.empty and "catalogNumber" in df_existing.columns:
                if df.empty:
                    df = df_existing
                else:
                    pre_count = len(df_existing)
                    df = pd_mod.concat([df_existing, df], ignore_index=True)
                    df = df.drop_duplicates(subset=["catalogNumber"], keep="last").reset_index(drop=True)
                    if logger:
                        logger.info(
                            f"Merged records with existing table: {pre_count} previous -> "
                            f"{len(df)} total unique vouchers by catalogNumber."
                        )
        except Exception as e:
            if logger:
                logger.warning(f"Could not read existing table at {output_path} for merging: {e}")

    df = df.copy()
    if not df.empty:
        if "scientificName" not in df.columns and "species_raw" in df.columns:
            df["scientificName"] = df["species_raw"]
        elif "species_raw" not in df.columns and "scientificName" in df.columns:
            df["species_raw"] = df["scientificName"]

        if "decimalLatitude" not in df.columns:
            df["decimalLatitude"] = df["latitude"] if "latitude" in df.columns else np_mod.nan
        elif "latitude" in df.columns:
            df["decimalLatitude"] = df["decimalLatitude"].where(df["decimalLatitude"].notna(), df["latitude"])
        if "decimalLongitude" not in df.columns:
            df["decimalLongitude"] = df["longitude"] if "longitude" in df.columns else np_mod.nan
        elif "longitude" in df.columns:
            df["decimalLongitude"] = df["decimalLongitude"].where(df["decimalLongitude"].notna(), df["longitude"])
        if "identifiedBy" not in df.columns:
            df["identifiedBy"] = ""

        dwc_contract_cols = [
            "catalogNumber", "scientificName", "decimalLatitude", "decimalLongitude",
            "eventDate", "identifiedBy", "determiner_tier", "image_path",
        ]
        for c in dwc_contract_cols:
            if c not in df.columns:
                df[c] = ""

        prev_tiers = df["determiner_tier"].astype(str).tolist()
        df = _harv.apply_corpus_authority_tiers(df, db_path=db_path)
        if "_raw_date_identified" in df.columns:
            df = df.drop(columns=["_raw_date_identified"])
        tiers_changed = df["determiner_tier"].astype(str).tolist() != prev_tiers

        if (
            tiers_changed
            or "exsiccatae_key" not in df.columns
            or "is_primary_duplicate" not in df.columns
            or "exsiccatae_taxonomic_conflict" not in df.columns
        ):
            df = _harv.stratify_duplicates(df, logger=logger, home_herbarium=home_herbarium)

    _allowed_extra_voucher_cols = {
        "locality", "doy", "species_standardized", "infraspecific_taxon",
        "exsiccatae_taxonomic_conflict", "exsiccatae_gold_taxon",
        "typified_basionym", "verbatim_taxon_concept", "concept_synonymy_status",
        "identification_year", "is_open_or_hybrid_concept",
        "is_centroid_or_inst_coord", "coord_quarantine_reason",
        "pixels_per_mm", "is_calibrated_mm", "scale_source",
        "qc_status", "failure_reason", "failure_details",
        "spatial_balance_weight", "countryCode", "gbifID", "occurrenceID",
        "minimumElevationInMeters", "elevation_m",
    }
    cols_to_export = [col for col in export_cols if col in df.columns]
    extra_cols = [
        c for c in df.columns
        if c not in cols_to_export and c in _allowed_extra_voucher_cols
    ]
    final_cols = cols_to_export + extra_cols

    df_export = df[final_cols] if not df.empty else pd_mod.DataFrame(columns=export_cols)

    temp_file = tempfile.NamedTemporaryFile(
        mode="w",
        delete=False,
        dir=output_path.parent,
        suffix=".tmp",
        encoding="utf-8",
    )
    temp_path = Path(temp_file.name)
    try:
        df_export.to_csv(temp_file, index=False, encoding="utf-8")
        temp_file.flush()
        os.fsync(temp_file.fileno())
        temp_file.close()
        temp_path.replace(output_path)
        if logger:
            logger.info(f"Atomically saved curated vouchers table ({len(df_export)} records) to: {output_path}")
        if db_path is not None or _harv.is_canonical_or_explicit_db_context(output_path):
            try:
                conn = _harv.get_db_connection(db_path=db_path, init_schema=True)
                try:
                    _harv.upsert_taxon_synonyms(
                        conn,
                        _harv._DYNAMIC_SYNONYM_TO_CANONICAL,
                        _harv._DYNAMIC_TAXON_AUTHOR_STRINGS,
                    )
                    if not df_export.empty:
                        sql_extra = (
                            "exsiccatae_taxonomic_conflict", "exsiccatae_gold_taxon",
                            "usda_symbol", "usda_county_status", "infraspecific_taxon",
                            "typified_basionym", "verbatim_taxon_concept",
                            "concept_synonymy_status", "identification_year",
                            "is_open_or_hybrid_concept", "georef_qc_status",
                            "is_centroid_or_inst_coord", "coord_quarantine_reason",
                            "spatial_balance_weight",
                        )
                        sql_cols = list(final_cols) + [
                            c for c in sql_extra if c in df.columns and c not in final_cols
                        ]
                        db_records = df[sql_cols].to_dict(orient="records")
                        for rec in db_records:
                            if not rec.get("species_standardized"):
                                rec["species_standardized"] = _harv.standardize_packera_taxon(
                                    rec.get("scientificName") or rec.get("species_raw")
                                )
                        _harv.upsert_vouchers_records(conn, db_records)
                finally:
                    conn.close()
            except Exception as db_exc:
                if logger:
                    logger.warning(f"SQLite synchronization warning for {output_path}: {db_exc}")
        return output_path
    except Exception as e:
        if temp_path.exists():
            try:
                temp_path.unlink()
            except Exception:
                pass
        if logger:
            logger.error(f"Failed to atomically write curated table to {output_path}: {e}")
        raise


def print_and_log_summary(
    df: pd.DataFrame,
    download_stats: Optional[Dict[str, Any]],
    logger: logging.Logger,
    usda_filter_stats: Optional[Dict[str, Any]] = None,
) -> None:
    """Generates a publication-grade summary of the harvested dataset."""
    total_records = len(df)
    logger.info("=" * 80)
    logger.info("                  PACKERA VOUCHER INGESTION & CURATION SUMMARY                  ")
    logger.info("=" * 80)
    logger.info(f"Total Quality-Filtered Specimen Vouchers: {total_records:,}")

    if total_records == 0:
        logger.warning("No records were retained. Check query parameters or network connection.")
        logger.info("=" * 80)
        return

    # Determiner Tier Breakdown
    tier_counts = df["determiner_tier"].value_counts()
    logger.info("\n--- TAXONOMIC DETERMINER AUTHORITY STRATIFICATION ---")
    for tier in ["Tier_1_Gold", "Tier_2_Silver", "Tier_3_Bronze"]:
        cnt = tier_counts.get(tier, 0)
        pct = (cnt / total_records) * 100.0
        logger.info(f"  * {tier:<15} : {cnt:>5} records ({pct:>5.1f}%)")

    # Species Breakdown
    species_counts = df["species_raw"].value_counts()
    logger.info("\n--- TAXON DISTRIBUTION (RAW DETERMINATIONS) ---")
    for sp, cnt in species_counts.head(8).items():
        pct = (cnt / total_records) * 100.0
        logger.info(f"  * {sp:<45} : {cnt:>5} records ({pct:>5.1f}%)")

    # Regional Ecological Groups Breakdown (EPA Level III Ecoregions)
    region_counts = df["regional_group"].fillna("NA").replace("", "NA").value_counts()
    logger.info("\n--- US EPA LEVEL III ECOREGIONS ---")
    for reg, cnt in region_counts.items():
        pct = (cnt / total_records) * 100.0
        logger.info(f"  * {str(reg):<30} : {cnt:>5} records ({pct:>5.1f}%)")

    # Herbarium Institutions (Top 10)
    inst_counts = df["institutionCode"].value_counts()
    logger.info("\n--- TOP HERBARIUM INSTITUTIONS ---")
    for inst, cnt in inst_counts.head(10).items():
        pct = (cnt / total_records) * 100.0
        logger.info(f"  * {inst:<15} : {cnt:>5} records ({pct:>5.1f}%)")

    # Exsiccatae & Duplicate Stratification Breakdown
    if "is_primary_duplicate" in df.columns:
        n_primary = int(df["is_primary_duplicate"].sum())
        n_duplicates = total_records - n_primary
        logger.info("\n--- EXSICCAE & DUPLICATE STRATIFICATION ---")
        logger.info(f"  * Primary Duplicate Vouchers  : {n_primary:>5} records ({(n_primary / total_records) * 100.0:>5.1f}%)")
        logger.info(f"  * Secondary Duplicate Sheets  : {n_duplicates:>5} records ({(n_duplicates / total_records) * 100.0:>5.1f}%)")
        if "duplicate_count" in df.columns:
            multi_sheets = int((df["duplicate_count"] > 1).sum())
            logger.info(f"  * Sheets from Multi-Sheet Events: {multi_sheets:>5} records")

    # USDA PLANTS County Range Filter Summary
    if usda_filter_stats and usda_filter_stats.get("evaluated", 0) > 0:
        logger.info("\n--- USDA PLANTS COUNTY RANGE FILTER ---")
        logger.info(f"  * Vouchers Evaluated          : {usda_filter_stats.get('evaluated', 0):>5}")
        logger.info(f"  * Documented In-Range         : {usda_filter_stats.get('in_range', 0):>5}")
        logger.info(f"  * Buffered / 1-Ring Adjacent  : {usda_filter_stats.get('buffered_adjacent', 0):>5}")
        logger.info(f"  * Exempt Retained (Type/Gold) : {usda_filter_stats.get('exempt_retained', 0):>5}")
        logger.info(f"  * Unmapped Taxon / Non-US     : {usda_filter_stats.get('unmapped_or_non_us', 0):>5}")
        logger.info(f"  * Dropped Out-of-Range        : {usda_filter_stats.get('dropped_out_of_range', 0):>5}")

    # Image Download & Quality Summary
    if download_stats:
        logger.info("\n--- SPECIMEN IMAGE DOWNLOAD & QUALITY STATUS ---")
        logger.info(f"  * Downloaded Successfully : {download_stats.get('success', 0):>5}")
        logger.info(f"  * Cached / Skipped        : {download_stats.get('skipped', 0):>5}")
        logger.info(f"  * Quality Filter Rejected : {download_stats.get('quality_rejected', 0):>5}")
        logger.info(f"  * Failed / Inaccessible   : {download_stats.get('failed', 0):>5}")
        if "median_mp" in download_stats:
            logger.info(f"  * Median Image Resolution : {download_stats.get('median_mp', 0.0):>5.2f} Megapixels")

    logger.info("=" * 80)


_CONFLICT_STOP_TOKENS: Set[str] = {
    "det", "conf", "rev", "annot", "by", "the", "and", "jr", "sr",
    "van", "von", "der", "den", "ter", "bar", "prof", "dr", "mr",
    "ms", "mrs", "herb", "univ", "bot", "ny", "ph", "ncu",
}
_NEGATED_TAXON_RE = re.compile(r"\b(?:not|non|nec|formerly|misapplied|excl\.?)\s+$", flags=re.IGNORECASE)
_MENTION_RE_CACHE: Dict[Tuple[str, ...], re.Pattern[str]] = {}


def _detect_verbatim_identification_conflict(
    rec: Optional[Dict[str, Any]],
    resolved_canonical: str,
    determiner_tier: str,
) -> Tuple[bool, str]:
    """
    Detect portal-lag identification conflicts where a voucher's structured `identifications`,
    `identificationRemarks` / `verbatimIdentificationHistory` (linked to `identifiedBy`), or
    `verbatimIdentification` points to a different canonical species than GBIF `scientificName`,
    while ignoring superseded `previousIdentifications` and `_DYNAMIC_SYNONYM_TO_CANONICAL` synonyms.
    """
    import scripts.core.harvester as _harv

    if not isinstance(rec, dict) or not resolved_canonical:
        return False, ""

    target_genus = (resolve_target_genus() or "Packera").strip().capitalize()
    allowed_genera = {target_genus, "Packera"} | set(_SYNONYM_GENERA)
    cache_key = tuple(sorted(allowed_genera))
    mention_re = _MENTION_RE_CACHE.get(cache_key)
    if mention_re is None:
        genus_alt = "|".join(sorted((re.escape(g) for g in allowed_genera), key=len, reverse=True))
        mention_re = re.compile(
            rf"\b({genus_alt})(?:\s*×\s*|\s+[x×]\s+|\s+)([a-z][a-z-]{{1,}})"
            r"(?:\s+(var|subsp|ssp|f)\.?\s+([a-z][a-z-]{1,}))?\b",
            flags=re.IGNORECASE,
        )
        _MENTION_RE_CACHE[cache_key] = mention_re

    def _norm_tx(raw_val: Any) -> Optional[str]:
        if raw_val is None or (isinstance(raw_val, float) and math.isnan(raw_val)):
            return None
        txt = re.sub(r"\s+", " ", str(raw_val).strip())
        if not txt:
            return None
        low = txt.lower()
        dyn_map = _harv._DYNAMIC_SYNONYM_TO_CANONICAL
        if low in dyn_map:
            return dyn_map[low]
        var_m = re.match(r"^([a-z]+\s+[a-z-]+\s+(?:var\.?|subsp\.?|ssp\.?|f\.?)\s+[a-z-]+)", low)
        if var_m:
            v_raw = var_m.group(1)
            v_dot = re.sub(r"\b(var|subsp|ssp|f)\b\.?", r"\1.", v_raw)
            if v_raw in dyn_map:
                return dyn_map[v_raw]
            if v_dot in dyn_map:
                return dyn_map[v_dot]
        cand_bin = _extract_canonical_binomial(txt, allowed_genera)
        if not cand_bin:
            return None
        if cand_bin.lower() in dyn_map:
            return dyn_map[cand_bin.lower()]
        if cand_bin.startswith(f"{target_genus} ") or cand_bin.startswith("Packera "):
            return cand_bin
        if _harv._is_gbif_species_mocked():
            std_m = _harv.standardize_packera_taxon(cand_bin)
            if std_m:
                return std_m
        return None

    def _extract_mentions(raw_text: Any) -> List[str]:
        if raw_text is None:
            return []
        txt = str(raw_text).strip()
        if not txt:
            return []
        found: List[str] = []
        for m in mention_re.finditer(txt):
            prefix = txt[max(0, m.start() - 28):m.start()]
            if _NEGATED_TAXON_RE.search(prefix):
                continue
            cand = m.group(0)
            tx = _norm_tx(cand) or _norm_tx(f"{m.group(1)} {m.group(2)}")
            if tx and tx not in found:
                found.append(tx)
        return found

    resolved_norm = _norm_tx(resolved_canonical) or _harv.standardize_packera_taxon(resolved_canonical) or ""
    if not resolved_norm:
        return False, ""

    is_gold = str(determiner_tier or "").strip() == "Tier_1_Gold"
    superseded_taxa: Set[str] = set()
    prev_raw = rec.get("previousIdentifications")
    if isinstance(prev_raw, list):
        for pv in prev_raw:
            superseded_taxa.update(_extract_mentions(pv))
    elif prev_raw:
        superseded_taxa.update(_extract_mentions(prev_raw))

    # 1. Structured identifications array: evaluate highest-priority current or latest dated determination
    idents = rec.get("identifications")
    if isinstance(idents, list) and idents:
        ranked_idents: List[Tuple[int, str, int, Dict[str, Any], Optional[str]]] = []
        for i_pos, entry in enumerate(idents):
            if not isinstance(entry, dict):
                continue
            curr_flag = 1 if entry.get("isCurrent") is True else (-1 if entry.get("isCurrent") is False else 0)
            date_s = str(entry.get("dateIdentified") or "").strip()
            raw_id_nm = entry.get("scientificName") or entry.get("taxonName") or entry.get("verbatimIdentification") or ""
            tx_id = _norm_tx(raw_id_nm) or (next(iter(_extract_mentions(raw_id_nm)), None))
            ranked_idents.append((curr_flag, date_s, i_pos, entry, tx_id))
        if ranked_idents:
            ranked_idents.sort(key=lambda t: (t[0], t[1], t[2]), reverse=True)
            top_curr_flag, top_date, _, top_entry, top_tx = ranked_idents[0]
            for cf, _, _, _, tx_other in ranked_idents[1:]:
                if tx_other and tx_other != resolved_norm and (cf < 0 or top_tx == resolved_norm):
                    superseded_taxa.add(tx_other)
            if top_tx and top_tx != resolved_norm and (top_curr_flag >= 0 or not any(r[4] == resolved_norm for r in ranked_idents)):
                top_det_tier, _, _ = _harv.parse_determiner_tier(
                    None,
                    top_entry.get("identifiedBy"),
                    rec.get("recordedBy"),
                    top_entry.get("identificationRemarks"),
                    rec.get("institutionCode"),
                    rec.get("locality"),
                    rec.get("habitat"),
                )
                gold_tx = top_tx if (is_gold or top_det_tier == "Tier_1_Gold") else ""
                return True, gold_tx

    # Extract surname tokens from identifiedBy for linking remarks/history to active determiner
    det_raw = str(rec.get("identifiedBy") or "").strip()
    det_surnames: List[str] = []
    if det_raw:
        for tok in re.findall(r"[A-Za-z][A-Za-z-]{1,}", det_raw):
            tl = tok.lower()
            if tl not in _CONFLICT_STOP_TOKENS and len(tl) >= 2:
                det_surnames.append(tl)

    # 2. identifiedBy-linked identificationRemarks / verbatimIdentificationHistory clauses
    if det_surnames:
        remark_clauses: List[str] = []
        for fld in ("verbatimIdentificationHistory", "identificationRemarks"):
            val = rec.get(fld)
            if val:
                remark_clauses.extend([c.strip() for c in re.split(r"[;|\n]+", str(val)) if c.strip()])
        linked_matches: List[str] = []
        for cl in remark_clauses:
            if any(re.search(rf"\b{re.escape(sn)}\b", cl, flags=re.IGNORECASE) for sn in det_surnames):
                cl_taxa = _extract_mentions(cl)
                if cl_taxa:
                    linked_matches.extend(cl_taxa)
        if linked_matches:
            latest_linked = linked_matches[-1]
            if latest_linked != resolved_norm:
                return True, (latest_linked if is_gold else "")

    # 3. Direct verbatimIdentification differing from resolved canonical (when not merely superseded history)
    v_raw = str(rec.get("verbatimIdentification") or "").strip()
    if v_raw:
        v_tx = _norm_tx(v_raw) or next(iter(_extract_mentions(v_raw)), None)
        if v_tx and v_tx != resolved_norm and v_tx not in superseded_taxa:
            return True, (v_tx if is_gold else "")

    return False, ""
