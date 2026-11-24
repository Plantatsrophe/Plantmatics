"""step01_eastern_cohort_scope.py - Geographic scope resolution, cis-Mississippi checks, and Eastern NA cohort gating."""

from __future__ import annotations

import math
import re
from typing import Any, Dict, List, Mapping, Optional, Set, Tuple, Union

try:
    import pandas as pd
except ImportError:
    pd = None  # type: ignore[assignment]

from scripts.core.config import (
    CIS_MISSISSIPPI_STATES_PROVINCES,
    DEFAULT_CIS_MISSISSIPPI_LON,
    DEFAULT_COUNTRIES,
    DEFAULT_GEOGRAPHIC_SCOPE,
    DEFAULT_MIN_CIS_MISSISSIPPI_FRAC,
    EXCLUDED_WESTERN_STATES,
    VALID_GEOGRAPHIC_SCOPES,
    WESTERN_LONGITUDE_THRESHOLD,
    normalize_countries_list,
)
from scripts.core.step01_determiner_surnames import _is_verified_type_status
from scripts.core.step01_outgroup_anchors import select_cascading_outgroup_anchors


def resolve_geographic_scope_settings(
    geographic_scope: Optional[str] = DEFAULT_GEOGRAPHIC_SCOPE,
    exclude_western: bool = True,
    include_western: bool = False,
    countries: Optional[Union[str, List[str]]] = None,
) -> Tuple[str, bool, List[str]]:
    """
    Resolves effective `(geographic_scope, exclude_western, countries)` settings:
      - `eastern_na` (default): `exclude_western=True` (unless `include_western=True` or
        `exclude_western=False`), `countries=["US", "CA"]` from `DEFAULT_COUNTRIES`
        (or `["US", "CA", "MX"]` when `include_western=True`).
      - `north_america`: automatically sets `exclude_western=False` and expands default
        `countries` to `["US", "CA", "MX"]` unless explicit `countries` were supplied.
      - `global`: automatically sets `exclude_western=False` and sets `countries=[]` (unrestricted
        global GBIF search) unless explicit `countries` were supplied.
    """
    scope = str(geographic_scope or DEFAULT_GEOGRAPHIC_SCOPE).strip().lower()
    if scope not in VALID_GEOGRAPHIC_SCOPES:
        raise ValueError(f"Unsupported geographic_scope={geographic_scope!r}; expected one of {sorted(VALID_GEOGRAPHIC_SCOPES)}.")

    eff_exclude_western = bool(exclude_western) and (not bool(include_western)) and (scope == "eastern_na")

    if countries is not None:
        eff_countries = normalize_countries_list(countries, default=[])
    elif scope == "global":
        eff_countries = []
    elif scope == "north_america" or bool(include_western):
        eff_countries = ["US", "CA", "MX"]
    else:
        eff_countries = normalize_countries_list(DEFAULT_COUNTRIES, default=["US", "CA"])

    return scope, eff_exclude_western, eff_countries


def is_excluded_western_region(
    state_province: Optional[str],
    lat: Optional[float] = None,
    lon: Optional[float] = None,
    western_longitude_threshold: float = WESTERN_LONGITUDE_THRESHOLD,
    exclude_western: bool = True,
    is_type_specimen: bool = False,
) -> bool:
    """
    Determines whether a record originates from a western state/province (`EXCLUDED_WESTERN_STATES`)
    or west of `western_longitude_threshold` (default `-104.05`), while exempting Eastern Canadian
    provinces, Washington D.C., and verified nomenclatural type specimens (`is_type_specimen=True`).
    """
    if not exclude_western:
        return False

    has_valid_state = False
    if state_province is not None and not (isinstance(state_province, float) and math.isnan(state_province)):
        if pd is None or not pd.isna(state_province):
            raw_state = str(state_province).strip()
            cleaned = re.sub(r"\((?:state|province|territory)\)", "", raw_state, flags=re.IGNORECASE).strip(" ._,-")
            upper_state = cleaned.upper()
            if upper_state and upper_state not in {"NONE", "NAN", "NA", "N/A", "NULL", "UNKNOWN", "UNSPECIFIED", "UNRECORDED"}:
                has_valid_state = True
                if upper_state in {
                    "WASHINGTON, D.C.", "WASHINGTON D.C.", "WASHINGTON DC", "DISTRICT OF COLUMBIA", "DC",
                    "ON", "ONTARIO", "QC", "QUEBEC", "QUÉBEC", "NB", "NEW BRUNSWICK",
                    "NS", "NOVA SCOTIA", "PE", "PEI", "PRINCE EDWARD ISLAND",
                    "NL", "NEWFOUNDLAND", "NEWFOUNDLAND AND LABRADOR",
                }:
                    return False

                if upper_state in EXCLUDED_WESTERN_STATES:
                    return True

    has_valid_lon = False
    if lon is not None:
        try:
            lon_val = float(lon)
            if math.isfinite(lon_val):
                has_valid_lon = True
                if lon_val < float(western_longitude_threshold):
                    return True
        except (ValueError, TypeError):
            pass

    if not has_valid_lon and not has_valid_state and not is_type_specimen:
        return True

    return False


def _clean_upper_state_province(state_province: Optional[Any]) -> str:
    if state_province is None or (isinstance(state_province, float) and math.isnan(state_province)):
        return ""
    if pd is not None and pd.isna(state_province):
        return ""
    raw_state = str(state_province).strip()
    cleaned = re.sub(r"\((?:state|province|territory)\)", "", raw_state, flags=re.IGNORECASE).strip(" ._,-")
    up = cleaned.upper()
    if up in {"NONE", "NAN", "NA", "N/A", "NULL", "UNKNOWN", "UNSPECIFIED", "UNRECORDED"}:
        return ""
    return up


def is_cis_mississippi_record(
    state_province: Optional[str],
    lat: Optional[float] = None,
    lon: Optional[float] = None,
    cis_longitude_threshold: float = DEFAULT_CIS_MISSISSIPPI_LON,
) -> bool:
    """
    Returns True strictly when an occurrence lies east of the Mississippi River
    (`CIS_MISSISSIPPI_STATES_PROVINCES`, or `lon >= cis_longitude_threshold` when state is missing).
    Trans-Mississippi states (MN, IA, MO, AR, LA, ND, SD, NE, KS, OK, TX, MB) return False so
    Cordilleran and Great Plains species cannot pass Eastern NA species-residency discovery.
    """
    up_state = _clean_upper_state_province(state_province)
    lon_val: Optional[float] = None
    if lon is not None:
        try:
            cand_lon = float(lon)
            if math.isfinite(cand_lon):
                lon_val = cand_lon
        except (ValueError, TypeError):
            lon_val = None

    if up_state:
        if up_state in CIS_MISSISSIPPI_STATES_PROVINCES:
            # Guard against grossly transposed western coordinates inside an eastern state label
            if lon_val is not None and lon_val < -92.5 and up_state not in {"ONTARIO", "ON"}:
                return False
            return True
        return False

    if lon_val is not None and float(cis_longitude_threshold) <= lon_val <= -50.0:
        if lat is not None:
            try:
                lat_f = float(lat)
                if math.isfinite(lat_f) and lat_f < 24.0:
                    return False
            except (ValueError, TypeError):
                pass
        return True
    return False


def evaluate_eastern_na_cohort_filter(
    df: "pd.DataFrame",
    standardize_taxon_fn: Any,
    *,
    western_longitude_threshold: float = WESTERN_LONGITUDE_THRESHOLD,
    transcontinental_longitude_threshold: Optional[float] = None,
    cis_longitude_threshold: float = DEFAULT_CIS_MISSISSIPPI_LON,
    min_cis_records: int = 10,
    min_cis_fraction: float = DEFAULT_MIN_CIS_MISSISSIPPI_FRAC,
    strict_species_whitelist: bool = False,
) -> Tuple["pd.Series", Dict[str, Any]]:
    """
    Two-stage Cis-Mississippi Eastern North America + Ecoregion-Bounded cohort filter:
      1. Species Gate: Retains genuinely cis-Mississippi Eastern Packera species from
         `resolve_eastern_resident_whitelist()` (YAML + SQLite TDWG/NatureServe) or
         meeting `cis_count >= eff_min_cis` and `cis_fraction >= min_cis_fraction`.
      2. Voucher Gate: Retains cis-Mississippi vouchers plus trans-Mississippi disjuncts
         up to `western_longitude_threshold` (default `-100.0W`), screens transcontinental
         boreal taxa by EPA/CEC Level III ecoregion via `evaluate_transcontinental_ecoregion_gate()`,
         and rescues verified nomenclatural types of retained Eastern species.
    """
    if df is None or df.empty:
        empty_mask = pd.Series(dtype=bool, index=df.index if df is not None else None)
        return empty_mask, {
            "retained_species": [],
            "dropped_species": [],
            "species_breakdown": {},
            "biogeographic_boundary_status": {},
        }

    from scripts.core.step01_biogeography_utils import (
        _canonical_taxon_equiv_key,
        evaluate_transcontinental_ecoregion_gate,
        resolve_eastern_resident_whitelist,
        resolve_transcontinental_taxa,
    )
    _ = transcontinental_longitude_threshold
    whitelist_keys = {_canonical_taxon_equiv_key(s): s for s in resolve_eastern_resident_whitelist()}
    transcontinental_keys = {_canonical_taxon_equiv_key(s) for s in resolve_transcontinental_taxa()}

    def _resolve_sp(row: Mapping[str, Any]) -> str:
        for col in ("species_standardized", "scientificName", "species_raw"):
            val = row.get(col)
            if val is not None and str(val).strip() and str(val).strip().lower() not in {"nan", "none"}:
                std = standardize_taxon_fn(val) if callable(standardize_taxon_fn) else None
                if std:
                    return str(std).strip()
        raw_fb = str(row.get("scientificName") or row.get("species_raw") or "").strip()
        return raw_fb

    lat_col = "decimalLatitude" if "decimalLatitude" in df.columns else ("latitude" if "latitude" in df.columns else None)
    lon_col = "decimalLongitude" if "decimalLongitude" in df.columns else ("longitude" if "longitude" in df.columns else None)

    n_rows = len(df)
    sp_list: List[str] = []
    is_type_list: List[bool] = []
    is_cis_list: List[bool] = []
    is_west_excl_list: List[bool] = []
    is_transcont_west_list: List[bool] = []
    bio_status_list: List[str] = []

    records = df.to_dict(orient="records")
    for rec in records:
        sp = _resolve_sp(rec)
        sp_list.append(sp)
        if "_is_type_rec" in rec and rec["_is_type_rec"] is not None and not pd.isna(rec["_is_type_rec"]):
            is_t = bool(rec["_is_type_rec"])
        else:
            is_t, _ = _is_verified_type_status(rec.get("type_status") or rec.get("_raw_type_status"))
        is_type_list.append(is_t)

        st_val = rec.get("stateProvince")
        lat_val = rec.get(lat_col) if lat_col else None
        lon_val = rec.get(lon_col) if lon_col else None

        cis_ok = is_cis_mississippi_record(
            st_val, lat=lat_val, lon=lon_val, cis_longitude_threshold=cis_longitude_threshold
        )
        is_cis_list.append(cis_ok)

        west_excl = is_excluded_western_region(
            st_val,
            lat=lat_val,
            lon=lon_val,
            western_longitude_threshold=western_longitude_threshold,
            exclude_western=True,
            is_type_specimen=is_t,
        )
        is_west_excl_list.append(west_excl)

        tc_west = False
        sp_key = _canonical_taxon_equiv_key(sp)
        if not is_t and sp_key in transcontinental_keys:
            tc_west, bio_status = evaluate_transcontinental_ecoregion_gate(
                st_val,
                lat=lat_val,
                lon=lon_val,
                regional_group=rec.get("regional_group"),
                is_type_specimen=is_t,
                western_longitude_threshold=western_longitude_threshold,
            )
        elif is_t:
            bio_status = "type_exempt"
        elif west_excl:
            bio_status = "excluded_plains_cordilleran"
        elif cis_ok:
            bio_status = "core_eastern"
        else:
            bio_status = "disjunct_western_extension"
        is_transcont_west_list.append(tc_west)
        bio_status_list.append(bio_status)

    eff_min_cis = int(min_cis_records) if n_rows >= 50 else min(3, max(1, int(min_cis_records)))

    sp_totals: Dict[str, int] = {}
    sp_non_type_totals: Dict[str, int] = {}
    sp_cis_counts: Dict[str, int] = {}
    for sp, is_t, is_cis in zip(sp_list, is_type_list, is_cis_list):
        if not sp:
            continue
        sp_totals[sp] = sp_totals.get(sp, 0) + 1
        if not is_t:
            sp_non_type_totals[sp] = sp_non_type_totals.get(sp, 0) + 1
            if is_cis:
                sp_cis_counts[sp] = sp_cis_counts.get(sp, 0) + 1

    retained_species: Set[str] = set()
    for sp, tot_n in sp_totals.items():
        cis_n = sp_cis_counts.get(sp, 0)
        nt_n = max(1, sp_non_type_totals.get(sp, 0))
        cis_frac = float(cis_n) / float(nt_n)
        sp_key = _canonical_taxon_equiv_key(sp)
        in_whitelist = sp_key in whitelist_keys
        if strict_species_whitelist:
            if in_whitelist:
                retained_species.add(sp)
        else:
            if in_whitelist and (cis_n >= 1 or tot_n == sp_totals[sp] - sp_non_type_totals.get(sp, 0) or n_rows < 50):
                if cis_n >= 1 or sp_non_type_totals.get(sp, 0) == 0:
                    retained_species.add(sp)
                elif cis_frac >= float(min_cis_fraction):
                    retained_species.add(sp)
                else:
                    if sp_key not in transcontinental_keys:
                        retained_species.add(sp)
            elif cis_n >= eff_min_cis and cis_frac >= float(min_cis_fraction):
                retained_species.add(sp)

    keep_flags: List[bool] = []
    species_breakdown: Dict[str, Dict[str, Any]] = {}
    for sp in sorted(sp_totals.keys()):
        species_breakdown[sp] = {
            "total": sp_totals[sp],
            "cis_mississippi": sp_cis_counts.get(sp, 0),
            "retained": 0,
            "disjunct_retained": 0,
            "dropped": 0,
            "resident_accepted": sp in retained_species,
        }

    bio_status_counts: Dict[str, int] = {}
    for sp, is_t, is_cis, west_excl, tc_west, b_st in zip(
        sp_list, is_type_list, is_cis_list, is_west_excl_list, is_transcont_west_list, bio_status_list
    ):
        if not sp or sp not in retained_species:
            keep_flags.append(False)
            bio_status_counts["excluded_non_resident_taxon"] = bio_status_counts.get("excluded_non_resident_taxon", 0) + 1
            if sp in species_breakdown:
                species_breakdown[sp]["dropped"] += 1
            continue
        if is_t:
            keep_flags.append(True)
            bio_status_counts["type_exempt"] = bio_status_counts.get("type_exempt", 0) + 1
            species_breakdown[sp]["retained"] += 1
            continue
        if west_excl or tc_west:
            keep_flags.append(False)
            bio_status_counts["excluded_plains_cordilleran"] = bio_status_counts.get("excluded_plains_cordilleran", 0) + 1
            species_breakdown[sp]["dropped"] += 1
            continue
        keep_flags.append(True)
        bio_status_counts[b_st] = bio_status_counts.get(b_st, 0) + 1
        species_breakdown[sp]["retained"] += 1
        if not is_cis:
            species_breakdown[sp]["disjunct_retained"] += 1

    keep_series = pd.Series(keep_flags, index=df.index, dtype=bool)
    dropped_species = sorted(set(sp_totals.keys()) - retained_species)
    stats = {
        "initial_rows": n_rows,
        "retained_rows": int(keep_series.sum()),
        "dropped_rows": int(n_rows - int(keep_series.sum())),
        "retained_species": sorted(retained_species),
        "dropped_species": dropped_species,
        "species_breakdown": species_breakdown,
        "biogeographic_boundary_status": bio_status_counts,
    }
    return keep_series, stats
