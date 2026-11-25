"""Step 01 single-record GBIF occurrence ingestion Gates 0 through 6."""

from __future__ import annotations

import datetime
import math
import re
from typing import Any, Dict, List, Optional, Set, Tuple

from scripts.core.epa_ecoregions import _US_STATE_ABBR_TO_FULL, _normalize_us_state_name
from scripts.core.step01_harvester_filters import (
    _FATAL_GBIF_GEO_ISSUES,
    _NA_STATE_PROVINCE_BBOXES,
    _SYNONYM_GENERA,
    _extract_canonical_binomial,
    _extract_infraspecific_taxon,
    _is_cultivated_occurrence,
    _normalize_ca_province_bbox_key,
)


def ingest_gbif_occurrence_record(
    harvester_self: Any,
    rec: Dict[str, Any],
    taxon: str,
    is_direct_genus_query: bool,
    current_year: int,
    us_spatial_graph: Dict[str, Any],
    all_curated_records: List[Dict[str, Any]],
    seen_catalog_keys: Set[Tuple[str, str]],
    catalog_to_inst: Dict[str, str],
    counters: Dict[str, int],
    require_type: bool = False,
    seen_gbif_keys: Optional[Set[str]] = None,
) -> None:
    """Validate, standardize, score, and append a single GBIF occurrence record."""
    import scripts.core.harvester as _harv

    rec_type_rank, _ = _harv.classify_nomenclatural_type_rank(rec.get("typeStatus"))
    is_rec_type = rec_type_rank >= 1
    if require_type and not is_rec_type:
        return

    counters["harvested"] = counters.get("harvested", 0) + 1
    early_gbif_key = str(rec.get("key") or "").strip()
    if early_gbif_key and seen_gbif_keys is not None and early_gbif_key in seen_gbif_keys:
        return

    # 0. Cultivated / Captive / Greenhouse Filter (exempting nomenclatural types)
    if not is_rec_type and _is_cultivated_occurrence(rec):
        return

    # 1. Geographic Coordinate, CoordinateCleaner Hygiene & Uncertainty Validation
    raw_issues = rec.get("issues") or []
    issue_set = {
        str(iss).strip().upper()
        for iss in (raw_issues if isinstance(raw_issues, (list, tuple, set)) else [raw_issues])
        if iss
    }
    georef_qc_status = "valid"
    lat = rec.get("decimalLatitude")
    lon = rec.get("decimalLongitude")
    lat_val: Optional[float] = None
    lon_val: Optional[float] = None
    if lat is not None and lon is not None:
        try:
            lf = float(lat)
            lnf = float(lon)
            if math.isfinite(lf) and math.isfinite(lnf):
                lat_val, lon_val = lf, lnf
        except (ValueError, TypeError):
            pass

    state_prov = rec.get("stateProvince")
    rec_country = str(rec.get("countryCode") or "").strip().upper()
    hit_cks: List[str] = []
    point_cls = _harv.Point
    np_mod = _harv.np

    if lat_val is not None and lon_val is not None:
        if issue_set & _FATAL_GBIF_GEO_ISSUES:
            lat_val, lon_val = None, None
            hit_cks = []
            georef_qc_status = "gbif_geo_issue_nullified"
        elif not (-90.0 <= lat_val <= 90.0 and -180.0 <= lon_val <= 180.0):
            lat_val, lon_val = None, None
            hit_cks = []
            georef_qc_status = "out_of_bounds_nullified"
        elif abs(lat_val) < 0.01 and abs(lon_val) < 0.01:
            lat_val, lon_val = None, None
            hit_cks = []
            georef_qc_status = "zero_coord_nullified"
        elif (
            lon_val > 0.0
            and (rec_country in {"US", "CA", "MX"} or (not rec_country and harvester_self.geographic_scope in {"eastern_na", "north_america"}))
        ):
            lat_val, lon_val = None, None
            hit_cks = []
            georef_qc_status = "hemisphere_sign_mismatch_nullified"
        elif state_prov:
            st_dwc_norm = _normalize_us_state_name(state_prov)
            checked_by_poly = False
            if st_dwc_norm and st_dwc_norm in _US_STATE_ABBR_TO_FULL.values() and point_cls is not None and us_spatial_graph.get("tree") is not None:
                pt_chk = point_cls(lon_val, lat_val)
                tree_c = us_spatial_graph["tree"]
                geoms_c = us_spatial_graph["geoms"]
                cks_c = us_spatial_graph["ck_list"]
                gid_c = us_spatial_graph["geom_id_map"]
                neigh_c = us_spatial_graph.get("neighbors") or {}
                hit_cks = []
                for h_it in tree_c.query(pt_chk):
                    h_idx = int(h_it) if (np_mod is not None and isinstance(h_it, (int, np_mod.integer))) or isinstance(h_it, int) else gid_c.get(id(h_it), -1)
                    if h_idx >= 0 and geoms_c[h_idx].covers(pt_chk):
                        hit_cks.append(cks_c[h_idx])
                if hit_cks:
                    checked_by_poly = True
                    ok_state = any(
                        ck_h.split("|", 1)[0] == st_dwc_norm
                        or any(nb.split("|", 1)[0] == st_dwc_norm for nb in neigh_c.get(ck_h, ()))
                        for ck_h in hit_cks
                    )
                    if not ok_state:
                        lat_val, lon_val = None, None
                        hit_cks = []
                        georef_qc_status = "state_polygon_mismatch_nullified"
            if not checked_by_poly and lat_val is not None and lon_val is not None:
                bbox_key = st_dwc_norm if (st_dwc_norm and st_dwc_norm in _NA_STATE_PROVINCE_BBOXES) else _normalize_ca_province_bbox_key(state_prov)
                if bbox_key and bbox_key in _NA_STATE_PROVINCE_BBOXES:
                    b_lat_min, b_lat_max, b_lon_min, b_lon_max = _NA_STATE_PROVINCE_BBOXES[bbox_key]
                    margin_deg = 0.5
                    if not (
                        (b_lat_min - margin_deg) <= lat_val <= (b_lat_max + margin_deg)
                        and (b_lon_min - margin_deg) <= lon_val <= (b_lon_max + margin_deg)
                    ):
                        lat_val, lon_val = None, None
                        hit_cks = []
                        georef_qc_status = "state_polygon_mismatch_nullified"
    else:
        georef_qc_status = "ungeoreferenced"

    is_west_rec = _harv.is_excluded_western_region(
        state_prov,
        lat=lat_val,
        lon=lon_val,
        western_longitude_threshold=harvester_self.western_longitude_threshold,
        exclude_western=True,
        is_type_specimen=is_rec_type,
    )
    has_any_state = bool(
        state_prov is not None
        and str(state_prov).strip()
        and str(state_prov).strip().upper() not in {"NONE", "NAN", "NA", "N/A", "NULL", "UNKNOWN", "UNSPECIFIED", "UNRECORDED"}
    )
    if harvester_self.exclude_western:
        if not harvester_self.harvest_full_continental_range and is_west_rec:
            counters["western_excluded"] = counters.get("western_excluded", 0) + 1
            return
        # Un-georeferenced non-types lacking any stateProvince cannot be placed biogeographically
        if lat_val is None and lon_val is None and not has_any_state and not is_rec_type:
            counters["western_excluded"] = counters.get("western_excluded", 0) + 1
            return

    # Preserve missing coordinateUncertaintyInMeters as None (NA) rather than fabricating 5000 m
    uncertainty_val: Optional[float] = None
    if lat_val is not None and lon_val is not None:
        uncertainty_raw = rec.get("coordinateUncertaintyInMeters")
        if uncertainty_raw is not None and str(uncertainty_raw).strip() != "":
            try:
                parsed_unc = float(uncertainty_raw)
                if math.isfinite(parsed_unc) and parsed_unc >= 0.0:
                    if parsed_unc > harvester_self.max_uncertainty_meters and not is_rec_type:
                        return
                    uncertainty_val = parsed_unc
            except (ValueError, TypeError):
                pass

    # 2. Darwin Core Temporal Validation (retaining month/year historical sheets for morphometrics
    #    with day=None -> doy=NA, and masking YYYY-01-01 digitization placeholders)
    raw_year = rec.get("year")
    raw_month = rec.get("month")
    raw_day = rec.get("day")
    raw_ev_str = str(rec.get("eventDate") or "").strip()

    year_val: Optional[int] = None
    month_val: Optional[int] = None
    day_val: Optional[int] = None
    try:
        if raw_year is not None and str(raw_year).strip() != "":
            year_val = int(float(raw_year))
    except (ValueError, TypeError):
        year_val = None
    if year_val is None and raw_ev_str:
        m_yr = re.match(r"^(\d{4})\b", raw_ev_str)
        if m_yr:
            year_val = int(m_yr.group(1))

    if year_val is not None and not (1700 <= year_val <= current_year):
        year_val = None
    if year_val is None and not is_rec_type:
        return

    try:
        if raw_month is not None and str(raw_month).strip() != "":
            m_cand = int(float(raw_month))
            if 1 <= m_cand <= 12:
                month_val = m_cand
    except (ValueError, TypeError):
        month_val = None
    try:
        if raw_day is not None and str(raw_day).strip() != "":
            d_cand = int(float(raw_day))
            if 1 <= d_cand <= 31:
                day_val = d_cand
    except (ValueError, TypeError):
        day_val = None

    if year_val is not None and month_val is not None and day_val is not None:
        try:
            valid_date = datetime.date(year_val, month_val, day_val)
            # Guard against YYYY-01-01 Symbiota/GBIF digitization placeholders unless verbatimEventDate explicitly names January
            verbatim_ev = str(rec.get("verbatimEventDate") or "").strip().lower()
            explicit_jan = bool(re.search(r"\bjan(?:uary)?\b", verbatim_ev))
            is_jan1_placeholder = (
                month_val == 1
                and day_val == 1
                and not explicit_jan
                and (
                    "RECORDED_DATE_UNLIKELY" in issue_set
                    or bool(re.fullmatch(r"\d{4}(?:-01-01(?:T00:00:00Z?)?)?", raw_ev_str))
                    or not raw_ev_str
                )
            )
            if is_jan1_placeholder:
                month_val, day_val = None, None
                event_date = f"{year_val:04d}"
            else:
                event_date = raw_ev_str or valid_date.isoformat()
        except ValueError:
            if not is_rec_type:
                return
            day_val = None
            event_date = raw_ev_str or f"{year_val:04d}-{month_val:02d}"
    elif year_val is not None and month_val is not None:
        day_val = None
        event_date = raw_ev_str or f"{year_val:04d}-{month_val:02d}"
    else:
        month_val, day_val = None, None
        event_date = raw_ev_str or (f"{year_val:04d}" if year_val is not None else "")

    # 3. High-Resolution Media Image Validation
    media_list = rec.get("media", [])
    image_url = _harv.extract_high_res_image_url(media_list)
    if not image_url:
        return

    # 4. Catalog Number & Herbarium Institution Namespace Protection
    raw_catalog = rec.get("catalogNumber")
    inst_code = str(rec.get("institutionCode") or rec.get("collectionCode") or "UNKNOWN_INST").strip() or "UNKNOWN_INST"
    inst_norm_up = inst_code.upper()
    gbif_key = str(rec.get("key", ""))

    if raw_catalog and str(raw_catalog).strip():
        raw_cat_str = str(raw_catalog).strip()
        clean_cat = _harv.sanitize_filename(raw_cat_str)
        # Prefix institutionCode onto bare numeric/non-alphabetic accession numbers (e.g. '12345' -> 'GA_12345')
        if not re.search(r"[A-Za-z]", raw_cat_str) and inst_norm_up not in _harv._INVALID_INSTITUTION_CODES:
            catalog_number = f"{_harv.sanitize_filename(inst_code)}_{clean_cat}"
        else:
            catalog_number = clean_cat
    else:
        catalog_number = f"{_harv.sanitize_filename(inst_code)}_{gbif_key}"

    prev_inst = catalog_to_inst.get(catalog_number)
    if prev_inst is not None and prev_inst != inst_norm_up and inst_norm_up not in _harv._INVALID_INSTITUTION_CODES:
        catalog_number = f"{_harv.sanitize_filename(inst_code)}_{catalog_number}"

    unique_key = (catalog_number, gbif_key)
    if unique_key in seen_catalog_keys:
        return
    seen_catalog_keys.add(unique_key)
    catalog_to_inst[catalog_number] = inst_norm_up

    relative_image_path = f"data/raw_vouchers/{catalog_number}.webp"
    raw_sci_val = str(rec.get("scientificName") or rec.get("species") or taxon).strip()
    infraspecific_val = _extract_infraspecific_taxon(raw_sci_val)

    # Dynamically map historical non-Packera genera (Senecio / Cineraria) AND occurrence records
    # explicitly flagged by GBIF as unaccepted synonyms (e.g. Packera tomentosa -> Packera dubia),
    # while protecting explicitly requested accepted Packera segregates.
    raw_packera_bin = _extract_canonical_binomial(raw_sci_val, {"Packera"})
    rec_tax_status = str(rec.get("taxonomicStatus") or rec.get("status") or "").strip().upper()
    is_rec_synonym = rec_tax_status in _harv._GBIF_SYNONYM_STATUSES
    is_protected_raw = bool(
        raw_packera_bin
        and raw_packera_bin.lower() in harvester_self._protected_packera_taxa
        and not (harvester_self.expand_synonyms and is_rec_synonym)
    )
    acc_from_rec: Optional[str] = None
    if raw_packera_bin is None or (is_rec_synonym and not is_protected_raw):
        for acc_fld in ("acceptedScientificName", "species"):
            cand_acc = _extract_canonical_binomial(rec.get(acc_fld), {"Packera"})
            if cand_acc:
                acc_from_rec = cand_acc
                break
        if acc_from_rec and (raw_packera_bin is None or acc_from_rec.lower() != raw_packera_bin.lower()):
            if acc_from_rec.lower() not in (harvester_self._protected_packera_taxa - {acc_from_rec.lower()}) or not is_protected_raw:
                syn_items_to_reg = [raw_sci_val] + ([raw_packera_bin] if raw_packera_bin else [])
                _harv._register_dynamic_synonyms(
                    acc_from_rec,
                    syn_items_to_reg,
                    [raw_sci_val],
                    allow_packera_species_remap=is_rec_synonym,
                )

    resolved_packera = _harv.standardize_packera_taxon(raw_sci_val) or acc_from_rec
    concept_info = _harv.parse_voucher_taxonomic_concept(
        rec=rec,
        resolved_species=resolved_packera,
        is_type_specimen=is_rec_type,
        synonym_genera=_SYNONYM_GENERA,
    )
    if not resolved_packera and concept_info.get("species_standardized_hint"):
        resolved_packera = (
            _harv.standardize_packera_taxon(concept_info["species_standardized_hint"])
            or concept_info["species_standardized_hint"]
        )
    if concept_info.get("infraspecific_taxon"):
        infraspecific_val = str(concept_info["infraspecific_taxon"])
    elif concept_info.get("is_open_or_hybrid_concept"):
        infraspecific_val = ""
    if is_direct_genus_query or harvester_self.is_genus_wide:
        if not resolved_packera or not re.match(r"^Packera\s+(?:[x×]\s*)?[a-z-]+$", resolved_packera):
            seen_catalog_keys.discard(unique_key)
            return
    is_open_or_hybrid = bool(concept_info.get("is_open_or_hybrid_concept"))
    if resolved_packera and (is_open_or_hybrid or not raw_sci_val.lower().startswith(resolved_packera.lower())):
        sci_name_val = resolved_packera
    else:
        sci_name_val = raw_sci_val
    species_raw_val = raw_sci_val
    typ_bas_val = str(concept_info.get("typified_basionym") or "").strip()
    if (
        is_rec_type
        and typ_bas_val
        and resolved_packera
        and raw_sci_val.lower().startswith(resolved_packera.lower())
        and not typ_bas_val.lower().startswith(resolved_packera.lower())
    ):
        species_raw_val = str(
            rec.get("verbatimScientificName")
            or rec.get("typifiedName")
            or typ_bas_val
        ).strip()

    raw_loc_val = rec.get("locality") or rec.get("verbatimLocality")
    raw_hab_val = rec.get("habitat")

    # 5. Initial Layer 1 + Layer 2 Determiner Tier Stratification
    determiner_tier, type_status, identified_by_clean = _harv.parse_determiner_tier(
        type_status_raw=rec.get("typeStatus"),
        identified_by_raw=rec.get("identifiedBy"),
        recorded_by_raw=rec.get("recordedBy"),
        history_raw=rec.get("verbatimIdentificationHistory"),
        institution_code_raw=inst_code,
        locality_raw=raw_loc_val,
        habitat_raw=raw_hab_val,
        dynamic_specialists=None,
        is_open_or_hybrid=is_open_or_hybrid,
    )

    # 6. North American (CEC / US EPA) Level III Ecoregion Polygon Assignment
    regional_group = _harv.infer_regional_group(
        lat=lat_val,
        lon=lon_val,
        state_province=state_prov,
        habitat=raw_hab_val,
        locality=raw_loc_val,
        county=rec.get("county"),
    )

    verbatim_conflict, verbatim_gold_taxon = _harv._detect_verbatim_identification_conflict(
        rec=rec,
        resolved_canonical=(resolved_packera or sci_name_val),
        determiner_tier=determiner_tier,
    )
    is_cent_q, coord_q_reason = _harv.evaluate_coordinate_centroid_quarantine(
        lat=lat_val,
        lon=lon_val,
        locality=str(raw_loc_val or ""),
        habitat=str(raw_hab_val or ""),
        county=str(rec.get("county") or ""),
        state_prov=str(state_prov or ""),
        uncertainty_m=uncertainty_val,
        issues=issue_set,
        is_type_specimen=bool(is_rec_type or _harv.classify_nomenclatural_type_rank(type_status)[0] >= 1),
        us_spatial_graph=us_spatial_graph,
        hit_county_keys=hit_cks,
    )
    curated_record = {
        "catalogNumber": catalog_number,
        "institutionCode": inst_code,
        "scientificName": sci_name_val,
        "species_raw": species_raw_val,
        "infraspecific_taxon": infraspecific_val,
        "typified_basionym": typ_bas_val,
        "verbatim_taxon_concept": str(concept_info.get("verbatim_taxon_concept") or ""),
        "concept_synonymy_status": str(concept_info.get("concept_synonymy_status") or "Accepted_Sensu_Stricto"),
        "identification_year": concept_info.get("identification_year"),
        "is_open_or_hybrid_concept": is_open_or_hybrid,
        "georef_qc_status": georef_qc_status,
        "is_centroid_or_inst_coord": bool(is_cent_q),
        "coord_quarantine_reason": str(coord_q_reason or ""),
        "recordedBy": str(rec.get("recordedBy") or "").strip(),
        "recordNumber": str(rec.get("recordNumber") or "").strip(),
        "identifiedBy": identified_by_clean,
        "determiner_tier": determiner_tier,
        "type_status": type_status,
        "county": rec.get("county") or "",
        "stateProvince": state_prov or "",
        "decimalLatitude": lat_val,
        "decimalLongitude": lon_val,
        "coordinateUncertainty": uncertainty_val,
        "year": year_val,
        "month": month_val,
        "day": day_val,
        "eventDate": event_date,
        "regional_group": regional_group,
        "image_path": relative_image_path,
        "exsiccatae_taxonomic_conflict": bool(verbatim_conflict),
        "exsiccatae_gold_taxon": verbatim_gold_taxon,
        "_image_url": image_url,
        "_is_western_rec": bool(is_west_rec),
        "_is_type_rec": bool(is_rec_type),
        "_raw_identified_by": str(rec.get("identifiedBy") or "").strip(),
        "_raw_date_identified": str(rec.get("dateIdentified") or "").strip(),
        "_raw_history": str(rec.get("verbatimIdentificationHistory") or "").strip(),
        "_raw_type_status": str(rec.get("typeStatus") or "").strip(),
        "_raw_locality": str(raw_loc_val or "").strip(),
        "_raw_habitat": str(raw_hab_val or "").strip(),
    }
    all_curated_records.append(curated_record)
    if gbif_key and seen_gbif_keys is not None:
        seen_gbif_keys.add(gbif_key)
    counters["retained"] = counters.get("retained", 0) + 1
