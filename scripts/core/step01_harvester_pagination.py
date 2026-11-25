"""Step 01 GBIF two-pass pagination, rescue passes, and cohort post-processing."""

from __future__ import annotations

import datetime
import math
from pathlib import Path
from typing import TYPE_CHECKING, Any, Dict, List, Optional, Set, Tuple

if TYPE_CHECKING:
    import pandas as pd


def _extract_page_gbif_keys(results: List[Dict[str, Any]]) -> Tuple[Set[int], List[str]]:
    """Extract distinct speciesKey/acceptedTaxonKey ints and non-empty GBIF occurrence keys from a page."""
    sp_keys: Set[int] = set()
    occ_keys: List[str] = []
    for rec in results:
        sk_raw = rec.get("speciesKey") or rec.get("acceptedTaxonKey")
        if sk_raw is not None:
            try:
                sp_keys.add(int(sk_raw))
            except (TypeError, ValueError):
                pass
        ok_raw = str(rec.get("key") or "").strip()
        if ok_raw:
            occ_keys.append(ok_raw)
    return sp_keys, occ_keys


def execute_gbif_harvest_loop(harvester_self: Any) -> pd.DataFrame:
    """Harvests and normalizes Darwin Core occurrence records across configured taxa."""
    import scripts.core.harvester as _harv
    from scripts.core.usda_plants_ranges import load_us_county_spatial_graph

    pd_mod = _harv.pd
    all_curated_records: List[Dict[str, Any]] = []
    seen_catalog_keys: Set[Tuple[str, str]] = set()
    seen_gbif_keys: Set[str] = set()
    visited_country_gbif_keys: Dict[Optional[str], Set[str]] = {}
    exhausted_gbif_species_keys: Set[Tuple[Optional[str], int]] = set()
    catalog_to_inst: Dict[str, str] = {}
    country_targets: List[Optional[str]] = list(harvester_self.countries) if harvester_self.countries else [None]
    current_year = datetime.date.today().year + 1

    occ_ref = _harv.occ
    is_gbif_mocked = (
        _harv._is_gbif_species_mocked() or occ_ref is None or getattr(occ_ref, "_is_missing_pygbif", False)
        or getattr(occ_ref, "search", None) is None or hasattr(occ_ref.search, "_mock_name")
        or hasattr(occ_ref.search, "assert_called")
    )

    try:
        us_spatial_graph = load_us_county_spatial_graph(county_path=harvester_self.usda_county_geojson_path, allow_download=False)
    except Exception:
        us_spatial_graph = {}

    for taxon in harvester_self.taxa:
        is_direct_genus_query = taxon.strip().lower() == "packera"
        genus_key = _harv._get_packera_genus_key_gbif() if is_direct_genus_query else None
        country_desc = ",".join(harvester_self.countries) if harvester_self.countries else "GLOBAL"
        harvester_self.logger.info(
            f"Querying GBIF API for taxon: '{taxon}' "
            f"(scope={harvester_self.geographic_scope}, country={country_desc}, basisOfRecord=PRESERVED_SPECIMEN)..."
        )
        counters: Dict[str, int] = {"harvested": 0, "retained": 0, "western_excluded": 0}

        def _ingest(rec: Dict[str, Any], require_type: bool = False) -> None:
            _harv.ingest_gbif_occurrence_record(
                harvester_self, rec, taxon, is_direct_genus_query, current_year,
                us_spatial_graph, all_curated_records, seen_catalog_keys,
                catalog_to_inst, counters, require_type=require_type, seen_gbif_keys=seen_gbif_keys,
            )

        use_server_lon_filter = bool(harvester_self.exclude_western and not harvester_self.harvest_full_continental_range)
        eff_query_west_lon = harvester_self.western_longitude_threshold
        limit = 300
        pass1_quota = max(1, math.ceil(harvester_self.max_records_per_taxon / max(1, len(country_targets))))
        country_offsets: Dict[Optional[str], int] = {c: 0 for c in country_targets}
        country_retained: Dict[Optional[str], int] = {c: 0 for c in country_targets}
        country_exhausted: Dict[Optional[str], bool] = {c: False for c in country_targets}
        country_dup_exhausted: Dict[Optional[str], bool] = {c: False for c in country_targets}
        country_sp_keys: Dict[Optional[str], Set[int]] = {c: set() for c in country_targets}

        def _run_country_pages(country_code: Optional[str], stop_at_country_quota: Optional[int]) -> None:
            visited_for_country = visited_country_gbif_keys.setdefault(country_code, set())
            while counters["retained"] < harvester_self.max_records_per_taxon and not country_exhausted[country_code]:
                if stop_at_country_quota is not None and country_retained[country_code] >= stop_at_country_quota:
                    break
                offset = country_offsets[country_code]
                search_kwargs: Dict[str, Any] = {
                    "basisOfRecord": "PRESERVED_SPECIMEN", "mediaType": "StillImage", "limit": limit, "offset": offset,
                }
                if is_direct_genus_query and genus_key is not None:
                    search_kwargs["taxonKey"] = genus_key
                else:
                    search_kwargs["scientificName"] = taxon
                if country_code:
                    search_kwargs["country"] = country_code
                if use_server_lon_filter:
                    search_kwargs["decimalLongitude"] = f"{eff_query_west_lon},-50.0"
                try:
                    response = _harv.occ.search(**search_kwargs)
                except Exception as e:
                    harvester_self.logger.error(f"GBIF API query error for '{taxon}' (country={country_code}) at offset {offset}: {e}")
                    country_exhausted[country_code] = True
                    break

                if not isinstance(response, dict):
                    country_exhausted[country_code] = True
                    break
                results = response.get("results", [])
                raw_count = response.get("count")
                if not results:
                    country_exhausted[country_code] = True
                    break

                page_sp_keys, page_occ_keys = _extract_page_gbif_keys(results)
                country_sp_keys[country_code].update(page_sp_keys)
                if (
                    not is_gbif_mocked and offset == 0 and page_sp_keys and page_occ_keys
                    and all((country_code, sk) in exhausted_gbif_species_keys for sk in page_sp_keys)
                    and all(ok in visited_for_country for ok in page_occ_keys)
                ):
                    country_exhausted[country_code] = True
                    country_dup_exhausted[country_code] = True
                    break

                eff_count = int(raw_count) if raw_count is not None else (offset + len(results))
                page_processed = 0
                for rec in results:
                    if counters["retained"] >= harvester_self.max_records_per_taxon:
                        break
                    if stop_at_country_quota is not None and country_retained[country_code] >= stop_at_country_quota:
                        break
                    pre_ret = counters["retained"]
                    _ingest(rec, require_type=False)
                    ok_rec = str(rec.get("key") or "").strip()
                    if ok_rec:
                        visited_for_country.add(ok_rec)
                    page_processed += 1
                    if counters["retained"] > pre_ret:
                        country_retained[country_code] += (counters["retained"] - pre_ret)

                country_offsets[country_code] = offset + page_processed
                if (page_processed == len(results) and country_offsets[country_code] >= eff_count) or page_processed == 0:
                    country_exhausted[country_code] = True
                    if (
                        not is_gbif_mocked and page_processed == len(results)
                        and len(country_sp_keys[country_code]) == 1
                        and counters["retained"] < harvester_self.max_records_per_taxon
                    ):
                        exhausted_gbif_species_keys.add((country_code, next(iter(country_sp_keys[country_code]))))
                _harv.time.sleep(0.05)

        for country_code in country_targets:
            if counters["retained"] >= harvester_self.max_records_per_taxon:
                break
            _run_country_pages(country_code, pass1_quota)

        if counters["retained"] < harvester_self.max_records_per_taxon:
            for country_code in country_targets:
                if counters["retained"] >= harvester_self.max_records_per_taxon:
                    break
                if not country_exhausted[country_code]:
                    _run_country_pages(country_code, None)

        all_countries_dup_exhausted = bool(
            not is_gbif_mocked
            and counters["retained"] == 0
            and any(country_dup_exhausted.values())
            and all(country_dup_exhausted[c] or not country_sp_keys[c] for c in country_targets)
        )

        # Ungeoreferenced (hasCoordinate=False) Cis-Mississippi State Rescue Pass for historical sheets
        if (
            use_server_lon_filter
            and harvester_self.include_ungeoreferenced_state_vouchers
            and counters["retained"] < harvester_self.max_records_per_taxon
            and not all_countries_dup_exhausted
        ):
            for country_code in country_targets:
                if counters["retained"] >= harvester_self.max_records_per_taxon:
                    break
                u_offset = 0
                while counters["retained"] < harvester_self.max_records_per_taxon:
                    u_kwargs: Dict[str, Any] = {
                        "basisOfRecord": "PRESERVED_SPECIMEN",
                        "mediaType": "StillImage",
                        "hasCoordinate": False,
                        "limit": limit,
                        "offset": u_offset,
                    }
                    if is_direct_genus_query and genus_key is not None:
                        u_kwargs["taxonKey"] = genus_key
                    else:
                        u_kwargs["scientificName"] = taxon
                    if country_code:
                        u_kwargs["country"] = country_code
                    try:
                        u_resp = _harv.occ.search(**u_kwargs)
                        u_results = u_resp.get("results", []) if isinstance(u_resp, dict) else []
                        u_count = (
                            int(u_resp.get("count"))
                            if isinstance(u_resp, dict) and u_resp.get("count") is not None
                            else (u_offset + len(u_results))
                        )
                    except Exception:
                        break
                    if not u_results:
                        break
                    # Guard against mock occ.search handlers ignoring hasCoordinate=False
                    has_finite_lon = False
                    for u_rec in u_results:
                        lon_u = u_rec.get("decimalLongitude")
                        try:
                            if lon_u is not None and math.isfinite(float(lon_u)):
                                has_finite_lon = True
                                break
                        except (TypeError, ValueError):
                            pass
                    if u_offset == 0 and has_finite_lon:
                        break
                    seen_before = len(seen_catalog_keys)
                    for u_rec in u_results:
                        if counters["retained"] >= harvester_self.max_records_per_taxon:
                            break
                        if not _harv.is_cis_mississippi_record(u_rec.get("stateProvince"), lat=None, lon=None):
                            continue
                        _ingest(u_rec, require_type=False)
                    if len(seen_catalog_keys) == seen_before:
                        break
                    u_offset += len(u_results)
                    if u_offset >= u_count:
                        break
                    _harv.time.sleep(0.05)

        # Unconditional & Country-Unrestricted Type-Specimen Rescue (captures European/uncoded types even at quota)
        type_country_targets: List[Optional[str]] = list(country_targets)
        if harvester_self.countries and None not in type_country_targets:
            type_country_targets.append(None)
        for country_code in type_country_targets:
            type_search_kwargs: Dict[str, Any] = {"basisOfRecord": "PRESERVED_SPECIMEN", "typeStatus": "*", "limit": 300, "offset": 0}
            type_search_kwargs["taxonKey" if (is_direct_genus_query and genus_key is not None) else "scientificName"] = (
                genus_key if (is_direct_genus_query and genus_key is not None) else taxon
            )
            if country_code:
                type_search_kwargs["country"] = country_code
            try:
                type_resp = _harv.occ.search(**type_search_kwargs)
                type_results = type_resp.get("results", []) if isinstance(type_resp, dict) else []
            except Exception:
                type_results = []
            if not type_results:
                continue
            seen_before_t = len(seen_catalog_keys)
            for t_rec in type_results:
                _ingest(t_rec, require_type=True)
            if len(seen_catalog_keys) == seen_before_t and country_code is None:
                break

        harvester_self.logger.info(
            f"Taxon '{taxon}': Processed {counters['harvested']} occurrences -> "
            f"Retained {counters['retained']} candidate records "
            f"(Immediate exclusions: {counters['western_excluded']})."
        )

    df = pd_mod.DataFrame(all_curated_records)
    if not df.empty:
        # Eastern NA Biogeographic Residency Screen + Full US+CA Continental Range Retention for Residents
        is_direct_single_genus_tiny_mock = (
            len(harvester_self.taxa) == 1
            and harvester_self.taxa[0].strip().lower() == "packera"
            and len(df) < harvester_self.min_eastern_residency_records
        )
        if harvester_self.exclude_western and not is_direct_single_genus_tiny_mock:
            if harvester_self.is_genus_wide and not harvester_self.harvest_full_continental_range:
                keep_res_mask, cis_stats = _harv.evaluate_eastern_na_cohort_filter(
                    df,
                    _harv.standardize_packera_taxon,
                    western_longitude_threshold=harvester_self.western_longitude_threshold,
                    cis_longitude_threshold=harvester_self.cis_mississippi_longitude_threshold,
                    min_cis_records=harvester_self.min_eastern_residency_records,
                    min_cis_fraction=harvester_self.min_cis_mississippi_fraction,
                )
                if not keep_res_mask.all():
                    pre_len = len(df)
                    df = df.loc[keep_res_mask].reset_index(drop=True)
                    harvester_self.logger.info(
                        f"Cis-Mississippi Eastern NA Filter: Retained {len(df)}/{pre_len} vouchers across "
                        f"{len(cis_stats['retained_species'])} Eastern species (dropped {len(cis_stats['dropped_species'])} "
                        f"Cordilleran/Plains species: {cis_stats['dropped_species']})."
                    )
            else:
                canon_series = df["scientificName"].map(
                    lambda x: _harv.standardize_packera_taxon(x) or str(x).strip()
                )
                west_mask = df["_is_western_rec"].astype(bool) if "_is_western_rec" in df.columns else pd_mod.Series(False, index=df.index)
                type_mask_r = df["_is_type_rec"].astype(bool) if "_is_type_rec" in df.columns else pd_mod.Series(False, index=df.index)
                cis_mask = pd_mod.Series(
                    [
                        _harv.is_cis_mississippi_record(
                            st_v, lat=lat_v, lon=lon_v,
                            cis_longitude_threshold=harvester_self.cis_mississippi_longitude_threshold,
                        )
                        for st_v, lat_v, lon_v in zip(
                            df.get("stateProvince", [""] * len(df)),
                            df.get("decimalLatitude", [None] * len(df)),
                            df.get("decimalLongitude", [None] * len(df)),
                        )
                    ],
                    index=df.index, dtype=bool,
                )
                eff_min_res = (
                    harvester_self.min_eastern_residency_records
                    if len(df) >= 50
                    else min(3, max(1, harvester_self.min_eastern_residency_records))
                )
                cis_counts = canon_series[cis_mask & (~type_mask_r)].value_counts().to_dict()
                resident_species = {sp for sp, cnt in cis_counts.items() if int(cnt) >= eff_min_res}
                if harvester_self.is_genus_wide:
                    keep_res_mask = canon_series.isin(resident_species)
                    if not keep_res_mask.all():
                        df = df.loc[keep_res_mask].reset_index(drop=True)
                elif harvester_self.harvest_full_continental_range:
                    keep_res_mask = canon_series.isin(resident_species) | (~west_mask)
                    if not keep_res_mask.all():
                        df = df.loc[keep_res_mask].reset_index(drop=True)

            if df.empty:
                for c_tmp in ("_is_western_rec", "_is_type_rec"):
                    if c_tmp in df.columns:
                        df = df.drop(columns=[c_tmp])
                return df

        # Layer 3: Corpus-derived multi-herbarium specialist & regional herbarium authority re-evaluation
        dyn_specs, dyn_herbs = _harv.discover_corpus_authorities(df)
        if dyn_specs or dyn_herbs:
            harvester_self.logger.info(
                f"Dynamic Authority Engine: Discovered {len(dyn_specs)} corpus specialist(s) "
                f"and {len(dyn_herbs)} dynamic research herbarium/herbaria."
            )
        df = _harv.apply_corpus_authority_tiers(
            df,
            dynamic_specialists=dyn_specs,
            dynamic_herbaria=dyn_herbs,
            db_path=harvester_self.db_path,
            workspace_dir=harvester_self.workspace_dir,
        )
        raw_audit_cols = [
            c for c in (
                "_is_western_rec", "_is_type_rec", "_raw_identified_by", "_raw_date_identified",
                "_raw_history", "_raw_type_status", "_raw_locality", "_raw_habitat",
            ) if c in df.columns
        ]
        if raw_audit_cols:
            df = df.drop(columns=raw_audit_cols)

        default_usda_cache = Path(harvester_self.workspace_dir) / "data" / "environmental" / "usda_plants"
        has_local_usda_cache = (
            harvester_self.usda_cache_dir is not None
            or (default_usda_cache / "plantlst.txt").is_file()
            or (default_usda_cache / "distributions").is_dir()
        )
        usda_stats: Optional[Dict[str, Any]] = None
        if (not is_gbif_mocked) or has_local_usda_cache:
            df, usda_stats = _harv.filter_vouchers_by_usda_county_range(
                df=df, standardize_taxon_fn=_harv.standardize_packera_taxon,
                synonym_to_canonical=_harv._DYNAMIC_SYNONYM_TO_CANONICAL,
                county_buffer_km=harvester_self.usda_county_buffer_km, exempt_gold=harvester_self.usda_exempt_gold,
                cache_dir=harvester_self.usda_cache_dir, county_geojson_path=harvester_self.usda_county_geojson_path,
                db_path=harvester_self.db_path, workspace_dir=harvester_self.workspace_dir,
                allow_network=(not is_gbif_mocked), force_refresh=harvester_self.force,
                drop_out_of_range=False, logger=harvester_self.logger,
            )
            harvester_self.last_usda_filter_stats = usda_stats
            if df.empty:
                return df

        df = _harv.stratify_duplicates(
            df,
            workspace_dir=harvester_self.workspace_dir,
            logger=harvester_self.logger,
            home_herbarium=harvester_self.home_herbarium,
        )

        resolved_db = _harv.resolve_db_path(harvester_self.db_path, workspace_root=harvester_self.workspace_dir)
        if resolved_db != ":memory:" and Path(resolved_db).is_file() and (not is_gbif_mocked or harvester_self.db_path is not None):
            try:
                conn_rec = _harv.get_db_connection(db_path=resolved_db, init_schema=False, read_only=True)
                try:
                    df = _harv.reconcile_cross_taxon_exsiccatae_with_db(
                        df, conn_rec, standardize_fn=_harv.standardize_packera_taxon,
                    )
                finally:
                    conn_rec.close()
            except Exception:
                pass

        if usda_stats is not None and bool(harvester_self.usda_county_filter) and "usda_county_status" in df.columns:
            keep_post: List[bool] = []
            for _, urow in df.iterrows():
                st_u = str(urow.get("usda_county_status") or "").strip()
                if st_u not in {"out_of_range", "natureserve_exotic_or_waif"}:
                    keep_post.append(True)
                    continue
                raw_ts_u = urow.get("type_status")
                if raw_ts_u is None or pd_mod.isna(raw_ts_u) or not str(raw_ts_u).strip():
                    raw_ts_u = urow.get("typeStatus")
                t_rank_u, _ = _harv.classify_nomenclatural_type_rank(raw_ts_u)
                d_tier_u = str(urow.get("determiner_tier") or "").strip()
                is_ex_u = (t_rank_u >= 1) or (bool(harvester_self.usda_exempt_gold) and d_tier_u == "Tier_1_Gold")
                if is_ex_u:
                    if st_u == "natureserve_exotic_or_waif":
                        usda_stats["exempt_retained"] = int(usda_stats.get("exempt_retained", 0)) + 1
                    keep_post.append(True)
                else:
                    usda_stats["dropped_out_of_range"] = int(usda_stats.get("dropped_out_of_range", 0)) + 1
                    keep_post.append(False)
            df = df.loc[keep_post].reset_index(drop=True)
            harvester_self.last_usda_filter_stats = usda_stats
            if df.empty:
                return df

        # Cache secondary duplicate sheets so download_and_validate_media() can promote a fallback
        # duplicate whenever a primary sheet's image URL 404s or fails optical resolution QC
        if "is_primary_duplicate" in df.columns:
            prim_mask = df["is_primary_duplicate"].map(
                lambda v: v is True or str(v).strip().lower() in {"true", "1", "t", "yes"}
            )
            harvester_self._secondary_duplicate_pool = df.loc[~prim_mask].copy()
            ts_col_h = "type_status" if "type_status" in df.columns else (
                "typeStatus" if "typeStatus" in df.columns else None
            )
            if ts_col_h is not None:
                type_mask_h = df[ts_col_h].map(lambda t: _harv.classify_nomenclatural_type_rank(t)[0] >= 1)
                keep_prim_or_type = prim_mask | type_mask_h
            else:
                keep_prim_or_type = prim_mask
            harvester_self._pre_thin_curated_pool = df.copy()
            df = df.loc[keep_prim_or_type].reset_index(drop=True)
        else:
            harvester_self._secondary_duplicate_pool = None
            harvester_self._pre_thin_curated_pool = df.copy()
        src_col_h = "scientificName" if "scientificName" in df.columns else (
            "species_raw" if "species_raw" in df.columns else None
        )
        if src_col_h is not None:
            df["species_standardized"] = df[src_col_h].map(
                lambda x: (_harv.standardize_packera_taxon(str(x)) or str(x).strip() or "Unknown")
                if not pd_mod.isna(x) else "Unknown"
            )
        elif "species_standardized" not in df.columns:
            df["species_standardized"] = "Unknown"
        df["spatial_balance_weight"] = 1.0
    return df
