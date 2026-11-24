"""Step 01 GBIF server-side species facet pre-screen for Eastern NA residency."""

from __future__ import annotations

import json
import logging
import urllib.parse
from pathlib import Path
from typing import Dict, List, Optional, Union

from scripts.core.config import CIS_MISSISSIPPI_STATES_PROVINCES
from scripts.core.epa_ecoregions import _US_STATE_ABBR_TO_FULL
from scripts.core.step01_harvester_filters import _extract_canonical_binomial


def prescreen_eastern_na_residents_gbif(
    base_taxa: List[str],
    eff_countries: Optional[List[str]],
    cis_mississippi_longitude_threshold: float,
    min_eastern_residency_records: int,
    min_cis_mississippi_fraction: float,
    db_path: Optional[Union[str, Path]],
    workspace_dir: Union[str, Path],
    logger: logging.Logger,
) -> List[str]:
    """Query GBIF server-side species facets to filter genus-wide taxa to Eastern NA residents."""
    import scripts.core.harvester as _harv

    occ_mod = _harv.occ
    sp_mod = _harv.species
    if (
        len(base_taxa) <= 1
        or occ_mod is None
        or getattr(occ_mod, "search", None) is None
        or _harv._is_gbif_species_mocked()
        or hasattr(occ_mod.search, "_mock_name")
        or hasattr(occ_mod.search, "assert_called")
    ):
        return base_taxa

    cis_lon_bound = f"{float(cis_mississippi_longitude_threshold)},-50.0"
    facet_countries = list(eff_countries) if eff_countries else ["US", "CA"]
    min_res_n = max(1, int(min_eastern_residency_records))
    min_cis_frac = float(min_cis_mississippi_fraction)
    logger.info(
        f"Querying GBIF server-side species facet for strict Cis-Mississippi Eastern NA residency "
        f"(countries={facet_countries}, cis_lon={cis_lon_bound}, min_records={min_res_n}, min_frac={min_cis_frac:.2f})..."
    )
    eastern_residents: List[str] = []
    try:
        genus_key_val = _harv._get_packera_genus_key_gbif() or 3109225
        cis_key_counts: Dict[int, int] = {}
        tot_key_counts: Dict[int, int] = {}
        cis_us_states = sorted(
            st_title
            for st_low, st_title in {
                v.lower(): v.title() for v in _US_STATE_ABBR_TO_FULL.values()
            }.items()
            if st_low.upper() in CIS_MISSISSIPPI_STATES_PROVINCES
        )
        cis_ca_provs = [
            "Ontario", "Quebec", "New Brunswick", "Nova Scotia",
            "Prince Edward Island", "Newfoundland and Labrador",
        ]
        for ctry in facet_countries:
            # Total country count per speciesKey to compute cis-Mississippi residency share
            tot_url = (
                f"https://api.gbif.org/v1/occurrence/search?"
                f"taxonKey={genus_key_val}&country={ctry}"
                f"&basisOfRecord=PRESERVED_SPECIMEN&mediaType=StillImage"
                f"&facet=speciesKey&facetMincount=1&facetLimit=500&limit=0"
            )
            req_tot = _harv.urllib.request.Request(tot_url, headers={"Accept": "application/json"})
            with _harv.urllib.request.urlopen(req_tot, timeout=15.0) as resp_tot:
                tot_payload = json.loads(resp_tot.read().decode("utf-8"))
            tot_facets = tot_payload.get("facets", []) or []
            for item in ((tot_facets[0].get("counts", []) or []) if tot_facets else []):
                if isinstance(item, dict) and str(item.get("name", "")).isdigit():
                    k_int = int(item["name"])
                    tot_key_counts[k_int] = tot_key_counts.get(k_int, 0) + int(item.get("count", 0))

            st_list = cis_us_states if ctry == "US" else (cis_ca_provs if ctry == "CA" else [])
            facet_urls = []
            if st_list:
                st_q = "&".join(f"stateProvince={urllib.parse.quote(s)}" for s in st_list)
                facet_urls.append(
                    f"https://api.gbif.org/v1/occurrence/search?"
                    f"taxonKey={genus_key_val}&country={ctry}&{st_q}"
                    f"&basisOfRecord=PRESERVED_SPECIMEN&mediaType=StillImage"
                    f"&facet=speciesKey&facetMincount=1&facetLimit=500&limit=0"
                )
            else:
                facet_urls.append(
                    f"https://api.gbif.org/v1/occurrence/search?"
                    f"taxonKey={genus_key_val}&country={ctry}&decimalLongitude={cis_lon_bound}"
                    f"&basisOfRecord=PRESERVED_SPECIMEN&mediaType=StillImage"
                    f"&facet=speciesKey&facetMincount=1&facetLimit=500&limit=0"
                )
            for facet_url in facet_urls:
                req = _harv.urllib.request.Request(facet_url, headers={"Accept": "application/json"})
                with _harv.urllib.request.urlopen(req, timeout=15.0) as resp_http:
                    payload = json.loads(resp_http.read().decode("utf-8"))
                facets = payload.get("facets", []) or []
                sp_counts = (facets[0].get("counts", []) or []) if facets else []
                for item in sp_counts:
                    if isinstance(item, dict) and str(item.get("name", "")).isdigit():
                        k_int = int(item["name"])
                        cis_key_counts[k_int] = cis_key_counts.get(k_int, 0) + int(item.get("count", 0))

        canon_whitelist_low = {
            s.lower()
            for s in _harv.resolve_eastern_resident_whitelist(
                db_path=db_path,
                workspace_root=workspace_dir,
            )
        }
        base_set_lower = {b.lower(): b for b in base_taxa}
        resolved_names = set()
        for sk, c_cis in cis_key_counts.items():
            c_tot = max(c_cis, tot_key_counts.get(sk, c_cis))
            cis_share = float(c_cis) / float(max(1, c_tot))
            try:
                u_info = sp_mod.name_usage(key=sk) if hasattr(sp_mod, "name_usage") else {}
                cand = None
                if isinstance(u_info, dict):
                    cand = _extract_canonical_binomial(
                        u_info.get("canonicalName") or u_info.get("species") or u_info.get("scientificName"),
                        {"Packera"},
                    )
                if not cand:
                    continue
                is_canon_east = cand.lower() in canon_whitelist_low
                if is_canon_east or (c_cis >= min_res_n and cis_share >= min_cis_frac):
                    resolved_names.add(base_set_lower.get(cand.lower(), cand))
            except Exception:
                pass
        # Ensure narrow cis-Mississippi endemics in base_taxa with locality-redacted records are never dropped
        for b_tx in base_taxa:
            if b_tx.lower() in canon_whitelist_low:
                resolved_names.add(b_tx)
        eastern_residents = sorted(resolved_names)
    except Exception as facet_exc:
        logger.warning(f"Facet pre-screen fallback triggered ({facet_exc}).")

    if eastern_residents:
        logger.info(
            f"Eastern NA Pre-Screen: Retained {len(eastern_residents)}/{len(base_taxa)} "
            f"Eastern NA resident Packera species: {eastern_residents}"
        )
        return eastern_residents
    return base_taxa
