#!/usr/bin/env python3
"""Unit tests for Step 01 systematics authority blocking, backbone lumping guards, and coordinate hygiene."""

import tempfile
import unittest
from pathlib import Path
from unittest.mock import MagicMock, patch

import pandas as pd
from PIL import Image


class TestStep1SystematicsAndDownstreamAuthorityAudit(unittest.TestCase):
    """Comprehensive regression tests for the Step 1 botanical systematics & downstream authority-tier review."""

    def test_authority_engine_checklist_blocking_surname_qualification_and_self_id_guard(self):
        """Verify flora/checklist blocking, initial-qualified common surnames, regional floristician Silver routing, and Layer 3 self-ID exclusion."""
        from scripts.core.taxonomy_authority import (
            parse_determiner_tier,
            discover_corpus_authorities,
        )

        # 1. Flora / checklist strings must NOT promote to Tier_1_Gold or Tier_2_Silver
        for flora_str in ("Weakley 2015", "Flora of Virginia (Weakley)", "FNA", "Kartesz", "BONAP", "sensu T.M. Barkley"):
            tier, _, _ = parse_determiner_tier(
                type_status_raw=None,
                identified_by_raw=flora_str,
                recorded_by_raw="Field Collector",
                history_raw=None,
                institution_code_raw="NCU",
                locality_raw="Rich calcareous seepage slope over dolostone along rocky bluff, elev. 420 m",
                habitat_raw="Calcareous seepage",
            )
            self.assertEqual(tier, "Tier_3_Bronze", f"Flora citation '{flora_str}' leaked above Bronze: {tier}")

        # 2. Unrelated collectors/annotators sharing common surnames (Fuller, Turner, Freeman, Weber) do NOT get Tier_1_Gold
        for non_spec in ("S. Fuller", "M. Turner", "J. Freeman", "K. Weber"):
            tier, _, _ = parse_determiner_tier(
                type_status_raw=None,
                identified_by_raw=non_spec,
                recorded_by_raw=non_spec,
                history_raw=None,
                institution_code_raw="UNKNOWN_HERB",
                locality_raw="Roadside",
                habitat_raw=None,
            )
            self.assertEqual(tier, "Tier_3_Bronze", f"Non-specialist '{non_spec}' falsely promoted to {tier}")

        # While initial-qualified core Packera specialists DO get Tier_1_Gold (and B.L. Turner routes Silver)
        for real_spec in ("T.M. Barkley", "Debra K. Trock", "R.R. Kowal", "J.F. Bain", "L.D. Estes", "C.C. Freeman"):
            tier, _, _ = parse_determiner_tier(
                type_status_raw=None,
                identified_by_raw=real_spec,
                recorded_by_raw="Field Collector",
                history_raw=None,
                institution_code_raw="UNKNOWN_HERB",
                locality_raw="Roadside",
                habitat_raw=None,
            )
            self.assertEqual(tier, "Tier_1_Gold", f"Core specialist '{real_spec}' failed Tier_1_Gold: {tier}")
        self.assertEqual(
            parse_determiner_tier(None, "B.L. Turner", "Field Collector", None, "UNKNOWN_HERB", "Roadside", None)[0],
            "Tier_2_Silver",
        )
        self.assertEqual(
            parse_determiner_tier(None, "B. Fuller", "Field Collector", None, "NCU", "Granite outcrop near Wake Forest", None)[0],
            "Tier_2_Silver",
        )
        self.assertEqual(
            parse_determiner_tier(None, "B. Fuller", "Field Collector", None, "UNKNOWN_HERB", "Roadside", None)[0],
            "Tier_3_Bronze",
        )

        # 3. Verbatim identification history does NOT override a later explicit generalist identifiedBy
        tier_reannot, _, _ = parse_determiner_tier(
            type_status_raw=None,
            identified_by_raw="Undergraduate Botany Class",
            recorded_by_raw="Field Collector",
            history_raw="Det. T.M. Barkley, 1978",
            institution_code_raw="NCU",
            locality_raw="Calcareous cliff ledge along river",
            habitat_raw="Limestone bluff",
        )
        self.assertEqual(tier_reannot, "Tier_3_Bronze")

        # 4. Regional floristician field collection with blank identifiedBy routes to Tier_2_Silver (not blind Tier_1_Gold)
        tier_reg_collector, _, _ = parse_determiner_tier(
            type_status_raw=None,
            identified_by_raw="",
            recorded_by_raw="Julian Campbell",
            history_raw=None,
            institution_code_raw="EKY",
            locality_raw="Limestone glade",
            habitat_raw=None,
        )
        self.assertEqual(tier_reg_collector, "Tier_2_Silver")

        # 5. Layer 3 self-identification across distributed duplicates (identifiedBy == recordedBy) must NOT promote collector
        dup_rows = [
            {
                "catalogNumber": f"DUP_{inst}_{i}",
                "institutionCode": inst,
                "scientificName": "Packera aurea" if i % 2 == 0 else "Packera anonyma",
                "identifiedBy": "H.E. Generalcollector",
                "recordedBy": "H.E. Generalcollector",
                "recordNumber": str(100 + i),
                "type_status": "None",
                "eventDate": "2018-05-10",
            }
            for i, inst in enumerate(["NCU", "NY", "MO", "US", "GH", "DUKE"])
        ]
        dyn_specs, _ = discover_corpus_authorities(pd.DataFrame(dup_rows))
        self.assertNotIn("generalcollector", dyn_specs)

    def test_no_packera_to_packera_backbone_lumping_and_infraspecific_persistence(self):
        """Verify GBIF Backbone only auto-maps historical generic transfers (Senecio/Cineraria -> Packera) and saves infraspecific_taxon."""
        import sqlite3
        from scripts.core import harvester as ch
        from scripts.core.harvester import VoucherHarvester, export_curated_table

        ch._DYNAMIC_SYNONYM_TO_CANONICAL.pop("packera pseudaurea", None)
        ch._GBIF_SYNONYM_CACHE.pop("packera pseudaurea", None)
        ch._DYNAMIC_SYNONYM_TO_CANONICAL.pop("senecio aureus", None)
        ch._GBIF_SYNONYM_CACHE.pop("senecio aureus", None)

        mock_sp = MagicMock()
        def _fake_bb(**kwargs):
            nm = str(kwargs.get("name", "")).lower()
            if "pseudaurea" in nm:
                return {"canonicalName": "Packera pseudaurea", "species": "Packera aurea", "accepted": "Packera aurea"}
            if "senecio aureus" in nm:
                return {"canonicalName": "Senecio aureus", "species": "Packera aurea", "accepted": "Packera aurea"}
            return {}
        mock_sp.name_backbone.side_effect = _fake_bb

        with patch("scripts.core.harvester.species", mock_sp):
            # Heterotypic Packera -> Packera backbone lump must NOT occur automatically
            self.assertEqual(ch.standardize_packera_taxon("Packera pseudaurea var. semicordata"), "Packera pseudaurea")
            # Historical Senecio -> Packera generic transfer DOES resolve cleanly
            self.assertEqual(ch.standardize_packera_taxon("Senecio aureus L."), "Packera aurea")

        def _fake_occ_var(**kwargs):
            if "typeStatus" in kwargs or kwargs.get("country") == "CA":
                return {"results": [], "count": 0}
            return {
                "results": [
                    {
                        "key": 55101,
                        "scientificName": "Packera pseudaurea (Rydb.) W.A.Weber & Á.Löve var. semicordata (Mack. & Bush) Trock & T.M.Barkley",
                        "infraspecificEpithet": "semicordata",
                        "taxonRank": "VARIETY",
                        "stateProvince": "Wisconsin",
                        "decimalLatitude": 43.1,
                        "decimalLongitude": -89.4,
                        "year": 1995,
                        "month": 5,
                        "recordedBy": "T.M. Barkley",
                        "identifiedBy": "Debra K. Trock",
                        "recordNumber": "2022-014",
                        "media": [{"identifier": "https://example.com/semicordata.jpg", "type": "StillImage", "format": "image/jpeg"}],
                        "catalogNumber": "12345",
                        "institutionCode": "WIS",
                    },
                    {
                        "key": 55102,
                        "scientificName": "Packera pseudaurea var. semicordata",
                        "stateProvince": "Minnesota",
                        "decimalLatitude": 44.9,
                        "decimalLongitude": -93.2,
                        "eventDate": "1988-01-01",
                        "year": 1988,
                        "month": 1,
                        "day": 1,
                        "recordedBy": "T.M. Barkley",
                        "identifiedBy": "Debra K. Trock",
                        "recordNumber": "2022-015",
                        "media": [{"identifier": "https://example.com/semicordata2.jpg", "type": "StillImage", "format": "image/jpeg"}],
                        "catalogNumber": "12345",
                        "institutionCode": "MIN",
                    },
                ],
                "count": 2,
            }

        with patch("scripts.core.harvester.occ.search", side_effect=_fake_occ_var):
            harv = VoucherHarvester(taxa=["Packera pseudaurea"], expand_synonyms=False)
            df_v = harv.harvest()

        self.assertEqual(len(df_v), 2)
        # Bare numeric catalogNumbers 12345 from WIS and MIN are namespaced instead of colliding!
        self.assertEqual(set(df_v["catalogNumber"]), {"WIS_12345", "MIN_12345"})
        self.assertEqual(set(df_v["infraspecific_taxon"]), {"var. semicordata"})
        # Missing coordinateUncertaintyInMeters stays NaN/None instead of fabricated 5000.0 m
        self.assertTrue(df_v["coordinateUncertainty"].isna().all())
        # Partial year+month (1995-05) retained, and 1988-01-01 placeholder masked to year-only (month/day NaN)
        min_row = df_v[df_v["catalogNumber"] == "MIN_12345"].iloc[0]
        self.assertTrue(pd.isna(min_row["month"]) and pd.isna(min_row["day"]))

        with tempfile.TemporaryDirectory() as tmp_d:
            tmp_p = Path(tmp_d)
            db_p = tmp_p / "var_test.sqlite"
            csv_p = tmp_p / "curated_vouchers.csv"
            export_curated_table(df_v, csv_p, db_path=db_p)
            conn = sqlite3.connect(str(db_p))
            try:
                rows_infra = dict(conn.execute("SELECT catalogNumber, infraspecific_taxon FROM vouchers").fetchall())
            finally:
                conn.close()
            self.assertEqual(rows_infra["WIS_12345"], "var. semicordata")
            self.assertEqual(rows_infra["MIN_12345"], "var. semicordata")

    def test_coordinate_cleaner_cultivated_filter_anonyma_usda_chesapeake_ecoregion_and_media_rescue(self):
        """Verify CoordinateCleaner rejection, cultivated filtering, Packera anonyma USDA symbol, Chesapeake VA ecoregion, and secondary duplicate media rescue."""
        from scripts.core.harvester import VoucherHarvester
        from scripts.core.usda_plants_ranges import resolve_usda_symbols_for_taxon
        from scripts.core.epa_ecoregions import _normalize_us_county_name, infer_regional_group

        # 1. Packera anonyma dynamic USDA symbol resolution (resolves via GBIF synonym Senecio smallii + algorithmic ICN generic transfer)
        with tempfile.TemporaryDirectory() as tmp_usda:
            p_lst = Path(tmp_usda) / "plantlst.txt"
            p_lst.write_text(
                '"Symbol","Synonym Symbol","Scientific Name with Author","Common Name","Family"\n'
                '"PAAN6","","Senecio smallii Britton","Small\'s ragwort","Asteraceae"\n'
                '"PACR11","","Senecio crawfordii (Britton) G.W. Douglas & Ruyle-Douglas","Crawford\'s ragwort","Asteraceae"\n',
                encoding="utf-8",
            )
            syms_craw = resolve_usda_symbols_for_taxon("Packera crawfordii", cache_dir=Path(tmp_usda), allow_network=False)
            self.assertEqual(syms_craw, ["PACR11"])
            mock_anon_sp = MagicMock()
            mock_anon_sp.name_backbone.return_value = {
                "matchType": "EXACT",
                "status": "ACCEPTED",
                "usageKey": 3109240,
                "canonicalName": "Packera anonyma",
                "species": "Packera anonyma",
            }
            mock_anon_sp.name_usage.return_value = {
                "results": [
                    {"canonicalName": "Senecio smallii", "scientificName": "Senecio smallii Britton", "taxonomicStatus": "SYNONYM"},
                ]
            }
            with patch("scripts.core.harvester.species", mock_anon_sp):
                syms_anon = resolve_usda_symbols_for_taxon("Packera anonyma", cache_dir=Path(tmp_usda), allow_network=False)
            self.assertEqual(syms_anon, ["PAAN6"])

        # 2. Chesapeake VA county name consolidation fix (not overwritten to 'norfolk'; South Norfolk -> 'chesapeake')
        self.assertEqual(_normalize_us_county_name("Chesapeake", "virginia"), "chesapeake")
        self.assertEqual(_normalize_us_county_name("South Norfolk", "virginia"), "chesapeake")
        # Eastern Canadian CEC Level III ecoregion assignment (even when CONUS-only us_eco_l3.geojson is cached on disk!)
        import json
        with tempfile.TemporaryDirectory() as tmp_eco_d:
            conus_only_geojson = Path(tmp_eco_d) / "us_eco_l3.geojson"
            conus_only_geojson.write_text(json.dumps({
                "type": "FeatureCollection",
                "features": [{
                    "type": "Feature",
                    "properties": {"US_L3CODE": "45", "US_L3NAME": "Piedmont", "ecoregion_label": "Piedmont"},
                    "geometry": {
                        "type": "Polygon",
                        "coordinates": [[[-80.0, 35.0], [-78.5, 35.0], [-78.5, 36.5], [-80.0, 36.5], [-80.0, 35.0]]],
                    },
                }],
            }))
            with patch("scripts.core.epa_ecoregions.DEFAULT_EPA_L3_GEOJSON", conus_only_geojson),                  patch("scripts.core.epa_ecoregions._ensure_epa_level3_cache_files", return_value=False):
                eco_ottawa = infer_regional_group(lat=45.42, lon=-75.69, state_province="Ontario")
                eco_erie = infer_regional_group(lat=42.98, lon=-81.25, state_province="Ontario")
                eco_nc = infer_regional_group(lat=35.91, lon=-79.05, state_province="North Carolina")
            self.assertEqual(eco_ottawa, "Eastern_Great_Lakes_Lowlands")
            self.assertEqual(eco_erie, "Lake_Erie_Lowland")
            self.assertEqual(eco_nc, "Piedmont")

        # 3. CoordinateCleaner + Cultivated occurrence filtering during harvest()
        def _fake_occ_hygiene(**kwargs):
            if "typeStatus" in kwargs or kwargs.get("country") == "CA":
                return {"results": [], "count": 0}
            return {
                "results": [
                    # Valid NC sheet
                    {
                        "key": 66001,
                        "scientificName": "Packera anonyma",
                        "stateProvince": "North Carolina",
                        "decimalLatitude": 35.9,
                        "decimalLongitude": -79.05,
                        "year": 2015,
                        "recordedBy": "R. Collector",
                        "recordNumber": "101",
                        "media": [{"identifier": "https://example.com/valid.jpg", "type": "StillImage", "format": "image/jpeg"}],
                        "catalogNumber": "NCU_VALID_1",
                        "institutionCode": "NCU",
                    },
                    # (0, 0) Gulf of Guinea error -> rejected
                    {
                        "key": 66002,
                        "scientificName": "Packera anonyma",
                        "stateProvince": "North Carolina",
                        "decimalLatitude": 0.0,
                        "decimalLongitude": 0.0,
                        "year": 2015,
                        "media": [{"identifier": "https://example.com/zero.jpg", "type": "StillImage", "format": "image/jpeg"}],
                        "catalogNumber": "BAD_ZERO_2",
                        "institutionCode": "NCU",
                    },
                    # Sign-inverted positive longitude in NA -> rejected
                    {
                        "key": 66003,
                        "scientificName": "Packera anonyma",
                        "stateProvince": "North Carolina",
                        "decimalLatitude": 35.9,
                        "decimalLongitude": 79.05,
                        "year": 2015,
                        "media": [{"identifier": "https://example.com/poslon.jpg", "type": "StillImage", "format": "image/jpeg"}],
                        "catalogNumber": "BAD_POSLON_3",
                        "institutionCode": "NCU",
                    },
                    # State-coordinate mismatch (says North Carolina, coordinates in Utah) -> rejected
                    {
                        "key": 66004,
                        "scientificName": "Packera anonyma",
                        "stateProvince": "North Carolina",
                        "decimalLatitude": 39.5,
                        "decimalLongitude": -111.5,
                        "year": 2015,
                        "media": [{"identifier": "https://example.com/utah.jpg", "type": "StillImage", "format": "image/jpeg"}],
                        "catalogNumber": "BAD_STATEMISMATCH_4",
                        "institutionCode": "NCU",
                    },
                    # Cultivated botanical garden sheet -> completely excluded
                    {
                        "key": 66005,
                        "scientificName": "Packera anonyma",
                        "stateProvince": "North Carolina",
                        "decimalLatitude": 35.9,
                        "decimalLongitude": -79.03,
                        "locality": "Cultivated in Sarah P. Duke Gardens ornamental bed",
                        "year": 2018,
                        "media": [{"identifier": "https://example.com/cult.jpg", "type": "StillImage", "format": "image/jpeg"}],
                        "catalogNumber": "BAD_CULT_5",
                        "institutionCode": "DUKE",
                    },
                    # (0, 0) error AND missing stateProvince -> completely excluded (cannot place biogeographically)
                    {
                        "key": 66006,
                        "scientificName": "Packera anonyma",
                        "stateProvince": "",
                        "decimalLatitude": 0.0,
                        "decimalLongitude": 0.0,
                        "year": 2015,
                        "media": [{"identifier": "https://example.com/zeronostate.jpg", "type": "StillImage", "format": "image/jpeg"}],
                        "catalogNumber": "BAD_ZERO_NOSTATE_6",
                        "institutionCode": "NCU",
                    },
                ],
                "count": 6,
            }

        with patch("scripts.core.harvester.occ.search", side_effect=_fake_occ_hygiene):
            harv_h = VoucherHarvester(taxa=["Packera anonyma"], expand_synonyms=False)
            df_h = harv_h.harvest().set_index("catalogNumber")
        self.assertNotIn("BAD_CULT_5", df_h.index)
        self.assertNotIn("BAD_ZERO_NOSTATE_6", df_h.index)
        self.assertAlmostEqual(float(df_h.loc["NCU_VALID_1", "decimalLatitude"]), 35.9, places=4)
        self.assertAlmostEqual(float(df_h.loc["NCU_VALID_1", "decimalLongitude"]), -79.05, places=4)
        self.assertEqual(df_h.loc["NCU_VALID_1", "georef_qc_status"], "valid")
        for bad_cat, exp_qc in [
            ("BAD_ZERO_2", "zero_coord_nullified"),
            ("BAD_POSLON_3", "hemisphere_sign_mismatch_nullified"),
            ("BAD_STATEMISMATCH_4", "state_polygon_mismatch_nullified"),
        ]:
            self.assertIn(bad_cat, df_h.index)
            self.assertTrue(pd.isna(df_h.loc[bad_cat, "decimalLatitude"]))
            self.assertTrue(pd.isna(df_h.loc[bad_cat, "decimalLongitude"]))
            self.assertEqual(df_h.loc[bad_cat, "georef_qc_status"], exp_qc)

        # 4. Secondary duplicate media rescue when primary URL fails download/QC
        with tempfile.TemporaryDirectory() as tmp_d:
            tmp_p = Path(tmp_d)
            good_src = tmp_p / "secondary_good.jpg"
            Image.new("RGB", (3200, 4000), color=(110, 130, 105)).save(good_src, format="JPEG", quality=95)

            harv_m = VoucherHarvester(
                taxa=["Packera anonyma"],
                workspace_dir=tmp_p,
                min_megapixels=4.0,
                min_file_size_kb=5,
                expand_synonyms=False,
            )
            harv_m.workspace_dir = tmp_p
            harv_m._secondary_duplicate_pool = pd.DataFrame([
                {
                    "catalogNumber": "MO_SEC_01",
                    "institutionCode": "MO",
                    "scientificName": "Packera anonyma",
                    "recordedBy": "T.M. Barkley",
                    "recordNumber": "500",
                    "eventDate": "1980-05-10",
                    "determiner_tier": "Tier_1_Gold",
                    "exsiccatae_key": "barkley_500_1980-05-10",
                    "image_path": "raw_vouchers/MO_SEC_01.webp",
                    "_image_url": "https://example.com/secondary_good.jpg",
                }
            ])
            df_prim = pd.DataFrame([
                {
                    "catalogNumber": "NCU_PRIM_01",
                    "institutionCode": "NCU",
                    "scientificName": "Packera anonyma",
                    "recordedBy": "T.M. Barkley",
                    "recordNumber": "500",
                    "eventDate": "1980-05-10",
                    "determiner_tier": "Tier_1_Gold",
                    "exsiccatae_key": "barkley_500_1980-05-10",
                    "image_path": "raw_vouchers/NCU_PRIM_01.webp",
                    "_image_url": "https://example.com/broken_primary.jpg",
                }
            ])

            async def _fake_dl_all(records_to_download, **kwargs):
                is_rescue = any("secondary_good" in str(u) for u, _ in records_to_download)
                if not is_rescue:
                    return {"success": 0, "skipped": 0, "failed": 1}
                dest = Path(records_to_download[0][1])
                dest.parent.mkdir(parents=True, exist_ok=True)
                Image.open(good_src).save(dest, format="WEBP", quality=92)
                return {"success": 1, "skipped": 0, "failed": 0}

            with patch("scripts.core.harvester.download_all_voucher_images", side_effect=_fake_dl_all):
                df_rescued, rescue_stats = harv_m.download_and_validate_media(df_prim)
            self.assertEqual(len(df_rescued), 1)
            self.assertEqual(df_rescued.iloc[0]["catalogNumber"], "MO_SEC_01")
            self.assertTrue(bool(df_rescued.iloc[0]["is_primary_duplicate"]))
            self.assertEqual(rescue_stats["quality_rejected"], 1)
            self.assertTrue((tmp_p / "raw_vouchers" / "MO_SEC_01.webp").exists())


if __name__ == "__main__":
    unittest.main()
