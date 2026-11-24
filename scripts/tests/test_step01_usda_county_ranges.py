#!/usr/bin/env python3
"""Unit tests for Step 01 USDA PLANTS symbol resolution, SQLite caching, and buffered county range filtering."""

import tempfile
import unittest
from pathlib import Path
from unittest.mock import MagicMock, patch

import pandas as pd

from scripts.core.harvester import VoucherHarvester


class TestPhase1UsdaCountyRanges(unittest.TestCase):
    """Tests for USDA PLANTS symbol resolution, CSV parsing, SQLite caching, and buffered county range filtering."""

    def test_usda_plants_symbol_resolution_csv_parsing_and_sqlite_caching(self):
        """Verify USDA PLANTS synonym symbol resolution (Packera dubia -> PATO4), CSV parsing, Excluded filtering, and SQLite caching."""
        from scripts.core.usda_plants_ranges import (
            fetch_or_load_usda_taxon_range,
            parse_usda_distribution_csv,
            reset_usda_network_state,
            resolve_usda_symbols_for_taxon,
        )

        raw_csv = (
            '"Symbol","Scientific Name","Common Name","State","County","FIPS Code","Native Status"\n'
            '"PATO4","Packera tomentosa (Michx.) C. Jeffrey","woolly ragwort","NC","Orange","37135","L48 (N)"\n'
            '"PATO4","Packera tomentosa (Michx.) C. Jeffrey","woolly ragwort","NC","Wake","37183","L48 (N)"\n'
            '"PATO4","Packera tomentosa (Michx.) C. Jeffrey","woolly ragwort","DE","","10","L48 (N)"\n'
            '"PATO4","Packera tomentosa (Michx.) C. Jeffrey","woolly ragwort","NY","Essex","36031","L48 (X) Excluded"\n'
        )
        parsed = parse_usda_distribution_csv(raw_csv, default_symbol="PATO4")
        self.assertEqual(parsed["counties"], {"north carolina|orange", "north carolina|wake"})
        self.assertEqual(parsed["state_only_states"], {"delaware"})
        self.assertNotIn("new york", parsed["documented_states"])

        with tempfile.TemporaryDirectory() as tmp_d:
            tmp_p = Path(tmp_d)
            cache_dir = tmp_p / "usda_cache"
            cache_dir.mkdir(parents=True, exist_ok=True)
            db_file = tmp_p / "usda_test.sqlite"

            # 1. Resolution via local plantlst.txt
            plantlst_sample = (
                '"Symbol","Synonym Symbol","Scientific Name with Author","Common Name","Family"\n'
                '"PATO4","","Packera tomentosa (Michx.) C. Jeffrey","woolly ragwort","Asteraceae"\n'
                '"PATO4","SETO2","Senecio tomentosus Michx.","","Asteraceae"\n'
                '"PASE16","","Packera serpenticola Boufford, Kartesz, S.H. Shi & Renchao Zhou","serpentine ragwort","Asteraceae"\n'
            )
            (cache_dir / "plantlst.txt").write_text(plantlst_sample, encoding="utf-8")
            mock_dubia_sp = MagicMock()
            mock_dubia_sp.name_backbone.return_value = {
                "matchType": "EXACT",
                "status": "ACCEPTED",
                "usageKey": 5387757,
                "canonicalName": "Packera dubia",
                "species": "Packera dubia",
            }
            mock_dubia_sp.name_usage.return_value = {
                "results": [
                    {"canonicalName": "Packera tomentosa", "scientificName": "Packera tomentosa (Michx.) C. Jeffrey", "taxonomicStatus": "SYNONYM"},
                    {"canonicalName": "Senecio tomentosus", "scientificName": "Senecio tomentosus Michx.", "taxonomicStatus": "SYNONYM"},
                ]
            }
            with patch("scripts.core.harvester.species", mock_dubia_sp):
                syms = resolve_usda_symbols_for_taxon(
                    "Packera dubia",
                    cache_dir=cache_dir,
                    allow_network=False,
                )
            self.assertEqual(syms, ["PATO4"])

            # 2. Fetch via mocked HTTP endpoint, verify CSV + SQLite persistence, then reload offline from SQLite
            reset_usda_network_state()
            with patch("scripts.core.harvester.species", mock_dubia_sp), \
                 patch("scripts.core.usda_plants_ranges._http_get_text", return_value=raw_csv) as mock_get:
                r_live = fetch_or_load_usda_taxon_range(
                    canonical_packera="Packera dubia",
                    cache_dir=cache_dir,
                    db_path=db_file,
                    workspace_dir=tmp_p,
                    allow_network=True,
                )
                self.assertEqual(mock_get.call_count, 1)
                self.assertEqual(r_live["symbols"], ["PATO4"])
                self.assertIn("north carolina|orange", r_live["counties"])
                self.assertIn("delaware", r_live["state_only_states"])
                self.assertTrue((cache_dir / "distributions" / "PATO4.csv").is_file())

            # Delete CSV cache and verify offline reload directly from SQLite usda_county_ranges
            (cache_dir / "distributions" / "PATO4.csv").unlink()
            r_sql = fetch_or_load_usda_taxon_range(
                canonical_packera="Packera dubia",
                cache_dir=cache_dir,
                db_path=db_file,
                workspace_dir=tmp_p,
                allow_network=False,
            )
            self.assertEqual(r_sql["symbols"], ["PATO4"])
            self.assertEqual(r_sql["counties"], {"north carolina|orange", "north carolina|wake"})
            self.assertEqual(r_sql["state_only_states"], {"delaware"})

            # 3. Verify USDA PLANTS 3.0 "State and Province,County Name,Distribution Status" CSV schema
            #    and two-step GET /api/PlantProfile -> POST /api/PlantProfile/getDownloadDistributionDocumentation flow
            v3_csv = (
                "Symbol,Scientific Name,Common Name,State and Province,County Name,FIPS Code,Distribution Status\n"
                "PATO4,Packera tomentosa (Michx.) C. Jeffrey,woolly ragwort,Alabama,Autauga County,01001,Native\n"
                "PATO4,Packera tomentosa (Michx.) C. Jeffrey,woolly ragwort,North Carolina,Durham County,37063,Native\n"
                "PATO4,Packera tomentosa (Michx.) C. Jeffrey,woolly ragwort,Michigan,Washtenaw County,26161,Excluded\n"
            )
            parsed_v3 = parse_usda_distribution_csv(v3_csv, default_symbol="PATO4")
            self.assertEqual(parsed_v3["counties"], {"alabama|autauga", "north carolina|durham"})

            db_v3 = tmp_p / "test_usda_v3.sqlite"
            cache_v3 = tmp_p / "usda_cache_v3"
            (cache_v3 / "plantlst.txt").parent.mkdir(parents=True, exist_ok=True)
            (cache_v3 / "plantlst.txt").write_text(plantlst_sample, encoding="utf-8")
            reset_usda_network_state()
            with patch("scripts.core.harvester.species", mock_dubia_sp), patch(
                "scripts.core.usda_plants_ranges._http_get_text",
                return_value='{"Id": 38530, "Symbol": "PATO4"}',
            ) as m_prof, patch(
                "scripts.core.usda_plants_ranges._http_post_text",
                return_value=v3_csv,
            ) as m_post:
                r_v3 = fetch_or_load_usda_taxon_range(
                    canonical_packera="Packera dubia",
                    cache_dir=cache_v3,
                    db_path=db_v3,
                    workspace_dir=tmp_p,
                    allow_network=True,
                )
                self.assertEqual(m_prof.call_count, 1)
                self.assertEqual(m_post.call_count, 1)
                post_payload_sent = m_post.call_args.kwargs.get("payload") or m_post.call_args[0][1]
                self.assertEqual(post_payload_sent["MasterId"], 38530)
                self.assertEqual(post_payload_sent["Symbol"], "PATO4")
                self.assertEqual(r_v3["counties"], {"alabama|autauga", "north carolina|durham"})

    def test_buffered_usda_county_filter_exemptions_and_harvester_integration(self):
        """Verify 1-ring/geodesic county buffer, Type & Gold exemptions, out-of-range misidentification removal, and 26-column CSV invariant."""
        with tempfile.TemporaryDirectory() as tmp_d:
            tmp_p = Path(tmp_d)
            usda_dir = tmp_p / "usda_plants"
            (usda_dir / "distributions").mkdir(parents=True, exist_ok=True)
            db_file = tmp_p / "harvester_usda.sqlite"
            out_csv = tmp_p / "curated_vouchers.csv"

            # Synthetic 3-county GeoJSON: Orange (documented), Durham (adjacent 1-ring neighbor ~18 km), and Essex NY (~850 km away)
            counties_geojson = {
                "type": "FeatureCollection",
                "features": [
                    {
                        "type": "Feature",
                        "properties": {"state_name": "north carolina", "county_norm": "orange", "centroid_lon": -79.10, "centroid_lat": 36.05},
                        "geometry": {"type": "Polygon", "coordinates": [[[-79.25, 35.90], [-78.95, 35.90], [-78.95, 36.20], [-79.25, 36.20], [-79.25, 35.90]]]},
                    },
                    {
                        "type": "Feature",
                        "properties": {"state_name": "north carolina", "county_norm": "durham", "centroid_lon": -78.85, "centroid_lat": 36.05},
                        "geometry": {"type": "Polygon", "coordinates": [[[-78.95, 35.90], [-78.70, 35.90], [-78.70, 36.20], [-78.95, 36.20], [-78.95, 35.90]]]},
                    },
                    {
                        "type": "Feature",
                        "properties": {"state_name": "new york", "county_norm": "essex", "centroid_lon": -73.90, "centroid_lat": 44.10},
                        "geometry": {"type": "Polygon", "coordinates": [[[-74.20, 43.80], [-73.50, 43.80], [-73.50, 44.40], [-74.20, 44.40], [-74.20, 43.80]]]},
                    },
                ],
            }
            cnty_geojson_path = tmp_p / "us_counties_test.geojson"
            import json
            cnty_geojson_path.write_text(json.dumps(counties_geojson), encoding="utf-8")

            (usda_dir / "plantlst.txt").write_text(
                '"Symbol","Synonym Symbol","Scientific Name with Author","Common Name","Family"\n'
                '"PATO4","","Packera tomentosa (Michx.) C. Jeffrey","woolly ragwort","Asteraceae"\n',
                encoding="utf-8",
            )
            (usda_dir / "distributions" / "PATO4.csv").write_text(
                '"Symbol","Scientific Name","Common Name","State","County","FIPS Code","Native Status"\n'
                '"PATO4","Packera tomentosa","woolly ragwort","NC","Orange","37135","L48 (N)"\n',
                encoding="utf-8",
            )

            mock_gbif_results = [
                # 1. Documented in-range county (Orange Co., NC) -> RETAIN (in_range)
                {
                    "catalogNumber": "NCU-INRANGE-01",
                    "institutionCode": "NCU",
                    "scientificName": "Packera dubia",
                    "recordedBy": "A. Collector",
                    "recordNumber": "101",
                    "identifiedBy": "Local Botanist",
                    "stateProvince": "North Carolina",
                    "county": "Orange",
                    "decimalLatitude": 36.05,
                    "decimalLongitude": -79.10,
                    "coordinateUncertaintyInMeters": 100.0,
                    "eventDate": "2021-04-20",
                    "year": 2021, "month": 4, "day": 20,
                    "media": [{"type": "StillImage", "identifier": "https://example.org/img1.jpg"}],
                },
                # 2. Adjacent county not in USDA checklist ("Swiss-cheese" gap: Durham Co., NC) -> RETAIN (buffered_adjacent)
                {
                    "catalogNumber": "DUKE-ADJ-02",
                    "institutionCode": "DUKE",
                    "scientificName": "Packera dubia",
                    "recordedBy": "B. Collector",
                    "recordNumber": "102",
                    "identifiedBy": "Local Botanist",
                    "stateProvince": "North Carolina",
                    "county": "Durham",
                    "decimalLatitude": 36.05,
                    "decimalLongitude": -78.85,
                    "coordinateUncertaintyInMeters": 100.0,
                    "eventDate": "2021-04-22",
                    "year": 2021, "month": 4, "day": 22,
                    "media": [{"type": "StillImage", "identifier": "https://example.org/img2.jpg"}],
                },
                # 3. Far out-of-range misidentification (Essex Co., NY ~850 km away, Bronze/Silver) -> DROP
                {
                    "catalogNumber": "NY-MISID-03",
                    "institutionCode": "NY",
                    "scientificName": "Packera dubia",
                    "recordedBy": "C. Collector",
                    "recordNumber": "103",
                    "identifiedBy": "Student Determiner",
                    "stateProvince": "New York",
                    "county": "Essex",
                    "decimalLatitude": 44.10,
                    "decimalLongitude": -73.90,
                    "coordinateUncertaintyInMeters": 250.0,
                    "eventDate": "2020-06-15",
                    "year": 2020, "month": 6, "day": 15,
                    "media": [{"type": "StillImage", "identifier": "https://example.org/img3.jpg"}],
                },
                # 4. Out-of-range Tier_1_Gold specialist determination (A.S. Weakley) -> EXEMPT RETAIN
                {
                    "catalogNumber": "NCU-GOLD-EXT-04",
                    "institutionCode": "NCU",
                    "scientificName": "Packera dubia",
                    "recordedBy": "D. Collector",
                    "recordNumber": "104",
                    "identifiedBy": "A.S. Weakley",
                    "stateProvince": "New York",
                    "county": "Essex",
                    "decimalLatitude": 44.12,
                    "decimalLongitude": -73.92,
                    "coordinateUncertaintyInMeters": 100.0,
                    "eventDate": "2022-05-30",
                    "year": 2022, "month": 5, "day": 30,
                    "media": [{"type": "StillImage", "identifier": "https://example.org/img4.jpg"}],
                },
                # 5. Out-of-range Holotype -> EXEMPT RETAIN
                {
                    "catalogNumber": "US-TYPE-05",
                    "institutionCode": "US",
                    "scientificName": "Packera dubia",
                    "typeStatus": "Holotype",
                    "recordedBy": "A. Michaux",
                    "recordNumber": "1",
                    "identifiedBy": "Unknown",
                    "stateProvince": "New York",
                    "county": "Essex",
                    "decimalLatitude": 44.05,
                    "decimalLongitude": -73.88,
                    "coordinateUncertaintyInMeters": 500.0,
                    "eventDate": "1803-05-10",
                    "year": 1803, "month": 5, "day": 10,
                    "media": [{"type": "StillImage", "identifier": "https://example.org/img5.jpg"}],
                },
            ]

            mock_occ = MagicMock()
            mock_occ.search.side_effect = lambda **kw: (
                {"results": [], "endOfRecords": True}
                if "typeStatus" in kw
                else {"results": mock_gbif_results, "endOfRecords": True}
            )

            mock_dubia_sp2 = MagicMock()
            mock_dubia_sp2.name_backbone.return_value = {
                "matchType": "EXACT",
                "status": "ACCEPTED",
                "usageKey": 5387757,
                "canonicalName": "Packera dubia",
                "species": "Packera dubia",
            }
            mock_dubia_sp2.name_usage.return_value = {
                "results": [
                    {"canonicalName": "Packera tomentosa", "scientificName": "Packera tomentosa (Michx.) C. Jeffrey", "taxonomicStatus": "SYNONYM"},
                ]
            }
            with patch("scripts.core.harvester.occ", mock_occ), patch("scripts.core.harvester.species", mock_dubia_sp2):
                harv = VoucherHarvester(
                    taxa=["Packera dubia"],
                    expand_synonyms=False,
                    output_csv=out_csv,
                    raw_dir=tmp_p / "raw",
                    workspace_dir=tmp_p,
                    db_path=db_file,
                    usda_county_filter=True,
                    usda_county_buffer_km=50.0,
                    usda_exempt_gold=True,
                    usda_cache_dir=usda_dir,
                    usda_county_geojson_path=cnty_geojson_path,
                )
                df_res = harv.harvest()
                harv.export(df_res)

            retained_cats = set(df_res["catalogNumber"].tolist())
            self.assertEqual(
                retained_cats,
                {"NCU-INRANGE-01", "DUKE-ADJ-02", "NCU-GOLD-EXT-04", "US-TYPE-05"},
            )
            self.assertNotIn("NY-MISID-03", retained_cats)
            self.assertEqual(harv.last_usda_filter_stats["dropped_out_of_range"], 1)
            self.assertEqual(harv.last_usda_filter_stats["buffered_adjacent"], 1)
            self.assertEqual(harv.last_usda_filter_stats["exempt_retained"], 2)

            # Verify SQLite-only USDA diagnostic columns never leak into curated_vouchers.csv
            csv_df = pd.read_csv(out_csv)
            from scripts.core.harvester import EXPORT_COLUMNS
            self.assertEqual(list(csv_df.columns[: len(EXPORT_COLUMNS)]), EXPORT_COLUMNS)
            self.assertNotIn("usda_symbol", csv_df.columns)
            self.assertNotIn("usda_county_status", csv_df.columns)

            import sqlite3
            conn = sqlite3.connect(str(db_file))
            try:
                rows_sql = dict(
                    conn.execute(
                        "SELECT catalogNumber, usda_county_status FROM vouchers;"
                    ).fetchall()
                )
                cnt_ranges = conn.execute(
                    "SELECT COUNT(*) FROM usda_county_ranges WHERE canonical_packera = 'Packera dubia';"
                ).fetchone()[0]
            finally:
                conn.close()
            self.assertEqual(rows_sql["NCU-INRANGE-01"], "in_range")
            self.assertEqual(rows_sql["DUKE-ADJ-02"], "buffered_adjacent")
            self.assertEqual(rows_sql["NCU-GOLD-EXT-04"], "out_of_range")
            self.assertEqual(rows_sql["US-TYPE-05"], "out_of_range")
            self.assertEqual(int(cnt_ranges), 1)


if __name__ == "__main__":
    unittest.main()
