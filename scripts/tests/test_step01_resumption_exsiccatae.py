#!/usr/bin/env python3
"""Unit tests for Step 01 download resumption, atomic CSV merge, and exsiccatae stratification."""

import tempfile
import unittest
from pathlib import Path
from unittest.mock import MagicMock, patch

import pandas as pd

from scripts.core.harvester import export_curated_table
from scripts.core.taxonomy_authority import (
    generate_exsiccatae_key,
    haversine_distance_meters,
    normalize_collector,
    normalize_event_date,
    normalize_record_number,
    stratify_duplicates,
)


class TestVoucherHarvesterResumptionLogic(unittest.TestCase):
    """Test suite for download resumption, atomic downloads, and CSV merging."""

    def setUp(self):
        self.temp_dir = tempfile.TemporaryDirectory()
        self.dir_path = Path(self.temp_dir.name)

    def tearDown(self):
        self.temp_dir.cleanup()

    def test_download_single_image_resumption_skips_existing(self):
        """Verify that an existing non-empty voucher image is retained and skipped when force=False."""
        import asyncio
        from scripts.core.harvester import download_single_image

        dest_file = self.dir_path / "NCU000123.jpg"
        dest_file.write_bytes(b"\xff\xd8\xff" + b"A" * 2048)

        mock_session = MagicMock()
        mock_logger = MagicMock()
        semaphore = asyncio.Semaphore(1)

        success, was_skipped = asyncio.run(
            download_single_image(
                session=mock_session,
                image_url="https://example.com/sheet.jpg",
                destination_path=dest_file,
                semaphore=semaphore,
                force=False,
                logger=mock_logger,
            )
        )

        self.assertTrue(success)
        self.assertTrue(was_skipped)
        mock_session.get.assert_not_called()
        mock_logger.info.assert_called()
        self.assertIn("Skipping download for existing voucher", mock_logger.info.call_args[0][0])

    def test_download_single_image_force_redownloads(self):
        """Verify that when force=True, existing files are re-downloaded and replaced atomically."""
        import asyncio
        from scripts.core.harvester import download_single_image

        dest_file = self.dir_path / "NCU000456.jpg"
        dest_file.write_bytes(b"old_data")
        new_image_data = b"\xff\xd8\xff" + b"X" * 2048

        # Mock aiohttp response context manager
        mock_response = MagicMock()
        mock_response.status = 200
        mock_response.headers = {"Content-Type": "image/jpeg"}

        async def mock_read():
            return new_image_data
        mock_response.read = mock_read

        class MockGetContext:
            async def __aenter__(self):
                return mock_response
            async def __aexit__(self, exc_type, exc_val, exc_tb):
                pass

        mock_session = MagicMock()
        mock_session.get.return_value = MockGetContext()
        semaphore = asyncio.Semaphore(1)

        success, was_skipped = asyncio.run(
            download_single_image(
                session=mock_session,
                image_url="https://example.com/sheet.jpg",
                destination_path=dest_file,
                semaphore=semaphore,
                force=True,
            )
        )

        self.assertTrue(success)
        self.assertFalse(was_skipped)
        mock_session.get.assert_called_once()
        self.assertEqual(dest_file.read_bytes(), new_image_data)

    def test_download_all_voucher_images_counters(self):
        """Verify download_all_voucher_images aggregates progress counters [Processed: X | Skipped: Y | Total: Z]."""
        import asyncio
        from scripts.core.harvester import download_all_voucher_images

        img1 = self.dir_path / "NCU_EXISTING.jpg"
        img1.write_bytes(b"\xff\xd8\xff" + b"Z" * 1500)
        img2 = self.dir_path / "NCU_PENDING.jpg"

        records = [
            ("https://example.com/img1.jpg", img1),
            ("https://example.com/img2.jpg", img2),
        ]

        mock_logger = MagicMock()
        with patch("scripts.core.harvester.download_single_image", return_value=(True, False)):
            stats = asyncio.run(
                download_all_voucher_images(
                    records_to_download=records,
                    concurrency_limit=2,
                    force=False,
                    logger=mock_logger,
                )
            )

        self.assertEqual(stats["total"], 2)
        self.assertEqual(stats["skipped"], 1)
        self.assertEqual(stats["success"], 1)
        mock_logger.info.assert_called()

    def test_export_curated_table_merges_without_duplicates(self):
        """Verify export_curated_table merges new records with existing CSV without duplicate catalogNumber rows."""
        csv_path = self.dir_path / "curated_vouchers.csv"

        df_initial = pd.DataFrame([
            {"catalogNumber": "NCU001", "scientificName": "Packera dubia", "year": 2020},
            {"catalogNumber": "NCU002", "scientificName": "Packera dubia", "year": 2021},
        ])
        export_curated_table(df_initial, csv_path)
        self.assertEqual(len(pd.read_csv(csv_path)), 2)

        # Second harvest includes NCU002 (duplicate) and NCU003 (new)
        df_new = pd.DataFrame([
            {"catalogNumber": "NCU002", "scientificName": "Packera dubia", "year": 2022},
            {"catalogNumber": "NCU003", "scientificName": "Packera dubia", "year": 2023},
        ])
        export_curated_table(df_new, csv_path)

        df_merged = pd.read_csv(csv_path)
        self.assertEqual(len(df_merged), 3, "Merged table must contain exactly 3 unique vouchers.")
        self.assertListEqual(list(df_merged["catalogNumber"]), ["NCU001", "NCU002", "NCU003"])
        # Should keep last updated value for NCU002
        ncu002_year = df_merged[df_merged["catalogNumber"] == "NCU002"]["year"].iloc[0]
        self.assertEqual(ncu002_year, 2022)


class TestExsiccataeFingerprintingAndStratification(unittest.TestCase):
    """Test suite for exsiccatae collection event fingerprinting and multi-tier duplicate stratification."""

    def test_normalize_collector(self):
        """Verify collector name normalization: lowercase, strip initials, extract primary surname."""
        test_cases = [
            ("J. Brandon Fuller", "fuller"),
            ("J.B. Fuller", "fuller"),
            ("Fuller, J. Brandon", "fuller"),
            ("Fuller, J. B.", "fuller"),
            ("Fuller", "fuller"),
            ("A. Cronquist", "cronquist"),
            ("R. R. Kowal & J. B. Fuller", "kowal"),
            ("Robert R. Kowal with J. Brandon Fuller", "kowal"),
            ("Radford, Ahles & Bell", "radford"),
            ("J. B.", "jb"),
            ("Unknown", "unknown"),
            ("", "unknown"),
            (None, "unknown"),
        ]
        for raw_name, expected in test_cases:
            with self.subTest(raw_name=raw_name):
                self.assertEqual(normalize_collector(raw_name), expected)

    def test_normalize_record_number(self):
        """Verify collection number normalization: strip prefixes/suffixes, extract numeric component."""
        test_cases = [
            ("#1042", "1042"),
            ("1042", "1042"),
            ("No. 1042", "1042"),
            ("No. 1042b", "1042b"),
            ("1042-A", "1042a"),
            ("2022-014", "2022-014"),
            ("2022-015", "2022-015"),
            ("Coll #5821", "5821"),
            ("s.n.", ""),
            ("s. n.", ""),
            ("sn", ""),
            ("unnumbered", ""),
            ("", ""),
            (None, ""),
        ]
        for raw_num, expected in test_cases:
            with self.subTest(raw_num=raw_num):
                self.assertEqual(normalize_record_number(raw_num), expected)

    def test_normalize_event_date(self):
        """Verify eventDate normalization to standard YYYY-MM-DD format."""
        test_cases = [
            ({"event_date_raw": "1984-05-12"}, "1984-05-12"),
            ({"event_date_raw": "1984-05-12T00:00:00Z"}, "1984-05-12"),
            ({"event_date_raw": "1984/05/12"}, "1984-05-12"),
            ({"year": 1984, "month": 5, "day": 12}, "1984-05-12"),
            ({"year": 1984}, "1984-00-00"),
            ({"event_date_raw": "1984"}, "1984-00-00"),
            ({}, "unknown_date"),
        ]
        for kwargs, expected in test_cases:
            with self.subTest(kwargs=kwargs):
                self.assertEqual(normalize_event_date(**kwargs), expected)

    def test_generate_exsiccatae_key(self):
        """Verify construction of standardized signature string f'{norm_collector}_{norm_number}_{norm_date}'."""
        key = generate_exsiccatae_key(
            collector="J. Brandon Fuller",
            number="#1042",
            date="1984-05-12",
        )
        self.assertEqual(key, "fuller_1042_1984-05-12")

    def test_haversine_distance_meters(self):
        """Verify haversine distance calculates accurate geodesic distances."""
        # Exact same point
        self.assertAlmostEqual(haversine_distance_meters(35.9123, -79.0512, 35.9123, -79.0512), 0.0, places=1)
        # Shift latitude by ~0.0003 deg (~33 m < 50 m)
        dist_near = haversine_distance_meters(35.9123, -79.0512, 35.9126, -79.0512)
        self.assertTrue(25.0 < dist_near < 40.0)
        # Shift latitude by ~0.001 deg (~111 m > 50 m)
        dist_far = haversine_distance_meters(35.9123, -79.0512, 35.9133, -79.0512)
        self.assertTrue(100.0 < dist_far < 125.0)

    def test_stratify_duplicates_numbered_exsiccatae(self):
        """Verify grouping of duplicate sheets sharing collector, collection number, and date."""
        df_exsiccatae = pd.DataFrame([
            {
                "catalogNumber": "NCU001",
                "institutionCode": "NCU",
                "recordedBy": "J. Brandon Fuller",
                "recordNumber": "1042",
                "eventDate": "1984-05-12",
                "determiner_tier": "Tier_2_Silver",
                "type_status": "None",
            },
            {
                "catalogNumber": "GA002",
                "institutionCode": "GA",
                "recordedBy": "J.B. Fuller",
                "recordNumber": "#1042",
                "eventDate": "1984-05-12",
                "determiner_tier": "Tier_3_Bronze",
                "type_status": "None",
            },
            {
                "catalogNumber": "US003",
                "institutionCode": "US",
                "recordedBy": "Fuller, J. Brandon",
                "recordNumber": "No. 1042",
                "eventDate": "1984-05-12",
                "determiner_tier": "Tier_3_Bronze",
                "type_status": "None",
            },
        ])
        strat_df = stratify_duplicates(df_exsiccatae)

        # All 3 records must share identical exsiccatae_key
        self.assertEqual(len(strat_df["exsiccatae_key"].unique()), 1)
        self.assertEqual(strat_df["exsiccatae_key"].iloc[0], "fuller_1042_1984-05-12")

        # duplicate_count must equal 3 for all records
        self.assertTrue((strat_df["duplicate_count"] == 3).all())

        # Exactly 1 primary duplicate, remaining are False
        self.assertEqual(strat_df["is_primary_duplicate"].sum(), 1)
        # Priority rule: NCU beats GA and US (Tier 2 vs Tier 3, and NCU is Home Herbarium)
        ncu_row = strat_df[strat_df["catalogNumber"] == "NCU001"].iloc[0]
        self.assertTrue(ncu_row["is_primary_duplicate"])
        self.assertFalse(strat_df[strat_df["catalogNumber"] == "GA002"]["is_primary_duplicate"].iloc[0])
        self.assertFalse(strat_df[strat_df["catalogNumber"] == "US003"]["is_primary_duplicate"].iloc[0])

    def test_stratify_duplicates_coordinate_proximity_fallback(self):
        """Verify unnumbered collections (< 50m) on identical dates are clustered as duplicates."""
        df_proximity = pd.DataFrame([
            {
                "catalogNumber": "SHEET_A",
                "institutionCode": "NCU",
                "recordedBy": "J. Brandon Fuller",
                "recordNumber": "s.n.",
                "eventDate": "1984-05-12",
                "decimalLatitude": 35.91230,
                "decimalLongitude": -79.05120,
                "determiner_tier": "Tier_2_Silver",
            },
            {
                "catalogNumber": "SHEET_B",
                "institutionCode": "GA",
                "recordedBy": "Fuller",
                "recordNumber": "",
                "eventDate": "1984-05-12",
                "decimalLatitude": 35.91245,  # ~17 m away (< 50 m)
                "decimalLongitude": -79.05125,
                "determiner_tier": "Tier_3_Bronze",
            },
            {
                "catalogNumber": "SHEET_C",
                "institutionCode": "MO",
                "recordedBy": "Fuller",
                "recordNumber": "s.n.",
                "eventDate": "1984-05-12",
                "decimalLatitude": 36.50000,  # ~65 km away (> 50 m)
                "decimalLongitude": -79.05120,
                "determiner_tier": "Tier_3_Bronze",
            },
        ])
        strat_df = stratify_duplicates(df_proximity, max_distance_meters=50.0)

        # SHEET_A and SHEET_B should share the same exsiccatae cluster key
        key_a = strat_df[strat_df["catalogNumber"] == "SHEET_A"]["exsiccatae_key"].iloc[0]
        key_b = strat_df[strat_df["catalogNumber"] == "SHEET_B"]["exsiccatae_key"].iloc[0]
        key_c = strat_df[strat_df["catalogNumber"] == "SHEET_C"]["exsiccatae_key"].iloc[0]

        self.assertEqual(key_a, key_b, "SHEET_A and SHEET_B must be clustered together within 50m proximity.")
        self.assertNotEqual(key_a, key_c, "SHEET_C (>50m away) must not be clustered with SHEET_A/B.")

        # Duplicate counts
        self.assertEqual(strat_df[strat_df["catalogNumber"] == "SHEET_A"]["duplicate_count"].iloc[0], 2)
        self.assertEqual(strat_df[strat_df["catalogNumber"] == "SHEET_B"]["duplicate_count"].iloc[0], 2)
        self.assertEqual(strat_df[strat_df["catalogNumber"] == "SHEET_C"]["duplicate_count"].iloc[0], 1)

        # Primary duplicate flag: NCU (SHEET_A) wins over GA (SHEET_B)
        self.assertTrue(strat_df[strat_df["catalogNumber"] == "SHEET_A"]["is_primary_duplicate"].iloc[0])
        self.assertFalse(strat_df[strat_df["catalogNumber"] == "SHEET_B"]["is_primary_duplicate"].iloc[0])
        # SHEET_C is singleton, so it is primary for its own group
        self.assertTrue(strat_df[strat_df["catalogNumber"] == "SHEET_C"]["is_primary_duplicate"].iloc[0])

    def test_institutional_priority_rule(self):
        """
        Verify institutional hierarchy:
        (1) Verified Type Status -> (2) Home Herbarium (NCU) -> (3) Tier 1 Specialist -> (4) Image resolution.
        """
        # Scenario 1: Verified Type status (GA) beats NCU (non-type)
        df_type_test = pd.DataFrame([
            {
                "catalogNumber": "NCU_NORM",
                "institutionCode": "NCU",
                "recordedBy": "Fuller",
                "recordNumber": "100",
                "eventDate": "1990-06-01",
                "determiner_tier": "Tier_1_Gold",
                "type_status": "None",
            },
            {
                "catalogNumber": "GA_ISOTYPE",
                "institutionCode": "GA",
                "recordedBy": "Fuller",
                "recordNumber": "100",
                "eventDate": "1990-06-01",
                "determiner_tier": "Tier_2_Silver",
                "type_status": "Isotype",
            },
        ])
        res1 = stratify_duplicates(df_type_test)
        self.assertTrue(res1[res1["catalogNumber"] == "GA_ISOTYPE"]["is_primary_duplicate"].iloc[0])
        self.assertFalse(res1[res1["catalogNumber"] == "NCU_NORM"]["is_primary_duplicate"].iloc[0])

        # Scenario 2: Home Herbarium NCU beats non-NCU when neither is a type
        df_ncu_test = pd.DataFrame([
            {
                "catalogNumber": "NCU_SHEET",
                "institutionCode": "NCU",
                "recordedBy": "Fuller",
                "recordNumber": "200",
                "eventDate": "1990-06-01",
                "determiner_tier": "Tier_2_Silver",
                "type_status": "None",
            },
            {
                "catalogNumber": "MO_SHEET",
                "institutionCode": "MO",
                "recordedBy": "Fuller",
                "recordNumber": "200",
                "eventDate": "1990-06-01",
                "determiner_tier": "Tier_2_Silver",
                "type_status": "None",
            },
        ])
        res2 = stratify_duplicates(df_ncu_test)
        self.assertTrue(res2[res2["catalogNumber"] == "NCU_SHEET"]["is_primary_duplicate"].iloc[0])
        self.assertFalse(res2[res2["catalogNumber"] == "MO_SHEET"]["is_primary_duplicate"].iloc[0])

        # Scenario 3: Tier 1 Specialist determination beats Tier 2/3 when neither is NCU and neither is Type
        df_tier1_test = pd.DataFrame([
            {
                "catalogNumber": "GA_TIER1",
                "institutionCode": "GA",
                "recordedBy": "Fuller",
                "recordNumber": "300",
                "eventDate": "1990-06-01",
                "determiner_tier": "Tier_1_Gold",
                "type_status": "None",
            },
            {
                "catalogNumber": "MO_TIER3",
                "institutionCode": "MO",
                "recordedBy": "Fuller",
                "recordNumber": "300",
                "eventDate": "1990-06-01",
                "determiner_tier": "Tier_3_Bronze",
                "type_status": "None",
            },
        ])
        res3 = stratify_duplicates(df_tier1_test)
        self.assertTrue(res3[res3["catalogNumber"] == "GA_TIER1"]["is_primary_duplicate"].iloc[0])
        self.assertFalse(res3[res3["catalogNumber"] == "MO_TIER3"]["is_primary_duplicate"].iloc[0])


if __name__ == "__main__":
    unittest.main()
