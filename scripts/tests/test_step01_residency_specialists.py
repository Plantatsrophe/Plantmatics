#!/usr/bin/env python3
"""Unit tests for Step 01 Eastern NA residency filtering, type rescue, duplicate dedup, and WebP migration."""

import tempfile
import unittest
from pathlib import Path
from unittest.mock import MagicMock, patch

import pandas as pd
from PIL import Image

from scripts.core.harvester import VoucherHarvester


class TestPhase1ResidencyAndWebpMigration(unittest.TestCase):
    """Tests for Eastern NA residency filtering, type rescue, secondary duplicate dedup, and WebP migration."""

    def test_genus_wide_eastern_na_residency_filter_and_type_rescue(self):
        """Verify genus-wide eastern_na residency check filters out western stray species (< 3 eastern records) and rescues un-georeferenced types."""
        mock_species = MagicMock()
        mock_species.name_backbone.return_value = {
            "usageKey": 3109225,
            "canonicalName": "Packera",
            "rank": "GENUS",
            "matchType": "EXACT",
        }
        mock_species.name_lookup.return_value = {
            "results": [
                {"canonicalName": "Packera aurea", "rank": "SPECIES", "taxonomicStatus": "ACCEPTED"},
                {"canonicalName": "Packera bolanderi", "rank": "SPECIES", "taxonomicStatus": "ACCEPTED"},
            ]
        }
        mock_species.name_usage.return_value = {"results": []}

        def _fake_occ_residency(**kwargs):
            sci = kwargs.get("scientificName")
            if kwargs.get("typeStatus") == "*":
                if sci == "Packera aurea":
                    # 19th-century Holotype lacking decimalLongitude
                    return {
                        "results": [
                            {
                                "key": 88001,
                                "scientificName": "Packera aurea (L.) Á.Löve & D.Löve",
                                "typeStatus": "Holotype",
                                "stateProvince": "Pennsylvania",
                                "decimalLatitude": None,
                                "decimalLongitude": None,
                                "year": 1814,
                                "recordedBy": "F. Pursh",
                                "recordNumber": "1",
                                "media": [{"identifier": "https://example.com/aurea_holo.jpg", "type": "StillImage", "format": "image/jpeg"}],
                                "catalogNumber": "PH_AUREA_HOLOTYPE",
                                "institutionCode": "PH",
                            }
                        ],
                        "count": 1,
                    }
                return {"results": [], "count": 0}

            ctry = kwargs.get("country")
            if sci == "Packera aurea":
                if ctry == "CA":
                    recs = [
                        {
                            "key": 8150,
                            "scientificName": "Packera aurea",
                            "stateProvince": "Ontario",
                            "decimalLatitude": 44.2,
                            "decimalLongitude": -79.5,
                            "coordinateUncertaintyInMeters": 50.0,
                            "year": 2020,
                            "month": 6,
                            "day": 12,
                            "recordedBy": "C. Collector",
                            "recordNumber": "701",
                            "media": [{"identifier": "https://example.com/aurea_on.jpg", "type": "StillImage", "format": "image/jpeg"}],
                            "catalogNumber": "CAN_AUREA_ON",
                            "institutionCode": "CAN",
                        }
                    ]
                    return {"results": recs, "count": len(recs)}
                # 2 US eastern records + 1 ON CA record = 3 Eastern NA records (>= 3 residency threshold),
                # plus 1 western Montana conspecific retained under full continental range harvesting!
                recs = [
                    {
                        "key": 8100 + i,
                        "scientificName": "Packera aurea",
                        "stateProvince": "North Carolina",
                        "decimalLatitude": 35.8 + i * 0.1,
                        "decimalLongitude": -79.0 - i * 0.1,
                        "coordinateUncertaintyInMeters": 50.0,
                        "year": 2020,
                        "month": 5,
                        "day": 10 + i,
                        "recordedBy": "A. Collector",
                        "recordNumber": str(500 + i),
                        "media": [{"identifier": f"https://example.com/aurea_{i}.jpg", "type": "StillImage", "format": "image/jpeg"}],
                        "catalogNumber": f"NCU_AUREA_{i}",
                        "institutionCode": "NCU",
                    }
                    for i in range(2)
                ] + [
                    {
                        "key": 8199,
                        "scientificName": "Packera aurea",
                        "stateProvince": "Montana",
                        "decimalLatitude": 46.8,
                        "decimalLongitude": -110.5,
                        "coordinateUncertaintyInMeters": 100.0,
                        "year": 2019,
                        "month": 7,
                        "day": 4,
                        "recordedBy": "W. Collector",
                        "recordNumber": "902",
                        "media": [{"identifier": "https://example.com/aurea_mt.jpg", "type": "StillImage", "format": "image/jpeg"}],
                        "catalogNumber": "MONT_AUREA_WEST",
                        "institutionCode": "MONT",
                    }
                ]
                return {"results": recs, "count": len(recs)}
            else:
                if ctry == "CA":
                    return {"results": [], "count": 0}
                # Western species with only 2 eastern stray records (< 3 threshold) + 2 Oregon records -> filtered out
                recs = [
                    {
                        "key": 8200 + i,
                        "scientificName": "Packera bolanderi",
                        "stateProvince": "Texas",
                        "decimalLatitude": 31.0 + i * 0.1,
                        "decimalLongitude": -97.0 - i * 0.1,
                        "coordinateUncertaintyInMeters": 50.0,
                        "year": 2020,
                        "month": 5,
                        "day": 10 + i,
                        "recordedBy": "B. Collector",
                        "recordNumber": str(600 + i),
                        "media": [{"identifier": f"https://example.com/bolanderi_{i}.jpg", "type": "StillImage", "format": "image/jpeg"}],
                        "catalogNumber": f"TEX_BOLANDERI_{i}",
                        "institutionCode": "TEX",
                    }
                    for i in range(2)
                ] + [
                    {
                        "key": 8250,
                        "scientificName": "Packera bolanderi",
                        "stateProvince": "Oregon",
                        "decimalLatitude": 44.0,
                        "decimalLongitude": -123.1,
                        "coordinateUncertaintyInMeters": 50.0,
                        "year": 2020,
                        "month": 6,
                        "day": 1,
                        "recordedBy": "P. Collector",
                        "recordNumber": "650",
                        "media": [{"identifier": "https://example.com/bolanderi_or.jpg", "type": "StillImage", "format": "image/jpeg"}],
                        "catalogNumber": "OSC_BOLANDERI_OR",
                        "institutionCode": "OSC",
                    }
                ]
                return {"results": recs, "count": len(recs)}

        with patch("scripts.core.harvester.species", mock_species), \
             patch("scripts.core.harvester.occ.search", side_effect=_fake_occ_residency):
            # Default eastern_na drops Montana (-110.5W) and non-cis-Mississippi P. bolanderi while keeping cis-Mississippi + Holotype
            harvester = VoucherHarvester(taxa=[], max_records_per_taxon=10, expand_synonyms=False)
            df = harvester.harvest()
            self.assertEqual(set(df["species_standardized"]), {"Packera aurea"})
            self.assertNotIn("Packera bolanderi", set(df["species_standardized"]))
            self.assertIn("PH_AUREA_HOLOTYPE", set(df["catalogNumber"]))
            self.assertIn("CAN_AUREA_ON", set(df["catalogNumber"]))
            self.assertNotIn("MONT_AUREA_WEST", set(df["catalogNumber"]))

            # Explicit harvest_full_continental_range=True opt-in still retains Montana conspecific
            harvester_full = VoucherHarvester(
                taxa=[],
                max_records_per_taxon=10,
                expand_synonyms=False,
                harvest_full_continental_range=True,
            )
            df_full = harvester_full.harvest()
            self.assertIn("MONT_AUREA_WEST", set(df_full["catalogNumber"]))

    def test_harvest_deduplicates_secondary_exsiccatae_and_migrates_jpg_to_webp(self):
        """Verify harvest() deduplicates secondary exsiccatae before download and migrate_existing_vouchers_to_webp converts .jpg to .webp."""
        from scripts.core.harvester import migrate_existing_vouchers_to_webp, resolve_voucher_image_path

        def _fake_occ_dups(**kwargs):
            if "typeStatus" in kwargs:
                return {"results": [], "count": 0}
            # Two duplicate sheets of the exact same collection event (Fuller #777 on 2021-05-01)
            return {
                "results": [
                    {
                        "key": 9101,
                        "scientificName": "Packera dubia",
                        "stateProvince": "North Carolina",
                        "decimalLatitude": 35.9,
                        "decimalLongitude": -79.0,
                        "coordinateUncertaintyInMeters": 50.0,
                        "year": 2021,
                        "month": 5,
                        "day": 1,
                        "recordedBy": "J. B. Fuller",
                        "recordNumber": "777",
                        "media": [{"identifier": "https://example.com/ncu_dup.jpg", "type": "StillImage", "format": "image/jpeg"}],
                        "catalogNumber": "NCU_PRIMARY_777",
                        "institutionCode": "NCU",
                    },
                    {
                        "key": 9102,
                        "scientificName": "Packera dubia",
                        "stateProvince": "North Carolina",
                        "decimalLatitude": 35.9,
                        "decimalLongitude": -79.0,
                        "coordinateUncertaintyInMeters": 50.0,
                        "year": 2021,
                        "month": 5,
                        "day": 1,
                        "recordedBy": "J. Brandon Fuller",
                        "recordNumber": "#777",
                        "media": [{"identifier": "https://example.com/ga_dup.jpg", "type": "StillImage", "format": "image/jpeg"}],
                        "catalogNumber": "GA_SECONDARY_777",
                        "institutionCode": "GA",
                    },
                ],
                "count": 2,
            }

        with patch("scripts.core.harvester.occ.search", side_effect=_fake_occ_dups):
            harvester = VoucherHarvester(taxa=["Packera dubia"], max_records_per_taxon=10, expand_synonyms=False)
            df_dedup = harvester.harvest()

        self.assertEqual(len(df_dedup), 1)
        self.assertEqual(df_dedup.iloc[0]["catalogNumber"], "NCU_PRIMARY_777")
        self.assertEqual(df_dedup.iloc[0]["image_path"], "data/raw_vouchers/NCU_PRIMARY_777.webp")

        with tempfile.TemporaryDirectory() as tmp_d:
            tmp_p = Path(tmp_d)
            raw_dir = tmp_p / "data" / "raw_vouchers"
            raw_dir.mkdir(parents=True, exist_ok=True)
            jpg_file = raw_dir / "NCU_LEGACY_001.jpg"
            Image.new("RGB", (600, 800), color=(90, 140, 100)).save(jpg_file, format="JPEG")

            csv_file = tmp_p / "curated_vouchers.csv"
            pd.DataFrame([
                {
                    "catalogNumber": "NCU_LEGACY_001",
                    "image_path": "data/raw_vouchers/NCU_LEGACY_001.jpg",
                }
            ]).to_csv(csv_file, index=False)

            migrated_n = migrate_existing_vouchers_to_webp(raw_dir, curated_csv=csv_file)
            self.assertEqual(migrated_n, 1)
            self.assertFalse(jpg_file.exists())
            webp_file = raw_dir / "NCU_LEGACY_001.webp"
            self.assertTrue(webp_file.is_file())
            self.assertGreater(webp_file.stat().st_size, 0)

            updated_csv = pd.read_csv(csv_file)
            self.assertEqual(updated_csv.iloc[0]["image_path"], "data/raw_vouchers/NCU_LEGACY_001.webp")
            resolved = resolve_voucher_image_path("data/raw_vouchers/NCU_LEGACY_001.jpg", workspace_dir=tmp_p)
            self.assertEqual(resolved, webp_file)


if __name__ == "__main__":
    unittest.main()
