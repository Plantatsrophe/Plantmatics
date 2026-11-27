#!/usr/bin/env python3
"""Unit tests for Step 01 spatial balancing/thinning and EPA Level III ecoregion assignment."""

import tempfile
import unittest
from pathlib import Path

import numpy as np
import pandas as pd

from scripts.core.epa_ecoregions import infer_regional_group
from scripts.core.taxonomy_authority import spatially_balance_and_thin_vouchers


class TestSpatialBalancingAndThinning(unittest.TestCase):
    """Test suite for Milestone 1/4 spatial grid thinning, state/region quotas, and inverse-density weights."""

    def test_spatially_balance_and_thin_vouchers_nc_cluster_vs_ar_spread(self):
        """
        Verify that 200 NC vouchers clustered in 2 counties and 15 AR vouchers across 10 counties
        are thinned by spatial_thin_km=15.0, max_per_state=25 to at most 25 NC vouchers while
        retaining all Tier_1_Gold vouchers and computing valid spatial_balance_weight in [0.20, 3.00].
        """
        rng = np.random.default_rng(42)
        records = []

        # 200 NC vouchers clustered across 2 counties (Wake & Franklin), including 6 Tier_1_Gold anchors
        # and 10 secondary duplicate sheets
        gold_catalog_numbers = set()
        for i in range(200):
            is_gold = i < 6
            cat = f"NC_VOUCHER_{i:03d}"
            if is_gold:
                gold_catalog_numbers.add(cat)
            county = "Wake" if i % 2 == 0 else "Franklin"
            # Spread across ~30 grid cells in Wake/Franklin so grid thinning leaves >25 before max_per_state=25 caps at 25
            base_lat = 35.78 + (i % 10) * 0.16
            base_lon = -78.64 + ((i // 10) % 6) * 0.18
            records.append({
                "catalogNumber": cat,
                "institutionCode": "NCU",
                "scientificName": "Packera anonyma",
                "species_raw": "Packera anonyma",
                "determiner_tier": "Tier_1_Gold" if is_gold else ("Tier_2_Silver" if i % 5 == 0 else "Tier_3_Bronze"),
                "type_status": "Holotype" if i == 0 else "None",
                "county": county,
                "stateProvince": "NC",
                "decimalLatitude": float(base_lat + rng.normal(0.0, 0.005)),
                "decimalLongitude": float(base_lon + rng.normal(0.0, 0.005)),
                "coordinateUncertainty": 50.0 if is_gold else 500.0,
                "regional_group": "Piedmont_Granite_Flatrocks",
                "is_primary_duplicate": False if i >= 190 else True,
            })

        # 15 AR vouchers spread across 10 counties, separated by >25 km each
        for j in range(15):
            cat = f"AR_VOUCHER_{j:03d}"
            if j == 0:
                gold_catalog_numbers.add(cat)
            records.append({
                "catalogNumber": cat,
                "institutionCode": "UARK",
                "scientificName": "Packera anonyma",
                "species_raw": "Packera anonyma",
                "determiner_tier": "Tier_1_Gold" if j == 0 else "Tier_2_Silver",
                "type_status": "None",
                "county": f"AR_County_{j % 10}",
                "stateProvince": "AR",
                "decimalLatitude": float(34.20 + (j % 5) * 0.35),
                "decimalLongitude": float(-94.10 + (j // 5) * 0.45),
                "coordinateUncertainty": 100.0,
                "regional_group": "Interior_Prairie_Midwest",
                "is_primary_duplicate": True,
            })

        df_in = pd.DataFrame(records)
        manifest_df = pd.DataFrame([
            {"catalogNumber": "NC_VOUCHER_010", "leaf_id": 1, "efa_qc_pass": True},
            {"catalogNumber": "NC_VOUCHER_010", "leaf_id": 2, "efa_qc_pass": True},
        ])

        thinned = spatially_balance_and_thin_vouchers(
            df_in,
            spatial_thin_km=15.0,
            max_per_state=25,
            primary_duplicates_only=True,
            preserve_all_gold=True,
            manifest_df=manifest_df,
        )

        # 1. Secondary duplicates must be excluded
        self.assertTrue(thinned["is_primary_duplicate"].all())
        for i in range(190, 200):
            self.assertNotIn(f"NC_VOUCHER_{i:03d}", set(thinned["catalogNumber"]))

        # 2. NC vouchers capped at max_per_state=25
        nc_retained = thinned[thinned["stateProvince"] == "NC"]
        ar_retained = thinned[thinned["stateProvince"] == "AR"]
        self.assertLessEqual(len(nc_retained), 25)
        self.assertEqual(len(nc_retained), 25)
        self.assertEqual(len(ar_retained), 15)

        # 3. Both NC counties represented evenly and 100% of Tier_1_Gold vouchers retained
        self.assertEqual(set(nc_retained["county"].unique()), {"Wake", "Franklin"})
        self.assertTrue(gold_catalog_numbers.issubset(set(thinned["catalogNumber"])))

        # 4. Continuous spatial_balance_weight bounded in [0.20, 3.00] and normalized around 1.0
        self.assertIn("spatial_balance_weight", thinned.columns)
        weights = thinned["spatial_balance_weight"].to_numpy(dtype=float)
        self.assertTrue(np.all(np.isfinite(weights)))
        self.assertTrue(np.all((weights >= 0.20) & (weights <= 3.00)))
        self.assertAlmostEqual(float(np.mean(weights)), 1.0, delta=0.15)

        # Under-represented, spatially dispersed AR vouchers should have higher mean weight than dense NC vouchers
        mean_w_ar = float(ar_retained["spatial_balance_weight"].mean())
        mean_w_nc = float(nc_retained["spatial_balance_weight"].mean())
        self.assertGreater(mean_w_ar, mean_w_nc)

    def test_max_per_region_and_cli_arguments(self):
        """Verify max_per_region quota balancing on Tier B thinning."""
        df_reg = pd.DataFrame([
            {
                "catalogNumber": f"REG_R1_{i}",
                "species_standardized": "Packera dubia",
                "determiner_tier": "Tier_2_Silver",
                "regional_group": "Piedmont",
                "stateProvince": "NC" if i % 2 == 0 else "VA",
                "decimalLatitude": 35.5 + i * 0.2,
                "decimalLongitude": -79.0 - i * 0.2,
                "is_primary_duplicate": True,
            }
            for i in range(6)
        ])
        thinned_reg = spatially_balance_and_thin_vouchers(df_reg, max_per_region=3, primary_duplicates_only=True)
        self.assertEqual(len(thinned_reg), 3)


class TestEPALevel3EcoregionAssignment(unittest.TestCase):
    """Test suite for official US EPA Level III Ecoregion polygon assignment & county-containment fallback."""

    def test_point_in_polygon_epa_level3_ecoregion_assignment(self):
        """Verify georeferenced points map strictly by polygon intersection across Atlantic, Gulf, TX/AR, and Interior ecoregions."""
        cases = [
            # (label, lat, lon, expected_ecoregion)
            ("Chapel Hill, NC (Piedmont)", 35.9132, -79.0558, "Piedmont"),
            ("Wilmington, NC (Middle Atlantic Coastal Plain)", 34.2257, -77.9447, "Middle_Atlantic_Coastal_Plain"),
            ("Dublin, GA (Southeastern Plains)", 32.5404, -82.9037, "Southeastern_Plains"),
            ("Gainesville, FL (Southern Coastal Plain)", 29.6516, -82.3248, "Southern_Coastal_Plain"),
            ("Nacogdoches, TX (South Central Plains)", 31.6035, -94.6555, "South_Central_Plains"),
            ("Camden, AR (South Central Plains)", 33.5846, -92.8343, "South_Central_Plains"),
            ("College Station, TX (East Central Texas Plains)", 30.6280, -96.3344, "East_Central_Texas_Plains"),
            ("Greenville, MS (Mississippi Alluvial Plain)", 33.4101, -91.0618, "Mississippi_Alluvial_Plain"),
            ("Asheville, NC (Blue Ridge)", 35.5951, -82.5515, "Blue_Ridge"),
            ("Knoxville, TN (Ridge and Valley)", 35.9606, -83.9207, "Ridge_and_Valley"),
            ("Sewanee, TN (Southwestern Appalachians)", 35.2031, -85.9211, "Southwestern_Appalachians"),
            ("Nashville, TN (Interior Plateau)", 36.1627, -86.7816, "Interior_Plateau"),
            ("Harrison, AR (Ozark Highlands)", 36.2298, -93.1077, "Ozark_Highlands"),
            ("Hot Springs, AR (Ouachita Mountains)", 34.5037, -93.0552, "Ouachita_Mountains"),
        ]
        for desc, lat, lon, expected in cases:
            with self.subTest(desc=desc):
                # Also pass misleading habitat/state text to prove text heuristics and state lumping are ignored
                got = infer_regional_group(
                    lat=lat,
                    lon=lon,
                    state_province="Florida",
                    habitat="sandhills granite flatrock",
                    locality="Blue Ridge mountains",
                )
                self.assertEqual(got, expected, f"Failed for {desc}: got {got!r}, expected {expected!r}")
                self.assertNotEqual(got, "Other_US")

        # Offshore / Null Island points outside all US EPA Level III polygons must return None (NA), never 'Other_US'
        self.assertIsNone(infer_regional_group(0.0, 0.0, state_province="Tennessee", county="Cheatham"))
        self.assertIsNone(infer_regional_group(35.16, -22.28, state_province="South Carolina", county="Greenville"))

    def test_county_containment_fallback_and_straddling_county_na(self):
        """Verify ungeoreferenced single-ecoregion counties resolve via centroid while straddling/missing counties return None (NA)."""
        # Single-ecoregion (>= 95% share) counties
        self.assertEqual(
            infer_regional_group(None, None, state_province="North Carolina", county="Orange County"),
            "Piedmont",
        )
        self.assertEqual(
            infer_regional_group(np.nan, np.nan, state_province="NC", county="New Hanover"),
            "Middle_Atlantic_Coastal_Plain",
        )
        self.assertEqual(
            infer_regional_group(None, None, state_province="TX", county="Nacogdoches Co."),
            "South_Central_Plains",
        )
        self.assertEqual(
            infer_regional_group(None, None, state_province="Arkansas", county="Baxter County"),
            "Ozark_Highlands",
        )
        self.assertEqual(
            infer_regional_group(None, None, state_province="AR", county="Garland"),
            "Ouachita_Mountains",
        )

        # Multi-ecoregion straddling counties (Fall Line & Blue Ridge Escarpment) must return None (NA)
        self.assertIsNone(infer_regional_group(None, None, state_province="North Carolina", county="Wake County"))
        self.assertIsNone(infer_regional_group(None, None, state_province="South Carolina", county="Richland"))
        self.assertIsNone(infer_regional_group(None, None, state_province="Georgia", county="Bibb County"))
        self.assertIsNone(infer_regional_group(None, None, state_province="North Carolina", county="Burke"))

        # Missing, ambiguous, or unmatched counties must return None (NA), never 'Other_US'
        self.assertIsNone(infer_regional_group(None, None, state_province="North Carolina", county=None))
        self.assertIsNone(infer_regional_group(None, None, state_province="Texas", county=""))
        self.assertIsNone(infer_regional_group(None, None, state_province="Tennessee", county="Sevier or Blount"))
        self.assertIsNone(infer_regional_group(None, None, state_province="Ontario", county="Bruce"))

    def test_deterministic_offline_fallback_and_curated_vouchers_integrity(self):
        """Verify deterministic offline polygon mode and zero 'Other_US' in curated_vouchers.csv."""
        with tempfile.TemporaryDirectory() as tmpdir:
            missing_eco = Path(tmpdir) / "nonexistent_eco.geojson"
            missing_cnty = Path(tmpdir) / "nonexistent_cnty.geojson"
            self.assertEqual(
                infer_regional_group(35.91, -79.05, ecoregion_path=missing_eco, county_path=missing_cnty),
                "Piedmont",
            )
            self.assertEqual(
                infer_regional_group(31.60, -94.65, ecoregion_path=missing_eco, county_path=missing_cnty),
                "South_Central_Plains",
            )
            self.assertEqual(
                infer_regional_group(None, None, state_province="NC", county="Orange", ecoregion_path=missing_eco, county_path=missing_cnty),
                "Piedmont",
            )
            self.assertIsNone(
                infer_regional_group(None, None, state_province="NC", county="Wake", ecoregion_path=missing_eco, county_path=missing_cnty),
            )
            # Offline fallback must never write synthetic files to disk
            self.assertFalse(missing_eco.exists())
            self.assertFalse(missing_cnty.exists())

        prod_csv = Path(__file__).resolve().parents[2] / "data" / "tables" / "curated_vouchers.csv"
        if prod_csv.exists():
            df = pd.read_csv(prod_csv, low_memory=False)
            vals = set(df["regional_group"].dropna().astype(str))
            self.assertNotIn("Other_US", vals)
            self.assertNotIn("Coastal_Plain_Sandhills", vals)
            self.assertNotIn("Piedmont_Granite_Flatrocks", vals)
            self.assertIn("Piedmont", vals)
            self.assertIn("Middle_Atlantic_Coastal_Plain", vals)
            self.assertIn("South_Central_Plains", vals)

    def test_habitat_stratified_thinning_preserves_azonal_outcrop_islands(self):
        """Verify 15 km thinning separates azonal granite outcrops from adjacent zonal sand matrix."""
        df_hab = pd.DataFrame([
            {"catalogNumber": "SAND_1", "species_standardized": "Packera dubia", "determiner_tier": "Tier_1_Gold",
             "decimalLatitude": 34.6650, "decimalLongitude": -80.5200, "coordinateUncertainty": 50.0, "lithic_outcrop_index": 0.05},
            {"catalogNumber": "ROCK_1", "species_standardized": "Packera dubia", "determiner_tier": "Tier_1_Gold",
             "decimalLatitude": 34.6690, "decimalLongitude": -80.5273, "coordinateUncertainty": 100.0, "lithic_outcrop_index": 0.82},
            {"catalogNumber": "ROCK_2", "species_standardized": "Packera dubia", "determiner_tier": "Tier_2_Silver",
             "decimalLatitude": 34.6705, "decimalLongitude": -80.5280, "coordinateUncertainty": 200.0, "lithic_outcrop_index": 0.78},
        ])
        res = spatially_balance_and_thin_vouchers(df_hab, spatial_thin_km=15.0, preserve_all_gold=False)
        kept = set(res["catalogNumber"].tolist())
        self.assertEqual(kept, {"SAND_1", "ROCK_1"})


if __name__ == "__main__":
    unittest.main()
