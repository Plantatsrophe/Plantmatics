#!/usr/bin/env python3
"""Unit tests for Step 01 dynamic synonym unification and Cis-Mississippi disjunct cohort pruning."""

import tempfile
import unittest
from pathlib import Path
from unittest.mock import MagicMock, patch

import pandas as pd

from scripts.tests._r_session import rscript_available, run_warm_r_cmd


class TestStep1SynonymUnificationAndCisDisjunct(unittest.TestCase):
    """Tests for dynamic synonym unification, R authority YAML config loading, and Cis-Mississippi cohort pruning."""

    def test_zero_hardcoded_synonyms_and_dynamic_packera_tomentosa_dubia_unification(self):
        """Verify zero hardcoded taxon names in source and 100% dynamic unification of Packera tomentosa into Packera dubia."""
        import subprocess
        import main as root_main
        from scripts.core import harvester as ch
        from scripts.core.harvester import VoucherHarvester

        project_root = Path(__file__).resolve().parents[2]
        harv_src = (project_root / "scripts/core/harvester.py").read_text(encoding="utf-8")
        usda_src = (project_root / "scripts/core/usda_plants_ranges.py").read_text(encoding="utf-8")
        morph_dir = project_root / "scripts" / "morphometrics"
        aw_r_paths = [morph_dir / "authority_weighting_utils.R", *sorted(morph_dir.glob("auth_*.R"))]
        aw_r_src = "\n".join(p.read_text(encoding="utf-8") for p in aw_r_paths if p.is_file())
        self.assertNotIn('"packera tomentosa"', harv_src.lower())
        self.assertNotIn('"senecio tomentosus"', harv_src.lower())
        self.assertNotIn("fallback_map", usda_src)
        self.assertNotIn('"packera dubia":', usda_src.lower())
        self.assertNotIn('grepl("dubia|tomentosa', aw_r_src)

        for k in ("packera dubia", "packera tomentosa", "senecio tomentosus", "cineraria dubia", "cineraria tomentosa", "senecio dubius"):
            ch._DYNAMIC_SYNONYM_TO_CANONICAL.pop(k, None)
            ch._GBIF_SYNONYM_CACHE.pop(k, None)
            ch._GBIF_BACKBONE_SYNONYM_FLAG.pop(k, None)

        mock_gbif_sp = MagicMock()
        def _gbif_bb(**kwargs):
            nm = str(kwargs.get("name") or kwargs.get("scientificName") or "").strip().lower()
            if "packera dubia" in nm:
                return {
                    "matchType": "EXACT",
                    "status": "ACCEPTED",
                    "usageKey": 5387757,
                    "canonicalName": "Packera dubia",
                    "scientificName": "Packera dubia (Spreng.) Trock & Mabb.",
                    "species": "Packera dubia",
                }
            if "packera tomentosa" in nm:
                return {
                    "matchType": "EXACT",
                    "status": "SYNONYM",
                    "usageKey": 3109265,
                    "acceptedUsageKey": 5387757,
                    "canonicalName": "Packera tomentosa",
                    "scientificName": "Packera tomentosa (Michx.) C. Jeffrey",
                    "accepted": "Packera dubia (Spreng.) Trock & Mabb.",
                    "species": "Packera dubia",
                }
            return {"matchType": "NONE"}

        def _gbif_nu(**kwargs):
            if kwargs.get("data") == "synonyms":
                return {
                    "results": [
                        {"canonicalName": "Packera tomentosa", "scientificName": "Packera tomentosa (Michx.) C. Jeffrey", "taxonomicStatus": "SYNONYM"},
                        {"canonicalName": "Senecio tomentosus", "scientificName": "Senecio tomentosus Michx.", "taxonomicStatus": "SYNONYM"},
                        {"canonicalName": "Cineraria dubia", "scientificName": "Cineraria dubia Spreng.", "taxonomicStatus": "HOMOTYPIC_SYNONYM"},
                        {"canonicalName": "Senecio dubius", "scientificName": "Senecio dubius (Spreng.)", "taxonomicStatus": "HOMOTYPIC_SYNONYM"},
                        {"canonicalName": "Packera aurea", "scientificName": "Packera aurea auct. non (L.)", "taxonomicStatus": "MISAPPLIED"},
                    ]
                }
            return {}

        mock_gbif_sp.name_backbone.side_effect = _gbif_bb
        mock_gbif_sp.name_usage.side_effect = _gbif_nu

        def _gbif_occ(**kwargs):
            if "typeStatus" in kwargs or kwargs.get("country") == "CA":
                return {"results": [], "count": 0}
            sci_q = str(kwargs.get("scientificName") or "")
            recs = []
            if sci_q == "Packera dubia":
                recs.append({
                    "key": 90101,
                    "scientificName": "Packera dubia (Spreng.) Trock & Mabb.",
                    "taxonomicStatus": "ACCEPTED",
                    "species": "Packera dubia",
                    "stateProvince": "North Carolina",
                    "decimalLatitude": 35.91,
                    "decimalLongitude": -79.05,
                    "year": 2021, "month": 4, "day": 18,
                    "recordedBy": "B. Fuller", "recordNumber": "401",
                    "identifiedBy": "A.S. Weakley",
                    "media": [{"identifier": "https://example.org/dubia1.jpg", "type": "StillImage"}],
                    "catalogNumber": "NCU_DUBIA_01",
                    "institutionCode": "NCU",
                })
            elif sci_q == "Packera tomentosa":
                recs.append({
                    "key": 90102,
                    "scientificName": "Packera tomentosa (Michx.) C. Jeffrey",
                    "taxonomicStatus": "SYNONYM",
                    "acceptedScientificName": "Packera dubia (Spreng.) Trock & Mabb.",
                    "species": "Packera dubia",
                    "stateProvince": "South Carolina",
                    "decimalLatitude": 34.00,
                    "decimalLongitude": -81.03,
                    "year": 1998, "month": 4, "day": 12,
                    "recordedBy": "J.B. Nelson", "recordNumber": "19200",
                    "identifiedBy": "T.M. Barkley",
                    "media": [{"identifier": "https://example.org/tomentosa1.jpg", "type": "StillImage"}],
                    "catalogNumber": "USCH_TOM_02",
                    "institutionCode": "USCH",
                })
            elif sci_q == "Senecio tomentosus":
                recs.append({
                    "key": 90103,
                    "scientificName": "Senecio tomentosus Michx.",
                    "taxonomicStatus": "SYNONYM",
                    "acceptedScientificName": "Packera dubia (Spreng.) Trock & Mabb.",
                    "species": "Packera dubia",
                    "stateProvince": "Virginia",
                    "decimalLatitude": 36.85,
                    "decimalLongitude": -76.29,
                    "year": 1968, "month": 5, "day": 3,
                    "recordedBy": "A.E. Radford", "recordNumber": "45012",
                    "identifiedBy": "A.E. Radford",
                    "media": [{"identifier": "https://example.org/sentom1.jpg", "type": "StillImage"}],
                    "catalogNumber": "NCU_SENTOM_03",
                    "institutionCode": "NCU",
                })
            return {"results": recs, "count": len(recs)}

        # Case A: Querying Packera dubia expands and combines Packera tomentosa + Senecio tomentosus into Packera dubia
        with patch("scripts.core.harvester.species", mock_gbif_sp), \
             patch("scripts.core.harvester.occ.search", side_effect=_gbif_occ):
            harv_dubia = VoucherHarvester(taxa=["Packera dubia"], expand_synonyms=True)
            self.assertIn("Packera tomentosa", harv_dubia.taxa)
            self.assertIn("Senecio tomentosus", harv_dubia.taxa)
            self.assertNotIn("Packera aurea", harv_dubia.taxa)  # MISAPPLIED rejected!
            df_combined = harv_dubia.harvest()

        self.assertEqual(set(df_combined["catalogNumber"]), {"NCU_DUBIA_01", "USCH_TOM_02", "NCU_SENTOM_03"})
        self.assertEqual(set(df_combined["species_standardized"]), {"Packera dubia"})
        tom_row = df_combined[df_combined["catalogNumber"] == "USCH_TOM_02"].iloc[0]
        self.assertEqual(tom_row["scientificName"], "Packera dubia")
        self.assertEqual(tom_row["species_raw"], "Packera tomentosa (Michx.) C. Jeffrey")

        # Case B: Querying unaccepted Packera tomentosa directly resolves canonical to Packera dubia and combines both
        for k in ("packera dubia", "packera tomentosa", "senecio tomentosus", "cineraria dubia", "senecio dubius"):
            ch._DYNAMIC_SYNONYM_TO_CANONICAL.pop(k, None)
            ch._GBIF_SYNONYM_CACHE.pop(k, None)
            ch._GBIF_BACKBONE_SYNONYM_FLAG.pop(k, None)

        with patch("scripts.core.harvester.species", mock_gbif_sp), \
             patch("scripts.core.harvester.occ.search", side_effect=_gbif_occ):
            harv_tom = VoucherHarvester(taxa=["Packera tomentosa"], expand_synonyms=True)
            df_from_tom = harv_tom.harvest()
            self.assertEqual(set(df_from_tom["catalogNumber"]), {"NCU_DUBIA_01", "USCH_TOM_02", "NCU_SENTOM_03"})
            self.assertEqual(set(df_from_tom["species_standardized"]), {"Packera dubia"})

        # Case C: Occurrence-level SYNONYM status combines Packera tomentosa into Packera dubia even when expand_synonyms=False
        for k in ("packera dubia", "packera tomentosa", "senecio tomentosus"):
            ch._DYNAMIC_SYNONYM_TO_CANONICAL.pop(k, None)
            ch._GBIF_SYNONYM_CACHE.pop(k, None)
            ch._GBIF_BACKBONE_SYNONYM_FLAG.pop(k, None)

        def _mixed_occ(**kwargs):
            if "typeStatus" in kwargs or kwargs.get("country") == "CA":
                return {"results": [], "count": 0}
            return {
                "results": [
                    {
                        "key": 90201,
                        "scientificName": "Packera dubia (Spreng.) Trock & Mabb.",
                        "taxonomicStatus": "ACCEPTED",
                        "species": "Packera dubia",
                        "stateProvince": "North Carolina",
                        "decimalLatitude": 35.91, "decimalLongitude": -79.05,
                        "year": 2021, "recordedBy": "B. Fuller", "recordNumber": "401",
                        "media": [{"identifier": "https://example.org/d.jpg", "type": "StillImage"}],
                        "catalogNumber": "NCU_D_01", "institutionCode": "NCU",
                    },
                    {
                        "key": 90202,
                        "scientificName": "Packera tomentosa (Michx.) C. Jeffrey",
                        "taxonomicStatus": "SYNONYM",
                        "acceptedScientificName": "Packera dubia (Spreng.) Trock & Mabb.",
                        "species": "Packera dubia",
                        "stateProvince": "South Carolina",
                        "decimalLatitude": 34.00, "decimalLongitude": -81.03,
                        "year": 1998, "recordedBy": "J.B. Nelson", "recordNumber": "19200",
                        "media": [{"identifier": "https://example.org/t.jpg", "type": "StillImage"}],
                        "catalogNumber": "USCH_T_02", "institutionCode": "USCH",
                    },
                ],
                "count": 2,
            }

        with patch("scripts.core.harvester.occ.search", side_effect=_mixed_occ):
            harv_occ_syn = VoucherHarvester(taxa=["Packera dubia"], expand_synonyms=False)
            df_occ_syn = harv_occ_syn.harvest()
            self.assertEqual(set(df_occ_syn["species_standardized"]), {"Packera dubia"})

        # Verify main.py forwards the 4 Step 1 flags when explicitly specified
        parser = root_main.build_parser()
        args_flags = parser.parse_args([
            "harvest",
            "--no-harvest-full-continental-range",
            "--spatial-thin-km", "12.5",
            "--max-per-state", "40",
            "--max-per-region", "80",
        ])
        cfg_stub = {
            "paths": {"raw_vouchers_dir": "/tmp/r", "curated_vouchers_csv": "/tmp/c.csv"},
            "harvesting": {
                "max_records_per_taxon": 100, "max_uncertainty_meters": 5000,
                "min_megapixels": 8.0, "min_file_size_kb": 500, "download_concurrency": 5,
            },
        }
        with patch("subprocess.run") as m_sub, patch("main.verify_file_exists"):
            root_main.run_harvest(args_flags, cfg_stub)
            sent_cmd = m_sub.call_args[0][0]
            self.assertIn("--no-harvest-full-continental-range", sent_cmd)
            self.assertNotIn("--spatial-thin-km", sent_cmd)
            self.assertNotIn("--max-per-state", sent_cmd)
            self.assertNotIn("--max-per-region", sent_cmd)

        # Verify R .load_authority_weight_config reads nested thresholds$authority_weights from YAML
        self.assertIn("cfg$thresholds$authority_weights", aw_r_src)
        self.assertIn("dyn_equiv", aw_r_src)
        if rscript_available():
            with tempfile.TemporaryDirectory() as tmp_r_d:
                custom_yaml = Path(tmp_r_d) / "cfg.yaml"
                custom_yaml.write_text(
                    "thresholds:\n  authority_weights:\n    tier_0_type: 3.5\n    silver_max_pool_share: 0.35\n",
                    encoding="utf-8",
                )
                r_out = run_warm_r_cmd(
                    [
                        "Rscript", "-e",
                        f"source('scripts/morphometrics/authority_weighting_utils.R'); "
                        f"rm(list=ls(envir=.AUTH_WEIGHT_CFG_ENV), envir=.AUTH_WEIGHT_CFG_ENV); "
                        f"on.exit(rm(list=ls(envir=.AUTH_WEIGHT_CFG_ENV), envir=.AUTH_WEIGHT_CFG_ENV), add=TRUE); "
                        f"w <- .load_authority_weight_config({repr(str(custom_yaml))}); "
                        f"stopifnot(abs(w$tier_0_type - 3.5) < 1e-6, abs(w$silver_max_pool_share - 0.35) < 1e-6); "
                        f"cat('R_CFG_OK\\n')",
                    ],
                    stage_key="step03",
                    cwd=project_root,
                )
                self.assertEqual(r_out.returncode, 0, r_out.stderr)
                self.assertIn("R_CFG_OK", r_out.stdout)

    def test_cis_mississippi_disjunct_filter_and_prune_eastern_cohort(self):
        """Verify Cis-Mississippi filter keeps P. dubia / P. obovata disjuncts up to -100W, drops P. cana & western P. paupercula, and prunes SQLite + CSVs."""
        import sqlite3
        from scripts.core.db import get_db_connection
        from scripts.data_prep.prune_eastern_cohort import prune_workspace_cohort

        with tempfile.TemporaryDirectory() as tmp_d:
            ws = Path(tmp_d)
            tables_dir = ws / "data" / "tables"
            tables_dir.mkdir(parents=True, exist_ok=True)
            db_p = ws / "data" / "plantmatics.sqlite"

            conn = get_db_connection(db_path=db_p, init_schema=True, workspace_root=ws)
            rows = [
                # 1. Packera dubia cis-Mississippi (NC) -> KEEP
                ("NCU_DUBIA_NC", "NCU", "Packera dubia", "Cineraria dubia", "North Carolina", 35.9, -79.0, "Tier_1_Gold", ""),
                # 2. Packera dubia trans-Mississippi disjunct in East Texas (-95.2W >= -100W) -> KEEP
                ("TEX_DUBIA_ETX", "TEX", "Packera dubia", "Senecio tomentosus", "Texas", 31.2, -95.2, "Tier_2_Silver", ""),
                # 3. Packera obovata Edwards Plateau TX disjunct (-98.6W >= -100W) -> KEEP
                ("BRIT_OBOVATA_TX", "BRIT", "Packera obovata", "Packera obovata", "Texas", 30.1, -98.6, "Tier_1_Gold", ""),
                # 4. Packera obovata West Texas stray (-103.9W < -100W) -> DROP
                ("SRSC_OBOVATA_WTX", "SRSC", "Packera obovata", "Packera obovata", "Texas", 30.6, -103.9, "Tier_3_Bronze", ""),
                # 5. Packera obovata Western Type Rescue (-104.5W Holotype) -> KEEP (nomenclatural type)
                ("NY_OBOVATA_TYPE", "NY", "Packera obovata", "Packera obovata", "New Mexico", 35.1, -104.5, "Tier_1_Gold", "Holotype"),
                # 6. Packera paupercula cis-Mississippi (NC) & Driftless IA (-91.5W) -> KEEP
                ("NCU_PAUP_NC", "NCU", "Packera paupercula", "Packera paupercula", "North Carolina", 36.1, -81.2, "Tier_1_Gold", ""),
                ("ISC_PAUP_IA", "ISC", "Packera paupercula", "Packera paupercula", "Iowa", 42.9, -91.5, "Tier_2_Silver", ""),
                # 7. Packera paupercula Western / Great Plains (SD -103.5W, TX -98.0W, BC -120.5W) -> DROP
                ("RM_PAUP_SD", "RM", "Packera paupercula", "Packera paupercula", "South Dakota", 43.8, -103.5, "Tier_2_Silver", ""),
                ("TEX_PAUP_TX", "TEX", "Packera paupercula", "Packera paupercula", "Texas", 30.4, -98.0, "Tier_3_Bronze", ""),
                ("UBC_PAUP_BC", "UBC", "Packera paupercula", "Packera paupercula", "British Columbia", 50.2, -120.5, "Tier_2_Silver", ""),
                # 8. Packera cana (Cordilleran/Plains with 1 stray IL sheet + 5 western sheets) -> DROP ENTIRE SPECIES
                ("ILL_CANA_STRAY", "ILL", "Packera cana", "Packera cana", "Illinois", 40.1, -89.2, "Tier_3_Bronze", ""),
                ("COLO_CANA_1", "COLO", "Packera cana", "Packera cana", "Colorado", 40.0, -105.3, "Tier_1_Gold", ""),
                ("COLO_CANA_2", "COLO", "Packera cana", "Packera cana", "Wyoming", 43.1, -108.2, "Tier_1_Gold", ""),
                ("MONT_CANA_3", "MONT", "Packera cana", "Packera cana", "Montana", 46.5, -112.0, "Tier_2_Silver", ""),
            ]
            with conn:
                for cat, inst, sn, sr, st, lat, lon, tier, tstat in rows:
                    conn.execute(
                        """
                        INSERT INTO vouchers (
                            catalogNumber, institutionCode, scientificName, species_raw, species_standardized,
                            stateProvince, decimalLatitude, decimalLongitude, determiner_tier, type_status
                        ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?);
                        """,
                        (cat, inst, sn, sr, sn, st, lat, lon, tier, tstat),
                    )
                    conn.execute(
                        "INSERT INTO plant_individuals (catalogNumber, plant_individual_id) VALUES (?, 0);",
                        (cat,),
                    )
                    conn.execute(
                        "INSERT INTO extracted_leaves (leaf_id, catalogNumber, plant_individual_id, leaf_zone, contour_xy_json) VALUES (?, ?, 0, 'basal', '[]');",
                        (f"{cat}_L1", cat),
                    )
            v_df = pd.read_sql_query("SELECT * FROM vouchers;", conn)
            v_df.to_csv(tables_dir / "curated_vouchers.csv", index=False)
            conn.close()

            res = prune_workspace_cohort(workspace_root=ws, db_path=db_p, vacuum=False, dry_run=False)
            expected_kept = {
                "NCU_DUBIA_NC",
                "TEX_DUBIA_ETX",
                "BRIT_OBOVATA_TX",
                "NY_OBOVATA_TYPE",
                "NCU_PAUP_NC",
                "ISC_PAUP_IA",
            }
            self.assertEqual(res["keep_catalog_count"], len(expected_kept))
            self.assertIn("Packera cana", res["dropped_species"])

            conn_chk = sqlite3.connect(str(db_p))
            try:
                rem_cats = {r[0] for r in conn_chk.execute("SELECT catalogNumber FROM vouchers;").fetchall()}
                rem_leaves = {r[0] for r in conn_chk.execute("SELECT catalogNumber FROM extracted_leaves;").fetchall()}
                self.assertEqual(rem_cats, expected_kept)
                self.assertEqual(rem_leaves, expected_kept)
            finally:
                conn_chk.close()

            pruned_csv = pd.read_csv(tables_dir / "curated_vouchers.csv")
            self.assertEqual(set(pruned_csv["catalogNumber"]), expected_kept)
            # Also verify missing extracted_leaf_manifest.csv was reconstructed automatically from SQLite!
            self.assertTrue((tables_dir / "extracted_leaf_manifest.csv").is_file())
            rebuilt_leaves = pd.read_csv(tables_dir / "extracted_leaf_manifest.csv")
            self.assertEqual(set(rebuilt_leaves["catalogNumber"]), expected_kept)


if __name__ == "__main__":
    unittest.main()
