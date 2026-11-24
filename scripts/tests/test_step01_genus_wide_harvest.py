#!/usr/bin/env python3
"""Unit tests for Step 01 genus-wide Packera discovery, explicit taxa restriction, and cross-stage CLI handoffs."""

import tempfile
import unittest
from pathlib import Path
from unittest.mock import MagicMock, patch

from scripts.core.harvester import VoucherHarvester


class TestGenusWideHarvestAndCliHandoffs(unittest.TestCase):
    """Test suite for genus-wide Packera harvest discovery, explicit taxa filtering, and cross-stage CLI handoffs."""

    def test_genus_wide_packera_harvest_when_taxa_empty(self):
        """Running VoucherHarvester with taxa=[] dynamically discovers and harvests genus-wide Packera species."""
        import scripts.core.config as core_cfg
        import scripts.core.harvester as core_harv
        from scripts.core.harvester import discover_genus_packera_species_gbif

        self.assertEqual(core_cfg.DEFAULT_TARGET_TAXA, [])
        cfg = core_cfg.PipelineConfig.from_yaml()
        self.assertEqual(cfg.taxa.target_species, [])
        self.assertEqual(cfg.taxa.synonyms, {})
        self.assertEqual(cfg.taxa.outgroups, [])
        self.assertFalse(hasattr(core_harv, "EASTERN_NA_PACKERA_TAXA"))
        self.assertFalse(hasattr(core_harv, "EASTERN_NA_GBIF_QUERY_TAXA"))

        # 1. Dynamic genus-wide species discovery via GBIF Species Backbone + name_lookup/name_usage
        mock_species = MagicMock()
        mock_species.name_backbone.return_value = {
            "usageKey": 3109225,
            "canonicalName": "Packera",
            "rank": "GENUS",
            "matchType": "EXACT",
        }
        mock_species.name_lookup.return_value = {
            "results": [
                {"canonicalName": "Packera bolanderi", "rank": "SPECIES", "taxonomicStatus": "ACCEPTED"},
                {"canonicalName": "Packera cardamine", "rank": "SPECIES", "taxonomicStatus": "ACCEPTED"},
            ]
        }
        mock_species.name_usage.return_value = {
            "results": [
                {"canonicalName": "Packera franciscana", "rank": "SPECIES", "taxonomicStatus": "ACCEPTED"},
                {"canonicalName": "Senecio bolanderi", "rank": "SPECIES", "taxonomicStatus": "SYNONYM"},
                {"canonicalName": "Packera dubitabilis", "rank": "SPECIES", "taxonomicStatus": "DOUBTFUL"},
            ]
        }

        queried_names = []

        def _fake_occ_by_species(**kwargs):
            if "typeStatus" in kwargs:
                return {"results": [], "count": 0}
            sci = kwargs.get("scientificName")
            queried_names.append(sci)
            sci_ord = len(dict.fromkeys(queried_names))
            recs = [
                {
                    "key": sci_ord * 10000 + idx + 1,
                    "scientificName": f"{sci} Auth.",
                    "stateProvince": "North Carolina",
                    "decimalLatitude": 35.80 + idx * 0.05,
                    "decimalLongitude": -79.50 - idx * 0.05,
                    "coordinateUncertaintyInMeters": 100.0,
                    "year": 2021,
                    "month": 5,
                    "day": 14 + idx,
                    "recordedBy": "R. Collector",
                    "recordNumber": str(sci_ord * 1000 + idx),
                    "media": [{"identifier": f"https://example.com/{sci.replace(' ', '_')}_{idx}.jpg", "type": "StillImage", "format": "image/jpeg"}],
                    "catalogNumber": f"CAT_{sci.replace(' ', '_').upper()}_{idx}",
                    "institutionCode": "TEX",
                }
                for idx in range(3)
            ]
            return {
                "results": recs,
                "count": len(recs),
            }

        with patch("scripts.core.harvester.species", mock_species), \
             patch("scripts.core.harvester.occ.search", side_effect=_fake_occ_by_species):
            discovered = discover_genus_packera_species_gbif()
            self.assertEqual(discovered, ["Packera bolanderi", "Packera cardamine", "Packera franciscana"])

            harvester = VoucherHarvester(taxa=[], max_records_per_taxon=5, expand_synonyms=False)
            self.assertEqual(harvester.taxa, ["Packera bolanderi", "Packera cardamine", "Packera franciscana"])
            df_genus = harvester.harvest()

        self.assertEqual(list(dict.fromkeys(queried_names)), ["Packera bolanderi", "Packera cardamine", "Packera franciscana"])
        self.assertEqual(len(df_genus), 9)
        self.assertTrue(all(str(p).endswith((".webp", ".jpg")) for p in df_genus["image_path"]))
        for unrequested in ("Packera dubia", "Packera anonyma", "Packera paupercula", "Packera plattensis"):
            self.assertNotIn(unrequested, queried_names)

        # 2. Direct genus-level GBIF occurrence query (taxa=["Packera"]) filters indets and resolves Senecio synonyms
        mock_direct_species = MagicMock()

        def _backbone_dispatch(**kwargs):
            if kwargs.get("rank") == "genus":
                return {"usageKey": 3109225, "canonicalName": "Packera", "matchType": "EXACT"}
            if "cardamine" in str(kwargs.get("name", "")).lower():
                return {
                    "usageKey": 777111,
                    "canonicalName": "Senecio cardamine",
                    "species": "Packera cardamine",
                    "accepted": "Packera cardamine (Greene) W.A.Weber & Á.Löve",
                }
            return {}

        mock_direct_species.name_backbone.side_effect = _backbone_dispatch
        mock_direct_species.name_usage.return_value = {"results": []}

        direct_search_calls = []

        def _fake_direct_genus_occ(**kwargs):
            direct_search_calls.append(dict(kwargs))
            return {
                "results": [
                    {
                        "key": 9001,
                        "scientificName": "Packera hartiana (A.Heller) W.A.Weber & Á.Löve",
                        "stateProvince": "Texas",
                        "decimalLatitude": 31.1,
                        "decimalLongitude": -95.1,
                        "coordinateUncertaintyInMeters": 50.0,
                        "year": 2020,
                        "month": 4,
                        "day": 20,
                        "media": [{"identifier": "https://example.com/hartiana.jpg", "type": "StillImage", "format": "image/jpeg"}],
                        "catalogNumber": "TEX_HART_01",
                        "institutionCode": "TEX",
                    },
                    {
                        "key": 9002,
                        "scientificName": "Senecio cardamine Greene",
                        "stateProvince": "Texas",
                        "decimalLatitude": 31.3,
                        "decimalLongitude": -95.3,
                        "coordinateUncertaintyInMeters": 50.0,
                        "year": 2020,
                        "month": 4,
                        "day": 22,
                        "media": [{"identifier": "https://example.com/cardamine.jpg", "type": "StillImage", "format": "image/jpeg"}],
                        "catalogNumber": "TEX_CARD_02",
                        "institutionCode": "TEX",
                    },
                    {
                        "key": 9003,
                        "scientificName": "Packera sp.",
                        "stateProvince": "Texas",
                        "decimalLatitude": 31.4,
                        "decimalLongitude": -95.4,
                        "coordinateUncertaintyInMeters": 50.0,
                        "year": 2020,
                        "month": 4,
                        "day": 25,
                        "media": [{"identifier": "https://example.com/indet.jpg", "type": "StillImage", "format": "image/jpeg"}],
                        "catalogNumber": "TEX_INDET_03",
                        "institutionCode": "TEX",
                    },
                ],
                "count": 3,
            }

        core_harv._DYNAMIC_SYNONYM_TO_CANONICAL.pop("senecio cardamine", None)
        core_harv._GBIF_SYNONYM_CACHE.pop("senecio cardamine", None)
        with patch("scripts.core.harvester.species", mock_direct_species), \
             patch("scripts.core.harvester.occ.search", side_effect=_fake_direct_genus_occ):
            harvester_direct = VoucherHarvester(taxa=["Packera"], max_records_per_taxon=10)
            df_direct = harvester_direct.harvest()

        self.assertEqual(direct_search_calls[0].get("taxonKey"), 3109225)
        self.assertEqual(set(df_direct["catalogNumber"]), {"TEX_HART_01", "TEX_CARD_02"})
        card_row = df_direct[df_direct["catalogNumber"] == "TEX_CARD_02"].iloc[0]
        self.assertEqual(card_row["scientificName"], "Packera cardamine")
        self.assertEqual(card_row["species_raw"], "Senecio cardamine Greene")

        # 3. Fallback when taxa=[] and species discovery returns empty: queries genus-level Packera directly
        mock_direct_species.name_lookup.return_value = {"results": []}
        direct_search_calls.clear()
        with patch("scripts.core.harvester.species", mock_direct_species), \
             patch("scripts.core.harvester.occ.search", side_effect=_fake_direct_genus_occ):
            harvester_fallback = VoucherHarvester(taxa=[], max_records_per_taxon=10)
            self.assertEqual(harvester_fallback.taxa, ["Packera"])
            df_fb = harvester_fallback.harvest()
        self.assertEqual(direct_search_calls[0].get("taxonKey"), 3109225)
        self.assertEqual(set(df_fb["catalogNumber"]), {"TEX_HART_01", "TEX_CARD_02"})

        # 4. Dynamic genus-wide standardize_packera_taxon edge cases (hybrid markers & open nomenclature)
        self.assertEqual(core_harv.standardize_packera_taxon("Packera ×memmingeri (Britton) Weakley"), "Packera memmingeri")
        self.assertEqual(core_harv.standardize_packera_taxon("Packera x memmingeri"), "Packera memmingeri")
        self.assertIsNone(core_harv.standardize_packera_taxon("Packera x"))
        self.assertIsNone(core_harv.standardize_packera_taxon("Packera ×"))
        self.assertIsNone(core_harv.standardize_packera_taxon("Packera sp."))

    def test_explicit_taxa_harvests_strictly_requested_species(self):
        """Running VoucherHarvester with explicit taxa harvests strictly those species without injecting Packera dubia."""
        queried_taxa = []

        def _fake_occ(**kwargs):
            if "typeStatus" in kwargs:
                return {"results": [], "count": 0}
            sci = kwargs.get("scientificName")
            queried_taxa.append(sci)
            return {
                "results": [
                    {
                        "key": 7101 if "aurea" in str(sci) else 7102,
                        "scientificName": sci,
                        "stateProvince": "North Carolina",
                        "decimalLatitude": 35.91,
                        "decimalLongitude": -79.05,
                        "coordinateUncertaintyInMeters": 50.0,
                        "year": 2022,
                        "month": 4,
                        "day": 18,
                        "media": [{"identifier": f"https://example.com/{sci.replace(' ', '_')}.jpg", "type": "StillImage", "format": "image/jpeg"}],
                        "catalogNumber": f"NCU_{sci.split()[-1].upper()}_01",
                        "institutionCode": "NCU",
                    }
                ],
                "count": 1,
            }

        with patch("scripts.core.harvester.occ.search", side_effect=_fake_occ):
            harvester = VoucherHarvester(
                taxa=["Packera aurea", "Packera obovata"],
                max_records_per_taxon=5,
                expand_synonyms=False,
            )
            self.assertEqual(harvester.taxa, ["Packera aurea", "Packera obovata"])
            df = harvester.harvest()

        self.assertEqual(list(dict.fromkeys(queried_taxa)), ["Packera aurea", "Packera obovata"])
        self.assertEqual(set(df["scientificName"]), {"Packera aurea", "Packera obovata"})
        for forbidden in ("Packera dubia", "Packera anonyma", "Packera paupercula", "Packera plattensis"):
            self.assertNotIn(forbidden, queried_taxa)
            self.assertNotIn(forbidden, set(df["scientificName"]))

        # Verify main.run_harvest omits --taxa when none supplied and passes strictly explicit --taxa
        import argparse
        import main as root_main

        cfg_default = root_main.PipelineConfig.from_yaml()
        with tempfile.TemporaryDirectory() as tmp_d:
            tmp_csv = Path(tmp_d) / "curated_vouchers.csv"
            tmp_csv.write_text("catalogNumber,scientificName\nNCU1,Packera aurea\n", encoding="utf-8")
            base_ns = dict(
                run_id=None,
                out_dir=Path(tmp_d) / "raw",
                output_csv=tmp_csv,
                max_records=10,
                max_uncertainty=1000.0,
                min_megapixels=8.0,
                min_file_size_kb=500.0,
                concurrency=4,
                download_images=False,
                force=False,
                verbose=False,
            )

            with patch.object(root_main.subprocess, "run") as mock_sub_empty:
                root_main.run_harvest(argparse.Namespace(taxa=None, **base_ns), cfg_default)
                cmd_empty = mock_sub_empty.call_args[0][0]
                self.assertNotIn("--taxa", cmd_empty)
                self.assertNotIn("--taxon", cmd_empty)

            with patch.object(root_main.subprocess, "run") as mock_sub_explicit:
                root_main.run_harvest(
                    argparse.Namespace(taxa=["Packera aurea", "Packera obovata"], **base_ns),
                    cfg_default,
                )
                cmd_explicit = mock_sub_explicit.call_args[0][0]
                self.assertTrue("--taxa" in cmd_explicit or "--taxon" in cmd_explicit)
                taxa_flag = "--taxa" if "--taxa" in cmd_explicit else "--taxon"
                taxa_idx = cmd_explicit.index(taxa_flag)
                self.assertEqual(cmd_explicit[taxa_idx + 1 : taxa_idx + 3], ["Packera aurea", "Packera obovata"])
                self.assertNotIn("Packera dubia", cmd_explicit)

            main_parser = root_main.build_parser()
            parsed_harvest = main_parser.parse_args(["harvest", "--taxon", "Packera aurea", "Packera obovata"])
            self.assertEqual(parsed_harvest.taxa, ["Packera aurea", "Packera obovata"])
            parsed_all_empty = main_parser.parse_args(["run-all"])
            self.assertIsNone(parsed_all_empty.taxa)

    def test_cross_stage_cli_table_handoffs_in_main_and_run_pipeline(self):
        """Verify main.py and run_pipeline.sh wire all cross-stage CLI table handoffs across Steps 04, 06, 07, and 08."""
        import os
        import subprocess
        import main as root_main

        def _arg_value(cmd_list, flag_name):
            self.assertIn(flag_name, cmd_list, f"Missing CLI flag {flag_name!r} in command: {cmd_list}")
            idx = cmd_list.index(flag_name)
            self.assertLess(idx + 1, len(cmd_list), f"Flag {flag_name!r} has no value in command: {cmd_list}")
            return cmd_list[idx + 1]

        with tempfile.TemporaryDirectory() as tmp_d:
            tmp_run = Path(tmp_d) / "run_handoff_audit"
            tables_dir = tmp_run / "data" / "tables"
            figures_dir = tmp_run / "outputs" / "figures"
            contours_dir = tmp_run / "data" / "contours"
            tables_dir.mkdir(parents=True, exist_ok=True)
            figures_dir.mkdir(parents=True, exist_ok=True)
            contours_dir.mkdir(parents=True, exist_ok=True)
            (contours_dir / "NCU001_leaf1.csv").write_text("x,y\n0,0\n1,1\n", encoding="utf-8")

            cfg_scoped = root_main.PipelineConfig.from_yaml().with_run_dir(tmp_run)
            parser = root_main.build_parser()

            # 1. Verify Step 03 -> Step 04 handoff in main.run_morphometrics()
            morph_args = parser.parse_args([
                "morphometrics",
                "--focal-taxon", "Packera aurea",
                "--intraspecific-mode",
            ])
            with patch.dict(os.environ, {"PLANTMATICS_RUN_DIR": "", "PACKERA_RUN_DIR": ""}, clear=False), \
                 patch.object(root_main.shutil, "which", return_value="/usr/bin/Rscript"), \
                 patch.object(root_main, "verify_file_exists"), \
                 patch.object(root_main, "verify_dir_has_files"), \
                 patch.object(root_main.subprocess, "run") as mock_morph_sub:
                root_main.run_morphometrics(morph_args, cfg_scoped)
                self.assertEqual(mock_morph_sub.call_count, 2)
                cmd_step04 = mock_morph_sub.call_args_list[1][0][0]

            self.assertTrue(any("04_gmm_morphotools.R" in tok for tok in cmd_step04))
            self.assertEqual(_arg_value(cmd_step04, "--reproductive"), str(tables_dir / "voucher_reproductive_metrics.csv"))
            self.assertEqual(_arg_value(cmd_step04, "--allometry-report"), str(tables_dir / "allometry_permanova_summary.csv"))
            self.assertEqual(_arg_value(cmd_step04, "--allometry-pairwise"), str(tables_dir / "allometry_pairwise_slopes.csv"))
            self.assertEqual(_arg_value(cmd_step04, "--allometry-plot"), str(figures_dir / "fig04a_foliar_allometry_audit.png"))
            self.assertEqual(_arg_value(cmd_step04, "--output-batch-audit"), str(tables_dir / "batch_variance_audit.csv"))
            self.assertEqual(_arg_value(cmd_step04, "--focal-taxon"), "Packera aurea")
            self.assertIn("--intraspecific-mode", cmd_step04)

            # 2. Verify Step 06 and Step 07 handoffs in main.run_synthesis()
            synth_args = parser.parse_args([
                "synthesis",
                "--focal-taxon", "Packera aurea",
                "--intraspecific-mode",
            ])
            with patch.dict(os.environ, {"PLANTMATICS_RUN_DIR": "", "PACKERA_RUN_DIR": ""}, clear=False), \
                 patch.object(root_main.shutil, "which", return_value="/usr/bin/Rscript"), \
                 patch.object(root_main, "verify_file_exists"), \
                 patch.object(root_main.subprocess, "run") as mock_synth_sub:
                root_main.run_synthesis(synth_args, cfg_scoped)
                self.assertEqual(mock_synth_sub.call_count, 2)
                cmd_step06 = mock_synth_sub.call_args_list[0][0][0]
                cmd_step07 = mock_synth_sub.call_args_list[1][0][0]

            # Step 05 assertions
            self.assertTrue(any("05_multimodal_spatial_rf.R" in tok for tok in cmd_step06))
            reports_dir = tmp_run / "outputs" / "reports"
            self.assertEqual(_arg_value(cmd_step06, "--pheno-states"), str(tables_dir / "voucher_phenological_states.csv"))
            self.assertEqual(_arg_value(cmd_step06, "--ssurgo-table"), str(tables_dir / "ssurgo_edaphic_validation.csv"))
            self.assertEqual(_arg_value(cmd_step06, "--contingency-table"), str(reports_dir / "micro_edaphic_outcrop_contingency.csv"))
            self.assertEqual(_arg_value(cmd_step06, "--pheno-anomalies"), str(tables_dir / "phenological_anomalies.csv"))
            self.assertEqual(_arg_value(cmd_step06, "--pheno-summary"), str(reports_dir / "phenological_anomaly_summary.csv"))
            self.assertEqual(_arg_value(cmd_step06, "--pheno-plot"), str(figures_dir / "fig05b_hopkins_circular_phenology.png"))
            self.assertEqual(_arg_value(cmd_step06, "--focal-taxon"), "Packera aurea")
            self.assertIn("--intraspecific-mode", cmd_step06)

            # Step 07 synthesis + complex delimitation assertions
            self.assertTrue(any("07_triage_dashboard_synthesis.R" in tok for tok in cmd_step07))
            self.assertEqual(_arg_value(cmd_step07, "--harmonics"), str(tables_dir / "leaf_efa_harmonics.csv"))
            self.assertEqual(_arg_value(cmd_step07, "--output-queue"), str(tables_dir / "triage_queue.csv"))
            self.assertEqual(_arg_value(cmd_step07, "--output-html"), str(reports_dir / "Triage_Dashboard_Interactive.html"))
            self.assertEqual(_arg_value(cmd_step07, "--output-complexes"), str(reports_dir / "empirical_taxon_complexes.csv"))
            self.assertEqual(_arg_value(cmd_step07, "--output-affinity"), str(reports_dir / "pairwise_taxon_affinity_matrix.csv"))
            self.assertEqual(_arg_value(cmd_step07, "--output-neighbors"), str(reports_dir / "focal_auto_neighbors.txt"))
            self.assertEqual(_arg_value(cmd_step07, "--output-complex-plot"), str(figures_dir / "fig07b_taxon_complex_delimitation_network.png"))
            self.assertEqual(_arg_value(cmd_step07, "--focal-taxon"), "Packera aurea")
            self.assertIn("--intraspecific-mode", cmd_step07)

        # 3. Verify run_pipeline.sh unified 7-stage ordering (Step 05 -> Step 06 genomics -> Step 07 synthesis)
        run_sh = Path(__file__).resolve().parents[2] / "run_pipeline.sh"
        sh_text = run_sh.read_text(encoding="utf-8")

        pos_step05 = sh_text.index("scripts/analysis/05_multimodal_spatial_rf.R")
        pos_step06 = sh_text.index("scripts/data_prep/06_ingest_genetic_data.py")
        pos_step07 = sh_text.index("scripts/analysis/07_triage_dashboard_synthesis.R")
        self.assertLess(pos_step05, pos_step06, "Step 06 genomics must follow Step 05 in run_pipeline.sh")
        self.assertLess(pos_step06, pos_step07, "Step 07 synthesis must run AFTER Step 06 in run_pipeline.sh")

        step07_block = sh_text[pos_step07:]
        self.assertIn('--multimodal-flags "${RUN_DIR}/data/tables/multimodal_conflict_flags.csv"', step07_block)
        self.assertIn('--output-complexes "${RUN_DIR}/outputs/reports/empirical_taxon_complexes.csv"', step07_block)
        self.assertIn('--output-complex-plot "${RUN_DIR}/outputs/figures/fig07b_taxon_complex_delimitation_network.png"', step07_block)

        proc = subprocess.run(["bash", "-n", str(run_sh)], capture_output=True, text=True, check=False)
        self.assertEqual(proc.returncode, 0, f"bash -n run_pipeline.sh failed: {proc.stderr}")


if __name__ == "__main__":
    unittest.main()
