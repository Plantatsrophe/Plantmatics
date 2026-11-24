#!/usr/bin/env python3
"""
Unit and Warm-R integration tests for public Packera Dryad/NCBI/SRA/BOLD downloads,
orthography and digit-prefix normalization, PhyTop node-tree lookup, focal tree
non-shadowing, and outgroup-exempt Newick MAD rogue-tip screening.
"""

from __future__ import annotations

import csv
import hashlib
import importlib
import io
import json
import os
import sqlite3
import sys
import tempfile
import unittest
import zipfile
from pathlib import Path
from unittest import mock

import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from scripts.core import config as cfg
from scripts.core.db import get_db_connection
from scripts.core.step01_biogeography_utils import _canonical_taxon_equiv_key, resolve_eastern_resident_whitelist
from scripts.data_prep import skim_offtarget_organellar as skim_mod
from scripts.data_prep import step06_public_cohort_ingestor as cohort_mod
from scripts.tests._r_session import rscript_available, run_warm_r_cmd

INGEST_MOD = importlib.import_module("scripts.data_prep.06_ingest_genetic_data")
PUBLIC_DL_MOD = importlib.import_module("scripts.data_prep.download_public_packera_sequences")


class TestStep06PublicBackboneDryad(unittest.TestCase):
    """Tests for Dryad/NCBI/SRA flat staging, orthography, PhyTop lookup, backbone polarization, and MAD QC."""

    def test_download_public_packera_sequences_flat_extraction_and_sha256(self) -> None:
        """Verify mocked Dryad (SHA-256 + nested zip), NCBI nuccore/SRA, and BOLD download flatly into subdirs."""
        class _MockHTTPResponse(io.BytesIO):
            def __enter__(self) -> "_MockHTTPResponse":
                return self

            def __exit__(self, *args: object) -> bool:
                self.close()
                return False

        def _zip_bytes(entries: dict[str, str]) -> bytes:
            buf = io.BytesIO()
            with zipfile.ZipFile(buf, "w", compression=zipfile.ZIP_DEFLATED) as zf:
                for arcname, content in entries.items():
                    zf.writestr(arcname, content.encode("utf-8"))
            return buf.getvalue()

        tree_zip = _zip_bytes({
            "nested/sub/trees/MO_Astral_lpp.tre": "((Packera_dubia:0.05,Packera_schweinitziana:0.05):0.10,Packera_aurea:0.20);\n",
            "nested/sub/trees/ASTRALPro_Astral_lpp.tre": "((Packera_dubia:0.06,Packera_schweinitziana:0.06):0.10,Packera_aurea:0.20);\n",
        })
        phytop_zip = _zip_bytes({
            "nested/phytop/MO_PhyTop.tsv": "taxon\tils\tih\nPackera dubia\t0.18\t0.06\nPackera schweinitziana\t0.14\t0.04\n",
        })
        align_zip = _zip_bytes({
            "nested/plastomes/18complete_Packera_plastomes.fasta": ">MK001|Packera_dubia\n" + ("ATGCTAGCTGATCGATCGTA" * 7) + "\n",
        })
        zip_payloads = {"101": ("speciestree_files.zip", tree_zip), "102": ("PhyTop_results.zip", phytop_zip), "103": ("alignment_files.zip", align_zip)}

        versions_json = json.dumps({"_embedded": {"stash:versions": [{"versionNumber": 1, "_links": {"self": {"href": "/api/v2/versions/999"}}}]}}).encode("utf-8")
        files_json = json.dumps({"_embedded": {"stash:files": [
            {"id": int(fid), "path": fname, "size": len(zbytes), "digestType": "sha-256", "digest": hashlib.sha256(zbytes).hexdigest(), "_links": {"stash:download": {"href": f"/api/v2/files/{fid}/download"}}}
            for fid, (fname, zbytes) in zip_payloads.items()
        ]}}).encode("utf-8")
        ncbi_nuccore_fasta = (
            ">NC_0001 Senecio tomentosus chloroplast, complete genome\n" + ("ACGTACGTACGTATCG" * 8) + "\n"
            ">MT_0002 Packera aurea internal transcribed spacer 1\n" + ("TGCATGCATGCATAGC" * 8) + "\n"
            ">JX_0003 Senecio tomentosus granule-bound starch synthase (waxy) gene\n" + ("GGCCAATTAAGGCCTT" * 8) + "\n"
        ).encode("utf-8")
        ena_tsv = (
            "run_accession\tsample_accession\tscientific_name\tlibrary_name\tfastq_ftp\tfastq_bytes\n"
            "SRR22900001\tSAMN01\tPackera dubia\tLIB1\tftp.sra.ebi.ac.uk/vol1/fastq/SRR229/001.fastq.gz\t1024000\n"
        ).encode("utf-8")
        sra_runinfo_csv = (
            "Run,BioProject,Sample,ScientificName,LibraryName,download_path,size_MB\n"
            "SRR22900002,PRJNA978568,SAMN02,Packera schweinitziana,LIB2,sra-pub-run-odp.s3.amazonaws.com/SRR22900002,12\n"
        ).encode("utf-8")
        bold_fasta = ">BOLD_01|Packera schweinitziana|ITS2\n" + ("ACGTGCATGCAT" * 11) + "\n"

        def _mock_open_url(req_or_url: object, data: object = None, timeout: object = 60) -> _MockHTTPResponse:
            _ = (data, timeout)
            url = req_or_url.full_url if hasattr(req_or_url, "full_url") else str(req_or_url)
            if "/datasets/" in url and url.endswith("/versions"):
                return _MockHTTPResponse(versions_json)
            if "/versions/999/files" in url:
                return _MockHTTPResponse(files_json)
            for fid, (_, zbytes) in zip_payloads.items():
                if f"/api/v2/files/{fid}/download" in url:
                    return _MockHTTPResponse(zbytes)
            if "esearch.fcgi" in url and "db=nuccore" in url:
                return _MockHTTPResponse(json.dumps({"esearchresult": {"count": "3", "idlist": ["101", "102", "103"]}}).encode("utf-8"))
            if "efetch.fcgi" in url and "db=nuccore" in url:
                return _MockHTTPResponse(ncbi_nuccore_fasta)
            if "ebi.ac.uk/ena/portal/api/filereport" in url:
                return _MockHTTPResponse(ena_tsv if "PRJNA907383" in url else b"")
            if "esearch.fcgi" in url and "db=sra" in url:
                return _MockHTTPResponse(json.dumps({"esearchresult": {"idlist": ["502"] if "PRJNA978568" in url else []}}).encode("utf-8"))
            if "efetch.fcgi" in url and "db=sra" in url:
                return _MockHTTPResponse(sra_runinfo_csv)
            if "portal.boldsystems.org/api/query" in url:
                return _MockHTTPResponse(bold_fasta.encode("utf-8"))
            raise AssertionError(f"Unexpected unmocked network URL in offline test: {url}")

        def _mock_opener_open(self_opener: object, fullurl: object, data: object = None, timeout: object = 60) -> _MockHTTPResponse:
            _ = self_opener
            return _mock_open_url(fullurl, data=data, timeout=timeout)

        with tempfile.TemporaryDirectory() as td:
            geno_dir = Path(td) / "genomics"
            with (
                mock.patch.dict(os.environ, {"DRYAD_API_TOKEN": "test_dryad_token", "NCBI_API_KEY": "test_ncbi_key"}),
                mock.patch("urllib.request.urlopen", side_effect=_mock_open_url),
                mock.patch("urllib.request.OpenerDirector.open", new=_mock_opener_open),
            ):
                summary = PUBLIC_DL_MOD.download_public_packera_sequences(genomics_dir=geno_dir)

            self.assertEqual(summary["dryad_status"], "downloaded")
            self.assertEqual(
                {p.name for p in geno_dir.iterdir() if p.is_dir()},
                {"public_nuclear_backbone", "public_phytop_conflict", "public_organellar_and_loci", "public_catalogs"},
            )
            self.assertEqual([p for p in geno_dir.iterdir() if p.is_file()], [])
            expected_rel_files = {
                "public_nuclear_backbone/public_packera_MO_Astral_lpp.tre",
                "public_nuclear_backbone/public_packera_ASTRALPro_Astral_lpp.tre",
                "public_phytop_conflict/public_packera_PhyTop_MO.tsv",
                "public_organellar_and_loci/public_packera_18complete_plastomes.fasta",
                "public_organellar_and_loci/public_packera_plastomes_and_cpdna.fasta",
                "public_organellar_and_loci/public_packera_its_ets.fasta",
                "public_organellar_and_loci/public_packera_waxy_nuclear.fasta",
                "public_organellar_and_loci/public_packera_bold_barcodes.fasta",
                "public_catalogs/public_packera_sra_runs.csv",
                "public_catalogs/public_packera_dryad_manifest.json",
                "public_catalogs/public_packera_sources_manifest.json",
            }
            actual_rel_files = {p.relative_to(geno_dir).as_posix() for p in geno_dir.rglob("*") if p.is_file()}
            self.assertTrue(expected_rel_files.issubset(actual_rel_files), f"Missing downloaded files: {expected_rel_files - actual_rel_files}")
            self.assertIn("Packera_dubia", (geno_dir / "public_organellar_and_loci" / "public_packera_plastomes_and_cpdna.fasta").read_text(encoding="utf-8"))
            self.assertIn("Packera_dubia", (geno_dir / "public_organellar_and_loci" / "public_packera_waxy_nuclear.fasta").read_text(encoding="utf-8"))
            sra_txt = (geno_dir / "public_catalogs" / "public_packera_sra_runs.csv").read_text(encoding="utf-8")
            self.assertIn("SRR22900001", sra_txt)
            self.assertIn("SRR22900002", sra_txt)
            self.assertIn("public_nuclear_backbone/public_packera_MO_Astral_lpp.tre", summary["files"])

    def test_public_reference_files_only_stay_noop_no_payload(self) -> None:
        with tempfile.TemporaryDirectory() as td:
            tmp = Path(td)
            geno_dir, out_dir = tmp / "genomics", tmp / "parsed_out"
            pub_rel_paths = (
                "public_nuclear_backbone/public_packera_MO_Astral_lpp.tre",
                "public_nuclear_backbone/public_packera_MO_concat_best.fasta",
                "public_organellar_and_loci/public_packera_18complete_plastomes.fasta",
                "public_catalogs/public_packera_sra_runs.csv",
                "public_catalogs/public_packera_dryad_manifest.json",
                "public_organellar_and_loci/public_packera_bold_barcodes.fasta",
            )
            for rel_p in pub_rel_paths:
                dst = geno_dir / rel_p
                dst.parent.mkdir(parents=True, exist_ok=True)
                dst.write_text("placeholder\n", encoding="utf-8")

            assets = INGEST_MOD.discover_genomic_inputs(geno_dir)
            self.assertFalse(assets.has_genetic_payload)
            self.assertEqual({p.name for p in assets.public_reference_files}, {Path(rp).name for rp in pub_rel_paths})
            self.assertEqual(INGEST_MOD.run_genomic_ingestion(genomics_dir=geno_dir, out_dir=out_dir)["status"], "skipped_empty")
            self.assertEqual(len([p for p in geno_dir.rglob("*") if p.is_file()]), len(pub_rel_paths))

    def test_genomic_naming_audit_parsers_orthography_and_outgroup_mad_qc(self) -> None:
        self.assertIn("Packera millefolium", cfg.CANONICAL_CIS_MISSISSIPPI_PACKERA_SPECIES)
        self.assertIn("Packera memmingeri", cfg.CANONICAL_CIS_MISSISSIPPI_PACKERA_SPECIES)
        self.assertEqual(_canonical_taxon_equiv_key("Packera millefolia"), "packera millefolium")
        self.assertIn("Packera millefolium", resolve_eastern_resident_whitelist())

        og = INGEST_MOD.parse_sample_code("EREHIENCS1M")
        self.assertEqual((og["taxon_code"], og["species_standardized"], og["population_id"], og["replicate_pair_id"], og["replicate_suffix"]), ("EREHIE", "Erechtites hieraciifolius", "NCS1", "EREHIENCS1", "M"))
        self.assertEqual(INGEST_MOD.parse_sample_code("EREHIE1-NCS1P")["species_standardized"], "Erechtites hieraciifolius")
        self.assertEqual(INGEST_MOD._normalize_event_date_part("3rd August 2020"), "2020-08-03")
        self.assertEqual(INGEST_MOD._normalize_event_date_part("August 3rd, 2020"), "2020-08-03")
        keys = INGEST_MOD._extract_meta_exsiccatae_keys({"collector": "Brandon Fuller", "Collector_ID": "21", "eventDate": "3rd August 2020"})
        self.assertTrue(any("21" in k and "2020-08-03" in k for k in keys))
        self.assertEqual(INGEST_MOD._canonical_taxon_equiv_key("Packera millefolia"), INGEST_MOD._canonical_taxon_equiv_key("Packera millefolium"))
        self.assertEqual(INGEST_MOD._canonical_taxon_equiv_key("Packera memmingeri"), INGEST_MOD._canonical_taxon_equiv_key("Packera \u00d7memmingeri"))

        # Outgroup-exempt Newick MAD screening: Erechtites hieraciifolius never pruned; multi-accession conspecific tips preserved
        nwk_qc = (
            "((Packera_dubia_acc1:0.020,Packera_dubia_acc2:0.028):0.02,"
            "((Packera_aurea:0.022,Packera_obovata:0.024):0.02,"
            "(Packera_anonyma:0.025,Packera_paupercula:0.026):0.02):0.02,"
            "Packera_rogue:0.950,Erechtites_hieraciifolius:0.950);\n"
        )
        tl, _sp, _cp, tip_hdrs = cohort_mod._parse_newick_tips_and_dists(nwk_qc)
        dub_hdrs = tip_hdrs.get("Packera dubia", [])
        self.assertEqual(set(dub_hdrs), {"Packera_dubia_acc1", "Packera_dubia_acc2"})
        kept = cohort_mod._mad_screen_tips(tl)
        kept_norm = {str(k).replace("_", " ").lower() for k in kept}
        self.assertTrue(any("erechtites" in k for k in kept_norm), f"Outgroup Erechtites must never be pruned: {kept}")
        self.assertFalse(any("rogue" in k for k in kept_norm), f"Rogue ingroup tip must be pruned: {kept}")
        self.assertTrue(any("dubia" in k for k in kept_norm))
        if any("acc1" in str(k).lower() for k in tl):
            self.assertIn("Packera_dubia_acc1", kept)
            self.assertIn("Packera_dubia_acc2", kept)

        with tempfile.TemporaryDirectory(prefix="packera_audit_nom_") as tmp:
            gdir = Path(tmp) / "genomics"
            gdir.mkdir(parents=True, exist_ok=True)
            fcsv = Path(tmp) / "Packera_Field_Specimen_Database.csv"
            fcsv.write_text("Collector_ID,Collector,EventDate,Taxon\n21,Brandon Fuller,3rd August 2020,Packera millefolia\n", encoding="utf-8")
            parsed_field = INGEST_MOD._parse_field_metadata([fcsv])
            self.assertIn("2020-08-03", json.dumps(parsed_field, default=str))
            self.assertIn("21", json.dumps(parsed_field, default=str))

            (gdir / "vcfref_packera_rm.txt.fasta").write_text(
                "".join(f">{q}\nATGCATGCATGCATGCATGC\n" for q in ("75949_5", "82071_300", "49275_9", "16157_6", "52797_44", "68356_89", "49053_5")),
                encoding="utf-8",
            )
            (gdir / "vcfref_packera_blast.txt").write_text(
                "Query= 75949_5\n>ref|NC_1| Packera chloroplast complete genome\n"
                "Query= 52797_44\n>ref|XM_1| 26S proteasome non-ATPase subunit\n"
                "Query= 68356_89\n>ref|XM_2| 40S ribosomal protein S3\n"
                "Query= 49053_5\n>ref|XM_3| probable plastid-lipid protein, chloroplastic-like\n"
                "Locus:\t82071_300\thas a significant blast hit:\t>gb|M1| Senecio chloroplast psbA gene\n"
                "Locus:\t16157_6\thas a significant blast hit:\t>gb|L1| Packera 28S ribosomal RNA gene\n"
                "49275_9\tref|NC_2|\t100\t122\t0\t0\t1\t122\t1\t122\t1e-50\t200\tPackera rbcL chloroplast gene\n",
                encoding="utf-8",
            )
            cp_s, its_s, _ = skim_mod._collect_organellar_seeds(gdir)
            self.assertTrue({"75949_5", "82071_300", "49275_9"}.issubset(cp_s.keys()))
            self.assertTrue({"49053_5", "52797_44", "68356_89"}.isdisjoint(cp_s.keys()))
            self.assertIn("16157_6", its_s)
            self.assertTrue({"52797_44", "68356_89"}.isdisjoint(its_s.keys()))
            self.assertEqual(skim_mod._extract_species_from_seed_key("AF161693.1|Packera_sp Senecio loratifolius"), "Packera loratifolia")
            by_id = {r["genetic_sample_id"]: r for r in skim_mod._detect_cytonuclear_discordance({"PACAURFK1M": "C1", "EREHIENCS1M": "C1"}, {})}
            self.assertEqual((by_id["EREHIENCS1M"]["nuclear_taxon"], by_id["EREHIENCS1M"]["cytonuclear_discordant"]), ("Erechtites hieraciifolius", 0))

    def test_public_backbone_digit_prefix_and_phytop_node_tree_lookup(self) -> None:
        """Verify Dryad file routing, .env auth defaults, and guarded Warm-R digit-prefix, PhyTop lookup, & MAD QC."""
        self.assertEqual(PUBLIC_DL_MOD._route_dryad_file(Path("PhyTop_results/MO_phytop.info.tsv"), "PhyTop_results.zip"), "public_packera_PhyTop_MO.tsv")
        self.assertEqual(PUBLIC_DL_MOD._route_dryad_file(Path("PhyTop_results/MO_phytop.tree"), "PhyTop_results.zip"), "public_packera_PhyTop_MO.tre")
        self.assertEqual(PUBLIC_DL_MOD._route_dryad_file(Path("PhyTop_results/ASTRAL-Pro_phytop.info.tsv"), "PhyTop_results.zip"), "public_packera_PhyTop_ASTRALPro.tsv")
        self.assertEqual(PUBLIC_DL_MOD._route_dryad_file(Path("RelTime_trees/S1outTest-HybPiper-df_nexus.tre"), "RelTime_trees.zip"), "public_packera_RelTime_S1outTest-HybPiper-df_nexus.tre")
        self.assertIsNone(PUBLIC_DL_MOD._route_dryad_file(Path("__MACOSX/PhyTop_results/._MO_phytop.tree"), "PhyTop_results.zip"))

        with tempfile.TemporaryDirectory() as td:
            tmp = Path(td)
            dotenv_file = tmp / ".env"
            dotenv_file.write_text("DRYAD_CLIENT_ID=dotenv_cid\nexport DRYAD_CLIENT_SECRET='dotenv_csec'\n", encoding="utf-8")
            with mock.patch.dict(os.environ, {"DRYAD_API_TOKEN": "keep_existing"}, clear=True):
                PUBLIC_DL_MOD._load_dotenv_defaults([dotenv_file])
                self.assertEqual((os.environ.get("DRYAD_API_TOKEN"), os.environ.get("DRYAD_CLIENT_ID"), os.environ.get("DRYAD_CLIENT_SECRET")), ("keep_existing", "dotenv_cid", "dotenv_csec"))

            if not rscript_available():
                self.skipTest("Rscript unavailable or PLANTMATICS_FAST_AUDIT=1")

            geno_dir, out_dir = tmp / "genomics", tmp / "out"
            geno_dir.mkdir(parents=True, exist_ok=True); out_dir.mkdir(parents=True, exist_ok=True)
            (geno_dir / "public_packera_MO_Astral_lpp.tre").write_text("(((Packera_2aurea:0.02,Packera_2obovata:0.02):0.05,Packera_dubia:0.07):0.15,Senecio_vulgaris:0.22);\n", encoding="utf-8")
            (geno_dir / "public_packera_PhyTop_MO.tre").write_text("(((Packera_2aurea:0.02,Packera_2obovata:0.02)N4:0.05,Packera_dubia:0.07)N3:0.15,Senecio_vulgaris:0.22);\n", encoding="utf-8")
            (geno_dir / "public_packera_PhyTop_MO.tsv").write_text("node\tILS_index\tIH_index\nN4\t0.92\t0.04\nN3\t0.85\t0.11\n", encoding="utf-8")

            r_script = tmp / "check_backbone_and_phytop.R"
            r_script.write_text(
                f"source({str(PROJECT_ROOT / 'scripts' / 'analysis' / 'step06_public_backbone.R')!r})\n"
                f"source({str(PROJECT_ROOT / 'scripts' / 'analysis' / 'genomic_triage_utils.R')!r})\n"
                "got <- .normalize_public_packera_taxon(c('Emilia 2fosbergii', 'AF161693.1|Packera_sp Senecio loratifolius', 'L33238.1|Senecio strepthanthifolius', 'Packera millefolia', 'Packera memmingeri', 'Packera \u00d7memmingeri', 'Packera sp.', 'Packera cf. aurea', 'Packera sp. nov.'))\n"
                "exp <- c('Emilia fosbergii', 'Packera loratifolia', 'Packera streptanthifolia', 'Packera millefolium', 'Packera xmemmingeri', 'Packera xmemmingeri', '', '', 'Packera sp. nov.')\n"
                "stopifnot(identical(got, exp))\n"
                "tr <- ape::read.tree(text = '((Packera_millefolium:1,Packera_xmemmingeri:1):1,(Packera_aurea:1,Packera_cana:3):1);')\n"
                "tr$tip.label <- .normalize_public_packera_taxon(tr$tip.label)\n"
                "bg <- select_backbone_guided_outgroup(c('Packera millefolia', 'Packera memmingeri', 'Packera aurea', 'Packera cana'), NULL, list(trees = list(MO = tr)))\n"
                "stopifnot(!is.null(bg), bg$outgroup == 'Packera cana', all(c('Packera millefolia', 'Packera memmingeri', 'Packera aurea') %in% bg$tree_order))\n"
                "stopifnot(identical(.align_packera_orthography(c('Packera millefolium', 'Packera xmemmingeri'), c('Packera millefolia', 'Packera memmingeri')), c('Packera millefolia', 'Packera memmingeri')))\n"
                "df <- data.frame(genetic_sample_id = c('S1', 'S2', 'S3', 'S4'), species_standardized = c('Packera aurea', 'Packera obovata', 'Packera dubia', 'Senecio vulgaris'), stringsAsFactors = FALSE)\n"
                f"res <- evaluate_public_backbone_and_organellar(df, {str(geno_dir)!r}, {str(out_dir)!r})\n"
                "cdf <- res$concordance_df\n"
                "stopifnot(all(cdf$in_mo_backbone), cdf$normalized_taxon[cdf$taxon == 'Senecio vulgaris'] == 'Senecio vulgaris')\n"
                "stopifnot(abs(cdf$public_phytop_ils_index[cdf$taxon == 'Packera aurea'] - 0.92) < 1e-4, abs(cdf$public_phytop_ih_index[cdf$taxon == 'Packera aurea'] - 0.04) < 1e-4)\n"
                "stopifnot(abs(cdf$public_phytop_ils_index[cdf$taxon == 'Packera obovata'] - 0.92) < 1e-4)\n"
                "stopifnot(abs(cdf$public_phytop_ils_index[cdf$taxon == 'Packera dubia'] - 0.85) < 1e-4, abs(cdf$public_phytop_ih_index[cdf$taxon == 'Packera dubia'] - 0.11) < 1e-4)\n"
                "tr_mad <- ape::read.tree(text = '((Packera_dubia_acc1:0.02,Packera_dubia_acc2:0.028):0.02,((Packera_aurea:0.022,Packera_obovata:0.024):0.02,(Packera_anonyma:0.025,Packera_paupercula:0.026):0.02):0.02,Packera_rogue:0.95,Erechtites_hieraciifolius:0.95);')\n"
                "qc_mad <- screen_public_tree_rogue_tips(tr_mad, tree_label = 'mad_test')\n"
                "stopifnot('Erechtites_hieraciifolius' %in% qc_mad$tree$tip.label, !('Packera_rogue' %in% qc_mad$tree$tip.label), all(c('Packera_dubia_acc1', 'Packera_dubia_acc2') %in% qc_mad$tree$tip.label))\n"
                "cat('PHYTOP_NODE_AND_MAD_OK\\n')\n",
                encoding="utf-8",
            )
            proc = run_warm_r_cmd(["Rscript", str(r_script)], stage_key="step06", cwd=PROJECT_ROOT, env={"PYTHONPATH": str(PROJECT_ROOT)}, skip_figures=True)
            self.assertEqual(proc.returncode, 0, f"R backbone/PhyTop/MAD check failed:\n{proc.stdout}\n{proc.stderr}")
            self.assertIn("PHYTOP_NODE_AND_MAD_OK", proc.stdout)

    def test_step06_public_backbone_outgroup_and_sqlite_columns(self) -> None:
        """Verify focal tree non-shadowing, public backbone outgroup polarization, and 6 SQLite columns."""
        if not rscript_available():
            self.skipTest("Rscript unavailable or PLANTMATICS_FAST_AUDIT=1")

        def _seq(seed_int: int, length: int = 120) -> str:
            return "".join(np.random.default_rng(seed_int).choice(list("ACGT"), size=length).tolist())

        cp_dub, cp_aur = _seq(303101), _seq(303102)
        its_dub, its_aur = _seq(303103), _seq(303104)
        waxy_dub, waxy_aur = _seq(303105), _seq(303106)

        with tempfile.TemporaryDirectory() as td:
            tmp = Path(td)
            geno_dir, parsed_dir, out_dir = tmp / "genomics", tmp / "parsed", tmp / "step06_out"
            tables_dir, figs_dir = tmp / "tables", tmp / "figures"
            (tmp / "data").mkdir(parents=True, exist_ok=True)
            for d in (geno_dir, parsed_dir, out_dir, tables_dir, figs_dir):
                d.mkdir(parents=True, exist_ok=True)
            db_path, run_id = tmp / "data" / "plantmatics.sqlite", "run_5563"

            sample_defs = [
                ("PACTOMSCG1M", "NCU_DUB_1", "Packera dubia", 35.12, -79.22, [1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 0.5, 0, 0, 0, 0, 0, 0.5, 0, 0, 0, 0, 0, 0, 0]),
                ("PACTOMVAL2M", "NCU_DUB_2", "Packera dubia", 35.14, -79.19, [1, 1, 1, 0.5, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0.5, 0, 0, 0, 0, 0, 0]),
                ("PACSCHNCB1M", "NCU_SCH_1", "Packera schweinitziana", 35.82, -82.91, [0.5, 0, 0, 0, 1, 1, 1, 1, 1, 0.5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]),
                ("PACSCHNCB2M", "NCU_SCH_2", "Packera schweinitziana", 35.85, -82.88, [0, 0.5, 0, 0, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]),
                ("PACMILVAL1M", "NCU_MIL_1", "Packera millefolium", 35.91, -82.11, [0, 0, 0, 0, 0.5, 0, 0, 0, 1, 1, 0, 0, 1, 1, 1, 1, 0.5, 0, 0, 0, 0, 0, 0, 0]),
                ("PACMILMIT1M", "NCU_MIL_2", "Packera millefolium", 35.93, -82.13, [0, 0, 0, 0, 0, 0.5, 0, 0, 1, 0.5, 0, 0, 1, 1, 1, 1, 0, 0.5, 0, 0, 0, 0, 0, 0]),
                ("PACAURFK1M", "NCU_AUR_1", "Packera aurea", 35.60, -81.90, [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 0.5]),
                ("PACAURFK2M", "NCU_AUR_2", "Packera aurea", 35.63, -81.87, [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 0.5, 1]),
            ]
            conn = get_db_connection(db_path=db_path, init_schema=True, workspace_root=tmp)
            try:
                conn.execute("INSERT INTO pipeline_runs (run_id, run_dir, intraspecific_mode) VALUES (?, ?, 0);", (run_id, str(tmp)))
                for sid, cat, sp, lat, lon, _ in sample_defs:
                    conn.execute("INSERT INTO vouchers (catalogNumber, scientificName, species_standardized, decimalLatitude, decimalLongitude, determiner_tier) VALUES (?, ?, ?, ?, ?, 'Tier_1_Gold');", (cat, sp, sp, lat, lon))
                    conn.execute("INSERT INTO plant_individuals (catalogNumber, plant_individual_id, rosette_centroid_x) VALUES (?, 0, 100.0);", (cat,))
                    conn.execute("INSERT INTO run_active_cohort (run_id, catalogNumber, plant_individual_id, species_standardized, determiner_tier, authority_weight) VALUES (?, ?, 0, ?, 'Tier_1_Gold', 2.0);", (run_id, cat, sp))
                    pub_sp = sp if sid in ("PACTOMSCG1M", "PACAURFK1M") else None
                    conn.execute("INSERT INTO genomic_samples (genetic_sample_id, catalogNumber, linkage_tier, species_standardized, geno_quality_pass, mapped_reads, public_plastome_match, public_its_match, waxy_lineage_match) VALUES (?, ?, 'G1_Direct_Voucher', ?, 1, 450000, ?, ?, ?);", (sid, cat, sp, pub_sp, pub_sp, pub_sp))
                conn.commit()
            finally:
                conn.close()

            loci = [f"L_{i:02d}" for i in range(24)]
            with open(parsed_dir / "allele_dosage_matrix.tsv", "w", encoding="utf-8", newline="") as fh:
                wr = csv.writer(fh, delimiter="\t"); wr.writerow(["genetic_sample_id"] + loci)
                for sid, _, _, _, _, vals in sample_defs:
                    wr.writerow([sid] + [f"{v:.4f}" for v in vals])
            with open(parsed_dir / "parsed_genomic_samples.csv", "w", encoding="utf-8", newline="") as fh:
                wr = csv.writer(fh)
                wr.writerow([
                    "genetic_sample_id", "project_id", "catalogNumber", "linkage_tier", "taxon_code", "taxon_raw",
                    "species_standardized", "population_id", "replicate_pair_id", "replicate_suffix",
                    "is_putative_novel_or_hybrid_query", "decimalLatitude", "decimalLongitude", "hq_reads",
                    "mapped_reads", "pct_mapped", "pct_repetitive", "n_missing_gt", "missing_genotype_rate",
                    "observed_het", "n_called_loci", "geno_quality_pass", "hdplot_retained_loci",
                    "inferred_cytotype", "nquire_delta_logl", "plastid_clade", "its_ribotype", "cytonuclear_discordant",
                ])
                for idx, (sid, cat, sp, lat, lon, _) in enumerate(sample_defs, start=1):
                    wr.writerow([sid, sid, cat, "G1_Direct_Voucher", sid[3:6], sp, sp, f"POP_{idx:02d}", sid[:-1], "M", "0", str(lat), str(lon), "450000", "410000", "91.1", "3.0", "0", "0.00", "0.20", "24", "1", "24", "2x", "1.1", f"cp_{sid[3:6]}", f"ITS_{sid[3:6]}", "0"])

            for rel_sd in ("snpsaurus_nextrad", "public_nuclear_backbone", "public_phytop_conflict", "public_organellar_and_loci"):
                (geno_dir / rel_sd).mkdir(parents=True, exist_ok=True)
            raxml_tree = geno_dir / "snpsaurus_nextrad" / "RAxML_bestTree.test"
            raxml_tree.write_text("(((PACTOMSCG1M:0.01,PACTOMVAL2M:0.01):0.04,(PACSCHNCB1M:0.01,PACSCHNCB2M:0.01):0.04):0.05,(PACMILVAL1M:0.01,PACMILMIT1M:0.01):0.10,(PACAURFK1M:0.01,PACAURFK2M:0.01):0.25);\n", encoding="utf-8")
            mo_tree = geno_dir / "public_nuclear_backbone" / "public_packera_MO_Astral_lpp.tre"
            mo_tree.write_text("(((Packera_dubia:0.04,Packera_schweinitziana:0.04):0.06,Packera_millefolium:0.10):0.15,Packera_aurea:0.25);\n", encoding="utf-8")
            (geno_dir / "public_nuclear_backbone" / "public_packera_ASTRALPro_Astral_lpp.tre").write_text("(((Packera_dubia:0.05,Packera_schweinitziana:0.05):0.06,Packera_millefolium:0.11):0.14,Packera_aurea:0.25);\n", encoding="utf-8")
            (geno_dir / "public_phytop_conflict" / "public_packera_PhyTop_MO.tsv").write_text("taxon\tils_index\tih_index\nPackera dubia\t0.21\t0.08\nPackera schweinitziana\t0.19\t0.05\nPackera millefolium\t0.14\t0.03\nPackera aurea\t0.11\t0.02\n", encoding="utf-8")
            (geno_dir / "public_organellar_and_loci" / "public_packera_18complete_plastomes.fasta").write_text(f">MK_CP_DUB Packera dubia\n{cp_dub}\n>MK_CP_AUR Packera aurea\n{cp_aur}\n", encoding="utf-8")
            (geno_dir / "public_organellar_and_loci" / "public_packera_its_ets.fasta").write_text(f">MT_ITS_DUB Packera dubia\n{its_dub}\n>MT_ITS_AUR Packera aurea\n{its_aur}\n", encoding="utf-8")
            (geno_dir / "public_organellar_and_loci" / "public_packera_waxy_nuclear.fasta").write_text(f">JX_WAXY_DUB Packera dubia\n{waxy_dub}\n>JX_WAXY_AUR Packera aurea\n{waxy_aur}\n", encoding="utf-8")
            (geno_dir / "public_organellar_and_loci" / "skimmed_plastid_consensus.fasta").write_text(f">PACTOMSCG1M|clade=cp_1\n{cp_dub}\n>PACAURFK1M|clade=cp_2\n{cp_aur}\n", encoding="utf-8")
            (geno_dir / "public_organellar_and_loci" / "skimmed_its_consensus.fasta").write_text(f">PACTOMSCG1M|ribotype=ITS_1\n{its_dub}\n>PACAURFK1M|ribotype=ITS_2\n{its_aur}\n", encoding="utf-8")
            (geno_dir / "public_organellar_and_loci" / "skimmed_waxy_consensus.fasta").write_text(f">PACTOMSCG1M|waxy=Waxy_1\n{waxy_dub}\n>PACAURFK1M|waxy=Waxy_2\n{waxy_aur}\n", encoding="utf-8")
            (parsed_dir / "genomic_manifest.json").write_text(json.dumps({"status": "ingested", "n_samples_total": len(sample_defs), "n_loci_retained": 24, "tree_files": [str(mo_tree), str(raxml_tree)]}), encoding="utf-8")

            proc = run_warm_r_cmd(
                [
                    "Rscript", str(PROJECT_ROOT / "scripts" / "analysis" / "06_genomic_phylogenomics.R"),
                    "--genomics-dir", str(geno_dir), "--raw-genomics-dir", str(geno_dir),
                    "--parsed-dir", str(parsed_dir), "--out-dir", str(out_dir), "--tables-dir", str(tables_dir),
                    "--figure-a", str(figs_dir / "fig06a.png"), "--figure-b", str(figs_dir / "fig06b.png"),
                    "--figure-c", str(figs_dir / "fig06c.png"), "--run-id", run_id, "--db", str(db_path),
                ],
                stage_key="step06", cwd=PROJECT_ROOT,
                env={"PLANTMATICS_DB_PATH": str(db_path), "PLANTMATICS_GENOMICS_DIR": str(geno_dir), "PYTHONPATH": str(PROJECT_ROOT)},
            )
            self.assertEqual(proc.returncode, 0, f"Step 06 public backbone run failed:\n{proc.stdout}\n{proc.stderr}")
            fig6c_pub = figs_dir / "fig06c.png"
            self.assertTrue(fig6c_pub.is_file() and fig6c_pub.stat().st_size > 10240)
            with open(fig6c_pub, "rb") as fh_p:
                self.assertEqual(fh_p.read(8), b"\x89PNG\r\n\x1a\n")
            self.assertEqual(list(figs_dir.glob("*.pdf")), [])

            with open(out_dir / "genomic_clade_cohesion_summary.csv", encoding="utf-8") as fh:
                self.assertGreaterEqual(len(list(csv.DictReader(fh))), 4, "Focal RAxML_bestTree.test must drive clade cohesion")
            with open(out_dir / "genomic_assignments.csv", encoding="utf-8") as fh:
                assign_rows = {r["genetic_sample_id"]: r for r in csv.DictReader(fh)}
            self.assertTrue(bool(assign_rows["PACTOMSCG1M"].get("nearest_tree_sister_taxon")))
            with open(out_dir / "introgressive_fbranch_summary.csv", encoding="utf-8") as fh:
                fb_rows = list(csv.DictReader(fh))
            self.assertTrue(len(fb_rows) > 0 and all(r["polarization_source"] == "public_nuclear_backbone" for r in fb_rows))
            for req_csv in ("public_backbone_concordance.csv", "public_organellar_placement.csv"):
                self.assertTrue((out_dir / req_csv).is_file() and (tables_dir / req_csv).is_file())

            with sqlite3.connect(str(db_path)) as conn_chk:
                gla_row = conn_chk.execute(
                    "SELECT public_backbone_clade, public_phytop_ils_index, public_phytop_ih_index FROM genetic_lineage_assignments WHERE run_id = ? AND genetic_sample_id = 'PACTOMSCG1M';",
                    (run_id,),
                ).fetchone()
                self.assertIsNotNone(gla_row)
                self.assertEqual(gla_row[0], "MO:Packera dubia+Packera schweinitziana")
                self.assertAlmostEqual(float(gla_row[1]), 0.21, places=4)
                self.assertAlmostEqual(float(gla_row[2]), 0.08, places=4)

                mv_row = conn_chk.execute(
                    "SELECT public_plastome_match, public_its_match, waxy_lineage_match, public_backbone_clade, public_phytop_ils_index, public_phytop_ih_index FROM v_integrative_triage_master WHERE run_id = ? AND catalogNumber = 'NCU_DUB_1';",
                    (run_id,),
                ).fetchone()
                self.assertIsNotNone(mv_row)
                self.assertEqual(mv_row[:4], ("Packera dubia", "Packera dubia", "Packera dubia", "MO:Packera dubia+Packera schweinitziana"))
                self.assertAlmostEqual(float(mv_row[4]), 0.21, places=4)
                self.assertAlmostEqual(float(mv_row[5]), 0.08, places=4)


if __name__ == "__main__":
    unittest.main()
