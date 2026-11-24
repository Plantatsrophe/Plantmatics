#!/usr/bin/env python3
"""
Batched Warm-R integration tests for Step 06 phylogenomics (`06_genomic_phylogenomics.R`),
out-of-sample Gower LORO-CV PCoA projection, ParaLoss-1272 weighted allele-sharing distance,
continuous-dosage SNPs2CF, `persist_step06_to_sqlite()` SQLite persistence, 300-DPI PNG
rendering (zero `.pdf` leaks), and Step 07 Phase 2 consensus & decision-matrix integration.
"""

from __future__ import annotations

import csv
import json
import sqlite3
import sys
import tempfile
import unittest
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from scripts.core.db import get_db_connection
from scripts.core.db_schema import PLANTMATICS_SCHEMA_USER_VERSION
from scripts.tests._r_session import rscript_available, run_warm_r_cmd


def _build_synthetic_vcf_text(
    samples: list[str],
    locus_specs: list[tuple[str, list[tuple[str, int, int, int]]]],
) -> str:
    lines = [
        "##fileformat=VCFv4.2",
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
        '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">',
        '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic Depths">',
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(samples),
    ]
    for pos, (lid, cells) in enumerate(locus_specs, start=1):
        fmt_cells = [f"{gt}:{dp}:{r},{a}" for gt, dp, r, a in cells]
        lines.append(f"contig_{pos}\t{100 + pos}\t{lid}\tA\tG\t99\tPASS\t.\tGT:DP:AD\t" + "\t".join(fmt_cells))
    return "\n".join(lines) + "\n"


class TestStep06PhylogenomicsRPipeline(unittest.TestCase):
    """Batched Warm-R verification of Step 06 phylogenomics, audit contracts, and Step 07 Phase 2 synthesis."""

    def test_batched_r_step06_and_step07_phase2_consensus_and_triage(self) -> None:
        if not rscript_available():
            self.skipTest("Rscript not found in PATH")

        # Verify LORO-CV uses fold-internal Gower out-of-sample projection (zero transductive X_dapc reuse)
        r_src_06 = (PROJECT_ROOT / "scripts" / "analysis" / "06_genomic_phylogenomics.R").read_text(encoding="utf-8")
        self.assertNotIn("X_dapc[ref_idx[tr_sub]", r_src_06, "LORO-CV must project test samples out-of-sample via fold-internal Gower PCoA")

        with tempfile.TemporaryDirectory() as td:
            tmp = Path(td)
            empty_geno, pop_geno, parsed_dir = tmp / "empty_geno", tmp / "pop_geno", tmp / "parsed"
            step06_empty_out, step06_out = tmp / "step06_empty_out", tmp / "step06_out"
            tables_dir, figs_dir, rep_dir = tmp / "tables", tmp / "figures", tmp / "reports"
            for d in (empty_geno, pop_geno, parsed_dir, step06_empty_out, step06_out, tables_dir, figs_dir, rep_dir):
                d.mkdir(parents=True, exist_ok=True)
            db_path = tmp / "packera_batch.sqlite"

            vouchers = [
                ("NCU00445952", 0, "Packera schweinitziana", "Packera schweinitziana", "Fuller", "NCU", "Tier_1_Gold", 35.82, -82.91, "NC", "Haywood", "2024-06-15", "Basal_Rosette", 6, 0, 1, 0.5, "Packera schweinitziana", 0.94, "Packera schweinitziana", 0.91, 166.0, 2.0, "Measured", 1, 0, 1, "Basal_Rosette", 0.4, 0.85, 18.0, 0, -0.8, 0.2, 0.1),
                ("NCU_SCH_2", 0, "Packera schweinitziana", "Packera schweinitziana", "Fuller", "NCU", "Tier_1_Gold", 35.85, -82.88, "NC", "Haywood", "2024-06-18", "Basal_Rosette", 5, 0, 1, 0.4, "Packera schweinitziana", 0.92, "Packera schweinitziana", 0.88, 169.0, 4.0, "Measured", 1, 0, 1, "Basal_Rosette", 0.5, 0.80, 17.5, 0, -0.7, 0.1, 0.0),
                ("NCU_DUB_1", 0, "Packera dubia", "Packera dubia", "Fuller", "NCU", "Tier_1_Gold", 35.12, -79.22, "NC", "Moore", "2024-04-20", "Basal_Rosette", 7, 0, 1, 0.2, "Packera dubia", 0.95, "Packera dubia", 0.90, 110.0, -1.0, "Measured", 1, 0, 1, "Basal_Rosette", 0.3, 0.90, 16.0, 0, 1.1, -0.1, 0.2),
                ("NCU_DUB_2", 0, "Packera dubia", "Packera dubia", "Fuller", "NCU", "Tier_1_Gold", 35.14, -79.19, "NC", "Moore", "2024-04-22", "Basal_Rosette", 6, 0, 1, 0.3, "Packera dubia", 0.93, "Packera dubia", 0.89, 112.0, 1.5, "Measured", 1, 0, 1, "Basal_Rosette", 0.4, 0.88, 15.5, 0, 1.0, 0.0, 0.1),
                ("NCU_MIL_1", 0, "Packera millefolium", "Packera millefolium", "Fuller", "NCU", "Tier_1_Gold", 35.91, -82.11, "NC", "Mitchell", "2024-05-20", "Basal_Rosette", 5, 0, 1, 0.1, "Packera millefolium", 0.96, "Packera millefolium", 0.93, 140.0, 0.5, "Measured", 1, 0, 1, "Basal_Rosette", 0.2, 0.91, 19.0, 0, 0.1, 1.2, -0.3),
                ("NCU_MIL_2", 0, "Packera millefolium", "Packera millefolium", "Fuller", "NCU", "Tier_1_Gold", 35.93, -82.13, "NC", "Mitchell", "2024-05-22", "Basal_Rosette", 6, 0, 1, 0.2, "Packera millefolium", 0.95, "Packera millefolium", 0.91, 142.0, -0.5, "Measured", 1, 0, 1, "Basal_Rosette", 0.3, 0.89, 18.5, 0, 0.2, 1.1, -0.2),
                ("NCU_F1_HYB", 0, "Packera dubia", "Packera dubia", "Fuller", "NCU", "Tier_2_Silver", 35.48, -81.05, "NC", "Burke", "2024-05-15", "Basal_Rosette", 4, 1, 1, 2.0, "Packera schweinitziana", 0.51, "Packera dubia", 0.55, 136.0, 12.0, "Measured", 1, 0, 1, "Basal_Rosette", 1.8, 0.35, 3.2, 0, 0.2, 0.1, 0.0),
                ("NCU_UNSEQ", 0, "Packera dubia", "Packera dubia", "Fuller", "NCU", "Tier_1_Gold", 35.10, -79.25, "NC", "Moore", "2024-04-18", "Basal_Rosette", 5, 0, 1, 0.1, "Packera dubia", 0.94, "Packera dubia", 0.92, 108.0, 0.0, "Measured", 1, 0, 1, "Basal_Rosette", 0.3, 0.85, 15.0, 0, 1.05, -0.05, 0.1),
            ]
            v_headers = [
                "catalogNumber", "plant_individual_id", "species_standardized", "species_raw",
                "recordedBy", "institutionCode", "determiner_tier", "decimalLatitude", "decimalLongitude",
                "stateProvince", "county", "eventDate", "leaf_zone", "basal_leaf_count",
                "is_morph_outlier", "gmm_cluster", "gmm_bayes_factor", "cda_predicted_taxon",
                "cda_posterior_prob", "edaphic_best_fit_taxon", "edaphic_prob", "hopkins_adjusted_doy",
                "delta_doy", "reproductive_trait_status", "is_anthesis", "is_passive", "is_basal_morph",
                "morph_tier_source", "type_edaphic_z_score", "open_set_novelty_p", "min_congener_mahal_d2",
                "cryptic_cluster_cohesion_pass", "PC1", "PC2", "PC3",
            ]
            vouchers_csv = tables_dir / "curated_vouchers.csv"
            with open(vouchers_csv, "w", encoding="utf-8", newline="") as fh:
                wr = csv.writer(fh); wr.writerow(v_headers); wr.writerows(vouchers)
            with open(tables_dir / "multimodal_conflict_flags.csv", "w", encoding="utf-8", newline="") as fh:
                wr = csv.writer(fh)
                wr.writerow(["catalogNumber", "plant_individual_id", "edaphic_best_fit_taxon", "edaphic_prob", "type_edaphic_z_score", "delta_doy", "hopkins_adjusted_doy"])
                for row in vouchers:
                    wr.writerow([row[0], row[1], row[19], row[20], row[28], row[22], row[21]])
            with open(tables_dir / "morphometric_flags.csv", "w", encoding="utf-8", newline="") as fh:
                wr = csv.writer(fh)
                wr.writerow(["catalogNumber", "plant_individual_id", "is_morph_outlier", "gmm_cluster", "gmm_bayes_factor", "cda_predicted_taxon", "cda_posterior_prob", "PC1", "PC2", "PC3"])
                for row in vouchers:
                    wr.writerow([row[0], row[1], row[14], row[15], row[16], row[17], row[18], row[32], row[33], row[34]])

            conn = get_db_connection(db_path=db_path, init_schema=True, workspace_root=tmp)
            try:
                for row in vouchers:
                    conn.execute(
                        "INSERT INTO vouchers (catalogNumber, recordedBy, institutionCode, scientificName, species_standardized, decimalLatitude, decimalLongitude, stateProvince, county, eventDate, determiner_tier) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?);",
                        (row[0], row[4], row[5], row[2], row[2], row[7], row[8], row[9], row[10], row[11], row[6]),
                    )
                    conn.execute("INSERT INTO plant_individuals (catalogNumber, plant_individual_id, rosette_centroid_x) VALUES (?, ?, ?);", (row[0], row[1], 100.0))
                for sp_c, r_cy, m2n, mpl in (("Packera dubia", "2x, 4x", 44, 2), ("Packera schweinitziana", "4x", 44, 4)):
                    conn.execute("INSERT OR REPLACE INTO taxon_cytology_summary (species_standardized, reported_cytotypes, modal_sporophytic_2n, modal_ploidy_level) VALUES (?, ?, ?, ?);", (sp_c, r_cy, m2n, mpl))
                conn.commit()
            finally:
                conn.close()

            s_defs = [
                ("NCU00445952", "NCU00445952", "G1_Direct_Voucher", "Packera schweinitziana", 0, 35.82, -82.91, "4x", "Clade_Sch", "ITS_Sch", 0, "SCH"),
                ("PACSCHNCB1M", "NCU_SCH_2", "G1_Direct_Voucher", "Packera schweinitziana", 0, 35.85, -82.88, "4x", "Clade_Sch", "ITS_Sch", 0, "SCH"),
                ("PACSCHNCB1P", "NCU_SCH_2", "G2_Population_Voucher_Proxy", "Packera schweinitziana", 0, 35.85, -82.88, "4x", "Clade_Sch", "ITS_Sch", 0, "SCH"),
                ("PACTOMSCG1M", "NCU_DUB_1", "G1_Direct_Voucher", "Packera dubia", 0, 35.12, -79.22, "2x", "Clade_Dub", "ITS_Dub", 0, "DUB"),
                ("PACTOMVAL2M", "NCU_DUB_2", "G1_Direct_Voucher", "Packera dubia", 0, 35.14, -79.19, "2x", "Clade_Dub", "ITS_Dub", 0, "DUB"),
                ("PACMILVAL1M", "NCU_MIL_1", "G1_Direct_Voucher", "Packera millefolium", 0, 35.91, -82.11, "4x", "Clade_Mil", "ITS_Mil", 0, "MIL"),
                ("PACMILMIT1M", "NCU_MIL_2", "G1_Direct_Voucher", "Packera millefolium", 0, 35.93, -82.13, "4x", "Clade_Mil", "ITS_Mil", 0, "MIL"),
                ("PACF1HYBM", "NCU_F1_HYB", "G1_Direct_Voucher", "Packera dubia", 1, 35.48, -81.05, "2x", "Clade_Sch", "ITS_Dub", 1, "F1_DUB_SCH"),
            ]
            loci = [f"L_{i:02d}" for i in range(30)]
            dos_rows = []
            for sid, *_, mode in s_defs:
                vals = [f"{(1.0 if li < 10 else 0.0) if mode == 'SCH' else ((1.0 if 10 <= li < 20 else 0.0) if mode == 'DUB' else ((1.0 if li >= 20 else 0.0) if mode == 'MIL' else (0.5 if li < 20 else 0.0))):.4f}" for li in range(30)]
                dos_rows.append([sid] + vals)
            with open(parsed_dir / "allele_dosage_matrix.tsv", "w", encoding="utf-8", newline="") as fh:
                wr = csv.writer(fh, delimiter="\t"); wr.writerow(["genetic_sample_id"] + loci); wr.writerows(dos_rows)
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
                for idx, (sid, cat, tier, sp, is_q, lat, lon, cyt, plast, its_r, cytonuc, _) in enumerate(s_defs, start=1):
                    wr.writerow([
                        sid, sid, cat, tier, sid[3:6], sp, sp, f"POP_{idx:02d}" if not sid.endswith(("M", "P")) else f"POP_{sid[:-1]}",
                        sid[:-1] if sid[-1] in ("M", "P") else sid, sid[-1] if sid[-1] in ("M", "P") else "",
                        str(is_q), str(lat), str(lon), "450000", "410000", "91.1", "3.2", "0", "0.00", "0.22", "30", "1",
                        "30", cyt, "5.2", plast, its_r, str(cytonuc),
                    ])
            with open(parsed_dir / "hdplot_locus_diagnostics.tsv", "w", encoding="utf-8", newline="") as fh:
                wr = csv.writer(fh, delimiter="\t")
                wr.writerow(["locus_id", "n_called_ref", "n_hets_ref", "H_obs", "het_ref_reads", "het_alt_reads", "het_total_reads", "D_read_ratio", "hdplot_keep", "filter_reason", "is_paraloss1272_single_copy", "locus_paralogy_weight"])
                for li, lid in enumerate(loci):
                    is_pl = 1 if li < 10 else 0
                    wr.writerow([lid, "8", "3", f"{0.18 + 0.01 * (li % 10):.3f}", "30", "28", "58", f"{-1.2 + 0.1 * (li % 10):.3f}", "1", "keep", str(is_pl), "1.5" if is_pl else ("0.5" if li >= 20 else "1.0")])
                wr.writerow(["L_PARA_1", "8", "7", "0.820", "45", "10", "55", "6.200", "0", "high_H_and_D", "0", "0.5"])
                wr.writerow(["L_PARA_2", "8", "6", "0.740", "12", "42", "54", "-5.400", "0", "high_H_and_D", "0", "0.5"])
            with open(parsed_dir / "locus_weights.tsv", "w", encoding="utf-8", newline="") as fh:
                wr = csv.writer(fh, delimiter="\t")
                wr.writerow(["locus_id", "paraloss_status", "locus_paralogy_weight"])
                for li, lid in enumerate(loci):
                    st_lbl, wt_s = ("single_copy", "1.5") if li < 10 else (("multi_copy", "0.5") if li >= 20 else ("unannotated", "1.0"))
                    wr.writerow([lid, st_lbl, wt_s])
            with open(parsed_dir / "cytonuclear_discordance.csv", "w", encoding="utf-8", newline="") as fh:
                wr = csv.writer(fh)
                wr.writerow(["genetic_sample_id", "nuclear_taxon", "plastid_clade", "its_ribotype", "cytonuclear_discordant", "public_plastome_match", "public_its_match", "waxy_lineage_match"])
                for sid, _, _, sp, _, _, _, _, plast, its_r, cytonuc, _ in s_defs:
                    wr.writerow([sid, sp, plast, its_r, str(cytonuc), "Packera schweinitziana" if sid == "PACF1HYBM" else sp, sp, sp])
            with open(parsed_dir / "genomic_manifest.json", "w", encoding="utf-8") as fh:
                json.dump({"status": "ingested", "n_samples_total": len(s_defs), "n_loci_retained": 30}, fh)
            (pop_geno / "packera_batch.vcf").write_text(
                _build_synthetic_vcf_text([r[0] for r in s_defs], [("L_00", [("1/1" if r[11] == "SCH" else "0/0", 30, 10, 20) for r in s_defs])]),
                encoding="utf-8",
            )
            (pop_geno / "RAxML_bestTree.batch").write_text(
                "(((NCU00445952:0.01,PACSCHNCB1M:0.01,PACSCHNCB1P:0.01):0.05,(PACTOMSCG1M:0.01,PACTOMVAL2M:0.01):0.05):0.02,PACF1HYBM:0.03,(PACMILVAL1M:0.01,PACMILMIT1M:0.01):0.08);\n",
                encoding="utf-8",
            )

            fig6a = figs_dir / "fig06a_genomic_pcoa_tree_concordance.png"
            fig6b = figs_dir / "fig06b_reticulation_gdi_landscape.png"
            fig6c = figs_dir / "fig06c_cytonuclear_paralog_cytotype_audit.png"
            fig7a_base, fig7a_4st = figs_dir / "fig07a_base_3stream.png", figs_dir / "fig07a_integrative_taxonomic_synthesis.png"
            fig7b = figs_dir / "fig07b_taxon_complex_delimitation_network.png"
            q_base_csv, q_4st_csv = tables_dir / "triage_queue_base3.csv", tables_dir / "triage_queue.csv"
            rep_base_md, rep_4st_md = rep_dir / "report_base3.md", rep_dir / "Packera_Taxonomic_Revision_Summary.md"
            cplx_csv, aff_csv, nbr_txt, html_out = rep_dir / "empirical_taxon_complexes.csv", rep_dir / "pairwise_taxon_affinity_matrix.csv", rep_dir / "focal_auto_neighbors.txt", rep_dir / "Triage_Dashboard_Interactive.html"

            # 1. Step 06 empty no-op check via warm R CLI
            proc_empty = run_warm_r_cmd(
                [
                    "Rscript", "scripts/analysis/06_genomic_phylogenomics.R",
                    "--genomics-dir", str(empty_geno), "--parsed-dir", str(tmp / "nonexistent_parsed"),
                    "--out-dir", str(step06_empty_out), "--tables-dir", str(tmp / "empty_tables"),
                    "--output-figures", str(tmp / "empty_figs"), "--run-id", "batch_smoke_empty",
                    "--db", str(db_path), "--demo",
                ],
                stage_key="step06", cwd=PROJECT_ROOT,
                env={"PLANTMATICS_DB_PATH": str(db_path), "PACKERA_SQLITE_DB": str(db_path), "PLANTMATICS_GENOMICS_DIR": str(empty_geno), "RUN_ID": "batch_smoke_empty"},
                skip_figures=False,
            )
            self.assertEqual(proc_empty.returncode, 0, f"Step 06 empty no-op failed:\n{proc_empty.stdout}\n{proc_empty.stderr}")
            self.assertFalse((step06_empty_out / "genomic_assignments.csv").exists())

            # 2. Step 06 populated Phase 2 execution via warm R CLI
            env_pop = {"PLANTMATICS_DB_PATH": str(db_path), "PACKERA_SQLITE_DB": str(db_path), "PLANTMATICS_GENOMICS_DIR": str(pop_geno), "RUN_ID": "batch_smoke"}
            proc_pop = run_warm_r_cmd(
                [
                    "Rscript", "scripts/analysis/06_genomic_phylogenomics.R",
                    "--genomics-dir", str(pop_geno), "--raw-genomics-dir", str(pop_geno),
                    "--parsed-dir", str(parsed_dir), "--out-dir", str(step06_out), "--tables-dir", str(tables_dir),
                    "--output-figures", str(figs_dir), "--figure-a", str(fig6a), "--figure-b", str(fig6b), "--figure-c", str(fig6c),
                    "--run-id", "batch_smoke", "--db", str(db_path), "--demo",
                ],
                stage_key="step06", cwd=PROJECT_ROOT, env=env_pop, skip_figures=False,
            )
            self.assertEqual(proc_pop.returncode, 0, f"Step 06 populated run failed:\n{proc_pop.stdout}\n{proc_pop.stderr}")

            # 3. Consolidated warm-R check: Phase 2 CSVs, ParaLoss weighting, continuous SNPs2CF, persist_step06_to_sqlite, Step 07 synthesis
            batch_r = tmp / "verify_batched_phase2.R"
            batch_r.write_text(f"""
Sys.setenv(PLANTMATICS_DB_PATH = {json.dumps(str(db_path))}, PACKERA_SQLITE_DB = {json.dumps(str(db_path))}, PLANTMATICS_GENOMICS_DIR = {json.dumps(str(pop_geno))}, RUN_ID = "batch_smoke")
source("scripts/core/db_utils.R")
source("scripts/analysis/step06_phylogenomics_utils.R")
source("scripts/analysis/step06_fbranch_quartets.R")

req_phase2_files <- c(
  "genomic_assignments.csv", "genomic_clade_cohesion_summary.csv", "genomic_model_metrics.csv",
  "genomic_multimodal_congruence_summary.csv", "coancestry_matrix_finerad.csv", "hybrid_index_q12_triangle.csv",
  "introgressive_fbranch_summary.csv", "quartet_concordance_factors_snps2cf.csv",
  "landscape_mmrr_feems_summary.csv", "feems_migration_edges.csv"
)
for (fn in req_phase2_files) {{
  fp <- file.path({json.dumps(str(step06_out))}, fn)
  if (!file.exists(fp)) stop(paste("Missing Step 06 Phase 2 file:", fp))
}}
stopifnot(file.exists({json.dumps(str(fig6a))}), file.info({json.dumps(str(fig6a))})$size > 5000)
stopifnot(file.exists({json.dumps(str(fig6b))}), file.info({json.dumps(str(fig6b))})$size > 5000)
stopifnot(file.exists({json.dumps(str(fig6c))}), file.info({json.dumps(str(fig6c))})$size > 10240)

ga <- utils::read.csv(file.path({json.dumps(str(step06_out))}, "genomic_assignments.csv"), stringsAsFactors = FALSE)
rownames(ga) <- ga$genetic_sample_id
stopifnot("NCU00445952" %in% rownames(ga), identical(ga["NCU00445952", "geno_predicted_taxon"], "Packera schweinitziana"))
stopifnot(identical(ga["PACF1HYBM", "hybrid_generation_class"], "F1_Hybrid"))
stopifnot(is.finite(ga["PACF1HYBM", "interspecific_ancestry_q12"]) && ga["PACF1HYBM", "interspecific_ancestry_q12"] >= 0.80)
stopifnot(ga["PACMILVAL1M", "hybrid_generation_class"] != "F1_Hybrid")
stopifnot(is.na(ga["PACMILVAL1M", "interspecific_ancestry_q12"]) || ga["PACMILVAL1M", "interspecific_ancestry_q12"] < 0.25)

mmrr_df <- utils::read.csv(file.path({json.dumps(str(step06_out))}, "landscape_mmrr_feems_summary.csv"), stringsAsFactors = FALSE)
stopifnot(isTRUE(as.logical(mmrr_df$small_n_guardrail_active[1L])))
stopifnot(grepl("OutFLANK", mmrr_df$excluded_large_n_methods[1L]), grepl("dadi", mmrr_df$excluded_large_n_methods[1L]), grepl("HZAR", mmrr_df$excluded_large_n_methods[1L]))

# Out-of-sample Gower PCoA projection (no transductive leakage) + ParaLoss-1272 locus_weights.tsv loading
pts_oos <- rbind(A1 = c(0, 0), A2 = c(0.1, 0), A3 = c(0, 0.1), B1 = c(4, 0), B2 = c(4.1, 0), B3 = c(4, 0.1), T_A = c(0.05, 0.05))
d_oos <- as.matrix(stats::dist(pts_oos))
g_oos <- .project_pcoa_gower_oos(d_oos[1:6, 1:6, drop = FALSE], d_oos[7L, 1:6, drop = FALSE], k_max = 3L)
stopifnot(nrow(g_oos$X_train) == 6L, nrow(g_oos$X_test) == 1L, all(is.finite(g_oos$X_test)))
cA <- colMeans(g_oos$X_train[1:3, , drop = FALSE]); cB <- colMeans(g_oos$X_train[4:6, , drop = FALSE])
stopifnot(sum((g_oos$X_test[1L, ] - cA)^2) < sum((g_oos$X_test[1L, ] - cB)^2))
lw_loaded <- .load_locus_weights_tsv({json.dumps(str(parsed_dir))}, c("L_00", "L_15", "L_25"))
stopifnot(abs(lw_loaded[["L_00"]] - 1.5) < 1e-6, abs(lw_loaded[["L_15"]] - 1.0) < 1e-6, abs(lw_loaded[["L_25"]] - 0.5) < 1e-6)

# ParaLoss-1272 weighted allele-sharing distance: 10 divergent single-copy (1.5) + 10 identical multi-copy (0.5) => 0.75 vs 0.50
dm_wt <- rbind(S1 = c(rep(0.0, 10L), rep(0.0, 10L)), S2 = c(rep(1.0, 10L), rep(0.0, 10L)))
w_vec <- c(rep(1.5, 10L), rep(0.5, 10L))
d_unw <- .compute_allele_sharing_distance(dm_wt)$dist[1L, 2L]
d_wtd <- .compute_allele_sharing_distance(dm_wt, w_vec)$dist[1L, 2L]
stopifnot(abs(d_unw - 0.50) < 1e-6, is.finite(d_wtd), abs(d_wtd - 0.75) < 1e-4)

# Continuous-dosage SNPs2CF: 0.5 heterozygotes must NOT round to 0
dos_het_cf <- cbind(
  matrix(rep(c(0.5, 0.5, 0.0, 0.0), 40L), nrow = 4L, ncol = 40L),
  matrix(rep(c(0.5, 0.0, 0.5, 0.0), 8L), nrow = 4L, ncol = 8L),
  matrix(rep(c(0.5, 0.0, 0.0, 0.5), 8L), nrow = 4L, ncol = 8L)
)
pass_cf <- data.frame(genetic_sample_id = c("QA", "QB", "QC", "QD"), species_standardized = c("TaxA", "TaxB", "TaxC", "TaxD"), stringsAsFactors = FALSE)
rownames(dos_het_cf) <- pass_cf$genetic_sample_id
cf_het <- compute_quartet_concordance_snps2cf(pass_cf, dos_het_cf)$cf_df
stopifnot(nrow(cf_het) == 1L, cf_het$n_inform_loci[1L] > 0.0, cf_het$cf_12_34[1L] > cf_het$cf_13_24[1L], identical(cf_het$major_topology[1L], "TaxA,TaxB|TaxC,TaxD"))

# persist_step06_to_sqlite alias & upsert of reconciled organellar + concept/quarantine columns into genomic_samples
stopifnot(exists("persist_step06_to_sqlite", mode = "function"), identical(persist_step06_to_sqlite, persist_step06_genomics_to_sqlite))
p_df <- data.frame(
  genetic_sample_id = "NCU00445952", catalogNumber = "NCU00445952", geno_predicted_taxon = "Packera schweinitziana",
  geno_posterior_prob = 0.98, cytonuclear_discordant = 1L, public_plastome_match = "Packera dubia",
  public_its_match = "Packera schweinitziana", waxy_lineage_match = "Packera schweinitziana",
  exsiccatae_key = "fuller|101|2024-06-15", genbank_accession = "OQ123456", coordinateUncertainty = 250.0,
  is_centroid_or_inst_coord = 1L, coord_quarantine_reason = "county_centroid", typified_basionym = "Senecio robbinsii",
  concept_synonymy_status = "heterotypic_synonym", is_open_or_hybrid_concept = 1L, stringsAsFactors = FALSE
)
persist_step06_to_sqlite(p_df, run_id = "batch_persist", db_path = {json.dumps(str(db_path))}, demo = FALSE)

# Step 07 synthesis: 3-stream vs 4-stream invariance + Phase 2 rules
source("scripts/analysis/07_triage_dashboard_synthesis.R")
opts_3 <- parse_args_robust(c("--vouchers", {json.dumps(str(vouchers_csv))}, "--morphometrics", {json.dumps(str(tables_dir / "morphometric_flags.csv"))}, "--multimodal-flags", {json.dumps(str(tables_dir / "multimodal_conflict_flags.csv"))}, "--output-queue", {json.dumps(str(q_base_csv))}, "--output-plot", {json.dumps(str(fig7a_base))}, "--output-report", {json.dumps(str(rep_base_md))}, "--genomic-dir", "", "--demo"))
q3 <- build_triage_queue(opts_3)
opts_4 <- parse_args_robust(c("--vouchers", {json.dumps(str(vouchers_csv))}, "--morphometrics", {json.dumps(str(tables_dir / "morphometric_flags.csv"))}, "--multimodal-flags", {json.dumps(str(tables_dir / "multimodal_conflict_flags.csv"))}, "--output-queue", {json.dumps(str(q_4st_csv))}, "--output-plot", {json.dumps(str(fig7a_4st))}, "--output-report", {json.dumps(str(rep_4st_md))}, "--output-html", {json.dumps(str(html_out))}, "--output-complexes", {json.dumps(str(cplx_csv))}, "--output-affinity", {json.dumps(str(aff_csv))}, "--output-neighbors", {json.dumps(str(nbr_txt))}, "--output-complex-plot", {json.dumps(str(fig7b))}, "--genomic-dir", {json.dumps(str(step06_out))}, "--demo"))
q4 <- build_triage_queue(opts_4)
render_synthesis_figure_png(q4, opts_4$output_plot, output_complex_plot = opts_4$output_complex_plot)
generate_taxonomic_revision_report(q4, opts_4)

u3 <- q3[q3$catalogNumber == "NCU_UNSEQ", , drop = FALSE]; u4 <- q4[q4$catalogNumber == "NCU_UNSEQ", , drop = FALSE]
stopifnot(nrow(u3) == 1L, nrow(u4) == 1L, identical(as.character(u3$taxonomic_treatment_subtype), as.character(u4$taxonomic_treatment_subtype)), abs(as.numeric(u3$probabilistic_consensus_score) - as.numeric(u4$probabilistic_consensus_score)) < 1e-6)
f1_row <- q4[q4$catalogNumber == "NCU_F1_HYB", , drop = FALSE]
stopifnot(nrow(f1_row) == 1L, identical(as.character(f1_row$taxonomic_treatment_subtype), "Hybrid_Introgressant"))
stopifnot(file.exists(file.path({json.dumps(str(tables_dir))}, "genomic_field_triage.csv")), file.exists({json.dumps(str(fig7a_4st))}), file.info({json.dumps(str(fig7a_4st))})$size > 5000)

source("scripts/figures/step06c_figures.R")
pub_conc_demo <- data.frame(taxon = c("Packera dubia", "Packera schweinitziana", "Packera millefolium"), public_phytop_ils_index = c(0.21, 0.19, 0.14), public_phytop_ih_index = c(0.08, 0.05, 0.03), stringsAsFactors = FALSE)
render_fig06c_cytonuclear_paralog_cytotype_png({json.dumps(str(fig6c))}, {json.dumps(str(parsed_dir))}, {json.dumps(str(step06_out))}, ga, pub_concordance_df = pub_conc_demo, db_path = {json.dumps(str(db_path))})
render_step07a_integrative_synthesis_png(q3, out_plot_path = {json.dumps(str(fig7a_base))})
q4_cyto <- q4; q4_cyto$cytology_linkage_tier <- ifelse(q4_cyto$catalogNumber == "NCU_DUB_1", "C1_Direct_Cytovoucher", NA_character_)
render_step07a_integrative_synthesis_png(q4_cyto, out_plot_path = {json.dumps(str(fig7a_4st))})
stopifnot(file.exists({json.dumps(str(fig6c))}), file.info({json.dumps(str(fig6c))})$size > 10240)
stopifnot(file.exists({json.dumps(str(fig7a_base))}), file.info({json.dumps(str(fig7a_base))})$size > 10240)
stopifnot(file.exists({json.dumps(str(fig7a_4st))}), file.info({json.dumps(str(fig7a_4st))})$size > 10240)
rep_txt <- paste(readLines({json.dumps(str(rep_4st_md))}, warn = FALSE), collapse = "\\n")
stopifnot(grepl("## 5\\\\. Independent Step 06 Genomic & Phylogenomic Evidence Stream", rep_txt))
stopifnot(length(list.files({json.dumps(str(figs_dir))}, pattern = "\\\\.pdf$", ignore.case = TRUE)) == 0L)

make_intra_probe <- function(prefix, gdi_val, hyb_cls, q12_val) {{
  df_k2 <- data.frame(
    catalogNumber = paste0(prefix, "_", 1:5), plant_individual_id = 0L, species_standardized = rep("Packera dubia", 5L),
    species_raw = rep("Packera dubia", 5L), recordedBy = c("Fuller", "Weakley", "Sorrie", "Poindexter", "Krings"),
    institutionCode = c("NCU", "DUKE", "US", "NY", "GH"), determiner_tier = rep("Tier_1_Gold", 5L),
    decimalLatitude = c(35.10, 35.45, 35.80, 36.15, 36.50), decimalLongitude = c(-79.10, -79.55, -80.05, -80.55, -81.05),
    stateProvince = rep("NC", 5L), county = c("Moore", "Chatham", "Orange", "Ashe", "Watauga"),
    eventDate = c("2020-04-15", "2021-04-18", "2022-04-20", "2023-04-22", "2024-04-25"), leaf_zone = rep("Basal_Rosette", 5L),
    basal_leaf_count = rep(6L, 5L), is_basal_morph = rep(TRUE, 5L), is_passive = rep(FALSE, 5L), morph_tier_source = rep("Basal_Rosette", 5L),
    is_morph_outlier = rep(TRUE, 5L), gmm_cluster = rep(2L, 5L), gmm_bayes_factor = rep(12.0, 5L), cda_predicted_taxon = rep("Packera dubia", 5L),
    cda_posterior_prob = rep(0.45, 5L), edaphic_best_fit_taxon = rep("Packera dubia", 5L), edaphic_prob = rep(0.20, 5L),
    type_edaphic_z_score = rep(3.2, 5L), hopkins_adjusted_doy = rep(155.0, 5L), delta_doy = rep(35.0, 5L),
    reproductive_trait_status = rep("Measured", 5L), is_anthesis = rep(TRUE, 5L), open_set_novelty_p = rep(0.01, 5L), open_set_df = rep(5L, 5L),
    min_congener_mahal_d2 = rep(28.0, 5L), cryptic_cluster_cohesion_pass = rep(TRUE, 5L),
    PC1 = c(2.4, 2.5, 2.3, 2.6, 2.4), PC2 = c(-1.1, -1.0, -1.2, -1.1, -0.9), PC3 = c(0.4, 0.3, 0.5, 0.4, 0.4),
    geno_predicted_taxon = rep("Packera dubia", 5L), geno_linkage_tier = rep("G1_Direct_Voucher", 5L), geno_quality_pass = rep(TRUE, 5L),
    geno_posterior_prob = rep(0.96, 5L), coalescent_gdi = rep(gdi_val, 5L), hybrid_generation_class = rep(hyb_cls, 5L),
    interspecific_ancestry_q12 = rep(q12_val, 5L), is_geno_admixed = rep(FALSE, 5L), stringsAsFactors = FALSE
  )
  df_k1 <- data.frame(
    catalogNumber = paste0(prefix, "_BASE_", 1:3), plant_individual_id = 0L, species_standardized = rep("Packera dubia", 3L),
    species_raw = rep("Packera dubia", 3L), recordedBy = c("Radford", "Bozeman", "Leonard"),
    institutionCode = c("NCU", "NCU", "DUKE"), determiner_tier = rep("Tier_1_Gold", 3L),
    decimalLatitude = c(35.05, 35.12, 35.18), decimalLongitude = c(-79.05, -79.12, -79.18),
    stateProvince = rep("NC", 3L), county = c("Moore", "Hoke", "Richmond"),
    eventDate = c("1960-04-10", "1968-04-12", "1972-04-14"), leaf_zone = rep("Basal_Rosette", 3L),
    basal_leaf_count = rep(6L, 3L), is_basal_morph = rep(TRUE, 3L), is_passive = rep(FALSE, 3L), morph_tier_source = rep("Basal_Rosette", 3L),
    is_morph_outlier = rep(FALSE, 3L), gmm_cluster = rep(1L, 3L), gmm_bayes_factor = rep(0.2, 3L), cda_predicted_taxon = rep("Packera dubia", 3L),
    cda_posterior_prob = rep(0.92, 3L), edaphic_best_fit_taxon = rep("Packera dubia", 3L), edaphic_prob = rep(0.88, 3L),
    type_edaphic_z_score = rep(0.4, 3L), hopkins_adjusted_doy = rep(115.0, 3L), delta_doy = rep(0.0, 3L),
    reproductive_trait_status = rep("Measured", 3L), is_anthesis = rep(TRUE, 3L), open_set_novelty_p = rep(0.85, 3L), open_set_df = rep(5L, 3L),
    min_congener_mahal_d2 = rep(16.0, 3L), cryptic_cluster_cohesion_pass = rep(FALSE, 3L),
    PC1 = c(-0.5, -0.4, -0.6), PC2 = c(0.1, 0.0, 0.2), PC3 = c(0.0, 0.1, -0.1),
    geno_predicted_taxon = rep("Packera dubia", 3L), geno_linkage_tier = rep("G1_Direct_Voucher", 3L), geno_quality_pass = rep(TRUE, 3L),
    geno_posterior_prob = rep(0.96, 3L), coalescent_gdi = rep(gdi_val, 3L), hybrid_generation_class = rep("Pure_Parent_1", 3L),
    interspecific_ancestry_q12 = rep(0.02, 3L), is_geno_admixed = rep(FALSE, 3L), stringsAsFactors = FALSE
  )
  rbind(df_k2, df_k1)
}}
intra_oversplit <- apply_taxonomic_decision_matrix(make_intra_probe("OVERSPLIT", 0.12, "Pure_Parent_1", 0.02), intraspecific_mode = TRUE, gmm_two_delta_bic = 18.0)
stopifnot(grepl("Coalescent Oversplit Block", intra_oversplit$synthesis_rationale[1L]))
stopifnot(intra_oversplit$taxonomic_treatment_subtype[1L] %in% c("Ecophenotypic_Plasticity", "Confirmed_Species"))
intra_cryptic <- apply_taxonomic_decision_matrix(make_intra_probe("CRYPTIC", 0.88, "Divergent_NonHybrid_Lineage", 0.04), intraspecific_mode = TRUE, gmm_two_delta_bic = 18.0)
stopifnot(identical(as.character(intra_cryptic$taxonomic_treatment_subtype[1L]), "Putative_Cryptic_Species"))
cat("BATCHED_PHASE2_R_OK\\n")
""", encoding="utf-8")

            proc = run_warm_r_cmd(["Rscript", str(batch_r)], stage_key="step07", cwd=PROJECT_ROOT, env=env_pop, skip_figures=False)
            self.assertEqual(proc.returncode, 0, f"Batched R verification failed:\nSTDOUT:\n{proc.stdout}\nSTDERR:\n{proc.stderr}")
            self.assertIn("BATCHED_PHASE2_R_OK", proc.stdout)

            # Verify PNG magic headers and zero .pdf leaks
            for png_p, min_sz in ((fig6a, 5000), (fig6b, 5000), (fig6c, 10240), (fig7a_base, 10240), (fig7a_4st, 10240)):
                self.assertTrue(png_p.is_file(), f"Missing rendered PNG: {png_p}")
                self.assertGreater(png_p.stat().st_size, min_sz, f"PNG too small: {png_p}")
                with open(png_p, "rb") as fh_png:
                    self.assertEqual(fh_png.read(8), b"\x89PNG\r\n\x1a\n", f"Invalid PNG magic header: {png_p}")
            self.assertEqual(list(figs_dir.glob("*.pdf")), [])

            # Verify persist_step06_to_sqlite saved all organellar + concept/quarantine columns & schema version 20261127
            with sqlite3.connect(str(db_path)) as conn_chk:
                self.assertEqual(conn_chk.execute("PRAGMA user_version;").fetchone()[0], PLANTMATICS_SCHEMA_USER_VERSION)
                self.assertEqual(PLANTMATICS_SCHEMA_USER_VERSION, 20261127)
                gs_row = conn_chk.execute(
                    """
                    SELECT cytonuclear_discordant, public_plastome_match, public_its_match, waxy_lineage_match,
                           exsiccatae_key, genbank_accession, coordinateUncertainty, is_centroid_or_inst_coord,
                           coord_quarantine_reason, typified_basionym, concept_synonymy_status, is_open_or_hybrid_concept
                    FROM genomic_samples WHERE genetic_sample_id = 'NCU00445952';
                    """
                ).fetchone()
                self.assertIsNotNone(gs_row)
                self.assertEqual(
                    (gs_row[0], gs_row[1], gs_row[2], gs_row[3], gs_row[4], gs_row[5], float(gs_row[6]), gs_row[7], gs_row[8], gs_row[9], gs_row[10], gs_row[11]),
                    (1, "Packera dubia", "Packera schweinitziana", "Packera schweinitziana", "fuller|101|2024-06-15", "OQ123456", 250.0, 1, "county_centroid", "Senecio robbinsii", "heterotypic_synonym", 1),
                )


if __name__ == "__main__":
    unittest.main()
