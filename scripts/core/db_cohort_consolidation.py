"""
Non-destructive SQLite healing of phantom pid=0 stubs and false multi-plant DBSCAN splits.
"""

from __future__ import annotations

import sqlite3
from typing import Any, Dict, List, Optional

from scripts.core.db_conn_pragmas import _ensure_dynamic_table_migrated

def consolidate_false_multi_plant_splits_db(conn: sqlite3.Connection) -> Dict[str, int]:
    """
    Non-destructively heal pre-existing Step 02 / preflight multi-plant artifacts in SQLite
    without requiring GPU re-segmentation:
      1. Remove all-NULL phantom `(catalogNumber, 0)` stubs in `plant_individuals` when a
         voucher already has segmented `plant_individual_id != 0` rows and zero `pid = 0`
         leaves, reproductive rows, or OTU morphometric records.
      2. Reunite false single-plant DBSCAN splits in `extracted_leaves` (cauline-only upper
         stem clusters, 1-leaf splayed outliers, and narrow same-column vertical splits)
         onto the primary basal `plant_individual_id`, coalescing per-individual capitula
         and phenology counts while preserving genuine multi-rosette sheets (`>= 2` basal
         leaves per independent rosette).
    """
    stats = {"phantom_zero_stubs_removed": 0, "false_split_individuals_merged": 0}
    _ensure_dynamic_table_migrated(conn, "voucher_reproductive_phenology")
    with conn:
        cur_stub = conn.execute(
            """
            DELETE FROM plant_individuals
            WHERE plant_individual_id = 0
              AND rosette_centroid_x IS NULL
              AND crown_x IS NULL
              AND petiole_length_px IS NULL
              AND petiole_length_mm IS NULL
              AND blade_length_px IS NULL
              AND blade_length_mm IS NULL
              AND EXISTS (
                  SELECT 1 FROM plant_individuals pi2
                  WHERE pi2.catalogNumber = plant_individuals.catalogNumber
                    AND pi2.plant_individual_id != 0
              )
              AND NOT EXISTS (
                  SELECT 1 FROM extracted_leaves el
                  WHERE el.catalogNumber = plant_individuals.catalogNumber
                    AND el.plant_individual_id = 0
              )
              AND NOT EXISTS (
                  SELECT 1 FROM voucher_reproductive_phenology vrp
                  WHERE vrp.catalogNumber = plant_individuals.catalogNumber
                    AND vrp.plant_individual_id = 0
              )
              AND NOT EXISTS (
                  SELECT 1 FROM otu_morphometrics om
                  WHERE om.catalogNumber = plant_individuals.catalogNumber
                    AND om.plant_individual_id = 0
              );
            """
        )
        stats["phantom_zero_stubs_removed"] = int(max(0, cur_stub.rowcount or 0))

        multi_Rows = conn.execute(
            """
            SELECT
                el.catalogNumber,
                el.plant_individual_id,
                COUNT(*) AS total_leaves,
                SUM(CASE WHEN LOWER(TRIM(COALESCE(el.leaf_zone, 'basal'))) != 'cauline' THEN 1 ELSE 0 END) AS basal_leaves,
                COALESCE(pi.crown_x, pi.rosette_centroid_x) AS anchor_x
            FROM extracted_leaves el
            LEFT JOIN plant_individuals pi
              ON pi.catalogNumber = el.catalogNumber
             AND pi.plant_individual_id = el.plant_individual_id
            WHERE el.catalogNumber IN (
                SELECT catalogNumber
                FROM extracted_leaves
                GROUP BY catalogNumber
                HAVING COUNT(DISTINCT plant_individual_id) > 1
            )
            GROUP BY el.catalogNumber, el.plant_individual_id
            ORDER BY el.catalogNumber, el.plant_individual_id;
            """
        ).fetchall()

        by_cat: Dict[str, List[Dict[str, Any]]] = {}
        for r in multi_Rows:
            by_cat.setdefault(str(r[0]), []).append(
                {
                    "pid": int(r[1]),
                    "total": int(r[2] or 0),
                    "basal": int(r[3] or 0),
                    "x": float(r[4]) if r[4] is not None else None,
                }
            )

        for cat, clusters in by_cat.items():
            if len(clusters) < 2:
                continue
            # True independent rosette hubs have >= 2 basal leaves
            rosette_hubs = [c for c in clusters if c["basal"] >= 2]
            if not rosette_hubs:
                rosette_hubs = [max(clusters, key=lambda c: (c["basal"], c["total"], -c["pid"]))]

            # Merge vertically stacked rosette clusters sitting in the same narrow shoot column (<= 320 px)
            kept_hubs: List[Dict[str, Any]] = []
            merge_map: Dict[int, int] = {}
            for hub in sorted(rosette_hubs, key=lambda c: (-c["basal"], -c["total"], c["pid"])):
                merged_to: Optional[int] = None
                if hub["x"] is not None:
                    for kh in kept_hubs:
                        if kh["x"] is not None and abs(hub["x"] - kh["x"]) <= 320.0:
                            merged_to = kh["pid"]
                            break
                if merged_to is not None:
                    merge_map[hub["pid"]] = merged_to
                else:
                    kept_hubs.append(hub)

            kept_pids = {h["pid"] for h in kept_hubs}
            for c in clusters:
                pid_c = c["pid"]
                if pid_c in kept_pids or pid_c in merge_map:
                    continue
                # Non-rosette (cauline-only) or 1-leaf stray cluster -> attach to nearest/primary hub
                if c["x"] is not None and any(h["x"] is not None for h in kept_hubs):
                    xh_cands = [h for h in kept_hubs if h["x"] is not None]
                    target_hub = min(xh_cands, key=lambda h: abs(float(c["x"]) - float(h["x"])))
                else:
                    target_hub = kept_hubs[0]
                merge_map[pid_c] = int(target_hub["pid"])

            # Also prune orphan leaf-less non-zero plant_individuals on this segmented sheet
            all_pi_pids = [
                int(row[0])
                for row in conn.execute(
                    "SELECT plant_individual_id FROM plant_individuals WHERE catalogNumber = ?;",
                    (cat,),
                ).fetchall()
            ]
            leaf_pids = {c["pid"] for c in clusters}
            primary_kept = int(kept_hubs[0]["pid"])
            for opid in all_pi_pids:
                if opid not in leaf_pids and opid != primary_kept:
                    merge_map[opid] = primary_kept

            if not merge_map:
                continue

            for src_pid, dst_pid in sorted(merge_map.items()):
                if src_pid == dst_pid:
                    continue
                conn.execute(
                    """
                    INSERT OR IGNORE INTO plant_individuals (catalogNumber, plant_individual_id)
                    VALUES (?, ?);
                    """,
                    (cat, dst_pid),
                )
                src_leaves = conn.execute(
                    "SELECT leaf_id, shape_id FROM extracted_leaves WHERE catalogNumber = ? AND plant_individual_id = ?;",
                    (cat, src_pid),
                ).fetchall()
                max_lid_row = conn.execute(
                    "SELECT COALESCE(MAX(leaf_id), 0) FROM extracted_leaves WHERE catalogNumber = ? AND plant_individual_id = ?;",
                    (cat, dst_pid),
                ).fetchone()
                next_lid = int((max_lid_row or [0])[0] or 0)
                existing_dst_lids = {
                    int(r[0])
                    for r in conn.execute(
                        "SELECT leaf_id FROM extracted_leaves WHERE catalogNumber = ? AND plant_individual_id = ?;",
                        (cat, dst_pid),
                    ).fetchall()
                }
                for raw_lid, raw_sid in src_leaves:
                    lid_cur = int(raw_lid)
                    new_lid = lid_cur
                    if new_lid in existing_dst_lids:
                        next_lid = max(next_lid + 1, lid_cur + 1)
                        new_lid = next_lid
                    existing_dst_lids.add(new_lid)
                    new_sid = f"{cat}_p{dst_pid}_leaf{new_lid}"
                    sid_taken = conn.execute(
                        """
                        SELECT 1 FROM extracted_leaves
                        WHERE shape_id = ?
                          AND NOT (catalogNumber = ? AND plant_individual_id = ? AND leaf_id = ?);
                        """,
                        (new_sid, cat, src_pid, lid_cur),
                    ).fetchone()
                    if sid_taken:
                        new_sid = f"{cat}_p{dst_pid}_leaf{new_lid}_m{src_pid}"
                    conn.execute(
                        """
                        UPDATE extracted_leaves
                        SET plant_individual_id = ?, leaf_id = ?, shape_id = ?
                        WHERE catalogNumber = ? AND plant_individual_id = ? AND leaf_id = ?;
                        """,
                        (dst_pid, new_lid, new_sid, cat, src_pid, lid_cur),
                    )

                # Coalesce reproductive/phenological counts onto dst_pid before removing src_pid
                vrp_src = conn.execute(
                    """
                    SELECT capitula_count, involucre_height_px, involucre_width_px,
                           involucre_height_mm, involucre_width_mm, AR_inv,
                           total_capitula, n_anthesis, n_bud, n_fruiting, n_senescent,
                           voucher_phenological_state, is_verified_anthesis,
                           reproductive_stage_confidence, dominant_pappus_ratio,
                           corymb_head_count_est, radiate_head_fraction
                    FROM voucher_reproductive_phenology
                    WHERE catalogNumber = ? AND plant_individual_id = ?;
                    """,
                    (cat, src_pid),
                ).fetchone()
                if vrp_src is not None:
                    vrp_dst = conn.execute(
                        """
                        SELECT capitula_count, involucre_height_px, involucre_width_px,
                               involucre_height_mm, involucre_width_mm, AR_inv,
                               total_capitula, n_anthesis, n_bud, n_fruiting, n_senescent,
                               voucher_phenological_state, is_verified_anthesis,
                               reproductive_stage_confidence, dominant_pappus_ratio,
                               corymb_head_count_est, radiate_head_fraction
                        FROM voucher_reproductive_phenology
                        WHERE catalogNumber = ? AND plant_individual_id = ?;
                        """,
                        (cat, dst_pid),
                    ).fetchone()
                    if vrp_dst is None:
                        conn.execute(
                            """
                            UPDATE voucher_reproductive_phenology
                            SET plant_individual_id = ?
                            WHERE catalogNumber = ? AND plant_individual_id = ?;
                            """,
                            (dst_pid, cat, src_pid),
                        )
                    else:
                        def _sum_opt(a: Any, b: Any) -> Optional[int]:
                            if a is None and b is None:
                                return None
                            return int(a or 0) + int(b or 0)

                        s_anth = _sum_opt(vrp_dst[7], vrp_src[7])
                        d_state = vrp_dst[11]
                        s_state = vrp_src[11]
                        merged_state = (
                            s_state
                            if (
                                (not d_state or str(d_state).lower() in {"sterile", "vegetative", "unknown"})
                                and s_state
                            )
                            else (d_state or s_state)
                        )
                        merged_ver = (
                            1
                            if (bool(vrp_dst[12]) or bool(vrp_src[12]) or (s_anth is not None and s_anth >= 1))
                            else (vrp_dst[12] if vrp_dst[12] is not None else vrp_src[12])
                        )
                        conn.execute(
                            """
                            UPDATE voucher_reproductive_phenology
                            SET capitula_count = ?,
                                involucre_height_px = COALESCE(involucre_height_px, ?),
                                involucre_width_px = COALESCE(involucre_width_px, ?),
                                involucre_height_mm = COALESCE(involucre_height_mm, ?),
                                involucre_width_mm = COALESCE(involucre_width_mm, ?),
                                AR_inv = COALESCE(AR_inv, ?),
                                total_capitula = ?,
                                n_anthesis = ?,
                                n_bud = ?,
                                n_fruiting = ?,
                                n_senescent = ?,
                                voucher_phenological_state = ?,
                                is_verified_anthesis = ?,
                                reproductive_stage_confidence = COALESCE(reproductive_stage_confidence, ?),
                                dominant_pappus_ratio = COALESCE(dominant_pappus_ratio, ?),
                                corymb_head_count_est = ?,
                                radiate_head_fraction = COALESCE(radiate_head_fraction, ?)
                            WHERE catalogNumber = ? AND plant_individual_id = ?;
                            """,
                            (
                                _sum_opt(vrp_dst[0], vrp_src[0]),
                                vrp_src[1],
                                vrp_src[2],
                                vrp_src[3],
                                vrp_src[4],
                                vrp_src[5],
                                _sum_opt(vrp_dst[6], vrp_src[6]),
                                s_anth,
                                _sum_opt(vrp_dst[8], vrp_src[8]),
                                _sum_opt(vrp_dst[9], vrp_src[9]),
                                _sum_opt(vrp_dst[10], vrp_src[10]),
                                merged_state,
                                merged_ver,
                                vrp_src[13],
                                vrp_src[14],
                                _sum_opt(vrp_dst[15], vrp_src[15]) or 0,
                                vrp_src[16],
                                cat,
                                dst_pid,
                            ),
                        )
                        conn.execute(
                            "DELETE FROM voucher_reproductive_phenology WHERE catalogNumber = ? AND plant_individual_id = ?;",
                            (cat, src_pid),
                        )

                # Coalesce linear scape/petiole metrics onto dst_pid when dst_pid lacked scape height
                conn.execute(
                    """
                    UPDATE plant_individuals
                    SET plant_height_px = COALESCE(
                            plant_height_px,
                            (SELECT s.plant_height_px FROM plant_individuals s WHERE s.catalogNumber = ? AND s.plant_individual_id = ?)
                        ),
                        plant_height_mm = COALESCE(
                            plant_height_mm,
                            (SELECT s.plant_height_mm FROM plant_individuals s WHERE s.catalogNumber = ? AND s.plant_individual_id = ?)
                        ),
                        plant_height_cm = COALESCE(
                            plant_height_cm,
                            (SELECT s.plant_height_cm FROM plant_individuals s WHERE s.catalogNumber = ? AND s.plant_individual_id = ?)
                        ),
                        scape_to_basal_leaf_ratio = COALESCE(
                            scape_to_basal_leaf_ratio,
                            (SELECT s.scape_to_basal_leaf_ratio FROM plant_individuals s WHERE s.catalogNumber = ? AND s.plant_individual_id = ?)
                        )
                    WHERE catalogNumber = ? AND plant_individual_id = ?;
                    """,
                    (cat, src_pid, cat, src_pid, cat, src_pid, cat, src_pid, cat, dst_pid),
                )

                for tier_b_tbl in ("leaf_efa_harmonics", "otu_morphometrics", "spatial_rf_predictions", "curation_triage", "run_active_cohort"):
                    conn.execute(
                        f"DELETE FROM {tier_b_tbl} WHERE catalogNumber = ? AND plant_individual_id = ?;",
                        (cat, src_pid),
                    )
                conn.execute(
                    "DELETE FROM plant_individuals WHERE catalogNumber = ? AND plant_individual_id = ?;",
                    (cat, src_pid),
                )
                stats["false_split_individuals_merged"] += 1

    return stats

__all__ = (
    "consolidate_false_multi_plant_splits_db",
)
