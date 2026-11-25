#!/usr/bin/env bash
# scripts/pipeline/sh/pipe_preflight_cohort_a.sh - First half + runner for preflight subset filtering.

source "${PROJECT_ROOT}/scripts/pipeline/sh/pipe_preflight_cohort_b.sh"

_pipe_preflight_py_part_a() {
  cat <<'__PLANTMATICS_PREFLIGHT_PART_A__'
import math, os, re, shutil, sys
from pathlib import Path
import pandas as pd

run_dir = Path(sys.argv[1]).resolve()
root_data = Path(sys.argv[2]).resolve()
focal_raw = sys.argv[3].strip()
exclude_csv = sys.argv[4].strip()
intraspecific_mode = sys.argv[5].lower() == 'true'
gold_anchors = int(sys.argv[6])
determiner_tier_raw = sys.argv[7].strip()
exclude_misidentified = sys.argv[8].lower() == 'true' if len(sys.argv) > 8 else False
triage_source_raw = sys.argv[9].strip() if len(sys.argv) > 9 else 'ALL'
auto_neighbors_flag = sys.argv[10].lower() == 'true' if len(sys.argv) > 10 else False
min_affinity_val = float(sys.argv[11]) if len(sys.argv) > 11 and sys.argv[11] else 0.03
spatial_thin_km = float(sys.argv[12]) if len(sys.argv) > 12 and sys.argv[12] else 0.0
max_per_state = int(sys.argv[13]) if len(sys.argv) > 13 and sys.argv[13].strip() else None
max_per_region = int(sys.argv[14]) if len(sys.argv) > 14 and sys.argv[14].strip() else None
no_spatial_weights = sys.argv[15].lower() == 'true' if len(sys.argv) > 15 else False
reuse_canonical = sys.argv[16].lower() == 'true' if len(sys.argv) > 16 else False
custom_taxa_csv = sys.argv[17].strip() if len(sys.argv) > 17 else ''
seed_data = Path(sys.argv[18]).resolve() if len(sys.argv) > 18 and sys.argv[18].strip() else root_data
start_step_val = int(sys.argv[19]) if len(sys.argv) > 19 and sys.argv[19].strip().isdigit() else 1
usda_filter_arg = sys.argv[20].strip().lower() if len(sys.argv) > 20 else ''
usda_buffer_arg = float(sys.argv[21]) if len(sys.argv) > 21 and sys.argv[21].strip() else 50.0
usda_exempt_gold_arg = (sys.argv[22].strip().lower() != 'false') if len(sys.argv) > 22 and sys.argv[22].strip() else True

proj_root_arg = Path(sys.argv[23]).resolve() if len(sys.argv) > 23 and sys.argv[23].strip() else root_data.parent
for _cand_root in ([proj_root_arg, root_data.parent, Path.cwd()]):
    if (_cand_root / 'scripts').is_dir() and str(_cand_root) not in sys.path:
        sys.path.insert(0, str(_cand_root))

from scripts.core.config import resolve_canonical_run_id
from scripts.core.harvester import standardize_packera_taxon
from scripts.core.taxonomy_authority import select_cascading_outgroup_anchors, spatially_balance_and_thin_vouchers

ALLOWED_TREATMENT_SUBTYPES = {
    'Confirmed_Species', 'Confirmed_Type_Control', 'Type_Synonymy_Candidate',
    'Hybrid_Introgressant', 'Subspecies_Ecotype', 'Ecophenotypic_Plasticity', 'Putative_Cryptic_Species',
}

tier_map = {
    'tier_1_gold': 'Tier_1_Gold', 'gold': 'Tier_1_Gold', 'tier_1': 'Tier_1_Gold', 'tier1': 'Tier_1_Gold', '1': 'Tier_1_Gold',
    'tier_2_silver': 'Tier_2_Silver', 'silver': 'Tier_2_Silver', 'tier_2': 'Tier_2_Silver', 'tier2': 'Tier_2_Silver', '2': 'Tier_2_Silver',
    'tier_3_bronze': 'Tier_3_Bronze', 'bronze': 'Tier_3_Bronze', 'tier_3': 'Tier_3_Bronze', 'tier3': 'Tier_3_Bronze', '3': 'Tier_3_Bronze',
}
target_tier = tier_map.get(determiner_tier_raw.lower(), determiner_tier_raw) if determiner_tier_raw else ''

custom_taxa = set()
if custom_taxa_csv:
    for item in custom_taxa_csv.split(','):
        item = item.strip()
        if item:
            std_item = standardize_packera_taxon(item)
            if std_item:
                custom_taxa.add(std_item)

exclude_taxa = set()
if exclude_csv:
    for item in exclude_csv.split(','):
        item = item.strip()
        if item:
            std_ex = standardize_packera_taxon(item)
            if std_ex:
                exclude_taxa.add(std_ex)

tables_dir = run_dir / 'data' / 'tables'
csv_path = tables_dir / 'curated_vouchers.csv'

def _resolve_preflight_taxon(sn, sr, ss):
    for cand in (sn, sr, ss):
        if pd.notna(cand) and str(cand).strip():
            std = standardize_packera_taxon(cand)
            if std:
                return std
    return None

focal_std_req = standardize_packera_taxon(focal_raw) if focal_raw else None
requested_taxa = set(custom_taxa)
if focal_std_req:
    requested_taxa.add(focal_std_req)

# Restore requested-taxon vouchers that exist in the canonical SQLite vouchers table but are missing
# from the curated_vouchers.csv snapshot (taxon-restricted harvests UPSERT into SQLite but are never
# promoted to data/tables/curated_vouchers.csv, so --skip-harvest runs would otherwise silently drop them).
if requested_taxa and os.environ.get('PLANTMATICS_DEMO_MODE') != '1' and os.environ.get('PACKERA_DEMO_MODE') != '1' and seed_data == root_data:
    from scripts.core.db import get_db_connection as _get_db_conn_sync, reconcile_curated_vouchers_csv_with_db
    sync_db_path = os.environ.get('PLANTMATICS_DB_PATH') or os.environ.get('PACKERA_DB_PATH') or str(root_data / 'plantmatics.sqlite')
    if Path(sync_db_path).is_file():
        conn_sync = _get_db_conn_sync(db_path=sync_db_path, init_schema=True, workspace_root=root_data.parent)
        try:
            n_synced = reconcile_curated_vouchers_csv_with_db(conn_sync, csv_path, taxa=requested_taxa, taxon_resolver=_resolve_preflight_taxon)
        finally:
            conn_sync.close()
        if n_synced:
            print(f'[Preflight: SQLite Voucher Sync] Restored {n_synced} vouchers for requested taxa {sorted(requested_taxa)} from SQLite vouchers that were missing from curated_vouchers.csv.')

df = pd.read_csv(csv_path, low_memory=False)
initial_voucher_count = len(df)

non_demo_pf = (os.environ.get('PLANTMATICS_DEMO_MODE') != '1' and os.environ.get('PACKERA_DEMO_MODE') != '1')
target_db_path = os.environ.get('PLANTMATICS_DB_PATH') or os.environ.get('PACKERA_DB_PATH') or str(root_data / 'plantmatics.sqlite')
if non_demo_pf and Path(target_db_path).is_file() and not df.empty and 'catalogNumber' in df.columns:
    try:
        from scripts.core.db import get_db_connection as _get_db_hyd
        _conn_hyd = _get_db_hyd(db_path=target_db_path, init_schema=False, read_only=True, workspace_root=root_data.parent)
        try:
            _vcols = {r[1] for r in _conn_hyd.execute('PRAGMA table_info(vouchers)').fetchall()}
            _hyd_cols = [c for c in ('usda_symbol', 'usda_county_status', 'exsiccatae_key', 'exsiccatae_taxonomic_conflict', 'exsiccatae_gold_taxon') if c in _vcols]
            if _hyd_cols:
                _rows_h = _conn_hyd.execute(
                    'SELECT catalogNumber, ' + ', '.join(_hyd_cols) + ' FROM vouchers'
                ).fetchall()
                _hmap = {str(r[0]).strip(): dict(zip(_hyd_cols, r[1:])) for r in _rows_h if r[0] is not None}
                _cats_df = df['catalogNumber'].astype(str).str.strip().tolist()
                for _hc in _hyd_cols:
                    if _hc not in df.columns:
                        df[_hc] = [_hmap.get(_c, {}).get(_hc) for _c in _cats_df]
                    else:
                        _cur_vals = df[_hc].tolist()
                        _new_vals = []
                        for _c, _cv in zip(_cats_df, _cur_vals):
                            _is_empty = _cv is None or (isinstance(_cv, float) and math.isnan(_cv)) or str(_cv).strip() in ('', 'None', 'nan')
                            _hv = _hmap.get(_c, {}).get(_hc)
                            _new_vals.append(_hv if (_is_empty and _hv is not None and str(_hv).strip() not in ('', 'None', 'nan')) else _cv)
                        df[_hc] = _new_vals
        finally:
            _conn_hyd.close()
    except Exception:
        pass

sn_vals = df['scientificName'] if 'scientificName' in df.columns else [None] * len(df)
sr_vals = df['species_raw'] if 'species_raw' in df.columns else [None] * len(df)
ss_vals = df['species_standardized'] if 'species_standardized' in df.columns else (df['species'] if 'species' in df.columns else [None] * len(df))
df['species_standardized'] = [_resolve_preflight_taxon(sn, sr, ss) for sn, sr, ss in zip(sn_vals, sr_vals, ss_vals)]
unres_mask = df['species_standardized'].isna() | df['species_standardized'].astype(str).str.strip().isin({'', 'None', 'nan'})
if unres_mask.any():
    n_unres = int(unres_mask.sum())
    df = df[~unres_mask].copy()
    print(f'[Preflight: Taxon Hygiene] Dropped {n_unres} vouchers with unresolvable or open-nomenclature taxon strings.')

from scripts.core.taxonomy_authority import apply_corpus_authority_tiers, evaluate_eastern_na_cohort_filter
df = apply_corpus_authority_tiers(df)

eff_geo_scope = (sys.argv[24] if len(sys.argv) > 24 else 'eastern_na').strip().lower() or 'eastern_na'
inc_west_pf = (sys.argv[25] if len(sys.argv) > 25 else 'false').strip().lower() == 'true'
has_geo_cols = (
    ('stateProvince' in df.columns and df['stateProvince'].astype(str).str.strip().replace({'nan': '', 'None': ''}).ne('').any())
    or ('decimalLongitude' in df.columns and pd.to_numeric(df['decimalLongitude'], errors='coerce').notna().any())
)
if eff_geo_scope == 'eastern_na' and not inc_west_pf and not custom_taxa and has_geo_cols:
    keep_east_pf, east_pf_stats = evaluate_eastern_na_cohort_filter(df, standardize_packera_taxon)
    if not keep_east_pf.all():
        before_east_pf = len(df)
        df = df.loc[keep_east_pf].reset_index(drop=True)
        ret_sp_n = len(east_pf_stats.get('retained_species', []))
        drp_sp_list = east_pf_stats.get('dropped_species', [])
        print(
            f'[Preflight: Cis-Mississippi Eastern NA Scope] Retained {len(df)} / {before_east_pf} vouchers '
            f'across {ret_sp_n} cis-Mississippi Eastern species '
            f'(dropped non-Eastern/western species: {drp_sp_list}).'
        )

if custom_taxa:
    before_custom_n = len(df)
    df = df[df['species_standardized'].isin(custom_taxa)].copy()
    print(f'[Preflight: Custom Taxa] Restricted strictly to {sorted(custom_taxa)}: retained {len(df)} / {before_custom_n} vouchers.')

# Fail fast when an explicitly requested taxon has zero vouchers (never silently analyze a reduced taxon set)
if requested_taxa:
    present_taxa_all = set(df['species_standardized'].dropna().astype(str).str.strip())
    missing_requested = sorted(t for t in requested_taxa if t not in present_taxa_all)
    if missing_requested:
        print(
            f'[ERROR] Requested taxa {missing_requested} have zero vouchers in {csv_path} or the SQLite vouchers table. '
            'Harvest them first (re-run without --skip-harvest using a new --run-id or --force), or remove them from --taxon / --focal-taxon.',
            file=sys.stderr,
        )
        sys.exit(1)

from scripts.core.db import ensure_pipeline_run, get_db_connection, populate_run_active_cohort, query_misidentified_vouchers_from_db

def _norm_pid(val):
    if val is None or (isinstance(val, float) and math.isnan(val)):
        return 0
    s = str(val).strip()
    if not s or s.lower() in ('nan', 'none', 'na', 'null'):
        return 0
    try:
        return int(float(s))
    except (TypeError, ValueError):
        return 0

g1_dna_cats = set()
if non_demo_pf and Path(target_db_path).is_file():
    try:
        _conn_g1 = get_db_connection(db_path=target_db_path, init_schema=False, read_only=True, workspace_root=root_data.parent)
        try:
            _tbls_g1 = {r[0] for r in _conn_g1.execute('SELECT name FROM sqlite_master WHERE type=?', ('table',)).fetchall()}
            if 'genomic_samples' in _tbls_g1:
                _g1_cols = {r[1] for r in _conn_g1.execute('PRAGMA table_info(genomic_samples)').fetchall()}
                if {'catalogNumber', 'linkage_tier', 'geno_quality_pass'}.issubset(_g1_cols):
                    _q_g1 = 'SELECT DISTINCT catalogNumber FROM genomic_samples WHERE linkage_tier IN (?, ?) AND COALESCE(geno_quality_pass, 0) = 1 AND TRIM(COALESCE(catalogNumber, ?)) != ?'
                    _direct_g1 = {str(r[0]).strip() for r in _conn_g1.execute(_q_g1, ('G1_Direct_Voucher', 'G1_Exsiccatae_Isoduplicate', '', '')).fetchall() if r[0] is not None and str(r[0]).strip()}
                    g1_dna_cats.update(_direct_g1)
                    if _direct_g1 and 'vouchers' in _tbls_g1:
                        _v_cols_g1 = {r[1] for r in _conn_g1.execute('PRAGMA table_info(vouchers)').fetchall()}
                        if {'catalogNumber', 'exsiccatae_key'}.issubset(_v_cols_g1):
                            _ph_g1 = ','.join('?' for _ in _direct_g1)
                            _exs_g1 = {str(r[0]).strip() for r in _conn_g1.execute(f'SELECT DISTINCT exsiccatae_key FROM vouchers WHERE catalogNumber IN ({_ph_g1}) AND TRIM(COALESCE(exsiccatae_key, ?)) != ?', tuple(sorted(_direct_g1)) + ('', '')).fetchall() if r[0] is not None and str(r[0]).strip()}
                            if _exs_g1:
                                _ph2_g1 = ','.join('?' for _ in _exs_g1)
                                _sib_g1 = {str(r[0]).strip() for r in _conn_g1.execute(f'SELECT DISTINCT catalogNumber FROM vouchers WHERE exsiccatae_key IN ({_ph2_g1}) AND TRIM(COALESCE(catalogNumber, ?)) != ?', tuple(sorted(_exs_g1)) + ('', '')).fetchall() if r[0] is not None and str(r[0]).strip()}
                                g1_dna_cats.update(_sib_g1)
        finally:
            _conn_g1.close()
    except Exception:
        pass
if g1_dna_cats and 'exsiccatae_key' in df.columns and 'catalogNumber' in df.columns:
    _df_cats = df['catalogNumber'].astype(str).str.strip()
    _df_exs = df['exsiccatae_key'].fillna('').astype(str).str.strip()
    _g1_exs_local = set(_df_exs[_df_cats.isin(g1_dna_cats)]) - {''}
    if _g1_exs_local:
        g1_dna_cats.update(c for c in _df_cats[_df_exs.isin(_g1_exs_local)] if c)

misidentified_cats = set()
flagged_ramets = set()
triaged_pids_by_cat = {}
misid_pids_by_cat = {}
if exclude_misidentified:
    tokens = [t.strip() for t in (triage_source_raw or 'ALL').split(',') if t.strip()]
    if not tokens:
        tokens = ['ALL']

    self_run_id = resolve_canonical_run_id(run_dir=run_dir, explicit_run_id=os.environ.get('RUN_ID') or None, workspace_root=root_data.parent)
    db_run_count = 0
    if non_demo_pf and Path(target_db_path).is_file():
        try:
            conn_mis = get_db_connection(db_path=target_db_path, init_schema=True, workspace_root=root_data.parent)
            try:
                _db_tbls = {r[0] for r in conn_mis.execute('SELECT name FROM sqlite_master WHERE type=?', ('table',)).fetchall()}
                if 'curation_triage' in _db_tbls:
                    _tri_cols = {r[1] for r in conn_mis.execute('PRAGMA table_info(curation_triage)').fetchall()}
                    if {'run_id', 'catalogNumber', 'taxonomic_treatment_subtype'}.issubset(_tri_cols):
                        _all_rids = [
                            str(r[0]).strip()
                            for r in conn_mis.execute('SELECT DISTINCT run_id FROM curation_triage WHERE TRIM(COALESCE(run_id, ?)) != ?', ('', '')).fetchall()
                            if r[0] is not None and str(r[0]).strip()
                        ]
                        _is_all_tok = any(t.upper() in ('ALL', '*') for t in tokens)
                        _tok_cands = set()
                        if not _is_all_tok:
                            for tok in tokens:
                                _tok_cands.add(tok)
                                _tok_cands.add(Path(tok).name)
                                _tn = resolve_canonical_run_id(explicit_run_id=tok, workspace_root=root_data.parent)
                                if _tn:
                                    _tok_cands.add(_tn)
                                _tp = Path(tok)
                                if _tp.name == 'triage_queue.csv' and len(_tp.parts) >= 4:
                                    _tok_cands.add(_tp.parents[2].name)
                        _matched_rids = []
                        for _rid in _all_rids:
                            if self_run_id and (_rid == self_run_id or _rid.startswith(self_run_id + '/')):
                                continue
                            if _is_all_tok:
                                _matched_rids.append(_rid)
                            elif any(_rid == _tc or _rid.startswith(_tc + '/') or _rid.endswith('/' + _tc) for _tc in _tok_cands if _tc):
                                _matched_rids.append(_rid)
                        db_run_count = len(_matched_rids)
                        if _matched_rids:
                            _pid_expr = 'COALESCE(plant_individual_id, 0)' if 'plant_individual_id' in _tri_cols else '0'
                            _ph_rids = ','.join('?' for _ in _matched_rids)
                            for _cat_r, _pid_r, _sub_r in conn_mis.execute(
                                f'SELECT catalogNumber, {_pid_expr}, taxonomic_treatment_subtype FROM curation_triage WHERE run_id IN ({_ph_rids})',
                                tuple(_matched_rids),
                            ).fetchall():
                                _cs = '' if _cat_r is None else str(_cat_r).strip()
                                if not _cs:
                                    continue
                                _pi = _norm_pid(_pid_r)
                                _ss = '' if _sub_r is None else str(_sub_r).strip()
                                triaged_pids_by_cat.setdefault(_cs, set()).add(_pi)
                                if _ss == 'Misidentification_Reassign' or (_ss != '' and _ss not in ALLOWED_TREATMENT_SUBTYPES):
                                    misid_pids_by_cat.setdefault(_cs, set()).add(_pi)
                                    flagged_ramets.add((_cs, _pi))
            finally:
                conn_mis.close()
        except Exception:
            db_run_count = 0

    runs_candidates = [root_data.parent / 'runs', Path.cwd() / 'runs', run_dir.parent]
    runs_base = next(
        (rc for rc in runs_candidates if rc.is_dir() and any(rc.glob('**/data/tables/triage_queue.csv'))),
        next((rc for rc in runs_candidates if rc.is_dir()), run_dir.parent),
    )
    triage_files = []
    seen_resolved = set()
    self_triage = (run_dir / 'data' / 'tables' / 'triage_queue.csv').resolve()

    for tok in tokens:
        if tok.upper() in ('ALL', '*'):
            for p in sorted(runs_base.glob('**/data/tables/triage_queue.csv')):
                rp = p.resolve()
                if rp != self_triage and rp not in seen_resolved and p.is_file() and p.stat().st_size > 0:
                    seen_resolved.add(rp)
                    triage_files.append(p)
        else:
            tok_norm = resolve_canonical_run_id(explicit_run_id=tok, workspace_root=root_data.parent)
            cand_paths = [Path(tok), Path(tok) / 'data' / 'tables' / 'triage_queue.csv', runs_base / tok / 'data' / 'tables' / 'triage_queue.csv', root_data.parent / 'runs' / tok / 'data' / 'tables' / 'triage_queue.csv', Path.cwd() / 'runs' / tok / 'data' / 'tables' / 'triage_queue.csv']
            if tok_norm and tok_norm != tok:
                cand_paths.extend([runs_base / tok_norm / 'data' / 'tables' / 'triage_queue.csv', root_data.parent / 'runs' / tok_norm / 'data' / 'tables' / 'triage_queue.csv', Path.cwd() / 'runs' / tok_norm / 'data' / 'tables' / 'triage_queue.csv'])
            matched = next((cp for cp in cand_paths if cp.is_file() and cp.stat().st_size > 0), None)
            if matched is None and db_run_count == 0:
                print(
                    f'[ERROR] --triage-source {tok!r} could not be resolved to a SQLite curation_triage run or triage_queue.csv file.',
                    file=sys.stderr,
                )
                sys.exit(1)
            if matched is not None:
                rp = matched.resolve()
                if rp not in seen_resolved:
                    seen_resolved.add(rp)
                    triage_files.append(matched)

    if db_run_count == 0 and not triage_files:
        print(
            f'[ERROR] --exclude-misidentified requested (triage_source={triage_source_raw!r}), but no SQLite curation_triage rows or runs/*/data/tables/triage_queue.csv files were found.',
            file=sys.stderr,
        )
        sys.exit(1)

    for tf in triage_files:
        tdf_tri = pd.read_csv(tf, low_memory=False)
        sub_col = 'taxonomic_treatment_subtype' if 'taxonomic_treatment_subtype' in tdf_tri.columns else None
        if sub_col is None or 'catalogNumber' not in tdf_tri.columns:
            print(
                f'[ERROR] Triage report {tf} is missing required columns (catalogNumber, taxonomic_treatment_subtype).',
                file=sys.stderr,
            )
            sys.exit(1)
        sub_series = tdf_tri[sub_col].fillna('').astype(str).str.strip()
        cat_series = tdf_tri['catalogNumber'].fillna('').astype(str).str.strip()
        pid_series = (
            tdf_tri['plant_individual_id'].map(_norm_pid)
            if 'plant_individual_id' in tdf_tri.columns
            else pd.Series(0, index=tdf_tri.index)
        )
        is_misid = (sub_series == 'Misidentification_Reassign') | (
            ~sub_series.isin(ALLOWED_TREATMENT_SUBTYPES) & (sub_series != '')
        )
        for c_val, p_val, bad_flag in zip(cat_series, pid_series, is_misid):
            if not c_val:
                continue
            triaged_pids_by_cat.setdefault(c_val, set()).add(int(p_val))
            if bool(bad_flag):
                misid_pids_by_cat.setdefault(c_val, set()).add(int(p_val))
                flagged_ramets.add((c_val, int(p_val)))

    # Evict a sheet from curated_vouchers.csv ONLY when 100% of its triaged ramets were misidentified
    for c_val, flagged_pids in misid_pids_by_cat.items():
        if flagged_pids and not (triaged_pids_by_cat.get(c_val, set()) - flagged_pids):
            misidentified_cats.add(c_val)

    before_misid_n = len(df)
    df = df[~df['catalogNumber'].astype(str).str.strip().isin(misidentified_cats)].copy()
    print(
        f'[Preflight: Exclude Misidentified] Source={triage_source_raw!r} (SQLite runs={db_run_count}, CSV reports={len(triage_files)}): '
        f'evicted {before_misid_n - len(df)} sheets where 100% of triaged ramets were flagged '
        f'({len(misidentified_cats)} full sheets, {len(flagged_ramets)} flagged (catalogNumber, plant_individual_id) ramets; '
        f'retaining mixed multi-ramet sheets and filtering only flagged ramets downstream): '
        f'retained {len(df)} / {before_misid_n} sheets.'
    )
__PLANTMATICS_PREFLIGHT_PART_A__
}

_pipe_run_preflight_subset_filter() {
  local _root_data="${1:-}"
  if [[ -z "${_root_data}" ]]; then
    _root_data="${PROJECT_ROOT}/data"
  fi
  echo "[Preflight] Balancing/filtering ${RUN_DIR}/data/tables/ for taxa='${CUSTOM_TAXA_CSV:-}' focal='${FOCAL_TAXON:-}' exclude='${EXCLUDE_TAXON_CSV:-}' determiner_tier='${DETERMINER_TIER:-}' exclude_misidentified='${EXCLUDE_MISIDENTIFIED:-false}' spatial_thin_km='${SPATIAL_THIN_KM:-0}' max_per_state='${MAX_PER_STATE:-}' max_per_region='${MAX_PER_REGION:-}' no_spatial_weights='${NO_SPATIAL_WEIGHTS:-false}'..."
  local _py_src
  _py_src="$(_pipe_preflight_py_part_a)"$'\n'"$(_pipe_preflight_py_part_b)"
  "${PYTHON_MAIN}" -c "${_py_src}" \
    "${RUN_DIR}" "${_root_data}" "${FOCAL_TAXON:-}" "${EXCLUDE_TAXON_CSV:-}" \
    "${INTRASPECIFIC_MODE:-false}" "${GOLD_ANCHORS:-30}" "${DETERMINER_TIER:-}" \
    "${EXCLUDE_MISIDENTIFIED:-false}" "${TRIAGE_SOURCE:-ALL}" "${AUTO_NEIGHBORS:-false}" \
    "${MIN_AFFINITY:-0.03}" "${SPATIAL_THIN_KM:-0}" "${MAX_PER_STATE:-}" "${MAX_PER_REGION:-}" \
    "${NO_SPATIAL_WEIGHTS:-false}" "${REUSE_CANONICAL_DATA:-false}" "${CUSTOM_TAXA_CSV:-}" \
    "${SEED_DATA_DIR:-}" "${START_STEP:-1}" "${USDA_COUNTY_FILTER:-}" \
    "${USDA_COUNTY_BUFFER_KM:-50.0}" "${USDA_EXEMPT_GOLD:-true}" \
    "${PROJECT_ROOT}" "${GEOGRAPHIC_SCOPE:-eastern_na}" "${INCLUDE_WESTERN:-false}"
}

run_preflight_subset_filter() {
  _pipe_run_preflight_subset_filter "$@"
}
