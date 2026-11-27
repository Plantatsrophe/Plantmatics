#!/usr/bin/env bash
# scripts/pipeline/sh/pipe_preflight_cohort_b.sh - Second half of preflight cohort subsetting & SQLite sync.

_pipe_preflight_py_part_b() {
  cat <<'__PLANTMATICS_PREFLIGHT_PART_B__'
if target_tier and 'determiner_tier' in df.columns:
    before_tier_n = len(df)
    tier_vals = df['determiner_tier'].astype(str).str.strip()
    is_g1_df = (
        df['catalogNumber'].astype(str).str.strip().isin(g1_dna_cats)
        if (g1_dna_cats and 'catalogNumber' in df.columns)
        else False
    )
    if target_tier in ('Tier_1_Gold', 'Tier_1_2_Gold_Silver'):
        is_type_col = df['is_type_specimen'].fillna(False).astype(bool) if 'is_type_specimen' in df.columns else False
        allow_t = ['Tier_1_Gold', 'Tier_0_Type'] + (['Tier_2_Silver'] if target_tier == 'Tier_1_2_Gold_Silver' else [])
        df = df[tier_vals.isin(allow_t) | is_type_col | is_g1_df].copy()
    else:
        df = df[tier_vals == target_tier].copy()
    print(f'[Preflight: Determiner Tier] Filtered to determiner_tier == {target_tier!r} (+ types and G1 DNA vouchers when Tier_1_Gold): retained {len(df)} / {before_tier_n} vouchers.')

if exclude_taxa:
    before_n = len(df)
    df = df[~df['species_standardized'].isin(exclude_taxa)].copy()
    retained_taxa = sorted([str(t).strip() for t in df['species_standardized'].dropna().unique().tolist() if str(t).strip()])
    print(f'[Preflight: Exclude Taxa] Excluded {sorted([str(t) for t in exclude_taxa if t])}: retained {len(df)} / {before_n} vouchers across {retained_taxa}.')

default_usda_cache_pf = root_data / 'environmental' / 'usda_plants'
has_local_usda_pf = (
    (default_usda_cache_pf / 'plantlst.txt').is_file()
    or (default_usda_cache_pf / 'distributions').is_dir()
)
syn_map_pf = {}
has_sqlite_usda_pf = False
if non_demo_pf and Path(target_db_path).is_file():
    try:
        conn_syn = get_db_connection(db_path=target_db_path, init_schema=False, read_only=True, workspace_root=root_data.parent)
        try:
            from scripts.core.db import load_taxon_synonyms
            syn_map_pf = load_taxon_synonyms(conn_syn)
            _u_cnt = conn_syn.execute('SELECT COUNT(*) FROM sqlite_master WHERE type=? AND name=?', ('table', 'usda_county_ranges')).fetchone()[0]
            if _u_cnt:
                has_sqlite_usda_pf = bool(conn_syn.execute('SELECT 1 FROM usda_county_ranges LIMIT 1').fetchone())
        finally:
            conn_syn.close()
    except Exception:
        syn_map_pf = {}

if (reuse_canonical or start_step_val >= 2) and (usda_filter_arg == 'true' or has_local_usda_pf or has_sqlite_usda_pf):
    from scripts.core.usda_plants_ranges import filter_vouchers_by_usda_county_range
    before_usda_n = len(df)
    _drop_flag = (usda_filter_arg == 'true')
    df_usda, usda_pf_stats = filter_vouchers_by_usda_county_range(
        df=df,
        standardize_taxon_fn=standardize_packera_taxon,
        synonym_to_canonical=syn_map_pf,
        county_buffer_km=usda_buffer_arg,
        exempt_gold=usda_exempt_gold_arg,
        db_path=target_db_path if non_demo_pf else None,
        workspace_dir=root_data.parent,
        allow_network=(_drop_flag and non_demo_pf),
        drop_out_of_range=_drop_flag,
    )
    df = df_usda
    print(
        f'[Preflight: USDA County Range Audit/Filter (drop_out_of_range={_drop_flag})] '
        f'Retained {len(df)} / {before_usda_n} vouchers '
        f'(buffer={usda_buffer_arg:.1f} km, exempt_gold={usda_exempt_gold_arg}; stats={usda_pf_stats}).'
    )

if non_demo_pf and Path(target_db_path).is_file() and not df.empty:
    try:
        from scripts.core.taxonomy_authority import reconcile_cross_taxon_exsiccatae_with_db
        conn_rec_pf = get_db_connection(db_path=target_db_path, init_schema=False, workspace_root=root_data.parent)
        try:
            if 'exsiccatae_key' in df.columns:
                df = reconcile_cross_taxon_exsiccatae_with_db(df, conn_rec_pf, standardize_fn=standardize_packera_taxon)
            # Persist refreshed USDA range and exsiccatae conflict fields back onto SQLite vouchers
            _v_cols_upd = {r[1] for r in conn_rec_pf.execute('PRAGMA table_info(vouchers)').fetchall()}
            _sync_cols = [c for c in ('usda_symbol', 'usda_county_status', 'exsiccatae_taxonomic_conflict', 'exsiccatae_gold_taxon') if c in _v_cols_upd and c in df.columns]
            if _sync_cols and 'catalogNumber' in df.columns:
                _q2 = chr(39) * 2
                _set_expr = ', '.join(f'{c} = COALESCE(NULLIF(?, {_q2}), {c})' if c != 'exsiccatae_taxonomic_conflict' else f'{c} = MAX(COALESCE({c}, 0), ?)' for c in _sync_cols)
                _upd_rows = []
                for _, _ur in df.iterrows():
                    _cat_u = str(_ur.get('catalogNumber') or '').strip()
                    if not _cat_u:
                        continue
                    _vals_u = []
                    for _sc in _sync_cols:
                        _v = _ur.get(_sc)
                        if _sc == 'exsiccatae_taxonomic_conflict':
                            _vals_u.append(1 if (_v is True or str(_v).strip().lower() in ('true', '1', 't', 'yes')) else 0)
                        else:
                            _sv = '' if (_v is None or (isinstance(_v, float) and math.isnan(_v))) else str(_v).strip()
                            _vals_u.append('' if _sv in ('None', 'nan') else _sv)
                    _vals_u.append(_cat_u)
                    _upd_rows.append(tuple(_vals_u))
                if _upd_rows:
                    conn_rec_pf.executemany(
                        f'UPDATE vouchers SET {_set_expr} WHERE catalogNumber = ?',
                        _upd_rows,
                    )
                    conn_rec_pf.commit()
        finally:
            conn_rec_pf.close()
    except Exception:
        pass

if focal_raw:
    focal_taxon = standardize_packera_taxon(focal_raw)
    focal_df = df[df['species_standardized'] == focal_taxon].copy()
    if intraspecific_mode:
        subset_df = focal_df
        print(f'[Strategy 2: Pure Intraspecific] Retained {len(subset_df)} {focal_taxon} vouchers (0 outgroup anchors).')
    else:
        present_taxa = sorted([t for t in df['species_standardized'].dropna().unique().tolist() if str(t).startswith('Packera ')])
        base_outgroups = present_taxa
        outgroup_taxa = [t for t in base_outgroups if t != focal_taxon and t not in exclude_taxa]

        if auto_neighbors_flag:
            discovered_neighbors = []
            runs_search_root = next(
                (
                    cand for cand in (root_data.parent / 'runs', Path.cwd() / 'runs')
                    if cand.is_dir() and cand in run_dir.parents
                ),
                None,
            )
            if runs_search_root is None:
                runs_search_root = next(
                    (p for p in run_dir.parents if p.name == 'runs'),
                    run_dir.parent,
                )
            ts_reports = None
            if triage_source_raw and triage_source_raw.upper() not in ('ALL', '*'):
                ts_norm = resolve_canonical_run_id(explicit_run_id=triage_source_raw, workspace_root=root_data.parent)
                ts_cands = [
                    runs_search_root / ts_norm / 'outputs' / 'reports',
                    Path(triage_source_raw) / 'outputs' / 'reports',
                    root_data.parent / 'runs' / ts_norm / 'outputs' / 'reports',
                ]
                ts_reports = next((p for p in ts_cands if p.is_dir()), None)
            if ts_reports is not None:
                neigh_candidates = [
                    run_dir / 'outputs' / 'reports' / 'focal_auto_neighbors.txt',
                    ts_reports / 'focal_auto_neighbors.txt',
                ]
                comp_candidates = [
                    run_dir / 'outputs' / 'reports' / 'empirical_taxon_complexes.csv',
                    ts_reports / 'empirical_taxon_complexes.csv',
                ]
                aff_candidates = [
                    run_dir / 'outputs' / 'reports' / 'pairwise_taxon_affinity_matrix.csv',
                    ts_reports / 'pairwise_taxon_affinity_matrix.csv',
                ]
            else:
                neigh_candidates = [
                    run_dir / 'outputs' / 'reports' / 'focal_auto_neighbors.txt',
                    runs_search_root / 'latest' / 'outputs' / 'reports' / 'focal_auto_neighbors.txt',
                ] + sorted(runs_search_root.glob('**/outputs/reports/focal_auto_neighbors.txt'), reverse=True)
                comp_candidates = [
                    run_dir / 'outputs' / 'reports' / 'empirical_taxon_complexes.csv',
                    runs_search_root / 'latest' / 'outputs' / 'reports' / 'empirical_taxon_complexes.csv',
                ] + sorted(runs_search_root.glob('**/outputs/reports/empirical_taxon_complexes.csv'), reverse=True)
                aff_candidates = [
                    run_dir / 'outputs' / 'reports' / 'pairwise_taxon_affinity_matrix.csv',
                    runs_search_root / 'latest' / 'outputs' / 'reports' / 'pairwise_taxon_affinity_matrix.csv',
                ] + sorted(runs_search_root.glob('**/outputs/reports/pairwise_taxon_affinity_matrix.csv'), reverse=True)

            comp_file = next((p for p in comp_candidates if p.is_file() and p.stat().st_size > 0), None)
            aff_file = next((p for p in aff_candidates if p.is_file() and p.stat().st_size > 0), None)
            if comp_file is not None:
                cdf = pd.read_csv(comp_file, low_memory=False)
                sp_c_col = 'species_standardized' if 'species_standardized' in cdf.columns else ('taxon' if 'taxon' in cdf.columns else None)
                if sp_c_col is not None and 'consensus_complex_id' in cdf.columns:
                    cdf['taxon_std'] = cdf[sp_c_col].apply(standardize_packera_taxon)
                    focal_rows = cdf[cdf['taxon_std'] == focal_taxon]
                    if not focal_rows.empty:
                        f_cid = focal_rows['consensus_complex_id'].iloc[0]
                        same_c = [
                            col_std for col_std in cdf[(cdf['consensus_complex_id'] == f_cid) & (cdf['taxon_std'] != focal_taxon)]['taxon_std'].tolist()
                            if isinstance(col_std, str) and col_std.startswith('Packera ') and col_std != focal_taxon
                        ]
                        discovered_neighbors.extend(same_c)
            if aff_file is not None:
                adf = pd.read_csv(aff_file, low_memory=False)
                sp_a_col = 'species_standardized' if 'species_standardized' in adf.columns else ('taxon' if 'taxon' in adf.columns else None)
                if sp_a_col is not None:
                    adf['taxon_std'] = adf[sp_a_col].apply(standardize_packera_taxon)
                    frow = adf[adf['taxon_std'] == focal_taxon]
                    if not frow.empty:
                        for col in adf.columns:
                            if col in ('species_standardized', 'taxon', 'taxon_std'):
                                continue
                            col_std = standardize_packera_taxon(col)
                            if isinstance(col_std, str) and col_std.startswith('Packera ') and col_std != focal_taxon:
                                try:
                                    val = float(frow[col].iloc[0])
                                    if val >= min_affinity_val:
                                        discovered_neighbors.append(col_std)
                                except Exception:
                                    pass
            if not discovered_neighbors:
                nf = next((p for p in neigh_candidates if p.is_file() and p.stat().st_size > 0), None)
                if nf is not None:
                    discovered_neighbors = [
                        col_std
                        for line in nf.read_text(encoding='utf-8').splitlines()
                        if line.strip()
                        for col_std in (standardize_packera_taxon(line.strip()),)
                        if isinstance(col_std, str) and col_std.startswith('Packera ') and col_std != focal_taxon
                    ]
            if discovered_neighbors:
                disc_set = {t for t in discovered_neighbors if isinstance(t, str) and t.startswith('Packera ')} - exclude_taxa - {focal_taxon}
                if disc_set:
                    outgroup_taxa = [t for t in outgroup_taxa if t in disc_set] or sorted(disc_set)
                    print(f'[Auto-Neighbors] Empirical complex & affinity neighbors for {focal_taxon}: {outgroup_taxa}')

        outgroup_dfs, outgroup_counts, outgroup_tiers = select_cascading_outgroup_anchors(
            df,
            outgroup_taxa,
            max_anchors=gold_anchors,
            workspace_dir=root_data.parent,
            focal_df=focal_df,
        )
        subset_df = pd.concat([focal_df] + outgroup_dfs, ignore_index=True) if outgroup_dfs else focal_df
        print(
            f'[Strategy 1: Focal + Outgroup Cascading Anchors] Focal {focal_taxon}: {len(focal_df)} | '
            f'Outgroup anchors (Gold->Silver->Bronze): {outgroup_counts} | Tier breakdown: {outgroup_tiers} '
            f'(Total: {len(subset_df)})'
        )
else:
    subset_df = df

if target_tier and 'determiner_tier' in subset_df.columns:
    sub_tiers = subset_df['determiner_tier'].astype(str).str.strip()
    is_g1_sub = (
        subset_df['catalogNumber'].astype(str).str.strip().isin(g1_dna_cats)
        if (g1_dna_cats and 'catalogNumber' in subset_df.columns)
        else False
    )
    if target_tier in ('Tier_1_Gold', 'Tier_1_2_Gold_Silver'):
        sub_types = subset_df['is_type_specimen'].fillna(False).astype(bool) if 'is_type_specimen' in subset_df.columns else False
        allow_sub = ['Tier_1_Gold', 'Tier_0_Type'] + (['Tier_2_Silver'] if target_tier == 'Tier_1_2_Gold_Silver' else [])
        subset_df = subset_df[sub_tiers.isin(allow_sub) | sub_types | is_g1_sub].copy()
    else:
        subset_df = subset_df[sub_tiers == target_tier].copy()
if exclude_misidentified and misidentified_cats:
    subset_df = subset_df[~subset_df['catalogNumber'].astype(str).str.strip().isin(misidentified_cats)].copy()

manifest_path = tables_dir / 'extracted_leaf_manifest.csv'
manifest_df = pd.read_csv(manifest_path, low_memory=False) if manifest_path.exists() and manifest_path.stat().st_size > 0 else None
if manifest_df is not None and flagged_ramets and {'catalogNumber', 'plant_individual_id'}.issubset(manifest_df.columns):
    _mc = manifest_df['catalogNumber'].fillna('').astype(str).str.strip()
    _mp = manifest_df['plant_individual_id'].map(_norm_pid)
    _keep_m = [((c_m, p_m) not in flagged_ramets) for c_m, p_m in zip(_mc, _mp)]
    manifest_df = manifest_df.loc[_keep_m].copy()

pre_balance_n = len(subset_df)
if spatial_thin_km > 0.0 and 'catalogNumber' in subset_df.columns:
    _hab_map = {}
    for _tf_h in (triage_files if exclude_misidentified else []):
        for _fn_h in ('multimodal_conflict_flags.csv', 'ssurgo_edaphic_validation.csv'):
            _cand_h = _tf_h.parent / _fn_h
            if _cand_h.is_file() and _cand_h.stat().st_size > 0:
                try:
                    _hdf = pd.read_csv(_cand_h, low_memory=False)
                    if 'catalogNumber' in _hdf.columns:
                        for _, _hr in _hdf.iterrows():
                            _hc = str(_hr.get('catalogNumber') or '').strip()
                            if not _hc or _hc in _hab_map:
                                continue
                            _er = str(_hr.get('ecotype_regime') or '').strip()
                            if _er.startswith('Ecotype_'):
                                _hab_map[_hc] = _er
                            elif str(_hr.get('is_rock_outcrop') or '').strip().lower() in ('true', '1', 't', 'yes'):
                                _hab_map[_hc] = 'azonal_outcrop'
                            elif 'lithic_outcrop_index' in _hdf.columns and pd.notna(_hr.get('lithic_outcrop_index')):
                                _hab_map[_hc] = 'azonal_outcrop' if float(_hr.get('lithic_outcrop_index')) >= 0.45 else 'zonal_matrix'
                except Exception:
                    pass
    if _hab_map:
        subset_df['_thin_substrate_stratum'] = subset_df['catalogNumber'].astype(str).str.strip().map(lambda c: _hab_map.get(c, ''))
subset_df = spatially_balance_and_thin_vouchers(
    subset_df,
    spatial_thin_km=spatial_thin_km,
    max_per_state=max_per_state,
    max_per_region=max_per_region,
    primary_duplicates_only=True,
    preserve_all_gold=(not (target_tier in ('Tier_1_Gold', 'Tier_1_2_Gold_Silver') and spatial_thin_km > 0.0)),
    preserve_catalog_numbers=sorted(g1_dna_cats) if g1_dna_cats else None,
    manifest_df=manifest_df,
    workspace_dir=root_data.parent,
)
if '_thin_substrate_stratum' in subset_df.columns:
    subset_df = subset_df.drop(columns=['_thin_substrate_stratum'])
if no_spatial_weights:
    subset_df['spatial_balance_weight'] = 1.0

has_spatial_thin_or_quota = (
    spatial_thin_km > 0.0
    or (max_per_state is not None and max_per_state > 0)
    or (max_per_region is not None and max_per_region > 0)
)
if has_spatial_thin_or_quota or len(subset_df) < pre_balance_n:
    print(
        f'[Preflight: Spatial Balancing & Deduplication] Retained {len(subset_df)} / {pre_balance_n} vouchers '
        f'(spatial_thin_km={spatial_thin_km}, max_per_state={max_per_state}, max_per_region={max_per_region}).'
    )
if has_spatial_thin_or_quota and not subset_df.empty:
    if 'stateProvince' in subset_df.columns:
        st_summary = subset_df['stateProvince'].fillna('UNKNOWN').astype(str).value_counts().to_dict()
        print(f'  Retained counts per stateProvince: {st_summary}')
    if 'regional_group' in subset_df.columns:
        reg_summary = subset_df['regional_group'].fillna('NA').astype(str).value_counts().to_dict()
        print(f'  Retained counts per regional_group: {reg_summary}')

def _write_atomic(df_obj, target_path):
    tmp_p = target_path.with_suffix(target_path.suffix + '.tmp')
    with open(tmp_p, 'w', encoding='utf-8', newline='') as f:
        df_obj.to_csv(f, index=False)
        f.flush()
        os.fsync(f.fileno())
    tmp_p.replace(target_path)

_write_atomic(subset_df, csv_path)
keep_cats = set(subset_df['catalogNumber'].astype(str).str.strip())

has_subset_filter = bool(
    custom_taxa_csv or focal_raw or exclude_csv or target_tier or exclude_misidentified or has_spatial_thin_or_quota or len(subset_df) < initial_voucher_count
)

if has_subset_filter:
    # Filter Step 02 and Stage 3 run-scoped tables in ${RUN_DIR}/data/tables/ by retained catalogNumbers
    # and strictly prune flagged (catalogNumber, plant_individual_id) ramets when present.
    _core_tbls = [
        'extracted_leaf_manifest.csv',
        'failed_qc_vouchers.csv',
        'voucher_reproductive_metrics.csv',
        'voucher_phenological_states.csv',
    ]
    _extra_tbls = sorted(
        p.name for p in tables_dir.glob('*.csv')
        if p.name not in _core_tbls and p.name != 'curated_vouchers.csv'
    )
    for tbl_name in (_core_tbls + _extra_tbls):
        tp = tables_dir / tbl_name
        if tp.exists() and tp.stat().st_size > 0:
            tdf = pd.read_csv(tp, low_memory=False)
            if 'catalogNumber' in tdf.columns:
                _t_cats = tdf['catalogNumber'].fillna('').astype(str).str.strip()
                keep_row_mask = _t_cats.isin(keep_cats)
                if flagged_ramets and 'plant_individual_id' in tdf.columns:
                    _t_pids = tdf['plant_individual_id'].map(_norm_pid)
                    _ramet_keep = pd.Series(
                        [((c_r, p_r) not in flagged_ramets) for c_r, p_r in zip(_t_cats, _t_pids)],
                        index=tdf.index,
                    )
                    keep_row_mask = keep_row_mask & _ramet_keep
                tdf_sub = tdf.loc[keep_row_mask].copy()
                _write_atomic(tdf_sub, tp)
                print(f'  Filtered {tbl_name}: {len(tdf_sub)} / {len(tdf)} rows')

# Register pipeline_runs and populate run_active_cohort in SQLite (skipping canonical DB under --demo)
if os.environ.get('PLANTMATICS_DEMO_MODE') != '1' and os.environ.get('PACKERA_DEMO_MODE') != '1':
    run_id_val = resolve_canonical_run_id(
        run_dir=run_dir,
        explicit_run_id=os.environ.get('RUN_ID') or None,
        workspace_root=root_data.parent,
    )
    conn_coh = get_db_connection(db_path=target_db_path, init_schema=True, workspace_root=root_data.parent)
    try:
        ensure_pipeline_run(
            conn_coh,
            run_id=run_id_val,
            run_dir=run_dir,
            focal_taxon=standardize_packera_taxon(focal_raw) if focal_raw else None,
            intraspecific_mode=intraspecific_mode,
        )
        n_cohort = populate_run_active_cohort(
            conn_coh,
            run_id=run_id_val,
            voucher_rows=subset_df.to_dict(orient='records'),
            clear_existing=True,
        )
        if exclude_misidentified and flagged_ramets:
            _del_params = [(run_id_val, c_d, int(p_d)) for (c_d, p_d) in sorted(flagged_ramets)]
            _coh_tbls = {r[0] for r in conn_coh.execute('SELECT name FROM sqlite_master WHERE type=?', ('table',)).fetchall()}
            for _tier_b_tbl in ('run_active_cohort', 'leaf_efa_harmonics', 'otu_morphometrics', 'spatial_rf_predictions', 'curation_triage'):
                if _tier_b_tbl in _coh_tbls:
                    _tb_cols = {r[1] for r in conn_coh.execute(f'PRAGMA table_info({_tier_b_tbl})').fetchall()}
                    if {'run_id', 'catalogNumber', 'plant_individual_id'}.issubset(_tb_cols):
                        conn_coh.executemany(
                            f'DELETE FROM {_tier_b_tbl} WHERE run_id = ? AND catalogNumber = ? AND COALESCE(plant_individual_id, 0) = ?',
                            _del_params,
                        )
            conn_coh.commit()
            n_cohort = conn_coh.execute(
                'SELECT COUNT(*) FROM run_active_cohort WHERE run_id = ?',
                (run_id_val,),
            ).fetchone()[0]
        print(f'[Preflight: SQLite Cohort] Registered run_id={run_id_val!r} with {n_cohort} active (catalogNumber, plant_individual_id) rows.')
    finally:
        conn_coh.close()

# Use directory-level symlinks to upstream/canonical cropped_patches/capitula
# whenever reusing Step 01-02 assets or starting at Step >= 3 (preserving non-empty run-local patch directories)
if reuse_canonical or start_step_val >= 3:
    run_patches_root = run_dir / 'data' / 'cropped_patches'
    if run_patches_root.is_symlink():
        run_patches_root.unlink()
    run_patches_root.mkdir(parents=True, exist_ok=True)
    for sub in ('capitula',):
        seed_sub = seed_data / 'cropped_patches' / sub
        root_sub = root_data / 'cropped_patches' / sub
        src_sub = seed_sub if (seed_sub.is_dir() and any(seed_sub.iterdir())) else root_sub
        if not src_sub.exists():
            continue
        dst_sub = run_patches_root / sub
        if dst_sub.resolve() == src_sub.resolve() and not dst_sub.is_symlink():
            continue
        if dst_sub.is_symlink() or dst_sub.is_file():
            dst_sub.unlink()
            os.symlink(src_sub, dst_sub)
        elif not dst_sub.exists():
            os.symlink(src_sub, dst_sub)
        elif dst_sub.is_dir() and next(dst_sub.iterdir(), None) is None:
            dst_sub.rmdir()
            os.symlink(src_sub, dst_sub)
__PLANTMATICS_PREFLIGHT_PART_B__
}
