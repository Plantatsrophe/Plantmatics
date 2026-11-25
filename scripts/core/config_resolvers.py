"""Botanical ontology constants, specialist patterns, defaults, and runtime compute/genomics resolvers."""

from __future__ import annotations

import os
from pathlib import Path
import sys
from typing import Any, Dict, List, Optional, Set, Tuple

import yaml

from scripts.core.config_pipeline_cls import PipelineConfig
from scripts.core.config_run_context import DEFAULT_CONFIG_PATH

# Botanical ontology mappings for Packera phenotyping
CLASS_NAMES: List[str] = [
    "basal_leaf_blade", "leaf_petiole", "cauline_leaf",
    "cauline_stem", "root_rhizome", "basal_rosette_clump", "capitulum",
]
CLASS_MAP: Dict[str, int] = {name: idx for idx, name in enumerate(CLASS_NAMES)}
CLASS_COLORS_BGR: Dict[int, Tuple[int, int, int]] = {
    0: (0, 200, 0), 1: (50, 255, 150), 2: (0, 165, 255),
    3: (0, 215, 255), 4: (50, 100, 200), 5: (0, 140, 70), 6: (230, 30, 230),
}

# Layer 1A: Core post-1960 Packera / Senecioneae monographers and segregate authors (Tier_1_Gold)
CORE_PACKERA_MONOGRAPHER_PATTERNS: List[str] = [
    r"\bBarkley\b", r"\bT\.?\s*M\.?\s*Barkley\b", r"\bTheodore\s+M\.?\s+Barkley\b",
    r"\bTrock\b", r"\bD\.?\s*K\.?\s*Trock\b", r"\bDebra\s+K\.?\s+Trock\b",
    r"\bKowal\b", r"\bR\.?\s*R\.?\s*Kowal\b", r"\bRobert\s+R\.?\s+Kowal\b",
    r"\bWeakley\b", r"\bA\.?\s*S\.?\s*Weakley\b", r"\bAlan\s+S\.?\s+Weakley\b",
    r"\bBain\b", r"\bJ\.?\s*F\.?\s*Bain\b", r"\bJohn\s+F\.?\s+Bain\b",
    r"\bMahoney\b", r"\bA\.?\s*M\.?\s*Mahoney\b", r"\bAlison\s+M\.?\s+Mahoney\b",
    r"\bSoderholm\b", r"\bJ\.?\s*Soderholm\b",
    r"\bJ\.?\s*G\.?\s*Packer\b|\bJohn\s+G\.?\s+Packer\b|\bPacker\s*,\s*(?:J\.?\s*G\.?|John)\b",
    r"\bMandel\b", r"\bJ\.?\s*R\.?\s*Mandel\b", r"\bJennifer\s+R\.?\s+Mandel\b",
    r"\bMoore[-\s]?Pollard\b",
    r"\b(?:E\.?\s*R\.?|Erika\s+R\.?)\s+Moore\b|\bMoore\s*,\s*(?:E\.?\s*R\.?|Erika)\b",
    r"\bBoufford\b", r"\bD\.?\s*E\.?\s*Boufford\b",
    r"\bGramling\b", r"\bA\.?\s*Gramling\b",
    r"\bC\.?\s*C\.?\s*Freeman\b|\bCraig\s+(?:C\.?\s+)?Freeman\b|\bFreeman\s*,\s*(?:C\.?\s*C\.?|Craig)\b",
    r"\bW\.?\s*A\.?\s*Weber\b|\bWilliam\s+A\.?\s+Weber\b|\bWeber\s*,\s*(?:W\.?\s*A\.?|William)\b",
    r"\bUttal\b", r"\bL\.?\s*J\.?\s*Uttal\b", r"\bLeonard\s+J\.?\s+Uttal\b",
    r"\bNesom\b", r"\bG\.?\s*L\.?\s*Nesom\b", r"\bGuy\s+L\.?\s+Nesom\b",
    r"\bJ\.?\s*L\.?\s*Golden\b|\bJudith\s+L\.?\s+Golden\b|\bGolden\s*,\s*(?:J\.?\s*L\.?|Judith)\b",
    r"\bL\.?\s*D\.?\s*Estes\b|\bDwayne\s+Estes\b|\bEstes\s*,\s*(?:L\.?\s*D\.?|Dwayne)\b",
    r"\bPoindexter\b", r"\bD\.?\s*B\.?\s*Poindexter\b",
]

# Layer 1B: Regional flora authorities and historical Senecioneae authors (Tier_2_Silver)
REGIONAL_AND_HISTORICAL_SPECIALIST_PATTERNS: List[str] = [
    r"\bC\.?\s*Jeffrey\b|\bCharles\s+Jeffrey\b|\bJeffrey\s*,\s*(?:C\.?|Charles)\b",
    r"\bMabberley\b", r"\bD\.?\s*J\.?\s*Mabberley\b",
    r"\bCronquist\b", r"\bA\.?\s*Cronquist\b", r"\bArthur\s+(?:J\.?\s+)?Cronquist\b",
    r"\bGreenman\b", r"\bJ\.?\s*M\.?\s*Greenman\b", r"\bJesse\s+M(?:ore)?\.?\s+Greenman\b",
    r"(?:^|[\s|,;(])[ÁA]\.?\s*L[öo]ve\b", r"\bD\.?\s*L[öo]ve\b",
    r"\b[ÁA]skell\s+L[öo]ve\b", r"\bDoris\s+L[öo]ve\b",
    r"\bSpellenberg\b", r"\bR\.?\s*(?:W\.?\s*)?Spellenberg\b", r"\bRichard\s+Spellenberg\b",
    r"\bPippen\b", r"\bR\.?\s*W\.?\s*Pippen\b", r"\bRichard\s+W\.?\s+Pippen\b",
    r"\bJudziewicz\b", r"\bE\.?\s*J\.?\s*Judziewicz\b",
    r"\bO['’]?Kennon\b", r"\bR\.?\s*J\.?\s*O['’]?Kennon\b",
    r"\bS\.?\s*L\.?\s*Welsh\b|\bStanley\s+L\.?\s+Welsh\b|\bWelsh\s*,\s*(?:S\.?\s*L\.?|Stanley)\b",
    r"\bN\.?\s*D\.?\s*Atwood\b|\bAtwood\s*,\s*(?:N\.?\s*D\.?|Duane)\b",
    r"\bR\.?\s*D\.?\s*Dorn\b|\bRobert\s+D\.?\s+Dorn\b|\bDorn\s*,\s*(?:R\.?\s*D\.?|Robert)\b",
    r"\bO['’]?Kane\b", r"\bS\.?\s*L\.?\s*O['’]?Kane\b",
    r"\bYeatts\b", r"\bL\.?\s*Yeatts\b",
    r"\bRuyle[-\s]?Douglas\b",
    r"\b(?:G\.?\s*W\.?|George\s+W\.?)\s+Douglas\b|\bDouglas\s*,\s*(?:G\.?\s*W\.?|George)\b",
    r"\bBoivin\b", r"\bB\.?\s*Boivin\b",
    r"\bFernald\b", r"\bM\.?\s*L\.?\s*Fernald\b", r"\bMerritt\s+L(?:yndon)?\.?\s+Fernald\b",
    r"\bRydberg\b", r"\bP\.?\s*A\.?\s*Rydberg\b", r"\bPer\s+Axel\s+Rydberg\b",
    r"\bN\.?\s*L\.?\s*Britton\b|\bNathaniel\s+L(?:ord)?\.?\s+Britton\b|\bBritton\s*,\s*(?:N\.?\s*L\.?|Nathaniel)\b",
    r"\b(?:E\.?\s*L\.?|Edward\s+Lee)\s+Greene\b|\bGreene\s*,\s*(?:E\.?\s*L\.?|Edward)\b",
    r"\b(?:A\.?|Asa)\s+Gray\b|\bGray\s*,\s*(?:A\.?|Asa)\b",
    r"\b(?:J\.?\s*K\.?|John\s+Kunkel)\s+Small\b|\bSmall\s*,\s*(?:J\.?\s*K\.?|John)\b",
    r"\b(?:Alph\.?|Alphonso)\s+Wood\b|\bWood\s*,\s*(?:Alph\.?|Alphonso)\b",
    r"\b(?:S\.?\s*F\.?|Sidney\s+Fay)\s+Blake\b|\bBlake\s*,\s*(?:S\.?\s*F\.?|sidney)\b",
    r"\b(?:S\.?|Sereno)\s+Watson\b|\bWatson\s*,\s*(?:S\.?|Sereno)\b",
    r"\bA\.?\s*A\.?\s*Heller\b|\bAmos\s+A(?:rthur)?\.?\s+Heller\b|\bHeller\s*,\s*(?:A\.?\s*A\.?|Amos)\b",
    r"\bKral\b", r"\bR\.?\s*Kral\b", r"\bRobert\s+Kral\b",
    r"\bGodfrey\b", r"\bR\.?\s*K\.?\s*Godfrey\b", r"\bRobert\s+K\.?\s+Godfrey\b",
    r"\bSorrie\b", r"\bB\.?\s*A\.?\s*Sorrie\b", r"\bBruce\s+A\.?\s+Sorrie\b",
    r"\bStrother\b", r"\bJ\.?\s*L\.?\s*Strother\b", r"\bJohn\s+L\.?\s+Strother\b",
    r"\bB\.?\s*L\.?\s*Turner\b|\bBillie\s+L(?:ee)?\.?\s+Turner\b|\bTurner\s*,\s*(?:B\.?\s*L\.?|Billie)\b",
    r"\bBrouillet\b", r"\bL\.?\s*Brouillet\b",
    r"\bSemple\b", r"\bJ\.?\s*C\.?\s*Semple\b",
    r"\bSchilling\b", r"\bE\.?\s*E\.?\s*Schilling\b",
    r"\bFloden\b", r"\bA\.?\s*J\.?\s*Floden\b",
    r"\bWieboldt\b", r"\bT\.?\s*F\.?\s*Wieboldt\b",
    r"\bMcMillan\b", r"\bP\.?\s*D\.?\s*McMillan\b",
    r"\bBeauchamp\b", r"\bR\.?\s*M\.?\s*Beauchamp\b",
    r"\b(?:J\.?\s*J\.?\s*N\.?|Julian(?:\s+J\.?\s*N\.?)?)\s+Campbell\b|\bCampbell\s*,\s*(?:J\.?\s*J\.?\s*N\.?|Julian)\b",
]

SPECIALIST_PATTERNS: List[str] = (
    CORE_PACKERA_MONOGRAPHER_PATTERNS + REGIONAL_AND_HISTORICAL_SPECIALIST_PATTERNS
)
MAJOR_HERBARIA_CODES: Set[str] = {
    "NCU", "DUKE", "NCSC", "GA", "US", "NY", "BRIT", "MO", "WIS", "VDB", "FLAS", "TEX", "LL", "TENN", "F",
    "GH", "PH", "KANU", "MSC", "COLO", "RM", "UC", "JEPS", "CAS", "DAO", "CAN", "MEXU",
    "USCH", "CLEMS", "VPI", "WILLI", "WVA", "CM", "MICH", "MIN", "OS", "ILL", "ILLS", "IND",
    "NEBC", "VT", "TRT", "MT", "QFA", "ACAD", "UNB", "UARK", "ARK", "LSU", "MISS", "UNA",
    "APSC", "EKY", "KY", "MU", "PAC", "PENN", "CONN", "NYS", "BH",
}
PRIMARY_TYPE_STATUSES: Set[str] = {
    "HOLOTYPE", "LECTOTYPE", "NEOTYPE", "EPITYPE", "ISOTYPE", "ISOLECTOTYPE", "ISONEOTYPE"
}
SECONDARY_TYPE_STATUSES: Set[str] = {
    "SYNTYPE", "ISOSYNTYPE", "PARATYPE", "ISOPARATYPE", "TYPE", "COTYPE"
}
VALID_TYPE_STATUSES: Set[str] = PRIMARY_TYPE_STATUSES | SECONDARY_TYPE_STATUSES
EXCLUDED_WESTERN_STATES: Set[str] = {
    "CO", "COLORADO", "NM", "NEW MEXICO", "WY", "WYOMING", "MT", "MONTANA",
    "UT", "UTAH", "AZ", "ARIZONA", "NV", "NEVADA", "ID", "IDAHO",
    "WA", "WASHINGTON", "OR", "OREGON", "CA", "CALIFORNIA", "AK", "ALASKA", "HI", "HAWAII",
    "BC", "BRITISH COLUMBIA", "AB", "ALBERTA", "SK", "SASKATCHEWAN",
    "YT", "YUKON", "YUKON TERRITORY", "NT", "NORTHWEST TERRITORIES", "NU", "NUNAVUT",
}

# Strictly cis-Mississippi US states + D.C. and Eastern Canadian provinces (excluding trans-Mississippi Plains/Highlands)
CIS_MISSISSIPPI_STATES_PROVINCES: Set[str] = {
    "AL", "ALABAMA", "CT", "CONNECTICUT", "DC", "DISTRICT OF COLUMBIA",
    "WASHINGTON, D.C.", "WASHINGTON D.C.", "WASHINGTON DC",
    "DE", "DELAWARE", "FL", "FLORIDA", "GA", "GEORGIA", "IL", "ILLINOIS",
    "IN", "INDIANA", "KY", "KENTUCKY", "MA", "MASSACHUSETTS", "MD", "MARYLAND",
    "ME", "MAINE", "MI", "MICHIGAN", "MS", "MISSISSIPPI", "NC", "NORTH CAROLINA",
    "NH", "NEW HAMPSHIRE", "NJ", "NEW JERSEY", "NY", "NEW YORK", "OH", "OHIO",
    "PA", "PENNSYLVANIA", "RI", "RHODE ISLAND", "SC", "SOUTH CAROLINA",
    "TN", "TENNESSEE", "VA", "VIRGINIA", "VT", "VERMONT", "WI", "WISCONSIN",
    "WV", "WEST VIRGINIA",
    "ON", "ONTARIO", "QC", "QUEBEC", "QUÉBEC", "NB", "NEW BRUNSWICK",
    "NS", "NOVA SCOTIA", "PE", "PEI", "PRINCE EDWARD ISLAND",
    "NL", "NEWFOUNDLAND", "NEWFOUNDLAND AND LABRADOR",
}

# Great Plains states/provinces west of the immediate Mississippi-border tier (MN, IA, MO, AR, LA)
GREAT_PLAINS_STRADDLING_STATES: Set[str] = {
    "ND", "NORTH DAKOTA", "SD", "SOUTH DAKOTA", "NE", "NEBRASKA",
    "KS", "KANSAS", "OK", "OKLAHOMA", "TX", "TEXAS", "MB", "MANITOBA",
}

# Module-level defaults instantiated dynamically from configuration YAML
_DEFAULT_CONFIG = PipelineConfig.from_yaml()
CANONICAL_CIS_MISSISSIPPI_PACKERA_SPECIES: Set[str] = set(_DEFAULT_CONFIG.taxa.eastern_resident_whitelist)
TRANSCONTINENTAL_BOREAL_PACKERA_SPECIES: Set[str] = set(_DEFAULT_CONFIG.taxa.transcontinental_boreal_taxa)
DEFAULT_WORKSPACE = _DEFAULT_CONFIG.paths.workspace_root
DEFAULT_DB_PATH = _DEFAULT_CONFIG.paths.db_path
DEFAULT_RAW_DIR = _DEFAULT_CONFIG.paths.raw_vouchers_dir
DEFAULT_CURATED_CSV = _DEFAULT_CONFIG.paths.curated_vouchers_csv
DEFAULT_OUTPUT_CSV = _DEFAULT_CONFIG.paths.curated_vouchers_csv
DEFAULT_SUMMARY_LOG = _DEFAULT_CONFIG.paths.summary_log
DEFAULT_QUARANTINE_DIR = _DEFAULT_CONFIG.paths.quarantine_dir
DEFAULT_MIN_MEGAPIXELS = _DEFAULT_CONFIG.thresholds.min_megapixels
DEFAULT_MIN_FILE_SIZE_KB = _DEFAULT_CONFIG.thresholds.min_file_size_kb
DEFAULT_MIN_SHARPNESS_LAPLACIAN = _DEFAULT_CONFIG.thresholds.min_sharpness_laplacian
DEFAULT_TARGET_TAXA: List[str] = []
WESTERN_LONGITUDE_THRESHOLD = _DEFAULT_CONFIG.thresholds.western_longitude_threshold
DEFAULT_CIS_MISSISSIPPI_LON = _DEFAULT_CONFIG.harvesting.cis_mississippi_longitude_threshold
DEFAULT_MIN_CIS_MISSISSIPPI_FRAC = _DEFAULT_CONFIG.harvesting.min_cis_mississippi_fraction
DEFAULT_TRANSCONTINENTAL_LON = _DEFAULT_CONFIG.harvesting.transcontinental_longitude_threshold
DEFAULT_GEOGRAPHIC_SCOPE = _DEFAULT_CONFIG.harvesting.geographic_scope
DEFAULT_COUNTRIES = _DEFAULT_CONFIG.harvesting.countries
DEFAULT_HOME_HERBARIUM = _DEFAULT_CONFIG.harvesting.home_herbarium
DEFAULT_USDA_COUNTY_FILTER = _DEFAULT_CONFIG.harvesting.usda_county_filter
DEFAULT_USDA_COUNTY_BUFFER_KM = _DEFAULT_CONFIG.harvesting.usda_county_buffer_km
DEFAULT_USDA_EXEMPT_GOLD = _DEFAULT_CONFIG.harvesting.usda_exempt_gold
DEFAULT_INCLUDE_UNGEOREFERENCED_STATE_VOUCHERS: bool = (
    _DEFAULT_CONFIG.harvesting.include_ungeoreferenced_state_vouchers
)


def resolve_cpu_threads(
    cli_val: Optional[int] = None,
    cfg: Optional[PipelineConfig] = None,
    default: Optional[int] = None,
) -> Optional[int]:
    """
    Resolve CPU worker thread count via priority chain:
    CLI (--threads / -j) -> ENV (PLANTMATICS_MAX_WORKERS / PACKERA_MAX_WORKERS) -> config.yaml (compute.threads) -> default.
    """
    if cli_val is not None:
        try:
            v = int(cli_val)
            if v >= 1:
                return v
        except (TypeError, ValueError):
            pass
    env_val = (
        os.environ.get("PLANTMATICS_MAX_WORKERS", "").strip()
        or os.environ.get("PACKERA_MAX_WORKERS", "").strip()
    )
    if env_val:
        try:
            v = int(env_val)
            if v >= 1:
                return v
        except (TypeError, ValueError):
            pass
    _cfg_mod = sys.modules.get("scripts.core.config")
    active_cfg = cfg if cfg is not None else getattr(_cfg_mod, "_DEFAULT_CONFIG", _DEFAULT_CONFIG)
    if active_cfg.compute.threads is not None and int(active_cfg.compute.threads) >= 1:
        return int(active_cfg.compute.threads)
    return default


def resolve_max_workers(
    default: Optional[int] = None,
    cap: Optional[int] = None,
    cfg: Optional[PipelineConfig] = None,
) -> int:
    """Resolve bounded CPU worker count, collapsing to 1 inside forked workers to block nested forks."""
    if os.environ.get("PLANTMATICS_IN_FORK_WORKER", "").strip() == "1":
        return 1
    env_val = (
        os.environ.get("PLANTMATICS_MAX_WORKERS", "").strip()
        or os.environ.get("PACKERA_MAX_WORKERS", "").strip()
    )
    resolved: Optional[int] = None
    if env_val:
        try:
            parsed = int(env_val)
            if parsed >= 1:
                resolved = parsed
        except ValueError:
            resolved = None
    if resolved is None:
        resolved = resolve_cpu_threads(None, cfg=cfg, default=default)
    if resolved is None:
        resolved = os.cpu_count() or 1
    out = max(1, int(resolved))
    if cap is not None and int(cap) >= 1:
        out = min(out, int(cap))
    return out


def resolve_gpu_batch_size(
    cli_val: Optional[int] = None,
    cfg: Optional[PipelineConfig] = None,
    default_cuda: int = 64,
    default_cpu: int = 32,
    use_cuda: bool = True,
    default: Optional[int] = None,
) -> int:
    """
    Resolve GPU/inference batch size via priority chain:
    CLI (--gpu-batch-size / --batch-size) -> ENV (PLANTMATICS_GPU_BATCH_SIZE / PACKERA_GPU_BATCH_SIZE) ->
    config.yaml (compute.gpu_batch_size) -> stage-optimal default.
    """
    if cli_val is not None:
        try:
            v = int(cli_val)
            if v >= 1:
                return v
        except (TypeError, ValueError):
            pass
    env_val = (
        os.environ.get("PLANTMATICS_GPU_BATCH_SIZE", "").strip()
        or os.environ.get("PACKERA_GPU_BATCH_SIZE", "").strip()
    )
    if env_val:
        try:
            v = int(env_val)
            if v >= 1:
                return v
        except (TypeError, ValueError):
            pass
    _cfg_mod = sys.modules.get("scripts.core.config")
    active_cfg = cfg if cfg is not None else getattr(_cfg_mod, "_DEFAULT_CONFIG", _DEFAULT_CONFIG)
    if active_cfg.compute.gpu_batch_size is not None and int(active_cfg.compute.gpu_batch_size) >= 1:
        return int(active_cfg.compute.gpu_batch_size)
    if default is not None:
        return max(1, int(default))
    return max(1, int(default_cuda) if use_cuda else int(default_cpu))


def _resolve_taxa_string(
    explicit: Optional[str],
    env_key: str,
    cfg: Optional[PipelineConfig],
    attr_name: str,
    default: str,
) -> str:
    """Resolve a target taxonomic rank via explicit -> env -> config/YAML -> default."""
    if explicit is not None and str(explicit).strip():
        return str(explicit).strip()
    env_val = os.environ.get(env_key, "").strip()
    if env_val:
        return env_val
    if cfg is not None:
        cfg_val = str(getattr(cfg.taxa, attr_name, "") or "").strip()
        if cfg_val:
            return cfg_val
    _cfg_mod = sys.modules.get("scripts.core.config")
    _default_cfg_path = getattr(_cfg_mod, "DEFAULT_CONFIG_PATH", DEFAULT_CONFIG_PATH)
    if _default_cfg_path.is_file():
        try:
            with _default_cfg_path.open("r", encoding="utf-8") as fh:
                raw_cfg = yaml.safe_load(fh) or {}
            yaml_val = str(((raw_cfg.get("taxa") or {}).get(attr_name, "")) or "").strip()
            if yaml_val:
                return yaml_val
        except (OSError, yaml.YAMLError):
            pass
    _def_cfg = getattr(_cfg_mod, "_DEFAULT_CONFIG", _DEFAULT_CONFIG)
    default_cfg_val = str(getattr(_def_cfg.taxa, attr_name, "") or "").strip()
    return default_cfg_val or default


def resolve_target_genus(
    cfg: Optional[PipelineConfig] = None,
    explicit_genus: Optional[str] = None,
) -> str:
    """Resolve the active target genus with explicit, environment, YAML, and default fallback."""
    return _resolve_taxa_string(
        explicit=explicit_genus,
        env_key="PLANTMATICS_TARGET_GENUS",
        cfg=cfg,
        attr_name="target_genus",
        default="Packera",
    )


def resolve_target_family(
    cfg: Optional[PipelineConfig] = None,
    explicit_family: Optional[str] = None,
) -> str:
    """Resolve the active target family with explicit, environment, YAML, and default fallback."""
    return _resolve_taxa_string(
        explicit=explicit_family,
        env_key="PLANTMATICS_TARGET_FAMILY",
        cfg=cfg,
        attr_name="target_family",
        default="Asteraceae",
    )


_DEFAULT_GENOMICS_CONFIG: Dict[str, Any] = {
    "min_hq_reads": 50000,
    "max_missing_rate": 0.80,
    "paralog_ho_threshold": 0.60,
    "hdplot_d_threshold": 4.0,
    "nquire_min_het_loci": 60,
    "nquire_strong_delta_bic": 25.0,
    "nquire_min_delta_logl_per_locus": 0.15,
    "default_fallback_ploidy": 2,
    "cytonuclear_min_containment": 0.80,
    "cytonuclear_min_shared_kmers": 80,
    "coalescent_block_len": 122,
    "coalescent_n_boot": 250,
}


def resolve_genomics_config(
    cfg: Optional[PipelineConfig] = None,
    config_path: Optional[Path] = None,
) -> Dict[str, Any]:
    """Resolve Step 06 genomic QC, cytotype, organellar, and coalescent settings from YAML/defaults."""
    out = dict(_DEFAULT_GENOMICS_CONFIG)
    raw_sec: Any = cfg.get("genomics") if cfg is not None else None
    if not isinstance(raw_sec, dict):
        _cfg_mod = sys.modules.get("scripts.core.config")
        _default_cfg_path = getattr(_cfg_mod, "DEFAULT_CONFIG_PATH", DEFAULT_CONFIG_PATH)
        yaml_p = Path(config_path) if config_path is not None else _default_cfg_path
        if yaml_p.is_file():
            try:
                loaded = yaml.safe_load(yaml_p.read_text(encoding="utf-8")) or {}
                raw_sec = loaded.get("genomics")
            except (OSError, yaml.YAMLError):
                raw_sec = None
    if isinstance(raw_sec, dict):
        for k, default_v in _DEFAULT_GENOMICS_CONFIG.items():
            if k in raw_sec and raw_sec[k] is not None:
                try:
                    out[k] = type(default_v)(raw_sec[k])
                except (TypeError, ValueError):
                    pass
    return out
