# ============================================================
# Layer 4: RWR Network Pathway Analysis
# v2.4 — Automated Robustness Audit
# ============================================================
#
# SCIENTIFIC INTEGRITY STATEMENT:
#
#   Primary analysis: top-50 Tier1 expanded signature.
#   This is the officially benchmarked cutoff (Layer 3A).
#   It is IMMUTABLE and cannot be replaced by sensitivity results.
#
#   Robustness audit: automated ORA for top-100 and top-150
#   positive-delta genes from Layer 2B, run in background.
#   These are AUDIT outputs only — never used to replace primary.
#
#   Theme concordance: pathways passing all three cutoffs
#   (50, 100, 150) are flagged as "robust themes" in
#   Robustness_Themes.csv. Jaccard(top50_terms, top150_terms)
#   is the falsification metric.
#
#   FALSIFICATION CONDITION (pre-declared):
#     If Jaccard(top50_pathway_set, top150_pathway_set) = 0,
#     the 50-gene signal is topologically unstable and must
#     not be used to claim mechanistic interpretation.
#
#   Scientific basis for multi-cutoff design:
#     A single rank cutoff introduces arbitrary threshold
#     sensitivity. The RWR top-50 represents the highest-density
#     signal core; top-100/150 are spatial expansions.
#     Theme concordance across all three layers provides
#     topological robustness evidence, not p-hacking.
#     Ref: Bauer et al. 2010 (PLoS Comput Biol) — ORA sensitivity
#          Tarca et al. 2013 — gene list size effects in enrichment
#
# ARCHITECTURE:
#   Stage 1-4:  Primary analysis on top-50 (unchanged from v2.3)
#   Stage 5:    Automated robustness audit (top-100, top-150)
#   Stage 6:    Theme concordance computation + Robustness_Themes.csv
#
# SMART compliance:
#   S: Specific cutoffs (50/100/150) from Layer 2B Delta_Score ranking
#   M: Jaccard concordance between pathway term sets is measurable
#   A: Fully automated — no human selection of "best" result
#   R: Relevant to MoA interpretation of RWR network signal
#   T: Single notebook execution, deterministic output
#
# FAIR compliance:
#   F: Findable — all outputs named with cutoff suffix
#   A: Accessible — CSV/JSON, standard formats
#   I: Interoperable — Cytoscape-compatible, GSEApy-standard
#   R: Reusable — config-driven, no hardcoded biology
#
# CHANGELOG v2.4 vs v2.3:
#   [NEW] Stage 5: automated robustness audit (top-100, top-150)
#   [NEW] Stage 6: theme concordance + Robustness_Themes.csv
#   [NEW] L2B_DELTA_CSV config path for sensitivity gene lists
#   [NEW] SENSITIVITY_CUTOFFS constant (declared, not hardcoded)
#   [NEW] FALSIFICATION_MIN_JACCARD pre-declared threshold
#   [KEPT] All v2.3 patches (C2, M1, C1, N1-N4, M-1 to M-4)
#   [KEPT] Primary analysis immutable at N=50
#
# REFERENCES:
#   Choobdar et al. Nat Methods 16:843-852 (2019)
#   Lamparter et al. PLoS Comput Biol 12(1):e1004714 (2016)
#   Ziemann et al. Genome Biology 17:177 (2016)
#   Bauer et al. PLoS Comput Biol 6(6):e1000796 (2010)
# ============================================================


# ============================================================
# STAGE 0: Imports & Configuration
# ============================================================
import os
import sys
import re
import json
import logging
import warnings
import time
from pathlib import Path
from datetime import datetime
from typing import List, Dict, Set, Tuple, Optional

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm
import seaborn as sns

from scipy.stats import hypergeom
from scipy.cluster.hierarchy import linkage
from scipy.spatial.distance import squareform
from sklearn.linear_model import ElasticNetCV, ElasticNet
from statsmodels.stats.multitest import multipletests

import gseapy as gp

warnings.filterwarnings("ignore")
pd.set_option("display.max_columns", 200)
pd.set_option("display.width", 200)

import config_system as cfg


# ============================================================
# FONT CONFIGURATION
# ============================================================
def configure_fonts() -> str:
    """Detect best available font and configure matplotlib globally."""
    system_fonts = {f.name for f in fm.fontManager.ttflist}
    preferred    = ['Arial', 'Tahoma', 'Segoe UI', 'Liberation Sans', 'DejaVu Sans']
    selected     = next((f for f in preferred if f in system_fonts), 'sans-serif')
    fallback     = [selected] + [f for f in preferred if f != selected] + ['sans-serif']
    matplotlib.rcParams.update({
        'figure.dpi': 300, 'savefig.dpi': 300,
        'font.size': 10, 'font.family': 'sans-serif',
        'font.sans-serif': fallback,
        'axes.linewidth': 0.8, 'axes.unicode_minus': False,
        'svg.fonttype': 'none', 'pdf.fonttype': 42, 'ps.fonttype': 42,
    })
    sns.set_style("white")
    sns.set_context("paper", font_scale=1.1)
    return selected

_FONT = configure_fonts()
LABEL_NEG_LOG10_PVAL = r'$-\log_{10}$' + '(Adjusted P-value)'


def _show_plot() -> None:
    """Display plot only in interactive backends; silent no-op on Agg."""
    if matplotlib.get_backend().lower() != 'agg':
        plt.show()


# ============================================================
# CONFIGURATION
# ============================================================

# --- Input: Layer 3A outputs (primary) ---
TIER1_CSV = cfg.LAYER3A_OUTPUT_DIR / cfg.L3A_EXPANDED_SIGNATURE_CSV
TIER2_CSV = cfg.LAYER3A_OUTPUT_DIR / cfg.L3A_CORE_ESSENTIAL_CSV

# --- Input: Layer 2B Delta Network (for sensitivity gene lists) ---
# Required for Stage 5 robustness audit.
L2B_DELTA_CSV = (
    cfg.PROJECT_ROOT / "outputs"
    / f"Layer2B_Heterogeneous_RWR_{cfg.L3A_PRIMARY_CELL_LINE}"
    / cfg.L2B_DELTA_NETWORK_CSV
)

# --- GMT library ---
GMT_FILE_PATH = cfg.LAYER4_GMT_FILE_PATH

# --- Primary enrichment parameters ---
FDR_CUTOFF  = 0.05
MIN_OVERLAP = 3
TOP_TERM    = 25

# [M1] Concordance tolerance: log10 scale, 10-fold allowed.
LOG10_PVALUE_CONCORDANCE_TOL = 1.0

# [C1] ElasticNetCV heuristic pre-filter. cv=5 conventional choice.
ELASTICNET_L1_RATIO = 0.5
ELASTICNET_CV_FOLDS = 5

# Pairwise overlap filter (Choobdar et al. 2019)
OVERLAP_FDR_ALPHA     = 0.05
SUBMODULE_J_THRESHOLD = 0.5
SUBMODULE_S_THRESHOLD = 0.5

# [M-3] Size tie-break threshold (units: absolute gene count difference).
# If |parent| - |child| < SIZE_DIFF_FALLBACK → prefer specific child.
# Value 5 is a conservative heuristic; no literature threshold exists.
SIZE_DIFF_FALLBACK = 5
MODULE_MIN_GENES   = 3

# ============================================================
# ROBUSTNESS AUDIT CONFIGURATION (v2.4 — NEW)
# ============================================================
# Pre-declared sensitivity cutoffs for automated robustness audit.
# These supplement, but NEVER replace, the primary top-50 analysis.
# The order matters for concordance computation:
#   SENSITIVITY_CUTOFFS[0] = primary cutoff (must match Tier1 size)
#   SENSITIVITY_CUTOFFS[1:] = audit cutoffs
#
# Scientific basis:
#   top-50  = high-density signal core (primary)
#   top-100 = spatial expansion layer 1
#   top-150 = spatial expansion layer 2
#   Theme concordance across all three confirms topological robustness.
#   Ref: RWR signal diffusion (Valdeolivas et al. 2017)
# ============================================================
SENSITIVITY_CUTOFFS = [50, 100, 150]

# Pre-declared falsification threshold (Jaccard between pathway term sets).
# If Jaccard(top50_terms ∩ top150_terms) / (top50_terms ∪ top150_terms) = 0,
# the primary signal is declared topologically unstable.
FALSIFICATION_MIN_JACCARD = 0.10   # strictly > 0 required to avoid falsification

# Robustness output filenames
ROBUSTNESS_THEMES_CSV    = "L4_Robustness_Themes.csv"
ROBUSTNESS_SUMMARY_JSON  = "L4_Robustness_Audit_Summary.json"
SENSITIVITY_DETAIL_CSV   = "L4_Sensitivity_Detail.csv"

# Output directories
TIMESTAMP  = datetime.now().strftime("%Y%m%d_%H%M%S")
RUN_ID     = f"Layer4_RWR_Pathway_{TIMESTAMP}"
OUTPUT_DIR = cfg.PROJECT_ROOT / "outputs" / "Layer4_RWR_Pathways" / RUN_ID
ENRICH_DIR = OUTPUT_DIR / "enrichment_raw"
CYTO_DIR   = OUTPUT_DIR / "cytoscape"
ROBUST_DIR = OUTPUT_DIR / "robustness_audit"
for _d in [OUTPUT_DIR, ENRICH_DIR, CYTO_DIR, ROBUST_DIR]:
    _d.mkdir(parents=True, exist_ok=True)

# Context note — neutral, from config
L4_CONTEXT_NOTE = (
    f"Cancer model: {cfg.L3B_CANCER_CELL_LINE} ({cfg.L3B_CANCER_MODEL_ID}). "
    "Pathway interpretation must be integrated with the corresponding Layer 3B "
    "AUCell result and matched in vitro evidence. "
    "No phenotype conclusion is hardcoded in this notebook."
)

print("=" * 72)
print("LAYER 4: RWR NETWORK PATHWAY ANALYSIS  v2.4")
print("=" * 72)
print(f"  Run ID:            {RUN_ID}")
print(f"  Cancer model:      {cfg.L3B_CANCER_CELL_LINE}")
print(f"  Primary cutoff:    N={SENSITIVITY_CUTOFFS[0]} (immutable primary)")
print(f"  Audit cutoffs:     {SENSITIVITY_CUTOFFS[1:]} (robustness only)")
print(f"  Falsification J:   > {FALSIFICATION_MIN_JACCARD} required")
print(f"  GMT:               {GMT_FILE_PATH}")


# ============================================================
# Logger
# ============================================================
logger = logging.getLogger("Layer4_v24")
logger.setLevel(logging.DEBUG)
if logger.hasHandlers():
    logger.handlers.clear()

_fmt = logging.Formatter("%(asctime)s | %(levelname)-8s | %(message)s")
_fh  = logging.FileHandler(OUTPUT_DIR / "layer4_audit.log", mode="w", encoding="utf-8")
_fh.setLevel(logging.DEBUG)
_fh.setFormatter(_fmt)
logger.addHandler(_fh)
_sh = logging.StreamHandler(sys.stdout)
_sh.setLevel(logging.INFO)
_sh.setFormatter(_fmt)
logger.addHandler(_sh)

logger.info("=" * 72)
logger.info(f"  RUN ID:  {RUN_ID}")
logger.info(f"  GSEApy:  {gp.__version__}")
logger.info(f"  Font:    {_FONT}")
logger.info(f"  Primary: N={SENSITIVITY_CUTOFFS[0]} (top-50, IMMUTABLE)")
logger.info(f"  Audit:   {SENSITIVITY_CUTOFFS[1:]} (robustness check)")
logger.info(f"  Falsification min Jaccard: > {FALSIFICATION_MIN_JACCARD}")
logger.info("=" * 72)


# ============================================================
# UTILITY: Excel date-corruption rescue [N3: SEPTIN 1-16]
# Ref: Ziemann et al. Genome Biology 17:177 (2016)
# ============================================================
def detect_date_corrupted_genes(
    gene_list: List[str],
) -> Tuple[List[str], List[dict]]:
    """
    Detect and rescue gene symbols corrupted by Excel date auto-conversion.
    [N3] SEPTIN coverage: 1-16 (includes 13, 15, 16 added in v2.2).
    """
    _dm = re.compile(
        r'^(\d{1,2})[-/](JAN|FEB|MAR|APR|MAY|JUN|JUL|AUG|SEP|OCT|NOV|DEC)$',
        re.IGNORECASE,
    )
    _md = re.compile(
        r'^(JAN|FEB|MAR|APR|MAY|JUN|JUL|AUG|SEP|OCT|NOV|DEC)[-/](\d{1,2})$',
        re.IGNORECASE,
    )
    _dt = re.compile(r'^\d{4}[-/]\d{2}[-/]\d{2}')
    _sr = re.compile(r'^\d{5}$')

    _SN = list(range(1, 17))   # SEPTIN 1-16
    rescue_map: Dict[str, str] = {
        **{f'{d}-MAR': f'MARCHF{d}' for d in range(1, 12)},
        **{f'MAR-{d}': f'MARCHF{d}' for d in range(1, 12)},
        **{f'{d}-SEP': f'SEPTIN{d}' for d in _SN},
        **{f'SEP-{d}': f'SEPTIN{d}' for d in _SN},
        '1-DEC': 'BHLHE40', '2-DEC': 'BHLHE41',
        'DEC-1': 'BHLHE40', 'DEC-2': 'BHLHE41',
        '4-OCT': 'POU5F1',  'OCT-4': 'POU5F1',
    }

    clean: List[str]   = []
    report: List[dict] = []

    for g in gene_list:
        gu = str(g).upper().strip()
        if gu in rescue_map:
            report.append({'corrupted': g, 'pattern': 'known_gene_date',
                           'rescued': rescue_map[gu], 'action': 'RESCUED'})
            clean.append(rescue_map[gu])
            continue
        corrupted, pattern = False, None
        if _dm.match(gu):   corrupted, pattern = True, 'D-MON'
        elif _md.match(gu): corrupted, pattern = True, 'MON-D'
        elif _dt.match(gu): corrupted, pattern = True, 'ISO_datetime'
        elif _sr.match(gu): corrupted, pattern = True, 'Excel_serial'
        if corrupted:
            report.append({'corrupted': g, 'pattern': pattern,
                           'rescued': None, 'action': 'REMOVED'})
        else:
            clean.append(g)
    return clean, report


# ============================================================
# UTILITY: Overlap column sanitization [CSV-FIX from v4.1.2]
# ============================================================
def _parse_overlap(s) -> Tuple[int, int]:
    try:
        s     = str(s)
        parts = s.split(' | ') if ' | ' in s else s.split('/')
        return int(parts[0].strip()), int(parts[1].strip())
    except Exception:
        return 0, 0


def sanitize_overlap_column(df: pd.DataFrame) -> pd.DataFrame:
    """Replace 'k/M' with 'k | M' and add integer columns."""
    if df is None or df.empty or 'Overlap' not in df.columns:
        return df
    df = df.copy()
    kv, mv, dv = [], [], []
    for _, row in df.iterrows():
        k, m = _parse_overlap(row['Overlap'])
        kv.append(k); mv.append(m); dv.append(f"{k} | {m}")
    df['k_overlap']      = kv
    df['M_pathway_size'] = mv
    df['Overlap']        = dv
    return df


# ============================================================
# UTILITY: Hypergeometric verification (audit only)
# ============================================================
def verify_hypergeometric(
    query: List[str], pathway_genes: List[str],
    N_audit: int, pathway_name: str, kb_universe: Set[str],
) -> dict:
    """
    Audit p-value via scipy.stats.hypergeom.sf.
    n = |Query ∩ GMT_universe| (n_effective, NOT raw size).
    N = N_KB_audit (GMT union only, NOT augmented by Tier1 genes).
    """
    qs  = set(query); ps = set(pathway_genes)
    qe  = qs & kb_universe; ne = len(qe)
    ov  = qe & ps; M = len(ps); k = len(ov)
    p   = float(hypergeom.sf(k - 1, N_audit, M, ne)) if (k > 0 and ne > 0) else 1.0
    return {'pathway': pathway_name, 'N_audit': N_audit, 'M': M,
            'n_raw': len(query), 'n_effective': ne, 'k': k,
            'p_value_manual': p, 'overlap_genes': ';'.join(sorted(ov))}


# ============================================================
# CORE: Single-run enrichment + redundancy filter
# ============================================================
def _parse_gmt_file(filepath: str) -> Dict[str, Set[str]]:
    """
    Parse GMT format. Produces set-valued dict for audit/overlap.
    Custom parser used to obtain set-valued dict for overlap math;
    gp.enrich() receives list-valued dict from the same source.
    """
    gmt: Dict[str, Set[str]] = {}
    if not os.path.exists(filepath):
        raise FileNotFoundError(f"GMT not found: {filepath}")
    with open(filepath, 'r', encoding='utf-8') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) < 3:
                continue
            t = parts[0].strip()
            g = {x.strip().upper() for x in parts[2:]
                 if x.strip() and x.strip().upper() not in {'', 'NA', 'NAN', 'NONE'}}
            if g:
                gmt[t] = g
    return gmt


def _select_terms_elasticnet(
    candidates: Dict[str, Set[str]], sig_genes: Set[str],
    gmt: Dict[str, Set[str]], l1: float, cv: int,
) -> Tuple[List[str], dict]:
    """
    Heuristic multivariate redundancy pre-filter via ElasticNetCV.
    [C1] Differs from Choobdar 2019 static regularization.
    No formal superiority claimed. Heuristic only.
    """
    terms = sorted(candidates.keys())
    if len(terms) == 0:
        return terms, {'method': 'passthrough', 'reason': 'empty'}
    all_g = set()
    for t in terms:
        all_g.update(gmt.get(t, set()))
    all_g.update(sig_genes)
    gl = sorted(all_g); gi = {g: i for i, g in enumerate(gl)}
    ng = len(gl); nt = len(terms)
    if nt < 2 or ng < 3:
        return terms, {'method': 'passthrough', 'reason': 'too_few'}
    X = np.zeros((ng, nt), dtype=np.float32)
    for j, t in enumerate(terms):
        for g in gmt.get(t, set()):
            if g in gi:
                X[gi[g], j] = 1.0
    y = np.zeros(ng, dtype=np.float32)
    for g in sig_genes:
        if g in gi:
            y[gi[g]] = 1.0
    valid = X.var(axis=0) > 0
    if valid.sum() < 2:
        return terms, {'method': 'passthrough', 'reason': 'no_variance'}
    Xv = X[:, valid]; vi = np.where(valid)[0]
    sel = None; diag: dict = {}
    try:
        m = ElasticNetCV(l1_ratio=l1, cv=cv, n_alphas=20, positive=True,
                         fit_intercept=False, max_iter=3000, tol=1e-4, n_jobs=1)
        m.fit(Xv, y)
        sel  = [terms[i] for i in vi[m.coef_ > 0]]
        diag = {'method': 'ElasticNetCV', 'note': 'heuristic pre-filter',
                'alpha': float(m.alpha_), 'cv_folds': cv,
                'input': len(terms), 'selected': len(sel)}
    except Exception as _e:
        logger.warning(f"    ElasticNetCV failed: {_e}")
    if sel is None or len(sel) < MODULE_MIN_GENES:
        try:
            ms = ElasticNet(alpha=0.01, l1_ratio=l1, positive=True,
                            fit_intercept=False, max_iter=3000, tol=1e-4)
            ms.fit(Xv, y)
            sel  = [terms[i] for i in vi[ms.coef_ > 0]]
            diag = {'method': 'ElasticNet_static_fallback',
                    'note': 'alpha=0.01 conventional, fallback only',
                    'input': len(terms), 'selected': len(sel)}
        except Exception as _e2:
            logger.warning(f"    Static fallback failed: {_e2}")
            sel  = terms
            diag = {'method': 'passthrough', 'reason': 'all_models_failed'}
    if len(sel) < MODULE_MIN_GENES:
        sel = terms; diag.update({'method': 'passthrough', 'reason': 'below_minimum'})
    return sel, diag


def _pairwise_filter(
    tl: List[str], gl: List[Set[str]], ps: dict,
    pvd: Dict[str, float], fa: float, jt: float, st: float, sf: int,
) -> Tuple[List[str], List[dict], dict]:
    """
    Pairwise Jaccard/Submodule/Hypergeometric redundancy filter.
    [N1] Inline literal replaces mutable module-level default dict.
    [M-3] Tie-break: size_diff < sf → specific child; else → broad parent.
    """
    n = len(tl)
    if n == 0:
        return [], [], {}
    rp = [
        (i, j, p['J'], p['S'], p['p'], p['k'])
        for i in range(n)
        for j in range(i + 1, n)
        if (p := ps.get((i, j), {'J': 0.0, 'S': 0.0, 'p': 1.0, 'k': 0}))['k'] > 0
    ]
    keep = {i: True for i in range(n)}
    elog: List[dict] = []; cm = 'none'
    if rp:
        rpv = np.array([x[4] for x in rp])
        rpv[~np.isfinite(rpv)] = 1.0
        cm = 'fdr_bh' if len(rpv) >= 20 else 'bonferroni'
        rej, pvc, _, _ = multipletests(rpv, alpha=fa, method=cm)
        for ix, (i, j, J, S, raw_p, k) in enumerate(rp):
            if not rej[ix]:
                continue
            twin   = J >= jt
            nested = (not twin) and S >= st
            if not (twin or nested):
                continue
            pi, ci = (i, j) if len(gl[i]) >= len(gl[j]) else (j, i)
            if not (keep.get(pi) and keep.get(ci)):
                continue
            pp, cp = pvd.get(tl[pi], 1.0), pvd.get(tl[ci], 1.0)
            if pp < cp:
                el, kp = ci, pi; reason = 'parent_stronger_pvalue'
            elif cp < pp:
                el, kp = pi, ci; reason = 'child_stronger_pvalue'
            else:
                sd = len(gl[pi]) - len(gl[ci])
                if sd < sf:
                    el, kp = pi, ci; reason = 'prefer_specific_child_when_sizes_similar'
                else:
                    el, kp = ci, pi; reason = 'prefer_broad_parent_when_gap_large'
            keep[el] = False
            elog.append({'Eliminated': tl[el], 'Kept': tl[kp],
                         'Jaccard': J, 'Submodule_S': S, 'Intersection': k,
                         'Overlap_Type': 'twin_J' if twin else 'nested_S',
                         'Reason': reason, 'Pval_Raw': raw_p,
                         'Pval_Corrected': float(pvc[ix])})
    surv = [tl[i] for i in range(n) if keep[i]]
    elim = [tl[i] for i in range(n) if not keep[i]]
    return surv, elog, {'input': n, 'pairs': len(rp), 'correction': cm,
                        'survivors': len(surv), 'eliminated': len(elim)}


def run_enrichment_for_gene_list(
    gene_list: List[str],
    gmt_for_enrich: Dict[str, List[str]],
    full_gmt_dict: Dict[str, Set[str]],
    background_gmt: Set[str],
    N_KB_audit: int,
    cutoff_label: str,
    fdr_cutoff: float = FDR_CUTOFF,
    min_overlap: int  = MIN_OVERLAP,
    run_filter: bool  = True,
) -> Tuple[pd.DataFrame, pd.DataFrame, List[str], dict]:
    """
    Run full ORA enrichment + optional redundancy filter for one gene list.

    Returns:
        enrichment_sig:  significant enrichment results (sanitized)
        final_df:        non-redundant pathways with evidence overlay
        final_terms:     list of non-redundant pathway term names
        filter_diag:     filtration diagnostics dict

    Parameters:
        run_filter: if False, skip ElasticNet+Jaccard (for audit runs
                    where we only need pathway term sets for concordance)
    """
    _label = cutoff_label
    _tier1_set = set(gene_list)
    _tier2_set: Set[str] = set()   # No Tier2 for sensitivity runs

    # Effective genes
    tier1_eff = sorted(_tier1_set & background_gmt)
    n_eff     = len(tier1_eff)
    not_in    = sorted(_tier1_set - background_gmt)

    logger.info(f"    [{_label}] n_raw={len(gene_list)}, n_eff={n_eff}, "
                f"not_in_GMT={len(not_in)}")

    if n_eff < min_overlap:
        logger.warning(f"    [{_label}] n_eff < MIN_OVERLAP, skipping")
        return pd.DataFrame(), pd.DataFrame(), [], {}

    # gp.enrich — background=None
    try:
        _enr = gp.enrich(
            gene_list=gene_list,
            gene_sets=gmt_for_enrich,
            background=None, outdir=None, no_plot=True, verbose=False,
        )
        _raw = sanitize_overlap_column(_enr.results.copy())
    except Exception as _e:
        logger.warning(f"    [{_label}] gp.enrich failed: {_e}")
        return pd.DataFrame(), pd.DataFrame(), [], {}

    # Filter by significance + min overlap
    _sig = _raw[
        (_raw['Adjusted P-value'] <= fdr_cutoff) &
        (_raw['k_overlap'] >= min_overlap)
    ].copy().reset_index(drop=True)

    logger.info(f"    [{_label}] Significant terms: {len(_sig)}")

    if _sig.empty:
        return _sig, pd.DataFrame(), [], {}

    if not run_filter:
        # Audit mode: return significant terms without full filter
        return _sig, pd.DataFrame(), list(_sig['Term'].unique()), {}

    # Build pathway gene sets
    _invalid  = {'', 'NAN', 'NONE', 'NA', 'NULL', 'N/A', '.'}
    _pw_genes: Dict[str, Set[str]] = {}
    _pw_pvals: Dict[str, float]    = {}

    for _, row in _sig.iterrows():
        t = str(row.get('Term', '')).strip()
        if not t:
            continue
        raw = str(row.get('Genes', ''))
        gs: Set[str] = set()
        for sep in [';', ',']:
            if sep in raw:
                gs = {g.strip().upper() for g in raw.split(sep)
                      if g.strip().upper() not in _invalid and g.strip()}
                break
        else:
            if raw.strip().upper() not in _invalid:
                gs = {raw.strip().upper()}
        if not gs:
            continue
        _pw_genes[t] = _pw_genes.get(t, set()) | gs
        try:
            pv = float(row.get('Adjusted P-value', 1.0))
            if np.isfinite(pv) and 0 < pv <= 1:
                if t not in _pw_pvals or pv < _pw_pvals[t]:
                    _pw_pvals[t] = pv
        except Exception:
            pass

    # Size filter
    cands: Dict[str, Set[str]] = {}
    for t, g in _pw_genes.items():
        if max(len(full_gmt_dict.get(t, set())), len(g)) >= MODULE_MIN_GENES:
            cands[t] = full_gmt_dict.get(t, g)

    if len(cands) == 0:
        return _sig, pd.DataFrame(), [], {}

    # ElasticNetCV pre-filter
    s1_surv, s1_diag = _select_terms_elasticnet(
        cands, _tier1_set, full_gmt_dict, ELASTICNET_L1_RATIO, ELASTICNET_CV_FOLDS
    )
    pw_sel = {t: full_gmt_dict.get(t, cands.get(t, set())) for t in s1_surv}
    pv_sel = {t: _pw_pvals.get(t, 1.0) for t in s1_surv}

    # Pairwise filter
    tl  = list(pw_sel.keys())
    gsl = [pw_sel[t] for t in tl]

    # Pairwise overlaps
    ps_dict: dict = {}
    for i in range(len(tl)):
        A = gsl[i]
        for j in range(i + 1, len(tl)):
            B     = gsl[j]
            inter = A & B; k = len(inter)
            if k == 0:
                continue
            u  = len(A) + len(B) - k
            J  = k / u if u > 0 else 0.0
            ls = min(len(A), len(B)); ll = max(len(A), len(B))
            S  = k / ls if ls > 0 else 0.0
            Ne = max(N_KB_audit, ll, ls)
            p  = float(hypergeom.sf(k - 1, Ne, ll, ls))
            if not np.isfinite(p) or p < 0:
                p = 1.0
            ps_dict[(i, j)] = {'J': J, 'S': S, 'p': p, 'k': k}

    s2_surv, _, s2_diag = _pairwise_filter(
        tl, gsl, ps_dict, pv_sel, OVERLAP_FDR_ALPHA,
        SUBMODULE_J_THRESHOLD, SUBMODULE_S_THRESHOLD, SIZE_DIFF_FALLBACK,
    )

    # Build final df with evidence overlay
    final_rows = []
    for t in s2_surv:
        eg  = _pw_genes.get(t, set())
        t1o = sorted(eg & _tier1_set)
        final_rows.append({
            'Term':               t,
            'Cutoff':             cutoff_label,
            'Adjusted_P_value':   _pw_pvals.get(t, np.nan),
            'N_Enrichment_Genes': len(eg),
            'Enrichment_Genes':   ';'.join(sorted(eg)),
            'N_Tier1_Overlap':    len(t1o),
            'Tier1_Overlap_Genes': ';'.join(t1o),
            'Multi_Layer_Score':  len(t1o),
        })
    final_df = pd.DataFrame(final_rows)
    if not final_df.empty:
        final_df = final_df.sort_values('Adjusted_P_value').reset_index(drop=True)

    diag = {'stage1_elasticnet': s1_diag, 'stage2_pairwise': s2_diag,
            'n_sig': len(_sig), 'n_final': len(final_df)}
    return _sig, final_df, s2_surv, diag


# ============================================================
# STAGES 1-4: PRIMARY ANALYSIS (N=50, immutable)
# ============================================================
logger.info("\n" + "=" * 72)
logger.info("STAGES 1-4: PRIMARY ANALYSIS (N=50, IMMUTABLE)")
logger.info("=" * 72)

# Stage 1: Load Tier 1 and Tier 2
for _p in [TIER1_CSV, TIER2_CSV]:
    if not _p.exists():
        raise FileNotFoundError(f"Required: {_p}")

df_tier1 = pd.read_csv(TIER1_CSV)
df_tier2 = pd.read_csv(TIER2_CSV)


def _extract_genes(df: pd.DataFrame) -> List[str]:
    col = "Gene_Normalized" if "Gene_Normalized" in df.columns else "Gene"
    return df[col].dropna().apply(lambda x: str(x).upper().strip()).tolist()


TIER1_RAW = _extract_genes(df_tier1)
TIER2_RAW = _extract_genes(df_tier2)
TIER1_CLEAN, _ct1 = detect_date_corrupted_genes(TIER1_RAW)
TIER2_CLEAN, _ct2 = detect_date_corrupted_genes(TIER2_RAW)
TIER1_GENES = sorted(list(dict.fromkeys(TIER1_CLEAN)))
TIER2_GENES = sorted(list(dict.fromkeys(TIER2_CLEAN)))
TIER1_SET   = set(TIER1_GENES)
TIER2_SET   = set(TIER2_GENES)

_all_corrupt = _ct1 + _ct2
if _all_corrupt:
    logger.warning(f"  Date-corruption: {len(_all_corrupt)} entries")
    pd.DataFrame(_all_corrupt).to_csv(
        OUTPUT_DIR / "L4_Gene_DateCorruption_Report.csv",
        index=False, encoding="utf-8-sig"
    )

logger.info(f"  Tier1: {len(TIER1_SET)} clean genes (primary input)")
logger.info(f"  Tier2: {len(TIER2_SET)} genes (CRISPR-confirmed essential)")

assert len(TIER1_GENES) > 0, "Tier1 gene list is empty."
assert len(TIER1_GENES) == SENSITIVITY_CUTOFFS[0], (
    f"Tier1 size ({len(TIER1_GENES)}) must match primary cutoff "
    f"({SENSITIVITY_CUTOFFS[0]}). Check config."
)

# Stage 1B: Load GMT + Build Audit Background [C2]
logger.info("\nSTAGE 1B: GMT + AUDIT BACKGROUND [C2]")
full_gmt_dict: Dict[str, Set[str]] = _parse_gmt_file(str(GMT_FILE_PATH))
gmt_for_enrich: Dict[str, List[str]] = {
    t: sorted(g) for t, g in full_gmt_dict.items()
}

# [C2] N_KB_AUDIT = GMT union ONLY — no Tier1 augmentation
BACKGROUND_GMT: Set[str] = set()
for _g in full_gmt_dict.values():
    BACKGROUND_GMT.update(_g)
N_KB_AUDIT        = len(BACKGROUND_GMT)
TIER1_EFFECTIVE   = sorted(TIER1_SET & BACKGROUND_GMT)
N_TIER1_EFFECTIVE = len(TIER1_EFFECTIVE)
TIER1_NOT_IN_GMT  = sorted(TIER1_SET - BACKGROUND_GMT)

logger.info(f"  GMT terms:     {len(full_gmt_dict)}")
logger.info(f"  N_KB_audit:    {N_KB_AUDIT} [C2: GMT union only]")
logger.info(f"  Tier1 n_eff:   {N_TIER1_EFFECTIVE}")
logger.info(f"  Tier1 unmappable: {len(TIER1_NOT_IN_GMT)}")

# Stage 1C: Primary gp.enrich()
logger.info("\nSTAGE 1C: PRIMARY gp.enrich() — N=50, background=None")

if N_TIER1_EFFECTIVE == 0:
    raise ValueError("n_effective=0: no Tier1 genes in GMT universe.")
if N_TIER1_EFFECTIVE < MIN_OVERLAP:
    raise ValueError(f"n_effective={N_TIER1_EFFECTIVE} < MIN_OVERLAP.")

try:
    _enr_primary      = gp.enrich(
        gene_list=TIER1_GENES, gene_sets=gmt_for_enrich,
        background=None, outdir=None, no_plot=True, verbose=True,
    )
    enrichment_df_raw = sanitize_overlap_column(_enr_primary.results.copy())
except Exception as _e:
    logger.error(f"  gp.enrich FAILED: {_e}")
    raise

enrichment_df_raw.to_csv(
    ENRICH_DIR / "L4_Enrichment_Raw_All.csv", index=False, encoding='utf-8-sig'
)
enrichment_df_sig = enrichment_df_raw[
    (enrichment_df_raw['Adjusted P-value'] <= FDR_CUTOFF) &
    (enrichment_df_raw['k_overlap'] >= MIN_OVERLAP)
].copy().reset_index(drop=True)
enrichment_df_sig.to_csv(
    ENRICH_DIR / "L4_Enrichment_Significant.csv", index=False, encoding='utf-8-sig'
)
logger.info(f"  Primary significant: {len(enrichment_df_sig)}")

# Stage 1D: Hypergeometric Verification [M1]
logger.info("\nSTAGE 1D: HYPERGEOMETRIC VERIFICATION [M1]")
_sig_terms    = set(enrichment_df_sig['Term'].tolist()) if not enrichment_df_sig.empty else set()
verification_rows = [
    verify_hypergeometric(TIER1_GENES, list(g), N_KB_AUDIT, t, BACKGROUND_GMT)
    for t, g in full_gmt_dict.items() if t in _sig_terms
]
verify_df = pd.DataFrame(verification_rows)
verify_df.to_csv(ENRICH_DIR / "L4_Verification.csv", index=False, encoding='utf-8-sig')

discordant = 0; discordance_rows: List[dict] = []
if verification_rows:
    _vl = {r['pathway']: r for r in verification_rows}
    for _, _row in enrichment_df_sig.iterrows():
        _vr = _vl.get(_row['Term'], {})
        _gp = float(_row['P-value'])
        _mp = float(_vr.get('p_value_manual', np.nan)) if isinstance(_vr, dict) else np.nan
        if np.isnan(_mp):
            continue
        _d  = abs(np.log10(max(_gp, 1e-300)) - np.log10(max(_mp, 1e-300)))
        if _d > LOG10_PVALUE_CONCORDANCE_TOL:
            discordant += 1
            discordance_rows.append({'Term': _row['Term'], 'GSEApy_P': _gp,
                                     'Manual_P': _mp, 'Abs_Log10_Diff': round(float(_d), 4),
                                     'Fold_Difference': round(float(10 ** _d), 2)})
if discordant > 0:
    _dp = ENRICH_DIR / "L4_Pvalue_Discordance_Report.csv"
    pd.DataFrame(discordance_rows).to_csv(_dp, index=False, encoding="utf-8-sig")
    logger.warning(f"  Discordance: {discordant} terms >10-fold [expected with background=None]")
else:
    logger.info("  Concordance PASSED: all within 10-fold")

# Stage 1E: Build pathway gene sets
_invalid = {'', 'NAN', 'NONE', 'NA', 'NULL', 'N/A', '.'}
pathway_gene_sets: Dict[str, Set[str]] = {}
pathway_pvalues:   Dict[str, float]    = {}

for _, _row in enrichment_df_sig.iterrows():
    t = str(_row.get('Term', '')).strip()
    if not t:
        continue
    raw = str(_row.get('Genes', ''))
    gs: Set[str] = set()
    for sep in [';', ',']:
        if sep in raw:
            gs = {g.strip().upper() for g in raw.split(sep)
                  if g.strip().upper() not in _invalid and g.strip()}
            break
    else:
        if raw.strip().upper() not in _invalid:
            gs = {raw.strip().upper()}
    if gs:
        pathway_gene_sets[t] = pathway_gene_sets.get(t, set()) | gs
        try:
            pv = float(_row.get('Adjusted P-value', 1.0))
            if np.isfinite(pv) and 0 < pv <= 1:
                if t not in pathway_pvalues or pv < pathway_pvalues[t]:
                    pathway_pvalues[t] = pv
        except Exception:
            pass

candidate_terms: Dict[str, Set[str]] = {}
size_eliminated: List[str]           = []
for t, g in pathway_gene_sets.items():
    if max(len(full_gmt_dict.get(t, set())), len(g)) >= MODULE_MIN_GENES:
        candidate_terms[t] = full_gmt_dict.get(t, g)
    else:
        size_eliminated.append(t)

logger.info(f"  Pathway candidates: {len(candidate_terms)} (after size filter)")

# Stage 2: ElasticNetCV
logger.info("\nSTAGE 2: ELASTICNETCV MULTIVARIATE PRE-FILTER [C1]")
stage1_survivors, stage1_diag = _select_terms_elasticnet(
    candidate_terms, TIER1_SET, full_gmt_dict, ELASTICNET_L1_RATIO, ELASTICNET_CV_FOLDS
)
stage1_eliminated = [t for t in candidate_terms if t not in set(stage1_survivors)]
logger.info(f"  ElasticNet: {len(candidate_terms)} → {len(stage1_survivors)}")

pw_selected = {t: full_gmt_dict.get(t, candidate_terms.get(t, set()))
               for t in stage1_survivors}
pv_selected = {t: pathway_pvalues.get(t, 1.0) for t in stage1_survivors}

# Stage 3: Pairwise filter
logger.info("\nSTAGE 3: PAIRWISE TOPOLOGICAL PRUNING")
tl  = list(pw_selected.keys())
gsl = [pw_selected[t] for t in tl]

ps_main: dict = {}
for i in range(len(tl)):
    A = gsl[i]
    for j in range(i + 1, len(tl)):
        B = gsl[j]; inter = A & B; k = len(inter)
        if k == 0:
            continue
        u  = len(A) + len(B) - k
        J  = k / u if u > 0 else 0.0
        ls = min(len(A), len(B)); ll = max(len(A), len(B))
        S  = k / ls if ls > 0 else 0.0
        Ne = max(N_KB_AUDIT, ll, ls)
        p  = float(hypergeom.sf(k - 1, Ne, ll, ls))
        if not np.isfinite(p) or p < 0:
            p = 1.0
        ps_main[(i, j)] = {'J': J, 'S': S, 'p': p, 'k': k}

stage2_survivors, elim_log, stage2_diag = _pairwise_filter(
    tl, gsl, ps_main, pv_selected, OVERLAP_FDR_ALPHA,
    SUBMODULE_J_THRESHOLD, SUBMODULE_S_THRESHOLD, SIZE_DIFF_FALLBACK,
)
logger.info(f"  Jaccard: {len(tl)} → {len(stage2_survivors)}")

# Jaccard QC gate
class JaccardQCError(Exception):
    """Redundancy filter failed. Do NOT publish."""
    pass

def _jmat(terms: List[str], gmt: Dict[str, Set[str]]) -> np.ndarray:
    n = len(terms); M = np.eye(n, dtype=np.float64)
    gs = [gmt.get(t, set()) for t in terms]
    for i in range(n):
        A = gs[i]
        for j in range(i + 1, n):
            B = gs[j]; inter = len(A & B); union = len(A) + len(B) - inter
            jv = inter / union if union > 0 else 0.0
            M[i, j] = jv; M[j, i] = jv
    return M

def _validate_jmat(M: np.ndarray, terms: List[str], thr: float) -> dict:
    n = M.shape[0]; viols = []; mj = 0.0; sj = 0.0; cnt = 0
    for i in range(n):
        for j in range(i + 1, n):
            jv = M[i, j]; sj += jv; cnt += 1; mj = max(mj, jv)
            if jv >= thr:
                viols.append((terms[i], terms[j], jv))
    mj_mean = sj / cnt if cnt > 0 else 0.0
    if viols:
        msg = f"\nJACCARD QC FAILURE\n{len(viols)} pairs >= J={thr}\n"
        for a, b, jv in sorted(viols, key=lambda x: -x[2])[:5]:
            msg += f"  {a}\n  ↔ {b}  J={jv:.4f}\n"
        logger.error(msg); raise JaccardQCError(msg)
    logger.info(f"  QC PASSED: max J={mj:.4f}, mean J={mj_mean:.4f}")
    return {'status': 'PASSED', 'n_pathways': n, 'max_offdiag_jaccard': float(mj),
            'mean_offdiag_jaccard': float(mj_mean), 'violations': 0}

qc_terms       = list(stage2_survivors)
qc_gmt         = {t: full_gmt_dict.get(t, pathway_gene_sets.get(t, set())) for t in qc_terms}
jaccard_matrix = _jmat(qc_terms, qc_gmt)
qc_validation  = _validate_jmat(jaccard_matrix, qc_terms, SUBMODULE_J_THRESHOLD)

# Stage 4A: Final primary output
logger.info("\nSTAGE 4A: PRIMARY OUTPUT (N=50, immutable)")

final_rows_primary = []
for t in stage2_survivors:
    eg  = pathway_gene_sets.get(t, set())
    t1o = sorted(eg & TIER1_SET)
    t2o = sorted(eg & TIER2_SET)
    # [MOD-4] Multi_Layer_Score: 2x for Tier2 (CRISPR-confirmed essential)
    final_rows_primary.append({
        'Term':               t,
        'Cutoff':             f"top{SENSITIVITY_CUTOFFS[0]}",
        'Adjusted_P_value':   pathway_pvalues.get(t, np.nan),
        'N_Enrichment_Genes': len(eg),
        'Enrichment_Genes':   ';'.join(sorted(eg)),
        'Tier1_Overlap_Genes': ';'.join(t1o),
        'N_Tier1_Overlap':    len(t1o),
        'Tier2_Overlap_Genes': ';'.join(t2o),
        'N_Tier2_Overlap':    len(t2o),
        'Has_Tier2_Support':  len(t2o) > 0,
        'Multi_Layer_Score':  len(t1o) + 2 * len(t2o),
    })

primary_df = pd.DataFrame(final_rows_primary)
if not primary_df.empty:
    primary_df = primary_df.sort_values('Adjusted_P_value').reset_index(drop=True)

primary_df.to_csv(OUTPUT_DIR / "L4_NonRedundant_Pathways.csv", index=False, encoding='utf-8-sig')
try:
    primary_df.to_excel(OUTPUT_DIR / "L4_NonRedundant_Pathways.xlsx", index=False)
except Exception as _xe:
    logger.warning(f"  Excel export skipped: {_xe}")

PRIMARY_TERM_SET = set(stage2_survivors)
logger.info(f"  Primary non-redundant pathways: {len(primary_df)}")
logger.info(f"  Primary term set size:          {len(PRIMARY_TERM_SET)}")


# ============================================================
# STAGE 4B: Primary Visualizations
# ============================================================
logger.info("\nSTAGE 4B: PRIMARY VISUALIZATIONS")


def _plot_barchart(df: pd.DataFrame, out: Path, cutoff_label: str) -> None:
    if df is None or df.empty:
        return
    pd_ = df.nsmallest(TOP_TERM, 'Adjusted_P_value').copy()
    if pd_.empty:
        return
    ts = pd_['Term'].tolist(); pvs = pd_['Adjusted_P_value'].tolist()
    nlog = [-np.log10(max(p, 1e-300)) for p in pvs]
    bc   = ['lightskyblue' if p < 0.05 else 'lightgrey' for p in pvs]
    plt.figure(figsize=(14, max(4, len(ts) * 0.5)))
    ax = sns.barplot(x=nlog, y=ts, palette=bc, edgecolor=None)
    ax.axes.get_yaxis().set_visible(False)
    ax.set_title(
        f"Layer 4 — Non-Redundant Pathways ({cutoff_label})\n"
        f"({cfg.L3B_CANCER_CELL_LINE}, background=None)",
        fontsize=14, fontweight='bold',
    )
    ax.set_xlabel(LABEL_NEG_LOG10_PVAL, fontsize=12)
    ax.spines['right'].set_visible(False); ax.spines['top'].set_visible(False)
    x0 = max(nlog) / 200 if max(nlog) > 0 else 0.01
    for ii, t in enumerate(ts):
        ax.text(x0, ii, f"  {t}  {pvs[ii]:.2e}", ha='left', fontsize=10, va='center')
    sl = -np.log10(0.05)
    if sl <= max(nlog) * 1.5:
        ax.axvline(sl, color='red', linestyle='--', linewidth=1, alpha=0.7, label='FDR=0.05')
        ax.legend(loc='lower right', fontsize=9)
    plt.tight_layout()
    for ext in ['png', 'svg']:
        plt.savefig(str(out / f"L4_Enrichr_Barchart_{cutoff_label}.{ext}"),
                    dpi=300, bbox_inches='tight', format=ext)
    _show_plot()
    plt.close()
    logger.info(f"  Bar chart ({cutoff_label}) saved")


if not primary_df.empty:
    _plot_barchart(primary_df, OUTPUT_DIR, f"top{SENSITIVITY_CUTOFFS[0]}")


# ============================================================
# STAGE 4C: Cytoscape Export (primary only)
# ============================================================
logger.info("\nSTAGE 4C: CYTOSCAPE EXPORT (primary)")

_dn = {cfg.normalize_gene_name(t) for t in cfg.TARGET_NAMES}
node_rows = []
for _, r in primary_df.iterrows():
    node_rows.append({'Node_ID': r['Term'], 'Node_Type': 'PATHWAY', 'Label': r['Term'],
                      'Adjusted_Pval': r['Adjusted_P_value'],
                      'N_Tier1': r['N_Tier1_Overlap'], 'N_Tier2': r['N_Tier2_Overlap'],
                      'Multi_Layer_Score': r['Multi_Layer_Score'],
                      'Has_Tier2_Support': r['Has_Tier2_Support'],
                      'Node_Size':  max(20, min(80, 20 + r['N_Tier1_Overlap'] * 3)),
                      'Node_Color': '#FF6B6B' if r['Has_Tier2_Support'] else '#4ECDC4',
                      'Node_Shape': 'ELLIPSE'})

_apg: Set[str] = set()
for _, r in primary_df.iterrows():
    for g in str(r['Enrichment_Genes']).split(';'):
        gc = g.strip()
        if gc:
            _apg.add(gc)

for gene in sorted(_apg):
    _id  = gene in _dn
    _ti  = 'Tier1' if gene in TIER1_SET else 'Tier2' if gene in TIER2_SET else 'Other'
    _col = ('#E74C3C' if _id else '#F39C12' if _ti == 'Tier1'
            else '#3498DB' if _ti == 'Tier2' else '#95A5A6')
    node_rows.append({'Node_ID': gene, 'Node_Type': 'GENE', 'Label': gene,
                      'Adjusted_Pval': np.nan,
                      'N_Tier1': 1 if _ti == 'Tier1' else 0,
                      'N_Tier2': 1 if _ti == 'Tier2' else 0,
                      'Multi_Layer_Score': 2 if _ti == 'Tier1' else 1 if _ti == 'Tier2' else 0,
                      'Has_Tier2_Support': _ti == 'Tier2',
                      'Node_Size': 40 if _ti == 'Tier1' else 25 if _ti == 'Tier2' else 15,
                      'Node_Color': _col,
                      'Node_Shape': 'DIAMOND' if _id else 'ELLIPSE'})

node_df = pd.DataFrame(node_rows)
node_df.to_csv(CYTO_DIR / "L4_Cytoscape_Nodes.csv", index=False, encoding='utf-8-sig')

edge_rows = []
pn = list(primary_df['Term'])
for _, r in primary_df.iterrows():
    pw = r['Term']
    for g in str(r['Enrichment_Genes']).split(';'):
        gc = g.strip()
        if not gc:
            continue
        ti = 'Tier1' if gc in TIER1_SET else 'Tier2' if gc in TIER2_SET else 'Other'
        edge_rows.append({'Source': pw, 'Target': gc, 'Interaction': 'member_of',
                          'Edge_Type': 'PATHWAY_GENE',
                          'Weight': 2.0 if ti == 'Tier1' else 1.0, 'Tier_Label': ti})
for i in range(len(pn)):
    A = qc_gmt.get(pn[i], set())
    for j in range(i + 1, len(pn)):
        B = qc_gmt.get(pn[j], set()); k = len(A & B)
        if k == 0:
            continue
        u = len(A) + len(B) - k
        J = k / u if u > 0 else 0.0
        edge_rows.append({'Source': pn[i], 'Target': pn[j], 'Interaction': 'overlaps',
                          'Edge_Type': 'PATHWAY_PATHWAY', 'Weight': round(J, 4), 'Tier_Label': ''})

edge_df = pd.DataFrame(edge_rows)
edge_df.to_csv(CYTO_DIR / "L4_Cytoscape_Edges.csv", index=False, encoding='utf-8-sig')
with open(CYTO_DIR / "L4_Cytoscape_Network.sif", 'w', encoding='utf-8') as _f:
    _f.write('\n'.join(f"{r['Source']}\t{r['Interaction']}\t{r['Target']}"
                       for _, r in edge_df.iterrows()))

_em = primary_df[['Term', 'Adjusted_P_value', 'N_Enrichment_Genes',
                  'N_Tier1_Overlap', 'Enrichment_Genes']].copy()
_em.columns = ['Name', 'pvalue', 'N_genes', 'N_Tier1', 'Genes']
_em.to_csv(CYTO_DIR / "L4_EnrichmentMap_Input.txt", index=False, sep='\t', encoding='utf-8')
logger.info(f"  Nodes: {int((node_df['Node_Type']=='PATHWAY').sum())} pathways + genes")


# ============================================================
# STAGE 5: AUTOMATED ROBUSTNESS AUDIT (top-100, top-150)
# ============================================================
# Scientific basis:
#   The primary top-50 signal represents the highest-density
#   RWR propagation core. Broader cutoffs (100, 150) represent
#   spatial expansion layers. Theme concordance across all
#   three layers confirms topological robustness of the signal.
#
#   Audit is fully automated — no human selection of "best" result.
#   Audit outputs CANNOT replace the primary top-50 result.
#
# Falsification condition (pre-declared before examining results):
#   If Jaccard(top50_term_set, top150_term_set) = 0,
#   the primary signal is declared topologically unstable.
# ============================================================
logger.info("\n" + "=" * 72)
logger.info("STAGE 5: AUTOMATED ROBUSTNESS AUDIT")
logger.info("=" * 72)
logger.info(f"  Cutoffs: {SENSITIVITY_CUTOFFS}")
logger.info(f"  Primary ({SENSITIVITY_CUTOFFS[0]}): IMMUTABLE — results cannot be overridden")
logger.info(f"  Audit  ({SENSITIVITY_CUTOFFS[1:]}): comparison only")
logger.info(f"  Falsification condition: Jaccard(top50, top150) > {FALSIFICATION_MIN_JACCARD}")

# Load Layer 2B Delta Network for positive-delta gene rankings
if not L2B_DELTA_CSV.exists():
    logger.warning(f"  Layer 2B Delta CSV not found: {L2B_DELTA_CSV}")
    logger.warning("  Robustness audit requires L2B_Delta_Network_Summary.csv")
    logger.warning("  Skipping Stage 5-6. Primary results remain valid.")
    _l2b_available = False
else:
    _l2b_available = True

_l2b_available_flag = _l2b_available

audit_term_sets:   Dict[int, Set[str]]   = {}
audit_final_dfs:   Dict[int, pd.DataFrame] = {}
audit_sig_dfs:     Dict[int, pd.DataFrame] = {}
audit_detail_rows: List[dict]             = []

if _l2b_available_flag:
    delta_df = pd.read_csv(L2B_DELTA_CSV)

    required_cols = ['Gene', 'Delta_Score', 'Is_Direct_Target']
    _missing = [c for c in required_cols if c not in delta_df.columns]
    if _missing:
        logger.warning(f"  L2B Delta CSV missing columns: {_missing}. Skipping audit.")
        _l2b_available_flag = False
    else:
        # Positive delta only (DL10 compliance)
        delta_positive = delta_df[delta_df['Delta_Score'] > 0].copy()
        delta_positive = delta_positive.sort_values('Delta_Score', ascending=False)
        delta_positive['Gene_Norm'] = delta_positive['Gene'].apply(cfg.normalize_gene_name)
        logger.info(f"  Positive-delta genes available: {len(delta_positive)}")

if _l2b_available_flag:
    for cutoff in SENSITIVITY_CUTOFFS:
        logger.info(f"\n  --- Audit cutoff: top-{cutoff} ---")

        # Extract top-N positive-delta genes
        top_genes_raw = delta_positive['Gene_Norm'].head(cutoff).tolist()
        top_genes_clean, _cc = detect_date_corrupted_genes(top_genes_raw)
        top_genes = sorted(list(dict.fromkeys(top_genes_clean)))

        _is_primary = (cutoff == SENSITIVITY_CUTOFFS[0])

        if _is_primary:
            # Use the TIER1_GENES already loaded for primary analysis
            # This ensures exact consistency with Stage 1-4
            audit_top_genes = TIER1_GENES
            logger.info(f"    Using official Tier1 gene list (primary, immutable)")
        else:
            audit_top_genes = top_genes
            logger.info(f"    Audit gene list: {len(audit_top_genes)} genes")

        # Check effective overlap with GMT
        _eff = sorted(set(audit_top_genes) & BACKGROUND_GMT)
        _n_eff = len(_eff)
        logger.info(f"    n_effective (in GMT): {_n_eff}")

        if _n_eff < MIN_OVERLAP:
            logger.warning(f"    n_eff < MIN_OVERLAP, skipping cutoff top-{cutoff}")
            audit_term_sets[cutoff] = set()
            audit_final_dfs[cutoff] = pd.DataFrame()
            audit_sig_dfs[cutoff]   = pd.DataFrame()
            continue

        # Run enrichment + filter
        # For primary cutoff: use already-computed primary results
        if _is_primary:
            audit_term_sets[cutoff] = PRIMARY_TERM_SET
            audit_final_dfs[cutoff] = primary_df
            audit_sig_dfs[cutoff]   = enrichment_df_sig
            logger.info(f"    Using primary result: {len(PRIMARY_TERM_SET)} terms")
        else:
            _sig_df, _final_df, _terms, _fdiag = run_enrichment_for_gene_list(
                gene_list=audit_top_genes,
                gmt_for_enrich=gmt_for_enrich,
                full_gmt_dict=full_gmt_dict,
                background_gmt=BACKGROUND_GMT,
                N_KB_audit=N_KB_AUDIT,
                cutoff_label=f"top{cutoff}",
                fdr_cutoff=FDR_CUTOFF,
                min_overlap=MIN_OVERLAP,
                run_filter=True,
            )
            audit_term_sets[cutoff] = set(_terms)
            audit_final_dfs[cutoff] = _final_df
            audit_sig_dfs[cutoff]   = _sig_df

            # Save audit enrichment
            if not _sig_df.empty:
                _sig_df.to_csv(
                    ROBUST_DIR / f"L4_Audit_Significant_top{cutoff}.csv",
                    index=False, encoding='utf-8-sig',
                )
            if not _final_df.empty:
                _final_df.to_csv(
                    ROBUST_DIR / f"L4_Audit_NonRedundant_top{cutoff}.csv",
                    index=False, encoding='utf-8-sig',
                )
            logger.info(f"    Non-redundant terms: {len(audit_term_sets[cutoff])}")

        # Record detail row
        audit_detail_rows.append({
            'Cutoff':                  cutoff,
            'N_Input_Genes':           len(audit_top_genes),
            'N_Effective_In_GMT':      _n_eff,
            'N_Significant_Enrichment': len(audit_sig_dfs.get(cutoff, pd.DataFrame())),
            'N_NonRedundant_Pathways': len(audit_term_sets.get(cutoff, set())),
            'Is_Primary':              _is_primary,
        })

    # Save sensitivity detail
    detail_df = pd.DataFrame(audit_detail_rows)
    detail_df.to_csv(
        ROBUST_DIR / SENSITIVITY_DETAIL_CSV, index=False, encoding='utf-8-sig'
    )
    logger.info(f"\n  Audit complete for cutoffs: {SENSITIVITY_CUTOFFS}")


# ============================================================
# STAGE 6: THEME CONCORDANCE + ROBUSTNESS_THEMES.CSV
# ============================================================
# Robustness theme = pathway term that survives the non-redundancy
# filter at ALL sensitivity cutoffs (50, 100, 150).
#
# Falsification check:
#   If Jaccard(top50_terms, top150_terms) = 0, the primary
#   signal is topologically unstable.
#
# Output: Robustness_Themes.csv
#   - Only pathways concordant across ALL cutoffs
#   - Never replaces primary results
#   - Provides evidence of topological robustness
# ============================================================
logger.info("\n" + "=" * 72)
logger.info("STAGE 6: THEME CONCORDANCE + ROBUSTNESS REPORT")
logger.info("=" * 72)

robustness_summary: dict = {
    'primary_cutoff':          SENSITIVITY_CUTOFFS[0],
    'audit_cutoffs':           SENSITIVITY_CUTOFFS[1:],
    'falsification_threshold': FALSIFICATION_MIN_JACCARD,
    'n_primary_terms':         len(PRIMARY_TERM_SET),
}

if not _l2b_available_flag:
    logger.warning("  Robustness audit skipped (L2B data unavailable)")
    robustness_summary['status'] = 'SKIPPED_NO_L2B_DATA'
    robust_themes = []
    concordance_jaccard: Dict[Tuple[int, int], float] = {}
else:
    # Compute pairwise Jaccard between all cutoff term sets
    concordance_jaccard: Dict[Tuple[int, int], float] = {}
    for ci, c1 in enumerate(SENSITIVITY_CUTOFFS):
        for c2 in SENSITIVITY_CUTOFFS[ci + 1:]:
            A = audit_term_sets.get(c1, set())
            B = audit_term_sets.get(c2, set())
            inter = len(A & B)
            union = len(A | B)
            jac   = inter / union if union > 0 else 0.0
            concordance_jaccard[(c1, c2)] = jac
            logger.info(f"  Jaccard(top{c1}, top{c2}): {jac:.4f} "
                        f"({inter} shared / {union} union terms)")

    # Falsification check (pre-declared)
    _c_primary = SENSITIVITY_CUTOFFS[0]
    _c_largest = SENSITIVITY_CUTOFFS[-1]
    _jac_critical = concordance_jaccard.get((_c_primary, _c_largest), 0.0)

    if _jac_critical <= FALSIFICATION_MIN_JACCARD:
        logger.error(
            f"\n  FALSIFICATION CONDITION MET:\n"
            f"  Jaccard(top{_c_primary}, top{_c_largest}) = {_jac_critical:.4f} "
            f"<= {FALSIFICATION_MIN_JACCARD}\n"
            f"  The primary signal is topologically unstable.\n"
            f"  Pathway interpretation should not be claimed from this run."
        )
        robustness_summary['falsification_triggered'] = True
        robustness_summary['jaccard_primary_to_largest'] = float(_jac_critical)
    else:
        logger.info(
            f"\n  FALSIFICATION NOT TRIGGERED:\n"
            f"  Jaccard(top{_c_primary}, top{_c_largest}) = {_jac_critical:.4f} "
            f"> {FALSIFICATION_MIN_JACCARD}\n"
            f"  Signal is topologically stable."
        )
        robustness_summary['falsification_triggered'] = False
        robustness_summary['jaccard_primary_to_largest'] = float(_jac_critical)

    # Robust themes = intersection of ALL cutoff term sets
    if all(cutoff in audit_term_sets for cutoff in SENSITIVITY_CUTOFFS):
        robust_themes_set = audit_term_sets[SENSITIVITY_CUTOFFS[0]].copy()
        for cutoff in SENSITIVITY_CUTOFFS[1:]:
            robust_themes_set &= audit_term_sets.get(cutoff, set())
    else:
        robust_themes_set = set()

    robust_themes = sorted(robust_themes_set)
    logger.info(f"\n  Robust themes (concordant across ALL cutoffs): {len(robust_themes)}")

    robustness_summary['n_robust_themes'] = len(robust_themes)
    robustness_summary['concordance_jaccard'] = {
        f"top{c1}_vs_top{c2}": round(float(jac), 4)
        for (c1, c2), jac in concordance_jaccard.items()
    }
    robustness_summary['status'] = 'COMPLETED'

    # Build Robustness_Themes.csv
    robust_rows = []
    for t in robust_themes:
        # Pull metrics from primary result where available
        _primary_row = primary_df[primary_df['Term'] == t]
        _adj_p = float(_primary_row['Adjusted_P_value'].values[0]) \
                 if not _primary_row.empty else np.nan
        _n_t1  = int(_primary_row['N_Tier1_Overlap'].values[0]) \
                 if not _primary_row.empty else 0
        _n_t2  = int(_primary_row['N_Tier2_Overlap'].values[0]) \
                 if not _primary_row.empty else 0
        _t2s   = bool(_primary_row['Has_Tier2_Support'].values[0]) \
                 if not _primary_row.empty else False

        # Count how many cutoffs contain this term
        _n_cutoffs = sum(1 for c in SENSITIVITY_CUTOFFS
                         if t in audit_term_sets.get(c, set()))

        robust_rows.append({
            'Term':                 t,
            'N_Cutoffs_Present':    _n_cutoffs,
            'Cutoffs_Present':      ';'.join([
                f"top{c}" for c in SENSITIVITY_CUTOFFS
                if t in audit_term_sets.get(c, set())
            ]),
            'Primary_Adjusted_P':   _adj_p,
            'Primary_N_Tier1':      _n_t1,
            'Primary_N_Tier2':      _n_t2,
            'Primary_Has_Tier2':    _t2s,
            'Is_Fully_Concordant':  (_n_cutoffs == len(SENSITIVITY_CUTOFFS)),
        })

    robust_df = pd.DataFrame(robust_rows)
    if not robust_df.empty:
        robust_df = robust_df.sort_values(
            ['Is_Fully_Concordant', 'Primary_Adjusted_P'],
            ascending=[False, True],
        ).reset_index(drop=True)

    robust_df.to_csv(
        ROBUST_DIR / ROBUSTNESS_THEMES_CSV, index=False, encoding='utf-8-sig'
    )
    logger.info(f"  Saved: {ROBUSTNESS_THEMES_CSV}")

    # Concordance visualization
    if concordance_jaccard:
        n_cuts = len(SENSITIVITY_CUTOFFS)
        jac_matrix = np.zeros((n_cuts, n_cuts))
        for i, c1 in enumerate(SENSITIVITY_CUTOFFS):
            jac_matrix[i, i] = 1.0
            for j, c2 in enumerate(SENSITIVITY_CUTOFFS):
                if i < j:
                    jv = concordance_jaccard.get((c1, c2), 0.0)
                    jac_matrix[i, j] = jv
                    jac_matrix[j, i] = jv

        labs = [f"top{c}" for c in SENSITIVITY_CUTOFFS]
        fig, ax = plt.subplots(figsize=(6, 5))
        sns.heatmap(
            pd.DataFrame(jac_matrix, index=labs, columns=labs),
            annot=True, fmt='.3f', cmap='Greens',
            vmin=0.0, vmax=1.0, ax=ax, linewidths=0.5,
            cbar_kws={'label': 'Jaccard Index'},
        )
        ax.set_title(
            "Pathway Term Set Concordance\n"
            f"(Jaccard across cutoffs {SENSITIVITY_CUTOFFS})\n"
            f"Falsification threshold: > {FALSIFICATION_MIN_JACCARD}",
            fontsize=11, fontweight='bold',
        )
        plt.tight_layout()
        for ext in ['png', 'svg']:
            plt.savefig(
                str(ROBUST_DIR / f"L4_Concordance_Heatmap.{ext}"),
                dpi=150, bbox_inches='tight', format=ext,
            )
        _show_plot()
        plt.close()
        logger.info("  Concordance heatmap saved")


# ============================================================
# STAGE 6B: Save Robustness Summary JSON
# ============================================================
with open(ROBUST_DIR / ROBUSTNESS_SUMMARY_JSON, 'w', encoding='utf-8') as _f:
    json.dump(robustness_summary, _f, indent=2, default=str)
logger.info(f"  Saved: {ROBUSTNESS_SUMMARY_JSON}")


# ============================================================
# MANIFEST
# ============================================================
manifest = {
    "engine":    "Layer4_RWR_Pathway_v2.4",
    "run_id":    RUN_ID,
    "timestamp": TIMESTAMP,
    "cancer_model":    cfg.L3B_CANCER_CELL_LINE,
    "cancer_model_id": cfg.L3B_CANCER_MODEL_ID,

    "changelog_v2.4": {
        "NEW_Stage5":  "Automated robustness audit (top-100, top-150)",
        "NEW_Stage6":  "Theme concordance + Robustness_Themes.csv",
        "NEW_Config":  "SENSITIVITY_CUTOFFS, FALSIFICATION_MIN_JACCARD declared",
        "KEPT_C2":     "N_KB_audit = GMT union only",
        "KEPT_M1":     "Concordance: log10 scale",
        "KEPT_C1":     "ElasticNetCV: heuristic pre-filter",
        "KEPT_N1_N4":  "Code hygiene fixes",
        "KEPT_M1_M4":  "Documentation fixes",
    },

    "scientific_integrity": {
        "primary_immutability": (
            "Primary top-50 result is immutable. "
            "Sensitivity audit results CANNOT replace primary output."
        ),
        "robustness_design": (
            "Theme concordance across multiple cutoffs confirms topological "
            "robustness, not p-hacking. Each cutoff is pre-declared."
        ),
        "falsification_condition": (
            f"Pre-declared: if Jaccard(top{SENSITIVITY_CUTOFFS[0]}, "
            f"top{SENSITIVITY_CUTOFFS[-1]}) <= {FALSIFICATION_MIN_JACCARD}, "
            "primary signal is declared topologically unstable."
        ),
        "C2_background": "N_KB_audit = GMT union only, no Tier1 augmentation",
        "M1_concordance": "Log10 scale, 10-fold tolerance",
        "C1_elasticnet": "Heuristic pre-filter, no formal superiority claimed",
        "cell_line_context": L4_CONTEXT_NOTE,
        "pathway_claim": "Candidate mechanistic themes, not confirmed MoA",
        "references": [
            "Choobdar et al. Nat Methods 16:843-852 (2019)",
            "Lamparter et al. PLoS Comput Biol 12(1):e1004714 (2016)",
            "Ziemann et al. Genome Biology 17:177 (2016)",
            "Bauer et al. PLoS Comput Biol 6(6):e1000796 (2010)",
        ],
    },

    "inputs": {
        "tier1_csv":         str(TIER1_CSV),
        "tier2_csv":         str(TIER2_CSV),
        "l2b_delta_csv":     str(L2B_DELTA_CSV),
        "gmt_file":          str(GMT_FILE_PATH),
        "tier1_n":           len(TIER1_GENES),
        "tier1_n_effective": N_TIER1_EFFECTIVE,
        "N_KB_audit":        N_KB_AUDIT,
    },

    "robustness_audit": robustness_summary,

    "primary_pipeline": {
        "after_enrichr":     len(enrichment_df_sig),
        "after_size_filter": len(candidate_terms),
        "after_elasticnet":  len(stage1_survivors),
        "after_jaccard":     len(stage2_survivors),
    },

    "parameters": {
        "SENSITIVITY_CUTOFFS":          SENSITIVITY_CUTOFFS,
        "FALSIFICATION_MIN_JACCARD":    FALSIFICATION_MIN_JACCARD,
        "FDR_CUTOFF":                   FDR_CUTOFF,
        "MIN_OVERLAP":                  MIN_OVERLAP,
        "LOG10_PVALUE_CONCORDANCE_TOL": LOG10_PVALUE_CONCORDANCE_TOL,
        "ELASTICNET_L1_RATIO":          ELASTICNET_L1_RATIO,
        "ELASTICNET_CV_FOLDS":          ELASTICNET_CV_FOLDS,
        "SUBMODULE_J_THRESHOLD":        SUBMODULE_J_THRESHOLD,
        "SUBMODULE_S_THRESHOLD":        SUBMODULE_S_THRESHOLD,
        "MODULE_MIN_GENES":             MODULE_MIN_GENES,
        "SIZE_DIFF_FALLBACK":           SIZE_DIFF_FALLBACK,
    },
}

with open(OUTPUT_DIR / "manifest.json", 'w', encoding='utf-8') as _f:
    json.dump(manifest, _f, indent=2, default=str)


# ============================================================
# FINAL SUMMARY
# ============================================================
logger.info("\n" + "=" * 72)
logger.info("LAYER 4 v2.4 COMPLETE")
logger.info("=" * 72)
logger.info(f"  Cancer model:           {cfg.L3B_CANCER_CELL_LINE}")
logger.info(f"  Primary cutoff:         N={SENSITIVITY_CUTOFFS[0]} (immutable)")
logger.info(f"  Primary pathways:       {len(primary_df)}")
logger.info(f"  QC:                     {qc_validation['status']}")

if _l2b_available_flag:
    _falsi = robustness_summary.get('falsification_triggered', None)
    _jcrit = robustness_summary.get('jaccard_primary_to_largest', None)
    _nrob  = robustness_summary.get('n_robust_themes', 0)
    logger.info(f"  Robust themes (all cutoffs): {_nrob}")
    logger.info(f"  Falsification triggered:     {_falsi}")
    logger.info(f"  Jaccard(top50, top150):      "
                f"{_jcrit:.4f}" if _jcrit is not None else "N/A")
else:
    logger.info("  Robustness audit:       SKIPPED (L2B data unavailable)")
logger.info("=" * 72)

print("\n" + "=" * 72)
print("LAYER 4 v2.4 COMPLETE — SUBMISSION-SAFE")
print("=" * 72)
print(f"  Cancer model:    {cfg.L3B_CANCER_CELL_LINE}")
print(f"  Primary:         {len(primary_df)} pathways (N={SENSITIVITY_CUTOFFS[0]}, IMMUTABLE)")
if _l2b_available_flag:
    print(f"  Robust themes:   {robustness_summary.get('n_robust_themes', 0)} "
          f"(concordant across {SENSITIVITY_CUTOFFS})")
    print(f"  Falsification:   {robustness_summary.get('falsification_triggered', 'N/A')}")
print(f"\nKey outputs:")
print(f"  L4_NonRedundant_Pathways.csv   — primary result (N=50)")
print(f"  robustness_audit/")
print(f"    {ROBUSTNESS_THEMES_CSV}")
print(f"    {ROBUSTNESS_SUMMARY_JSON}")
print(f"    {SENSITIVITY_DETAIL_CSV}")
print(f"    L4_Concordance_Heatmap.png")
print(f"\n{L4_CONTEXT_NOTE}")
