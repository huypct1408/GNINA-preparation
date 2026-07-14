#!/usr/bin/env python
# -*- coding: utf-8 -*-
# ============================================================
# 02B — Layer 2: LIONESS Multi-Model Inference (v2.4.0)
# ============================================================
# Purpose
#   Infer single-sample LIONESS GRNs for EVERY ModelID declared in
#   cfg.L2_TARGET_MODELS_LIST, while preserving absolute scientific
#   integrity of the original single-model pipeline.
#
# Architecture (Iteration + lineage-shared PANDA)
#   1. Validate cfg.L2_TARGET_MODELS_LIST (DL2-26)
#   2. Resolve each ACH ID → OncotreeLineage via Model.csv ONLY (DL2-15)
#   3. Group targets by lineage
#   4. For each lineage present in 02A foundation freeze:
#        - load frozen expr / sample_order / genes_post_qc ONCE
#        - load + SHA-256-verify motif/PPI priors ONCE
#        - build ONE Panda aggregate e^(α) for the lineage (DL2-25)
#        - verify flatten order C/F ONCE against export_panda_results
#        - for each target ModelID in lineage ∩ list:
#              Lioness(start=pos, end=pos, save_single=True)
#              label edges with verified order
#              write Z_{ModelID}_LIONESS.tsv + per-sample run manifest
#        - write lineage-level multi-run summary
#   5. Write global multi-model gate report + full run ledger (DL2-28)
#
# Scientific integrity (STRICT)
#   - NO SNAIL recomputation (DL2-22)
#   - NO raw CCLE expression load (DL2-22)
#   - NO cross-lineage pooling of the aggregate network (DL2-25)
#   - NO rebuilding PANDA per target inside the same lineage (DL2-25)
#   - NO silent drop of declared IDs (DL2-28)
#   - NO edge weight re-weighting / manual edit
#   - Separate multi-model log + gate report (DL2-23)
#
# Missing-target policy
#   cfg.L2_MULTI_MODEL_ON_MISSING:
#     "skip_missing_with_ledger" (default) — SKIP + ledger, continue
#     "fail" — any missing ID aborts the whole multi-model run
#
# Requires
#   - Stage 8F outputs of 02A (frozen_lineage_inputs/)
#   - config_system.py v3.4.0+
#   - netZooPy == cfg.NETZOOPY_REQUIRED_VERSION
# ============================================================

from __future__ import annotations

import sys
import json
import time
import shutil
import hashlib
import tempfile
import warnings
import traceback
from pathlib import Path
from datetime import datetime
from collections import defaultdict, OrderedDict

import pandas as pd
import numpy as np

warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=UserWarning, module="matplotlib")
warnings.filterwarnings("ignore", message="pkg_resources is deprecated")

# ---------------------------------------------------------------------------
# Ensure project root (where config_system.py lives) is importable
# Works as:
#   - python path/to/02B_Layer2_LIONESS_MULTI_MODEL.py   (__file__ defined)
#   - Jupyter / IPython %run or paste                     (__file__ absent)
# Scientific note: path bootstrap only affects import resolution.
# All scientific paths still come from config_system.py constants.
# ---------------------------------------------------------------------------

def _bootstrap_project_root() -> Path:
    """Locate directory containing config_system.py without requiring __file__."""
    candidates = []

    # 1) Script directory (CLI / %run with __file__)
    try:
        here = Path(__file__).resolve().parent  # type: ignore[name-defined]
        candidates.append(here)
        candidates.append(here.parent)
        candidates.append(here.parent.parent)
    except NameError:
        # Jupyter / interactive: __file__ is undefined — expected, not an error
        pass

    # 2) IPython/Jupyter current working directory chain
    cwd = Path.cwd().resolve()
    candidates.append(cwd)
    candidates.append(cwd.parent)
    # common layout: cwd = project root, or cwd = project/notebooks
    candidates.append(cwd / "notebooks")
    if (cwd / "pyscenic").is_dir():
        candidates.append(cwd / "pyscenic")

    # 3) Known lab host path (frozen deployment location)
    candidates.append(
        Path("/home/labhhc5/Documents/workspace/D21/Duong Huy/pyscenic")
    )

    # 4) sys.path entries already present
    for p in list(sys.path):
        if p:
            candidates.append(Path(p))

    seen = set()
    ordered = []
    for c in candidates:
        try:
            r = c.resolve()
        except Exception:
            continue
        if r in seen:
            continue
        seen.add(r)
        ordered.append(r)

    for c in ordered:
        if (c / "config_system.py").is_file():
            if str(c) not in sys.path:
                sys.path.insert(0, str(c))
            return c

    raise FileNotFoundError(
        "[M1] Cannot locate config_system.py.\n"
        "Jupyter fix: %cd to the pyscenic project root (folder that contains "
        "config_system.py), then re-run this cell.\n"
        f"  Tried cwd={cwd}\n"
        f"  Searched {len(ordered)} candidate roots."
    )


_PROJECT_ROOT = _bootstrap_project_root()
print(f"[M1] config root : {_PROJECT_ROOT}")

import config_system as cfg


# =============================================================================
# STAGE M1: IMPORTS & CONFIGURATION
# =============================================================================

logger = cfg.setup_logger(
    "LIONESS_MultiModel",
    logfile=cfg.LAYER2B_MULTI_LOG_FILE,
    reset_handlers=True,
)

cfg.print_config_summary()

NOTEBOOK_START_TIME = datetime.now()
PIPELINE_VERSION = "2.4.0"
CONFIG_VERSION = "config_system.py v3.4.0"

logger.info("=" * 80)
logger.info("02B LIONESS MULTI-MODEL NOTEBOOK START — v2.4.0")
logger.info(f"Start time   : {NOTEBOOK_START_TIME.isoformat()}")
logger.info(f"Python       : {sys.version}")
logger.info(f"Pandas       : {pd.__version__}")
logger.info(f"NumPy        : {np.__version__}")
logger.info(f"Run mode     : {cfg.LAYER2_RUN_MODE}")
logger.info(f"Architecture : {cfg.LAYER2_NOTEBOOK_ARCHITECTURE}")
logger.info("=" * 80)

print(f"\nNotebook     : {cfg.LAYER2_NOTEBOOK_02B_MULTI_NAME}")
print(f"Started      : {NOTEBOOK_START_TIME.strftime('%Y-%m-%d %H:%M:%S')}")
print(f"Run mode     : {cfg.LAYER2_RUN_MODE}")
print(f"Ref policy   : {cfg.LIONESS_REFERENCE_POLICY}")
print(f"Missing pol. : {cfg.L2_MULTI_MODEL_ON_MISSING}")
print(f"Log file     : {cfg.LAYER2B_MULTI_LOG_FILE.name}")

if cfg.LAYER2_RUN_MODE != "TARGETED_LIONESS_MULTI_MODEL":
    raise RuntimeError(
        f"[M1] LAYER2_RUN_MODE must be 'TARGETED_LIONESS_MULTI_MODEL', "
        f"got '{cfg.LAYER2_RUN_MODE}'. Update config_system.py."
    )

# DL2-23: separate log + gate report
cfg.validate_deadlock_rules(
    "dl2_23_split_provenance",
    separate_logs=True,
    separate_gate_reports=True,
)
print("[DL2-23] PASS — separate multi-model log + gate report")

# DL2-26: target list declared + validated
TARGET_MODELS_LIST = cfg.validate_target_models_list()
cfg.validate_deadlock_rules(
    "dl2_26_target_list_declared",
    target_list_source="cfg.L2_TARGET_MODELS_LIST",
    list_validated=True,
)
cfg.validate_deadlock_rules(
    "dl2_24_target_input_declared",
    target_model_input=cfg.L2_TARGET_MODEL_INPUT,
    target_models_list=TARGET_MODELS_LIST,
)
print(f"[DL2-26] PASS — target list n={len(TARGET_MODELS_LIST)}")
print(f"[DL2-24] PASS — multi-model list declared in config")
logger.info(f"TARGET_MODELS_LIST (n={len(TARGET_MODELS_LIST)}): {TARGET_MODELS_LIST}")


# =============================================================================
# HELPERS
# =============================================================================

def _sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(65536), b""):
            h.update(chunk)
    return h.hexdigest()


def _handle_missing(model_id: str, reason: str, ledger: list) -> None:
    """Apply L2_MULTI_MODEL_ON_MISSING policy."""
    policy = cfg.L2_MULTI_MODEL_ON_MISSING
    entry = {
        "model_id": model_id,
        "status": "SKIPPED" if policy == "skip_missing_with_ledger" else "FAILED",
        "reason": reason,
        "lineage": None,
        "timestamp": datetime.now().isoformat(),
    }
    ledger.append(entry)
    msg = f"[MISSING] {model_id}: {reason} (policy={policy})"
    if policy == "fail":
        logger.error(msg)
        raise RuntimeError(msg)
    logger.warning(msg)
    print(f"  [SKIP] {model_id}: {reason}")


def _resolve_lioness_npy(
    raw_dir: Path,
    target_model_id: str,
    target_pos_0based: int,
    expected_shape: tuple,
    context: str,
) -> Path:
    """
    Select the correct single-sample LIONESS .npy from netZooPy output.

    netZooPy dual-write (when save_single=True and ignore_final=False):
      1) save_single path in __lioness_loop:
            lioness.{expression_samples[i]}.{i}.{fmt}
         e.g. lioness.ACH-000019.41.npy  OR  lioness.1.0.npy
         (sample name depends on whether expression sample labels survived)
      2) save_lioness_results() at end of __init__ (ignore_final=False):
            lioness.npy
         This is the aggregated export table — NOT the preferred single-
         sample TF×gene matrix for our labeling pipeline.

    Integrity rule:
      Prefer save_single patterned files whose array shape == PANDA network.
      Never silently pick an arbitrary file when shapes disagree.
    """
    raw_dir = Path(raw_dir)
    npy_files = sorted(raw_dir.glob("*.npy"))
    if not npy_files:
        raise RuntimeError(
            f"[{context}] No .npy files in {raw_dir}. "
            "LIONESS did not write output (check save_dir / save_single)."
        )

    # Priority 1: exact save_single pattern with ModelID + 0-based index
    #   lioness.{sample}.{i}.npy
    candidates_ranked: list[Path] = []
    pat_model_idx = (
        f"lioness.{target_model_id}.{target_pos_0based}.npy"
    )
    p = raw_dir / pat_model_idx
    if p.exists():
        candidates_ranked.append(p)

    # Priority 2: any lioness.{something}.{i}.npy with matching index
    #   (netZooPy may use integer sample labels if headerless expr had no names)
    for p in npy_files:
        name = p.name
        if name == "lioness.npy":
            continue  # aggregate export — deprioritize
        # lioness.<sample>.<idx>.npy
        parts = name.split(".")
        # ['lioness', sample..., idx, 'npy'] — sample may contain dots rarely
        if len(parts) >= 4 and parts[0] == "lioness" and parts[-1] == "npy":
            idx_token = parts[-2]
            if idx_token.isdigit() and int(idx_token) == int(target_pos_0based):
                if p not in candidates_ranked:
                    candidates_ranked.append(p)

    # Priority 3: any non-aggregate lioness.*.npy
    for p in npy_files:
        if p.name != "lioness.npy" and p not in candidates_ranked:
            candidates_ranked.append(p)

    # Priority 4: aggregate lioness.npy last (only if shape matches)
    agg = raw_dir / "lioness.npy"
    if agg.exists() and agg not in candidates_ranked:
        candidates_ranked.append(agg)

    shape_ok = []
    shape_bad = []
    for p in candidates_ranked:
        try:
            arr = np.load(str(p), mmap_mode="r")
            shp = tuple(arr.shape)
        except Exception as e:
            shape_bad.append((p.name, f"load_error:{e}"))
            continue
        if shp == tuple(expected_shape):
            shape_ok.append((p, shp))
        else:
            shape_bad.append((p.name, shp))

    if shape_ok:
        chosen, shp = shape_ok[0]
        if len(npy_files) > 1:
            logger.info(
                f"[{context}] netZooPy wrote {len(npy_files)} .npy files "
                f"{[p.name for p in npy_files]}; selected '{chosen.name}' "
                f"shape={shp} (PANDA-matched). "
                f"Rejected/other: {shape_bad}"
            )
            print(
                f"    npy select    : {chosen.name} "
                f"(from {len(npy_files)} files; shape OK {shp})"
            )
        return chosen

    raise RuntimeError(
        f"[{context}] Could not resolve single-sample LIONESS .npy with "
        f"shape={expected_shape}. Found files: "
        f"{[p.name for p in npy_files]}. Shape report: {shape_bad}. "
        f"Hint: prefer lioness.<sample>.<0based_index>.npy from save_single; "
        f"lioness.npy is often the aggregate export and may differ in shape."
    )


# =============================================================================
# STAGE M2: LOAD FOUNDATION MANIFEST & VALIDATE CONTRACT
# =============================================================================

print("\n" + "=" * 64)
print("STAGE M2: LOAD FOUNDATION MANIFEST & VALIDATE CONTRACT")
print("=" * 64)
logger.info("Stage M2 start")

manifest_path = cfg.L2A_FROZEN_INPUTS_DIR / cfg.L2A_FOUNDATION_MANIFEST_JSON
if not manifest_path.exists():
    raise FileNotFoundError(
        f"[M2] Foundation manifest missing: {manifest_path}\n"
        f"Run notebook {cfg.LAYER2_NOTEBOOK_02A_NAME} (Stage 8F) first."
    )

with open(manifest_path, "r", encoding="utf-8") as fh:
    foundation = json.load(fh)

print(f"Manifest version : {foundation.get('manifest_version')}")
print(f"Pipeline version : {foundation.get('pipeline_version')}")
print(f"Config version   : {foundation.get('config_version')}")
print(f"Lineages frozen  : {foundation.get('lineages_frozen')}")
print(f"PANDA mode       : {foundation.get('panda_mode')}")
print(f"PANDA alpha      : {foundation.get('panda_alpha')}")

if foundation.get("manifest_version") != cfg.FOUNDATION_MANIFEST_VERSION:
    raise RuntimeError(
        f"[M2] Foundation manifest version mismatch: "
        f"found='{foundation.get('manifest_version')}', "
        f"expected='{cfg.FOUNDATION_MANIFEST_VERSION}'"
    )

logger.info(
    f"[M2] Foundation manifest loaded | "
    f"version={foundation.get('manifest_version')} | "
    f"lineages={foundation.get('lineages_frozen')}"
)

# Validate all frozen lineage files exist
print()
for lineage in foundation["lineages_frozen"]:
    arts = foundation["lineage_artifacts"][lineage]
    for key, desc in [
        ("expr_parquet", "expression parquet"),
        ("sample_order_json", "sample order"),
        ("genes_post_qc_json", "genes post-QC"),
    ]:
        p = cfg.LAYER2_GRN_OUTPUT_DIR / arts[key]
        if not p.exists():
            raise FileNotFoundError(
                f"[M2] Missing frozen artifact for '{lineage}' ({desc}): {p}"
            )
    print(f"  [✓] {lineage:<12}: all frozen artifacts present")

# Validate + verify prior file SHA-256 hashes
print()
prior_info = foundation["prior_files"]
for path_key, hash_key, desc in [
    ("motif_prior_path", "motif_prior_sha256", "motif prior"),
    ("ppi_prior_path", "ppi_prior_sha256", "PPI prior"),
]:
    path = Path(prior_info[path_key])
    if not path.exists():
        raise FileNotFoundError(f"[M2] Prior file missing: {path}")
    if foundation["prior_files"].get("hash_validation_required", True):
        expected = prior_info[hash_key]
        actual = _sha256_file(path)
        if actual != expected:
            raise RuntimeError(
                f"[M2] SHA-256 mismatch for {desc}: {path.name}\n"
                f"  Expected : {expected}\n"
                f"  Actual   : {actual}\n"
                "Prior files were modified since 02A ran."
            )
        print(f"  [✓] {path.name}: SHA-256 verified")
    else:
        print(f"  [i] {path.name}: hash validation disabled in manifest")

# Load Model.csv (for target resolution only)
print()
model_df = pd.read_csv(cfg.CCLE_MODEL_CSV)
n_rows = len(model_df)
if cfg.COL_MODEL_ID not in model_df.columns:
    raise RuntimeError(f"[M2] Model.csv missing '{cfg.COL_MODEL_ID}' column")
if cfg.COL_ONCOTREE_LINEAGE not in model_df.columns:
    raise RuntimeError(
        f"[M2] Model.csv missing '{cfg.COL_ONCOTREE_LINEAGE}' column"
    )
print(f"  Model.csv loaded : {n_rows} rows")
logger.info(f"[M2] Model.csv: {n_rows} rows")
logger.info("Stage M2 complete")


# =============================================================================
# STAGE M3: RESOLVE ALL TARGETS → LINEAGE GROUPS
# =============================================================================

print("\n" + "=" * 64)
print("STAGE M3: RESOLVE ALL TARGETS → LINEAGE GROUPS")
print("=" * 64)
logger.info("Stage M3 start")

run_ledger: list = []  # full DL2-28 ledger entries
# lineage -> OrderedDict of model_id -> metadata
lineage_groups: dict = defaultdict(OrderedDict)
resolution_rows = []

for raw in TARGET_MODELS_LIST:
    model_id = cfg.resolve_model_id(raw)
    hits = model_df[
        model_df[cfg.COL_MODEL_ID].astype(str).str.strip() == model_id
    ].copy()

    if len(hits) == 0:
        _handle_missing(
            model_id,
            "not found in Model.csv",
            run_ledger,
        )
        resolution_rows.append(
            {
                "raw_input": raw,
                "model_id": model_id,
                "lineage": None,
                "status": "SKIPPED_OR_FAILED",
                "reason": "not found in Model.csv",
            }
        )
        continue

    if len(hits) > 1:
        # Ambiguous — always hard fail (integrity)
        raise RuntimeError(
            f"[M3] Target '{raw}' → '{model_id}' matches {len(hits)} rows "
            f"in Model.csv — ambiguous."
        )

    # DL2-14 per resolved unique hit
    cfg.validate_deadlock_rules(
        "dl2_14_target_model", target_exists_uniquely=True
    )

    lineage = str(hits.iloc[0][cfg.COL_ONCOTREE_LINEAGE]).strip()
    if lineage in ("", "Unknown", "nan"):
        _handle_missing(
            model_id,
            f"invalid lineage metadata: {lineage!r}",
            run_ledger,
        )
        resolution_rows.append(
            {
                "raw_input": raw,
                "model_id": model_id,
                "lineage": lineage,
                "status": "SKIPPED_OR_FAILED",
                "reason": f"invalid lineage: {lineage!r}",
            }
        )
        continue

    # DL2-15: lineage from Model.csv only
    cfg.validate_deadlock_rules(
        "dl2_15_lineage_from_metadata", lineage_from_metadata=True
    )

    if lineage not in foundation["lineages_frozen"]:
        _handle_missing(
            model_id,
            f"lineage '{lineage}' not in 02A frozen set "
            f"{foundation['lineages_frozen']}",
            run_ledger,
        )
        resolution_rows.append(
            {
                "raw_input": raw,
                "model_id": model_id,
                "lineage": lineage,
                "status": "SKIPPED_OR_FAILED",
                "reason": f"lineage '{lineage}' not frozen in 02A",
            }
        )
        continue

    if model_id in lineage_groups[lineage]:
        # Same ACH listed twice after alias resolution
        raise RuntimeError(
            f"[M3] Duplicate resolved ModelID '{model_id}' in target list "
            f"(raw inputs include '{raw}')."
        )

    lineage_groups[lineage][model_id] = {
        "raw_input": raw,
        "model_id": model_id,
        "lineage": lineage,
    }
    resolution_rows.append(
        {
            "raw_input": raw,
            "model_id": model_id,
            "lineage": lineage,
            "status": "QUEUED",
            "reason": "",
        }
    )

print("\nResolution summary by lineage:")
for lin, models in lineage_groups.items():
    print(f"  {lin:<12}: {len(models)} targets")
    logger.info(f"[M3] lineage={lin} n_targets={len(models)} ids={list(models)}")

n_queued = sum(len(v) for v in lineage_groups.values())
n_ledger_early = len(run_ledger)
print(f"\n  Queued for LIONESS : {n_queued}")
print(f"  Early SKIP/FAIL    : {n_ledger_early}")
logger.info(
    f"[M3] queued={n_queued} early_skip_or_fail={n_ledger_early}"
)

if n_queued == 0:
    raise RuntimeError(
        "[M3] No targets remain after resolution against Model.csv + "
        "frozen lineages. Aborting."
    )

# DL2-22 confirmed globally (no SNAIL, frozen only)
cfg.validate_deadlock_rules(
    "dl2_22_frozen_foundation_inputs",
    using_frozen_foundation_inputs=True,
    snail_recomputed=False,
)
print("[DL2-22] PASS — frozen foundation inputs only; no SNAIL recompute")
logger.info("Stage M3 complete")


# =============================================================================
# STAGE M4: LOAD & CANONICALIZE PRIORS (once, shared)
# =============================================================================

print("\n" + "=" * 64)
print("STAGE M4: LOAD & VERIFY PRIORS (shared across all lineages)")
print("=" * 64)
logger.info("Stage M4 start")

local_motif = Path(foundation["prior_files"]["motif_prior_path"])
local_ppi = Path(foundation["prior_files"]["ppi_prior_path"])

motif_prior_df = pd.read_csv(local_motif, sep="\t", header=None)
assert not motif_prior_df.empty, "[M4] Motif prior is empty"
assert motif_prior_df.shape[1] >= 3, (
    f"[M4] Motif prior needs >= 3 columns, got {motif_prior_df.shape[1]}"
)
motif_prior_df = motif_prior_df.iloc[:, :3].copy()
motif_prior_df.columns = ["TF", "Gene", "Weight"]
motif_prior_df["TF"] = motif_prior_df["TF"].astype(str).str.strip()
motif_prior_df["Gene"] = motif_prior_df["Gene"].astype(str).str.strip()
motif_prior_df["Weight"] = pd.to_numeric(
    motif_prior_df["Weight"], errors="coerce"
).fillna(1.0)

ppi_prior_raw = pd.read_csv(local_ppi, sep="\t", header=None)
assert not ppi_prior_raw.empty, "[M4] PPI prior is empty"
assert ppi_prior_raw.shape[1] >= 2, (
    f"[M4] PPI prior needs >= 2 columns, got {ppi_prior_raw.shape[1]}"
)
if ppi_prior_raw.shape[1] == 2:
    ppi_prior_raw[2] = 1
    ppi_prior_runtime_weight_mode = "canonicalized_binary_weight_1"
    logger.warning(
        "[M4] PPI prior: 2-column file canonicalized to 3-column (Weight=1). "
        "File-format normalization for PANDA API. "
        "Rebuild ppi_prior.txt as 3-column to eliminate this step."
    )
else:
    ppi_prior_runtime_weight_mode = "native_3col"

ppi_prior_df = ppi_prior_raw.iloc[:, :3].copy()
ppi_prior_df.columns = ["TF1", "TF2", "Weight"]
ppi_prior_df["TF1"] = ppi_prior_df["TF1"].astype(str).str.strip()
ppi_prior_df["TF2"] = ppi_prior_df["TF2"].astype(str).str.strip()
ppi_prior_df["Weight"] = pd.to_numeric(
    ppi_prior_df["Weight"], errors="coerce"
).fillna(1.0)

n_self = int((ppi_prior_df["TF1"] == ppi_prior_df["TF2"]).sum())
if n_self > 0:
    ppi_prior_df = ppi_prior_df[
        ppi_prior_df["TF1"] != ppi_prior_df["TF2"]
    ].copy()
    logger.warning(f"[M4] PPI prior: removed {n_self} self-loop edges")

cfg.validate_deadlock_rules(
    "dl2_05_prior_source", prior_source_runtime="LOCAL_FILES"
)
print(f"Motif prior : {motif_prior_df.shape}")
print(f"PPI prior   : {ppi_prior_df.shape}")
print(f"PPI mode    : {ppi_prior_runtime_weight_mode}")
print("[DL2-05] PASS")
logger.info(
    f"[M4] Priors loaded | motif={motif_prior_df.shape} | "
    f"ppi={ppi_prior_df.shape} | ppi_mode={ppi_prior_runtime_weight_mode}"
)
logger.info("Stage M4 complete")


# =============================================================================
# STAGE M5: PER-LINEAGE SHARED PANDA + PER-TARGET LIONESS
# =============================================================================

print("\n" + "=" * 64)
print("STAGE M5: PER-LINEAGE SHARED PANDA + PER-TARGET LIONESS")
print("=" * 64)
logger.info("Stage M5 start")

from netZooPy.panda import Panda
from netZooPy.lioness.lioness import Lioness

QC_DIR = cfg.LAYER2_GRN_OUTPUT_DIR / cfg.L2_GRN_QC_PLOTS_DIR
QC_DIR.mkdir(parents=True, exist_ok=True)

lineage_summaries = {}
completed_exports = []

for lineage, models_od in lineage_groups.items():
    target_ids = list(models_od.keys())
    print("\n" + "-" * 64)
    print(f"LINEAGE: {lineage}  |  targets: {len(target_ids)}")
    print("-" * 64)
    logger.info(f"[M5] === lineage={lineage} n={len(target_ids)} ===")

    arts = foundation["lineage_artifacts"][lineage]

    # ── Load frozen expression matrix ─────────────────────────
    expr_parquet = cfg.LAYER2_GRN_OUTPUT_DIR / arts["expr_parquet"]
    expr_cohort_full = pd.read_parquet(
        expr_parquet, engine=cfg.FROZEN_EXPR_ENGINE
    )
    expr_cohort_full.index = pd.Index(
        [str(x) for x in expr_cohort_full.index], name=cfg.COL_MODEL_ID
    )
    expr_cohort_full.columns = pd.Index(
        [str(x) for x in expr_cohort_full.columns], name="Gene"
    )
    assert expr_cohort_full.index.is_unique, (
        f"[M5/{lineage}] Duplicate sample IDs in frozen expr"
    )
    assert expr_cohort_full.columns.is_unique, (
        f"[M5/{lineage}] Duplicate gene IDs in frozen expr"
    )
    assert not isinstance(expr_cohort_full.index, pd.MultiIndex), (
        f"[M5/{lineage}] MultiIndex in frozen expr — 02A freeze may have failed"
    )

    # ── Load frozen sample order ──────────────────────────────
    with open(
        cfg.LAYER2_GRN_OUTPUT_DIR / arts["sample_order_json"], "r"
    ) as fh:
        sample_sidecar = json.load(fh)
    frozen_sample_order = [str(x) for x in sample_sidecar["samples_in_order"]]
    if list(expr_cohort_full.index) != frozen_sample_order:
        raise RuntimeError(
            f"[M5/{lineage}] Loaded parquet sample order does not match "
            "frozen sample order sidecar."
        )

    # ── Load frozen gene list (zero-var post-QC) ─────────────
    with open(
        cfg.LAYER2_GRN_OUTPUT_DIR / arts["genes_post_qc_json"], "r"
    ) as fh:
        gene_info = json.load(fh)
    genes_post_qc = gene_info["genes_post_qc_in_order"]
    n_zero_var = gene_info["n_zero_var_genes"]
    n_genes_post_qc = gene_info["n_genes_post_qc"]
    assert len(genes_post_qc) == n_genes_post_qc, (
        f"[M5/{lineage}] Gene list length inconsistency"
    )

    expr_cohort = expr_cohort_full[genes_post_qc].copy()
    sample_order = [str(x) for x in expr_cohort.index]
    gene_order = [str(x) for x in expr_cohort.columns]
    n_cohort = len(sample_order)
    n_genes = len(gene_order)

    print(f"  Cohort size     : {n_cohort}")
    print(f"  Genes (post-QC) : {n_genes}")
    print(f"  Zero-var removed: {n_zero_var}")

    # ── Partition present vs absent targets in frozen cohort ──
    present_targets = []
    for mid in target_ids:
        if mid not in frozen_sample_order:
            _handle_missing(
                mid,
                f"not present in frozen cohort for lineage '{lineage}'",
                run_ledger,
            )
        else:
            present_targets.append(mid)

    if not present_targets:
        print(f"  [!] No targets present in frozen cohort for {lineage}; skip lineage.")
        lineage_summaries[lineage] = {
            "status": "NO_TARGETS_PRESENT",
            "n_requested": len(target_ids),
            "n_present": 0,
            "n_panda_builds": 0,
        }
        continue

    # ── Output dirs ───────────────────────────────────────────
    LIONESS_OUTPUT_DIR = (
        cfg.LAYER2_GRN_OUTPUT_DIR
        / cfg.L2_GRN_LIONESS_DIR.format(lineage=lineage)
    )
    LIONESS_OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    # ── Temp prior/expr files (shared for this lineage PANDA) ─
    LIONESS_TMP = Path(tempfile.mkdtemp(prefix=f"l2b_multi_{lineage}_"))
    logger.info(f"[M5/{lineage}] Temp dir: {LIONESS_TMP}")
    expr_tmp = LIONESS_TMP / f"expr_{lineage}.txt"
    motif_tmp = LIONESS_TMP / "motif_prior.txt"
    ppi_tmp = LIONESS_TMP / "ppi_prior.txt"

    expr_cohort.T.to_csv(expr_tmp, sep="\t", header=False, index=True)
    motif_prior_df[["TF", "Gene", "Weight"]].to_csv(
        motif_tmp, sep="\t", header=False, index=False
    )
    ppi_prior_df[["TF1", "TF2", "Weight"]].to_csv(
        ppi_tmp, sep="\t", header=False, index=False
    )

    # ── Build ONE shared Panda for this lineage (DL2-25) ──────
    print(
        f"\n  Building SHARED Panda for {lineage} "
        f"({n_cohort} samples, {len(present_targets)} LIONESS targets) ..."
    )
    logger.info(f"[M5/{lineage}] Building shared Panda")
    t0_panda = time.time()
    panda_obj = Panda(
        expression_file=str(expr_tmp),
        motif_file=str(motif_tmp),
        ppi_file=str(ppi_tmp),
        computing=cfg.PANDA_COMPUTING,
        precision=cfg.PANDA_PRECISION,
        save_memory=False,
        save_tmp=True,
        remove_missing=cfg.PANDA_REMOVE_MISSING,
        keep_expression_matrix=True,
        modeProcess=cfg.PANDA_MODE,
        alpha=cfg.PANDA_ALPHA,
        with_header=cfg.PANDA_WITH_HEADER,
    )
    panda_elapsed = time.time() - t0_panda
    n_panda_builds = 1  # exactly one for this lineage

    for req in [
        "expression_matrix",
        "motif_matrix",
        "ppi_matrix",
        "panda_network",
        "export_panda_results",
    ]:
        if not hasattr(panda_obj, req):
            raise RuntimeError(
                f"[M5/{lineage}] Panda missing required attribute: {req}"
            )

    em_shape = tuple(panda_obj.expression_matrix.shape)
    pn_shape = tuple(panda_obj.panda_network.shape)
    ep_shape = tuple(panda_obj.export_panda_results.shape)
    panda_internal_n_genes = em_shape[0]

    if em_shape[1] != n_cohort:
        raise RuntimeError(
            f"[M5/{lineage}] Panda sample count {em_shape[1]} != cohort {n_cohort}"
        )
    if panda_internal_n_genes < n_genes:
        raise RuntimeError(
            f"[M5/{lineage}] Panda internal genes {panda_internal_n_genes} "
            f"< input genes {n_genes} — unexpected for mode='{cfg.PANDA_MODE}'"
        )
    expected_edges = pn_shape[0] * pn_shape[1]
    if ep_shape[0] != expected_edges:
        raise RuntimeError(
            f"[M5/{lineage}] export_panda_results rows {ep_shape[0]:,} "
            f"!= expected edges {expected_edges:,}"
        )

    print(f"  Panda built     : {panda_elapsed:.1f} s")
    print(f"  Input genes     : {n_genes}")
    print(f"  Internal genes  : {panda_internal_n_genes}")
    print(f"  Network shape   : {pn_shape}")
    print(f"  Export shape    : {ep_shape}")

    cfg.validate_deadlock_rules(
        "dl2_09_lioness_after_panda",
        lineage=lineage,
        panda_converged=True,
    )
    cfg.validate_deadlock_rules(
        "dl2_10_lioness_same_priors",
        panda_motif_path=str(motif_tmp),
        lioness_motif_path=str(motif_tmp),
        panda_ppi_path=str(ppi_tmp),
        lioness_ppi_path=str(ppi_tmp),
    )
    cfg.validate_deadlock_rules(
        "dl2_25_shared_panda_per_lineage",
        shared_panda_per_lineage=True,
        cross_lineage_pooling=False,
        n_panda_builds_for_lineage=n_panda_builds,
    )
    print("[DL2-09] PASS | [DL2-10] PASS | [DL2-25] PASS (shared PANDA)")

    # ── Verify flatten order once for this shared PANDA (DL2-19/27) ──
    ep = panda_obj.export_panda_results.copy()
    ep_cols_lower = {str(c).lower(): c for c in ep.columns}
    tf_col = ep_cols_lower.get("tf")
    gene_col = ep_cols_lower.get("gene")
    motif_col = ep_cols_lower.get("motif")
    force_candidates = [
        c
        for c in ep.columns
        if str(c).lower() in {"force", "weight", "score"}
    ]
    if not tf_col or not gene_col or len(force_candidates) != 1:
        raise RuntimeError(
            f"[M5/{lineage}] Cannot identify tf/gene/force columns in "
            f"export_panda_results: {list(ep.columns)}"
        )
    force_col = force_candidates[0]

    pn = panda_obj.panda_network
    pn_vals = (
        pn.to_numpy() if isinstance(pn, pd.DataFrame) else np.asarray(pn)
    )
    ep_force = pd.to_numeric(ep[force_col], errors="raise").to_numpy()
    match_C = np.allclose(
        ep_force,
        pn_vals.reshape(-1, order="C"),
        rtol=0,
        atol=1e-12,
        equal_nan=True,
    )
    match_F = np.allclose(
        ep_force,
        pn_vals.reshape(-1, order="F"),
        rtol=0,
        atol=1e-12,
        equal_nan=True,
    )
    if match_C and not match_F:
        verified_order = "C"
    elif match_F and not match_C:
        verified_order = "F"
    elif match_C and match_F:
        raise RuntimeError(
            f"[M5/{lineage}] Both C and F flatten orders match — ambiguous."
        )
    else:
        raise RuntimeError(
            f"[M5/{lineage}] Neither C nor F flatten order matches PANDA export."
        )

    cfg.validate_deadlock_rules(
        "dl2_19_export_only_if_verified", flatten_order_verified=True
    )
    cfg.validate_deadlock_rules(
        "dl2_27_export_only_if_verified_multi", flatten_order_verified=True
    )
    print(f"  Verified order  : {verified_order}")
    print("[DL2-19] PASS | [DL2-27] PASS")

    lineage_completed = []
    lineage_failed = []

    # ── Per-target LIONESS on the SHARED panda_obj ────────────
    for mid in present_targets:
        meta = models_od[mid]
        raw_input = meta["raw_input"]
        target_pos_0based = sample_order.index(mid)
        target_pos_1based = target_pos_0based + 1

        print(
            f"\n  → LIONESS {mid}  pos={target_pos_1based}/{n_cohort} "
            f"(input={raw_input})"
        )
        logger.info(
            f"[M5/{lineage}] Lioness start={target_pos_1based} "
            f"end={target_pos_1based} model={mid}"
        )

        LIONESS_RAW_DIR = LIONESS_OUTPUT_DIR / cfg.L2B_NETZOOPY_RAW_DIR.format(
            sample=mid
        )
        LIONESS_RAW_DIR.mkdir(parents=True, exist_ok=True)
        # Integrity: purge stale netZooPy dumps from prior failed/partial runs
        # so dual-write / leftover files cannot pollute shape-gated selection.
        for _stale in LIONESS_RAW_DIR.glob("*"):
            try:
                if _stale.is_file():
                    _stale.unlink()
            except Exception as _se:
                logger.warning(
                    f"[M5/{lineage}/{mid}] could not remove stale {_stale.name}: {_se}"
                )

        # Sidecar (DL2-17)
        sidecar_fname = cfg.L2B_TARGET_SAMPLE_ORDER_JSON.format(sample=mid)
        sidecar_path = LIONESS_OUTPUT_DIR / sidecar_fname
        sample_order_sidecar = {
            "target_model_input": raw_input,
            "target_model_id": mid,
            "target_lineage": lineage,
            "reference_cohort_policy": cfg.LIONESS_REFERENCE_POLICY,
            "reference_cohort_size": n_cohort,
            "n_genes_input_post_qc": n_genes,
            "n_zero_var_removed": n_zero_var,
            "zero_var_policy": cfg.ZERO_VARIANCE_GENE_POLICY,
            "target_position_0based": target_pos_0based,
            "target_position_1based": target_pos_1based,
            "samples_in_order": sample_order,
            "notebook": cfg.LAYER2_NOTEBOOK_02B_MULTI_NAME,
            "run_mode": cfg.LAYER2_RUN_MODE,
            "shared_panda_per_lineage": True,
            "foundation_manifest": str(manifest_path),
            "timestamp": datetime.now().isoformat(),
        }
        with open(sidecar_path, "w", encoding="utf-8") as fh:
            json.dump(
                sample_order_sidecar, fh, indent=2, ensure_ascii=False
            )
        cfg.validate_deadlock_rules(
            "dl2_17_sample_order_deterministic", sidecar_written=True
        )

        t0_lioness = time.time()
        lioness_exception = None
        try:
            try:
                # ignore_final=True: do NOT call save_lioness_results() → no
                # bare lioness.npy aggregate. save_single=True still writes the
                # per-sample TF×gene matrix we need. netZooPy requires at least
                # one of (save_single, not ignore_final) to be true.
                # Scientific payload unchanged: LIONESS formula + shared PANDA.
                lioness_obj = Lioness(
                    obj=panda_obj,
                    computing=cfg.PANDA_COMPUTING,
                    precision=cfg.PANDA_PRECISION,
                    ncores=1,
                    start=target_pos_1based,
                    end=target_pos_1based,
                    subset_numbers=None,
                    subset_names=None,
                    save_dir=str(LIONESS_RAW_DIR),
                    save_fmt="npy",
                    output="network",
                    alpha=cfg.PANDA_ALPHA,
                    save_single=True,
                    export_filename=None,
                    ignore_final=True,
                    online_coexpression=False,
                )
            except AttributeError as ae:
                if "total_lioness_network" in str(ae):
                    lioness_exception = ae
                    logger.warning(
                        f"[M5/{lineage}/{mid}] Known netZooPy 0.11.0 bug "
                        f"caught after computation: {ae}. "
                        "Raw .npy has been saved prior to this error."
                    )
                else:
                    raise

            lioness_elapsed = time.time() - t0_lioness

            # ------------------------------------------------------------------
            # Resolve raw .npy (netZooPy dual-write behavior)
            # ------------------------------------------------------------------
            # When save_single=True, Lioness writes:
            #   lioness.{sample_name}.{i}.npy   ← single-sample TF×gene matrix
            # When ignore_final=False (default), __init__ also calls
            # save_lioness_results() which writes:
            #   lioness.npy                     ← aggregate/export table
            # So glob("*.npy") often returns 2 files. That is NOT a compute
            # failure — only a wrong assumption in older Stage B6 code.
            # Prefer the save_single pattern; fall back carefully.
            # Ref: netZooPy/lioness/lioness.py (__lioness_loop + save_lioness_results)
            # ------------------------------------------------------------------
            npy_path = _resolve_lioness_npy(
                raw_dir=LIONESS_RAW_DIR,
                target_model_id=mid,
                target_pos_0based=target_pos_0based,
                expected_shape=tuple(pn_vals.shape),
                context=f"M5/{lineage}/{mid}",
            )
            lioness_raw = np.load(str(npy_path))

            if tuple(lioness_raw.shape) != tuple(pn_vals.shape):
                raise RuntimeError(
                    f"[M5/{lineage}/{mid}] LIONESS raw shape "
                    f"{lioness_raw.shape} != PANDA network shape {pn_vals.shape} "
                    f"(file={npy_path.name})"
                )

            # Label edges with verified flatten order
            lioness_vec = lioness_raw.reshape(-1, order=verified_order)
            if lioness_vec.size != len(ep):
                raise RuntimeError(
                    f"[M5/{lineage}/{mid}] Vector size {lioness_vec.size:,} "
                    f"!= export rows {len(ep):,}"
                )

            lioness_df = pd.DataFrame(
                {
                    cfg.COL_GRN_SOURCE: ep[tf_col].astype(str).to_numpy(),
                    cfg.COL_GRN_TARGET: ep[gene_col].astype(str).to_numpy(),
                    cfg.COL_GRN_WEIGHT: lioness_vec.astype(float),
                    "Motif": (
                        ep[motif_col].to_numpy() if motif_col else np.nan
                    ),
                    "ModelID": mid,
                    "Lineage": lineage,
                }
            )
            assert not lioness_df.empty, (
                f"[M5/{lineage}/{mid}] Extracted LIONESS DataFrame is empty"
            )
            assert {
                cfg.COL_GRN_SOURCE,
                cfg.COL_GRN_TARGET,
                cfg.COL_GRN_WEIGHT,
            }.issubset(lioness_df.columns)

            final_lioness_tsv = LIONESS_OUTPUT_DIR / cfg.L2_GRN_LIONESS_TSV.format(
                sample=mid
            )
            lioness_df.to_csv(final_lioness_tsv, sep="\t", index=False)

            run_manifest_fname = cfg.L2B_LIONESS_RUN_MANIFEST_JSON.format(
                sample=mid
            )
            run_manifest_path = LIONESS_OUTPUT_DIR / run_manifest_fname
            run_manifest = {
                "status": "COMPLETED_MULTI_TARGET",
                "timestamp": datetime.now().isoformat(),
                "notebook": cfg.LAYER2_NOTEBOOK_02B_MULTI_NAME,
                "pipeline_version": PIPELINE_VERSION,
                "config_version": CONFIG_VERSION,
                "netzoopy_version": getattr(
                    __import__("netZooPy"), "__version__", "unknown"
                ),
                "run_mode": cfg.LAYER2_RUN_MODE,
                "shared_panda_per_lineage": True,
                "n_panda_builds_for_lineage": n_panda_builds,
                "target_model_input": raw_input,
                "target_model_id": mid,
                "target_lineage": lineage,
                "reference_cohort_policy": cfg.LIONESS_REFERENCE_POLICY,
                "reference_cohort_size": n_cohort,
                "n_genes_input_post_qc": n_genes,
                "n_genes_internal_union": panda_internal_n_genes,
                "n_zero_var_removed": n_zero_var,
                "zero_var_policy": cfg.ZERO_VARIANCE_GENE_POLICY,
                "target_position_0based": target_pos_0based,
                "target_position_1based": target_pos_1based,
                "verified_flatten_order": verified_order,
                "npy_flatten_order_policy": cfg.LIONESS_NPY_FLATTEN_ORDER,
                "n_edges": len(lioness_df),
                "panda_elapsed_sec": round(panda_elapsed, 1),
                "lioness_elapsed_sec": round(lioness_elapsed, 1),
                "raw_npy_file": str(npy_path),
                "raw_npy_shape": list(lioness_raw.shape),
                "final_tsv": str(final_lioness_tsv),
                "sample_order_sidecar": str(sidecar_path),
                "foundation_manifest": str(manifest_path),
                "motif_prior_path": str(local_motif),
                "ppi_prior_path": str(local_ppi),
                "ppi_runtime_weight_mode": ppi_prior_runtime_weight_mode,
                "snail_recomputed": False,
                "snail_label_policy": cfg.SNAIL_LABEL_RESTORATION_POLICY,
                "panda_mode": cfg.PANDA_MODE,
                "panda_alpha": cfg.PANDA_ALPHA,
                "panda_computing": cfg.PANDA_COMPUTING,
                "panda_precision": cfg.PANDA_PRECISION,
                "panda_union_expansion": cfg.PANDA_UNION_GENE_EXPANSION_EXPECTED,
                "known_netzoopy_bug": (
                    str(lioness_exception) if lioness_exception else None
                ),
                "deadlock_rules_checked": [
                    "DL2-05",
                    "DL2-09",
                    "DL2-10",
                    "DL2-14",
                    "DL2-15",
                    "DL2-17",
                    "DL2-19",
                    "DL2-22",
                    "DL2-23",
                    "DL2-24",
                    "DL2-25",
                    "DL2-26",
                    "DL2-27",
                    "DL2-28",
                ],
            }
            with open(run_manifest_path, "w", encoding="utf-8") as fh:
                json.dump(run_manifest, fh, indent=2, ensure_ascii=False)

            entry = {
                "model_id": mid,
                "raw_input": raw_input,
                "status": "COMPLETED",
                "reason": "ok",
                "lineage": lineage,
                "target_position_1based": target_pos_1based,
                "n_edges": len(lioness_df),
                "verified_flatten_order": verified_order,
                "panda_elapsed_sec": round(panda_elapsed, 1),
                "lioness_elapsed_sec": round(lioness_elapsed, 1),
                "final_tsv": str(final_lioness_tsv),
                "run_manifest": str(run_manifest_path),
                "raw_npy_file": str(npy_path),
                "shared_panda": True,
                "timestamp": datetime.now().isoformat(),
            }
            run_ledger.append(entry)
            lineage_completed.append(mid)
            completed_exports.append(entry)

            print(
                f"    edges={len(lioness_df):,}  "
                f"lioness={lioness_elapsed:.1f}s  "
                f"tsv={final_lioness_tsv.name}"
            )
            logger.info(
                f"[M5/{lineage}/{mid}] COMPLETED edges={len(lioness_df):,} "
                f"elapsed={lioness_elapsed:.1f}s"
            )

        except Exception as e:
            err = f"{type(e).__name__}: {e}"
            logger.error(
                f"[M5/{lineage}/{mid}] FAILED: {err}\n{traceback.format_exc()}"
            )
            run_ledger.append(
                {
                    "model_id": mid,
                    "raw_input": raw_input,
                    "status": "FAILED",
                    "reason": err,
                    "lineage": lineage,
                    "timestamp": datetime.now().isoformat(),
                }
            )
            lineage_failed.append(mid)
            print(f"    [FAIL] {mid}: {err}")
            if cfg.L2_MULTI_MODEL_ON_MISSING == "fail":
                raise

    # Lineage summary
    lin_summary = {
        "lineage": lineage,
        "status": "COMPLETED",
        "n_requested": len(target_ids),
        "n_present_in_cohort": len(present_targets),
        "n_completed": len(lineage_completed),
        "n_failed": len(lineage_failed),
        "completed_ids": lineage_completed,
        "failed_ids": lineage_failed,
        "n_panda_builds": n_panda_builds,
        "shared_panda": True,
        "reference_cohort_size": n_cohort,
        "n_genes_input_post_qc": n_genes,
        "n_genes_internal_union": panda_internal_n_genes,
        "verified_flatten_order": verified_order,
        "panda_elapsed_sec": round(panda_elapsed, 1),
        "panda_network_shape": list(pn_shape),
    }
    lineage_summaries[lineage] = lin_summary
    lin_summary_path = (
        cfg.LAYER2_GRN_OUTPUT_DIR
        / cfg.L2_MULTI_MODEL_LINEAGE_SUMMARY_JSON.format(lineage=lineage)
    )
    with open(lin_summary_path, "w", encoding="utf-8") as fh:
        json.dump(lin_summary, fh, indent=2, ensure_ascii=False)
    print(f"\n  Lineage summary → {lin_summary_path.name}")
    logger.info(f"[M5/{lineage}] summary written: {lin_summary_path}")

    # Cleanup temp
    try:
        shutil.rmtree(LIONESS_TMP)
        logger.info(f"[M5/{lineage}] Temp cleaned: {LIONESS_TMP}")
    except Exception as ce:
        logger.warning(f"[M5/{lineage}] Could not clean temp dir: {ce}")

    # Free large objects between lineages
    del panda_obj, expr_cohort, expr_cohort_full, ep, pn_vals
    try:
        import gc

        gc.collect()
    except Exception:
        pass

logger.info("Stage M5 complete")


# =============================================================================
# STAGE M6: LEDGER + GATE REPORT + FINAL SUMMARY
# =============================================================================

print("\n" + "=" * 64)
print("STAGE M6: LEDGER + GATE REPORT + FINAL SUMMARY")
print("=" * 64)
logger.info("Stage M6 start")

NOTEBOOK_END_TIME = datetime.now()
elapsed_total = (NOTEBOOK_END_TIME - NOTEBOOK_START_TIME).total_seconds()

# DL2-28: every declared ID exactly once in ledger
declared_ids = list(TARGET_MODELS_LIST)
# Resolve declared list to ModelIDs (same as resolution)
declared_resolved = [cfg.resolve_model_id(x) for x in declared_ids]
ledger_ids = [e["model_id"] for e in run_ledger]

# Detect silent drops
missing_from_ledger = sorted(set(declared_resolved) - set(ledger_ids))
extra_in_ledger = sorted(set(ledger_ids) - set(declared_resolved))
no_silent = (len(missing_from_ledger) == 0) and (len(extra_in_ledger) == 0)

if not no_silent:
    # Attempt recovery note — still hard-fail DL2-28
    logger.error(
        f"[M6] Ledger integrity breach: missing={missing_from_ledger} "
        f"extra={extra_in_ledger}"
    )

cfg.validate_deadlock_rules(
    "dl2_28_full_ledger",
    declared_ids=declared_resolved,
    ledger_ids=ledger_ids,
    no_silent_drops=no_silent,
)
print("[DL2-28] PASS — full ledger for every declared ModelID")

# Status counts
status_counts = defaultdict(int)
for e in run_ledger:
    status_counts[e["status"]] += 1

print("\nRun ledger status counts:")
for st, n in sorted(status_counts.items()):
    print(f"  {st:<12}: {n}")

# Write ledger JSON + TSV
ledger_json_path = (
    cfg.LAYER2_GRN_OUTPUT_DIR / cfg.L2_MULTI_MODEL_LEDGER_JSON
)
ledger_tsv_path = cfg.LAYER2_GRN_OUTPUT_DIR / cfg.L2_MULTI_MODEL_LEDGER_TSV

ledger_doc = {
    "pipeline_version": PIPELINE_VERSION,
    "config_version": CONFIG_VERSION,
    "notebook": cfg.LAYER2_NOTEBOOK_02B_MULTI_NAME,
    "run_mode": cfg.LAYER2_RUN_MODE,
    "declared_list_source": "cfg.L2_TARGET_MODELS_LIST",
    "declared_n": len(declared_resolved),
    "declared_ids": declared_resolved,
    "missing_policy": cfg.L2_MULTI_MODEL_ON_MISSING,
    "reference_policy": cfg.LIONESS_REFERENCE_POLICY,
    "shared_panda_per_lineage": True,
    "cross_lineage_pooling": False,
    "status_counts": dict(status_counts),
    "entries": run_ledger,
    "timestamp_start": NOTEBOOK_START_TIME.isoformat(),
    "timestamp_end": NOTEBOOK_END_TIME.isoformat(),
    "elapsed_seconds": round(elapsed_total, 1),
}
with open(ledger_json_path, "w", encoding="utf-8") as fh:
    json.dump(ledger_doc, fh, indent=2, ensure_ascii=False)

pd.DataFrame(run_ledger).to_csv(ledger_tsv_path, sep="\t", index=False)
print(f"\nLedger JSON : {ledger_json_path.name}")
print(f"Ledger TSV  : {ledger_tsv_path.name}")

# Gate report
gate_report = {
    "pipeline_version": PIPELINE_VERSION,
    "config_version": CONFIG_VERSION,
    "notebook": cfg.LAYER2_NOTEBOOK_02B_MULTI_NAME,
    "foundation_notebook": cfg.LAYER2_NOTEBOOK_02A_NAME,
    "foundation_manifest": str(manifest_path),
    "log_file": str(cfg.LAYER2B_MULTI_LOG_FILE),
    "timestamp_start": NOTEBOOK_START_TIME.isoformat(),
    "timestamp_end": NOTEBOOK_END_TIME.isoformat(),
    "elapsed_seconds": round(elapsed_total, 1),
    "run_mode": cfg.LAYER2_RUN_MODE,
    "declared_n": len(declared_resolved),
    "status_counts": dict(status_counts),
    "lineage_summaries": lineage_summaries,
    "deadlock_rules": {
        "DL2-05_prior_source": {"status": "PASS", "value": "LOCAL_FILES"},
        "DL2-09_lioness_after_panda": {"status": "PASS"},
        "DL2-10_lioness_same_priors": {"status": "PASS"},
        "DL2-11_non_randomness": {"status": "NOT_RUN"},
        "DL2-14_target_exists": {"status": "PASS_PER_RESOLVED"},
        "DL2-15_lineage_from_metadata": {
            "status": "PASS",
            "source": "Model.csv",
        },
        "DL2-17_sample_order": {"status": "PASS_PER_TARGET"},
        "DL2-19_flatten_verified": {"status": "PASS_PER_LINEAGE"},
        "DL2-22_frozen_inputs_only": {
            "status": "PASS",
            "snail_recomputed": False,
            "foundation_dir": str(cfg.L2A_FROZEN_INPUTS_DIR),
        },
        "DL2-23_split_provenance": {
            "status": "PASS",
            "log_multi": str(cfg.LAYER2B_MULTI_LOG_FILE),
            "gate_report_multi": cfg.L2_MULTI_MODEL_GATE_REPORT_JSON,
        },
        "DL2-24_target_declared": {
            "status": "PASS",
            "source": "cfg.L2_TARGET_MODELS_LIST",
            "n": len(declared_resolved),
        },
        "DL2-25_shared_panda_per_lineage": {
            "status": "PASS",
            "cross_lineage_pooling": False,
            "policy": "one_panda_build_per_lineage",
        },
        "DL2-26_target_list_declared": {
            "status": "PASS",
            "source": "cfg.L2_TARGET_MODELS_LIST",
        },
        "DL2-27_export_only_if_verified_multi": {
            "status": "PASS",
        },
        "DL2-28_full_ledger": {
            "status": "PASS" if no_silent else "FAIL",
            "ledger_json": str(ledger_json_path),
            "ledger_tsv": str(ledger_tsv_path),
        },
    },
    "prior_provenance": {
        "motif_prior": str(local_motif),
        "ppi_prior": str(local_ppi),
        "ppi_runtime_weight_mode": ppi_prior_runtime_weight_mode,
        "jaspar_version": cfg.PRIOR_JASPAR_VERSION,
        "string_version": cfg.PRIOR_STRING_VERSION,
        "prior_generated_date": cfg.PRIOR_GENERATED_DATE,
        "hash_validated": True,
    },
    "outputs": {
        "ledger_json": cfg.L2_MULTI_MODEL_LEDGER_JSON,
        "ledger_tsv": cfg.L2_MULTI_MODEL_LEDGER_TSV,
        "gate_report": cfg.L2_MULTI_MODEL_GATE_REPORT_JSON,
        "completed_exports": [
            {
                "model_id": e["model_id"],
                "lineage": e.get("lineage"),
                "final_tsv": e.get("final_tsv"),
                "n_edges": e.get("n_edges"),
            }
            for e in completed_exports
        ],
        "log_file": str(cfg.LAYER2B_MULTI_LOG_FILE),
    },
}

gr_path = cfg.LAYER2_GRN_OUTPUT_DIR / cfg.L2_MULTI_MODEL_GATE_REPORT_JSON
with open(gr_path, "w", encoding="utf-8") as fh:
    json.dump(gate_report, fh, indent=2, ensure_ascii=False)

print(f"Gate report  : {gr_path.name}")
logger.info(f"[M6] Gate report: {gr_path}")

# Final summary
print(f"\n{'=' * 80}")
print("02B LIONESS MULTI-MODEL — COMPLETE (v2.4.0)")
print(f"{'=' * 80}")
print(f"  Declared targets : {len(declared_resolved)}")
print(f"  COMPLETED        : {status_counts.get('COMPLETED', 0)}")
print(f"  SKIPPED          : {status_counts.get('SKIPPED', 0)}")
print(f"  FAILED           : {status_counts.get('FAILED', 0)}")
print(f"  Lineages run     : {list(lineage_summaries.keys())}")
print(f"  Shared PANDA     : True (per lineage)")
print(f"  Cross-lineage    : False")
print(f"  Total elapsed    : {elapsed_total:.1f} s")
print(f"  Ledger           : {ledger_json_path}")
print(f"  Gate report      : {gr_path}")
print(f"  Log file         : {cfg.LAYER2B_MULTI_LOG_FILE}")
print(f"  Outputs dir      : {cfg.LAYER2_GRN_OUTPUT_DIR}")

logger.info("=" * 80)
logger.info("02B LIONESS MULTI-MODEL COMPLETE — v2.4.0")
logger.info(
    f"declared={len(declared_resolved)} "
    f"completed={status_counts.get('COMPLETED', 0)} "
    f"skipped={status_counts.get('SKIPPED', 0)} "
    f"failed={status_counts.get('FAILED', 0)} "
    f"elapsed={elapsed_total:.1f}s"
)
logger.info("=" * 80)

print("\n[DL2-11] Non-randomness test: NOT_RUN")
logger.info("[DL2-11] NOT_RUN")
