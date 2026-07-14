# %% [markdown]
# # 03 — Layer 3: RWR-MultiXrank Batch (single-file Jupyter)
# **Science:** v3.3.0 (geom-mean strict · DL3_07 CLASS_SHARED · 9c/9d L4 handoff)  
# **Architecture:** v4.0.0 batch over N LIONESS samples (e.g. 35/70 from L2 multi-model)  
# **Engine:** MultiXrank 0.3 — https://github.com/anthbapt/multixrank  
# **Integrity:** DL3_01…DL3_07 · frozen L3 PPI · no Stage-9b `all_results` / max_score
#
# ## Run in Jupyter
# ```python
# %cd "/home/labhhc5/Documents/workspace/D21/Duong Huy/pyscenic"
# # edit USER KNOBS in Stage U below, then:
# %run -i notebooks/03_Layer3_RWR_BATCH_JUPYTER.py
# ```

# %% [markdown]
# ## Stage U — User knobs (edit before run)

# %%
from __future__ import annotations

# ---------------------------------------------------------------------------
# USER KNOBS
# ---------------------------------------------------------------------------
# Sample discovery:
#   "l2_multimodel_ledger" — COMPLETED rows from L2 multi-model ledger (default)
#   "explicit_list"        — use MODELS + LINEAGE
#   "glob"                 — discover Z_*_LIONESS.tsv under LIONESS_{LINEAGE}/
SOURCE = "l2_multimodel_ledger"

LINEAGE = "Breast"  # None = all lineages in ledger / glob

# Smoke test: set e.g. ["ACH-000019"] then None for full ledger batch
MODELS = None  # list[str] | None

KEEP_WORKDIR = False  # True → snapshot MultiXrank workdir per sample
RUN_LAMBDA_BENCHMARK = True  # Stage 10 on first completed sample only
BATCH_ID = None  # None → auto timestamp batch_YYYYMMDDTHHMMSS

# Optional absolute project root if cwd is wrong
PROJECT_ROOT_OVERRIDE = None
# PROJECT_ROOT_OVERRIDE = "/home/labhhc5/Documents/workspace/D21/Duong Huy/pyscenic"

# %% [markdown]
# ## Stage 0 — Bootstrap + imports + config

# %%
import hashlib
import json
import os
import re
import shutil
import sys
import tempfile
import warnings
from concurrent.futures import ProcessPoolExecutor  # reserved; default sequential
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Any, Iterable

import numpy as np
import pandas as pd
import yaml
from scipy.stats import spearmanr


def bootstrap_project_root(override: str | None = None) -> Path:
    """Locate config_system.py without requiring __file__ (Jupyter-safe)."""
    candidates: list[Path] = []
    if override:
        candidates.append(Path(override))
    try:
        here = Path(__file__).resolve().parent  # type: ignore[name-defined]
        candidates += [here, here.parent, here.parent.parent]
    except NameError:
        pass
    cwd = Path.cwd().resolve()
    candidates += [cwd, cwd.parent]
    candidates.append(
        Path("/home/labhhc5/Documents/workspace/D21/Duong Huy/pyscenic")
    )
    for p in list(sys.path):
        if p:
            candidates.append(Path(p))
    seen: set[Path] = set()
    for c in candidates:
        try:
            r = c.resolve()
        except Exception:
            continue
        if r in seen:
            continue
        seen.add(r)
        if (r / "config_system.py").is_file():
            if str(r) not in sys.path:
                sys.path.insert(0, str(r))
            os.chdir(r)
            return r
    raise FileNotFoundError(
        "Cannot find config_system.py.\n"
        "In Jupyter: %cd to the pyscenic project root, then re-run.\n"
        f"cwd={Path.cwd()}"
    )


_PROJECT_ROOT = bootstrap_project_root(PROJECT_ROOT_OVERRIDE)
print(f"[M0] project root: {_PROJECT_ROOT}")

import config_system as cfg

import multixrank

_MX_VER = getattr(multixrank, "__version__", "unknown")
if _MX_VER != cfg.L3_MULTIXRANK_REQUIRED_VERSION:
    warnings.warn(
        f"MultiXrank {_MX_VER} != required {cfg.L3_MULTIXRANK_REQUIRED_VERSION}",
        stacklevel=1,
    )

logger = cfg.setup_logger(
    name="Layer3_RWR_Batch",
    logfile=getattr(cfg, "LAYER3_BATCH_LOG_FILE", cfg.LAYER3_LOG_FILE),
    reset_handlers=True,
)

NOTEBOOK_NAME = getattr(
    cfg, "LAYER3_NOTEBOOK_BATCH_NAME", "03_Layer3_RWR_BATCH_JUPYTER.py"
)
SCIENCE_VERSION = "v3.3.0"
ARCH_VERSION = "v4.0.0-monolith"

logger.info("=" * 72)
logger.info("LAYER 3 BATCH — science %s | arch %s", SCIENCE_VERSION, ARCH_VERSION)
logger.info("MultiXrank: %s | notebook: %s", _MX_VER, NOTEBOOK_NAME)
logger.info("cwd: %s", Path.cwd())
print(f"Stage 0: OK | MultiXrank {_MX_VER} | science {SCIENCE_VERSION}")

# %% [markdown]
# ## Stage 0b — Domain helpers (integrity core)

# %%


def normalize_gene_name(name: str | None) -> str:
    if name is None:
        return ""
    return str(name).replace("-", "").strip().upper()


def sanitize_node_name(value: str) -> str:
    txt = re.sub(r"[^A-Za-z0-9_.-]+", "_", str(value).strip())
    txt = re.sub(r"_+", "_", txt).strip("_")
    return txt or "CLASS_SHARED"


def sha256_of_file(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(65536), b""):
            h.update(chunk)
    return h.hexdigest()


def geometric_mean_strict(values: np.ndarray) -> np.ndarray:
    """Any layer score <= 0 ⇒ global score 0 (v3.3.0 multiplex penalty)."""
    arr = np.asarray(values, dtype=np.float64)
    if arr.ndim != 2:
        raise ValueError(f"Expected 2D, got {arr.shape}")
    out = np.zeros(arr.shape[0], dtype=np.float64)
    ok = np.all(arr > 0, axis=1)
    if np.any(ok):
        out[ok] = np.exp(np.mean(np.log(arr[ok]), axis=1))
    return out


def layer_name_to_role(layer_name: str) -> str:
    txt = str(layer_name).upper()
    if "PPI" in txt or txt in {"1", "01"}:
        return "PPI"
    if "GRN" in txt or txt in {"2", "11"}:
        return "GRN"
    raise ValueError(f"Unrecognized layer label: {layer_name!r}")


def build_global_gene_ranking_from_per_layer(
    per_layer_df: pd.DataFrame,
) -> pd.DataFrame:
    """Long gene×layer → wide + rwr_score_geom_mean (DL3_05 preserves long table)."""
    required = {"gene_symbol", "layer", "rwr_score"}
    missing = required - set(per_layer_df.columns)
    if missing:
        raise ValueError(f"per_layer_df missing: {sorted(missing)}")

    tmp = per_layer_df.copy()
    tmp["_role"] = tmp["layer"].map(layer_name_to_role)
    pivot = tmp.pivot_table(
        index="gene_symbol",
        columns="_role",
        values="rwr_score",
        aggfunc="first",
    ).reset_index()
    pivot.columns.name = None
    for role in ("PPI", "GRN"):
        if role not in pivot.columns:
            pivot[role] = 0.0
    pivot["PPI"] = pivot["PPI"].fillna(0.0).astype(float)
    pivot["GRN"] = pivot["GRN"].fillna(0.0).astype(float)
    scores = pivot[["PPI", "GRN"]].to_numpy(dtype=np.float64)
    pivot["rwr_score_geom_mean"] = geometric_mean_strict(scores)
    pivot["layer_support_min"] = scores.min(axis=1)
    pivot["dominant_layer"] = np.select(
        [pivot["PPI"] > pivot["GRN"], pivot["GRN"] > pivot["PPI"]],
        ["PPI_layer", "GRN_layer"],
        default="TIE",
    )
    pivot = pivot.rename(
        columns={"PPI": "score_PPI_layer", "GRN": "score_GRN_layer"}
    )
    pivot = pivot.sort_values(
        by=[
            "rwr_score_geom_mean",
            "layer_support_min",
            "score_PPI_layer",
            "score_GRN_layer",
            "gene_symbol",
        ],
        ascending=[False, False, False, False, True],
    ).reset_index(drop=True)
    pivot.insert(0, "rank_geom_mean", np.arange(1, len(pivot) + 1))
    pivot["original_rank"] = pivot["rank_geom_mean"]
    return pivot


def transform_grn_weights(
    grn_df: pd.DataFrame,
    source_col: str,
    target_col: str,
    weight_col: str,
) -> tuple[pd.DataFrame, dict]:
    w = grn_df[weight_col]
    stats = {
        "n_edges_before": int(len(grn_df)),
        "n_negative": int((w < 0).sum()),
        "n_zero": int((w == 0).sum()),
        "n_positive": int((w > 0).sum()),
        "weight_min_before": float(w.min()) if len(grn_df) else None,
        "weight_max_before": float(w.max()) if len(grn_df) else None,
        "transform_method": cfg.L3_GRN_WEIGHT_TRANSFORM,
    }
    out = grn_df[grn_df[weight_col] > 0].copy()
    stats["n_edges_after"] = int(len(out))
    stats["n_edges_removed"] = stats["n_edges_before"] - stats["n_edges_after"]
    stats["pct_removed"] = round(
        100 * stats["n_edges_removed"] / max(stats["n_edges_before"], 1), 2
    )
    if len(out) == 0:
        raise ValueError("GRN has 0 edges after clip_nonpositive.")
    stats["weight_min_after"] = float(out[weight_col].min())
    stats["weight_max_after"] = float(out[weight_col].max())
    return out, stats


def write_network_file(edges: list[tuple], path: Path, weighted: bool) -> None:
    lines = []
    for e in edges:
        if weighted:
            lines.append(f"{e[0]}\t{e[1]}\t{float(e[2])}")
        else:
            lines.append(f"{e[0]}\t{e[1]}")
    path.write_text("\n".join(lines), encoding="utf-8")


def write_drug_monoplex(drug_ids: list[str], path: Path) -> None:
    if not drug_ids:
        raise ValueError("drug_ids empty")
    path.write_text(
        "\n".join(f"{d}\t{d}\t1.0" for d in drug_ids) + "\n", encoding="utf-8"
    )


def write_seed_file(seed_genes: list[str], path: Path) -> None:
    if not seed_genes:
        raise ValueError("seed_genes empty")
    path.write_text("\n".join(seed_genes) + "\n", encoding="utf-8")


def build_multixrank_config(
    wdir: Path,
    seed_path: Path,
    ppi_path: Path,
    grn_path: Path,
    drug_path: Path,
    bipartite_path: Path,
    lambda_val: float,
) -> dict:
    def rel(p: Path) -> str:
        return str(p.relative_to(wdir))

    return {
        "seed": rel(seed_path),
        "r": cfg.L3_RWR_RESTART_PROB,
        "eta": list(cfg.L3_ETA_DEFAULT),
        "lamb": [
            [lambda_val, 1 - lambda_val],
            [1 - lambda_val, lambda_val],
        ],
        "multiplex": {
            cfg.L3_GENE_MULTIPLEX_ID: {
                "layers": [rel(ppi_path), rel(grn_path)],
                "delta": cfg.L3_DELTA_GENE,
                "graph_type": [
                    cfg.L3_GRAPH_TYPE_WEIGHTED_UNDIRECTED,
                    cfg.L3_GRAPH_TYPE_WEIGHTED_DIRECTED,
                ],
                "tau": list(cfg.L3_TAU_GENE_DEFAULT),
            },
            cfg.L3_DRUG_MULTIPLEX_ID: {
                "layers": [rel(drug_path)],
                "delta": cfg.L3_DELTA_DRUG,
                "graph_type": [cfg.L3_GRAPH_TYPE_WEIGHTED_UNDIRECTED],
                "tau": list(cfg.L3_TAU_DRUG_DEFAULT),
            },
        },
        "bipartite": {
            rel(bipartite_path): {
                "source": cfg.L3_GENE_MULTIPLEX_ID,
                "target": cfg.L3_DRUG_MULTIPLEX_ID,
                "graph_type": cfg.L3_GRAPH_TYPE_UNWEIGHTED_UNDIRECTED,
            }
        },
    }


def run_rwr(wdir: Path, config_path: Path) -> pd.DataFrame:
    obj = multixrank.Multixrank(config=str(config_path), wdir=str(wdir))
    result = obj.random_walk_rank()
    expected = ["multiplex", "node", "layer", "score"]
    if list(result.columns) != expected:
        raise RuntimeError(
            f"MultiXrank schema changed: {list(result.columns)} != {expected}"
        )
    return result


def extract_per_layer_scores(
    ranking_df: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    gene_mask = ranking_df["multiplex"].astype(str) == str(cfg.L3_GENE_MULTIPLEX_ID)
    drug_mask = ranking_df["multiplex"].astype(str) == str(cfg.L3_DRUG_MULTIPLEX_ID)
    gene_df = (
        ranking_df.loc[gene_mask, ["node", "layer", "score"]]
        .copy()
        .rename(columns={"node": "gene_symbol", "score": "rwr_score"})
        .reset_index(drop=True)
    )
    drug_df = (
        ranking_df.loc[drug_mask, ["node", "layer", "score"]]
        .copy()
        .rename(columns={"node": "drug_id", "score": "rwr_score"})
        .sort_values("rwr_score", ascending=False)
        .reset_index(drop=True)
    )
    n_neg = int((gene_df["rwr_score"] < 0).sum())
    if n_neg:
        logger.warning("INTEGRITY: %d negative gene scores", n_neg)
    total = float(ranking_df["score"].sum())
    if abs(total - 1.0) > 0.01:
        logger.warning("INTEGRITY: score sum=%.6f (expected ~1)", total)
    return gene_df, drug_df


def load_housekeeping_gmt(gmt_path: Path) -> tuple[set[str], dict]:
    if not gmt_path.exists():
        raise FileNotFoundError(f"Housekeeping GMT not found: {gmt_path}")
    genes: set[str] = set()
    n_valid = 0
    with open(gmt_path, "r", encoding="utf-8") as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            n_valid += 1
            for g in parts[2:]:
                ng = normalize_gene_name(g)
                if ng:
                    genes.add(ng)
    if not genes:
        raise ValueError(f"GMT empty: {gmt_path}")
    meta = {
        "file": str(gmt_path),
        "file_name": gmt_path.name,
        "sha256": sha256_of_file(gmt_path),
        "n_valid_lines": n_valid,
        "n_unique_hk_genes_norm": len(genes),
    }
    return genes, meta


def build_layer4_ready_nohk(
    global_ranking: pd.DataFrame,
    hk_set: set[str],
    top_n: int,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, dict]:
    lb = global_ranking.copy()
    lb["gene_symbol"] = lb["gene_symbol"].astype(str).str.strip()
    lb["gene_symbol_norm"] = lb["gene_symbol"].map(normalize_gene_name)
    if lb["gene_symbol_norm"].duplicated().any():
        raise ValueError("Normalization collision in global ranking.")
    lb["is_housekeeping"] = lb["gene_symbol_norm"].isin(hk_set)
    hk_rows = lb[lb["is_housekeeping"]].copy()
    non_hk = lb[~lb["is_housekeeping"]].copy()
    non_hk = non_hk.sort_values(
        ["rwr_score_geom_mean", "layer_support_min", "rank_geom_mean"],
        ascending=[False, False, True],
    ).reset_index(drop=True)
    non_hk["rank_noHK"] = np.arange(1, len(non_hk) + 1)

    cols = [
        "rank_noHK",
        "rank_geom_mean",
        "gene_symbol",
        "gene_symbol_norm",
        "score_GRN_layer",
        "score_PPI_layer",
        "rwr_score_geom_mean",
        "layer_support_min",
        "dominant_layer",
    ]
    ora = non_hk.head(top_n)[cols].rename(columns={"rank_geom_mean": "original_rank"})
    full_nohk = non_hk[cols].rename(columns={"rank_geom_mean": "original_rank"})
    hk_export = hk_rows[
        [
            "rank_geom_mean",
            "gene_symbol",
            "gene_symbol_norm",
            "score_GRN_layer",
            "score_PPI_layer",
            "rwr_score_geom_mean",
            "layer_support_min",
            "dominant_layer",
        ]
    ].rename(columns={"rank_geom_mean": "original_rank"})
    qc = {
        "n_total": int(len(lb)),
        "n_hk": int(len(hk_rows)),
        "n_non_hk": int(len(non_hk)),
        "n_exported_top": int(len(ora)),
        "n_hk_in_top100": int((hk_rows["rank_geom_mean"] <= 100).sum())
        if len(hk_rows)
        else 0,
    }
    return ora, full_nohk, hk_export, qc


def build_gene_universe_artifact(
    non_hk_full: pd.DataFrame, hk_export: pd.DataFrame
) -> tuple[pd.DataFrame, str]:
    a = non_hk_full[["gene_symbol", "gene_symbol_norm"]].copy()
    a["source_partition"] = "NON_HK"
    b = hk_export[["gene_symbol", "gene_symbol_norm"]].copy()
    b["source_partition"] = "HK"
    bg = pd.concat([a, b], ignore_index=True)
    if bg["gene_symbol_norm"].duplicated().any():
        raise ValueError("Duplicate genes across NON_HK/HK partitions.")
    n_diff = int(
        (bg["gene_symbol"].astype(str) != bg["gene_symbol_norm"].astype(str)).sum()
    )
    canonical = "gene_symbol" if n_diff == 0 else "gene_symbol_norm"
    bg = bg.sort_values([canonical, "source_partition"]).reset_index(drop=True)
    return bg, canonical


def verify_seed_equivalence(
    valid_drug_seed_map: dict[str, list[str]],
    representative_drug_id: str,
) -> dict:
    """DL3_07 structural audit: identical effective seed lists."""
    ref = tuple(sorted(valid_drug_seed_map[representative_drug_id]))
    bad = {
        d: sorted(set(ref) ^ set(s))
        for d, s in valid_drug_seed_map.items()
        if tuple(sorted(s)) != ref
    }
    ok = len(bad) == 0
    if not ok:
        raise ValueError(
            "DL3_07: compounds not mathematically equivalent for single-run RWR. "
            f"Non-identical seed sets (examples): {list(bad)[:5]}"
        )
    return {
        "verification_method": "structural_equivalence_audit",
        "seed_representation_in_rwr": "unweighted seed gene list + bipartite edges",
        "n_compounds_verified": len(valid_drug_seed_map),
        "representative_drug_id": representative_drug_id,
        "seed_sets_identical": True,
        "shared_seed_genes": list(ref),
        "mathematically_equivalent": True,
        "same_ppi_network": True,
        "same_grn_network": True,
        "same_hyperparameters": True,
        "interpretation": (
            "All valid compounds share the same effective seed gene list after "
            "universe intersection; one CLASS_SHARED RWR is exact under current impl."
        ),
    }


def sample_run_id(
    model_id: str, grn_sha: str | None, ppi_sha: str, lambda_val: float
) -> str:
    raw = (
        f"{model_id}|{grn_sha}|{ppi_sha}|{lambda_val}|"
        f"{cfg.L3_RWR_RESTART_PROB}|{cfg.L3_DELTA_GENE}|geom_mean|{ARCH_VERSION}"
    )
    return hashlib.sha256(raw.encode()).hexdigest()[:16]


print("Stage 0b: helpers defined.")

# %% [markdown]
# ## Stage 1 — Resource check

# %%
logger.info("STAGE 1: Resource check")
l3_resources = cfg.check_layer3_resources()
critical = {
    "multixrank_available",
    "yaml_available",
    "layer1_output_dir",
    "layer1_p0_long_csv",
    "layer2_grn_output_dir",
    "local_ppi_prior_l3_tsv",
    "local_ppi_prior_l3_manifest",
    "local_ppi_prior_l3_sha256_ok",
}
print("Layer 3 resource check:")
all_ok = True
for k, v in l3_resources.items():
    flag = "✓" if v is True else ("✗" if v is False else " ")
    crit = " ← CRITICAL" if k in critical and v is False else ""
    print(f"  [{flag}] {k}: {v}{crit}")
    if k in critical and v is False:
        all_ok = False
if not all_ok:
    raise RuntimeError("Critical Layer 3 resources missing or SHA mismatch.")
print("Stage 1: PASSED")

# %% [markdown]
# ## Stage 2 — Load Layer 1 seeds (once)

# %%
logger.info("STAGE 2: Load Layer 1 seeds")
p0_path = cfg.LAYER1_OUTPUT_DIR / cfg.L1_P0_VECTOR_LONG_CSV
if not p0_path.exists():
    raise FileNotFoundError(p0_path)
p0_long = pd.read_csv(p0_path)
P0_DRUG_COL = getattr(cfg, "COL_LIGAND_ID", "lig_id")
P0_GENE_COL = getattr(cfg, "COL_TARGET", "target")
P0_IN_TD_COL = "in_T_d"
for col in (P0_DRUG_COL, P0_GENE_COL):
    if col not in p0_long.columns:
        raise ValueError(f"P0 missing {col}: {p0_long.columns.tolist()}")
if P0_IN_TD_COL in p0_long.columns:
    p0_long = p0_long[p0_long[P0_IN_TD_COL] == True].copy()  # noqa: E712
p0_long[P0_GENE_COL] = p0_long[P0_GENE_COL].map(normalize_gene_name)
cfg.validate_deadlock_rules("dl3_01_p0_no_reweight", reweighted=False)

drug_seed_map: dict[str, list[str]] = {}
for drug_id, grp in p0_long.groupby(P0_DRUG_COL):
    seeds = sorted(grp[P0_GENE_COL].dropna().unique().tolist())
    if seeds:
        drug_seed_map[str(drug_id)] = seeds
if not drug_seed_map:
    raise ValueError("No seeds from Layer 1.")
print(f"Stage 2: {len(drug_seed_map)} drugs with seeds")

# %% [markdown]
# ## Stage 3 — Load frozen L3 PPI (once, DL3_06)

# %%
logger.info("STAGE 3: Load L3 PPI")
ppi_path = Path(cfg.LOCAL_PPI_PRIOR_L3_TSV)
ppi_manifest_path = Path(cfg.LOCAL_PPI_PRIOR_L3_MANIFEST_JSON)
if not ppi_path.exists() or not ppi_manifest_path.exists():
    raise FileNotFoundError("L3 PPI artifact or manifest missing.")

actual_sha256 = sha256_of_file(ppi_path)
sha256_ok = actual_sha256 == cfg.LOCAL_PPI_PRIOR_L3_EXPECTED_SHA256
if not sha256_ok:
    raise ValueError(
        f"SHA-256 mismatch for {ppi_path.name}\n"
        f"  expected {cfg.LOCAL_PPI_PRIOR_L3_EXPECTED_SHA256}\n"
        f"  actual   {actual_sha256}"
    )
with open(ppi_manifest_path, "r", encoding="utf-8") as fh:
    ppi_manifest = json.load(fh)
params = ppi_manifest.get("parameters", {})
tf_restricted = params.get("tf_universe_restriction", True)
if tf_restricted is not False:
    raise ValueError(f"tf_universe_restriction must be False, got {tf_restricted!r}")
for key, expected in [
    ("string_version", cfg.L3_PPI_STRING_VERSION),
    ("network_type", cfg.L3_PPI_NETWORK_TYPE),
]:
    if str(params.get(key)) != str(expected):
        raise ValueError(f"PPI manifest {key}: {params.get(key)!r} != {expected!r}")
if int(params.get("min_combined_score", -1)) != cfg.L3_PPI_MIN_STRING_SCORE:
    raise ValueError("PPI min_combined_score mismatch")

cfg.validate_deadlock_rules(
    "dl3_06_layer3_ppi_source",
    ppi_path_used=str(ppi_path),
    manifest_validated=True,
    tf_restricted=False,
    sha256_ok=True,
)

ppi_raw = pd.read_csv(ppi_path, sep="\t")
missing = [c for c in cfg.L3_PPI_REQUIRED_COLUMNS if c not in ppi_raw.columns]
if missing:
    raise ValueError(f"PPI missing columns {missing}")
ppi_raw = ppi_raw[cfg.L3_PPI_REQUIRED_COLUMNS].copy()
ppi_raw.columns = ["gene_a", "gene_b", "weight"]
ppi_raw["gene_a"] = ppi_raw["gene_a"].map(normalize_gene_name)
ppi_raw["gene_b"] = ppi_raw["gene_b"].map(normalize_gene_name)
ppi_raw["weight"] = pd.to_numeric(ppi_raw["weight"], errors="raise")
ppi_raw = ppi_raw[ppi_raw["gene_a"] != ppi_raw["gene_b"]].copy()
if (ppi_raw["weight"] <= 0).any():
    raise ValueError("Non-positive PPI weights.")
ppi_genes = sorted(set(ppi_raw["gene_a"]) | set(ppi_raw["gene_b"]))
ppi_gene_set = set(ppi_genes)
ppi_edges = list(
    ppi_raw[["gene_a", "gene_b", "weight"]].itertuples(index=False, name=None)
)
print(
    f"Stage 3: PPI edges={len(ppi_edges):,} genes={len(ppi_genes):,} "
    f"SHA={actual_sha256[:16]}... OK"
)

# %% [markdown]
# ## Stage S — Resolve sample list (batch axis = LIONESS cell lines)

# %%
logger.info("STAGE S: Resolve samples")


@dataclass(frozen=True)
class SampleSpec:
    model_id: str
    lineage: str
    grn_tsv: Path
    grn_sha256: str | None = None


def resolve_samples(
    source: str,
    lineage: str | None,
    models: list[str] | None,
) -> list[SampleSpec]:
    source = source or getattr(cfg, "L3_BATCH_SOURCE", "l2_multimodel_ledger")
    samples: list[SampleSpec] = []

    if source == "explicit_list":
        ids = models or list(getattr(cfg, "L3_SAMPLE_MODELS_LIST", None) or [])
        if not ids or not lineage:
            raise ValueError("explicit_list needs MODELS/L3_SAMPLE_MODELS_LIST + LINEAGE")
        for mid in ids:
            tsv = (
                cfg.LAYER2_GRN_OUTPUT_DIR
                / cfg.L2_GRN_LIONESS_DIR.format(lineage=lineage)
                / cfg.L2_GRN_LIONESS_TSV.format(sample=mid)
            )
            samples.append(
                SampleSpec(
                    model_id=mid,
                    lineage=lineage,
                    grn_tsv=tsv,
                    grn_sha256=sha256_of_file(tsv) if tsv.exists() else None,
                )
            )
        return samples

    if source == "glob":
        root = cfg.LAYER2_GRN_OUTPUT_DIR
        dirs = (
            [root / cfg.L2_GRN_LIONESS_DIR.format(lineage=lineage)]
            if lineage
            else sorted(root.glob("LIONESS_*"))
        )
        for d in dirs:
            if not d.is_dir():
                continue
            lin = d.name.replace("LIONESS_", "", 1)
            for tsv in sorted(d.glob("Z_*_LIONESS.tsv")):
                mid = tsv.name[2 : -len("_LIONESS.tsv")]
                if models and mid not in models:
                    continue
                samples.append(
                    SampleSpec(
                        model_id=mid,
                        lineage=lin,
                        grn_tsv=tsv,
                        grn_sha256=sha256_of_file(tsv),
                    )
                )
        return samples

    # l2_multimodel_ledger
    ledger_path = cfg.LAYER2_GRN_OUTPUT_DIR / getattr(
        cfg, "L2_MULTI_MODEL_LEDGER_JSON", "L2B_LIONESS_MultiModel_Run_Ledger.json"
    )
    if not ledger_path.exists():
        # fallback TSV
        ledger_path = cfg.LAYER2_GRN_OUTPUT_DIR / getattr(
            cfg, "L2_MULTI_MODEL_LEDGER_TSV", "L2B_LIONESS_MultiModel_Run_Ledger.tsv"
        )
    if not ledger_path.exists():
        raise FileNotFoundError(
            f"L2 multi-model ledger not found under {cfg.LAYER2_GRN_OUTPUT_DIR}. "
            "Run 02B multi-model first, or use SOURCE='glob' / 'explicit_list'."
        )
    if ledger_path.suffix == ".json":
        with open(ledger_path, "r", encoding="utf-8") as fh:
            doc = json.load(fh)
        rows = doc.get("entries", doc if isinstance(doc, list) else [])
        df = pd.DataFrame(rows)
    else:
        df = pd.read_csv(ledger_path, sep="\t")
    if "status" in df.columns:
        df = df[df["status"].astype(str).str.upper() == "COMPLETED"].copy()
    for _, r in df.iterrows():
        mid = str(r.get("model_id", "")).strip()
        lin = str(r.get("lineage", "")).strip() or "UNKNOWN"
        if lineage and lin != lineage:
            continue
        if models and mid not in models:
            continue
        tsv = r.get("final_tsv") or r.get("grn_tsv")
        if not mid or not tsv:
            continue
        tsv_path = Path(str(tsv))
        if not tsv_path.is_absolute():
            tsv_path = cfg.LAYER2_GRN_OUTPUT_DIR / tsv_path
        samples.append(
            SampleSpec(
                model_id=mid,
                lineage=lin,
                grn_tsv=tsv_path,
                grn_sha256=sha256_of_file(tsv_path) if tsv_path.exists() else None,
            )
        )
    return samples


SAMPLES = resolve_samples(SOURCE, LINEAGE, MODELS)
if not SAMPLES:
    raise RuntimeError("No samples resolved.")
missing_grn = [s for s in SAMPLES if not s.grn_tsv.exists()]
print(f"Stage S: {len(SAMPLES)} samples resolved | missing GRN files: {len(missing_grn)}")
for s in SAMPLES[:8]:
    print(f"  [{'✓' if s.grn_tsv.exists() else '✗'}] {s.lineage}/{s.model_id}")
if len(SAMPLES) > 8:
    print(f"  ... +{len(SAMPLES) - 8} more")
if missing_grn:
    raise FileNotFoundError(
        f"Missing GRN e.g. {missing_grn[0].grn_tsv} — fix L2 outputs first."
    )

# %% [markdown]
# ## Stage B — Batch loop: one CLASS_SHARED RWR per sample

# %%
logger.info("STAGE B: Batch RWR over %d samples", len(SAMPLES))

batch_id = BATCH_ID or datetime.now().strftime("batch_%Y%m%dT%H%M%S")
out_root = cfg.LAYER3_OUTPUT_DIR / batch_id
out_root.mkdir(parents=True, exist_ok=True)
print(f"Batch output root: {out_root}")

lambda_val = float(cfg.L3_LAMBDA_DEFAULT)
hk_set, hk_meta = load_housekeeping_gmt(Path(cfg.L3_HK_GMT_PATH))
print(f"HK genes loaded: {len(hk_set):,} from {hk_meta['file_name']}")

ledger_rows: list[dict[str, Any]] = []
first_completed_global: pd.DataFrame | None = None
first_completed_seeds: list[str] | None = None
first_completed_node: str | None = None
shared_ppi_for_bench: Path | None = None
shared_grn_for_bench: Path | None = None
bench_tmpdir: Path | None = None

batch_t0 = datetime.now()

for idx, sample in enumerate(SAMPLES, 1):
    print(f"\n[{idx}/{len(SAMPLES)}] {sample.lineage}/{sample.model_id}")
    logger.info("[%d/%d] %s", idx, len(SAMPLES), sample.model_id)
    row: dict[str, Any] = {
        "model_id": sample.model_id,
        "lineage": sample.lineage,
        "status": "FAILED",
        "reason": "",
        "sample_run_id": "",
        "final_top100": "",
        "n_universe": 0,
        "dl3_07_class_shared": False,
    }
    t_sample = datetime.now()
    try:
        # --- GRN ---
        grn_df = pd.read_csv(
            sample.grn_tsv,
            sep="\t",
            usecols=[cfg.COL_GRN_SOURCE, cfg.COL_GRN_TARGET, cfg.COL_GRN_WEIGHT],
        )
        grn_df[cfg.COL_GRN_SOURCE] = grn_df[cfg.COL_GRN_SOURCE].map(
            normalize_gene_name
        )
        grn_df[cfg.COL_GRN_TARGET] = grn_df[cfg.COL_GRN_TARGET].map(
            normalize_gene_name
        )
        grn_df[cfg.COL_GRN_WEIGHT] = pd.to_numeric(
            grn_df[cfg.COL_GRN_WEIGHT], errors="raise"
        )
        n_grn_raw = len(grn_df)
        grn_df = grn_df[
            grn_df[cfg.COL_GRN_SOURCE].isin(ppi_gene_set)
            & grn_df[cfg.COL_GRN_TARGET].isin(ppi_gene_set)
        ].copy()
        cfg.validate_deadlock_rules(
            "dl3_03_frozen_layer2",
            grn_recomputed=False,
            using_frozen_layer2=True,
        )
        grn_clipped, clip_stats = transform_grn_weights(
            grn_df,
            cfg.COL_GRN_SOURCE,
            cfg.COL_GRN_TARGET,
            cfg.COL_GRN_WEIGHT,
        )
        cfg.validate_deadlock_rules(
            "dl3_02_grn_clip_applied",
            provenance={"grn_weight_transform": clip_stats},
        )

        grn_genes = set(grn_clipped[cfg.COL_GRN_SOURCE]) | set(
            grn_clipped[cfg.COL_GRN_TARGET]
        )
        gene_universe = sorted(ppi_gene_set & grn_genes)
        gene_universe_set = set(gene_universe)
        if not gene_universe:
            raise ValueError("Empty gene universe.")

        # --- seeds in universe ---
        valid_drug_seed_map: dict[str, list[str]] = {}
        for drug_id, seeds in drug_seed_map.items():
            valid = sorted(set(seeds) & gene_universe_set)
            if valid:
                valid_drug_seed_map[drug_id] = valid
        if not valid_drug_seed_map:
            raise ValueError("No drugs retain valid seeds in universe.")

        rep_drug = sorted(valid_drug_seed_map.keys())[0]
        eq_report = verify_seed_equivalence(valid_drug_seed_map, rep_drug)
        seed_genes = list(eq_report["shared_seed_genes"])
        class_node = sanitize_node_name(
            f"{cfg.L3_CLASS_SHARED_NODE_PREFIX}__{sample.lineage}__{sample.model_id}"
        )
        cfg.validate_deadlock_rules(
            "dl3_07_single_representative_run",
            mathematically_equivalent=True,
            synthetic_class_node=True,
            n_runs_executed=1,
        )

        ppi_f = [
            (a, b, w)
            for a, b, w in ppi_edges
            if a in gene_universe_set and b in gene_universe_set
        ]
        grn_f = list(
            grn_clipped[
                grn_clipped[cfg.COL_GRN_SOURCE].isin(gene_universe_set)
                & grn_clipped[cfg.COL_GRN_TARGET].isin(gene_universe_set)
            ][
                [cfg.COL_GRN_SOURCE, cfg.COL_GRN_TARGET, cfg.COL_GRN_WEIGHT]
            ].itertuples(index=False, name=None)
        )
        if not ppi_f or not grn_f:
            raise ValueError("Empty PPI/GRN after universe filter.")

        # --- MultiXrank workdir ---
        wdir = Path(tempfile.mkdtemp(prefix=f"l3_{sample.model_id}_"))
        try:
            net = wdir / "network"
            seeds_dir = wdir / "seeds"
            net.mkdir(parents=True)
            seeds_dir.mkdir()
            ppi_file = net / "PPI.tsv"
            grn_file = net / "GRN.tsv"
            drug_file = net / "drug.tsv"
            bip_file = net / "bipartite.tsv"
            seed_file = seeds_dir / "seeds.txt"
            write_network_file(ppi_f, ppi_file, weighted=True)
            write_network_file(grn_f, grn_file, weighted=True)
            write_drug_monoplex([class_node], drug_file)
            cfg.validate_deadlock_rules(
                "dl3_04_drug_monoplex",
                drug_monoplex_policy=cfg.L3_DRUG_MONOPLEX_POLICY,
            )
            write_network_file(
                [(g, class_node) for g in seed_genes], bip_file, weighted=False
            )
            write_seed_file(seed_genes, seed_file)
            conf = build_multixrank_config(
                wdir, seed_file, ppi_file, grn_file, drug_file, bip_file, lambda_val
            )
            cfg_path = wdir / "config.yml"
            with open(cfg_path, "w", encoding="utf-8") as fh:
                yaml.dump(conf, fh, default_flow_style=False)

            t0 = datetime.now()
            ranking_df = run_rwr(wdir, cfg_path)
            rwr_sec = (datetime.now() - t0).total_seconds()
            gene_long, drug_df = extract_per_layer_scores(ranking_df)
            cfg.validate_deadlock_rules("dl3_05_no_aggregation", gene_df=gene_long)

            global_ranking = build_global_gene_ranking_from_per_layer(gene_long)
            ora, full_nohk, hk_export, hk_qc = build_layer4_ready_nohk(
                global_ranking, hk_set, int(cfg.L3_ORA_TOP_N_NO_HK)
            )
            if len(full_nohk) < int(cfg.L3_ORA_MIN_NONHK):
                raise ValueError(
                    f"Only {len(full_nohk)} non-HK genes < L3_ORA_MIN_NONHK"
                )
            bg_df, canonical = build_gene_universe_artifact(full_nohk, hk_export)
            bg_list = bg_df[canonical].astype(str).tolist()
            # query ⊂ background
            if set(ora["gene_symbol_norm"]) - set(bg_df["gene_symbol_norm"]):
                raise ValueError("top100_noHK not subset of gene universe.")

            run_id = sample_run_id(
                sample.model_id, sample.grn_sha256, actual_sha256, lambda_val
            )
            out_dir = (
                out_root
                / "CLASS_SHARED"
                / sample.lineage
                / sample.model_id
            )
            out_dir.mkdir(parents=True, exist_ok=True)
            prefix = f"lambda{lambda_val}"

            paths = {}
            p = out_dir / f"{prefix}_CLASS_SHARED_full_per_layer.tsv"
            gene_long.to_csv(p, sep="\t", index=False)
            paths["full_per_layer"] = p

            p = out_dir / f"{prefix}_CLASS_SHARED_global_geomean_all_genes.tsv"
            global_ranking.to_csv(p, sep="\t", index=False)
            paths["global_all"] = p

            p = out_dir / f"{prefix}_CLASS_SHARED_top100.tsv"
            global_ranking.head(100).to_csv(p, sep="\t", index=False)

            p = out_dir / f"{prefix}_CLASS_SHARED_top100_noHK.tsv"
            ora.to_csv(p, sep="\t", index=False)
            paths["top100_noHK"] = p

            p = out_dir / f"{prefix}_CLASS_SHARED_full_ranking_noHK.tsv"
            full_nohk.to_csv(p, sep="\t", index=False)

            p = out_dir / f"{prefix}_CLASS_SHARED_HK_in_ranking.tsv"
            hk_export.to_csv(p, sep="\t", index=False)

            p = out_dir / f"{prefix}_CLASS_SHARED_gene_universe.txt"
            p.write_text("\n".join(bg_list) + "\n", encoding="utf-8")
            paths["universe_txt"] = p

            p = out_dir / f"{prefix}_CLASS_SHARED_gene_universe.tsv"
            bg_df.to_csv(p, sep="\t", index=False)

            uni_man = {
                "stage": "9d_export_gene_universe_artifact",
                "timestamp": datetime.now().isoformat(),
                "n_genes": len(bg_list),
                "canonical_id_field": canonical,
                "policy": "explicit_layer3_gene_universe",
                "background_definition": "NON_HK UNION HK",
                "sample_run_id": run_id,
            }
            p = out_dir / f"{prefix}_CLASS_SHARED_gene_universe_manifest.json"
            p.write_text(json.dumps(uni_man, indent=2), encoding="utf-8")

            prov = {
                "sample_run_id": run_id,
                "pipeline": "Layer3_RWR_MultiXrank_batch_monolith",
                "science_version": SCIENCE_VERSION,
                "architecture_version": ARCH_VERSION,
                "notebook": NOTEBOOK_NAME,
                "timestamp": datetime.now().isoformat(),
                "model_id": sample.model_id,
                "lineage": sample.lineage,
                "grn_tsv": str(sample.grn_tsv),
                "grn_sha256": sample.grn_sha256,
                "ppi_sha256": actual_sha256,
                "ppi_path": str(ppi_path),
                "n_grn_raw": n_grn_raw,
                "grn_weight_transform": clip_stats,
                "n_gene_universe": len(gene_universe),
                "n_ppi_edges": len(ppi_f),
                "n_grn_edges": len(grn_f),
                "dl3_07": eq_report,
                "synthetic_class_node": class_node,
                "n_rwr_runs": 1,
                "rwr_runtime_seconds": rwr_sec,
                "lambda": lambda_val,
                "global_score_method": cfg.L3_GLOBAL_GENE_SCORE_METHOD,
                "global_score_reference": cfg.L3_GLOBAL_GENE_SCORE_REFERENCE,
                "housekeeping": hk_meta,
                "hk_filter_qc": hk_qc,
                "rwr_parameters": {
                    "r": cfg.L3_RWR_RESTART_PROB,
                    "delta_gene": cfg.L3_DELTA_GENE,
                    "delta_drug": cfg.L3_DELTA_DRUG,
                    "eta": cfg.L3_ETA_DEFAULT,
                    "tau_gene": cfg.L3_TAU_GENE_DEFAULT,
                    "tau_drug": cfg.L3_TAU_DRUG_DEFAULT,
                },
                "multixrank_version": _MX_VER,
                "top10": global_ranking.head(10)["gene_symbol"].tolist(),
            }
            p = out_dir / f"{prefix}_CLASS_SHARED_provenance.json"
            p.write_text(json.dumps(prov, indent=2, default=str), encoding="utf-8")

            if KEEP_WORKDIR:
                snap = out_dir / "workdir_snapshot"
                if snap.exists():
                    shutil.rmtree(snap)
                shutil.copytree(wdir, snap)

            # keep network files for optional lambda bench (first sample)
            if first_completed_global is None:
                first_completed_global = global_ranking
                first_completed_seeds = seed_genes
                first_completed_node = class_node
                bench_tmpdir = Path(tempfile.mkdtemp(prefix="l3_bench_shared_"))
                shared_ppi_for_bench = bench_tmpdir / "PPI.tsv"
                shared_grn_for_bench = bench_tmpdir / "GRN.tsv"
                shutil.copy2(ppi_file, shared_ppi_for_bench)
                shutil.copy2(grn_file, shared_grn_for_bench)

        finally:
            if not KEEP_WORKDIR:
                shutil.rmtree(wdir, ignore_errors=True)

        elapsed = (datetime.now() - t_sample).total_seconds()
        row.update(
            {
                "status": "COMPLETED",
                "reason": "ok",
                "sample_run_id": run_id,
                "final_top100": str(paths["top100_noHK"]),
                "n_universe": len(bg_list),
                "n_grn_edges_post_clip": clip_stats["n_edges_after"],
                "dl3_07_class_shared": True,
                "synthetic_class_node": class_node,
                "rwr_runtime_sec": round(rwr_sec, 1),
                "elapsed_sec": round(elapsed, 1),
                "output_dir": str(out_dir),
                "top_gene": global_ranking.iloc[0]["gene_symbol"],
            }
        )
        print(
            f"  COMPLETED  top={row['top_gene']}  "
            f"rwr={rwr_sec:.1f}s  universe={len(bg_list):,}"
        )
    except Exception as e:
        row["reason"] = f"{type(e).__name__}: {e}"
        logger.error("[%s] FAILED: %s", sample.model_id, row["reason"])
        print(f"  FAILED: {row['reason']}")
    ledger_rows.append(row)

# %% [markdown]
# ## Stage 10 — Lambda sensitivity (optional, first completed sample)

# %%
lambda_bench_summary: dict[str, Any] = {}
if (
    RUN_LAMBDA_BENCHMARK
    and first_completed_global is not None
    and shared_ppi_for_bench
    and shared_grn_for_bench
    and first_completed_seeds
    and first_completed_node
):
    logger.info("STAGE 10: lambda benchmark (geom-mean)")
    print("\nStage 10: Lambda sensitivity (first completed sample)")

    def jaccard(a: set, b: set) -> float:
        return len(a & b) / len(a | b) if (a | b) else 1.0

    rankings: dict[float, pd.Series] = {}
    for lv in cfg.L3_LAMBDA_BENCH_VALUES:
        try:
            bm = Path(tempfile.mkdtemp()) / f"bm_{lv}"
            net = bm / "network"
            sd = bm / "seeds"
            net.mkdir(parents=True)
            sd.mkdir()
            shutil.copy2(shared_ppi_for_bench, net / "PPI.tsv")
            shutil.copy2(shared_grn_for_bench, net / "GRN.tsv")
            write_drug_monoplex([first_completed_node], net / "drug.tsv")
            write_network_file(
                [(g, first_completed_node) for g in first_completed_seeds],
                net / "bipartite.tsv",
                weighted=False,
            )
            write_seed_file(first_completed_seeds, sd / "seeds.txt")
            conf = build_multixrank_config(
                bm,
                sd / "seeds.txt",
                net / "PPI.tsv",
                net / "GRN.tsv",
                net / "drug.tsv",
                net / "bipartite.tsv",
                float(lv),
            )
            cp = bm / "config.yml"
            with open(cp, "w", encoding="utf-8") as fh:
                yaml.dump(conf, fh, default_flow_style=False)
            res = run_rwr(bm, cp)
            gdf, _ = extract_per_layer_scores(res)
            series = (
                build_global_gene_ranking_from_per_layer(gdf)
                .set_index("gene_symbol")["rwr_score_geom_mean"]
                .sort_values(ascending=False)
            )
            rankings[float(lv)] = series
            shutil.rmtree(bm, ignore_errors=True)
        except Exception as exc:
            logger.warning("lambda bench λ=%s failed: %s", lv, exc)

    rows_bm = []
    lvs = sorted(rankings)
    for i in range(len(lvs)):
        for j in range(i + 1, len(lvs)):
            ri, rj = rankings[lvs[i]], rankings[lvs[j]]
            top_i = set(ri.head(50).index)
            top_j = set(rj.head(50).index)
            common = sorted(set(ri.index) & set(rj.index))
            sp = (
                spearmanr(ri.loc[common], rj.loc[common])[0]
                if len(common) >= 3
                else float("nan")
            )
            rows_bm.append(
                {
                    "lambda_a": lvs[i],
                    "lambda_b": lvs[j],
                    "score_method": cfg.L3_GLOBAL_GENE_SCORE_METHOD,
                    "jaccard_top50": round(jaccard(top_i, top_j), 4),
                    "spearman": round(sp, 4) if not np.isnan(sp) else None,
                }
            )
    bm_df = pd.DataFrame(rows_bm)
    if len(bm_df):
        bm_path = out_root / cfg.L3_LAMBDA_BENCHMARK_TSV
        bm_df.to_csv(bm_path, sep="\t", index=False)
        jac_min = float(bm_df["jaccard_top50"].min())
        sp_min = (
            float(bm_df["spearman"].dropna().min())
            if bm_df["spearman"].notna().any()
            else float("nan")
        )
        lambda_bench_summary = {
            "jaccard_min": jac_min,
            "spearman_min": sp_min if not np.isnan(sp_min) else None,
            "jaccard_stable": jac_min >= cfg.L3_JACCARD_THRESHOLD,
            "spearman_stable": (not np.isnan(sp_min))
            and sp_min >= cfg.L3_SPEARMAN_THRESHOLD,
            "file": str(bm_path),
        }
        print(bm_df.to_string(index=False))
        print(f"  Jaccard min={jac_min:.4f}  Spearman min={sp_min}")

if bench_tmpdir:
    shutil.rmtree(bench_tmpdir, ignore_errors=True)

# %% [markdown]
# ## Stage G — Batch ledger, FAIR catalog, gate report

# %%
logger.info("STAGE G: Ledger + FAIR + gate")

n_ok = sum(1 for r in ledger_rows if r["status"] == "COMPLETED")
n_fail = sum(1 for r in ledger_rows if r["status"] == "FAILED")
elapsed_total = (datetime.now() - batch_t0).total_seconds()

ledger_tsv = out_root / getattr(
    cfg, "L3_BATCH_LEDGER_TSV", "L3_Batch_Ledger.tsv"
)
ledger_json = out_root / getattr(
    cfg, "L3_BATCH_LEDGER_JSON", "L3_Batch_Ledger.json"
)
pd.DataFrame(ledger_rows).to_csv(ledger_tsv, sep="\t", index=False)
ledger_json.write_text(
    json.dumps(
        {
            "batch_id": batch_id,
            "n": len(ledger_rows),
            "entries": ledger_rows,
            "timestamp": datetime.now().isoformat(),
        },
        indent=2,
        default=str,
    ),
    encoding="utf-8",
)

catalog = {
    "@context": "https://www.go-fair.org/fair-principles/",
    "batch_id": batch_id,
    "timestamp": datetime.now().isoformat(),
    "science_version": SCIENCE_VERSION,
    "architecture_version": ARCH_VERSION,
    "notebook": NOTEBOOK_NAME,
    "multixrank_version": _MX_VER,
    "global_score_method": cfg.L3_GLOBAL_GENE_SCORE_METHOD,
    "dl3_07": True,
    "ppi_sha256": actual_sha256,
    "rwr_parameters": {
        "r": cfg.L3_RWR_RESTART_PROB,
        "delta_gene": cfg.L3_DELTA_GENE,
        "delta_drug": cfg.L3_DELTA_DRUG,
        "lambda": lambda_val,
        "eta": cfg.L3_ETA_DEFAULT,
    },
    "samples": [
        {
            "model_id": r["model_id"],
            "lineage": r.get("lineage"),
            "status": r["status"],
            "sample_run_id": r.get("sample_run_id"),
            "top100_noHK": r.get("final_top100"),
        }
        for r in ledger_rows
    ],
    "fair_notes": {
        "F1": "sample_run_id content hash",
        "F2": "provenance.json per sample",
        "R1.2": "L1 P0 + L2 GRN SHA + L3 PPI SHA + RWR params",
        "R1.3": "MultiXrank community workdir + DepMap/STRING IDs",
    },
}
catalog_path = out_root / getattr(
    cfg, "L3_BATCH_CATALOG_JSON", "L3_Batch_Catalog.json"
)
catalog_path.write_text(json.dumps(catalog, indent=2, default=str), encoding="utf-8")

gate = {
    "layer": "Layer3_RWR_MultiXrank_batch",
    "science_version": SCIENCE_VERSION,
    "architecture_version": ARCH_VERSION,
    "notebook": NOTEBOOK_NAME,
    "batch_id": batch_id,
    "timestamp_start": batch_t0.isoformat(),
    "timestamp_end": datetime.now().isoformat(),
    "elapsed_seconds": elapsed_total,
    "multixrank_version": _MX_VER,
    "multixrank_version_ok": _MX_VER == cfg.L3_MULTIXRANK_REQUIRED_VERSION,
    "source": SOURCE,
    "lineage_filter": LINEAGE,
    "n_samples": len(SAMPLES),
    "n_completed": n_ok,
    "n_failed": n_fail,
    "status": "PASS"
    if n_fail == 0 and n_ok > 0
    else ("PARTIAL" if n_ok > 0 else "FAIL"),
    "ppi_sha256": actual_sha256,
    "ppi_sha256_ok": sha256_ok,
    "global_gene_score_method": cfg.L3_GLOBAL_GENE_SCORE_METHOD,
    "dl3_07_single_run_per_sample": True,
    "lambda_benchmark": lambda_bench_summary,
    "deadlock_rules_checked": [
        "dl3_01_p0_no_reweight",
        "dl3_02_grn_clip_applied",
        "dl3_03_frozen_layer2",
        "dl3_04_drug_monoplex",
        "dl3_05_no_aggregation",
        "dl3_06_layer3_ppi_source",
        "dl3_07_single_representative_run",
    ],
    "outputs": {
        "root": str(out_root),
        "ledger_tsv": str(ledger_tsv),
        "ledger_json": str(ledger_json),
        "catalog": str(catalog_path),
    },
    "sample_ids": [s.model_id for s in SAMPLES],
}
gate_path = out_root / getattr(
    cfg, "L3_BATCH_GATE_REPORT_JSON", "L3_Batch_Gate_Report.json"
)
# also write canonical name at LAYER3_OUTPUT_DIR for discoverability
gate_path.write_text(json.dumps(gate, indent=2, default=str), encoding="utf-8")
(cfg.LAYER3_OUTPUT_DIR / "L3_RWR_Gate_Report.json").write_text(
    json.dumps(gate, indent=2, default=str), encoding="utf-8"
)

print(f"\n{'=' * 72}")
print("LAYER 3 BATCH COMPLETE")
print(f"{'=' * 72}")
print(f"  Status:     {gate['status']}")
print(f"  Batch:      {batch_id}")
print(f"  Completed:  {n_ok}/{len(SAMPLES)}")
print(f"  Failed:     {n_fail}")
print(f"  Elapsed:    {elapsed_total:.1f} s")
print(f"  Root:       {out_root}")
print(f"  Ledger:     {ledger_tsv.name}")
print(f"  Catalog:    {catalog_path.name}")
print(f"  Gate:       {gate_path.name}")
print(f"{'=' * 72}")
if n_ok:
    print("\nCompleted samples (first 10):")
    for r in [x for x in ledger_rows if x["status"] == "COMPLETED"][:10]:
        print(
            f"  {r['model_id']}: top={r.get('top_gene')}  "
            f"universe={r.get('n_universe')}"
        )

# expose for interactive follow-up cells
GATE_REPORT = gate
LEDGER_ROWS = ledger_rows
BATCH_OUT_ROOT = out_root
