# Layer 2A: SCENIC GRN Inference Pipelined
# ============================================================
# STAGE 1: IMPORTS & CONFIGURATION
# ============================================================

import pandas as pd
import numpy as np
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

# pySCENIC imports
from arboreto.algo import grnboost2
from ctxcore.rnkdb import FeatherRankingDatabase as RankingDatabase
from pyscenic.utils import modules_from_adjacencies
from pyscenic.prune import prune2df, df2regulons
from pyscenic.aucell import aucell

# Import central config
import config_system as cfg

# Setup logger
logger = cfg.setup_logger("Layer2A")

# Print config summary
cfg.print_config_summary()

# Check SCENIC resources
print("\n" + "="*64)
print("SCENIC Resource Check:")
print("="*64)
resources = cfg.check_scenic_resources()
for key, available in resources.items():
    status = "OK" if available else "MISSING"
    print(f"  {key}: {status}")

# Create output directory
cfg.LAYER2A_OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
print(f"\nOutput directory: {cfg.LAYER2A_OUTPUT_DIR}")
# ============================================================
# STAGE 2: LOAD CCLE DATA & FILTER BY LINEAGE
# ============================================================

logger.info("Loading CCLE Model metadata...")

# Load Model.csv with cell line metadata
model_df = pd.read_csv(cfg.CCLE_MODEL_CSV)
print(f"Model.csv loaded: {len(model_df)} cell lines")
print(f"Columns: {list(model_df.columns)[:10]}...")

# Check lineage distribution
print("\nLineage distribution (OncotreeLineage):")
lineage_counts = model_df[cfg.COL_ONCOTREE_LINEAGE].value_counts()
for lineage in cfg.LINEAGES_OF_INTEREST:
    count = lineage_counts.get(lineage, 0)
    status = "OK" if count >= cfg.MIN_SAMPLES_PER_LINEAGE else f"LOW (<{cfg.MIN_SAMPLES_PER_LINEAGE})"
    print(f"  {lineage}: {count} samples [{status}]")

# Load TPM expression data
logger.info("Loading TPM expression data...")
if not cfg.CCLE_TPM_EXPRESSION_CSV.exists():
    raise FileNotFoundError(
        f"TPM expression file not found: {cfg.CCLE_TPM_EXPRESSION_CSV}\n"
        f"Please download from DepMap 25Q3:\n"
        f"https://depmap.org/portal/data_page/?tab=allData&releasename=DepMap%20Public%2025Q3\n"
        f"File: OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv"
    )

# Đọc file, lấy cột số thứ tự vô danh ngoài cùng làm index (index_col=0)
tpm_df = pd.read_csv(cfg.CCLE_TPM_EXPRESSION_CSV, index_col=0)
print(f"\nRaw TPM loaded: {tpm_df.shape[0]} rows x {tpm_df.shape[1]} columns")

# ============================================================
# STAGE 3: PREPROCESSING
# ============================================================

logger.info("Preprocessing expression data...")

# ---------------------------------------------------------
# [HOTFIX DEPMAP 25Q3]: Chuẩn hóa cấu trúc ma trận
# ---------------------------------------------------------
# 1. Lọc chỉ giữ lại các bản ghi mặc định (tránh trùng lặp)
if 'IsDefaultEntryForModel' in tpm_df.columns:
    n_before = len(tpm_df)
    tpm_df = tpm_df[tpm_df['IsDefaultEntryForModel'] == 'Yes']
    print(f"Filtered Default Entries: {n_before} -> {len(tpm_df)} rows")

# 2. Đặt ModelID làm Index chuẩn để có thể khớp với Model.csv
if 'ModelID' in tpm_df.columns:
    tpm_df = tpm_df.set_index('ModelID')

# 3. Dọn dẹp rác metadata để ma trận CHỈ CÒN chứa dữ liệu Gen
meta_cols = ['SequencingID', 'IsDefaultEntryForModel', 'ModelConditionID', 'IsDefaultEntryForMC']
cols_to_drop = [c for c in meta_cols if c in tpm_df.columns]
if cols_to_drop:
    tpm_df = tpm_df.drop(columns=cols_to_drop)

print(f"Processed TPM Matrix (Genes only): {tpm_df.shape[0]} samples x {tpm_df.shape[1]} genes")
# ---------------------------------------------------------

# CRITICAL: Strip Entrez IDs from gene names (DL5 compliance)
# Column format: "GENE (ENTREZ_ID)" -> "GENE"
print("\nBefore stripping Entrez IDs:")
print(f"  Example columns: {list(tpm_df.columns)[:3]}")

tpm_df.columns = [c.split(' (')[0] for c in tpm_df.columns]

print("After stripping Entrez IDs:")
print(f"  Example columns: {list(tpm_df.columns)[:3]}")

# Validate DL5 compliance
cfg.validate_deadlock_rules("gene_name_format", stripped_entrez=True)
logger.info("DL5 compliance: Gene names stripped of Entrez IDs")

# Handle NaN values
nan_count_before = tpm_df.isna().sum().sum()
tpm_df = tpm_df.fillna(0)
print(f"\nNaN handling: {nan_count_before} NaN values replaced with 0")

# Load TF list
logger.info("Loading TF list...")
tf_list = pd.read_csv(cfg.SCENIC_TF_LIST, header=None)[0].tolist()
print(f"TF list loaded: {len(tf_list)} transcription factors")

# Intersect TFs with expression data
tf_in_data = [tf for tf in tf_list if tf in tpm_df.columns]
print(f"TFs present in expression data: {len(tf_in_data)} / {len(tf_list)}")

# Verify matrix orientation (DL4 compliance)
# GRNBoost2 expects: Rows=Samples (ACH-xxxx), Cols=Genes (TSPAN6, TNMD...)
cfg.validate_deadlock_rules("grnboost2_orientation", transposed=False)
logger.info("DL4 compliance: Matrix orientation correct (Rows=Samples, Cols=Genes)")
# ============================================================
# STAGE 4: GRNBOOST2 INFERENCE (PER LINEAGE)
# ============================================================

logger.info("Starting GRNBoost2 inference...")

# Store results
all_adjacencies = {}
all_expr_subsets = {}  # HOTFIX: Store lineage-specific expression for Stage 5 (avoid data leakage)
lineage_stats = []

for lineage in cfg.LINEAGES_OF_INTEREST:
    print(f"\n{'='*64}")
    print(f"Processing lineage: {lineage}")
    print(f"{'='*64}")
    
    # Get samples for this lineage
    lineage_samples = model_df[
        model_df[cfg.COL_ONCOTREE_LINEAGE] == lineage
    ][cfg.COL_MODEL_ID].tolist()
    
    # Intersect with expression data index
    samples_in_expr = [s for s in lineage_samples if s in tpm_df.index]
    
    print(f"  Samples in Model.csv: {len(lineage_samples)}")
    print(f"  Samples in expression data: {len(samples_in_expr)}")
    
    if len(samples_in_expr) < cfg.MIN_SAMPLES_PER_LINEAGE:
        logger.warning(f"Skipping {lineage}: only {len(samples_in_expr)} samples (min={cfg.MIN_SAMPLES_PER_LINEAGE})")
        lineage_stats.append({
            "lineage": lineage,
            "n_samples": len(samples_in_expr),
            "status": "SKIPPED",
            "n_edges": 0
        })
        continue
    
    # Subset expression matrix
    expr_subset = tpm_df.loc[samples_in_expr]
    print(f"  Expression matrix: {expr_subset.shape[0]} samples x {expr_subset.shape[1]} genes")
    
    # Run GRNBoost2
    logger.info(f"Running GRNBoost2 for {lineage}...")
    print(f"  GRNBoost2 parameters:")
    print(f"    n_workers: {cfg.GRNBOOST2_N_WORKERS}")
    print(f"    seed: {cfg.GRNBOOST2_SEED}")
    
    adjacencies = grnboost2(
        expression_data=expr_subset,
        tf_names=tf_in_data,
        seed=cfg.GRNBOOST2_SEED,
        client_or_address='custom_multiprocessing',
        num_workers=cfg.GRNBOOST2_N_WORKERS,
        verbose=cfg.GRNBOOST2_VERBOSE
    )
    
    # IMPORTANT: Keep original column names (TF, target, importance) for modules_from_adjacencies()
    # pySCENIC's modules_from_adjacencies() expects hardcoded column names: TF, target, importance
    # Column renaming to Source/Target/Weight deferred to Stage 6 export for FAIR compliance
    
    # Add lineage column (modules_from_adjacencies ignores extra columns)
    adjacencies['Lineage'] = lineage
    
    print(f"  GRNBoost2 output: {len(adjacencies)} edges")
    print(f"  Top 5 edges (original column names):")
    print(adjacencies.head())
    
    # Save lineage-specific adjacencies with FAIR-compliant column names
    adj_export = adjacencies.rename(columns={
        cfg.COL_GRNBOOST2_TF: cfg.COL_GRN_SOURCE,
        cfg.COL_GRNBOOST2_TARGET: cfg.COL_GRN_TARGET,
        cfg.COL_GRNBOOST2_IMPORTANCE: cfg.COL_GRN_WEIGHT
    })
    output_path = cfg.LAYER2A_OUTPUT_DIR / cfg.L2A_ADJACENCIES_CSV.format(lineage=lineage)
    adj_export.to_csv(output_path, index=False)
    logger.info(f"Saved adjacencies to: {output_path}")
    
    # Store for aggregation
    all_adjacencies[lineage] = adjacencies
    all_expr_subsets[lineage] = expr_subset  # HOTFIX: Store for cisTarget in Stage 5
    lineage_stats.append({
        "lineage": lineage,
        "n_samples": len(samples_in_expr),
        "status": "COMPLETED",
        "n_edges": len(adjacencies)
    })

print(f"\n{'='*64}")
print("GRNBoost2 Summary:")
print(f"{'='*64}")
stats_df = pd.DataFrame(lineage_stats)
print(stats_df.to_string(index=False))
# ============================================================
# STAGE 5: CISTARGET PRUNING (OPTIONAL - GRACEFUL SKIP)
# ============================================================

logger.info("Checking cisTarget resources...")

# Check if motif databases are available
motif_db_available = (
    cfg.SCENIC_MOTIF_DB_DIR.exists() and 
    str(cfg.SCENIC_MOTIF_DB_DIR) != "<USER_FILL_IN_PATH>"
)
motif_annot_available = (
    cfg.SCENIC_MOTIF_ANNOTATIONS.exists() and 
    str(cfg.SCENIC_MOTIF_ANNOTATIONS) != "<USER_FILL_IN_PATH>"
)

if motif_db_available and motif_annot_available:
    logger.info("cisTarget resources available - running motif pruning...")
    
    # Load ranking databases
    db_fnames = list(cfg.SCENIC_MOTIF_DB_DIR.glob("*.feather"))
    print(f"Found {len(db_fnames)} ranking databases")
    
    dbs = [RankingDatabase(fname=str(f), name=f.stem) for f in db_fnames]
    
    # Process each lineage
    all_regulons = {}
    for lineage, adjacencies in all_adjacencies.items():
        print(f"\nPruning {lineage} adjacencies...")
        
        # Create modules from adjacencies
        # HOTFIX: Use lineage-specific expression (not global tpm_df) to avoid data leakage
        modules = modules_from_adjacencies(adjacencies, all_expr_subsets[lineage])
        
        # Prune using cisTarget
        df_regulons = prune2df(
            dbs,
            modules,
            str(cfg.SCENIC_MOTIF_ANNOTATIONS),
            num_workers=cfg.GRNBOOST2_N_WORKERS,
            nes_threshold=cfg.CISTARGET_NES_THRESHOLD,
            rank_threshold=cfg.CISTARGET_RANK_THRESHOLD,
            auc_threshold=cfg.CISTARGET_AUC_THRESHOLD,
            motif_similarity_fdr=cfg.CISTARGET_MOTIF_SIMILARITY_FDR
        )
        
        # Convert to regulons
        regulons = df2regulons(df_regulons)
        print(f"  {len(regulons)} regulons identified")
        
        # Save regulons
        regulon_output = cfg.LAYER2A_OUTPUT_DIR / cfg.L2A_REGULONS_CSV.format(lineage=lineage)
        df_regulons.to_csv(regulon_output, index=False)
        
        all_regulons[lineage] = regulons
        
else:
    logger.warning(
        "cisTarget resources not available - SKIPPING motif pruning.\n"
        "Raw GRNBoost2 adjacencies will be used as preliminary GRN.\n"
        "\n"
        "To enable cisTarget pruning, download from Aerts Lab:\n"
        "  - Motif databases: https://resources.aertslab.org/cistarget/databases/\n"
        "  - Motif annotations: https://resources.aertslab.org/cistarget/motif2tf/\n"
        "\n"
        "Then update config_system.py:\n"
        "  SCENIC_MOTIF_DB_DIR = Path(r'path/to/databases')\n"
        "  SCENIC_MOTIF_ANNOTATIONS = Path(r'path/to/motifs.tbl')"
    )
    print("\nProceeding with raw GRNBoost2 adjacencies...")
# ============================================================
# STAGE 6: EXPORT MASTER FILES (ADJACENCIES + REGULONS)
# ============================================================
# HOTFIX: Separate raw GRNBoost2 adjacencies from cisTarget-pruned regulons
# - Master Adjacencies: Raw co-expression edges (always exported)
# - Master Regulons: Motif-pruned validated edges with NES scores (if cisTarget ran)
#
# SMART Compliance:
#   S - Export two distinct output files with clear naming
#   M - Include NES scores for confidence ranking
#   A - Handle graceful skip when cisTarget unavailable
#   R - Correct terminology (adjacencies vs regulons)
#   T - Immediate export after inference
#
# FAIR Compliance:
#   F - Persistent filenames defined in config_system.py
#   A - Human-readable CSV format
#   I - Standard pySCENIC output columns
#   R - Full provenance (lineage, NES, context)
# ============================================================

logger.info("Aggregating all lineage GRNs into master files...")

# --- Part A: Export Master Adjacencies (raw GRNBoost2 output) ---
if all_adjacencies:
    master_adj_df = pd.concat(all_adjacencies.values(), ignore_index=True)
    # Sort by original column name (importance) before renaming
    master_adj_df = master_adj_df.sort_values(cfg.COL_GRNBOOST2_IMPORTANCE, ascending=False)
    
    # Rename columns for FAIR-compliant export: TF->Source, target->Target, importance->Weight
    master_adj_df = master_adj_df.rename(columns={
        cfg.COL_GRNBOOST2_TF: cfg.COL_GRN_SOURCE,
        cfg.COL_GRNBOOST2_TARGET: cfg.COL_GRN_TARGET,
        cfg.COL_GRNBOOST2_IMPORTANCE: cfg.COL_GRN_WEIGHT
    })
    
    master_adj_output = cfg.LAYER2A_OUTPUT_DIR / cfg.L2A_MASTER_ADJACENCIES_CSV
    master_adj_df.to_csv(master_adj_output, index=False)
    
    print(f"\n{'='*64}")
    print("MASTER ADJACENCIES (raw GRNBoost2)")
    print(f"{'='*64}")
    print(f"Saved to: {master_adj_output}")
    print(f"Total edges: {len(master_adj_df)}")
    print(f"Unique TFs: {master_adj_df[cfg.COL_GRN_SOURCE].nunique()}")
    print(f"Unique targets: {master_adj_df[cfg.COL_GRN_TARGET].nunique()}")
    print(f"Lineages: {master_adj_df['Lineage'].unique().tolist()}")
    
    print("\nTop 10 edges by weight:")
    print(master_adj_df.head(10).to_string(index=False))
else:
    logger.error("No adjacencies to export - all lineages skipped.")

# --- Part B: Export Master Regulons (cisTarget pruned, if available) ---
# Check if all_regulons exists and has data (set in Stage 5)
if 'all_regulons' in dir() and all_regulons:
    print(f"\n{'='*64}")
    print("MASTER REGULONS (cisTarget pruned)")
    print(f"{'='*64}")
    
    # Build regulon dataframe with NES scores
    regulon_records = []
    for lineage, regulons in all_regulons.items():
        for regulon in regulons:
            tf = regulon.transcription_factor
            # Get context (activator/repressor) - frozenset, convert to string
            context_str = ','.join(sorted(regulon.context)) if regulon.context else 'unknown'
            # NES score may be stored as 'score' attribute
            nes_score = getattr(regulon, 'score', None)
            
            for target, weight in regulon.gene2weight.items():
                regulon_records.append({
                    cfg.COL_GRN_SOURCE: tf,
                    cfg.COL_GRN_TARGET: target,
                    cfg.COL_GRN_WEIGHT: weight,
                    'Context': context_str,
                    'NES': nes_score,
                    'Lineage': lineage
                })
    
    if regulon_records:
        master_reg_df = pd.DataFrame(regulon_records)
        master_reg_df = master_reg_df.sort_values(cfg.COL_GRN_WEIGHT, ascending=False)
        
        master_reg_output = cfg.LAYER2A_OUTPUT_DIR / cfg.L2A_MASTER_REGULONS_CSV
        master_reg_df.to_csv(master_reg_output, index=False)
        
        print(f"Saved to: {master_reg_output}")
        print(f"Total edges: {len(master_reg_df)}")
        print(f"Unique TFs: {master_reg_df[cfg.COL_GRN_SOURCE].nunique()}")
        print(f"Unique targets: {master_reg_df[cfg.COL_GRN_TARGET].nunique()}")
        
        print("\nTop 10 regulon edges:")
        print(master_reg_df.head(10).to_string(index=False))
    else:
        logger.warning("Regulons exist but no edges extracted.")
else:
    print(f"\n{'='*64}")
    print("MASTER REGULONS: SKIPPED")
    print(f"{'='*64}")
    print("cisTarget was not run (motif databases unavailable).")
    print("Only raw adjacencies have been exported.")
    print("Re-run with motif databases to generate pruned regulons.")

# --- Summary ---
print(f"\n{'='*64}")
print("Layer 2A SCENIC GRN Inference - COMPLETE")
print(f"{'='*64}")
print(f"Output directory: {cfg.LAYER2A_OUTPUT_DIR}")
print(f"\nGenerated files:")
for f in cfg.LAYER2A_OUTPUT_DIR.glob("*.csv"):
    size_kb = f.stat().st_size / 1024
    print(f"  {f.name}: {size_kb:.1f} KB")
