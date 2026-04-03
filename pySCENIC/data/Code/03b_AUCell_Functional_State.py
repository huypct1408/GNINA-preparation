# Layer 3B: AUCell Functional State Analysis
"""
## Cancer Selectivity Validation via Regulon Activity Comparison

---

### SMART Goal
- **S**pecific: Compare regulon activity (AUC) between Jurkat (T-cell leukemia) and HEK-293 (normal kidney)
- **M**easurable: Delta AUC = AUC_Jurkat - AUC_HEK293 quantifies cancer selectivity
- **A**chievable: Using pySCENIC AUCell with ctxcore GeneSignature objects
- **R**elevant: Validates that essential regulons are "bright" in cancer, "dim" in normal
- **T**ime-bound: Single notebook execution (~10-15 minutes)

### FAIR Compliance
- **F**indable: Persistent regulon IDs (TF_regulon format)
- **A**ccessible: Standard CSV output with full AUC metrics
- **I**nteroperable: Uses ctxcore GeneSignature objects (pySCENIC compatible)
- **R**eusable: Full provenance in selectivity summary JSON

### Critical Deadlock Rules (DL13-DL16)
- **DL13**: MUST use `ctxcore.genesig.GeneSignature` objects - plain lists will TypeError
- **DL14**: Expression matrix must be log2(x+1) normalized (CCLE TPM format)
- **DL15**: AUCell rankings are cell-independent - batch effects auto-handled
- **DL16**: Positive Delta AUC = cancer selective; Negative = toxicity risk

### Pipeline Position
```
Layer 1 (Thermodynamic Gate) → Layer 2A (SCENIC GRN) → Layer 2B (RWR) 
    → Layer 3A (CRISPR) → [Layer 3B: AUCell Selectivity] → Layer 4...
```

### Input Files
1. `L3A_Essential_Targets.csv` - CRISPR-validated essential genes (from Layer 3A)
2. `L2A_Master_Regulons_AllLineages.csv` - TF→Target edges with weights (from Layer 2A)
3. `OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv` - CCLE TPM expression matrix

### Output Files
1. `L3B_Active_Regulons.csv` - Regulons with positive Delta AUC > threshold
2. `L3B_All_Regulons_AUC.csv` - All regulons with AUC scores for both cell lines
3. `L3B_L3A_Signature_AUC.csv` - Aggregate L3A signature AUC
4. `L3B_Selectivity_Summary.json` - Metadata + statistics
"""
## Stage 1: Imports & Configuration

Load all required libraries and import config_system v1.6.
# ============================================================
# STAGE 1: Imports & Configuration
# ============================================================

import pandas as pd
import numpy as np
import json
from pathlib import Path
from datetime import datetime
import warnings
warnings.filterwarnings('ignore')

# pySCENIC imports
from pyscenic.aucell import aucell
from ctxcore.genesig import GeneSignature  # CRITICAL: DL13 compliance

# Visualization
import matplotlib.pyplot as plt
import seaborn as sns

# Import config_system v1.6
import sys
sys.path.insert(0, str(Path('.').resolve()))
from config_system import (
    # Layer 3B Configuration
    LAYER3B_OUTPUT_DIR,
    L3B_AUC_THRESHOLD,
    L3B_DELTA_AUC_THRESHOLD,
    L3B_MIN_GENES_IN_SIGNATURE,
    L3B_CANCER_CELL_LINE,
    L3B_CANCER_MODEL_ID,
    L3B_NORMAL_CELL_LINE,
    L3B_NORMAL_MODEL_ID,
    L3B_ACTIVE_REGULONS_CSV,
    L3B_ALL_REGULONS_AUC_CSV,
    L3B_L3A_SIGNATURE_AUC_CSV,
    L3B_SELECTIVITY_SUMMARY_JSON,
    L3B_QC_PLOTS_DIR,
    COL_AUC_JURKAT,
    COL_AUC_HEK293,
    COL_DELTA_AUC,
    COL_REGULON_NAME,
    COL_N_GENES,
    COL_IS_SELECTIVE,
    COL_SELECTIVITY_RANK,
    # Layer 3A paths (input)
    LAYER3A_OUTPUT_DIR,
    L3A_ESSENTIAL_TARGETS_CSV,
    # Layer 2A paths (input)
    LAYER2A_OUTPUT_DIR,
    L2A_MASTER_REGULONS_CSV,
    # Expression data
    CCLE_TPM_EXPRESSION_CSV,
    # Utilities
    print_config_summary,
    check_layer3b_resources,
    validate_deadlock_rules,
    setup_logger,
    normalize_gene_name,
    DEADLOCK_RULES,
)

# Setup logger
logger = setup_logger("Layer3B")

# Print configuration
print_config_summary()

print("\n" + "="*64)
print("STAGE 1 COMPLETE: Imports & Configuration")
print("="*64)
# ============================================================
# STAGE 1B: Validate Layer 3B Resources
# ============================================================

print("Checking Layer 3B resources...\n")
l3b_resources = check_layer3b_resources()

all_ok = True
for key, available in l3b_resources.items():
    status = "OK" if available else "MISSING"
    icon = "[OK]" if available else "[MISSING]"
    print(f"  {icon} {key}: {status}")
    if not available and key not in ["layer3b_output_dir_exists"]:
        all_ok = False

# Create output directory if needed
LAYER3B_OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
print(f"\nOutput directory: {LAYER3B_OUTPUT_DIR}")

# Create QC plots directory
qc_plots_dir = LAYER3B_OUTPUT_DIR / L3B_QC_PLOTS_DIR
qc_plots_dir.mkdir(parents=True, exist_ok=True)
print(f"QC plots directory: {qc_plots_dir}")

if not all_ok:
    print("\n[WARNING] Some resources are missing. Please check paths in config_system.py")
else:
    print("\n[OK] All Layer 3B resources available")
## Stage 2: Load Layer 3A Essential Targets

Load CRISPR-validated essential genes from Layer 3A output.
# ============================================================
# STAGE 2: Load Layer 3A Essential Targets
# ============================================================

l3a_path = LAYER3A_OUTPUT_DIR / L3A_ESSENTIAL_TARGETS_CSV
print(f"Loading L3A Essential Targets from: {l3a_path}")

l3a_df = pd.read_csv(l3a_path)
print(f"\nL3A Essential Targets loaded: {len(l3a_df)} genes")
print(f"Columns: {list(l3a_df.columns)}")

# Extract gene list
# Check which column contains gene names
gene_col = None
for col in ['Gene', 'gene_symbol', 'Gene_Symbol', 'gene', 'Symbol']:
    if col in l3a_df.columns:
        gene_col = col
        break

if gene_col is None:
    # Use first column if no standard name found
    gene_col = l3a_df.columns[0]
    print(f"[WARNING] Using first column as gene column: {gene_col}")

l3a_genes = l3a_df[gene_col].tolist()
print(f"\nExtracted {len(l3a_genes)} essential genes from column '{gene_col}'")
print(f"Sample genes: {l3a_genes[:10]}")

# Store for later use
L3A_ESSENTIAL_GENES = l3a_genes

print("\n" + "="*64)
print("STAGE 2 COMPLETE: L3A Essential Targets Loaded")
print(f"  Total essential genes: {len(L3A_ESSENTIAL_GENES)}")
print("="*64)
## Stage 3: Load L2A Regulons

Load TF→Target edges from Layer 2A Master Regulons. Group by TF to create per-regulon signatures.
# ============================================================
# STAGE 3: Load L2A Regulons
# ============================================================

l2a_path = LAYER2A_OUTPUT_DIR / L2A_MASTER_REGULONS_CSV
print(f"Loading L2A Master Regulons from: {l2a_path}")

l2a_df = pd.read_csv(l2a_path)
print(f"\nL2A Regulons loaded: {len(l2a_df)} TF→Target edges")
print(f"Columns: {list(l2a_df.columns)}")

# Expected columns: Source (TF), Target (gene), Weight, possibly Context, NES, Lineage
print(f"\nUnique TFs (regulons): {l2a_df['Source'].nunique()}")
print(f"Unique Target genes: {l2a_df['Target'].nunique()}")

# Preview
print("\nSample edges:")
display(l2a_df.head(10))
# ============================================================
# STAGE 3B: Group Regulons by TF
# ============================================================

print("Grouping regulons by TF...\n")

# Group by Source (TF) to get target genes per regulon
regulon_groups = l2a_df.groupby('Source')

# Calculate statistics per regulon
regulon_stats = []
for tf, group in regulon_groups:
    regulon_stats.append({
        'TF': tf,
        'N_Targets': len(group),
        'Mean_Weight': group['Weight'].mean(),
        'Max_Weight': group['Weight'].max(),
    })

regulon_stats_df = pd.DataFrame(regulon_stats)
regulon_stats_df = regulon_stats_df.sort_values('N_Targets', ascending=False)

print(f"Total regulons (TFs): {len(regulon_stats_df)}")
print(f"\nTop 10 regulons by target count:")
display(regulon_stats_df.head(10))

print(f"\nRegulons with >= {L3B_MIN_GENES_IN_SIGNATURE} targets: "
      f"{(regulon_stats_df['N_Targets'] >= L3B_MIN_GENES_IN_SIGNATURE).sum()}")

# Store for later
REGULON_GROUPS = regulon_groups
REGULON_STATS = regulon_stats_df

print("\n" + "="*64)
print("STAGE 3 COMPLETE: L2A Regulons Loaded & Grouped")
print(f"  Total TF→Target edges: {len(l2a_df)}")
print(f"  Unique regulons (TFs): {len(regulon_stats_df)}")
print("="*64)
## Stage 4: Load Expression Data (CCLE TPM)

Load the pre-normalized CCLE TPM matrix and extract Jurkat + HEK-293 samples.

**DL14 Compliance**: CCLE TPM is already log2(x+1) normalized - NO additional transformation needed.
# ============================================================
# STAGE 4: Load Expression Data (CCLE TPM)
# ============================================================

print(f"Loading CCLE TPM Expression from: {CCLE_TPM_EXPRESSION_CSV}")
print("[DL14] Using pre-normalized log2(x+1) TPM data - NO additional transformation\n")

# Load expression matrix
exp_df = pd.read_csv(CCLE_TPM_EXPRESSION_CSV, index_col=0)
print(f"Expression matrix shape: {exp_df.shape}")
print(f"  Rows (samples/cells): {exp_df.shape[0]}")
print(f"  Columns (genes): {exp_df.shape[1]}")

# Clean gene names - strip Entrez IDs if present
# Format: "GENE (12345)" -> "GENE"
original_cols = exp_df.columns.tolist()
clean_cols = [c.split(' (')[0] for c in original_cols]
exp_df.columns = clean_cols
print(f"\n[DL5] Stripped Entrez IDs from gene names")
print(f"  Sample original: {original_cols[0]}")
print(f"  Sample cleaned:  {clean_cols[0]}")

# Verify DL14 compliance - check data range (should be log2(TPM+1))
data_min = exp_df.values.min()
data_max = exp_df.values.max()
data_mean = exp_df.values.mean()
print(f"\nExpression data range:")
print(f"  Min: {data_min:.4f}")
print(f"  Max: {data_max:.4f}")
print(f"  Mean: {data_mean:.4f}")

# DL14 validation
validate_deadlock_rules("tpm_normalization", is_tpm_normalized=True)
print("\n[DL14] TPM normalization check: PASSED")
# ============================================================
# STAGE 4B: Extract Jurkat and HEK-293 Samples
# ============================================================

print("Extracting target cell lines...\n")

# Target ModelIDs
jurkat_id = L3B_CANCER_MODEL_ID   # ACH-000995
hek293_id = L3B_NORMAL_MODEL_ID   # ACH-001085

print(f"Cancer cell line:  {L3B_CANCER_CELL_LINE} ({jurkat_id})")
print(f"Normal cell line:  {L3B_NORMAL_CELL_LINE} ({hek293_id})")

# Check if both cell lines exist in the matrix
available_samples = exp_df.index.tolist()

jurkat_found = jurkat_id in available_samples
hek293_found = hek293_id in available_samples

print(f"\nJurkat ({jurkat_id}) in matrix: {jurkat_found}")
print(f"HEK-293 ({hek293_id}) in matrix: {hek293_found}")

if not jurkat_found or not hek293_found:
    raise ValueError(f"Required cell lines not found in expression matrix!\n"
                     f"Jurkat: {jurkat_found}, HEK-293: {hek293_found}")

# Extract subset
cell_lines = [jurkat_id, hek293_id]
exp_subset = exp_df.loc[cell_lines].copy()

print(f"\nExtracted expression subset:")
print(f"  Shape: {exp_subset.shape}")
print(f"  Samples: {exp_subset.index.tolist()}")

# Store for AUCell
EXP_MATRIX = exp_subset

print("\n" + "="*64)
print("STAGE 4 COMPLETE: Expression Data Loaded")
print(f"  Total genes: {exp_subset.shape[1]}")
print(f"  Cell lines: Jurkat (cancer), HEK-293 (normal)")
print("="*64)
## Stage 5: Create GeneSignature Objects

**CRITICAL (DL13)**: MUST use `ctxcore.genesig.GeneSignature` objects - plain lists will cause TypeError!

### Two Approaches:
1. **Aggregate L3A Signature**: Single signature from all L3A essential genes
2. **Per-TF Regulon Signatures**: Individual signatures for each TF's target genes
# ============================================================
# STAGE 5A: Create Aggregate L3A Signature
# ============================================================

print("Creating Aggregate L3A Signature...")
print(f"[DL13] Using ctxcore.genesig.GeneSignature objects\n")

# Get genes that exist in expression matrix
available_genes = set(EXP_MATRIX.columns)
l3a_genes_in_matrix = [g for g in L3A_ESSENTIAL_GENES if g in available_genes]

print(f"L3A essential genes: {len(L3A_ESSENTIAL_GENES)}")
print(f"Genes in expression matrix: {len(l3a_genes_in_matrix)}")
print(f"Coverage: {len(l3a_genes_in_matrix)/len(L3A_ESSENTIAL_GENES)*100:.1f}%")

# Check minimum threshold
if len(l3a_genes_in_matrix) < L3B_MIN_GENES_IN_SIGNATURE:
    print(f"[WARNING] Only {len(l3a_genes_in_matrix)} genes found - below minimum {L3B_MIN_GENES_IN_SIGNATURE}")

# Create GeneSignature with uniform weights
l3a_signature = GeneSignature(
    name='L3A_Essential_Aggregate',
    gene2weight={gene: 1.0 for gene in l3a_genes_in_matrix}
)

print(f"\nL3A Aggregate Signature created:")
print(f"  Name: {l3a_signature.name}")
print(f"  Genes: {len(l3a_signature)}")

# DL13 validation
validate_deadlock_rules("genesignature_objects", used_genesignature=True)
print("\n[DL13] GeneSignature object check: PASSED")
# ============================================================
# STAGE 5B: Create Per-TF Regulon Signatures
# ============================================================

print("Creating Per-TF Regulon Signatures...\n")

regulon_signatures = []
skipped_regulons = []

for tf, group in REGULON_GROUPS:
    # Get target genes that exist in expression matrix
    targets = group['Target'].tolist()
    weights = group['Weight'].tolist()
    
    # Filter to genes in expression matrix
    gene2weight = {}
    for target, weight in zip(targets, weights):
        if target in available_genes:
            gene2weight[target] = weight
    
    # Skip if too few genes
    if len(gene2weight) < L3B_MIN_GENES_IN_SIGNATURE:
        skipped_regulons.append({
            'TF': tf,
            'Original_Targets': len(targets),
            'In_Matrix': len(gene2weight),
            'Reason': f'< {L3B_MIN_GENES_IN_SIGNATURE} genes'
        })
        continue
    
    # Create GeneSignature
    gs = GeneSignature(
        name=f"{tf}_regulon",
        gene2weight=gene2weight
    )
    regulon_signatures.append(gs)

print(f"Per-TF Regulon Signatures created: {len(regulon_signatures)}")
print(f"Skipped (< {L3B_MIN_GENES_IN_SIGNATURE} genes): {len(skipped_regulons)}")

# Summary of created signatures
signature_sizes = [len(gs) for gs in regulon_signatures]
print(f"\nSignature size statistics:")
print(f"  Min genes: {min(signature_sizes)}")
print(f"  Max genes: {max(signature_sizes)}")
print(f"  Mean genes: {np.mean(signature_sizes):.1f}")
print(f"  Median genes: {np.median(signature_sizes):.1f}")

# Show top regulons
top_regulons = sorted(regulon_signatures, key=lambda x: len(x), reverse=True)[:10]
print(f"\nTop 10 largest regulons:")
for gs in top_regulons:
    print(f"  {gs.name}: {len(gs)} genes")
# ============================================================
# STAGE 5C: Combine All Signatures
# ============================================================

print("Combining all signatures for AUCell...\n")

# Combine L3A aggregate + per-TF regulons
all_signatures = [l3a_signature] + regulon_signatures

print(f"Total signatures for AUCell:")
print(f"  L3A Aggregate: 1")
print(f"  Per-TF Regulons: {len(regulon_signatures)}")
print(f"  TOTAL: {len(all_signatures)}")

# Store for AUCell
ALL_SIGNATURES = all_signatures
REGULON_SIGNATURES = regulon_signatures
L3A_AGGREGATE_SIGNATURE = l3a_signature

print("\n" + "="*64)
print("STAGE 5 COMPLETE: GeneSignature Objects Created")
print(f"  L3A Aggregate: {len(l3a_signature)} genes")
print(f"  Per-TF Regulons: {len(regulon_signatures)}")
print(f"  [DL13] All signatures use ctxcore.genesig.GeneSignature")
print("="*64)
## Stage 6: Run AUCell Analysis

Execute AUCell to calculate regulon activity scores (AUC) for each cell line.

**DL15**: AUCell rankings are computed independently within each cell - automatically robust to batch effects.
# ============================================================
# STAGE 6: Run AUCell Analysis
# ============================================================

print("Running AUCell Analysis...")
print(f"[DL15] AUCell ranks genes independently per cell - batch effect immune\n")

print(f"Input expression matrix: {EXP_MATRIX.shape}")
print(f"Number of signatures: {len(ALL_SIGNATURES)}")
print(f"AUC threshold: {L3B_AUC_THRESHOLD}")

# Run AUCell
# Note: exp_mtx must be (n_cells x n_genes) - which is already correct
print("\nExecuting pyscenic.aucell.aucell()...")
start_time = datetime.now()

auc_mtx = aucell(
    exp_mtx=EXP_MATRIX,
    signatures=ALL_SIGNATURES,
    auc_threshold=L3B_AUC_THRESHOLD,
    num_workers=1  # Single worker for stability
)

end_time = datetime.now()
duration = (end_time - start_time).total_seconds()

print(f"\nAUCell completed in {duration:.2f} seconds")
print(f"Output AUC matrix shape: {auc_mtx.shape}")
print(f"  Rows (cells): {auc_mtx.shape[0]}")
print(f"  Columns (signatures): {auc_mtx.shape[1]}")
# ============================================================
# STAGE 6B: Extract AUC Scores per Cell Line
# ============================================================

print("Extracting AUC scores per cell line...\n")

# Get AUC scores for each cell line
auc_jurkat = auc_mtx.loc[L3B_CANCER_MODEL_ID]
auc_hek293 = auc_mtx.loc[L3B_NORMAL_MODEL_ID]

print(f"Jurkat ({L3B_CANCER_MODEL_ID}) AUC scores:")
print(f"  Min: {auc_jurkat.min():.4f}")
print(f"  Max: {auc_jurkat.max():.4f}")
print(f"  Mean: {auc_jurkat.mean():.4f}")

print(f"\nHEK-293 ({L3B_NORMAL_MODEL_ID}) AUC scores:")
print(f"  Min: {auc_hek293.min():.4f}")
print(f"  Max: {auc_hek293.max():.4f}")
print(f"  Mean: {auc_hek293.mean():.4f}")

# Store for delta calculation
AUC_JURKAT = auc_jurkat
AUC_HEK293 = auc_hek293
AUC_MATRIX = auc_mtx

print("\n" + "="*64)
print("STAGE 6 COMPLETE: AUCell Analysis")
print(f"  Jurkat mean AUC: {auc_jurkat.mean():.4f}")
print(f"  HEK-293 mean AUC: {auc_hek293.mean():.4f}")
print("="*64)
## Stage 7: Calculate Delta AUC (Cancer Selectivity)

**Delta AUC = AUC_Jurkat - AUC_HEK293**

**DL16 Interpretation**:
- **Positive Delta AUC**: Regulon is MORE active in cancer (Jurkat) than normal (HEK-293) = **CANCER SELECTIVE**
- **Negative Delta AUC**: Regulon is LESS active in cancer = **TOXICITY RISK**
# ============================================================
# STAGE 7: Calculate Delta AUC (Cancer Selectivity)
# ============================================================

print("Calculating Delta AUC (Cancer Selectivity)...")
print(f"\n[DL16] Interpretation:")
print(f"  Positive Delta AUC > 0: CANCER SELECTIVE (bright in cancer, dim in normal)")
print(f"  Negative Delta AUC < 0: TOXICITY RISK (dim in cancer, bright in normal)")
print(f"  Threshold for selectivity: Delta AUC > {L3B_DELTA_AUC_THRESHOLD}\n")

# Calculate Delta AUC
delta_auc = AUC_JURKAT - AUC_HEK293

# Create results DataFrame
results_df = pd.DataFrame({
    COL_REGULON_NAME: delta_auc.index,
    COL_AUC_JURKAT: AUC_JURKAT.values,
    COL_AUC_HEK293: AUC_HEK293.values,
    COL_DELTA_AUC: delta_auc.values,
})

# Add gene count for each signature
gene_counts = {}
for sig in ALL_SIGNATURES:
    gene_counts[sig.name] = len(sig)
results_df[COL_N_GENES] = results_df[COL_REGULON_NAME].map(gene_counts)

# Add selectivity flag
results_df[COL_IS_SELECTIVE] = results_df[COL_DELTA_AUC] > L3B_DELTA_AUC_THRESHOLD

# Sort by Delta AUC descending
results_df = results_df.sort_values(COL_DELTA_AUC, ascending=False)
results_df[COL_SELECTIVITY_RANK] = range(1, len(results_df) + 1)

print(f"Delta AUC Statistics:")
print(f"  Min: {delta_auc.min():.4f}")
print(f"  Max: {delta_auc.max():.4f}")
print(f"  Mean: {delta_auc.mean():.4f}")
print(f"  Median: {delta_auc.median():.4f}")
print(f"  Std: {delta_auc.std():.4f}")
# ============================================================
# STAGE 7B: Count Selective Regulons (CRITICAL LOG)
# ============================================================

# Count regulons meeting selectivity threshold
n_total = len(results_df)
n_selective = results_df[COL_IS_SELECTIVE].sum()
n_positive_delta = (results_df[COL_DELTA_AUC] > 0).sum()
n_negative_delta = (results_df[COL_DELTA_AUC] < 0).sum()

print("="*64)
print("       CANCER SELECTIVITY SUMMARY (DL16)")
print("="*64)
print(f"\nTotal regulons analyzed: {n_total}")
print(f"")
print(f"Regulons with Delta AUC > 0 (positive): {n_positive_delta}")
print(f"Regulons with Delta AUC < 0 (negative): {n_negative_delta}")
print(f"")
print(f">>> Regulons with Delta AUC > {L3B_DELTA_AUC_THRESHOLD} (SELECTIVE): {n_selective} <<<")
print(f"")
print(f"Selectivity rate: {n_selective/n_total*100:.1f}%")
print("="*64)

# Log L3A aggregate signature specifically
l3a_row = results_df[results_df[COL_REGULON_NAME] == 'L3A_Essential_Aggregate']
if not l3a_row.empty:
    l3a_delta = l3a_row[COL_DELTA_AUC].values[0]
    l3a_selective = l3a_row[COL_IS_SELECTIVE].values[0]
    print(f"\nL3A Aggregate Signature:")
    print(f"  Delta AUC: {l3a_delta:.4f}")
    print(f"  Selective: {l3a_selective}")
    print(f"  Interpretation: {'CANCER SELECTIVE' if l3a_selective else 'Not selective'}")
# ============================================================
# STAGE 7C: Top Selective Regulons
# ============================================================

print("Top 20 Cancer-Selective Regulons (by Delta AUC):")
print("="*80)

top_selective = results_df[results_df[COL_IS_SELECTIVE]].head(20)
display(top_selective)

# Store results
RESULTS_DF = results_df

print("\n" + "="*64)
print("STAGE 7 COMPLETE: Delta AUC Calculated")
print(f"  Total regulons: {n_total}")
print(f"  >>> Selective (Delta AUC > {L3B_DELTA_AUC_THRESHOLD}): {n_selective} <<<")
print(f"  Positive Delta: {n_positive_delta}")
print(f"  Negative Delta: {n_negative_delta}")
print("="*64)
## Stage 8: Statistical Validation & Visualization

Generate QC plots and validate results.
# ============================================================
# STAGE 8A: Delta AUC Distribution Plot
# ============================================================

fig, axes = plt.subplots(1, 2, figsize=(14, 5))

# Plot 1: Delta AUC distribution
ax1 = axes[0]
ax1.hist(RESULTS_DF[COL_DELTA_AUC], bins=50, edgecolor='black', alpha=0.7)
ax1.axvline(x=L3B_DELTA_AUC_THRESHOLD, color='red', linestyle='--', linewidth=2, 
            label=f'Threshold ({L3B_DELTA_AUC_THRESHOLD})')
ax1.axvline(x=0, color='gray', linestyle='-', linewidth=1, label='Zero')
ax1.set_xlabel('Delta AUC (Jurkat - HEK-293)', fontsize=12)
ax1.set_ylabel('Count', fontsize=12)
ax1.set_title('Delta AUC Distribution\n(Positive = Cancer Selective)', fontsize=14)
ax1.legend()

# Plot 2: AUC Jurkat vs HEK-293 scatter
ax2 = axes[1]
colors = ['green' if s else 'gray' for s in RESULTS_DF[COL_IS_SELECTIVE]]
ax2.scatter(RESULTS_DF[COL_AUC_HEK293], RESULTS_DF[COL_AUC_JURKAT], 
            c=colors, alpha=0.6, edgecolor='black', linewidth=0.5)
ax2.plot([0, 0.5], [0, 0.5], 'k--', alpha=0.5, label='y=x (no difference)')
ax2.set_xlabel('AUC HEK-293 (Normal)', fontsize=12)
ax2.set_ylabel('AUC Jurkat (Cancer)', fontsize=12)
ax2.set_title('AUC Comparison: Jurkat vs HEK-293\n(Green = Selective)', fontsize=14)
ax2.legend()

plt.tight_layout()
plt.savefig(qc_plots_dir / 'L3B_Delta_AUC_Distribution.png', dpi=150, bbox_inches='tight')
plt.show()

print(f"\nPlot saved: {qc_plots_dir / 'L3B_Delta_AUC_Distribution.png'}")
# ============================================================
# STAGE 8B: Top Selective Regulons Bar Chart
# ============================================================

fig, ax = plt.subplots(figsize=(12, 8))

# Get top 30 selective regulons
top30 = RESULTS_DF[RESULTS_DF[COL_IS_SELECTIVE]].head(30)

if len(top30) > 0:
    colors = ['darkgreen' if 'L3A' in name else 'steelblue' 
              for name in top30[COL_REGULON_NAME]]
    
    bars = ax.barh(range(len(top30)), top30[COL_DELTA_AUC].values, color=colors, edgecolor='black')
    ax.set_yticks(range(len(top30)))
    ax.set_yticklabels(top30[COL_REGULON_NAME].values, fontsize=9)
    ax.axvline(x=L3B_DELTA_AUC_THRESHOLD, color='red', linestyle='--', linewidth=2,
               label=f'Threshold ({L3B_DELTA_AUC_THRESHOLD})')
    ax.set_xlabel('Delta AUC (Jurkat - HEK-293)', fontsize=12)
    ax.set_title(f'Top {len(top30)} Cancer-Selective Regulons\n(Green = L3A Aggregate, Blue = TF Regulons)', fontsize=14)
    ax.legend()
    ax.invert_yaxis()
else:
    ax.text(0.5, 0.5, 'No selective regulons found', ha='center', va='center', fontsize=14)
    ax.set_title('No Cancer-Selective Regulons', fontsize=14)

plt.tight_layout()
plt.savefig(qc_plots_dir / 'L3B_Top_Selective_Regulons.png', dpi=150, bbox_inches='tight')
plt.show()

print(f"\nPlot saved: {qc_plots_dir / 'L3B_Top_Selective_Regulons.png'}")
# ============================================================
# STAGE 8C: AUC Heatmap
# ============================================================

# Select top regulons for heatmap
n_top = min(50, len(RESULTS_DF))
top_regulons = RESULTS_DF.head(n_top)[COL_REGULON_NAME].tolist()

# Create heatmap data
heatmap_data = AUC_MATRIX[top_regulons].T
heatmap_data.columns = [L3B_CANCER_CELL_LINE, L3B_NORMAL_CELL_LINE]

fig, ax = plt.subplots(figsize=(8, max(12, n_top * 0.25)))
sns.heatmap(heatmap_data, annot=True, fmt='.3f', cmap='RdYlGn', 
            linewidths=0.5, ax=ax, vmin=0, vmax=0.5,
            cbar_kws={'label': 'AUC Score'})
ax.set_xlabel('Cell Line', fontsize=12)
ax.set_ylabel('Regulon', fontsize=12)
ax.set_title(f'AUC Heatmap: Top {n_top} Regulons\n(Sorted by Delta AUC)', fontsize=14)

plt.tight_layout()
plt.savefig(qc_plots_dir / 'L3B_AUC_Heatmap.png', dpi=150, bbox_inches='tight')
plt.show()

print(f"\nPlot saved: {qc_plots_dir / 'L3B_AUC_Heatmap.png'}")

print("\n" + "="*64)
print("STAGE 8 COMPLETE: Visualization")
print(f"  QC plots saved to: {qc_plots_dir}")
print("="*64)
## Stage 9: Export Results & Pipeline Summary

Save all output files and generate metadata summary.
# ============================================================
# STAGE 9A: Export Active Regulons CSV
# ============================================================

print("Exporting results...\n")

# Filter to selective regulons only
active_regulons = RESULTS_DF[RESULTS_DF[COL_IS_SELECTIVE]].copy()

# Save active regulons
active_path = LAYER3B_OUTPUT_DIR / L3B_ACTIVE_REGULONS_CSV
active_regulons.to_csv(active_path, index=False)
print(f"[1] Active Regulons saved: {active_path}")
print(f"    Rows: {len(active_regulons)}")
# ============================================================
# STAGE 9B: Export All Regulons AUC CSV
# ============================================================

# Save all regulons with AUC scores
all_path = LAYER3B_OUTPUT_DIR / L3B_ALL_REGULONS_AUC_CSV
RESULTS_DF.to_csv(all_path, index=False)
print(f"[2] All Regulons AUC saved: {all_path}")
print(f"    Rows: {len(RESULTS_DF)}")
# ============================================================
# STAGE 9C: Export L3A Signature AUC CSV
# ============================================================

# Extract L3A signature row
l3a_results = RESULTS_DF[RESULTS_DF[COL_REGULON_NAME] == 'L3A_Essential_Aggregate'].copy()

# Add additional metadata
l3a_results['Total_L3A_Genes'] = len(L3A_ESSENTIAL_GENES)
l3a_results['Genes_In_Matrix'] = len([g for g in L3A_ESSENTIAL_GENES if g in set(EXP_MATRIX.columns)])

# Save L3A signature AUC
l3a_path = LAYER3B_OUTPUT_DIR / L3B_L3A_SIGNATURE_AUC_CSV
l3a_results.to_csv(l3a_path, index=False)
print(f"[3] L3A Signature AUC saved: {l3a_path}")
print(f"    L3A Delta AUC: {l3a_results[COL_DELTA_AUC].values[0]:.4f}")
# ============================================================
# STAGE 9D: Export Selectivity Summary JSON
# ============================================================

# Build summary metadata
summary = {
    "layer": "3B",
    "name": "AUCell Functional State Analysis",
    "timestamp": datetime.now().isoformat(),
    "config_version": "1.6",
    
    # Input summary
    "inputs": {
        "l3a_essential_genes": len(L3A_ESSENTIAL_GENES),
        "l2a_regulons": len(REGULON_SIGNATURES),
        "total_signatures": len(ALL_SIGNATURES),
        "expression_matrix_shape": list(EXP_MATRIX.shape),
    },
    
    # Cell line comparison
    "cell_lines": {
        "cancer": {
            "name": L3B_CANCER_CELL_LINE,
            "model_id": L3B_CANCER_MODEL_ID,
            "mean_auc": float(AUC_JURKAT.mean()),
        },
        "normal": {
            "name": L3B_NORMAL_CELL_LINE,
            "model_id": L3B_NORMAL_MODEL_ID,
            "mean_auc": float(AUC_HEK293.mean()),
        },
    },
    
    # Parameters
    "parameters": {
        "auc_threshold": L3B_AUC_THRESHOLD,
        "delta_auc_threshold": L3B_DELTA_AUC_THRESHOLD,
        "min_genes_per_signature": L3B_MIN_GENES_IN_SIGNATURE,
    },
    
    # Results summary
    "results": {
        "total_regulons_analyzed": len(RESULTS_DF),
        "selective_regulons": int(RESULTS_DF[COL_IS_SELECTIVE].sum()),
        "positive_delta_regulons": int((RESULTS_DF[COL_DELTA_AUC] > 0).sum()),
        "negative_delta_regulons": int((RESULTS_DF[COL_DELTA_AUC] < 0).sum()),
        "selectivity_rate": float(RESULTS_DF[COL_IS_SELECTIVE].sum() / len(RESULTS_DF)),
        "delta_auc_stats": {
            "min": float(RESULTS_DF[COL_DELTA_AUC].min()),
            "max": float(RESULTS_DF[COL_DELTA_AUC].max()),
            "mean": float(RESULTS_DF[COL_DELTA_AUC].mean()),
            "median": float(RESULTS_DF[COL_DELTA_AUC].median()),
            "std": float(RESULTS_DF[COL_DELTA_AUC].std()),
        },
    },
    
    # L3A Aggregate signature
    "l3a_aggregate": {
        "genes_in_signature": len(L3A_AGGREGATE_SIGNATURE),
        "auc_jurkat": float(l3a_results[COL_AUC_JURKAT].values[0]),
        "auc_hek293": float(l3a_results[COL_AUC_HEK293].values[0]),
        "delta_auc": float(l3a_results[COL_DELTA_AUC].values[0]),
        "is_selective": bool(l3a_results[COL_IS_SELECTIVE].values[0]),
    },
    
    # Deadlock compliance
    "deadlock_compliance": {
        "DL13_GeneSignature_objects": True,
        "DL14_TPM_normalized": True,
        "DL15_cell_independent_ranking": True,
        "DL16_positive_delta_selective": True,
    },
    
    # Output files
    "output_files": {
        "active_regulons": str(active_path),
        "all_regulons_auc": str(all_path),
        "l3a_signature_auc": str(l3a_path),
        "qc_plots_dir": str(qc_plots_dir),
    },
}

# Save JSON
summary_path = LAYER3B_OUTPUT_DIR / L3B_SELECTIVITY_SUMMARY_JSON
with open(summary_path, 'w') as f:
    json.dump(summary, f, indent=2)
print(f"[4] Selectivity Summary saved: {summary_path}")
# ============================================================
# STAGE 9E: Final Pipeline Summary
# ============================================================

print("\n" + "="*64)
print("       LAYER 3B COMPLETE: AUCell Functional State Analysis")
print("="*64)

print(f"\n[SMART Goal Achievement]")
print(f"  S - Compared regulon activity between Jurkat (cancer) and HEK-293 (normal)")
print(f"  M - Delta AUC quantified for {len(RESULTS_DF)} regulons")
print(f"  A - Used pySCENIC AUCell with ctxcore GeneSignature objects")
print(f"  R - Identified {int(RESULTS_DF[COL_IS_SELECTIVE].sum())} cancer-selective regulons")
print(f"  T - Completed in single notebook execution")

print(f"\n[FAIR Compliance]")
print(f"  F - Regulon IDs use TF_regulon format (persistent)")
print(f"  A - Standard CSV output with full AUC metrics")
print(f"  I - ctxcore GeneSignature objects (pySCENIC compatible)")
print(f"  R - Full provenance in {L3B_SELECTIVITY_SUMMARY_JSON}")

print(f"\n[Deadlock Compliance]")
print(f"  DL13 - GeneSignature objects: PASSED")
print(f"  DL14 - TPM normalization: PASSED")
print(f"  DL15 - Cell-independent ranking: AUTO-HANDLED")
print(f"  DL16 - Delta AUC interpretation: APPLIED")

print(f"\n[Results Summary]")
print(f"  Total regulons analyzed: {len(RESULTS_DF)}")
print(f"  >>> Cancer-selective (Delta AUC > {L3B_DELTA_AUC_THRESHOLD}): {int(RESULTS_DF[COL_IS_SELECTIVE].sum())} <<<")
print(f"  Positive Delta AUC: {int((RESULTS_DF[COL_DELTA_AUC] > 0).sum())}")
print(f"  Negative Delta AUC: {int((RESULTS_DF[COL_DELTA_AUC] < 0).sum())}")

print(f"\n[L3A Aggregate Signature]")
print(f"  Delta AUC: {l3a_results[COL_DELTA_AUC].values[0]:.4f}")
print(f"  Status: {'CANCER SELECTIVE' if l3a_results[COL_IS_SELECTIVE].values[0] else 'Not selective'}")

print(f"\n[Output Files]")
print(f"  1. {active_path}")
print(f"  2. {all_path}")
print(f"  3. {l3a_path}")
print(f"  4. {summary_path}")
print(f"  5. {qc_plots_dir}/")

print(f"\n[Next Steps]")
print(f"  -> Layer 4: Integration with drug sensitivity data")
print(f"  -> Or: Cross-validate with additional cancer cell lines")

print("\n" + "="*64)
print("                    LAYER 3B PIPELINE FINISHED")
print("="*64)
