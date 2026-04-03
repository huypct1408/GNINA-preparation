# Layer 3A: CRISPR Essentiality Validation (Physical Death)

## Stage 1: Imports & Configuration
# ============================================================
# Stage 1: Imports & Configuration
# ============================================================
import pandas as pd
import numpy as np
import json
from pathlib import Path
from datetime import datetime
import matplotlib.pyplot as plt
import seaborn as sns
import warnings
warnings.filterwarnings('ignore')

# Import central configuration
import config_system as cfg

# Print configuration summary
cfg.print_config_summary()
# ============================================================
# Check Layer 3A Resources
# ============================================================
print("\n" + "="*64)
print("LAYER 3A RESOURCE CHECK")
print("="*64)

l3a_resources = cfg.check_layer3a_resources()
all_ok = True

for key, available in l3a_resources.items():
    status = "OK" if available else "MISSING"
    symbol = "[+]" if available else "[X]"
    print(f"  {symbol} {key}: {status}")
    if not available and key in ['crispr_gene_dependency', 'layer2b_delta_network_jurkat']:
        all_ok = False

if not all_ok:
    print("\n*** CRITICAL: Missing required resources. Cannot proceed. ***")
else:
    print("\n[OK] All required resources available.")
# ============================================================
# Define Paths
# ============================================================
# Input paths
L2B_OUTPUT_DIR = cfg.PROJECT_ROOT / "outputs" / "Layer2B_Heterogeneous_RWR_Jurkat"
L2B_DELTA_CSV = L2B_OUTPUT_DIR / cfg.L2B_DELTA_NETWORK_CSV
CRISPR_CSV = cfg.CRISPR_GENE_DEPENDENCY_CSV

# Output paths
OUTPUT_DIR = cfg.LAYER3A_OUTPUT_DIR
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

QC_PLOTS_DIR = OUTPUT_DIR / cfg.L3A_QC_PLOTS_DIR
QC_PLOTS_DIR.mkdir(parents=True, exist_ok=True)

print(f"L2B Delta Network:  {L2B_DELTA_CSV}")
print(f"CRISPR Data:        {CRISPR_CSV}")
print(f"Output Directory:   {OUTPUT_DIR}")
print(f"QC Plots:           {QC_PLOTS_DIR}")
# ============================================================
# Initialize Logger
# ============================================================
logger = cfg.setup_logger("Layer3A")
logger.info("Layer 3A: CRISPR Essentiality Validation initialized")
logger.info(f"Primary cell line: {cfg.L3A_PRIMARY_CELL_LINE} ({cfg.L3A_PRIMARY_MODEL_ID})")
logger.info(f"Cancer cell lines: {list(cfg.LAYER3A_CANCER_CELL_LINES.keys())}")
logger.info(f"P(dep) threshold: > {cfg.CRISPR_PDEP_THRESHOLD}")
---
## Stage 2: Load L2B Delta Network (Filter Positive Delta, Top 60)
# ============================================================
# Stage 2: Load L2B Delta Network Summary
# ============================================================
logger.info("Loading L2B Delta Network Summary...")

df_delta = pd.read_csv(L2B_DELTA_CSV)
# After loading df_delta, add validation:
required_l2b_cols = ['Gene', 'Delta_Score', 'RWR_Active', 'RWR_Inactive', 'Is_Direct_Target']
missing_cols = [c for c in required_l2b_cols if c not in df_delta.columns]
if missing_cols:
    raise ValueError(
        f"L2B Delta Network missing required columns: {missing_cols}\n"
        f"Available columns: {list(df_delta.columns)}\n"
        f"Ensure Layer 2B completed successfully before running Layer 3A."
    )
logger.info(f"Loaded {len(df_delta):,} genes from Delta Network")

# Display summary
print("\nDelta Network Summary:")
print(f"  Total genes:        {len(df_delta):,}")
print(f"  Columns:            {list(df_delta.columns)}")
print(f"  Direct targets:     {df_delta['Is_Direct_Target'].sum()}")
print(f"\nDelta Score Statistics:")
print(f"  Min:    {df_delta['Delta_Score'].min():.6f}")
print(f"  Max:    {df_delta['Delta_Score'].max():.6f}")
print(f"  Mean:   {df_delta['Delta_Score'].mean():.6f}")
print(f"  Median: {df_delta['Delta_Score'].median():.6f}")

# Count positive vs negative delta
n_positive = (df_delta['Delta_Score'] > 0).sum()
n_negative = (df_delta['Delta_Score'] < 0).sum()
n_zero = (df_delta['Delta_Score'] == 0).sum()
print(f"\nDelta Score Distribution:")
print(f"  Positive (> 0):   {n_positive:,} ({100*n_positive/len(df_delta):.1f}%)")
print(f"  Negative (< 0):   {n_negative:,} ({100*n_negative/len(df_delta):.1f}%)")
print(f"  Zero (= 0):       {n_zero:,}")
# ============================================================
# DL10 Compliance: Filter Positive Delta ONLY
# ============================================================
logger.info("Applying DL10: Filtering positive delta scores only...")

# CRITICAL: Only positive delta (Mechanism of Killing from ACTIVE drugs)
df_positive = df_delta[df_delta['Delta_Score'] > 0].copy()
logger.info(f"Filtered to {len(df_positive):,} genes with positive delta")

# Validate DL10
cfg.validate_deadlock_rules("delta_positive_only", used_positive_only=True)
print("\n[DL10 PASS] Only positive Delta_Score values used.")
# ============================================================
# Select Top 60 by Delta Score (Descending)
# ============================================================
logger.info(f"Selecting Top {cfg.L3A_TOP_DELTA_GENES} genes by Delta_Score...")

# Sort by Delta_Score descending and take top N
df_top60 = df_positive.nlargest(cfg.L3A_TOP_DELTA_GENES, 'Delta_Score').copy()
df_top60 = df_top60.reset_index(drop=True)

# Add Delta Rank column
df_top60[cfg.COL_DELTA_RANK] = range(1, len(df_top60) + 1)

logger.info(f"Top {len(df_top60)} genes selected for CRISPR validation")

# Display top 10
print(f"\nTop 10 Positive Delta Genes (of {len(df_top60)}):")
print(df_top60[['Gene', 'Delta_Score', 'RWR_Active', 'RWR_Inactive', 'Is_Direct_Target', cfg.COL_DELTA_RANK]].head(10).to_string(index=False))
# ============================================================
# Normalize Gene Names for Matching
# ============================================================
logger.info("Normalizing gene names for cross-database matching...")

# Apply normalization (remove hyphens, uppercase)
df_top60[cfg.COL_GENE_NORMALIZED] = df_top60['Gene'].apply(cfg.normalize_gene_name)

# Check for any changes
changed = df_top60[df_top60['Gene'] != df_top60[cfg.COL_GENE_NORMALIZED]]
if len(changed) > 0:
    print(f"\nGene names normalized ({len(changed)} changes):")
    print(changed[['Gene', cfg.COL_GENE_NORMALIZED]].to_string(index=False))
else:
    print("\nNo gene name normalization needed.")
---
## Stage 3: Load CRISPR Dependency Data
# ============================================================
# Stage 3: Load CRISPR Gene Dependency Data
# ============================================================
logger.info("Loading CRISPR Gene Dependency data...")

# Load CRISPR data (large file - may take a moment)
df_crispr = pd.read_csv(CRISPR_CSV, index_col=0)
logger.info(f"Loaded CRISPR data: {df_crispr.shape[0]:,} cell lines x {df_crispr.shape[1]:,} genes")

print(f"\nCRISPR Gene Dependency Matrix:")
print(f"  Cell lines (rows): {df_crispr.shape[0]:,}")
print(f"  Genes (columns):   {df_crispr.shape[1]:,}")
print(f"  P(dep) range:      {df_crispr.min().min():.4f} to {df_crispr.max().max():.4f}")
# ============================================================
# Parse Gene Names from CRISPR Columns
# ============================================================
# CRISPR columns are formatted as "SYMBOL (EntrezID)" - need to extract SYMBOL
logger.info("Parsing gene symbols from CRISPR column headers...")

# Create gene symbol mapping
crispr_gene_map = {}
for col in df_crispr.columns:
    # Extract gene symbol (everything before the space and parentheses)
    if ' (' in col:
        gene_symbol = col.split(' (')[0]
    else:
        gene_symbol = col
    # Normalize for matching
    gene_normalized = cfg.normalize_gene_name(gene_symbol)
    crispr_gene_map[gene_normalized] = col

print(f"Parsed {len(crispr_gene_map):,} unique gene symbols from CRISPR data")
print(f"\nExample mappings:")
for i, (norm, orig) in enumerate(list(crispr_gene_map.items())[:5]):
    print(f"  {norm} -> {orig}")
# ============================================================
# DL12 Compliance: Verify HEK-293 Excluded
# ============================================================
logger.info("Validating DL12: HEK-293 exclusion...")

# Get cell lines to query (excluding HEK-293)
cancer_cell_lines = cfg.LAYER3A_CANCER_CELL_LINES
cell_lines_to_query = list(cancer_cell_lines.keys())

# Validate DL12
cfg.validate_deadlock_rules("exclude_hek293", cell_lines=cell_lines_to_query)
print(f"\n[DL12 PASS] HEK-293 excluded from CRISPR query.")
print(f"Cancer cell lines to query: {cell_lines_to_query}")
# ============================================================
# Verify Cell Lines Exist in CRISPR Data
# ============================================================
logger.info("Verifying cell lines exist in CRISPR data...")

cell_line_status = {}
for name, model_id in cancer_cell_lines.items():
    exists = model_id in df_crispr.index
    cell_line_status[name] = {
        'model_id': model_id,
        'exists': exists
    }
    status = "FOUND" if exists else "NOT FOUND"
    print(f"  {name} ({model_id}): {status}")

# Check if all cell lines found
all_found = all(s['exists'] for s in cell_line_status.values())
if not all_found:
    logger.warning("Some cell lines not found in CRISPR data!")
else:
    logger.info("All cancer cell lines found in CRISPR data.")
---
## Stage 4: Cross-Reference & Filter (5 Cancer Cell Lines)
# ============================================================
# Stage 4: Query CRISPR P(dep) for Top 50 Genes
# ============================================================
logger.info("Querying CRISPR P(dep) for Top 50 genes...")

# Initialize results DataFrame
results = df_top60.copy()

# Track gene coverage
genes_found = []
genes_not_found = []

# Query P(dep) for each gene across all cancer cell lines
for name, model_id in cancer_cell_lines.items():
    col_name = f"{cfg.COL_P_DEP_PREFIX}{name}"
    results[col_name] = np.nan
    
    if model_id not in df_crispr.index:
        logger.warning(f"Cell line {name} ({model_id}) not in CRISPR data")
        continue
    
    for idx, row in results.iterrows():
        gene_norm = row[cfg.COL_GENE_NORMALIZED]
        
        if gene_norm in crispr_gene_map:
            crispr_col = crispr_gene_map[gene_norm]
            p_dep = df_crispr.loc[model_id, crispr_col]
            results.at[idx, col_name] = p_dep
            if gene_norm not in genes_found:
                genes_found.append(gene_norm)
        else:
            if gene_norm not in genes_not_found:
                genes_not_found.append(gene_norm)

logger.info(f"Genes found in CRISPR: {len(genes_found)}/{len(results)}")
logger.info(f"Genes NOT found: {len(genes_not_found)}")

print(f"\nCRISPR Coverage:")
print(f"  Found:     {len(genes_found)} genes")
print(f"  Not found: {len(genes_not_found)} genes")
if genes_not_found:
    print(f"  Missing:   {genes_not_found}")
# ============================================================
# Add CRISPR Coverage Status Column
# ============================================================
# Determine coverage status for each gene
p_dep_cols = [f"{cfg.COL_P_DEP_PREFIX}{name}" for name in cancer_cell_lines.keys()]

def get_coverage_status(row):
    """Determine if gene was found in CRISPR data."""
    p_dep_values = [row[col] for col in p_dep_cols if col in row.index]
    if any(pd.notna(v) for v in p_dep_values):
        return cfg.CRISPR_FOUND
    return cfg.CRISPR_NOT_FOUND

results[cfg.COL_CRISPR_COVERAGE] = results.apply(get_coverage_status, axis=1)

# Count coverage
n_found = (results[cfg.COL_CRISPR_COVERAGE] == cfg.CRISPR_FOUND).sum()
n_not_found = (results[cfg.COL_CRISPR_COVERAGE] == cfg.CRISPR_NOT_FOUND).sum()
print(f"\nCRISPR Coverage Summary:")
print(f"  FOUND:        {n_found} genes")
print(f"  NOT_IN_CRISPR: {n_not_found} genes")
# ============================================================
# Calculate Max P(dep) and Essential Status
# ============================================================
logger.info("Calculating Max P(dep) and essential status...")

# Calculate max P(dep) across all cancer cell lines
results[cfg.COL_MAX_P_DEP] = results[p_dep_cols].max(axis=1)

# Count how many cell lines show essentiality
def count_essential_lines(row):
    """Count cell lines where P(dep) > threshold."""
    count = 0
    for col in p_dep_cols:
        if pd.notna(row[col]) and row[col] > cfg.CRISPR_PDEP_THRESHOLD:
            count += 1
    return count

results[cfg.COL_ESSENTIAL_IN_N_LINES] = results.apply(count_essential_lines, axis=1)

# DL11 Compliance: Essential = P(dep) > 0.8 in ANY cancer cell line
results[cfg.COL_IS_ESSENTIAL] = results[cfg.COL_MAX_P_DEP] > cfg.CRISPR_PDEP_THRESHOLD

# Validate DL11
cfg.validate_deadlock_rules("pdep_threshold", threshold=cfg.CRISPR_PDEP_THRESHOLD)
print(f"\n[DL11 PASS] P(dep) > {cfg.CRISPR_PDEP_THRESHOLD} threshold applied.")
# ============================================================
# DL9 Compliance: Verify Cell Line Match
# ============================================================
logger.info("Validating DL9: Cell line match...")

# L2B was run with Jurkat (from directory name)
l2b_cell_line = "Jurkat"  # Derived from Layer2B_Heterogeneous_RWR_Jurkat
l3a_primary = cfg.L3A_PRIMARY_CELL_LINE

# Validate DL9
cfg.validate_deadlock_rules("cell_line_match", 
                           l2b_cell_line=l2b_cell_line, 
                           l3a_primary_cell_line=l3a_primary)
print(f"\n[DL9 PASS] L2B cell line ({l2b_cell_line}) matches L3A primary ({l3a_primary}).")
# ============================================================
# Summary Statistics
# ============================================================
print("\n" + "="*64)
print("LAYER 3A CRISPR ESSENTIALITY SUMMARY")
print("="*64)

n_essential = results[cfg.COL_IS_ESSENTIAL].sum()
n_found_genes = (results[cfg.COL_CRISPR_COVERAGE] == cfg.CRISPR_FOUND).sum()

print(f"\nInput:")
print(f"  Top positive-delta genes queried: {len(results)}")
print(f"  Cancer cell lines queried:        {len(cancer_cell_lines)}")

print(f"\nCRISPR Coverage:")
print(f"  Genes found in CRISPR:     {n_found_genes}/{len(results)} ({100*n_found_genes/len(results):.1f}%)")
print(f"  Genes NOT in CRISPR:       {len(results) - n_found_genes}")

print(f"\nEssentiality (P(dep) > {cfg.CRISPR_PDEP_THRESHOLD}):")
print(f"  Essential genes:           {n_essential}")
print(f"  Non-essential genes:       {n_found_genes - n_essential}")
print(f"  Validation rate:           {100*n_essential/n_found_genes:.1f}% of CRISPR-covered genes" if n_found_genes > 0 else "  N/A")

# Per cell line breakdown
print(f"\nPer Cell Line Breakdown:")
for name in cancer_cell_lines.keys():
    col = f"{cfg.COL_P_DEP_PREFIX}{name}"
    n_essential_line = (results[col] > cfg.CRISPR_PDEP_THRESHOLD).sum()
    mean_pdep = results[col].mean()
    primary_marker = " [PRIMARY]" if name == cfg.L3A_PRIMARY_CELL_LINE else ""
    print(f"  {name}{primary_marker}: {n_essential_line} essential genes (mean P(dep) = {mean_pdep:.3f})")
# ============================================================
# Display Essential Genes
# ============================================================
print("\n" + "="*64)
print("ESSENTIAL GENES (P(dep) > 0.8 in any cancer line)")
print("="*64)

df_essential = results[results[cfg.COL_IS_ESSENTIAL]].copy()

if len(df_essential) > 0:
    display_cols = ['Gene', 'Delta_Score', cfg.COL_MAX_P_DEP, cfg.COL_ESSENTIAL_IN_N_LINES, 'Is_Direct_Target']
    # Add cell-line specific columns
    for name in cancer_cell_lines.keys():
        col = f"{cfg.COL_P_DEP_PREFIX}{name}"
        display_cols.append(col)
    
    print(df_essential[display_cols].to_string(index=False))
else:
    print("\nNo genes met the P(dep) > 0.8 threshold.")
    print("\nTop 5 genes by Max P(dep):")
    display_cols = ['Gene', 'Delta_Score', cfg.COL_MAX_P_DEP, cfg.COL_ESSENTIAL_IN_N_LINES]
    print(results.nlargest(5, cfg.COL_MAX_P_DEP)[display_cols].to_string(index=False))
---
## Stage 5: Generate Outputs & Visualizations
# ============================================================
# Stage 5A: Save Essential Targets CSV
# ============================================================
logger.info("Saving essential targets...")

# Filter to essential genes only
df_essential_out = results[results[cfg.COL_IS_ESSENTIAL]].copy()

# Reorder columns for output
output_cols = [
    'Gene', cfg.COL_GENE_NORMALIZED, cfg.COL_DELTA_RANK, 'Delta_Score',
    'RWR_Active', 'RWR_Inactive', 'Is_Direct_Target',
    cfg.COL_MAX_P_DEP, cfg.COL_ESSENTIAL_IN_N_LINES, cfg.COL_IS_ESSENTIAL,
    cfg.COL_CRISPR_COVERAGE
]
# Add cell-line specific P(dep) columns
for name in cancer_cell_lines.keys():
    output_cols.append(f"{cfg.COL_P_DEP_PREFIX}{name}")

df_essential_out = df_essential_out[output_cols]

# Save
essential_path = OUTPUT_DIR / cfg.L3A_ESSENTIAL_TARGETS_CSV
df_essential_out.to_csv(essential_path, index=False)
logger.info(f"Saved {len(df_essential_out)} essential targets to {essential_path}")
print(f"\n[SAVED] Essential targets: {essential_path}")
print(f"        {len(df_essential_out)} genes")
# ============================================================
# Stage 5B: Save All Candidates CSV
# ============================================================
logger.info("Saving all candidates...")

# Save all Top 50 with P(dep) values
df_all_out = results[output_cols].copy()

# Save
all_path = OUTPUT_DIR / cfg.L3A_ALL_CANDIDATES_CSV
df_all_out.to_csv(all_path, index=False)
logger.info(f"Saved {len(df_all_out)} candidates to {all_path}")
print(f"\n[SAVED] All candidates: {all_path}")
print(f"        {len(df_all_out)} genes")
# ============================================================
# Stage 5C: Save Validation Summary JSON
# ============================================================
logger.info("Generating validation summary...")

# Build summary JSON
validation_summary = {
    "metadata": {
        "layer": "3A",
        "name": "CRISPR Essentiality Validation",
        "timestamp": datetime.now().isoformat(),
        "config_version": "1.5",
        "primary_cell_line": cfg.L3A_PRIMARY_CELL_LINE,
        "primary_model_id": cfg.L3A_PRIMARY_MODEL_ID,
    },
    "parameters": {
        "top_delta_genes": cfg.L3A_TOP_DELTA_GENES,
        "delta_filter": "positive_only (DL10)",
        "pdep_threshold": cfg.CRISPR_PDEP_THRESHOLD,
        "cancer_cell_lines": list(cancer_cell_lines.keys()),
        "hek293_excluded": True,
    },
    "input_files": {
        "l2b_delta_network": str(L2B_DELTA_CSV),
        "crispr_gene_dependency": str(CRISPR_CSV),
    },
    "statistics": {
        "total_genes_queried": len(results),
        "genes_in_crispr": int(n_found_genes),
        "genes_not_in_crispr": len(results) - int(n_found_genes),
        "essential_genes": int(n_essential),
        "non_essential_genes": int(n_found_genes) - int(n_essential),
        "validation_rate": float(n_essential / n_found_genes) if n_found_genes > 0 else 0.0,
    },
    "per_cell_line": {},
    "deadlock_compliance": {
        "DL9_cell_line_match": True,
        "DL10_delta_positive_only": True,
        "DL11_pdep_threshold": True,
        "DL12_hek293_excluded": True,
    },
    "output_files": {
        "essential_targets": str(essential_path),
        "all_candidates": str(all_path),
    },
    "essential_gene_list": df_essential_out['Gene'].tolist() if len(df_essential_out) > 0 else [],
    "genes_not_in_crispr": genes_not_found,
}

# Add per-cell-line stats
for name in cancer_cell_lines.keys():
    col = f"{cfg.COL_P_DEP_PREFIX}{name}"
    validation_summary["per_cell_line"][name] = {
        "model_id": cancer_cell_lines[name],
        "essential_count": int((results[col] > cfg.CRISPR_PDEP_THRESHOLD).sum()),
        "mean_pdep": float(results[col].mean()) if pd.notna(results[col].mean()) else None,
        "max_pdep": float(results[col].max()) if pd.notna(results[col].max()) else None,
    }

# Save JSON
summary_path = OUTPUT_DIR / cfg.L3A_VALIDATION_SUMMARY_JSON
with open(summary_path, 'w') as f:
    json.dump(validation_summary, f, indent=2)

logger.info(f"Saved validation summary to {summary_path}")
print(f"\n[SAVED] Validation summary: {summary_path}")
# ============================================================
# Stage 5D: Visualization - P(dep) Distribution
# ============================================================
logger.info("Generating QC visualizations...")

fig, axes = plt.subplots(2, 2, figsize=(14, 10))

# Plot 1: Max P(dep) distribution
ax1 = axes[0, 0]
data_for_hist = results[results[cfg.COL_CRISPR_COVERAGE] == cfg.CRISPR_FOUND][cfg.COL_MAX_P_DEP].dropna()
ax1.hist(data_for_hist, bins=20, edgecolor='black', alpha=0.7, color='steelblue')
ax1.axvline(x=cfg.CRISPR_PDEP_THRESHOLD, color='red', linestyle='--', linewidth=2, label=f'Threshold ({cfg.CRISPR_PDEP_THRESHOLD})')
ax1.set_xlabel('Max P(dep) across cancer lines')
ax1.set_ylabel('Number of genes')
ax1.set_title('Distribution of Max P(dependency)')
ax1.legend()

# Plot 2: P(dep) by cell line (boxplot)
ax2 = axes[0, 1]
pdep_data = []
pdep_labels = []
for name in cancer_cell_lines.keys():
    col = f"{cfg.COL_P_DEP_PREFIX}{name}"
    values = results[col].dropna().values
    if len(values) > 0:
        pdep_data.append(values)
        pdep_labels.append(name)

if pdep_data:
    bp = ax2.boxplot(pdep_data, labels=pdep_labels, patch_artist=True)
    colors = ['#FF9999' if l == cfg.L3A_PRIMARY_CELL_LINE else '#99CCFF' for l in pdep_labels]
    for patch, color in zip(bp['boxes'], colors):
        patch.set_facecolor(color)
    ax2.axhline(y=cfg.CRISPR_PDEP_THRESHOLD, color='red', linestyle='--', linewidth=2)
    ax2.set_ylabel('P(dependency)')
    ax2.set_title('P(dep) Distribution by Cell Line')
    ax2.tick_params(axis='x', rotation=45)

# Plot 3: Delta Score vs Max P(dep) scatter
ax3 = axes[1, 0]
df_plot = results[results[cfg.COL_CRISPR_COVERAGE] == cfg.CRISPR_FOUND].copy()
colors = ['red' if x else 'blue' for x in df_plot[cfg.COL_IS_ESSENTIAL]]
ax3.scatter(df_plot['Delta_Score'], df_plot[cfg.COL_MAX_P_DEP], c=colors, alpha=0.7, edgecolors='black')
ax3.axhline(y=cfg.CRISPR_PDEP_THRESHOLD, color='red', linestyle='--', linewidth=2, label='Essential threshold')
ax3.set_xlabel('Delta Score (RWR_Active - RWR_Inactive)')
ax3.set_ylabel('Max P(dependency)')
ax3.set_title('Delta Score vs CRISPR Essentiality')
ax3.legend(['Essential threshold', 'Essential', 'Non-essential'])

# Annotate top genes
for _, row in df_plot.nlargest(3, cfg.COL_MAX_P_DEP).iterrows():
    ax3.annotate(row['Gene'], (row['Delta_Score'], row[cfg.COL_MAX_P_DEP]), 
                fontsize=8, ha='left')

# Plot 4: Essential genes per cell line (bar chart)
ax4 = axes[1, 1]
cell_line_counts = []
for name in cancer_cell_lines.keys():
    col = f"{cfg.COL_P_DEP_PREFIX}{name}"
    count = (results[col] > cfg.CRISPR_PDEP_THRESHOLD).sum()
    cell_line_counts.append(count)

bar_colors = ['#FF6B6B' if n == cfg.L3A_PRIMARY_CELL_LINE else '#4ECDC4' for n in cancer_cell_lines.keys()]
bars = ax4.bar(list(cancer_cell_lines.keys()), cell_line_counts, color=bar_colors, edgecolor='black')
ax4.set_ylabel('Number of Essential Genes')
ax4.set_title(f'Essential Genes per Cell Line (P(dep) > {cfg.CRISPR_PDEP_THRESHOLD})')
ax4.tick_params(axis='x', rotation=45)

# Add count labels on bars
for bar, count in zip(bars, cell_line_counts):
    ax4.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.1, 
            str(count), ha='center', va='bottom', fontsize=10)

plt.tight_layout()

# Save figure
fig_path = QC_PLOTS_DIR / "L3A_CRISPR_Essentiality_Overview.png"
plt.savefig(fig_path, dpi=150, bbox_inches='tight')
logger.info(f"Saved QC plot to {fig_path}")
print(f"\n[SAVED] QC plot: {fig_path}")
plt.show()
# ============================================================
# Stage 5E: Heatmap - P(dep) across Cell Lines
# ============================================================
logger.info("Generating P(dep) heatmap...")

# Prepare data for heatmap (top 30 genes by Max P(dep))
df_heatmap = results[results[cfg.COL_CRISPR_COVERAGE] == cfg.CRISPR_FOUND].nlargest(30, cfg.COL_MAX_P_DEP).copy()

if len(df_heatmap) > 0:
    # Extract P(dep) columns
    heatmap_data = df_heatmap.set_index('Gene')[p_dep_cols]
    heatmap_data.columns = [col.replace(cfg.COL_P_DEP_PREFIX, '') for col in heatmap_data.columns]
    
    # Create heatmap
    fig, ax = plt.subplots(figsize=(10, 12))
    sns.heatmap(heatmap_data, annot=True, fmt='.2f', cmap='RdYlBu_r', 
                center=cfg.CRISPR_PDEP_THRESHOLD, vmin=0, vmax=1,
                linewidths=0.5, ax=ax, cbar_kws={'label': 'P(dependency)'})
    ax.set_title(f'P(dependency) Heatmap - Top 30 Genes by Max P(dep)\n(Threshold: {cfg.CRISPR_PDEP_THRESHOLD})')
    ax.set_xlabel('Cancer Cell Line')
    ax.set_ylabel('Gene')
    
    plt.tight_layout()
    
    # Save figure
    heatmap_path = QC_PLOTS_DIR / "L3A_Pdep_Heatmap.png"
    plt.savefig(heatmap_path, dpi=150, bbox_inches='tight')
    logger.info(f"Saved heatmap to {heatmap_path}")
    print(f"\n[SAVED] Heatmap: {heatmap_path}")
    plt.show()
else:
    print("\nNo genes with CRISPR data available for heatmap.")
---
## Stage 6: Pipeline Summary & Next Steps
# ============================================================
# Stage 6: Final Pipeline Summary
# ============================================================
print("\n" + "="*64)
print("LAYER 3A CRISPR ESSENTIALITY - PIPELINE COMPLETE")
print("="*64)

print(f"\n[INPUT]")
print(f"  L2B Delta Network: {L2B_DELTA_CSV.name}")
print(f"  CRISPR Data:       {CRISPR_CSV.name}")
print(f"  Top genes queried: {len(results)} (positive delta only)")
print(f"  Cell lines:        {', '.join(cancer_cell_lines.keys())}")

print(f"\n[OUTPUT]")
print(f"  Essential targets: {essential_path.name} ({len(df_essential_out)} genes)")
print(f"  All candidates:    {all_path.name} ({len(df_all_out)} genes)")
print(f"  Summary JSON:      {summary_path.name}")
print(f"  QC Plots:          {QC_PLOTS_DIR}")

print(f"\n[RESULTS]")
print(f"  CRISPR coverage:   {n_found_genes}/{len(results)} genes ({100*n_found_genes/len(results):.1f}%)")
print(f"  Essential genes:   {n_essential} (P(dep) > {cfg.CRISPR_PDEP_THRESHOLD})")
print(f"  Validation rate:   {100*n_essential/n_found_genes:.1f}%" if n_found_genes > 0 else "  N/A")

print(f"\n[DEADLOCK COMPLIANCE]")
print(f"  DL9  (Cell line match):      PASS")
print(f"  DL10 (Delta positive only):  PASS")
print(f"  DL11 (P(dep) threshold):     PASS")
print(f"  DL12 (HEK-293 excluded):     PASS")

if len(df_essential_out) > 0:
    print(f"\n[ESSENTIAL GENES]")
    for _, row in df_essential_out.iterrows():
        direct = "*" if row['Is_Direct_Target'] else ""
        print(f"  {direct}{row['Gene']}: Max P(dep) = {row[cfg.COL_MAX_P_DEP]:.3f} (in {row[cfg.COL_ESSENTIAL_IN_N_LINES]} cell lines)")
    print("  (* = Direct drug target)")

print(f"\n[NEXT STEPS]")
print(f"  -> Layer 3B: AUCell Selectivity Analysis")
print(f"     - Compare regulon activity: Cancer vs HEK-293")
print(f"     - Identify cancer-selective essential genes")
print(f"     - Generate final MoA target list")

print("\n" + "="*64)
logger.info("Layer 3A pipeline completed successfully.")
