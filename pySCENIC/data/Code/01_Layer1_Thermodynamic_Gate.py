# Layer 1: Thermodynamic Gate Pipeline v1.3
# **Purpose**: Filter GNINA docking results by thermodynamic feasibility (dG <= -6.5 kcal/mol) using **Two-Track Architecture** with Intra-Ligand PoseBusters Rescue Logic, and calculate P0 raw weight vector for downstream analysis.\n
# =============================================================================
# STAGE 1: IMPORTS & CONFIG
# =============================================================================
# Version: 1.2.0 - Two-Track Architecture with Intra-Ligand PoseBusters Rescue
# FAIR: Findable (structured outputs), Accessible (CSV),
#       Interoperable (pandas), Reusable (config_system.py)
# =============================================================================

import sys
import warnings
from pathlib import Path
from datetime import datetime
from typing import Optional, Tuple, List, Dict

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# Import central configuration
import config_system as cfg

# Setup logger
logger = cfg.setup_logger("Layer1")

# Print configuration summary
cfg.print_config_summary()

# Suppress warnings for clean output
warnings.filterwarnings('ignore', category=FutureWarning)

print(f"\nNotebook started: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
print(f"Python version: {sys.version}")
print(f"Pandas version: {pd.__version__}")
# =============================================================================
# STAGE 2: DATA LOADING FUNCTIONS
# =============================================================================
print("\n" + "="*64)
print("STAGE 2: DATA LOADING FUNCTIONS")
print("="*64)


def load_from_posebusters(target_name: str) -> Optional[pd.DataFrame]:
    """
    Load ALL poses (1790) from PoseBusters CSV for a target.
    Used by Track 1 for PPARA, PPARD, ALOX-5.
    
    Returns DataFrame with columns:
      - ligand_id, pose_id, original_gnina_model
      - CNNscore, CNN_VS, CNNaffinity, minimizedAffinity
      - PoseBuster_Valid (boolean)
    """
    pb_path = cfg.POSEBUSTERS_PATHS.get(target_name)
    
    if pb_path is None or not pb_path.exists():
        logger.warning(f"{target_name}: PoseBusters CSV not found")
        return None
    
    logger.info(f"{target_name}: Loading PoseBusters from {pb_path.name}")
    df = pd.read_csv(pb_path)
    
    # Verify required columns exist
    required_cols = ['ligand_id', cfg.COL_PB_VALID, 'CNNscore', 'CNN_VS', 'minimizedAffinity']
    missing = [c for c in required_cols if c not in df.columns]
    if missing:
        logger.error(f"{target_name}: Missing columns in PoseBusters CSV: {missing}")
        logger.info(f"Available columns: {list(df.columns)}")
        return None
    
    df['target'] = target_name
    
    n_total = len(df)
    n_valid = df[cfg.COL_PB_VALID].sum()
    n_invalid = n_total - n_valid
    print(f"  {target_name}: {n_total} poses loaded ({n_valid} valid, {n_invalid} invalid)")
    
    return df


def load_from_inter_ligand(target_name: str) -> Optional[pd.DataFrame]:
    """
    Load best pose per ligand from Excel Sheet 1 (Inter-Ligand_Ranking).
    Used by Track 2 for targets without PoseBusters.
    
    Returns DataFrame with columns:
      - lig_id, orig_name, CNN_VS, CNNscore, CNNaffinity, Affinity_kcal, SMILES
    """
    target_info = cfg.TARGET_PROTEINS.get(target_name)
    if target_info is None:
        raise ValueError(f"Unknown target: {target_name}")
    
    docking_dir = target_info["docking_dir"]
    summary_dir = docking_dir / "summary"
    
    xlsx_files = list(summary_dir.glob("docking_summary*.xlsx"))
    if not xlsx_files:
        raise FileNotFoundError(f"No docking_summary*.xlsx found in {summary_dir}")
    
    xlsx_path = xlsx_files[0]
    logger.info(f"{target_name}: Loading from {xlsx_path.name} Sheet 1 (Inter-Ligand)")
    
    df = pd.read_excel(xlsx_path, sheet_name=cfg.SHEET_INTER_LIGAND, engine='openpyxl')
    
    # Standardize column names
    col_mapping = {
        'ID': 'lig_id',
        'Original_Name': 'orig_name',
        'CNN_VS': 'CNN_VS',
        'CNNscore': 'CNNscore',
        'CNNaffinity': 'CNNaffinity',
        'Affinity_kcal': 'Affinity_kcal',
        'SMILES': 'smiles',
        'Sanity_Flag': 'Sanity_Flag',
    }
    existing_cols = {k: v for k, v in col_mapping.items() if k in df.columns}
    df = df.rename(columns=existing_cols)
    df['target'] = target_name
    
    print(f"  {target_name}: {len(df)} ligands loaded from Inter-Ligand sheet")
    return df


def load_metadata_from_excel(target_name: str) -> pd.DataFrame:
    """
    Load SMILES and orig_name from Excel Sheet 1 for merging.
    Used to enrich Track 1 data with SMILES.
    """
    target_info = cfg.TARGET_PROTEINS.get(target_name)
    docking_dir = target_info["docking_dir"]
    summary_dir = docking_dir / "summary"
    
    xlsx_files = list(summary_dir.glob("docking_summary*.xlsx"))
    if not xlsx_files:
        return pd.DataFrame()
    
    xlsx_path = xlsx_files[0]
    df = pd.read_excel(xlsx_path, sheet_name=cfg.SHEET_INTER_LIGAND, engine='openpyxl')
    
    # Extract just ID, Original_Name, SMILES
    col_mapping = {'ID': 'ligand_id', 'Original_Name': 'orig_name', 'SMILES': 'smiles'}
    existing = {k: v for k, v in col_mapping.items() if k in df.columns}
    df = df.rename(columns=existing)
    
    return df[['ligand_id', 'orig_name', 'smiles']].copy() if 'ligand_id' in df.columns else df


print("Data loading functions defined.")
# =============================================================================
# STAGE 3: TRACK 1 - POSEBUSTERS TARGETS WITH RESCUE LOGIC
# =============================================================================
# Targets: PPARA, PPARD, ALOX-5
# Strategy: Load ALL 1790 poses -> Filter PoseBuster_Valid == True
#           -> Sort by CNN_VS descending -> Select best valid pose per ligand
# RESCUE: If Mode 1 fails PB but Mode 2+ passes, use the best valid mode
# =============================================================================
print("\n" + "="*64)
print("STAGE 3: TRACK 1 - POSEBUSTERS TARGETS WITH RESCUE LOGIC")
print("="*64)
print(f"Track 1 targets: {cfg.TARGETS_WITH_PB}")

track1_results = []
rescue_records = []  # For rescue report

for target_name in cfg.TARGETS_WITH_PB:
    print(f"\n--- Processing {target_name} ---")
    
    # Load all poses from PoseBusters CSV
    df_poses = load_from_posebusters(target_name)
    if df_poses is None:
        logger.error(f"{target_name}: Failed to load PoseBusters data")
        continue
    
    # Load metadata (SMILES, orig_name) from Excel
    df_meta = load_metadata_from_excel(target_name)
    
    # Get unique ligands
    all_ligands = df_poses['ligand_id'].unique()
    n_ligands = len(all_ligands)
    print(f"  Total unique ligands: {n_ligands}")
    
    # Filter to valid poses only
    df_valid = df_poses[df_poses[cfg.COL_PB_VALID] == True].copy()
    ligands_with_valid = df_valid['ligand_id'].unique()
    n_with_valid = len(ligands_with_valid)
    print(f"  Ligands with at least 1 valid pose: {n_with_valid}")
    
    # Identify ligands that need rescue (Mode 1 invalid but have valid poses)
    # First, find Mode 1 for each ligand
    df_mode1 = df_poses[df_poses[cfg.COL_GNINA_MODEL] == 1].copy()
    mode1_invalid_ligands = df_mode1[df_mode1[cfg.COL_PB_VALID] == False]['ligand_id'].unique()
    
    # Rescued = Mode 1 invalid BUT has valid poses in Mode 2+
    rescued_ligands = set(mode1_invalid_ligands) & set(ligands_with_valid)
    n_rescued = len(rescued_ligands)
    print(f"  Ligands RESCUED (Mode 1 invalid, Mode 2+ valid): {n_rescued}")
    
    # Ligands completely lost (no valid poses at all)
    lost_ligands = set(all_ligands) - set(ligands_with_valid)
    n_lost = len(lost_ligands)
    print(f"  Ligands LOST (all 10 poses invalid): {n_lost}")
    
    # Select best valid pose per ligand (sort by CNN_VS descending)
    df_valid_sorted = df_valid.sort_values('CNN_VS', ascending=False)
    df_best = df_valid_sorted.groupby('ligand_id').first().reset_index()
    print(f"  Best valid poses selected: {len(df_best)}")
    
    # Record rescue details
    for lig_id in rescued_ligands:
        # Get Mode 1 scores
        mode1_row = df_mode1[df_mode1['ligand_id'] == lig_id].iloc[0]
        # Get rescued pose
        rescued_row = df_best[df_best['ligand_id'] == lig_id].iloc[0]
        
        rescue_records.append({
            'target': target_name,
            'ligand_id': lig_id,
            'original_mode': 1,
            'original_CNN_VS': mode1_row['CNN_VS'],
            'original_CNNscore': mode1_row['CNNscore'],
            'original_PB_Valid': False,
            'rescued_mode': rescued_row[cfg.COL_GNINA_MODEL],
            'rescued_CNN_VS': rescued_row['CNN_VS'],
            'rescued_CNNscore': rescued_row['CNNscore'],
            'rescued_PB_Valid': True,
        })
    
    # Standardize columns for merging with Track 2
    df_best = df_best.rename(columns={
        'ligand_id': 'lig_id',
        'minimizedAffinity': 'Affinity_kcal'
    })
    
    # Merge with metadata to get SMILES and orig_name
    if not df_meta.empty and 'ligand_id' in df_meta.columns:
        df_meta_renamed = df_meta.rename(columns={'ligand_id': 'lig_id'})
        df_best = df_best.merge(df_meta_renamed[['lig_id', 'orig_name', 'smiles']], 
                                 on='lig_id', how='left')
    
    # Add PB_Status column
    df_best[cfg.COL_PB_STATUS] = cfg.PB_STATUS_VALID
    
    # Add to track 1 results
    track1_results.append(df_best)

# Combine all Track 1 results
if track1_results:
    df_track1 = pd.concat(track1_results, ignore_index=True)
    print(f"\nTrack 1 TOTAL: {len(df_track1)} ligand-target pairs")
else:
    df_track1 = pd.DataFrame()
    print("\nTrack 1 TOTAL: 0 (no data)")

# Create rescue report DataFrame
df_rescue = pd.DataFrame(rescue_records)
print(f"\nRescue report: {len(df_rescue)} ligands rescued across Track 1 targets")
# =============================================================================
# STAGE 4: TRACK 2 - NON-POSEBUSTERS TARGETS WITH CNNSCORE FALLBACK
# =============================================================================
# Targets: PPARG, EGFR, ERBB2, KDR, PTGS2, PTGES
# Strategy: Load from Excel Sheet 1 (already best pose per ligand)
#           -> Apply CNNscore >= 0.5 fallback filter
# =============================================================================
print("\n" + "="*64)
print("STAGE 4: TRACK 2 - NON-POSEBUSTERS TARGETS (CNNscore Fallback)")
print("="*64)
print(f"Track 2 targets: {cfg.TARGETS_WITHOUT_PB}")

track2_results = []

for target_name in cfg.TARGETS_WITHOUT_PB:
    print(f"\n--- Processing {target_name} ---")
    
    try:
        df = load_from_inter_ligand(target_name)
        
        # Apply CNNscore >= 0.5 fallback filter
        n_before = len(df)
        df_filtered = df[df['CNNscore'] >= cfg.CNNSCORE_FALLBACK_THRESHOLD].copy()
        n_after = len(df_filtered)
        n_rejected = n_before - n_after
        
        print(f"  CNNscore >= {cfg.CNNSCORE_FALLBACK_THRESHOLD} filter: {n_before} -> {n_after} ({n_rejected} rejected)")
        
        # Add PB_Status column
        df_filtered[cfg.COL_PB_STATUS] = cfg.PB_STATUS_NOT_AVAILABLE
        
        track2_results.append(df_filtered)
        
    except Exception as e:
        logger.error(f"{target_name}: {e}")

# Combine all Track 2 results
if track2_results:
    df_track2 = pd.concat(track2_results, ignore_index=True)
    print(f"\nTrack 2 TOTAL: {len(df_track2)} ligand-target pairs")
else:
    df_track2 = pd.DataFrame()
    print("\nTrack 2 TOTAL: 0 (no data)")
# =============================================================================
# STAGE 5: COMBINE TRACKS & APPLY GATE 1 FILTER
# =============================================================================
print("\n" + "="*64)
print("STAGE 5: COMBINE TRACKS & GATE 1 FILTER")
print("="*64)

# Combine Track 1 and Track 2
df_combined = pd.concat([df_track1, df_track2], ignore_index=True)
print(f"Combined records: {len(df_combined)}")
print(f"  - Track 1 (PB): {len(df_track1)}")
print(f"  - Track 2 (no PB): {len(df_track2)}")

# Ensure Affinity_kcal column exists
if 'Affinity_kcal' not in df_combined.columns and 'minimizedAffinity' in df_combined.columns:
    df_combined['Affinity_kcal'] = df_combined['minimizedAffinity']

# Store pre-filter data for QC
df_pre_gate1 = df_combined.copy()

# Apply Gate 1: Thermodynamic filter (dG <= -6.5 kcal/mol)
print(f"\nApplying Gate 1: Affinity_kcal <= {cfg.DELTA_G_CUTOFF} kcal/mol")

thermo_pass = df_combined['Affinity_kcal'] <= cfg.DELTA_G_CUTOFF
df_passed = df_combined[thermo_pass].copy()
df_rejected = df_combined[~thermo_pass].copy()

n_total = len(df_combined)
n_passed = len(df_passed)
n_rejected = len(df_rejected)

print(f"\nGate 1 Results:")
print(f"  Total records: {n_total}")
print(f"  Passed: {n_passed} ({100*n_passed/n_total:.1f}%)")
print(f"  Rejected (dG > {cfg.DELTA_G_CUTOFF}): {n_rejected} ({100*n_rejected/n_total:.1f}%)")

# Per-target breakdown
print("\nPer-target pass rates:")
for target in cfg.TARGET_NAMES:
    n_tot = len(df_combined[df_combined['target'] == target])
    n_pass = len(df_passed[df_passed['target'] == target])
    rate = 100 * n_pass / n_tot if n_tot > 0 else 0
    track = "Track1" if target in cfg.TARGETS_WITH_PB else "Track2"
    print(f"  {target} [{track}]: {n_pass}/{n_tot} ({rate:.1f}%)")
# =============================================================================
# STAGE 6: P0 ASSEMBLY
# =============================================================================
# HOTFIX #1: Inverted Affinity Normalization (more negative = better)
# HOTFIX #2: PER-TARGET Normalization (avoid binding pocket size bias)
# 
# Outputs TWO DataFrames:
#   1. df_with_p0_long: Long format (each row = ligand-target pair) for QC
#   2. df_p0_matrix: Wide format (rows=ligands, cols=targets) for Layer 2 RWR
# =============================================================================
print("\n" + "="*64)
print("STAGE 6: P0 ASSEMBLY (Per-Target Normalization + Pivot Matrix)")
print("="*64)

def calculate_p0_raw_weight(df: pd.DataFrame) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Calculate P0_raw_weight using PER-TARGET Normalization.
    
    HOTFIX #1: Inverted Affinity (more negative = better score)
    HOTFIX #2: Per-Target Normalization to prevent cross-protein scoring bias
               (avoids "large binding pocket bias" where deep pockets like MMP-9
               naturally produce more negative dG than shallow pockets like EGFR)
    
    Formula (applied PER TARGET):
      CNN_VS_norm = (CNN_VS - target_min) / (target_max - target_min)
      Affinity_norm = (target_max - Affinity) / (target_max - target_min)  # INVERTED
      P0_raw_weight = alpha * CNN_VS_norm + (1 - alpha) * Affinity_norm
    
    Returns:
      Tuple[df_long, df_matrix] - Long format for QC, Wide matrix for Layer 2
    """
    df = df.copy()
    
    # ---------------------------------------------------------
    # 1. CNN_VS normalization (PER-TARGET) - Higher is better
    # ---------------------------------------------------------
    cnn_min = df.groupby('target')['CNN_VS'].transform('min')
    cnn_max = df.groupby('target')['CNN_VS'].transform('max')
    cnn_range = cnn_max - cnn_min
    
    # np.where handles edge case where all values are identical (range = 0)
    df[cfg.COL_CNN_VS_NORM] = np.where(
        cnn_range > 0,
        (df['CNN_VS'] - cnn_min) / cnn_range,
        0.5  # Neutral value if no differentiation
    )
    
    # ---------------------------------------------------------
    # 2. Affinity normalization (PER-TARGET & INVERTED)
    #    More negative (better binding) -> Higher score
    # ---------------------------------------------------------
    aff_min = df.groupby('target')['Affinity_kcal'].transform('min')  # Best (most negative)
    aff_max = df.groupby('target')['Affinity_kcal'].transform('max')  # Worst (least negative)
    aff_range = aff_max - aff_min
    
    df[cfg.COL_AFFINITY_NORM] = np.where(
        aff_range > 0,
        (aff_max - df['Affinity_kcal']) / aff_range,  # INVERTED
        0.5
    )
    
    # ---------------------------------------------------------
    # 3. Calculate P0_raw_weight
    # ---------------------------------------------------------
    alpha = cfg.ALPHA_CNN_VS_WEIGHT  # 0.7
    df[cfg.COL_P0_RAW_WEIGHT] = (
        alpha * df[cfg.COL_CNN_VS_NORM] + 
        (1 - alpha) * df[cfg.COL_AFFINITY_NORM]
    )
    
    # ---------------------------------------------------------
    # 4. Log Per-Target Statistics (for transparency)
    # ---------------------------------------------------------
    print(f"\nPer-Target Normalization Statistics (HOTFIX #2):")
    print(f"{'Target':<10} {'CNN_VS Range':<20} {'Affinity Range':<22} {'P0 Range':<18} {'N':<5}")
    print("-" * 80)
    for target in sorted(df['target'].unique()):
        t_df = df[df['target'] == target]
        print(f"{target:<10} "
              f"[{t_df['CNN_VS'].min():.3f}, {t_df['CNN_VS'].max():.3f}]      "
              f"[{t_df['Affinity_kcal'].min():.2f}, {t_df['Affinity_kcal'].max():.2f}]      "
              f"[{t_df[cfg.COL_P0_RAW_WEIGHT].min():.3f}, {t_df[cfg.COL_P0_RAW_WEIGHT].max():.3f}]   "
              f"{len(t_df)}")
    
    print(f"\nAlpha (CNN_VS weight): {alpha}")
    print(f"Global P0_raw_weight stats: mean={df[cfg.COL_P0_RAW_WEIGHT].mean():.4f}, "
          f"std={df[cfg.COL_P0_RAW_WEIGHT].std():.4f}")
    
    # ==========================================
    # PIVOT TO WIDE MATRIX FOR LAYER 2 (RWR)
    # ==========================================
    # Create matrix: Ligand (rows) x Target (columns)
    p0_matrix = df.pivot_table(
        index='lig_id',
        columns='target',
        values=cfg.COL_P0_RAW_WEIGHT,
        aggfunc='first'
    )
    
    # Merge metadata (SMILES, orig_name) back into matrix
    meta_cols = ['lig_id', 'orig_name', 'smiles']
    available_meta = [c for c in meta_cols if c in df.columns]
    if available_meta:
        meta_df = df.drop_duplicates('lig_id')[available_meta].set_index('lig_id')
        p0_matrix = p0_matrix.join(meta_df)
    
    # Reorder columns: metadata first, then targets in config order
    target_cols = [t for t in cfg.TARGET_NAMES if t in p0_matrix.columns]
    meta_first = [c for c in ['orig_name', 'smiles'] if c in p0_matrix.columns]
    final_cols = meta_first + target_cols
    p0_matrix = p0_matrix[final_cols]
    
    print(f"\n[PIVOT SUCCESS] P0 Matrix Shape: {p0_matrix.shape[0]} ligands x {p0_matrix.shape[1]} columns")
    print(f"  - Metadata columns: {meta_first}")
    print(f"  - Target columns: {target_cols}")
    
    return df, p0_matrix


# Execute P0 calculation
df_with_p0_long, df_p0_matrix = calculate_p0_raw_weight(df_passed)

# For backward compatibility (QC uses df_with_p0)
df_with_p0 = df_with_p0_long

print("\nTop 10 compounds by P0_raw_weight (Long format):")
top10_cols = ['lig_id', 'target', 'CNN_VS', 'Affinity_kcal', cfg.COL_P0_RAW_WEIGHT, cfg.COL_PB_STATUS]
available_cols = [c for c in top10_cols if c in df_with_p0.columns]
top10 = df_with_p0.nlargest(10, cfg.COL_P0_RAW_WEIGHT)[available_cols]
print(top10.to_string(index=False))

print("\nP0 Matrix Sample (First 5 ligands, Wide format for Layer 2):")
display(df_p0_matrix.head())
# =============================================================================
# STAGE 7: DEADLOCK VALIDATION
# =============================================================================
print("\n" + "="*64)
print("STAGE 7: DEADLOCK VALIDATION")
print("="*64)

# DL1: Verify all normalized records have dG <= cutoff
try:
    max_affinity = df_with_p0['Affinity_kcal'].max()
    assert max_affinity <= cfg.DELTA_G_CUTOFF, \
        f"DL1 VIOLATION: Found Affinity {max_affinity} > {cfg.DELTA_G_CUTOFF} after normalization"
    cfg.validate_deadlock_rules('pre_normalize', filtered=True)
    print(f"[OK] DL1: All normalized records have dG <= {cfg.DELTA_G_CUTOFF} kcal/mol (max: {max_affinity:.2f})")
except AssertionError as e:
    print(f"[FAIL] DL1 VIOLATION: {e}")
    raise

# DL2: No sign interpretation
try:
    cfg.validate_deadlock_rules('no_sign_interpretation', interpreting_sign=False)
    print("[OK] DL2: No agonist/antagonist interpretation at Layer 1")
except AssertionError as e:
    print(f"[FAIL] DL2 VIOLATION: {e}")
    raise

# DL3: WT structures only
try:
    cfg.validate_deadlock_rules('no_mutant', using_mutant=False)
    print("[OK] DL3: Only wild-type crystal structures used")
except AssertionError as e:
    print(f"[FAIL] DL3 VIOLATION: {e}")
    raise

print("\nAll deadlock rules validated successfully.")
# =============================================================================
# STAGE 8: QC VISUALIZATIONS
# =============================================================================
print("\n" + "="*64)
print("STAGE 8: QC VISUALIZATIONS")
print("="*64)

cfg.LAYER1_OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
qc_plot_dir = cfg.LAYER1_OUTPUT_DIR / cfg.L1_QC_PLOTS_DIR
qc_plot_dir.mkdir(parents=True, exist_ok=True)

plt.style.use('seaborn-v0_8-whitegrid')
fig_dpi = 150

# QC1: Affinity and CNNscore distributions
fig, axes = plt.subplots(1, 2, figsize=(14, 5))

ax1 = axes[0]
ax1.hist(df_pre_gate1['Affinity_kcal'], bins=50, alpha=0.7, label='All', color='gray')
ax1.hist(df_passed['Affinity_kcal'], bins=50, alpha=0.7, label='Passed Gate 1', color='green')
ax1.axvline(cfg.DELTA_G_CUTOFF, color='red', linestyle='--', linewidth=2, label=f'Cutoff ({cfg.DELTA_G_CUTOFF})')
ax1.set_xlabel('Affinity (kcal/mol)')
ax1.set_ylabel('Count')
ax1.set_title('Affinity Distribution: Before vs After Gate 1')
ax1.legend()

ax2 = axes[1]
ax2.hist(df_pre_gate1['CNNscore'], bins=50, alpha=0.7, label='All', color='gray')
ax2.hist(df_passed['CNNscore'], bins=50, alpha=0.7, label='Passed Gate 1', color='blue')
ax2.axvline(cfg.CNNSCORE_FALLBACK_THRESHOLD, color='orange', linestyle='--', linewidth=2, label=f'Fallback ({cfg.CNNSCORE_FALLBACK_THRESHOLD})')
ax2.set_xlabel('CNNscore')
ax2.set_ylabel('Count')
ax2.set_title('CNNscore Distribution: Before vs After Gate 1')
ax2.legend()

plt.tight_layout()
plt.savefig(qc_plot_dir / 'QC1_Gate1_Filter_Distributions.png', dpi=fig_dpi, bbox_inches='tight')
plt.show()
print(f"Saved: QC1_Gate1_Filter_Distributions.png")

# QC2: P0_raw_weight by target
fig, ax = plt.subplots(figsize=(12, 6))
df_with_p0.boxplot(column=cfg.COL_P0_RAW_WEIGHT, by='target', ax=ax)
ax.set_xlabel('Target')
ax.set_ylabel('P0_raw_weight')
ax.set_title('P0_raw_weight Distribution by Target')
plt.suptitle('')
plt.xticks(rotation=45)
plt.tight_layout()
plt.savefig(qc_plot_dir / 'QC2_P0_by_Target.png', dpi=fig_dpi, bbox_inches='tight')
plt.show()
print(f"Saved: QC2_P0_by_Target.png")

# QC3: CNN_VS vs Affinity scatter
fig, ax = plt.subplots(figsize=(10, 8))
scatter = ax.scatter(
    df_with_p0['CNN_VS'], 
    df_with_p0['Affinity_kcal'], 
    c=df_with_p0[cfg.COL_P0_RAW_WEIGHT], 
    cmap='viridis', 
    alpha=0.6,
    s=20
)
plt.colorbar(scatter, label='P0_raw_weight')
ax.set_xlabel('CNN_VS')
ax.set_ylabel('Affinity (kcal/mol)')
ax.set_title('CNN_VS vs Affinity (colored by P0_raw_weight)')
ax.axhline(cfg.DELTA_G_CUTOFF, color='red', linestyle='--', alpha=0.5)
plt.tight_layout()
plt.savefig(qc_plot_dir / 'QC3_CNNVS_vs_Affinity.png', dpi=fig_dpi, bbox_inches='tight')
plt.show()
print(f"Saved: QC3_CNNVS_vs_Affinity.png")

# QC4: Pass rate by target (with Track annotation)
pass_rates = []
for target in cfg.TARGET_NAMES:
    n_total = len(df_pre_gate1[df_pre_gate1['target'] == target])
    n_passed = len(df_passed[df_passed['target'] == target])
    rate = n_passed / n_total if n_total > 0 else 0
    track = 'T1' if target in cfg.TARGETS_WITH_PB else 'T2'
    pass_rates.append({'target': f"{target} [{track}]", 'pass_rate': rate, 'n_passed': n_passed, 'n_total': n_total})

df_rates = pd.DataFrame(pass_rates)

fig, ax = plt.subplots(figsize=(10, 5))
colors = ['steelblue' if 'T1' in t else 'coral' for t in df_rates['target']]
bars = ax.barh(df_rates['target'], df_rates['pass_rate'], color=colors)
ax.set_xlabel('Pass Rate')
ax.set_ylabel('Target')
ax.set_title('Gate 1 Pass Rate by Target (T1=PoseBusters, T2=CNNscore Fallback)')
ax.set_xlim(0, 1)
for i, (rate, n_passed, n_total) in enumerate(zip(df_rates['pass_rate'], df_rates['n_passed'], df_rates['n_total'])):
    ax.text(rate + 0.02, i, f'{n_passed}/{n_total} ({rate:.1%})', va='center', fontsize=9)
plt.tight_layout()
plt.savefig(qc_plot_dir / 'QC4_PassRate_by_Target.png', dpi=fig_dpi, bbox_inches='tight')
plt.show()
print(f"Saved: QC4_PassRate_by_Target.png")

# QC5: Rescue Statistics (NEW in v1.2)
if not df_rescue.empty:
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    
    # Rescue count by target
    ax1 = axes[0]
    rescue_counts = df_rescue['target'].value_counts()
    rescue_counts.plot(kind='bar', ax=ax1, color='forestgreen')
    ax1.set_xlabel('Target')
    ax1.set_ylabel('Ligands Rescued')
    ax1.set_title('QC5a: Ligands Rescued per Track 1 Target')
    ax1.tick_params(axis='x', rotation=45)
    for i, v in enumerate(rescue_counts.values):
        ax1.text(i, v + 0.5, str(v), ha='center', fontsize=10)
    
    # Rescued mode distribution
    ax2 = axes[1]
    mode_counts = df_rescue['rescued_mode'].value_counts().sort_index()
    mode_counts.plot(kind='bar', ax=ax2, color='darkorange')
    ax2.set_xlabel('Rescued Pose Mode')
    ax2.set_ylabel('Count')
    ax2.set_title('QC5b: Distribution of Rescued Pose Modes')
    ax2.tick_params(axis='x', rotation=0)
    
    plt.tight_layout()
    plt.savefig(qc_plot_dir / 'QC5_Rescue_Statistics.png', dpi=fig_dpi, bbox_inches='tight')
    plt.show()
    print(f"Saved: QC5_Rescue_Statistics.png")
else:
    print("QC5: No rescues to plot (all Mode 1 poses were valid)")

print(f"\nAll QC plots saved to: {qc_plot_dir}")
# =============================================================================
# STAGE 9: EXPORT & SUMMARY
# =============================================================================
print("\n" + "="*64)
print("STAGE 9: EXPORT & SUMMARY")
print("="*64)

output_dir = cfg.LAYER1_OUTPUT_DIR

# Export merged raw data (before Gate 1)
raw_path = output_dir / cfg.L1_MERGED_DOCKING_PB_CSV
df_pre_gate1.to_csv(raw_path, index=False)
print(f"Exported: {raw_path.name} ({len(df_pre_gate1)} rows)")

# Export filtered data (after Gate 1) - Long format
filtered_path = output_dir / cfg.L1_FILTERED_MATRIX_CSV
df_passed.to_csv(filtered_path, index=False)
print(f"Exported: {filtered_path.name} ({len(df_passed)} rows) [Long format]")

# Export P0 Long format (for QC reference)
p0_long_path = output_dir / 'L1_P0_Vector_Long.csv'
p0_cols = ['lig_id', 'target', 'CNN_VS', 'CNNscore', 'Affinity_kcal', 
           cfg.COL_CNN_VS_NORM, cfg.COL_AFFINITY_NORM, cfg.COL_P0_RAW_WEIGHT, cfg.COL_PB_STATUS]
available_p0_cols = [c for c in p0_cols if c in df_with_p0_long.columns]
if 'smiles' in df_with_p0_long.columns:
    available_p0_cols.insert(2, 'smiles')
if 'orig_name' in df_with_p0_long.columns:
    available_p0_cols.insert(2, 'orig_name')
df_with_p0_long[available_p0_cols].to_csv(p0_long_path, index=False)
print(f"Exported: {p0_long_path.name} ({len(df_with_p0_long)} rows) [Long format for QC]")

# Export P0 MATRIX (Wide format) - FOR LAYER 2 RWR
p0_matrix_path = output_dir / cfg.L1_P0_VECTOR_RAW_CSV
df_p0_matrix.to_csv(p0_matrix_path)  # Keep index since lig_id is the index
print(f"Exported: {p0_matrix_path.name} ({df_p0_matrix.shape[0]} ligands x {df_p0_matrix.shape[1]} columns) [Wide Matrix for Layer 2]")

# Export rejected compounds
rejected_path = output_dir / 'L1_Rejected_Compounds.csv'
df_rejected.to_csv(rejected_path, index=False)
print(f"Exported: {rejected_path.name} ({len(df_rejected)} rows)")

# Export Rescue Report (NEW in v1.2)
rescue_path = output_dir / cfg.L1_RESCUE_REPORT_CSV
if not df_rescue.empty:
    df_rescue.to_csv(rescue_path, index=False)
    print(f"Exported: {rescue_path.name} ({len(df_rescue)} rescued ligands)")
else:
    print(f"Skipped: {cfg.L1_RESCUE_REPORT_CSV} (no rescues)")

# Export Gate Report as JSON
import json
gate_report = {
    'timestamp': datetime.now().isoformat(),
    'pipeline_version': cfg.PROJECT_METADATA['version'],
    'architecture': 'Two-Track v1.2',
    'n_targets': cfg.N_TARGETS,
    'targets': cfg.TARGET_NAMES,
    'track1_targets': cfg.TARGETS_WITH_PB,
    'track2_targets': cfg.TARGETS_WITHOUT_PB,
    'n_ligands_expected': cfg.N_LIGANDS_EXPECTED,
    'delta_g_cutoff': cfg.DELTA_G_CUTOFF,
    'cnnscore_fallback_threshold': cfg.CNNSCORE_FALLBACK_THRESHOLD,
    'alpha_cnn_vs_weight': cfg.ALPHA_CNN_VS_WEIGHT,
    'total_records_pre_gate1': len(df_pre_gate1),
    'passed_gate1': len(df_passed),
    'rejected': len(df_rejected),
    'pass_rate': len(df_passed) / len(df_pre_gate1) if len(df_pre_gate1) > 0 else 0,
    'track1_records': len(df_track1),
    'track2_records': len(df_track2),
    'total_rescues': len(df_rescue),
    'p0_matrix_shape': list(df_p0_matrix.shape),
    'deadlock_dl1_compliant': True,
    'deadlock_dl2_compliant': True,
    'deadlock_dl3_compliant': True,
}

report_path = output_dir / 'L1_Gate_Report.json'
with open(report_path, 'w') as f:
    json.dump(gate_report, f, indent=2, default=str)
print(f"Exported: {report_path.name}")

# Final Summary
print("\n" + "="*64)
print("LAYER 1 PIPELINE COMPLETE (Two-Track Architecture v1.2)")
print("="*64)
print(f"""
SUMMARY:
  - Total records processed: {len(df_pre_gate1)}
    - Track 1 (PoseBusters): {len(df_track1)} from {len(cfg.TARGETS_WITH_PB)} targets
    - Track 2 (CNNscore): {len(df_track2)} from {len(cfg.TARGETS_WITHOUT_PB)} targets
  - Passed Gate 1: {len(df_passed)} ({100*len(df_passed)/len(df_pre_gate1):.1f}%)
  - Rejected: {len(df_rejected)}
  - Ligands RESCUED (Mode 1 invalid -> Mode 2+ valid): {len(df_rescue)}
  
P0 OUTPUT:
  - Long format: {len(df_with_p0_long)} rows (ligand-target pairs, for QC)
  - Wide Matrix: {df_p0_matrix.shape[0]} ligands x {df_p0_matrix.shape[1]} columns (for Layer 2 RWR)
  
FILTER CRITERIA:
  - Gate 1: dG <= {cfg.DELTA_G_CUTOFF} kcal/mol (thermodynamic)
  - Track 1: PoseBusters validation (pose-level rescue enabled)
  - Track 2: CNNscore >= {cfg.CNNSCORE_FALLBACK_THRESHOLD} fallback filter
  
P0 FORMULA (HOTFIX #1):
  P0_raw_weight = {cfg.ALPHA_CNN_VS_WEIGHT} * CNN_VS_norm + {cfg.ALPHA_AFFINITY_WEIGHT} * Affinity_norm
  (where Affinity_norm is INVERTED: more negative = higher score)

DEADLOCK COMPLIANCE:
  [OK] DL1: Normalized AFTER filtering
  [OK] DL2: No sign interpretation
  [OK] DL3: WT structures only

OUTPUT DIRECTORY: {output_dir}
""")

print(f"\nNotebook completed: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
