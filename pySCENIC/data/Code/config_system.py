# ============================================================
# config_system.py - Central Configuration Hub v1.2
# ============================================================
# Purpose: Single source of truth for multi-layer analysis pipeline
#          (Layer 1 -> Layer N). All notebooks import from here
#          to ensure consistency (FAIR: Interoperable).
#
# Design principles:
#   - SMART: Each constant has specific, measurable meaning
#   - FAIR:  Full metadata, reusable
#   - No hardcoded values in notebooks
#
# History:
#   v1.0 - Initial for Layer 1 (Thermodynamic Gate)
#   v1.1 - Added LAYER1_OUTPUT_DIR, PoseBusters paths, logging config,
#          CNNscore fallback threshold, P0_raw_weight constants
#   v1.2 - Two-Track Architecture: Added pose-level columns, sheet indices,
#          track definitions (TARGETS_WITH_PB, TARGETS_WITHOUT_PB),
#          rescue report path, fixed COL_PB_VALID column name
# ============================================================

from pathlib import Path
import os
import logging

# ============================================================
# SYSTEM PATHS
# ============================================================
PROJECT_ROOT = Path(__file__).resolve().parent
DOCKING_PARENT_DIR = Path(r"D:\khoa_luan\(gnina) docking results")

# ============================================================
# LAYER 1 OUTPUT DIRECTORIES (NEW v1.1)
# ============================================================
LAYER1_OUTPUT_DIR = PROJECT_ROOT / "outputs" / "Layer1_Thermodynamic_Gate"

# ============================================================
# POSEBUSTERS VALIDATION PATHS (NEW v1.1)
# ============================================================
POSEBUSTERS_PARENT_DIR = DOCKING_PARENT_DIR / "PoseBuster_results"

# Mapping: target_name -> path to posebusters CSV (None if not available)
POSEBUSTERS_PATHS = {
    "PPARA": POSEBUSTERS_PARENT_DIR / "ppara_7bq2" / "summary" / "posebusters_master_results.csv",
    "PPARD": POSEBUSTERS_PARENT_DIR / "ppard_7wgn" / "summary" / "posebusters_master_results.csv",
    "PPARG": None,  # PB_NOT_AVAILABLE
    "EGFR": None,   # PB_NOT_AVAILABLE
    "ERBB2": None,  # PB_NOT_AVAILABLE
    "KDR": None,    # PB_NOT_AVAILABLE
    "PTGS2": None,  # PB_NOT_AVAILABLE
    "ALOX-5": POSEBUSTERS_PARENT_DIR / "alox5_6n2w" / "summary" / "posebusters_master_results.csv",
    "PTGES": None,  # PB_NOT_AVAILABLE
}

# ============================================================
# TWO-TRACK ARCHITECTURE (NEW v1.2)
# ============================================================
# Track 1: Targets with PoseBusters validation (pose-level rescue logic)
# Track 2: Targets without PoseBusters (CNNscore fallback at ligand level)
TARGETS_WITH_PB = ["PPARA", "PPARD", "ALOX-5"]
TARGETS_WITHOUT_PB = ["PPARG", "EGFR", "ERBB2", "KDR", "PTGS2", "PTGES"]

# Excel sheet indices
SHEET_INTRA_LIGAND = 0  # Intra-Ligand_Poses (wide format, P1_*, P2_*, P3_*)
SHEET_INTER_LIGAND = 1  # Inter-Ligand_Ranking (179 rows, best pose per ligand)

# ============================================================
# 9 TARGET PROTEINS
# ============================================================
# 9 targets: PPAR family, ErbB family, COX-2
# Each has WT PDB ID (wild-type crystal structure)
#
# DEADLOCK #3: Only use WT crystal structures.
#    DO NOT simulate mutant proteins at Layer 1.
# ============================================================
TARGET_PROTEINS = {
    # -- PPAR family (Peroxisome Proliferator-Activated Receptors) --
    "PPARA": {
        "pdb_id": "7BQ2",
        "family": "PPAR",
        "flex_residues": "A:273",
        "docking_dir": DOCKING_PARENT_DIR / "docking_results_ppara_7bq2",
        "description": "PPARa",
    },
    "PPARD": {
        "pdb_id": "7WGN",
        "family": "PPAR",
        "flex_residues": "A:293 A:303 A:312 A:417",
        "docking_dir": DOCKING_PARENT_DIR / "docking_results_ppard_7wgn",
        "description": "PPARd",
    },
    "PPARG": {
        "pdb_id": "9F7W",
        "family": "PPAR",
        "flex_residues": "A:259 A:262 A:263 A:283",
        "docking_dir": DOCKING_PARENT_DIR / "docking_results_pparg_9f7w",
        "description": "PPARg",
    },
    # -- ErbB family (Epidermal growth factor receptor family) --
    "EGFR": {
        "pdb_id": "1XKK",
        "family": "ErbB",
        "flex_residues": "A:745 A:776 A:777 A:790 A:856",
        "docking_dir": DOCKING_PARENT_DIR / "docking_results_egfr_1xkk",
        "description": "Epidermal growth factor receptor",
    },
    "ERBB2": {
        "pdb_id": "7PCD",
        "family": "ErbB",
        "flex_residues": "A:774 A:783 A:862 A:864",
        "docking_dir": DOCKING_PARENT_DIR / "docking_results_erbb2_7pcd",
        "description": "Receptor tyrosine-protein kinase erbB-2",
    },
    "KDR": {
        "pdb_id": "5EW3",
        "family": "Tyr protein kinase family",
        "flex_residues": "A:885 A:916 A:1026",
        "docking_dir": DOCKING_PARENT_DIR / "docking_results_kdr_5ew3",
        "description": "Kinase Insert Domain Receptor - VEGF",
    },
    # -- COX-2 --
    "PTGS2": {
        "pdb_id": "5KIR",
        "family": "COX",
        "flex_residues": "A:90 A:385 A:434 A:513 A:523",
        "docking_dir": DOCKING_PARENT_DIR / "docking_results_ptgs2_5kir",
        "description": "COX-2 - Cyclooxygenase-2",
    },
    # -- ALOX-5 --
    "ALOX-5": {
        "pdb_id": "6N2W",
        "family": "LOX",
        "flex_residues": "B:359 B:363 B:368 B:596 B:599 B:607",
        "docking_dir": DOCKING_PARENT_DIR / "docking_results_alox5_6n2w",
        "description": "Arachidonate lipoxygenases",
    },
    # -- PTGES --
    "PTGES": {
        "pdb_id": "5TL9",
        "family": "MAPEG",
        "flex_residues": "B:28,B:32,B:52,B:53,A:130,A:134",
        "docking_dir": DOCKING_PARENT_DIR / "docking_results_ptges_5tl9",
        "description": "Prostaglandin synthases",
    },
}

# Target names list (for iteration)
TARGET_NAMES = list(TARGET_PROTEINS.keys())
N_TARGETS = len(TARGET_NAMES)

# ============================================================
# LIGAND COUNT
# ============================================================
N_LIGANDS_EXPECTED = 179  # 179 benzyl ether derivatives

# ============================================================
# THERMODYNAMIC THRESHOLD - LAYER 1 (Thermodynamic Gate)
# ============================================================
# Gibbs free energy (dG) threshold for interaction filtering
#
# Scientific basis:
#   dG >= -6.5 kcal/mol -> binding too weak, no biological significance
#   at lab conditions (25C, 1 atm).
#   Equivalent to Kd ~ 17 uM - outside typical drug-like range (Kd < 10 uM).
#
# DEADLOCK #1 (Normalization Paradox):
#   ABSOLUTELY DO NOT normalize P0 vector if compounds with
#   dG > -6.5 kcal/mol have not been filtered out.
#   Filter FIRST, normalize LATER.
# ============================================================
DELTA_G_CUTOFF = -6.5  # kcal/mol - Gate 1 threshold

# Positive affinity = repulsive = physically meaningless
DELTA_G_REPULSIVE = 0.0  # kcal/mol

# ============================================================
# CNNSCORE FALLBACK THRESHOLD (NEW v1.1)
# ============================================================
# Fallback CNNscore threshold for targets missing PoseBusters validation
#
# Scientific basis:
#   CNNscore is probability that pose reflects natural binding state
#   CNNscore >= 0.5 -> likelihood RMSD < 2A is positive
#   Used as proxy when PoseBusters unavailable for geometry validation
#
# Application logic:
#   - PB_VALID: Pass through (geometry validated)
#   - PB_NOT_AVAILABLE: Require CNNscore >= 0.5
#   - PB_INVALID: Reject completely
# ============================================================
CNNSCORE_FALLBACK_THRESHOLD = 0.5

# ============================================================
# WEIGHT MIXING - CNN_VS vs AFFINITY (Layer 1)
# ============================================================
# Alpha weight for linear combination function:
#   P0_raw_weight = alpha * CNN_VS_norm + (1 - alpha) * Affinity_norm
#
# Basis:
#   alpha = 0.7 -> prioritize CNN_VS (statistical ranking from
#   neural network) over Affinity (force field energy).
#   CNN_VS proven superior in virtual screening by GNINA team
#   (CACHE Challenge #1, first place).
#
# Normalization formulas:
#   CNN_VS_norm = (CNN_VS - min) / (max - min)        # Higher is better
#   Affinity_norm = (max - Affinity) / (max - min)    # More negative is better (inverted)
# ============================================================
ALPHA_CNN_VS_WEIGHT = 0.7
ALPHA_AFFINITY_WEIGHT = 1.0 - ALPHA_CNN_VS_WEIGHT  # 0.3

# ============================================================
# GNINA DATA COLUMNS (Column Mappings)
# ============================================================
# Column names from GNINA v2.6.2 SDF/Excel output
COL_MINIMIZED_AFFINITY = "minimizedAffinity"
COL_CNN_SCORE = "CNNscore"
COL_CNN_AFFINITY = "CNNaffinity"
COL_CNN_VS = "CNN_VS"
COL_VINARDO = "vinardo"

# Metadata columns
COL_LIGAND_ID = "lig_id"
COL_LIGAND_NAME = "orig_name"
COL_SMILES = "smiles"
COL_TARGET = "target"

# Computed columns (NEW v1.1)
COL_CNN_VS_NORM = "CNN_VS_norm"
COL_AFFINITY_NORM = "Affinity_norm"
COL_P0_RAW_WEIGHT = "P0_raw_weight"

# ============================================================
# POSEBUSTERS COLUMN CONSTANTS (UPDATED v1.2)
# ============================================================
# CRITICAL FIX: Column name is "PoseBuster_Valid" NOT "all_valid"
COL_PB_VALID = "PoseBuster_Valid"    # Boolean from posebusters_master_results.csv
COL_PB_STATUS = "PB_Status"          # New column added to merged DataFrame

# Pose-level columns from PoseBusters CSV (NEW v1.2)
COL_POSE_ID = "pose_id"              # e.g., "pose_000" to "pose_009"
COL_GNINA_MODEL = "original_gnina_model"  # 1-10

# Enum values for PB_Status
PB_STATUS_VALID = "PB_VALID"
PB_STATUS_INVALID = "PB_INVALID"
PB_STATUS_NOT_AVAILABLE = "PB_NOT_AVAILABLE"

# ============================================================
# FLAG LABELS - synchronized with GNINA v2.6.2
# ============================================================
FLAG_POSITIVE_AFFINITY = "POSITIVE_AFFINITY"
FLAG_POOR_AFFINITY = "POOR_AFFINITY"
FLAG_LOW_CNNSCORE = "LOW_CNNSCORE"  # NEW v1.1: CNNscore < 0.5 when PB missing
FLAG_PB_INVALID = "PB_INVALID"      # NEW v1.1: PoseBusters failed
FLAG_CLEAN = "CLEAN"

# Flag column name
COL_FLAG = "Sanity_Flag"

# ============================================================
# NORMALIZATION PARAMETERS
# ============================================================
NORMALIZATION_METHOD = "min_max"  # "min_max" | "z_score" | "robust"

# ============================================================
# OUTPUT FILENAMES
# ============================================================
# Layer 1 outputs
L1_RAW_MATRIX_CSV = "L1_Raw_Docking_Matrix.csv"
L1_FILTERED_MATRIX_CSV = "L1_Filtered_Matrix_PostGate.csv"
L1_P0_RAW_WEIGHTS_CSV = "L1_P0_Raw_Weights.csv"
L1_GATE_REPORT_CSV = "L1_Gate_Report.csv"
L1_QC_PLOTS_DIR = "L1_QC_Plots"

# Extended output filenames (NEW v1.1)
L1_MERGED_DOCKING_PB_CSV = "L1_Merged_Docking_PoseBusters.csv"
L1_P0_VECTOR_RAW_CSV = "L1_P0_Vector_Raw.csv"           # BEFORE normalization
L1_P0_VECTOR_FILTERED_CSV = "L1_P0_Vector_Filtered.csv"  # AFTER Gate 1
L1_P0_VECTOR_NORMALIZED_CSV = "L1_P0_Vector_Normalized.csv"  # AFTER normalization

# Rescue report (NEW v1.2)
L1_RESCUE_REPORT_CSV = "L1_Rescue_Report.csv"

# ============================================================
# LOGGING CONFIGURATION (NEW v1.1)
# ============================================================
LOG_FORMAT = "%(asctime)s | %(levelname)-8s | %(name)s | %(message)s"
LOG_DATEFMT = "%Y-%m-%d %H:%M:%S"
LOG_LEVEL = logging.INFO


def setup_logger(name: str = "Layer1") -> logging.Logger:
    """
    Create standard logger for notebook/module.
    
    SMART:
      S - Create logger with standard format
      M - Log level is measurable (INFO/DEBUG/WARNING/ERROR)
      A - Uses Python standard logging module
      R - Supports debugging and audit trail
      T - Instant setup
    """
    logger = logging.getLogger(name)
    if not logger.handlers:
        handler = logging.StreamHandler()
        handler.setFormatter(logging.Formatter(LOG_FORMAT, LOG_DATEFMT))
        logger.addHandler(handler)
    logger.setLevel(LOG_LEVEL)
    return logger


# ============================================================
# CNN_VS CLARIFICATION (NEW v1.1)
# ============================================================
CNN_VS_NOTE = (
    "CNN_VS is a RELATIVE RANKING metric - NO absolute threshold applies. "
    "The previous 0.70 threshold was WRONG. Only use ranking to select top-N."
)

# ============================================================
# PROJECT METADATA (FAIR Compliant)
# ============================================================
PROJECT_METADATA = {
    "project_name": "Benzyl Ether Multi-Target Docking Campaign",
    "version": "1.2.0",
    "pipeline_version": "GNINA_v2.6.2",
    "n_ligands": N_LIGANDS_EXPECTED,
    "n_targets": N_TARGETS,
    "targets": TARGET_NAMES,
    "docking_engine": "GNINA (CNN-based scoring)",
    "scoring_functions": [
        "minimizedAffinity (Vina forcefield)",
        "CNNscore (pose probability)",
        "CNNaffinity (CNN-predicted affinity)",
        "CNN_VS (virtual screening rank score)",
    ],
    "flex_docking": True,
    "exhaustiveness": 64,
    "num_modes": 10,
    "cnn_scoring": "rescore",
    "pose_sort_order": "CNNscore",
    "seed": 42,
    "author": None,
    "institution": None,
    "license": "CC-BY-4.0",
    "fair_principles": {
        "findable": "Persistent ligand IDs (LIG_XXXX), structured directories",
        "accessible": "Local filesystem, exportable CSV/XLSX",
        "interoperable": "Standard SDF/PDB/CSV formats, RDKit compatible",
        "reusable": "Full provenance in STATUS.txt, META.txt, command.txt",
    },
}

# ============================================================
# DEADLOCK RULES - CHECKLIST
# ============================================================
DEADLOCK_RULES = {
    "DL1_NORMALIZE_BEFORE_FILTER": (
        "ABSOLUTELY DO NOT normalize P0 vector if compounds with "
        f"dG > {DELTA_G_CUTOFF} kcal/mol have not been filtered out. "
        "Filter FIRST, normalize LATER."
    ),
    "DL2_SIGN_INTERPRETATION": (
        "ABSOLUTELY DO NOT decode agonist/antagonist (+/-) "
        "at Layer 1. Kinetic signs belong to Layer 2+."
    ),
    "DL3_MUTANT_SIMULATION": (
        "ABSOLUTELY DO NOT simulate mutant proteins at Layer 1. "
        "Only use wild-type (WT) crystal structures."
    ),
}

# ============================================================
# UTILITY FUNCTIONS
# ============================================================
def print_config_summary():
    """Print current configuration summary - call at start of each notebook."""
    print("=" * 64)
    print("CONFIG SYSTEM - Central Configuration Hub v1.2")
    print("=" * 64)
    print(f"Project root:       {PROJECT_ROOT}")
    print(f"Docking data:       {DOCKING_PARENT_DIR}")
    print(f"Layer 1 output:     {LAYER1_OUTPUT_DIR}")
    print(f"Targets:            {N_TARGETS} ({', '.join(TARGET_NAMES)})")
    print(f"  Track 1 (PB):     {', '.join(TARGETS_WITH_PB)}")
    print(f"  Track 2 (no PB):  {', '.join(TARGETS_WITHOUT_PB)}")
    print(f"Ligands expected:   {N_LIGANDS_EXPECTED}")
    print(f"dG cutoff:          {DELTA_G_CUTOFF} kcal/mol")
    print(f"CNNscore fallback:  {CNNSCORE_FALLBACK_THRESHOLD}")
    print(f"alpha (CNN_VS):     {ALPHA_CNN_VS_WEIGHT}")
    print(f"Normalization:      {NORMALIZATION_METHOD}")
    print(f"Deadlock rules:     {len(DEADLOCK_RULES)} active")
    print("=" * 64)


def validate_deadlock_rules(step: str, **context):
    """
    Validate deadlock rules at each step.
    
    SMART:
      S - Check specific rule by step name
      M - Return True/False or raise AssertionError
      A - Simple assertion-based validation
      R - Critical for pipeline integrity
      T - O(1) per check
    """
    if step == "pre_normalize":
        assert context.get("filtered", False), (
            f"DEADLOCK VIOLATION - DL1:\n"
            f"{DEADLOCK_RULES['DL1_NORMALIZE_BEFORE_FILTER']}"
        )
    elif step == "no_sign_interpretation":
        assert not context.get("interpreting_sign", False), (
            f"DEADLOCK VIOLATION - DL2:\n"
            f"{DEADLOCK_RULES['DL2_SIGN_INTERPRETATION']}"
        )
    elif step == "no_mutant":
        assert not context.get("using_mutant", False), (
            f"DEADLOCK VIOLATION - DL3:\n"
            f"{DEADLOCK_RULES['DL3_MUTANT_SIMULATION']}"
        )


# ============================================================
# Self-test on import
# ============================================================
if __name__ == "__main__":
    print_config_summary()
    print("\nconfig_system.py v1.2 loaded successfully.")
    print(f"   Deadlock rules: {list(DEADLOCK_RULES.keys())}")
    pb_configured = sum(1 for v in POSEBUSTERS_PATHS.values() if v is not None)
    print(f"   PoseBusters paths configured: {pb_configured}/{len(POSEBUSTERS_PATHS)}")
    print(f"   Track 1 targets: {TARGETS_WITH_PB}")
    print(f"   Track 2 targets: {TARGETS_WITHOUT_PB}")
