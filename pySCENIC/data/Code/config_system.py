# ============================================================
# config_system.py - Central Configuration Hub v1.4
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
#   v1.3 - Layer 2A: SCENIC GRN Inference configuration
#          - CCLE/DepMap data paths
#          - pySCENIC resource paths (TF list, motif DBs)
#          - Lineage definitions (OncotreeLineage level)
#          - GRNBoost2 parameters
#   v1.4 - Layer 2B: Heterogeneous Topology RWR configuration
#          - STRING PPI data paths
#          - RWR parameters (alpha, pseudo_count)
#          - Target cell line configuration (MCF-7)
#          - Active/Inactive percentile thresholds
#          - Delta-Network analysis parameters
#          - normalize_gene_name() utility function
#          - Deadlock rules DL6, DL7, DL8
# ============================================================

from pathlib import Path
import os
import logging

# ============================================================
# SYSTEM PATHS
# ============================================================
PROJECT_ROOT = Path(__file__).resolve().parent
DOCKING_PARENT_DIR = Path(r"/home/labhhc5/Documents/workspace/D21/Duong Huy/gnina_project")

# ============================================================
# LAYER 1 OUTPUT DIRECTORIES (NEW v1.1)
# ============================================================
LAYER1_OUTPUT_DIR = PROJECT_ROOT / "outputs" / "Layer1_Thermodynamic_Gate"

# ============================================================
# LAYER 2A OUTPUT DIRECTORIES (NEW v1.3)
# ============================================================
LAYER2A_OUTPUT_DIR = PROJECT_ROOT / "outputs" / "Layer2A_SCENIC_GRN"

# ============================================================
# POSEBUSTERS VALIDATION PATHS (NEW v1.1)
# ============================================================
POSEBUSTERS_PARENT_DIR = Path(r'/home/labhhc5/Documents/workspace/D21/Duong Huy/posebusters')

# Mapping: target_name -> path to posebusters CSV (None if not available)
POSEBUSTERS_PATHS = {
    "PPARA": POSEBUSTERS_PARENT_DIR / "ppara_7bq2" / "summary" / "posebusters_master_results.csv",
    "PPARD": POSEBUSTERS_PARENT_DIR / "ppard_7wgn" / "summary" / "posebusters_master_results.csv",
    "PPARG": POSEBUSTERS_PARENT_DIR / "pparg_9f7w" / "summary" / "posebusters_master_results.csv",  # PB_NOT_AVAILABLE -> None
    "EGFR": POSEBUSTERS_PARENT_DIR / "egfr_1xkk" / "summary" / "posebusters_master_results.csv",   # PB_NOT_AVAILABLE
    "ERBB2": POSEBUSTERS_PARENT_DIR / "erbb2_7pcd" / "summary" / "posebusters_master_results.csv",  # PB_NOT_AVAILABLE
    "KDR": POSEBUSTERS_PARENT_DIR / "kdr_5ew3" / "summary" / "posebusters_master_results.csv",   # PB_NOT_AVAILABLE
    "PTGS2": POSEBUSTERS_PARENT_DIR / "ptgs2_5kir" / "summary" / "posebusters_master_results.csv",  # PB_NOT_AVAILABLE
    "ALOX-5": POSEBUSTERS_PARENT_DIR / "alox5_6n2w" / "summary" / "posebusters_master_results.csv",
    "PTGES": POSEBUSTERS_PARENT_DIR / "ptges_5tl9" / "summary" / "posebusters_master_results.csv",  # PB_NOT_AVAILABLE
}

# ============================================================
# TWO-TRACK ARCHITECTURE (NEW v1.2)
# ============================================================
# Track 1: Targets with PoseBusters validation (pose-level rescue logic)
# Track 2: Targets without PoseBusters (CNNscore fallback at ligand level)
TARGETS_WITH_PB = ["PPARA", "PPARD", "ALOX-5", "PPARG", "EGFR", "ERBB2", "KDR", "PTGS2", "PTGES"]
TARGETS_WITHOUT_PB = ""

# Excel sheet indices
SHEET_INTRA_LIGAND = 0  # Intra-Ligand_Poses (wide format, P1_*, P2_*, P3_*)
SHEET_INTER_LIGAND = 1  # Inter-Ligand_Ranking (179 rows, best pose per ligand)

# ============================================================
# CCLE / DEPMAP DATA PATHS (NEW v1.3)
# ============================================================
# CCLE (Cancer Cell Line Encyclopedia) from DepMap portal
# https://depmap.org/portal/data_page/?tab=allData&releasename=DepMap%20Public%2025Q3
#
# Scientific basis:
#   - TPM log2(x+1) is GOLD STANDARD for GRNBoost2 input
#   - NO additional transformation needed (already normalized)
#   - Matrix orientation: Rows=Samples, Cols=Genes (GRNBoost2 compatible)
# ============================================================
CCLE_DATA_DIR = Path(r"/home/labhhc5/Documents/workspace/D21/Duong Huy/pyscenic/pyscenic_data")

# Model metadata with OncotreeLineage column for filtering
CCLE_MODEL_CSV = CCLE_DATA_DIR / "Model.csv"

# Expression data - USER MUST DOWNLOAD from DepMap 25Q3
# File: OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv
CCLE_TPM_EXPRESSION_CSV = CCLE_DATA_DIR / "OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv"

# ============================================================
# PYSCENIC RESOURCE PATHS (NEW v1.3)
# ============================================================
# pySCENIC requires:
#   1. TF list (transcription factors)
#   2. Motif databases (.feather) - for cisTarget pruning
#   3. Motif annotations (.tbl) - TF-motif mappings
#
# Download from Aerts Lab:
#   https://resources.aertslab.org/cistarget/databases/
#   https://resources.aertslab.org/cistarget/motif2tf/
# ============================================================
SCENIC_RESOURCE_DIR = CCLE_DATA_DIR

# Human TF list - EXISTS (1839 TFs)
SCENIC_TF_LIST = SCENIC_RESOURCE_DIR / "hs_hgnc_tfs.txt"

# Motif databases directory - USER MUST DOWNLOAD
# Expected files: hg38_*.feather (e.g., hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather và hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather), https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg38/refseq_r80/mc_v10_clust/gene_based/
SCENIC_MOTIF_DB_DIR = Path(r"/home/labhhc5/Documents/workspace/D21/Duong Huy/pyscenic/pyscenic_data/scenic_motifs")

# Motif annotations - USER MUST DOWNLOAD
# Expected file: 	motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl (https://resources.aertslab.org/cistarget/motif2tf/)
SCENIC_MOTIF_ANNOTATIONS = Path(r"/home/labhhc5/Documents/workspace/D21/Duong Huy/pyscenic/pyscenic_data/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl")

# ============================================================
# LINEAGE DEFINITIONS (NEW v1.3)
# ============================================================
# OncotreeLineage level filtering for GRNBoost2 statistical power
# Minimum 40-50+ samples per lineage required for robust inference
#
# Architecture Decision: Use OncotreeLineage (NOT OncotreePrimaryDisease)
#   - OncotreePrimaryDisease too granular (often <20 samples)
#   - OncotreeLineage provides biological coherence + sample size
#
# Cell line mappings:
#   MCF-7, MDA-MB-231 -> Breast
#   A549             -> Lung
#   Jurkat           -> Lymphoid (T-cell Leukemia)
#   SW480            -> Bowel (Colorectal)
#   HEK-293          -> Kidney (reference, non-malignant)
# ============================================================
LINEAGES_OF_INTEREST = [
    "Breast",
    "Lung",
    "Lymphoid",
    "Bowel",
    "Kidney",
]

# Column name in Model.csv for lineage filtering
COL_ONCOTREE_LINEAGE = "OncotreeLineage"
COL_MODEL_ID = "ModelID"  # Sample identifier column

# Minimum samples per lineage for GRNBoost2
MIN_SAMPLES_PER_LINEAGE = 40

# ============================================================
# GRNBOOST2 PARAMETERS (NEW v1.3)
# ============================================================
# GRNBoost2 is the gradient boosting implementation of GENIE3
# Optimized for single-cell / bulk RNA-seq GRN inference
#
# Scientific basis:
#   - Uses Dask for distributed computing
#   - client_or_address='custom_multiprocessing' for local execution
#   - n_workers=4 recommended for RAM safety (16GB+ system)
#   - seed=42 for reproducibility
# ============================================================
GRNBOOST2_N_WORKERS = 4          # Number of parallel workers
GRNBOOST2_SEED = 42              # Random seed for reproducibility
GRNBOOST2_VERBOSE = True         # Print progress

# Output column names (after renaming)
COL_GRN_SOURCE = "Source"        # TF (transcription factor)
COL_GRN_TARGET = "Target"        # Target gene
COL_GRN_WEIGHT = "Weight"        # Importance score

# Raw GRNBoost2 output columns (before renaming)
COL_GRNBOOST2_TF = "TF"
COL_GRNBOOST2_TARGET = "target"
COL_GRNBOOST2_IMPORTANCE = "importance"

# ============================================================
# CISTARGET PARAMETERS (NEW v1.3)
# ============================================================
# cisTarget prunes GRNBoost2 adjacencies using motif enrichment
# Requires motif databases and annotations from Aerts Lab
#
# Graceful skip: If motif DB not available, skip cisTarget step
#   - Output raw GRNBoost2 adjacencies as preliminary GRN
#   - User can run cisTarget later when resources available
# ============================================================
CISTARGET_NES_THRESHOLD = 3.0    # Normalized Enrichment Score threshold
CISTARGET_RANK_THRESHOLD = 5000  # Top 5% of ranked genes
CISTARGET_AUC_THRESHOLD = 0.05   # AUC threshold for motif recovery
CISTARGET_MOTIF_SIMILARITY_FDR = 0.001  # FDR for motif similarity

# ============================================================
# LAYER 2A OUTPUT FILENAMES (NEW v1.3)
# ============================================================
L2A_ADJACENCIES_CSV = "L2A_GRNBoost2_Adjacencies_{lineage}.csv"
L2A_REGULONS_CSV = "L2A_Regulons_{lineage}.csv"
L2A_MASTER_ADJACENCIES_CSV = "L2A_Master_Adjacencies_AllLineages.csv"  # Raw GRNBoost2
L2A_MASTER_REGULONS_CSV = "L2A_Master_Regulons_AllLineages.csv"        # cisTarget pruned
L2A_AUC_MATRIX_CSV = "L2A_AUC_Matrix_{lineage}.csv"
L2A_QC_PLOTS_DIR = "L2A_QC_Plots"

# ============================================================
# LAYER 2B OUTPUT DIRECTORIES (NEW v1.4)
# ============================================================
LAYER2B_OUTPUT_DIR = PROJECT_ROOT / "outputs" / "Layer2B_Heterogeneous_RWR"

# ============================================================
# STRING PPI DATA PATHS (NEW v1.4)
# ============================================================
# STRING Database v12.0 - Human protein-protein interactions
# Download from: https://string-db.org/cgi/download
#
# Scientific basis:
#   - physical.links contains direct binding interactions only
#   - Score range 0-1000 (use >= 700 for high confidence)
#   - More specific than full links (avoids co-expression/text-mining noise)
# ============================================================
STRING_DATA_DIR = Path(r"/home/labhhc5/Documents/workspace/D21/Duong Huy/pyscenic/pyscenic_data")
STRING_PHYSICAL_LINKS_FILE = STRING_DATA_DIR / "9606.protein.physical.links.v12.0.txt.gz"
STRING_LINKS_FILE = STRING_DATA_DIR / "9606.protein.links.v12.0.txt.gz"
STRING_ALIASES_FILE = STRING_DATA_DIR / "9606.protein.aliases.v12.0.txt.gz"
STRING_INFO_FILE = STRING_DATA_DIR / "9606.protein.info.v12.0.txt.gz"

# STRING Configuration
STRING_SPECIES_PREFIX = "9606."
STRING_MIN_CONFIDENCE = 700  # High confidence threshold (0-1000 scale)

# ============================================================
# RWR PARAMETERS (NEW v1.4)
# ============================================================
# Random Walk with Restart via NetworkX PageRank
#
# Mathematical basis (NetworkX implementation):
#   P_{t+1} = alpha * W^T * P_t + (1-alpha) * P_0
#
# NetworkX convention:
#   alpha = DAMPING FACTOR (probability to continue walking)
#   (1 - alpha) = RESTART PROBABILITY (teleport back to seed nodes)
#   Setting alpha = 0.3 means: 30% random walk, 70% restart!
#
# Scientific basis:
#   Higher restart (lower alpha) -> signal stays tightly around drug targets.
#   Lower restart (higher alpha) -> signal diffuses broadly into the hairball.
#   Kohler et al. recommends a 0.7 restart probability -> We set alpha = 0.3
# ============================================================
RWR_ALPHA = 0.3              # Damping factor (yields 70% restart probability)
RWR_MAX_ITER = 100           # Maximum iterations for convergence
RWR_TOL = 1e-6               # Convergence tolerance
RWR_PSEUDO_COUNT = 0.001     # Epsilon for P0 calibration (prevents zero-division)

# ============================================================
# TARGET CELL LINE CONFIGURATION (NEW v1.4)
# ============================================================
# MCF-7: Human breast adenocarcinoma cell line
# Used for P0 calibration with tissue-specific expression
#
# Scientific basis:
#   - MCF-7 is estrogen receptor positive (ER+)
#   - Well-characterized transcriptome in CCLE/DepMap
#   - Representative of Breast lineage from Layer 2A
# ============================================================
TARGET_CELL_LINE = "MCF7" # Jurkat -> HEK-293 -> A549 -> SW480 -> MCF-7 và MDA-MB-231 (multi-lineage analysis)
TARGET_CELL_LINE_MODEL_ID = "ACH-000019"

# Alternative cell lines for multi-lineage analysis
CELL_LINE_MODELS = {
    "MCF7": "ACH-000019",       # Breast - ER+
    "MDA-MB-231": "ACH-000768", # Breast - Triple negative
    "A549": "ACH-000681",       # Lung - NSCLC
    "Jurkat": "ACH-000995",     # Lymphoid - T-cell leukemia
    "SW480": "ACH-000842",      # Bowel - Colorectal
    "HEK293": "ACH-001085",     # Kidney - Reference (non-cancer)
}

# ============================================================
# ACTIVE/INACTIVE THRESHOLD FOR DELTA-NETWORK (NEW v1.4)
# ============================================================
# Definition of "Active" vs "Inactive" compounds for comparative analysis
#
# Scientific basis:
#   - Active = Top 10% by CNN_VS (highest predicted binding affinity)
#   - Inactive = Bottom 10% by CNN_VS (weakest predicted binding)
#   - Delta = RWR_Active - RWR_Inactive eliminates housekeeping hubs
#   - Only drug-specific hubs survive delta filtering
# ============================================================
ACTIVE_PERCENTILE = 90       # Top 10% (100 - 90 = 10%)
INACTIVE_PERCENTILE = 10     # Bottom 10%

# ============================================================
# LAYER 2B OUTPUT FILENAMES (NEW v1.4)
# ============================================================
L2B_TOP_HUB_GENES = 50       # Number of top hub genes to extract per ligand
L2B_HUB_GENES_CSV = "L2B_Top50_Hub_Genes_{ligand_id}.csv"
L2B_RWR_SCORES_CSV = "L2B_RWR_Full_Scores_{ligand_id}.csv"
L2B_DELTA_NETWORK_CSV = "L2B_Delta_Network_Summary.csv"
L2B_MASTER_RWR_CSV = "L2B_Master_RWR_AllLigands.csv"
L2B_GRAPH_STATS_JSON = "L2B_Graph_Statistics.json"
L2B_GRAPH_GRAPHML = "L2B_Heterogeneous_Graph.graphml"

# ============================================================
# LAYER 2B COLUMN NAMES (NEW v1.4)
# ============================================================
COL_RWR_SCORE = "RWR_Score"
COL_RWR_RANK = "RWR_Rank"
COL_DELTA_SCORE = "Delta_Score"
COL_IS_DIRECT_TARGET = "Is_Direct_Target"
COL_EDGE_TYPE = "Edge_Type"
COL_GENE = "Gene"
COL_GENE_NORMALIZED = "Gene_Normalized"

# Edge type identifiers
EDGE_TYPE_SCENIC = "SCENIC_GRN"    # TF -> Target (one-way)
EDGE_TYPE_STRING = "STRING_PPI"    # Protein <-> Protein (bidirectional)

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
    "version": "1.4.0",
    "pipeline_version": "GNINA_v2.6.2 + pySCENIC + NetworkX RWR",
    "n_ligands": N_LIGANDS_EXPECTED,
    "n_targets": N_TARGETS,
    "targets": TARGET_NAMES,
    "docking_engine": "GNINA (CNN-based scoring)",
    "grn_engine": "pySCENIC (GRNBoost2 + cisTarget)",
    "rwr_engine": "NetworkX PageRank (RWR with alpha=0.7)",
    "ppi_database": "STRING v12.0 (physical.links, score >= 700)",
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
    # Layer 1 Rules
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
    # Layer 2A Rules
    "DL4_GRNBOOST2_ORIENTATION": (
        "DO NOT transpose expression matrix for GRNBoost2. "
        "CCLE TPM format (Rows=Samples, Cols=Genes) is already correct."
    ),
    "DL5_GENE_NAME_FORMAT": (
        "MUST strip Entrez IDs from gene names before GRNBoost2. "
        "Use: df.columns = [c.split(' (')[0] for c in df.columns]"
    ),
    # Layer 2B Rules (NEW v1.4)
    "DL6_DIRECTED_GRN": (
        "ABSOLUTELY DO NOT use Undirected Graph for transcription data. "
        "SCENIC edges must be one-way: TF -> Target. "
        "Signal cannot flow backwards from gene to transcription factor. "
        "This violates biological causality (Algebraic Conflict)."
    ),
    "DL7_P0_CALIBRATION": (
        "ABSOLUTELY DO NOT directly load CNN_VS into P0 vector. "
        "MUST multiply by expression level: P0,i = (CNN_VS_i * E_i) + epsilon. "
        "Raw CNN_VS alone creates Competitive Reservoir Paradox - "
        "a drug cannot influence a gene that isn't expressed."
    ),
    "DL8_DELTA_NETWORK": (
        "ABSOLUTELY DO NOT conclude MoA from single compound analysis. "
        "MUST run Delta-Network comparison: Active vs Inactive substances. "
        "Delta_Score = RWR_Active - RWR_Inactive. "
        "If graphs don't diverge, discard hypothesis (Against Confirmation Bias)."
    ),
}

# ============================================================
# UTILITY FUNCTIONS
# ============================================================

def normalize_gene_name(name: str) -> str:
    """
    Normalize gene names for cross-database matching.
    
    Architecture Decision (Lead Architect Approved):
        Auto hyphen-removal + uppercase normalization
        Example: ALOX-5 -> ALOX5, Alox5 -> ALOX5
    
    Scientific basis:
        - STRING uses ENSP IDs mapped to HGNC symbols (no hyphens)
        - CCLE/DepMap uses HGNC symbols (no hyphens)
        - Layer 1 may use variant naming (e.g., ALOX-5)
        - Normalization ensures consistent matching across layers
    
    Args:
        name: Gene name (may contain hyphens, mixed case)
    
    Returns:
        Normalized gene name (uppercase, no hyphens)
    """
    if name is None:
        return ""
    return name.replace("-", "").upper()


def print_config_summary():
    """Print current configuration summary - call at start of each notebook."""
    print("=" * 64)
    print("CONFIG SYSTEM - Central Configuration Hub v1.4")
    print("=" * 64)
    print(f"Project root:       {PROJECT_ROOT}")
    print(f"Docking data:       {DOCKING_PARENT_DIR}")
    print(f"Layer 1 output:     {LAYER1_OUTPUT_DIR}")
    print(f"Layer 2A output:    {LAYER2A_OUTPUT_DIR}")
    print(f"Layer 2B output:    {LAYER2B_OUTPUT_DIR}")
    print(f"Targets:            {N_TARGETS} ({', '.join(TARGET_NAMES)})")
    print(f"  Track 1 (PB):     {', '.join(TARGETS_WITH_PB)}")
    print(f"  Track 2 (no PB):  {', '.join(TARGETS_WITHOUT_PB)}")
    print(f"Ligands expected:   {N_LIGANDS_EXPECTED}")
    print(f"dG cutoff:          {DELTA_G_CUTOFF} kcal/mol")
    print(f"CNNscore fallback:  {CNNSCORE_FALLBACK_THRESHOLD}")
    print(f"alpha (CNN_VS):     {ALPHA_CNN_VS_WEIGHT}")
    print(f"Normalization:      {NORMALIZATION_METHOD}")
    print(f"Deadlock rules:     {len(DEADLOCK_RULES)} active")
    print("-" * 64)
    print("Layer 2A (SCENIC GRN):")
    print(f"  CCLE data dir:    {CCLE_DATA_DIR}")
    print(f"  Lineages:         {', '.join(LINEAGES_OF_INTEREST)}")
    print(f"  GRNBoost2 workers:{GRNBOOST2_N_WORKERS}")
    print(f"  Min samples:      {MIN_SAMPLES_PER_LINEAGE}")
    print("-" * 64)
    print("Layer 2B (Heterogeneous RWR):")
    print(f"  STRING data dir:  {STRING_DATA_DIR}")
    print(f"  STRING threshold: {STRING_MIN_CONFIDENCE} (high confidence)")
    print(f"  RWR alpha:        {RWR_ALPHA} (restart probability)")
    print(f"  Target cell line: {TARGET_CELL_LINE} ({TARGET_CELL_LINE_MODEL_ID})")
    print(f"  Active/Inactive:  Top {100-ACTIVE_PERCENTILE}% / Bottom {INACTIVE_PERCENTILE}%")
    print(f"  Top hub genes:    {L2B_TOP_HUB_GENES}")
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
    
    Layer 2B Steps:
      - directed_grn: Ensure SCENIC edges are directed (DL6)
      - p0_calibration: Ensure P0 = CNN_VS * Expression (DL7)
      - delta_network: Ensure Active vs Inactive comparison (DL8)
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
    elif step == "grnboost2_orientation":
        assert not context.get("transposed", False), (
            f"DEADLOCK VIOLATION - DL4:\n"
            f"{DEADLOCK_RULES['DL4_GRNBOOST2_ORIENTATION']}"
        )
    elif step == "gene_name_format":
        assert context.get("stripped_entrez", False), (
            f"DEADLOCK VIOLATION - DL5:\n"
            f"{DEADLOCK_RULES['DL5_GENE_NAME_FORMAT']}"
        )
    # Layer 2B Rules (NEW v1.4)
    elif step == "directed_grn":
        assert context.get("is_directed", False), (
            f"DEADLOCK VIOLATION - DL6:\n"
            f"{DEADLOCK_RULES['DL6_DIRECTED_GRN']}"
        )
    elif step == "p0_calibration":
        assert context.get("multiplied_by_expression", False), (
            f"DEADLOCK VIOLATION - DL7:\n"
            f"{DEADLOCK_RULES['DL7_P0_CALIBRATION']}"
        )
    elif step == "delta_network":
        assert context.get("compared_active_inactive", False), (
            f"DEADLOCK VIOLATION - DL8:\n"
            f"{DEADLOCK_RULES['DL8_DELTA_NETWORK']}"
        )


def check_scenic_resources() -> dict:
    """
    Check availability of pySCENIC resources.
    
    Returns:
        dict with keys: tf_list, motif_db, motif_annotations
        Values: True if available, False otherwise
    """
    return {
        "tf_list": SCENIC_TF_LIST.exists(),
        "motif_db": SCENIC_MOTIF_DB_DIR.exists() and str(SCENIC_MOTIF_DB_DIR) != "<USER_FILL_IN_PATH>",
        "motif_annotations": SCENIC_MOTIF_ANNOTATIONS.exists() and str(SCENIC_MOTIF_ANNOTATIONS) != "<USER_FILL_IN_PATH>",
        "ccle_model": CCLE_MODEL_CSV.exists(),
        "ccle_tpm": CCLE_TPM_EXPRESSION_CSV.exists(),
    }


def check_layer2b_resources() -> dict:
    """
    Check availability of Layer 2B resources.
    
    Returns:
        dict with keys for each required resource
        Values: True if available, False otherwise
    
    Resources checked:
        - STRING physical links (required)
        - STRING aliases (required for gene name mapping)
        - Layer 1 P0 Vector (required)
        - Layer 2A SCENIC GRN (optional - graceful degradation)
        - CCLE expression (required for P0 calibration)
    """
    # Layer 1 P0 Vector path
    l1_p0_path = LAYER1_OUTPUT_DIR / "L1_P0_Vector_Long.csv"
    
    # Layer 2A SCENIC GRN path (optional)
    l2a_grn_path = LAYER2A_OUTPUT_DIR / L2A_MASTER_REGULONS_CSV
    
    return {
        "string_physical_links": STRING_PHYSICAL_LINKS_FILE.exists(),
        "string_aliases": STRING_ALIASES_FILE.exists(),
        "string_info": STRING_INFO_FILE.exists(),
        "layer1_p0_vector": l1_p0_path.exists(),
        "layer2a_scenic_grn": l2a_grn_path.exists(),  # Optional - graceful degradation
        "ccle_model": CCLE_MODEL_CSV.exists(),
        "ccle_expression": Path(CCLE_DATA_DIR / "OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv").exists(),
    }


# ============================================================
# Self-test on import
# ============================================================
if __name__ == "__main__":
    print_config_summary()
    print("\nconfig_system.py v1.4 loaded successfully.")
    print(f"   Deadlock rules: {list(DEADLOCK_RULES.keys())}")
    pb_configured = sum(1 for v in POSEBUSTERS_PATHS.values() if v is not None)
    print(f"   PoseBusters paths configured: {pb_configured}/{len(POSEBUSTERS_PATHS)}")
    print(f"   Track 1 targets: {TARGETS_WITH_PB}")
    print(f"   Track 2 targets: {TARGETS_WITHOUT_PB}")
    print("\n   SCENIC Resources (Layer 2A):")
    resources = check_scenic_resources()
    for key, available in resources.items():
        status = "OK" if available else "MISSING"
        print(f"      {key}: {status}")
    print("\n   Layer 2B Resources:")
    l2b_resources = check_layer2b_resources()
    for key, available in l2b_resources.items():
        status = "OK" if available else "MISSING"
        if key == "layer2a_scenic_grn" and not available:
            status = "MISSING (will use STRING-only graph)"
        print(f"      {key}: {status}")
    print("\n   Gene Name Normalization Test:")
    test_cases = [("ALOX-5", "ALOX5"), ("Alox5", "ALOX5"), ("PPARG", "PPARG")]
    for original, expected in test_cases:
        result = normalize_gene_name(original)
        match = "PASS" if result == expected else "FAIL"
        print(f"      {original} -> {result} [{match}]")
