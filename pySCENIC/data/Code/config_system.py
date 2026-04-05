# ============================================================
# config_system.py - Central Configuration Hub v1.6
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
#          - High-Affinity/Low-Affinity percentile thresholds
#          - Delta-Network analysis parameters
#          - normalize_gene_name() utility function
#          - Deadlock rules DL6, DL7, DL8
#   v1.5 - Layer 3A: CRISPR Essentiality Validation configuration
#          - CRISPR DepMap data path (CRISPRGeneDependency.csv)
#          - 5 cancer cell lines (excluding HEK-293)
#          - P(dep) threshold configuration
#          - Positive delta filtering (DL10 compliance)
#          - Deadlock rules DL9, DL10, DL11, DL12
#          - check_layer3a_resources() function
#   v1.6 - Layer 3B: AUCell Functional State configuration
#          - AUCell threshold parameters (L3B_AUC_THRESHOLD, L3B_DELTA_AUC_THRESHOLD)
#          - GeneSignature requirements (ctxcore compliance)
#          - Jurkat (cancer) vs HEK-293 (normal) comparison
#          - Deadlock rules DL13, DL14, DL15, DL16
#          - check_layer3b_resources() function
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
# LAYER 3A OUTPUT DIRECTORIES (NEW v1.5)
# ============================================================
LAYER3A_OUTPUT_DIR = PROJECT_ROOT / "outputs" / "Layer3A_CRISPR_Validation"

# ============================================================
# LAYER 3B OUTPUT DIRECTORIES (NEW v1.6)
# ============================================================
LAYER3B_OUTPUT_DIR = PROJECT_ROOT / "outputs" / "Layer3B_AUCell_Selectivity"

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
# TARGET CELL LINE CONFIGURATION (NEW v1.4) (This is for layer 2B RWR analysis, different from the cell lines used for GRNBoost2 in layer 2A or layer 3A CRISPR validation)
# ============================================================
# MCF-7: Human breast adenocarcinoma cell line
# Used for P0 calibration with tissue-specific expression
#
# Scientific basis:
#   - MCF-7 is estrogen receptor positive (ER+)
#   - Well-characterized transcriptome in CCLE/DepMap
#   - Representative of Breast lineage from Layer 2A
# ============================================================
TARGET_CELL_LINE = "A549" # Jurkat -> HEK-293 -> A549 -> SW480 -> MCF-7 và MDA-MB-231 (multi-lineage analysis)
TARGET_CELL_LINE_MODEL_ID = "ACH-000681"

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
# HIGH-AFFINITY/LOW-AFFINITY THRESHOLD FOR DELTA-NETWORK (NEW v1.4)
# ============================================================
# Definition of "High-Affinity" vs "Low-Affinity" compounds for comparative analysis
#
# IMPORTANT TERMINOLOGY NOTE (v1.6.1):
#   CNN_VS is a PREDICTED binding score, NOT experimental biological activity.
#   We use "High-Affinity" / "Low-Affinity" instead of "Active" / "Inactive"
#   to avoid implying experimental validation (IC50/EC50) which does not exist.
#
# Scientific basis:
#   - High-Affinity = Top 10% by CNN_VS (highest predicted binding affinity)
#   - Low-Affinity = Bottom 10% by CNN_VS (weakest predicted binding)
#   - Delta = RWR_HighAffinity - RWR_LowAffinity eliminates housekeeping hubs
#   - Only drug-specific hubs survive delta filtering
#
# Cut-off justification (HTS hit rate ~1-10%):
#   - 5% (N~9): Too small, high variance from outliers
#   - 10% (N~18): Sweet spot - statistically meaningful, elite binders only
#   - 20% (N~36): Dilutes signal with mediocre binders
# ============================================================
ACTIVE_PERCENTILE = 90       # High-Affinity: Top 10% by CNN_VS (100 - 90 = 10%)
INACTIVE_PERCENTILE = 10     # Low-Affinity: Bottom 10% by CNN_VS

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
# CRISPR DEPMAP DATA PATH (NEW v1.5)
# ============================================================
# DepMap CRISPR Gene Dependency - P(dependency) probability scores
# Download from: https://depmap.org/portal/data_page/?tab=allData&releasename=DepMap%20Public%2025Q3
#
# Scientific basis:
#   - P(dep) is probability that gene knockout reduces cell viability
#   - Range 0.0 to 1.0 (0 = not essential, 1.0 = essential)
#   - P(dep) > 0.8 indicates high-confidence essential gene
#   - Matrix orientation: Rows=Samples (ModelID), Cols=Genes (Symbol + EntrezID)
#
# File format:
#   - Column 1: ModelID (e.g., ACH-000001)
#   - Columns 2-N: Gene names as "SYMBOL (EntrezID)" format
#   - Values: P(dependency) probability scores
# ============================================================
CRISPR_GENE_DEPENDENCY_CSV = Path(r"/home/labhhc5/Documents/workspace/D21/Duong Huy/pyscenic/pyscenic_data/CRISPRGeneDependency.csv")

# ============================================================
# LAYER 3A CRISPR CONFIGURATION (NEW v1.5)
# ============================================================
# CRISPR Essentiality Validation parameters
#
# SMART Goal:
#   S - Filter Top 50 positive delta genes against CRISPR P(dep)
#   M - Retain genes with P(dep) > 0.5 in any cancer cell line
#   A - Using existing CRISPRGeneDependency.csv (1,186 cell lines)
#   R - Identifies genes whose knockout causes cancer cell death
#   T - Single notebook execution (~5-10 minutes)
#
# FAIR Compliance:
#   F - Persistent gene IDs (HGNC symbols)
#   A - Standard CSV output format
#   I - Uses normalized gene names (compatible with STRING/CCLE)
#   R - Full provenance in validation summary JSON
# ============================================================

# P(dependency) threshold for essential genes
CRISPR_PDEP_THRESHOLD = 0.5  # P(dep) > 0.5 indicates essential gene

# Number of top positive delta genes to validate
L3A_TOP_DELTA_GENES = 20  
'''
"Phân tích độ nhạy (Sensitivity Analysis) o layer 3A chỉ ra rằng: 
Ngay tại ngưỡng Top 20, mạng lưới đã bao phủ trọn vẹn toàn bộ các đích tác dụng vật lý (4 Direct Targets). 
Đồng thời, mốc Top 20 mang lại tỷ lệ Tín hiệu/Nhiễu tốt nhất, với tỷ lệ nhiễu chỉ ở mức 75% (so với >85% ở các mốc mở rộng). 
Việc nới rộng ngưỡng cắt ra Top 40 hay Top 60 không mang lại thêm bất kỳ mỏ neo nào mới, mà chỉ làm pha loãng mạng lưới bởi các gen không thiết yếu. 
Do đó, Top 20 là ngưỡng cắt bảo thủ (conservative) và tối ưu nhất."
'''

# CRITICAL: Only use positive delta (DL10 compliance)
# Positive delta = genes elevated when ACTIVE drugs are present
# This isolates "Mechanism of Killing" from Active compounds only
L3A_DELTA_POSITIVE_ONLY = True

# ============================================================
# LAYER 3A CANCER CELL LINES (NEW v1.5)
# ============================================================
# 5 cancer cell lines for CRISPR validation
# HEK-293 EXCLUDED from CRISPR query (DL12 compliance)
# Selectivity against HEK-293 evaluated in Layer 3B via AUCell
#
# Scientific basis:
#   - CRISPR knockout in immortalized lines (HEK-293) artificially
#     inflates toxicity profiles
#   - Cancer selectivity must be evaluated using transcriptomic
#     comparison (AUCell), not CRISPR essentiality
# ============================================================
LAYER3A_CANCER_CELL_LINES = {
    "MCF7": "ACH-000019",        # Breast - ER+
    "MDAMB231": "ACH-000768",    # Breast - Triple negative (normalized name)
    "A549": "ACH-000681",        # Lung - NSCLC
    "Jurkat": "ACH-000995",      # Lymphoid - T-cell leukemia
    "SW480": "ACH-000842",       # Bowel - Colorectal
}
# NOTE: HEK293 (ACH-001085) intentionally EXCLUDED - handled in Layer 3B

# Primary cell line for validation (must match Layer 2B RWR analysis)
L3A_PRIMARY_CELL_LINE = "A549"
L3A_PRIMARY_MODEL_ID = "ACH-000681"

# ============================================================
# LAYER 3A OUTPUT FILENAMES (NEW v1.5)
# ============================================================
L3A_ESSENTIAL_TARGETS_CSV = "L3A_Essential_Targets.csv"
L3A_ALL_CANDIDATES_CSV = "L3A_All_Candidates.csv"
L3A_VALIDATION_SUMMARY_JSON = "L3A_Validation_Summary.json"
L3A_QC_PLOTS_DIR = "L3A_QC_Plots"

# ============================================================
# LAYER 3A COLUMN NAMES (NEW v1.5)
# ============================================================
COL_P_DEP = "P_dep"
COL_P_DEP_PREFIX = "P_dep_"  # Column prefix for cell line-specific P(dep)
COL_MAX_P_DEP = "Max_P_dep"
COL_ESSENTIAL_IN_N_LINES = "Essential_In_N_Lines"
COL_IS_ESSENTIAL = "Is_Essential"
COL_CRISPR_COVERAGE = "CRISPR_Coverage"
COL_DELTA_RANK = "Delta_Rank"

# CRISPR coverage status values
CRISPR_FOUND = "FOUND"
CRISPR_NOT_FOUND = "NOT_IN_CRISPR"

# ============================================================
# LAYER 3B AUCELL CONFIGURATION (NEW v1.6)
# ============================================================
# AUCell Functional State Analysis - Cancer Selectivity Validation
#
# SMART Goal:
#   S - Compare regulon activity (AUC) between Jurkat (cancer) and HEK-293 (normal)
#   M - Delta AUC = AUC_Jurkat - AUC_HEK293 quantifies cancer selectivity
#   A - Using pySCENIC AUCell with ctxcore GeneSignature objects
#   R - Validates that essential regulons are "bright" in cancer, "dim" in normal
#   T - Single notebook execution (~10-15 minutes)
#
# FAIR Compliance:
#   F - Persistent regulon IDs (TF_regulon format)
#   A - Standard CSV output with full AUC metrics
#   I - Uses ctxcore GeneSignature objects (pySCENIC compatible)
#   R - Full provenance in selectivity summary JSON
#
# Scientific Basis:
#   - AUCell ranks genes within each cell independently (DL15)
#   - Positive Delta AUC = regulon is MORE active in cancer cells
#   - Negative Delta AUC = regulon is LESS active in cancer = toxicity risk
#   - Threshold 0.05 Delta AUC provides robust selectivity signal
# ============================================================

# AUCell threshold for gene recovery curve
L3B_AUC_THRESHOLD = 0.05  # Default AUCell recovery curve threshold

# Delta AUC threshold for cancer selectivity
# Delta AUC > 0.05 = cancer selective (bright in cancer, dim in normal)
L3B_DELTA_AUC_THRESHOLD = 0.05

# Minimum genes required in a GeneSignature
L3B_MIN_GENES_IN_SIGNATURE = 5

# ============================================================
# LAYER 3B CELL LINE COMPARISON (NEW v1.6)
# ============================================================
# Jurkat (cancer) vs HEK-293 (normal) for selectivity calculation
#
# Scientific basis:
#   - Jurkat: T-cell leukemia (cancer) - Expected BRIGHT regulon activity
#   - HEK-293: Kidney non-malignant (normal) - Expected DIM regulon activity
#   - Delta AUC = AUC_Jurkat - AUC_HEK293
#   - Positive delta indicates cancer-selective regulon activation
# ============================================================
L3B_CANCER_CELL_LINE = "A549"
L3B_CANCER_MODEL_ID = "ACH-000681"  # Jurkat T-cell leukemia or other cancer cells

L3B_NORMAL_CELL_LINE = "HEK293"
L3B_NORMAL_MODEL_ID = "ACH-001085"  # HEK-293 non-malignant kidney

# ============================================================
# LAYER 3B OUTPUT FILENAMES (NEW v1.6)
# ============================================================
L3B_ACTIVE_REGULONS_CSV = "L3B_Active_Regulons.csv"          # Regulons with positive Delta AUC > threshold
L3B_ALL_REGULONS_AUC_CSV = "L3B_All_Regulons_AUC.csv"        # All regulons with AUC scores
L3B_L3A_SIGNATURE_AUC_CSV = "L3B_L3A_Signature_AUC.csv"      # Aggregate L3A signature AUC
L3B_SELECTIVITY_SUMMARY_JSON = "L3B_Selectivity_Summary.json"  # Metadata + statistics
L3B_QC_PLOTS_DIR = "L3B_QC_Plots"

# ============================================================
# LAYER 3B COLUMN NAMES (NEW v1.6)
# ============================================================
COL_AUC_JURKAT = "AUC_A549" # "AUC_Jurkat"
COL_AUC_HEK293 = "AUC_HEK293"
COL_DELTA_AUC = "Delta_AUC"
COL_REGULON_NAME = "Regulon_Name"
COL_N_GENES = "N_Genes"
COL_IS_SELECTIVE = "Is_Selective"
COL_SELECTIVITY_RANK = "Selectivity_Rank"

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
    "version": "1.6.0",
    "pipeline_version": "GNINA_v2.6.2 + pySCENIC + NetworkX RWR + CRISPR + AUCell",
    "n_ligands": N_LIGANDS_EXPECTED,
    "n_targets": N_TARGETS,
    "targets": TARGET_NAMES,
    "docking_engine": "GNINA (CNN-based scoring)",
    "grn_engine": "pySCENIC (GRNBoost2 + cisTarget)",
    "rwr_engine": "NetworkX PageRank (RWR with alpha=0.7)",
    "ppi_database": "STRING v12.0 (physical.links, score >= 700)",
    "crispr_database": "DepMap CRISPR Gene Dependency (P(dep) scores)",
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
        "MUST run Delta-Network comparison: High-Affinity vs Low-Affinity compounds. "
        "Delta_Score = RWR_HighAffinity - RWR_LowAffinity. "
        "If graphs don't diverge, discard hypothesis (Against Confirmation Bias). "
        "NOTE: 'High/Low-Affinity' refers to PREDICTED binding (CNN_VS), not experimental activity."
    ),
    # Layer 3A Rules (NEW v1.5)
    "DL9_CELL_LINE_MATCH": (
        "ABSOLUTELY DO NOT mix cell line contexts between layers. "
        "Layer 2B RWR cell line MUST match Layer 3A primary CRISPR query. "
        "Example: Jurkat RWR analysis -> Jurkat P(dep) as primary validation. "
        "Cross-cell-line validation is secondary/supplementary only."
    ),
    "DL10_DELTA_POSITIVE_ONLY": (
        "ABSOLUTELY DO NOT use absolute value of Delta_Score for L3A input. "
        "MUST filter Delta_Score > 0 (positive only) to isolate Mechanism of Killing. "
        "Negative delta genes = pathways elevated in INACTIVE (failing) drugs. "
        "Including negative delta contaminates essential gene list with non-MoA genes."
    ),
    "DL11_PDEP_THRESHOLD": (
        "ABSOLUTELY DO NOT claim gene is 'essential' with P(dep) <= 0.5. "
        "P(dep) > 0.5 is the Broad Institute's official threshold for CRISPR essentiality. (https://forum.depmap.org/t/dependency-threshold/498/2)"
        "Report all candidates but only label P(dep) > 0.5 as Is_Essential=True."
    ),
    "DL12_EXCLUDE_HEK293_CRISPR": (
        "ABSOLUTELY DO NOT include HEK-293 in Layer 3A CRISPR essentiality query. "
        "HEK-293 is immortalized non-cancer line - CRISPR knockout artificially "
        "inflates toxicity profiles due to viral immortalization machinery. "
        "Cancer selectivity vs HEK-293 MUST be evaluated via AUCell in Layer 3B."
    ),
    # Layer 3B Rules (NEW v1.6)
    "DL13_GENESIGNATURE_OBJECTS": (
        "ABSOLUTELY DO NOT use plain Python lists/arrays as AUCell signatures. "
        "MUST use ctxcore.genesig.GeneSignature objects. "
        "Plain lists will cause TypeError in pySCENIC AUCell function. "
        "Correct: GeneSignature(name='...', gene2weight={...})"
    ),
    "DL14_TPM_NORMALIZATION": (
        "ABSOLUTELY DO NOT use raw read counts for AUCell analysis. "
        "Expression matrix MUST be log2(x+1) normalized (CCLE TPM format). "
        "Raw counts create artificial variance that corrupts ranking."
    ),
    "DL15_AUCELL_CELL_INDEPENDENCE": (
        "AUCell rankings are computed INDEPENDENTLY within each cell. "
        "This makes cross-cell comparison robust to batch effects. "
        "Do NOT apply additional batch correction before AUCell - it auto-handles."
    ),
    "DL16_DELTA_AUC_INTERPRETATION": (
        "Positive Delta AUC (Jurkat - HEK293 > 0) = CANCER SELECTIVE. "
        "Regulon is BRIGHT in cancer, DIM in normal = therapeutic opportunity. "
        "Negative Delta AUC = regulon more active in normal = TOXICITY RISK. "
        "Only positive Delta AUC regulons are valid drug targets."
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
    print("CONFIG SYSTEM - Central Configuration Hub v1.6")
    print("=" * 64)
    print(f"Project root:       {PROJECT_ROOT}")
    print(f"Docking data:       {DOCKING_PARENT_DIR}")
    print(f"Layer 1 output:     {LAYER1_OUTPUT_DIR}")
    print(f"Layer 2A output:    {LAYER2A_OUTPUT_DIR}")
    print(f"Layer 2B output:    {LAYER2B_OUTPUT_DIR}")
    print(f"Layer 3A output:    {LAYER3A_OUTPUT_DIR}")
    print(f"Layer 3B output:    {LAYER3B_OUTPUT_DIR}")
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
    print(f"  High/Low-Affinity: Top {100-ACTIVE_PERCENTILE}% / Bottom {INACTIVE_PERCENTILE}% (by CNN_VS)")
    print(f"  Top hub genes:    {L2B_TOP_HUB_GENES}")
    print("-" * 64)
    print("Layer 3A (CRISPR Essentiality):")
    print(f"  CRISPR data:      {CRISPR_GENE_DEPENDENCY_CSV}")
    print(f"  P(dep) threshold: > {CRISPR_PDEP_THRESHOLD}")
    print(f"  Top delta genes:  {L3A_TOP_DELTA_GENES}")
    print(f"  Delta filter:     Positive only (DL10)")
    print(f"  Primary cell line:{L3A_PRIMARY_CELL_LINE} ({L3A_PRIMARY_MODEL_ID})")
    print(f"  Cancer cell lines:{', '.join(LAYER3A_CANCER_CELL_LINES.keys())}")
    print(f"  HEK-293:          EXCLUDED (DL12 - handled in L3B)")
    print("-" * 64)
    print("Layer 3B (AUCell Selectivity):")
    print(f"  AUC threshold:    {L3B_AUC_THRESHOLD}")
    print(f"  Delta AUC cutoff: > {L3B_DELTA_AUC_THRESHOLD} (cancer selective)")
    print(f"  Min genes/sig:    {L3B_MIN_GENES_IN_SIGNATURE}")
    print(f"  Cancer cell line: {L3B_CANCER_CELL_LINE} ({L3B_CANCER_MODEL_ID})")
    print(f"  Normal cell line: {L3B_NORMAL_CELL_LINE} ({L3B_NORMAL_MODEL_ID})")
    print(f"  Signature type:   ctxcore GeneSignature (DL13)")
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
      - p0_calibration: Ensure P0 = P0_raw_weight * Expression (DL7)
      - delta_network: Ensure High-Affinity vs Low-Affinity comparison (DL8)
    
    Layer 3A Steps (NEW v1.5):
      - cell_line_match: Ensure L2B cell line matches L3A primary (DL9)
      - delta_positive_only: Ensure only positive Delta_Score used (DL10)
      - pdep_threshold: Ensure P(dep) > 0.5 for essential label (DL11)
      - exclude_hek293: Ensure HEK-293 not in CRISPR query (DL12)
    
    Layer 3B Steps (NEW v1.6):
      - genesignature_objects: Ensure ctxcore GeneSignature used (DL13)
      - tpm_normalization: Ensure expression is log2(x+1) TPM (DL14)
      - aucell_independence: Informational - AUCell handles batch effects (DL15)
      - delta_auc_interpretation: Ensure positive Delta AUC for selectivity (DL16)
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
    # Layer 3A Rules (NEW v1.5)
    elif step == "cell_line_match":
        l2b_cell_line = context.get("l2b_cell_line", "")
        l3a_primary = context.get("l3a_primary_cell_line", L3A_PRIMARY_CELL_LINE)
        assert l2b_cell_line == l3a_primary, (
            f"DEADLOCK VIOLATION - DL9:\n"
            f"{DEADLOCK_RULES['DL9_CELL_LINE_MATCH']}\n"
            f"L2B cell line: {l2b_cell_line}, L3A primary: {l3a_primary}"
        )
    elif step == "delta_positive_only":
        assert context.get("used_positive_only", False), (
            f"DEADLOCK VIOLATION - DL10:\n"
            f"{DEADLOCK_RULES['DL10_DELTA_POSITIVE_ONLY']}"
        )
    elif step == "pdep_threshold":
        threshold = context.get("threshold", 0)
        assert threshold >= CRISPR_PDEP_THRESHOLD, (
            f"DEADLOCK VIOLATION - DL11:\n"
            f"{DEADLOCK_RULES['DL11_PDEP_THRESHOLD']}\n"
            f"Used threshold: {threshold}, Required: > {CRISPR_PDEP_THRESHOLD}"
        )
    elif step == "exclude_hek293":
        cell_lines_used = context.get("cell_lines", [])
        assert "HEK293" not in cell_lines_used and "HEK-293" not in cell_lines_used, (
            f"DEADLOCK VIOLATION - DL12:\n"
            f"{DEADLOCK_RULES['DL12_EXCLUDE_HEK293_CRISPR']}\n"
            f"Cell lines used: {cell_lines_used}"
        )
    # Layer 3B Rules (NEW v1.6)
    elif step == "genesignature_objects":
        assert context.get("used_genesignature", False), (
            f"DEADLOCK VIOLATION - DL13:\n"
            f"{DEADLOCK_RULES['DL13_GENESIGNATURE_OBJECTS']}"
        )
    elif step == "tpm_normalization":
        assert context.get("is_tpm_normalized", False), (
            f"DEADLOCK VIOLATION - DL14:\n"
            f"{DEADLOCK_RULES['DL14_TPM_NORMALIZATION']}"
        )
    elif step == "aucell_independence":
        # Informational check - AUCell auto-handles batch effects
        pass  # DL15 is informational, not a hard constraint
    elif step == "delta_auc_interpretation":
        delta_auc = context.get("delta_auc", 0)
        assert delta_auc > 0, (
            f"DEADLOCK VIOLATION - DL16:\n"
            f"{DEADLOCK_RULES['DL16_DELTA_AUC_INTERPRETATION']}\n"
            f"Delta AUC: {delta_auc} (negative = toxicity risk)"
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


def check_layer3a_resources() -> dict:
    """
    Check availability of Layer 3A CRISPR Essentiality resources.
    
    Returns:
        dict with keys for each required resource
        Values: True if available, False otherwise
    
    Resources checked:
        - CRISPR Gene Dependency CSV (required)
        - Layer 2B Delta Network Summary (required)
        - Layer 3A output directory (will be created if missing)
    
    DL9 Compliance:
        Also validates that L2B was run with matching cell line.
    """
    # Layer 2B Delta Network path (required input)
    # Note: Check for Jurkat-specific output directory
    l2b_delta_path = PROJECT_ROOT / "outputs" / "Layer2B_Heterogeneous_RWR_Jurkat" / L2B_DELTA_NETWORK_CSV
    
    # Alternative: check generic L2B output if Jurkat-specific doesn't exist
    l2b_generic_path = LAYER2B_OUTPUT_DIR / L2B_DELTA_NETWORK_CSV
    
    return {
        "crispr_gene_dependency": CRISPR_GENE_DEPENDENCY_CSV.exists(),
        "layer2b_delta_network_jurkat": l2b_delta_path.exists(),
        "layer2b_delta_network_generic": l2b_generic_path.exists(),
        "layer3a_output_dir_exists": LAYER3A_OUTPUT_DIR.exists(),
        "primary_cell_line_match": L3A_PRIMARY_CELL_LINE in LAYER3A_CANCER_CELL_LINES,
    }


def check_layer3b_resources() -> dict:
    """
    Check availability of Layer 3B AUCell Functional State resources.
    
    Returns:
        dict with keys for each required resource
        Values: True if available, False otherwise
    
    Resources checked:
        - L3A Essential Targets CSV (required - from Layer 3A)
        - L2A Master Regulons CSV (required - for per-TF signatures)
        - CCLE TPM Expression CSV (required - log2(x+1) normalized)
        - Layer 3B output directory (will be created if missing)
    
    DL13 Compliance:
        Validates that ctxcore is available for GeneSignature creation.
    DL14 Compliance:
        Validates that expression data is TPM normalized (CCLE format).
    """
    # Layer 3A Essential Targets (required input)
    l3a_essential_path = LAYER3A_OUTPUT_DIR / L3A_ESSENTIAL_TARGETS_CSV
    
    # Layer 2A Master Regulons (required for per-TF signatures)
    l2a_regulons_path = LAYER2A_OUTPUT_DIR / L2A_MASTER_REGULONS_CSV
    
    # CCLE TPM Expression (required)
    ccle_tpm_path = Path(CCLE_TPM_EXPRESSION_CSV)
    
    # Check ctxcore availability
    try:
        from ctxcore.genesig import GeneSignature
        ctxcore_available = True
    except ImportError:
        ctxcore_available = False
    
    return {
        "l3a_essential_targets": l3a_essential_path.exists(),
        "l2a_master_regulons": l2a_regulons_path.exists(),
        "ccle_tpm_expression": ccle_tpm_path.exists(),
        "layer3b_output_dir_exists": LAYER3B_OUTPUT_DIR.exists(),
        "ctxcore_available": ctxcore_available,
        "cancer_cell_line_configured": L3B_CANCER_MODEL_ID == "ACH-000995",  # Jurkat
        "normal_cell_line_configured": L3B_NORMAL_MODEL_ID == "ACH-001085",  # HEK-293
    }


# ============================================================
# Self-test on import
# ============================================================
if __name__ == "__main__":
    print_config_summary()
    print("\nconfig_system.py v1.6 loaded successfully.")
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
    print("\n   Layer 3A Resources:")
    l3a_resources = check_layer3a_resources()
    for key, available in l3a_resources.items():
        status = "OK" if available else "MISSING"
        if key == "layer3a_output_dir_exists" and not available:
            status = "MISSING (will be created)"
        if key == "layer2b_delta_network_generic" and not available:
            # Only show warning if Jurkat-specific also missing
            if not l3a_resources.get("layer2b_delta_network_jurkat", False):
                status = "MISSING (run Layer 2B first)"
        print(f"      {key}: {status}")
    print("\n   Layer 3B Resources:")
    l3b_resources = check_layer3b_resources()
    for key, available in l3b_resources.items():
        status = "OK" if available else "MISSING"
        if key == "layer3b_output_dir_exists" and not available:
            status = "MISSING (will be created)"
        if key == "ctxcore_available" and not available:
            status = "MISSING (install ctxcore: pip install ctxcore)"
        print(f"      {key}: {status}")
    print("\n   Gene Name Normalization Test:")
    test_cases = [("ALOX-5", "ALOX5"), ("Alox5", "ALOX5"), ("PPARG", "PPARG")]
    for original, expected in test_cases:
        result = normalize_gene_name(original)
        match = "PASS" if result == expected else "FAIL"
        print(f"      {original} -> {result} [{match}]")
