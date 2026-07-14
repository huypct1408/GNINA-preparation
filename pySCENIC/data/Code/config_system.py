# ============================================================
# config_system.py - Central Configuration Hub v3.4.0
# ============================================================
# Purpose: Single source of truth for:
#   - Layer 1: Thermodynamic Gate (PoseBusters-based core seed)
#   - Layer 2: GRN PANDA/LIONESS (LOCAL FROZEN PRIORS ONLY)
#   - Layer 3: RWR-MultiXrank (universal multilayer network)
#   - Layer 4: ORA + Redundancy Filter
#
# Scientific Integrity Policy (STRICT):
#   - NO hardcoded notebook parameters outside config
#   - NO on-the-fly prior generation inside any notebook
#   - NO heuristic fallback branches in runtime pipeline
#   - ALL priors must be precomputed, frozen, provenance-tracked
#   - Config stores decisions and constraints, not post-hoc justification
#
# Version history:
#   v2.0.0  Layer 1 unified PoseBusters architecture
#   v2.1.0  Layer 2 PANDA/LIONESS architecture
#   v2.1.1  cisTarget motif prior + STRING physical PPI prior
#   v2.1.2  STRICT LOCAL-FROZEN-PRIORS policy
#   v2.2.0  LIONESS TARGETED SINGLE MODEL mode
#   v2.3.0  SPLIT NOTEBOOK ARCHITECTURE (02A / 02B)
#   v3.0.0  LAYER 3 RWR-MultiXrank added (DL3_01..DL3_05)
#   v3.1.0  LAYER 3 PPI FROZEN ARTIFACT separated from Layer 2
#   v3.3.0  LAYER 3 GEOMETRIC-MEAN + SINGLE-RUN POLICY
#   v3.4.0  LAYER 2 MULTI-MODEL LIONESS
#            ADDED: LAYER2_RUN_MODE = TARGETED_LIONESS_MULTI_MODEL
#            ADDED: L2_TARGET_MODELS_LIST (70 ACH ModelIDs)
#            ADDED: L2_MULTI_MODEL_* policy constants
#            ADDED: DL2_25..DL2_28 multi-model deadlock rules
#            KEPT:  L2_TARGET_MODEL_INPUT as backward-compat alias
#            KEPT:  same_lineage_from_frozen_lineage_inputs policy
#            KEPT:  lineage-shared PANDA aggregate (one e^(α) per lineage)
#   v3.5.0  LAYER 3 BATCH (N LIONESS samples) + Jupyter monolith
#            ADDED: L3_BATCH_* discovery / ledger / FAIR catalog keys
#            ADDED: LAYER3_NOTEBOOK_BATCH_NAME, LAYER3_BATCH_LOG_FILE
#            ADDED: L3_SCIENCE_VERSION / L3_ARCHITECTURE_VERSION locks
#            ADDED: L3_BATCH_RUN_LAMBDA_BENCHMARK, L3_BATCH_KEEP_WORKDIR
#            KEPT:  DL3_01..07, geom-mean global score, CLASS_SHARED policy
#            NOTE:  Stage 9b multi-compound all_results path is FORBIDDEN
#                   under DL3_07 (use structural seed audit only)
#
# PPI artifact provenance:
#   Layer 2: TF-only PPI (for PANDA TF cooperativity prior)
#            → LOCAL_PPI_PRIOR_TXT (2-col, no header)
#   Layer 3: Gene-level PPI (for RWR-M multiplex layer)
#            → LOCAL_PPI_PRIOR_L3_TSV (3-col, with header)
#            Built by: 01b_build_ppi_prior_l3_from_string.py
#            SHA-256: 209941abc26f562a1ba9adc000c09b2fc12354da765abc7cab7b718f05567fbe
# ============================================================

from pathlib import Path
import logging

# ============================================================
# SYSTEM PATHS
# ============================================================

PROJECT_ROOT = Path(__file__).resolve().parent

# ============================================================
# LAYER 1 — OUTPUT DIRECTORY
# ============================================================

LAYER1_OUTPUT_DIR = PROJECT_ROOT / "outputs" / "Layer1_Consensus_Target_Prior"

# ============================================================
# LAYER 1 — INPUT ARTIFACT PATHS
# ============================================================
# These are FROZEN input artifacts.
# They must NOT be modified after being declared.
# Any change requires updating the manifest and SHA-256.
#
# L1_DATA_DIR: dedicated directory for Layer 1 frozen priors

L1_DATA_DIR = PROJECT_ROOT / "pyscenic_data"

# File 1: Consensus target list
# Format: 1 column, header = "Target", one gene symbol per line
# Inclusion rule: supported by >= L1_MIN_SUPPORT_METHODS algorithms
LOCAL_CONSENSUS_TARGET_LIST_TXT = (
    L1_DATA_DIR / "huymly_boron_consensus_target_ge3.txt"
    # Note: rename to "manual_consensus_targets_support_ge3.txt"
    # to revert to project-179 compound set
)

# File 2: Ligand (compound) IDs
# Format: 1 column, header = "ligand_id", one ID per line
LOCAL_CONSENSUS_LIGAND_IDS_TXT = (
    L1_DATA_DIR / "manual_ligand_ids.txt"
)

# File 3: Provenance manifest (JSON)
LOCAL_CONSENSUS_MANIFEST_JSON = (
    L1_DATA_DIR / "consensus_target_prior_manifest.json"
)

# SHA-256 checksums of frozen artifacts
# Leave as empty string "" if file not yet hashed.
LOCAL_CONSENSUS_TARGET_LIST_EXPECTED_SHA256 = ""
LOCAL_CONSENSUS_LIGAND_IDS_EXPECTED_SHA256 = ""

# ============================================================
# LAYER 1 — SEED MODE
# ============================================================
# CONSENSUS_GLOBAL_UNWEIGHTED:
#   T_d = T (global, same for all compounds)
#   p0(i|d) = 1/|T| if i in T, else 0
#   |T| computed at runtime from artifact, NOT hardcoded
#
# Do NOT change this value without re-running all downstream layers.

L1_SEED_MODE = "CONSENSUS_GLOBAL_UNWEIGHTED"
L1_SEED_FORMULA = "p0(i) = 1/|T| if i in T, else 0 (unweighted)"

# Minimum number of prediction methods that must support a target
# for it to be included in the consensus prior.
# This value describes the ARTIFACT that was pre-built —
# it is NOT enforced at runtime (runtime trusts the frozen artifact).
# It MUST match the manifest.
L1_MIN_SUPPORT_METHODS = 3

# Seed weighting policy
CORE_SEED_UNWEIGHTED = True  # MUST remain True for this mode
SEED_METHOD = "uniform_over_consensus_global_target_set"

# ============================================================
# LAYER 1 — GENE NORMALIZATION POLICY
# ============================================================
# Applied to the target list at load time only.
# Rule: uppercase, strip whitespace, remove hyphens.
# Example: "ALOX-5" -> "ALOX5"
# This must match normalize_gene_name() defined below.

L1_GENE_NORMALIZATION_POLICY = "uppercase_strip_remove_hyphens"

# ============================================================
# LAYER 1 — COLUMN NAMES (OUTPUT CONTRACT)
# ============================================================
# These column names are shared with Layer 3 (do NOT rename).

COL_LIGAND_ID = "lig_id"
COL_TARGET = "target"
COL_SMILES = "smiles"
COL_LIGAND_NAME = "orig_name"

# ============================================================
# LAYER 1 — OUTPUT FILENAMES
# ============================================================
# L1_P0_VECTOR_LONG_CSV is the canonical downstream contract.
# Layer 3 reads this file. Do NOT rename.

L1_P0_VECTOR_LONG_CSV = "L1_P0_Vector_Long.csv"
L1_P0_MATRIX_WIDE_CSV = "L1_P0_Matrix_Wide.csv"
L1_TARGET_AUDIT_CSV = "L1_Consensus_Target_Audit.csv"
L1_GATE_REPORT_JSON = "L1_Gate_Report.json"
L1_QC_PLOTS_DIR = "L1_QC_Plots"

# Kept for backward compatibility with Layer 3 checks
L1_P0_VECTOR_NORMALIZED_CSV = "L1_P0_Vector_Normalized.csv"
L1_TD_DISTRIBUTION_CSV = "L1_Td_Distribution.csv"

# ============================================================
# LAYER 1 — EXPECTED COUNTS (soft validation only)
# ============================================================
# These are INFORMATIONAL only — code reads from file and counts.
# Mismatch triggers a WARNING, not a hard failure,
# because the file is the ground truth.

N_LIGANDS_EXPECTED = 179
# N_TARGETS_EXPECTED is intentionally NOT hardcoded.
# Layer 1 will count from the artifact at runtime.

# ============================================================
# CCLE / DEPMAP DATA PATHS
# ============================================================

CCLE_DATA_DIR = Path(
    "/home/labhhc5/Documents/workspace/D21/Duong Huy/pyscenic/pyscenic_data"
)

CCLE_MODEL_CSV = CCLE_DATA_DIR / "Model.csv"
CCLE_TPM_EXPRESSION_CSV = (
    CCLE_DATA_DIR / "OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv"
)

# ============================================================
# LAYER 2 — DIRECTORIES
# ============================================================

LAYER2_GRN_DATA_DIR = Path(
    "/home/labhhc5/Documents/workspace/D21/Duong Huy/pyscenic"
)
LAYER2_GRN_OUTPUT_DIR = LAYER2_GRN_DATA_DIR / "outputs" / "Layer2_GRN"

# ============================================================
# LAYER 2 — NOTEBOOK ARCHITECTURE (v2.4.0 multi-model)
# ============================================================

LAYER2_NOTEBOOK_ARCHITECTURE = "SPLIT_02A_FOUNDATION__02B_TARGETED_MULTI"
LAYER2_NOTEBOOK_02A_NAME = "02A_Layer2_PANDA_FOUNDATION.ipynb"
LAYER2_NOTEBOOK_02B_NAME = "02B_Layer2_LIONESS_TARGETED.ipynb"
LAYER2_NOTEBOOK_02B_MULTI_NAME = "02B_Layer2_LIONESS_MULTI_MODEL.py"

# ============================================================
# LAYER 2 — LINEAGES
# ============================================================

LINEAGES_OF_INTEREST = ["Breast", "Lung", "Lymphoid", "Bowel", "Kidney"]
COL_ONCOTREE_LINEAGE = "OncotreeLineage"
COL_MODEL_ID = "ModelID"
MIN_SAMPLES_PANDA = 20

# ============================================================
# LAYER 2 — PANDA / SNAIL
# ============================================================

NETZOOPY_REQUIRED_VERSION = "0.11.0"
SCENIC_FORBIDDEN_ON_BULK = True
GLOBAL_QNORM_FORBIDDEN = True
SNAIL_NORMALIZATION_REQUIRED = True
GRN_BUILDER = "PANDA"
PANDA_MODE = "union"
PANDA_COMPUTING = "gpu"
PANDA_PRECISION = "double"
PANDA_SAVE_MEMORY = True
PANDA_SAVE_TMP = False
PANDA_REMOVE_MISSING = True
PANDA_KEEP_EXPRESSION_MATRIX_CORE = False
PANDA_KEEP_EXPRESSION_MATRIX_LIONESS = True
PANDA_WITH_HEADER = False
PANDA_ALPHA = 0.1
PANDA_WEIGHTS_UNMODIFIED = True

SNAIL_AGGREGATION = "median"
SNAIL_THRESHOLD = 0.25
SNAIL_CUTOFF = 0.15
SNAIL_LABEL_RESTORATION_POLICY = "restore_exact_input_labels_after_qsmooth"

# ============================================================
# LAYER 2 — LIONESS MULTI-MODEL MODE (v2.4.0 / config v3.4.0)
# ============================================================
# Scientific note:
#   LIONESS itself is not single-sample-only. The previous
#   TARGETED_LIONESS_SINGLE_MODEL mode was a resource wrapper.
#   Multi-model mode keeps the SAME lineage-frozen reference policy
#   and shares one PANDA aggregate e^(α) per lineage across all
#   requested targets in that lineage (DL2-25).
#
# Integrity constraints:
#   - Targets declared ONLY via L2_TARGET_MODELS_LIST (DL2-26)
#   - No SNAIL recompute in 02B (DL2-22)
#   - No cross-lineage pooling of the aggregate network
#   - Per-target export only after flatten-order verification (DL2-27)
#   - Every declared ID must appear in the run ledger (DL2-28)

LAYER2_RUN_MODE = "TARGETED_LIONESS_MULTI_MODEL"

# Canonical multi-model target list (ACH ModelIDs, as provided).
# Order is preserved for provenance; duplicates are rejected at runtime.
L2_TARGET_MODELS_LIST = [
    "ACH-000828",
    "ACH-000117",
    "ACH-000017",
    "ACH-000554",
    "ACH-000276",
    "ACH-000818",
    "ACH-000856",
    "ACH-000223",
    "ACH-000019",  # MCF7
    "ACH-000330",
    "ACH-000725",
    "ACH-000643",
    "ACH-000711",
    "ACH-000248",
    "ACH-000857",
    "ACH-002499",
    "ACH-000196",
    "ACH-000258",
    "ACH-000691",
    "ACH-000374",
    "ACH-000148",
    "ACH-001705",
    "ACH-001392",
    "ACH-000876",
    "ACH-000668",
    "ACH-002921",
    "ACH-001396",
    "ACH-001419",
    "ACH-001389",
    "ACH-001683",
    "ACH-000568",
    "ACH-000111",
    "ACH-000621",
    "ACH-000755",
    "ACH-000028",
    "ACH-000277",
    "ACH-001394",
    "ACH-000721",
    "ACH-000902",
    "ACH-000859",
    "ACH-000927",
    "ACH-000768",  # MDA-MB-231
    "ACH-000349",
    "ACH-001662",
    "ACH-001393",
    "ACH-000573",
    "ACH-000783",
    "ACH-000759",
    "ACH-000910",
    "ACH-000147",
    "ACH-000930",
    "ACH-000624",
    "ACH-000699",
    "ACH-000044",
    "ACH-000352",
    "ACH-001390",
    "ACH-001388",
    "ACH-001395",
    "ACH-000536",
    "ACH-000642",  # HMEL (Breast non-cancerous immortalized)
    "ACH-001819",
    "ACH-002401",
    "ACH-000212",
    "ACH-000849",
    "ACH-000097",
    "ACH-000288",
    "ACH-001820",
    "ACH-002399",
    "ACH-001391",
    "ACH-000934",
]

# Soft expected count — runtime counts the list; mismatch = hard fail
L2_TARGET_MODELS_LIST_EXPECTED_N = 70

# Backward-compatible single-model alias (first list element only).
# Single-model notebooks may still read this; multi-model MUST use the list.
L2_TARGET_MODEL_INPUT = L2_TARGET_MODELS_LIST[0]
L2_TARGET_CELL_LINE = L2_TARGET_MODEL_INPUT

# Optional alias map (name → ModelID). Multi-model list uses raw ACH IDs.
# Kept for single-model UX and documentation.
CELL_LINE_MODELS = {
    "MCF7": "ACH-000019",
    "MDA-MB-231": "ACH-000768",
    "A549": "ACH-000681",
    "Jurkat": "ACH-000995",
    "SW480": "ACH-000842",
    "HEK293": "ACH-001085",
    "SALE": "ACH-000064",   # Lung, Non-Cancerous Immortalized Lung Cells
    "HEKTE": "ACH-000049",  # Kidney, Non-Cancerous Immortalized Kidney Cells
    "HMEL": "ACH-000642",   # Breast, Non-Cancerous Immortalized Breast Cells
    "HCC1739BL": "ACH-001828",  # Lymphoid, B lymphoblastoid cells
}

# Reference / export policy (UNCHANGED scientific contract)
LIONESS_REFERENCE_POLICY = "same_lineage_from_frozen_lineage_inputs"
LIONESS_FAIL_IF_TARGET_ABSENT = True  # per-target in single-model mode
LIONESS_EXPORT_ONLY_TARGET = True
ZERO_VARIANCE_GENE_POLICY = "drop_variance_leq_0_exact"
PANDA_UNION_GENE_EXPANSION_EXPECTED = True
LIONESS_NPY_FLATTEN_ORDER = "verified_at_runtime_against_panda_export"
NETZOOPY_KNOWN_BUG_LIONESS_SINGLE_SAMPLE = (
    "AttributeError: 'Lioness' object has no attribute 'total_lioness_network'"
)

# Multi-model specific policies
# "skip_missing_with_ledger": record SKIP for IDs absent from frozen lineage
#   cohort / not in Model.csv lineage of interest, continue other targets.
# "fail": any missing ID aborts the entire multi-model run.
L2_MULTI_MODEL_ON_MISSING = "skip_missing_with_ledger"
L2_MULTI_MODEL_REQUIRE_SHARED_PANDA_PER_LINEAGE = True
L2_MULTI_MODEL_ALLOW_CROSS_LINEAGE_POOLING = False  # MUST remain False
L2_MULTI_MODEL_DEDUPLICATE_LIST = False  # duplicates → hard fail (integrity)
L2_MULTI_MODEL_GATE_REPORT_JSON = "L2B_LIONESS_MultiModel_Gate_Report.json"
L2_MULTI_MODEL_LEDGER_JSON = "L2B_LIONESS_MultiModel_Run_Ledger.json"
L2_MULTI_MODEL_LEDGER_TSV = "L2B_LIONESS_MultiModel_Run_Ledger.tsv"
L2_MULTI_MODEL_LINEAGE_SUMMARY_JSON = (
    "L2B_LIONESS_MultiModel_{lineage}_summary.json"
)

# ============================================================
# LAYER 2 — 02A FOUNDATION FROZEN INPUT CONTRACT
# ============================================================

L2A_FROZEN_INPUTS_DIRNAME = "frozen_lineage_inputs"
L2A_FROZEN_INPUTS_DIR = LAYER2_GRN_OUTPUT_DIR / L2A_FROZEN_INPUTS_DIRNAME
FOUNDATION_MANIFEST_VERSION = "8F_v1"
L2A_FOUNDATION_MANIFEST_JSON = "foundation_manifest.json"
L2A_EXPR_SNAIL_PARQUET = "{lineage}_expr_snail.parquet"
L2A_SAMPLE_ORDER_JSON = "{lineage}_sample_order.json"
L2A_GENES_POST_QC_JSON = "{lineage}_genes_post_qc.json"
FROZEN_EXPR_STORAGE_FORMAT = "parquet"
FROZEN_EXPR_ENGINE = "pyarrow"
FOUNDATION_PRIOR_HASH_VALIDATION_REQUIRED = True

# ============================================================
# LAYER 2 — 02B TARGETED LIONESS OUTPUT CONTRACT
# ============================================================

L2B_TARGET_SAMPLE_ORDER_JSON = "{sample}_sample_order.json"
L2B_LIONESS_RUN_MANIFEST_JSON = "{sample}_lioness_run_manifest.json"
L2B_NETZOOPY_RAW_DIR = "_netzoopy_raw_{sample}"

# ============================================================
# LAYER 2 — PRIOR POLICY (STRICT LOCAL FROZEN ONLY)
# ============================================================

LOCAL_FROZEN_PRIORS_ONLY = True
PRIOR_SOURCE_POLICY = "LOCAL_FROZEN_PRIORS_ONLY"

LOCAL_MOTIF_PRIOR_TXT = LAYER2_GRN_DATA_DIR / "motif_prior_500bp.txt"
LOCAL_MOTIF_PRIOR_10KB_TXT = LAYER2_GRN_DATA_DIR / "motif_prior_10kb.txt"
LOCAL_PPI_PRIOR_TXT = LAYER2_GRN_DATA_DIR / "ppi_prior.txt"
ALL_TFS_HG38_TXT = L1_DATA_DIR / "allTFs_hg38.txt"

LOCAL_MOTIF_PRIOR_MANIFEST_JSON = (
    LAYER2_GRN_DATA_DIR / "motif_prior_500bp_manifest.json"
)
LOCAL_MOTIF_PRIOR_10KB_MANIFEST_JSON = (
    LAYER2_GRN_DATA_DIR / "motif_prior_10kb_manifest.json"
)
LOCAL_PPI_PRIOR_MANIFEST_JSON = (
    LAYER2_GRN_DATA_DIR / "ppi_prior_manifest.json"
)

# ============================================================
# LAYER 2 — PRIOR PROVENANCE
# ============================================================

PRIOR_JASPAR_VERSION = "cisTarget_v10_clust_hg38_refseq_r80_JASPAR2020_basis"
PRIOR_STRING_VERSION = "STRINGv12.0_physical_links"
PRIOR_GENERATED_DATE = "2026-05-04"
MOTIF_PRIOR_CORE_WINDOW = "500bp_up_100bp_down"
MOTIF_PRIOR_SENS_WINDOW = "10kbp_up_10kbp_down"
MOTIF_PRIOR_WEIGHT_MODE = "binary_fixed_1"
MOTIF_PRIOR_TOP_N_PER_MOTIF = 5000
PPI_PRIOR_SCORE_CUTOFF = 700
PPI_PRIOR_NETWORK_TYPE = "physical_links_only"
STRING_SCORE_CUTOFF = PPI_PRIOR_SCORE_CUTOFF

# ============================================================
# LAYER 2 — GRN COLUMN NAMES
# ============================================================

COL_GRN_SOURCE = "Source"
COL_GRN_TARGET = "Target"
COL_GRN_WEIGHT = "Weight"
EDGE_TYPE_GRN = "PANDA_GRN"
EDGE_TYPE_STRING = "STRING_PPI"
EDGE_TYPE_MULTIPLEX = "MULTIPLEX"

# ============================================================
# LAYER 2 — OUTPUT FILENAMES
# ============================================================

L2_GRN_PANDA_TSV = "Z_{lineage}_PANDA.tsv"
L2_GRN_LIONESS_DIR = "LIONESS_{lineage}"
L2_GRN_LIONESS_TSV = "Z_{sample}_LIONESS.tsv"
L2_GRN_QC_SUMMARY_CSV = "L2_GRN_QC_Summary.csv"
L2_GRN_QC_PLOTS_DIR = "L2_GRN_QC_Plots"
L2_GRN_PRIOR_MANIFEST_JSON = "L2_GRN_Prior_Manifest.json"
L2_GRN_UNIVERSE_CSV = "L2_GRN_Gene_Universe.csv"
L2A_FOUNDATION_GATE_REPORT_JSON = "L2A_PANDA_Foundation_Gate_Report.json"
L2B_LIONESS_GATE_REPORT_JSON = "L2B_LIONESS_Targeted_Gate_Report.json"
L2_GRN_GATE_REPORT_JSON = "L2_GRN_Gate_Report.json"  # deprecated alias

# ============================================================
# LAYER 3 — DIRECTORIES
# ============================================================

LAYER3_OUTPUT_DIR = PROJECT_ROOT / "outputs" / "Layer3_RWR_MultiXrank"

# ============================================================
# LAYER 3 — MULTIXRANK VERSION LOCK
# ============================================================
# All parameters confirmed empirically against MultiXrank 0.3
# Toy data diagnostic runs: 2026-05-09
# Ref: Baptista et al. 2022, Nature Communications

L3_MULTIXRANK_REQUIRED_VERSION = "0.3"

# graph_type encoding confirmed from MultiplexLayer source & official docs:
"""
Code thật (commit 8738164) — chuẩn, https://github.com/anthbapt/multixrank/blob/873816443225851abbb5a02c5c66bb4c8d3e3bf6/multixrank/MultiplexLayer.py
Python

# weighted?
if self.graph_type[1] == '1':   # ← CHỮ SỐ THỨ HAI
    # đọc 3 cột, weight = float64

# directed?
if self.graph_type[0] == '1':   # ← CHỮ SỐ THỨ NHẤT
    networkx_graph_obj = networkx.DiGraph()
else:
    networkx_graph_obj = networkx.Graph()

Ý nghĩa là:
Index	Ý nghĩa (CODE)
graph_type[0]	directed? '0' = undirected, '1' = directed
graph_type[1]	weighted? '0' = unweighted (2 cột), '1' = weighted (3 cột)

Trong thiết kế code của MultiXrank, người ta thiết kế chữ số dầu tiên là có hướng/ vô hướng, còn chữ số sau biểu thị có trọng số/ không có trọng số
        """
L3_GRAPH_TYPE_WEIGHTED_UNDIRECTED = "01"  # PPI
L3_GRAPH_TYPE_WEIGHTED_DIRECTED = "11"    # GRN
L3_GRAPH_TYPE_UNWEIGHTED_UNDIRECTED = "00"  # bipartite gene-drug

# Confirmed: multiplex column in output dtype=object (strings '1','2')
L3_GENE_MULTIPLEX_ID = "1"
L3_DRUG_MULTIPLEX_ID = "2"

# ============================================================
# LAYER 3 — RWR PARAMETERS
# ============================================================

L3_RWR_RESTART_PROB = 0.7
L3_DELTA_GENE = 0.5
L3_DELTA_DRUG = 0.0
L3_LAMBDA_DEFAULT = 0.5
L3_LAMBDA_BENCH_VALUES = [0.2, 0.5, 0.8]
L3_ETA_DEFAULT = [0.5, 0.5]
L3_TAU_GENE_DEFAULT = [0.5, 0.5]
L3_TAU_DRUG_DEFAULT = [1.0]

# ============================================================
# LAYER 3 — GRN WEIGHT TRANSFORM POLICY
# ============================================================

L3_GRN_WEIGHT_TRANSFORM = "clip_nonpositive"
L3_GRN_WEIGHT_TRANSFORM_DOC = (
    "Remove edges with PANDA/LIONESS z-score <= 0 before RWR. "
    "Negative weights cause non-stochastic transition matrix. "
    "Confirmed necessary: MultiXrank 0.3 does not clip internally."
)

# ============================================================
# LAYER 3 — DRUG MONOPLEX POLICY
# ============================================================

L3_DRUG_MONOPLEX_POLICY = "self_loop_per_drug_node"
L3_DRUG_MONOPLEX_DOC = (
    "Drug monoplex written as self-loop per drug (D D 1.0). "
    "MultiXrank removes self-loops internally; drug nodes persist. "
    "Drug scores propagate via bipartite only. No drug-drug diffusion."
)

# ============================================================
# LAYER 3 — FROZEN PPI ARTIFACT (GENE-LEVEL STRING, v3.1.0)
# ============================================================
# CRITICAL DISTINCTION:
#   LOCAL_PPI_PRIOR_TXT       → Layer 2 PANDA only (TF-restricted, 2-col)
#   LOCAL_PPI_PRIOR_L3_TSV    → Layer 3 RWR only  (gene-level, 3-col weighted)
#
# Layer 2 PPI: TF–TF edges only (PANDA TF cooperativity prior)
# Layer 3 PPI: All gene–gene physical interactions (RWR diffusion backbone)
#              Includes kinases, receptors, enzymes NOT in motif TF set
#              e.g. EGFR, ERBB2, KDR, PTGS2, PTGES, ALOX5
#
# Built by: 01b_build_ppi_prior_l3_from_string.py
# SHA-256:  209941abc26f562a1ba9adc000c09b2fc12354da765abc7cab7b718f05567fbe
# Stats:    86519 gene pairs, 10746 unique genes

L3_STRING_PHYSICAL_LINKS = (
    CCLE_DATA_DIR / "9606.protein.physical.links.v12.0.txt.gz"
)
L3_STRING_INFO = CCLE_DATA_DIR / "9606.protein.info.v12.0.txt.gz"

# Frozen artifact paths
LOCAL_PPI_PRIOR_L3_TSV = LAYER2_GRN_DATA_DIR / "ppi_prior_l3_with_score.tsv"
LOCAL_PPI_PRIOR_L3_MANIFEST_JSON = (
    LAYER2_GRN_DATA_DIR / "ppi_prior_l3_manifest.json"
)

# Expected SHA-256 of artifact (from 01b build run 2026-05-09)
LOCAL_PPI_PRIOR_L3_EXPECTED_SHA256 = (
    "209941abc26f562a1ba9adc000c09b2fc12354da765abc7cab7b718f05567fbe"
)

# Policy constants (must match manifest generated by 01b script)
L3_PPI_STRING_VERSION = "STRINGv12.0"
L3_PPI_NETWORK_TYPE = "physical_links_only"
L3_PPI_MIN_STRING_SCORE = 700
L3_PPI_GENE_MAPPING = "ENSP_to_preferred_name_via_STRING_info"
L3_PPI_COLLAPSE_RULE = "max_score_per_gene_pair"
L3_PPI_TF_RESTRICTED = False  # MUST be False for Layer 3
L3_PPI_SOURCE_POLICY = "LOCAL_FROZEN_STRING_GENE_LEVEL_FOR_LAYER3_ONLY"

# Required schema for Layer 3 PPI artifact
# Header required; exact column names enforced at load time
L3_PPI_REQUIRED_COLUMNS = ["GeneA", "GeneB", "Score"]

# ============================================================
# LAYER 3 — OUTPUT FILENAMES
# ============================================================

L3_GENE_SCORES_PER_LAYER_TSV = "{prefix}_gene_scores_per_layer.tsv"
L3_GENE_TOP_N_PREVIEW_TSV = "{prefix}_gene_top{n}_PREVIEW_ONLY.tsv"
L3_DRUG_SCORES_TSV = "{prefix}_drug_scores.tsv"
L3_PROVENANCE_JSON = "{prefix}_provenance.json"
L3_LOG_FILE = LAYER3_OUTPUT_DIR / "layer3_execution.log"
L3_TOP_N_PREVIEW = 100
L3_LAMBDA_BENCHMARK_TSV = "lambda_sensitivity_benchmark.tsv"

# ============================================================
# LAYER 3 — STABILITY THRESHOLDS
# ============================================================

L3_JACCARD_THRESHOLD = 0.8
L3_SPEARMAN_THRESHOLD = 0.9

# ============================================================
# LAYER 3 — GLOBAL SCORE / SINGLE-RUN POLICY (v3.3.0)
# ============================================================

L3_GLOBAL_GENE_SCORE_METHOD = "geometric_mean_across_layers"
L3_GLOBAL_GENE_SCORE_REFERENCE = (
    "Valdeolivas et al. 2019, Bioinformatics, doi:10.1093/bioinformatics/bty543"
)
L3_CLASS_SHARED_NODE_PREFIX = "CLASS_SHARED"
L3_USE_SINGLE_REPRESENTATIVE_RUN = True
L3_SINGLE_RUN_EQUIVALENCE_METHOD = "structural_equivalence_audit"

# Housekeeping / ORA handoff policy
L3_HK_GMT_PATH = CCLE_DATA_DIR / "HOUNKPE_HOUSEKEEPING_GENES.v2026.1.Hs.gmt"
L3_ORA_TOP_N_NO_HK = 100
L3_ORA_MIN_NONHK = 20
L3_ORA_BACKGROUND_RECOMMENDATION = "full_layer3_gene_universe"

# ============================================================
# LAYER 3 — BATCH / MULTI-SAMPLE (architecture v4.0.0 monolith)
# ============================================================
# Canonical Jupyter entry (single file, all stages inline):
#   notebooks/03_Layer3_RWR_BATCH_JUPYTER.py
# Optional modular package (same science): layer3/
#
# Sample discovery for N LIONESS GRNs (e.g. 35 or 70 from 02B multi-model).
#
# L3_BATCH_SOURCE:
#   "l2_multimodel_ledger" — COMPLETED rows from L2 multi-model ledger (default)
#   "explicit_list"        — L3_SAMPLE_MODELS_LIST + L3_BATCH_LINEAGE_FILTER
#   "glob"                 — discover Z_*_LIONESS.tsv under LIONESS_{lineage}/
#
# Scientific locks (must match notebook SCIENCE_VERSION / ARCH_VERSION):
L3_SCIENCE_VERSION = "v3.3.0"
L3_ARCHITECTURE_VERSION = "v4.0.0-monolith"
LAYER3_NOTEBOOK_BATCH_NAME = "03_Layer3_RWR_BATCH_JUPYTER.py"
# Path uses PROJECT_ROOT (LOG_DIR is defined later in LOGGING section)
LAYER3_BATCH_LOG_FILE = (
    PROJECT_ROOT / "logs" / "03_Layer3_RWR_MultiXrank_BATCH.log"
)

L3_BATCH_SOURCE = "l2_multimodel_ledger"
L3_BATCH_LINEAGE_FILTER = "Breast"  # None = all lineages in ledger
# When source=explicit_list, set IDs here (None = unused)
L3_SAMPLE_MODELS_LIST = None  # e.g. ["ACH-000019", "ACH-000768", ...]
L3_BATCH_KEEP_WORKDIR = False
L3_BATCH_WORKERS = 1  # sequential default; use job arrays for horizontal scale
L3_BATCH_RUN_LAMBDA_BENCHMARK = True  # Stage 10 on first completed sample only

# Batch durable artifacts (under LAYER3_OUTPUT_DIR / batch_{id}/)
L3_BATCH_GATE_REPORT_JSON = "L3_Batch_Gate_Report.json"
L3_BATCH_LEDGER_JSON = "L3_Batch_Ledger.json"
L3_BATCH_LEDGER_TSV = "L3_Batch_Ledger.tsv"
L3_BATCH_CATALOG_JSON = "L3_Batch_Catalog.json"  # FAIR index

# Per-sample CLASS_SHARED filename contracts (Layer 4 globs these)
L3_CLASS_SHARED_TOP100_NOHK = "{prefix}_CLASS_SHARED_top100_noHK.tsv"
L3_CLASS_SHARED_GENE_UNIVERSE_TXT = "{prefix}_CLASS_SHARED_gene_universe.txt"
L3_CLASS_SHARED_GENE_UNIVERSE_TSV = "{prefix}_CLASS_SHARED_gene_universe.tsv"
L3_CLASS_SHARED_GENE_UNIVERSE_MANIFEST = (
    "{prefix}_CLASS_SHARED_gene_universe_manifest.json"
)
L3_CLASS_SHARED_FULL_PER_LAYER = "{prefix}_CLASS_SHARED_full_per_layer.tsv"
L3_CLASS_SHARED_GLOBAL_GEOMEAN = (
    "{prefix}_CLASS_SHARED_global_geomean_all_genes.tsv"
)
L3_CLASS_SHARED_PROVENANCE = "{prefix}_CLASS_SHARED_provenance.json"

# Forbidden under DL3_07 (legacy Stage 9b multi-compound score concordance
# that required all_results and max_score ranking). Structural seed audit only.
L3_FORBID_MULTI_COMPOUND_SCORE_CONCORDANCE = True
L3_FORBID_MAX_SCORE_AS_GLOBAL_RANKING = True

# ============================================================
# LAYER 4 — DIRECTORIES
# ============================================================

LAYER4_OUTPUT_DIR = PROJECT_ROOT / "outputs" / "Layer4_ORA"

# ============================================================
# LAYER 4 — INPUT ARTIFACT CONTRACT (from Layer 3 / Stage 9c-9d)
# ============================================================
# These are the durable file artifacts exported by Layer 3.
# Layer 4 MUST NOT depend on runtime Python variables from Layer 3.
#
# Observed real-data schemas (confirmed 2026-05-11):
#
#   top100_noHK.tsv columns:
#     rank_noHK, original_rank, gene_symbol, gene_symbol_norm,
#     score_GRN_layer, score_PPI_layer, rwr_score_geom_mean,
#     layer_support_min, dominant_layer
#
#   gene_universe.txt:
#     one gene symbol per line, 10286 genes total
#
#   gene_universe.tsv columns:
#     gene_symbol, gene_symbol_norm, source_partition
#
# Naming convention: files are prefixed with "lambda{value}_"
# e.g. "lambda0.5_CLASS_SHARED_top100_noHK.tsv"

L4_QUERY_FILENAME_PATTERN = "*_CLASS_SHARED_top100_noHK.tsv"
L4_BACKGROUND_TXT_PATTERN = "*_CLASS_SHARED_gene_universe.txt"
L4_BACKGROUND_TSV_PATTERN = "*_CLASS_SHARED_gene_universe.tsv"
L4_BACKGROUND_MANIFEST_PATTERN = "*_CLASS_SHARED_gene_universe_manifest.json"

# Required columns in query file (from observed real data)
L4_QUERY_REQUIRED_COLUMNS = [
    "rank_noHK",
    "original_rank",
    "gene_symbol",
    "gene_symbol_norm",
    "score_GRN_layer",
    "score_PPI_layer",
    "rwr_score_geom_mean",
    "layer_support_min",
    "dominant_layer",
]
L4_QUERY_SCORE_COLUMN = "rwr_score_geom_mean"
L4_QUERY_LAYER_SUPPORT_COLUMN = "layer_support_min"
# Column used as query gene list for ORA
L4_QUERY_GENE_COLUMN = "gene_symbol"

# Expected background size (from Stage 9d verified output)
L4_BACKGROUND_POLICY = "explicit_layer3_gene_universe"
L4_REQUIRE_BACKGROUND_MANIFEST = True
L4_REQUIRE_QUERY_SUBSET_BACKGROUND = True
L4_MIN_BACKGROUND_GENES_WARN = 1000

# ============================================================
# LAYER 4 — GMT LIBRARY
# ============================================================
# Observed real data (2026-05-11):
#   MERGED_Reference_Universe.gmt
#   n_terms = 4652
#   Format: TERM \t DESC/URL \t GENE1 \t GENE2 ...

L4_GMT_PATH = Path(
    "/home/labhhc5/Documents/workspace/D21/Duong Huy/pyscenic/"
    "pyscenic_data/MERGED_Reference_Universe.gmt"
)
L4_GMT_LABEL = "MERGED_Reference_Universe"

# ============================================================
# LAYER 4 — TARGET CELL LINE (must be declared explicitly)
# ============================================================
# When multiple CLASS_SHARED directories (multiple cell lines) exist,
# Layer 4 must know exactly which one to process.
# Change these two lines when switching cell line for ORA.

L4_TARGET_LINEAGE = "Breast"
L4_TARGET_MODEL = "ACH-000019"  # MCF-7

# ============================================================
# LAYER 4 — ORA ENGINE
# ============================================================
# Confirmed behavior from toy run + production run (gseapy 1.2.1):
#   N = len(background)
#   n = len(query ∩ background)
#   M = len(pathway ∩ background)
#   k = len(query ∩ pathway ∩ background)
#   Overlap column = "k/M" (string)
#   Genes column = ";" separated
#
# background=None uses N = |union(all GMT genes)| = 17071
# background=gene_universe uses N = 10286
# Layer 4 MUST use background=gene_universe (explicit)

L4_GSEAPY_REQUIRED_VERSION = "1.2.1"
# L4_BACKGROUND_POLICY already set above to explicit_layer3_gene_universe

# ORA result schema confirmed from gseapy 1.2.1
L4_GSEAPY_EXPECTED_COLUMNS = [
    "Gene_set", "Term", "Overlap", "P-value",
    "Adjusted P-value", "Odds Ratio", "Combined Score", "Genes",
]

# ============================================================
# LAYER 4 — ORA THRESHOLDS
# ============================================================

L4_FDR_CUTOFF = 0.05
L4_MIN_OVERLAP = 3
L4_TOP_TERM_DISPLAY = 25
L4_PVALUE_CONCORDANCE_TOL = 0.01

# ============================================================
# LAYER 4 — OVERLAP SANITIZATION (Excel date-corruption prevention)
# ============================================================
# GSEApy outputs Overlap as "k/M" (e.g. "12/11").
# Excel interprets "/" as date separator → "12-Nov".
# Sanitization: "k/M" → "k | M" + integer columns k_overlap, M_pathway_size.

L4_OVERLAP_SANITIZE = True

# ============================================================
# LAYER 4 — REDUNDANCY FILTER
# ============================================================
# Redundancy filter runs on SIGNIFICANT output (confirmed by user decision).
# References:
#   [1] Choobdar et al. Nat Methods 16, 843-852 (2019)
#   [2] Lamparter et al. PLoS Comput Biol 12(1):e1004714 (2016)

L4_REDUNDANCY_INPUT_POLICY = "SIGNIFICANT_only"
# ElasticNet parameters
L4_ELASTICNET_L1_RATIO = 0.5
L4_ELASTICNET_CV_FOLDS = 3
# Pairwise overlap filter parameters
L4_OVERLAP_FDR_ALPHA = 0.05
L4_SUBMODULE_S_THRESHOLD = 0.5
L4_SUBMODULE_J_THRESHOLD = 0.5
L4_SIZE_DIFF_FALLBACK = 5
# Module size bounds (Choobdar et al.)
L4_MODULE_MIN_GENES = 3
L4_MODULE_MAX_GENES = 100
# CRITICAL: Redundancy filter MUST use pathway gene sets restricted
# to the same background as ORA, not full GMT raw sets.
L4_REDUNDANCY_RESTRICT_TO_BACKGROUND = True

# ============================================================
# LAYER 4 — OUTPUT FILENAMES
# ============================================================

L4_ALL_RESULTS_TSV = "L4_ORA_ALL_results.tsv"
L4_SIGNIFICANT_TSV = "L4_ORA_SIGNIFICANT.tsv"
L4_VERIFICATION_TSV = "L4_ORA_VERIFICATION.tsv"
L4_ORA_MANIFEST_JSON = "L4_ORA_manifest.json"
L4_ORA_GATE_REPORT_JSON = "L4_ORA_Gate_Report.json"
# Redundancy filter outputs
L4_NONREDUNDANT_TSV = "L4_NON_REDUNDANT_PATHWAYS.tsv"
L4_ELIMINATED_TSV = "L4_ELIMINATED_PATHWAYS.tsv"
L4_OVERLAP_PAIRS_TSV = "L4_OVERLAP_PAIRS.tsv"
L4_FILTRATION_SUMMARY_TSV = "L4_FILTRATION_SUMMARY.tsv"
L4_REDUNDANCY_MANIFEST_JSON = "L4_REDUNDANCY_manifest.json"
# QC outputs
L4_JACCARD_HEATMAP_PNG = "L4_JACCARD_HEATMAP.png"
L4_JACCARD_HEATMAP_SVG = "L4_JACCARD_HEATMAP.svg"
L4_JACCARD_MATRIX_CSV = "L4_JACCARD_MATRIX.csv"

# ============================================================
# LOGGING
# ============================================================

LOG_FORMAT = "%(asctime)s | %(levelname)-8s | %(name)s | %(message)s"
LOG_DATEFMT = "%Y-%m-%d %H:%M:%S"
LOG_LEVEL = logging.INFO

LOG_DIR = PROJECT_ROOT / "logs"
LAYER1_LOG_FILE = LOG_DIR / "01_Layer1_Thermodynamic_Gate.log"
LAYER2A_LOG_FILE = LOG_DIR / "02A_Layer2_PANDA_FOUNDATION.log"
LAYER2B_LOG_FILE = LOG_DIR / "02B_Layer2_LIONESS_TARGETED.log"
LAYER2B_MULTI_LOG_FILE = LOG_DIR / "02B_Layer2_LIONESS_MULTI_MODEL.log"
LAYER3_LOG_FILE = LOG_DIR / "03_Layer3_RWR_MultiXrank.log"
# Keep batch log path consistent with LOG_DIR (re-bind after LOG_DIR exists)
LAYER3_BATCH_LOG_FILE = LOG_DIR / "03_Layer3_RWR_MultiXrank_BATCH.log"
LAYER2_GRN_LOG_FILE = LAYER2A_LOG_FILE  # deprecated alias
L4_LOG_FILE = LOG_DIR / "04_Layer4_ORA.log"


def setup_logger(
    name: str = "Pipeline",
    logfile=None,
    reset_handlers: bool = False,
) -> logging.Logger:
    logger = logging.getLogger(name)
    if reset_handlers:
        for h in list(logger.handlers):
            logger.removeHandler(h)
            try:
                h.close()
            except Exception:
                pass
    logger.setLevel(LOG_LEVEL)
    logger.propagate = False
    has_stream = any(
        isinstance(h, logging.StreamHandler)
        and not isinstance(h, logging.FileHandler)
        for h in logger.handlers
    )
    if not has_stream:
        sh = logging.StreamHandler()
        sh.setFormatter(logging.Formatter(LOG_FORMAT, LOG_DATEFMT))
        logger.addHandler(sh)
    if logfile is not None:
        logfile = Path(logfile)
        logfile.parent.mkdir(parents=True, exist_ok=True)
        existing = [
            h
            for h in logger.handlers
            if isinstance(h, logging.FileHandler)
            and Path(getattr(h, "baseFilename", "")) == logfile.resolve()
        ]
        if not existing:
            fh = logging.FileHandler(str(logfile), mode="w", encoding="utf-8")
            fh.setFormatter(logging.Formatter(LOG_FORMAT, LOG_DATEFMT))
            logger.addHandler(fh)
    return logger


# ============================================================
# DEADLOCK RULES
# ============================================================

DEADLOCK_RULES = {
    # ── Layer 1 ─────────────────────────────────────────────
    "DL1C_NO_RUNTIME_TARGET_DISCOVERY": (
        "Layer 1 MUST read target prior from a frozen artifact file. "
        "Runtime target discovery (docking scan, expression filter, etc.) "
        "is FORBIDDEN. The artifact must have a known SHA-256."
    ),
    "DL1C_NO_INLINE_GENE_LIST": (
        "The consensus target list MUST NOT be hardcoded inline in notebooks "
        "or config. It must be read from LOCAL_CONSENSUS_TARGET_LIST_TXT."
    ),
    "DL19C_UNWEIGHTED_SEED": (
        "Seed MUST remain unweighted: p0(i) = 1/|T| if i in T, else 0. "
        "Prediction confidence scores, network degree, or expression level "
        "MUST NOT be used to weight p0 at Layer 1 or downstream."
    ),
    "DL20C_SUPPORT_RULE_FIXED_AT_SOURCE": (
        f"Target inclusion rule is fixed: supported by >= "
        f"{L1_MIN_SUPPORT_METHODS} independent algorithms. "
        "This rule applies to the artifact creation step only. "
        "Layer 1 runtime trusts the frozen artifact as ground truth. "
        "The rule must be documented in the manifest."
    ),
    "DL21C_NO_DOWNSTREAM_REWEIGHTING": (
        "The P0 vector exported by Layer 1 MUST NOT be modified, "
        "re-weighted, or re-filtered downstream (Layer 2, 3, 4). "
        "Any modification must go through a new Layer 1 run."
    ),
    "DL22C_SEED_COLLAPSE_ACKNOWLEDGED": (
        "Because all compounds share the same global target set, "
        "seed collapse by design is expected and MUST be explicitly "
        "acknowledged in the gate report. It is NOT a bug."
    ),
    # ── Layer 2 ─────────────────────────────────────────────
    "DL2_01_SCENIC_FORBIDDEN": (
        "SCENIC/GRNBoost2 MUST NOT be used on bulk CCLE."
    ),
    "DL2_02_MIN_SAMPLES": (
        f"Each lineage MUST have >= {MIN_SAMPLES_PANDA} samples."
    ),
    "DL2_03_NO_GLOBAL_QNORM": (
        "Global quantile normalization is forbidden."
    ),
    "DL2_04_SNAIL_REQUIRED": (
        "SNAIL MUST be applied per lineage BEFORE PANDA (02A only). "
        "SNAIL MUST NOT be recomputed in notebook 02B."
    ),
    "DL2_05_LOCAL_FROZEN_PRIORS_ONLY": (
        "Layer 2 notebooks MUST use only LOCAL frozen priors."
    ),
    "DL2_06_FIXED_PRIOR_PROVENANCE": (
        f"Core motif prior: {MOTIF_PRIOR_CORE_WINDOW}, binary weight, "
        f"top_n={MOTIF_PRIOR_TOP_N_PER_MOTIF}. "
        f"PPI prior: STRING physical links, score>={PPI_PRIOR_SCORE_CUTOFF}."
    ),
    "DL2_07_GRN_BUILDER": (
        f"GRN builder MUST be {GRN_BUILDER} (netZooPy {NETZOOPY_REQUIRED_VERSION})."
    ),
    "DL2_08_PANDA_WEIGHTS_UNMODIFIED": (
        "PANDA force values MUST NOT be manually edited."
    ),
    "DL2_09_LIONESS_AFTER_PANDA": (
        "LIONESS MUST run only after successful PANDA convergence."
    ),
    "DL2_10_LIONESS_SAME_PRIORS": (
        "LIONESS MUST use identical motif and PPI files as PANDA."
    ),
    "DL2_11_NON_RANDOMNESS_TEST": (
        "Non-randomness test is NOT_IMPLEMENTED. Gate report MUST record NOT_RUN."
    ),
    "DL2_12_CONSISTENT_GENE_IDS": (
        "All gene IDs MUST be HGNC symbols; Entrez suffixes stripped."
    ),
    "DL2_13_NO_MANUAL_EDGE_EDIT": (
        "PANDA/LIONESS outputs MUST NOT be manually edited."
    ),
    "DL2_14_TARGET_MODEL_MUST_EXIST": (
        "Requested target ModelID MUST exist uniquely in Model.csv."
    ),
    "DL2_15_TARGET_LINEAGE_FROM_METADATA_ONLY": (
        "Target lineage MUST be derived from Model.csv metadata."
    ),
    "DL2_16_SINGLE_MODEL_REFERENCE_POLICY_FIXED": (
        "LIONESS reference cohort MUST follow predeclared fixed policy."
    ),
    "DL2_17_SAMPLE_ORDER_DETERMINISTIC": (
        "Expression sample order MUST be deterministic and saved to sidecar."
    ),
    "DL2_18_NO_ENUMERATION_LABELING": (
        "Sample labels MUST NOT be assigned by positional enumeration."
    ),
    "DL2_19_TARGET_EXPORT_ONLY_IF_VERIFIED": (
        "Target-specific LIONESS network exported only if flatten order verified."
    ),
    "DL2_20_RUNTIME_PRIOR_TRANSFORM_MINIMAL": (
        "Runtime prior canonicalization must be minimal, logged, and recorded."
    ),
    "DL2_21_GATE_REPORT_MUST_RECORD_REFERENCE_COHORT": (
        "Gate report MUST store reference cohort membership and flatten order."
    ),
    "DL2_22_FROZEN_FOUNDATION_INPUTS_ONLY": (
        "Notebook 02B MUST consume only frozen lineage inputs from 02A."
    ),
    "DL2_23_SPLIT_NOTEBOOK_PROVENANCE": (
        "Notebooks 02A and 02B MUST have separate log files and gate reports."
    ),
    "DL2_24_TARGET_INPUT_DECLARED_IN_CONFIG": (
        "Notebook 02B target MUST be declared via cfg.L2_TARGET_MODEL_INPUT "
        "or cfg.L2_TARGET_MODELS_LIST (multi-model)."
    ),
    # ── Layer 2 multi-model (v2.4.0 / config v3.4.0) ────────
    "DL2_25_SHARED_PANDA_PER_LINEAGE": (
        "In TARGETED_LIONESS_MULTI_MODEL mode, all targets belonging to the "
        "same OncotreeLineage MUST share one PANDA aggregate network e^(α) "
        "built once from the frozen lineage cohort. Rebuilding PANDA per "
        "target inside the same lineage is FORBIDDEN (would reintroduce "
        "numerical non-identity of the aggregate and false biological "
        "differences). Cross-lineage pooling of the aggregate is also FORBIDDEN."
    ),
    "DL2_26_TARGET_LIST_DECLARED_IN_CONFIG": (
        "Multi-model targets MUST come exclusively from "
        "cfg.L2_TARGET_MODELS_LIST. Runtime discovery, ad-hoc CLI injection "
        "without config update, or silent list mutation is FORBIDDEN."
    ),
    "DL2_27_EXPORT_ONLY_IF_VERIFIED_MULTI": (
        "Each per-target LIONESS TSV may be exported only after the lineage "
        "PANDA flatten order (C or F) has been verified against "
        "export_panda_results for that shared PANDA object."
    ),
    "DL2_28_FULL_LEDGER_FOR_DECLARED_LIST": (
        "Every ModelID in cfg.L2_TARGET_MODELS_LIST MUST appear exactly once "
        "in the multi-model run ledger with status COMPLETED, SKIPPED, or "
        "FAILED and an explicit reason. Silent drops are FORBIDDEN."
    ),
    # ── Layer 3 (v3.0.0) ────────────────────────────────────
    "DL3_01_P0_NO_REWEIGHT": (
        "Layer 3 MUST use Layer 1 core seed as-is (uniform 1/|T_d|). "
        "No re-weighting by docking score, CNN_VS, or expression downstream."
    ),
    "DL3_02_GRN_CLIP_NONPOSITIVE": (
        f"GRN edges with weight <= 0 MUST be removed before RWR. "
        f"Policy: {L3_GRN_WEIGHT_TRANSFORM}."
    ),
    "DL3_03_FROZEN_LAYER2_OUTPUTS": (
        "Layer 3 MUST consume only frozen Layer 2 outputs (PANDA/LIONESS TSV). "
        "No GRN recomputation in Layer 3 notebook."
    ),
    "DL3_04_DRUG_MONOPLEX_POLICY": (
        f"Drug monoplex policy: {L3_DRUG_MONOPLEX_POLICY}. "
        "No drug-drug similarity edges permitted."
    ),
    "DL3_05_NO_AGGREGATION_IN_LAYER3": (
        "Layer 3 MUST preserve raw per-layer RWR scores. Any global gene ranking "
        "computed downstream MUST remain explicit, provenance-tracked, and must NOT "
        "overwrite the per-layer outputs."
    ),
    # ── Layer 3 (v3.1.0 / v3.3.0) ───────────────────────────
    "DL3_06_LAYER3_PPI_MUST_USE_GENE_LEVEL_STRING": (
        "Layer 3 MUST use LOCAL_PPI_PRIOR_L3_TSV built from STRING v12.0 "
        "physical links only, score>=700, gene-level (NOT TF-restricted), "
        "3-column weighted format. "
        "Layer 2 TF-only PPI prior (LOCAL_PPI_PRIOR_TXT) is FORBIDDEN for Layer 3 RWR. "
        "Rationale: TF-only PPI excludes kinases, receptors, enzymes that are "
        "valid docking targets (EGFR, ERBB2, KDR, PTGS2, PTGES, ALOX5)."
    ),
    "DL3_07_SINGLE_REPRESENTATIVE_RUN": (
        "If all compounds are mathematically equivalent for RWR under the current "
        "implementation (same valid seed gene list after universe filtering, shared PPI, "
        "shared GRN, shared hyperparameters, and no downstream seed reweighting), Layer 3 "
        "MUST execute exactly one synthetic CLASS_SHARED run instead of redundant per-compound runs. "
        "This policy prevents computational waste and avoids implying compound-specific variation "
        "that does not exist."
    ),
    # ── Layer 4 (v4.0.0) ────────────────────────────────────
    "DL4_01_EXPLICIT_BACKGROUND": (
        "Layer 4 ORA MUST use background=gene_universe (explicit Layer 3 "
        "gene universe artifact). background=None is FORBIDDEN. "
        "Confirmed: background=None uses N=17071 (GMT union), "
        "background=gene_universe uses N=10286 (correct)."
    ),
    "DL4_02_QUERY_SUBSET_BACKGROUND": (
        "All query genes MUST exist in background. GSEApy silently drops "
        "query genes outside background without warning. Hard fail required."
    ),
    "DL4_03_PRESERVE_RANKING_ORDER": (
        "Query gene list MUST preserve Layer 3 ranking order from "
        "top100_noHK.tsv. Alphabetical sorting is FORBIDDEN for query."
    ),
    "DL4_04_SANITY_CHECK_FORMULA": (
        "Sanity check hypergeometric must use: "
        "N=len(background), n=len(query∩background), "
        "M=len(pathway∩background), k=len(query∩pathway∩background). "
        "Using full GMT pathway size for M is WRONG."
    ),
    "DL4_05_REDUNDANCY_SAME_UNIVERSE": (
        "Redundancy filter MUST use pathway gene sets restricted to the "
        "same Layer 3 gene universe as ORA. Full GMT raw sets are FORBIDDEN "
        "for pairwise Jaccard/hypergeometric/submodule calculations."
    ),
    "DL4_06_REDUNDANCY_ON_SIGNIFICANT": (
        "Redundancy filter MUST run on SIGNIFICANT results only, "
        "not on ALL_results."
    ),
}


# ============================================================
# UTILITY FUNCTIONS
# ============================================================


def normalize_gene_name(name: str) -> str:
    """Remove hyphens, uppercase. 'ALOX-5' → 'ALOX5'"""
    if name is None:
        return ""
    return str(name).replace("-", "").strip().upper()


def combine_edge_weights_bayes(w_a: float, w_b: float) -> float:
    """P(A∪B) = 1-(1-P(A))(1-P(B)). For downstream multiplex layers."""
    return 1.0 - (1.0 - float(w_a)) * (1.0 - float(w_b))


def strip_entrez_ids(columns) -> list:
    """'TP53 (7157)' → 'TP53'. Enforces DL2-12."""
    return [str(c).split(" (")[0].strip() for c in columns]


def verify_file_sha256(path: Path, expected_sha256: str) -> bool:
    """Verify SHA-256 hash of a file for provenance validation."""
    import hashlib

    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(65536), b""):
            h.update(chunk)
    actual = h.hexdigest()
    return actual == expected_sha256


def resolve_model_id(raw_input: str) -> str:
    """
    Resolve a user input to a DepMap ModelID.
    Accepts alias from CELL_LINE_MODELS or a raw ACH-* ModelID.
    """
    s = str(raw_input).strip()
    return CELL_LINE_MODELS.get(s, s)


def validate_target_models_list(models_list=None) -> list:
    """
    Validate L2_TARGET_MODELS_LIST integrity.
    Returns cleaned list of unique ModelIDs (order preserved).
    Hard-fails on empty list, wrong expected count, or duplicates
    when L2_MULTI_MODEL_DEDUPLICATE_LIST is False.
    """
    if models_list is None:
        models_list = list(L2_TARGET_MODELS_LIST)
    cleaned = [str(x).strip() for x in models_list if str(x).strip()]
    if not cleaned:
        raise ValueError(
            "DEADLOCK DL2-26: L2_TARGET_MODELS_LIST is empty."
        )
    if L2_TARGET_MODELS_LIST_EXPECTED_N is not None:
        if len(cleaned) != int(L2_TARGET_MODELS_LIST_EXPECTED_N):
            raise ValueError(
                f"DEADLOCK DL2-26: L2_TARGET_MODELS_LIST has {len(cleaned)} "
                f"entries, expected {L2_TARGET_MODELS_LIST_EXPECTED_N}."
            )
    seen = set()
    dups = []
    ordered_unique = []
    for m in cleaned:
        if m in seen:
            dups.append(m)
        else:
            seen.add(m)
            ordered_unique.append(m)
    if dups and not L2_MULTI_MODEL_DEDUPLICATE_LIST:
        raise ValueError(
            f"DEADLOCK DL2-26: duplicate ModelIDs in L2_TARGET_MODELS_LIST: "
            f"{sorted(set(dups))}"
        )
    if L2_MULTI_MODEL_ALLOW_CROSS_LINEAGE_POOLING:
        raise ValueError(
            "DEADLOCK DL2-25: L2_MULTI_MODEL_ALLOW_CROSS_LINEAGE_POOLING "
            "must be False."
        )
    return ordered_unique if L2_MULTI_MODEL_DEDUPLICATE_LIST else cleaned


# ============================================================
# DEADLOCK VALIDATION ENGINE
# ============================================================


def validate_deadlock_rules(step: str, context=None, **kwargs) -> None:
    """Validate deadlock rules at pipeline checkpoints."""
    if context is None:
        context = {}
    # Allow both context dict and kwargs (notebook convenience)
    if kwargs:
        context = {**context, **kwargs}

    # ── Layer 1 ─────────────────────────────────────────────
    if step == "pre_normalize":
        assert context.get("filtered", False), (
            "DEADLOCK DL1: normalize/filter ordering violation."
        )
    elif step == "dl1c_no_runtime_discovery":
        assert context.get("using_frozen_artifact", False), (
            f"DEADLOCK DL1C: {DEADLOCK_RULES['DL1C_NO_RUNTIME_TARGET_DISCOVERY']}"
        )
    elif step == "dl1c_no_inline_list":
        assert context.get("read_from_file", False), (
            f"DEADLOCK DL1C: {DEADLOCK_RULES['DL1C_NO_INLINE_GENE_LIST']}"
        )
    elif step == "dl19c_unweighted_seed":
        assert CORE_SEED_UNWEIGHTED, (
            "DEADLOCK DL19C: CORE_SEED_UNWEIGHTED must be True."
        )
        assert L1_SEED_MODE == "CONSENSUS_GLOBAL_UNWEIGHTED", (
            f"DEADLOCK DL19C: L1_SEED_MODE='{L1_SEED_MODE}' "
            "must be 'CONSENSUS_GLOBAL_UNWEIGHTED'."
        )
    elif step == "dl20c_support_rule_fixed":
        manifest_rule = context.get("manifest_support_rule", -1)
        assert manifest_rule >= L1_MIN_SUPPORT_METHODS, (
            f"DEADLOCK DL20C: manifest declares support rule "
            f"{manifest_rule} but config requires >= {L1_MIN_SUPPORT_METHODS}."
        )
    elif step == "dl21c_no_downstream_reweight":
        assert not context.get("p0_modified_downstream", False), (
            f"DEADLOCK DL21C: {DEADLOCK_RULES['DL21C_NO_DOWNSTREAM_REWEIGHTING']}"
        )
    elif step == "dl22c_seed_collapse_acknowledged":
        assert context.get("seed_collapse_acknowledged", False), (
            f"DEADLOCK DL22C: {DEADLOCK_RULES['DL22C_SEED_COLLAPSE_ACKNOWLEDGED']}"
        )

    # ── Layer 2 ─────────────────────────────────────────────
    elif step == "dl2_01_scenic_forbidden":
        assert SCENIC_FORBIDDEN_ON_BULK
        assert not context.get("scenic_used", False)
    elif step == "dl2_02_min_samples":
        lineage = context.get("lineage", "UNKNOWN")
        n_samples = context.get("n_samples", 0)
        assert n_samples >= MIN_SAMPLES_PANDA, (
            f"DEADLOCK DL2-02: '{lineage}' has {n_samples} < {MIN_SAMPLES_PANDA}"
        )
    elif step == "dl2_03_no_global_qnorm":
        assert GLOBAL_QNORM_FORBIDDEN
        assert not context.get("global_qnorm_applied", False)
    elif step == "dl2_04_snail_applied":
        lineage = context.get("lineage", "UNKNOWN")
        assert SNAIL_NORMALIZATION_REQUIRED
        assert context.get("snail_applied", False), (
            f"DEADLOCK DL2-04: SNAIL not applied for '{lineage}'"
        )
    elif step == "dl2_05_prior_source":
        psr = context.get("prior_source_runtime", "UNKNOWN")
        assert LOCAL_FROZEN_PRIORS_ONLY
        assert psr == "LOCAL_FILES", (
            f"DEADLOCK DL2-05: prior_source_runtime='{psr}'"
        )
    elif step == "dl2_06_prior_cutoffs":
        sc = context.get("string_cutoff", None)
        if sc is not None:
            assert sc >= PPI_PRIOR_SCORE_CUTOFF
    elif step == "dl2_07_grn_builder":
        builder = context.get("builder_used", "")
        assert builder == GRN_BUILDER
    elif step == "dl2_09_lioness_after_panda":
        assert context.get("panda_converged", False)
    elif step == "dl2_10_lioness_same_priors":
        pm = context.get("panda_motif_path", "")
        lm = context.get("lioness_motif_path", "")
        pp = context.get("panda_ppi_path", "")
        lp = context.get("lioness_ppi_path", "")
        assert pm == lm, "DEADLOCK DL2-10: motif mismatch"
        assert pp == lp, "DEADLOCK DL2-10: ppi mismatch"
    elif step == "dl2_12_gene_ids":
        assert context.get("entrez_stripped", False)
    elif step == "dl2_14_target_model":
        assert context.get("target_exists_uniquely", False)
    elif step == "dl2_15_lineage_from_metadata":
        assert context.get("lineage_from_metadata", False)
    elif step == "dl2_17_sample_order_deterministic":
        assert context.get("sidecar_written", False)
    elif step == "dl2_19_export_only_if_verified":
        assert context.get("flatten_order_verified", False)
    elif step == "dl2_22_frozen_foundation_inputs":
        assert context.get("using_frozen_foundation_inputs", False)
        assert not context.get("snail_recomputed", False)
    elif step == "dl2_23_split_provenance":
        assert context.get("separate_logs", False)
        assert context.get("separate_gate_reports", False)
    elif step == "dl2_24_target_input_declared":
        # Single-model: non-empty string; multi-model: non-empty list also OK
        single = str(context.get("target_model_input", "")).strip()
        multi = context.get("target_models_list", None)
        ok = bool(single) or (
            isinstance(multi, (list, tuple)) and len(multi) > 0
        )
        assert ok, (
            "DEADLOCK DL2-24: neither L2_TARGET_MODEL_INPUT nor "
            "L2_TARGET_MODELS_LIST declared."
        )

    # ── Layer 2 multi-model ─────────────────────────────────
    elif step == "dl2_25_shared_panda_per_lineage":
        assert L2_MULTI_MODEL_REQUIRE_SHARED_PANDA_PER_LINEAGE is True
        assert L2_MULTI_MODEL_ALLOW_CROSS_LINEAGE_POOLING is False
        assert context.get("shared_panda_per_lineage", False), (
            f"DEADLOCK DL2-25: {DEADLOCK_RULES['DL2_25_SHARED_PANDA_PER_LINEAGE']}"
        )
        assert not context.get("cross_lineage_pooling", False), (
            "DEADLOCK DL2-25: cross-lineage pooling detected."
        )
        n_panda = int(context.get("n_panda_builds_for_lineage", -1))
        assert n_panda == 1, (
            f"DEADLOCK DL2-25: expected 1 PANDA build per lineage, got {n_panda}."
        )
    elif step == "dl2_26_target_list_declared":
        assert LAYER2_RUN_MODE == "TARGETED_LIONESS_MULTI_MODEL"
        src = context.get("target_list_source", "")
        assert src == "cfg.L2_TARGET_MODELS_LIST", (
            f"DEADLOCK DL2-26: target_list_source='{src}'"
        )
        assert context.get("list_validated", False), (
            "DEADLOCK DL2-26: list not validated via validate_target_models_list()."
        )
    elif step == "dl2_27_export_only_if_verified_multi":
        assert context.get("flatten_order_verified", False), (
            f"DEADLOCK DL2-27: {DEADLOCK_RULES['DL2_27_EXPORT_ONLY_IF_VERIFIED_MULTI']}"
        )
    elif step == "dl2_28_full_ledger":
        declared = list(context.get("declared_ids", []))
        ledger_ids = list(context.get("ledger_ids", []))
        assert len(declared) > 0
        assert sorted(declared) == sorted(ledger_ids), (
            f"DEADLOCK DL2-28: ledger IDs mismatch declared list. "
            f"missing={sorted(set(declared) - set(ledger_ids))} "
            f"extra={sorted(set(ledger_ids) - set(declared))}"
        )
        assert context.get("no_silent_drops", False), (
            "DEADLOCK DL2-28: silent drops detected."
        )

    # ── Layer 3 (DL3_01–DL3_05) ─────────────────────────────
    elif step == "dl3_01_p0_no_reweight":
        assert not context.get("reweighted", False), (
            f"DEADLOCK DL3-01: {DEADLOCK_RULES['DL3_01_P0_NO_REWEIGHT']}"
        )
    elif step == "dl3_02_grn_clip_applied":
        prov = context.get("provenance", {})
        transform = prov.get("grn_weight_transform", {})
        assert transform.get("transform_method") == L3_GRN_WEIGHT_TRANSFORM, (
            f"DEADLOCK DL3-02: transform_method="
            f"'{transform.get('transform_method')}' "
            f"!= '{L3_GRN_WEIGHT_TRANSFORM}'"
        )
        assert transform.get("n_edges_after", -1) >= 0
    elif step == "dl3_03_frozen_layer2":
        assert not context.get("grn_recomputed", False)
        assert context.get("using_frozen_layer2", False)
    elif step == "dl3_04_drug_monoplex":
        policy = context.get("drug_monoplex_policy", "")
        assert policy == L3_DRUG_MONOPLEX_POLICY, (
            f"DEADLOCK DL3-04: policy='{policy}' != '{L3_DRUG_MONOPLEX_POLICY}'"
        )
    elif step == "dl3_05_no_aggregation":
        gene_df = context.get("gene_df", None)
        assert gene_df is not None
        assert "layer" in gene_df.columns, (
            "DEADLOCK DL3-05: 'layer' column missing — aggregation may have occurred."
        )

    # ── Layer 3 (DL3_06, v3.1.0) ────────────────────────────
    elif step == "dl3_06_layer3_ppi_source":
        ppi_path_used = context.get("ppi_path_used", "")
        assert str(ppi_path_used) == str(LOCAL_PPI_PRIOR_L3_TSV), (
            f"DEADLOCK DL3-06: ppi_path_used='{ppi_path_used}' "
            f"!= '{LOCAL_PPI_PRIOR_L3_TSV}'"
        )
        assert context.get("manifest_validated", False), (
            "DEADLOCK DL3-06: PPI L3 manifest not validated."
        )
        tf_restricted = context.get("tf_restricted", True)
        assert tf_restricted is False, (
            f"DEADLOCK DL3-06: tf_restricted={tf_restricted} must be False."
        )
        sha_ok = context.get("sha256_ok", False)
        assert sha_ok, (
            "DEADLOCK DL3-06: SHA-256 of PPI L3 artifact does not match expected. "
            "File may have been modified."
        )
    elif step == "dl3_07_single_representative_run":
        assert context.get("mathematically_equivalent", False), (
            f"DEADLOCK DL3-07: {DEADLOCK_RULES['DL3_07_SINGLE_REPRESENTATIVE_RUN']}"
        )
        assert context.get("synthetic_class_node", False), (
            "DEADLOCK DL3-07: synthetic CLASS_SHARED node required."
        )
        assert int(context.get("n_runs_executed", -1)) == 1, (
            "DEADLOCK DL3-07: exactly one RWR run must be executed."
        )

    # ── Layer 4 (DL4_01–DL4_06) ─────────────────────────────
    elif step == "dl4_01_explicit_background":
        bg_policy = context.get("background_policy", "")
        assert bg_policy == L4_BACKGROUND_POLICY, (
            f"DEADLOCK DL4-01: background_policy='{bg_policy}' "
            f"!= '{L4_BACKGROUND_POLICY}'"
        )
        assert context.get("background_is_none") is False, (
            "DEADLOCK DL4-01: background=None is FORBIDDEN."
        )
    elif step == "dl4_02_query_subset_background":
        outside = context.get("genes_outside_background", set())
        assert len(outside) == 0, (
            f"DEADLOCK DL4-02: {len(outside)} query genes not in background: "
            f"{sorted(list(outside))[:10]}"
        )
    elif step == "dl4_03_preserve_ranking_order":
        assert not context.get("query_was_sorted", False), (
            "DEADLOCK DL4-03: Query was alphabetically sorted. "
            "This destroys Layer 3 ranking."
        )
    elif step == "dl4_04_sanity_check_formula":
        assert context.get("M_is_background_restricted", False), (
            "DEADLOCK DL4-04: M must be pathway∩background, not full GMT size."
        )
    elif step == "dl4_05_redundancy_same_universe":
        assert context.get("redundancy_restricted_to_background", False), (
            "DEADLOCK DL4-05: Redundancy filter must use background-restricted "
            "pathway gene sets."
        )
    elif step == "dl4_06_redundancy_on_significant":
        assert context.get("input_is_significant_only", False), (
            "DEADLOCK DL4-06: Redundancy filter must run on SIGNIFICANT only."
        )
    else:
        import warnings as _w

        _w.warn(
            f"validate_deadlock_rules: unknown step '{step}'.",
            stacklevel=2,
        )


# ============================================================
# RESOURCE CHECKS
# ============================================================


def check_layer1_resources() -> dict:
    """Check Layer 1 (Consensus Prior) runtime dependencies."""
    results: dict = {}
    results["consensus_target_list_exists"] = (
        LOCAL_CONSENSUS_TARGET_LIST_TXT.exists()
    )
    results["consensus_ligand_ids_exists"] = (
        LOCAL_CONSENSUS_LIGAND_IDS_TXT.exists()
    )
    results["consensus_manifest_exists"] = (
        LOCAL_CONSENSUS_MANIFEST_JSON.exists()
    )
    if LOCAL_CONSENSUS_TARGET_LIST_EXPECTED_SHA256:
        results["consensus_target_list_sha256_ok"] = verify_file_sha256(
            LOCAL_CONSENSUS_TARGET_LIST_TXT,
            LOCAL_CONSENSUS_TARGET_LIST_EXPECTED_SHA256,
        )
    else:
        results["consensus_target_list_sha256_ok"] = None
    if LOCAL_CONSENSUS_LIGAND_IDS_EXPECTED_SHA256:
        results["consensus_ligand_ids_sha256_ok"] = verify_file_sha256(
            LOCAL_CONSENSUS_LIGAND_IDS_TXT,
            LOCAL_CONSENSUS_LIGAND_IDS_EXPECTED_SHA256,
        )
    else:
        results["consensus_ligand_ids_sha256_ok"] = None
    LAYER1_OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    results["layer1_output_dir"] = LAYER1_OUTPUT_DIR.exists()
    LOG_DIR.mkdir(parents=True, exist_ok=True)
    results["layer1_log_dir"] = LOG_DIR.exists()
    return results


def check_layer2_grn_resources() -> dict:
    results: dict = {}
    results["ccle_model_csv"] = CCLE_MODEL_CSV.exists()
    results["ccle_tpm_expression_csv"] = CCLE_TPM_EXPRESSION_CSV.exists()
    results["local_motif_prior_500bp"] = LOCAL_MOTIF_PRIOR_TXT.exists()
    results["local_motif_prior_10kb"] = LOCAL_MOTIF_PRIOR_10KB_TXT.exists()
    results["local_ppi_prior"] = LOCAL_PPI_PRIOR_TXT.exists()
    results["both_core_priors_exist"] = (
        LOCAL_MOTIF_PRIOR_TXT.exists() and LOCAL_PPI_PRIOR_TXT.exists()
    )
    results["motif_prior_manifest_500bp"] = (
        LOCAL_MOTIF_PRIOR_MANIFEST_JSON.exists()
    )
    results["ppi_prior_manifest"] = LOCAL_PPI_PRIOR_MANIFEST_JSON.exists()
    results["layer2_grn_data_dir"] = LAYER2_GRN_DATA_DIR.exists()
    results["layer2_grn_output_dir"] = LAYER2_GRN_OUTPUT_DIR.exists()
    try:
        import netZooPy as _nz

        results["netzoopy_available"] = True
        results["netzoopy_version"] = getattr(_nz, "__version__", "unknown")
    except ImportError:
        results["netzoopy_available"] = False
        results["netzoopy_version"] = "NOT_INSTALLED"
    try:
        from pysnail import qsmooth as _qs  # noqa: F401

        results["pysnail_available"] = True
    except ImportError:
        results["pysnail_available"] = False
    results["l1_p0_vector_long"] = (
        LAYER1_OUTPUT_DIR / L1_P0_VECTOR_LONG_CSV
    ).exists()
    results["l2_target_models_list_n"] = len(L2_TARGET_MODELS_LIST)
    results["l2_run_mode"] = LAYER2_RUN_MODE
    return results


def check_layer2b_resources() -> dict:
    results: dict = {}
    results["ccle_model_csv"] = CCLE_MODEL_CSV.exists()
    results["local_motif_prior_500bp"] = LOCAL_MOTIF_PRIOR_TXT.exists()
    results["local_ppi_prior"] = LOCAL_PPI_PRIOR_TXT.exists()
    results["frozen_inputs_dir"] = L2A_FROZEN_INPUTS_DIR.exists()
    _manifest = L2A_FROZEN_INPUTS_DIR / L2A_FOUNDATION_MANIFEST_JSON
    results["foundation_manifest"] = _manifest.exists()
    if results["foundation_manifest"]:
        try:
            import json as _json

            with open(_manifest, "r", encoding="utf-8") as _fh:
                _fm = _json.load(_fh)
            _frozen = _fm.get("lineages_frozen", [])
            results["lineages_frozen"] = _frozen
            for _l in _frozen:
                _arts = _fm.get("lineage_artifacts", {}).get(_l, {})
                _ep = LAYER2_GRN_OUTPUT_DIR / _arts.get("expr_parquet", "")
                _sp = LAYER2_GRN_OUTPUT_DIR / _arts.get("sample_order_json", "")
                _gp = LAYER2_GRN_OUTPUT_DIR / _arts.get("genes_post_qc_json", "")
                results[f"{_l}_expr_parquet"] = _ep.exists()
                results[f"{_l}_sample_order"] = _sp.exists()
                results[f"{_l}_genes_post_qc"] = _gp.exists()
        except Exception as _e:
            results["foundation_manifest_readable"] = False
            results["foundation_manifest_error"] = str(_e)
    else:
        results["lineages_frozen"] = []
    try:
        import netZooPy as _nz

        results["netzoopy_available"] = True
        results["netzoopy_version"] = getattr(_nz, "__version__", "unknown")
    except ImportError:
        results["netzoopy_available"] = False
        results["netzoopy_version"] = "NOT_INSTALLED"
    return results


def check_layer3_resources() -> dict:
    """Check Layer 3 runtime dependencies and upstream outputs."""
    results: dict = {}
    try:
        import multixrank as _mx

        results["multixrank_available"] = True
        _ver = getattr(_mx, "__version__", "unknown")
        results["multixrank_version"] = _ver
        results["multixrank_version_ok"] = (
            _ver == L3_MULTIXRANK_REQUIRED_VERSION
        )
    except ImportError:
        results["multixrank_available"] = False
        results["multixrank_version"] = "NOT_INSTALLED"
        results["multixrank_version_ok"] = False
    try:
        import yaml  # noqa: F401

        results["yaml_available"] = True
    except ImportError:
        results["yaml_available"] = False
    results["layer1_output_dir"] = LAYER1_OUTPUT_DIR.exists()
    results["layer1_p0_long_csv"] = (
        LAYER1_OUTPUT_DIR / L1_P0_VECTOR_LONG_CSV
    ).exists()
    results["layer2_grn_output_dir"] = LAYER2_GRN_OUTPUT_DIR.exists()
    results["l3_string_physical_links"] = L3_STRING_PHYSICAL_LINKS.exists()
    results["l3_string_info"] = L3_STRING_INFO.exists()
    results["local_ppi_prior_l3_tsv"] = LOCAL_PPI_PRIOR_L3_TSV.exists()
    results["local_ppi_prior_l3_manifest"] = (
        LOCAL_PPI_PRIOR_L3_MANIFEST_JSON.exists()
    )
    if results["local_ppi_prior_l3_tsv"]:
        sha_ok = verify_file_sha256(
            LOCAL_PPI_PRIOR_L3_TSV,
            LOCAL_PPI_PRIOR_L3_EXPECTED_SHA256,
        )
        results["local_ppi_prior_l3_sha256_ok"] = sha_ok
    else:
        results["local_ppi_prior_l3_sha256_ok"] = False
    LAYER3_OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    results["layer3_output_dir"] = LAYER3_OUTPUT_DIR.exists()
    LOG_DIR.mkdir(parents=True, exist_ok=True)
    results["layer3_log_dir"] = LOG_DIR.exists()
    return results


def check_layer4_resources(class_shared_dir: Path = None) -> dict:
    """Check Layer 4 runtime dependencies and upstream artifacts."""
    results: dict = {}
    try:
        import gseapy as _gp

        results["gseapy_available"] = True
        _ver = getattr(_gp, "__version__", "unknown")
        results["gseapy_version"] = _ver
        results["gseapy_version_ok"] = (_ver == L4_GSEAPY_REQUIRED_VERSION)
    except ImportError:
        results["gseapy_available"] = False
        results["gseapy_version"] = "NOT_INSTALLED"
        results["gseapy_version_ok"] = False
    try:
        import scipy as _sp

        results["scipy_available"] = True
        results["scipy_version"] = _sp.__version__
    except ImportError:
        results["scipy_available"] = False
    try:
        import sklearn as _sk

        results["sklearn_available"] = True
        results["sklearn_version"] = _sk.__version__
    except ImportError:
        results["sklearn_available"] = False
    results["gmt_file_exists"] = L4_GMT_PATH.exists()
    if class_shared_dir is not None:
        _csd = Path(class_shared_dir)
        results["class_shared_dir_exists"] = _csd.exists()
        _q = sorted(_csd.glob(L4_QUERY_FILENAME_PATTERN))
        results["query_file_found"] = len(_q) == 1
        if _q:
            results["query_file_name"] = _q[0].name
        _bg_txt = sorted(_csd.glob(L4_BACKGROUND_TXT_PATTERN))
        results["background_txt_found"] = len(_bg_txt) == 1
        if _bg_txt:
            results["background_txt_name"] = _bg_txt[0].name
        _bg_manifest = sorted(_csd.glob(L4_BACKGROUND_MANIFEST_PATTERN))
        results["background_manifest_found"] = len(_bg_manifest) == 1
    else:
        results["class_shared_dir_exists"] = None
        results["query_file_found"] = None
        results["background_txt_found"] = None
    LAYER4_OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    results["layer4_output_dir"] = LAYER4_OUTPUT_DIR.exists()
    return results


# ============================================================
# CONFIG SUMMARY PRINTER
# ============================================================


def print_config_summary() -> None:
    w = 66
    print("=" * w)
    print("CONFIG SYSTEM v3.5.0 — Layer 1 + Layer 2 + Layer 3 + Layer 4")
    print("=" * w)
    print(f"Project root          : {PROJECT_ROOT}")
    print(f"CCLE data dir         : {CCLE_DATA_DIR}")
    print(f"Layer 1 output        : {LAYER1_OUTPUT_DIR}")
    print(f"Layer 2 GRN output    : {LAYER2_GRN_OUTPUT_DIR}")
    print(f"Layer 3 output        : {LAYER3_OUTPUT_DIR}")
    print("-" * w)
    print("LAYER 1 — Thermodynamic Gate")
    print(f"  Seed mode           : UNWEIGHTED (p0 = 1/|T_d|)")
    print(f"  Deadlock            : DL19, DL20, DL21, DL22")
    print("-" * w)
    print("LAYER 2 — PANDA/LIONESS (v2.4.0 multi-model)")
    print(f"  Architecture        : {LAYER2_NOTEBOOK_ARCHITECTURE}")
    print(f"  Run mode            : {LAYER2_RUN_MODE}")
    print(f"  Target list n       : {len(L2_TARGET_MODELS_LIST)}")
    print(f"  Compat single alias : {L2_TARGET_MODEL_INPUT}")
    print(f"  Multi missing policy: {L2_MULTI_MODEL_ON_MISSING}")
    print(f"  Shared PANDA/lineage: {L2_MULTI_MODEL_REQUIRE_SHARED_PANDA_PER_LINEAGE}")
    print(f"  Cross-lineage pool  : {L2_MULTI_MODEL_ALLOW_CROSS_LINEAGE_POOLING}")
    print(f"  netZooPy required   : {NETZOOPY_REQUIRED_VERSION}")
    print(f"  Prior policy        : {PRIOR_SOURCE_POLICY}")
    print(f"  PPI prior (L2)      : {LOCAL_PPI_PRIOR_TXT.name} [TF-only, 2-col]")
    print(f"  Motif prior (core)  : {LOCAL_MOTIF_PRIOR_TXT.name}")
    print(f"  Frozen inputs dir   : {L2A_FROZEN_INPUTS_DIR}")
    print(f"  Log 02A             : {LAYER2A_LOG_FILE.name}")
    print(f"  Log 02B single      : {LAYER2B_LOG_FILE.name}")
    print(f"  Log 02B multi       : {LAYER2B_MULTI_LOG_FILE.name}")
    print(f"  Gate report 02A     : {L2A_FOUNDATION_GATE_REPORT_JSON}")
    print(f"  Gate report 02B     : {L2B_LIONESS_GATE_REPORT_JSON}")
    print(f"  Gate report multi   : {L2_MULTI_MODEL_GATE_REPORT_JSON}")
    print("-" * w)
    print("LAYER 3 — RWR-MultiXrank (v3.3.0)")
    print(f"  MultiXrank required : {L3_MULTIXRANK_REQUIRED_VERSION}")
    print(f"  r (restart prob)    : {L3_RWR_RESTART_PROB}")
    print(f"  delta_gene          : {L3_DELTA_GENE}")
    print(f"  delta_drug          : {L3_DELTA_DRUG}")
    print(f"  lambda_default      : {L3_LAMBDA_DEFAULT}")
    print(f"  lambda_bench        : {L3_LAMBDA_BENCH_VALUES}")
    print(f"  eta_default         : {L3_ETA_DEFAULT}")
    print(f"  tau_gene_default    : {L3_TAU_GENE_DEFAULT}")
    print(f"  tau_drug_default    : {L3_TAU_DRUG_DEFAULT}")
    print(f"  graph_type PPI      : {L3_GRAPH_TYPE_WEIGHTED_UNDIRECTED}")
    print(f"  graph_type GRN      : {L3_GRAPH_TYPE_WEIGHTED_DIRECTED}")
    print(f"  graph_type bipartite: {L3_GRAPH_TYPE_UNWEIGHTED_UNDIRECTED}")
    print(f"  GRN transform       : {L3_GRN_WEIGHT_TRANSFORM}")
    print(f"  Drug monoplex       : {L3_DRUG_MONOPLEX_POLICY}")
    print(f"  Top-N preview       : {L3_TOP_N_PREVIEW}")
    print(f"  Global gene score   : {L3_GLOBAL_GENE_SCORE_METHOD}")
    print(f"  CLASS_SHARED prefix : {L3_CLASS_SHARED_NODE_PREFIX}")
    print(f"  Single-run (DL3_07) : {L3_USE_SINGLE_REPRESENTATIVE_RUN}")
    print(f"  Science version     : {L3_SCIENCE_VERSION}")
    print(f"  Architecture        : {L3_ARCHITECTURE_VERSION}")
    print(f"  Batch source        : {L3_BATCH_SOURCE}")
    print(f"  Batch lineage filter: {L3_BATCH_LINEAGE_FILTER}")
    print(f"  Batch notebook      : {LAYER3_NOTEBOOK_BATCH_NAME}")
    print(f"  Jaccard threshold   : {L3_JACCARD_THRESHOLD}")
    print(f"  Spearman threshold  : {L3_SPEARMAN_THRESHOLD}")
    print(f"  Log                 : {LAYER3_LOG_FILE.name}")
    print(f"  Batch log           : {LAYER3_BATCH_LOG_FILE.name}")
    print("-" * w)
    print("LAYER 3 — PPI ARTIFACT (v3.1.0)")
    print(f"  Source policy       : {L3_PPI_SOURCE_POLICY}")
    print(f"  STRING version      : {L3_PPI_STRING_VERSION}")
    print(f"  Network type        : {L3_PPI_NETWORK_TYPE}")
    print(f"  Score cutoff        : {L3_PPI_MIN_STRING_SCORE}")
    print(f"  TF restricted       : {L3_PPI_TF_RESTRICTED}")
    print(f"  Gene mapping        : {L3_PPI_GENE_MAPPING}")
    print(f"  Collapse rule       : {L3_PPI_COLLAPSE_RULE}")
    print(f"  Required columns    : {L3_PPI_REQUIRED_COLUMNS}")
    print(f"  Artifact TSV        : {LOCAL_PPI_PRIOR_L3_TSV.name}")
    print(f"  Artifact manifest   : {LOCAL_PPI_PRIOR_L3_MANIFEST_JSON.name}")
    print(f"  Expected SHA-256    : {LOCAL_PPI_PRIOR_L3_EXPECTED_SHA256[:16]}...")
    print(f"  Raw STRING source   : {L3_STRING_PHYSICAL_LINKS.name}")
    print(f"  Deadlock            : DL3_01..DL3_07")
    print("-" * w)
    print(f"Total deadlock rules  : {len(DEADLOCK_RULES)}")
    print("=" * w)
    print("-" * w)
    print("LAYER 4 — ORA + Redundancy Filter (v4.0.0)")
    print(f"  Background policy   : {L4_BACKGROUND_POLICY}")
    print(f"  Query score column  : {L4_QUERY_SCORE_COLUMN}")
    print(f"  GMT library         : {L4_GMT_LABEL}")
    print(f"  GSEApy required     : {L4_GSEAPY_REQUIRED_VERSION}")
    print(f"  FDR cutoff          : {L4_FDR_CUTOFF}")
    print(f"  Min overlap         : {L4_MIN_OVERLAP}")
    print(f"  Redundancy input    : {L4_REDUNDANCY_INPUT_POLICY}")
    print(f"  Restrict to BG      : {L4_REDUNDANCY_RESTRICT_TO_BACKGROUND}")
    print(f"  ElasticNet L1 ratio : {L4_ELASTICNET_L1_RATIO}")
    print(f"  Jaccard QC threshold: {L4_SUBMODULE_J_THRESHOLD}")
    print(f"  Log                 : {L4_LOG_FILE.name}")
    print(f"  Deadlock            : DL4_01..DL4_06")


# ============================================================
# CONFIG INTEGRITY CHECK
# ============================================================


def _check_config_integrity() -> None:
    import warnings as _w

    # Layer 1
    if not CORE_SEED_UNWEIGHTED:
        _w.warn("[CONFIG] CORE_SEED_UNWEIGHTED must be True.", stacklevel=2)

    # Layer 2
    if PANDA_MODE != "union":
        _w.warn(
            f"[CONFIG] PANDA_MODE='{PANDA_MODE}' expected 'union'.",
            stacklevel=2,
        )
    if not LOCAL_FROZEN_PRIORS_ONLY:
        _w.warn(
            "[CONFIG] LOCAL_FROZEN_PRIORS_ONLY must be True.", stacklevel=2
        )
    if not LOCAL_MOTIF_PRIOR_TXT.exists():
        _w.warn(
            f"[CONFIG] Missing motif prior: {LOCAL_MOTIF_PRIOR_TXT}",
            stacklevel=2,
        )
    if not LOCAL_PPI_PRIOR_TXT.exists():
        _w.warn(
            f"[CONFIG] Missing L2 PPI prior: {LOCAL_PPI_PRIOR_TXT}",
            stacklevel=2,
        )
    if LAYER2_RUN_MODE not in {
        "TARGETED_LIONESS_SINGLE_MODEL",
        "TARGETED_LIONESS_MULTI_MODEL",
    }:
        _w.warn(
            f"[CONFIG L2] Unknown LAYER2_RUN_MODE='{LAYER2_RUN_MODE}'.",
            stacklevel=2,
        )
    if L2_MULTI_MODEL_ALLOW_CROSS_LINEAGE_POOLING:
        _w.warn(
            "[CONFIG L2] L2_MULTI_MODEL_ALLOW_CROSS_LINEAGE_POOLING must be False.",
            stacklevel=2,
        )
    if len(L2_TARGET_MODELS_LIST) != L2_TARGET_MODELS_LIST_EXPECTED_N:
        _w.warn(
            f"[CONFIG L2] L2_TARGET_MODELS_LIST n={len(L2_TARGET_MODELS_LIST)} "
            f"!= expected {L2_TARGET_MODELS_LIST_EXPECTED_N}.",
            stacklevel=2,
        )
    if len(L2_TARGET_MODELS_LIST) != len(set(L2_TARGET_MODELS_LIST)):
        _w.warn(
            "[CONFIG L2] L2_TARGET_MODELS_LIST contains duplicates.",
            stacklevel=2,
        )

    # Layer 3 — RWR parameters
    if not (0 < L3_RWR_RESTART_PROB < 1):
        _w.warn(
            f"[CONFIG L3] L3_RWR_RESTART_PROB={L3_RWR_RESTART_PROB} not in (0,1).",
            stacklevel=2,
        )
    if not (0 <= L3_DELTA_GENE <= 1):
        _w.warn(
            f"[CONFIG L3] L3_DELTA_GENE={L3_DELTA_GENE} not in [0,1].",
            stacklevel=2,
        )
    if len(L3_ETA_DEFAULT) != 2:
        _w.warn(
            "[CONFIG L3] L3_ETA_DEFAULT must have 2 elements.", stacklevel=2
        )
    if abs(sum(L3_ETA_DEFAULT) - 1.0) > 1e-9:
        _w.warn(
            f"[CONFIG L3] L3_ETA_DEFAULT={L3_ETA_DEFAULT} should sum to 1.0.",
            stacklevel=2,
        )
    if len(L3_TAU_GENE_DEFAULT) != 2:
        _w.warn(
            "[CONFIG L3] L3_TAU_GENE_DEFAULT must have 2 elements.",
            stacklevel=2,
        )
    if L3_GLOBAL_GENE_SCORE_METHOD != "geometric_mean_across_layers":
        _w.warn(
            f"[CONFIG L3] L3_GLOBAL_GENE_SCORE_METHOD="
            f"'{L3_GLOBAL_GENE_SCORE_METHOD}' should be "
            f"'geometric_mean_across_layers'.",
            stacklevel=2,
        )
    if not L3_USE_SINGLE_REPRESENTATIVE_RUN:
        _w.warn(
            "[CONFIG L3] L3_USE_SINGLE_REPRESENTATIVE_RUN is False — "
            "per-compound fan-out only valid if seeds truly differ (DL3_07).",
            stacklevel=2,
        )
    if L3_BATCH_SOURCE not in {
        "l2_multimodel_ledger",
        "explicit_list",
        "glob",
    }:
        _w.warn(
            f"[CONFIG L3] Unknown L3_BATCH_SOURCE='{L3_BATCH_SOURCE}'.",
            stacklevel=2,
        )
    if L3_FORBID_MAX_SCORE_AS_GLOBAL_RANKING and (
        L3_GLOBAL_GENE_SCORE_METHOD != "geometric_mean_across_layers"
    ):
        _w.warn(
            "[CONFIG L3] max_score global ranking is forbidden under v3.3.0+.",
            stacklevel=2,
        )
    if not L3_HK_GMT_PATH.exists():
        _w.warn(
            f"[CONFIG L3] Missing housekeeping GMT: {L3_HK_GMT_PATH}",
            stacklevel=2,
        )

    # Layer 3 — PPI artifact
    if not LOCAL_PPI_PRIOR_L3_TSV.exists():
        _w.warn(
            f"[CONFIG L3] Missing Layer 3 PPI artifact: {LOCAL_PPI_PRIOR_L3_TSV}. "
            f"Run 01b_build_ppi_prior_l3_from_string.py first.",
            stacklevel=2,
        )
    if not LOCAL_PPI_PRIOR_L3_MANIFEST_JSON.exists():
        _w.warn(
            f"[CONFIG L3] Missing Layer 3 PPI manifest: "
            f"{LOCAL_PPI_PRIOR_L3_MANIFEST_JSON}.",
            stacklevel=2,
        )
    if L3_PPI_TF_RESTRICTED is not False:
        _w.warn(
            "[CONFIG L3] L3_PPI_TF_RESTRICTED must be False for Layer 3.",
            stacklevel=2,
        )
    if LOCAL_PPI_PRIOR_L3_TSV.exists():
        sha_ok = verify_file_sha256(
            LOCAL_PPI_PRIOR_L3_TSV,
            LOCAL_PPI_PRIOR_L3_EXPECTED_SHA256,
        )
        if not sha_ok:
            _w.warn(
                f"[CONFIG L3] SHA-256 mismatch for {LOCAL_PPI_PRIOR_L3_TSV.name}. "
                f"File may have been modified since build.",
                stacklevel=2,
            )
    if L3_GRN_WEIGHT_TRANSFORM != "clip_nonpositive":
        _w.warn(
            f"[CONFIG L3] L3_GRN_WEIGHT_TRANSFORM='{L3_GRN_WEIGHT_TRANSFORM}' "
            f"should be 'clip_nonpositive'.",
            stacklevel=2,
        )

    # Layer 4
    if not L4_GMT_PATH.exists():
        _w.warn(f"[CONFIG L4] Missing GMT: {L4_GMT_PATH}", stacklevel=2)
    if L4_BACKGROUND_POLICY != "explicit_layer3_gene_universe":
        _w.warn(
            "[CONFIG L4] L4_BACKGROUND_POLICY must be "
            "'explicit_layer3_gene_universe'.",
            stacklevel=2,
        )
    if not L4_REDUNDANCY_RESTRICT_TO_BACKGROUND:
        _w.warn(
            "[CONFIG L4] L4_REDUNDANCY_RESTRICT_TO_BACKGROUND must be True.",
            stacklevel=2,
        )
    if "rwr_score_geom_mean" not in L4_QUERY_REQUIRED_COLUMNS:
        _w.warn(
            "[CONFIG L4] L4_QUERY_REQUIRED_COLUMNS must include "
            "'rwr_score_geom_mean'.",
            stacklevel=2,
        )
    if "layer_support_min" not in L4_QUERY_REQUIRED_COLUMNS:
        _w.warn(
            "[CONFIG L4] L4_QUERY_REQUIRED_COLUMNS must include "
            "'layer_support_min'.",
            stacklevel=2,
        )


_check_config_integrity()


# ============================================================
# Self-test
# ============================================================

if __name__ == "__main__":
    import pandas as pd

    print_config_summary()
    print(
        f"\nRunning deadlock rule self-tests ({len(DEADLOCK_RULES)} rules)..."
    )

    # Validate multi-model list itself
    _validated = validate_target_models_list()
    print(
        f"  validate_target_models_list                 → PASS "
        f"(n={len(_validated)})"
    )

    _tests = [
        ("dl19c_unweighted_seed", {}),
        ("dl2_01_scenic_forbidden", {"scenic_used": False}),
        ("dl2_02_min_samples", {"lineage": "Lymphoid", "n_samples": 70}),
        ("dl2_04_snail_applied", {"lineage": "Lymphoid", "snail_applied": True}),
        ("dl2_05_prior_source", {"prior_source_runtime": "LOCAL_FILES"}),
        ("dl2_14_target_model", {"target_exists_uniquely": True}),
        ("dl2_15_lineage_from_metadata", {"lineage_from_metadata": True}),
        ("dl2_17_sample_order_deterministic", {"sidecar_written": True}),
        ("dl2_19_export_only_if_verified", {"flatten_order_verified": True}),
        (
            "dl2_22_frozen_foundation_inputs",
            {
                "using_frozen_foundation_inputs": True,
                "snail_recomputed": False,
            },
        ),
        (
            "dl2_23_split_provenance",
            {
                "separate_logs": True,
                "separate_gate_reports": True,
            },
        ),
        (
            "dl2_24_target_input_declared",
            {
                "target_model_input": L2_TARGET_MODEL_INPUT,
                "target_models_list": L2_TARGET_MODELS_LIST,
            },
        ),
        (
            "dl2_25_shared_panda_per_lineage",
            {
                "shared_panda_per_lineage": True,
                "cross_lineage_pooling": False,
                "n_panda_builds_for_lineage": 1,
            },
        ),
        (
            "dl2_26_target_list_declared",
            {
                "target_list_source": "cfg.L2_TARGET_MODELS_LIST",
                "list_validated": True,
            },
        ),
        (
            "dl2_27_export_only_if_verified_multi",
            {"flatten_order_verified": True},
        ),
        (
            "dl2_28_full_ledger",
            {
                "declared_ids": list(L2_TARGET_MODELS_LIST),
                "ledger_ids": list(L2_TARGET_MODELS_LIST),
                "no_silent_drops": True,
            },
        ),
        # Layer 3
        ("dl3_01_p0_no_reweight", {"reweighted": False}),
        (
            "dl3_02_grn_clip_applied",
            {
                "provenance": {
                    "grn_weight_transform": {
                        "transform_method": "clip_nonpositive",
                        "n_edges_after": 268417,
                    }
                }
            },
        ),
        (
            "dl3_03_frozen_layer2",
            {
                "grn_recomputed": False,
                "using_frozen_layer2": True,
            },
        ),
        (
            "dl3_04_drug_monoplex",
            {
                "drug_monoplex_policy": "self_loop_per_drug_node",
            },
        ),
        (
            "dl3_05_no_aggregation",
            {
                "gene_df": pd.DataFrame(
                    {
                        "gene_symbol": ["PPARA", "PPARA", "PPARD", "PPARD"],
                        "layer": ["PPI", "GRN", "PPI", "GRN"],
                        "rwr_score": [0.22, 0.21, 0.20, 0.19],
                    }
                )
            },
        ),
        (
            "dl3_06_layer3_ppi_source",
            {
                "ppi_path_used": str(LOCAL_PPI_PRIOR_L3_TSV),
                "manifest_validated": True,
                "tf_restricted": False,
                "sha256_ok": True,
            },
        ),
        (
            "dl3_07_single_representative_run",
            {
                "mathematically_equivalent": True,
                "synthetic_class_node": True,
                "n_runs_executed": 1,
            },
        ),
        # Layer 4
        (
            "dl4_01_explicit_background",
            {
                "background_policy": "explicit_layer3_gene_universe",
                "background_is_none": False,
            },
        ),
        (
            "dl4_02_query_subset_background",
            {
                "genes_outside_background": set(),
            },
        ),
        (
            "dl4_03_preserve_ranking_order",
            {
                "query_was_sorted": False,
            },
        ),
        (
            "dl4_04_sanity_check_formula",
            {
                "M_is_background_restricted": True,
            },
        ),
        (
            "dl4_05_redundancy_same_universe",
            {
                "redundancy_restricted_to_background": True,
            },
        ),
        (
            "dl4_06_redundancy_on_significant",
            {
                "input_is_significant_only": True,
            },
        ),
    ]

    for _step, _ctx in _tests:
        try:
            validate_deadlock_rules(_step, **_ctx)
            print(f"  {_step:48s} → PASS")
        except AssertionError as _e:
            print(f"  {_step:48s} → FAIL: {str(_e)[:60]}")

    print("\ncheck_layer3_resources():")
    for k, v in check_layer3_resources().items():
        flag = "✓" if v is True else ("✗" if v is False else " ")
        print(f"  [{flag}] {k}: {v}")
