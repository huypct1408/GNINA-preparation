# ========================================
# CONFIGURATION LOADING
# ========================================

NOTEBOOK_DIR = Path.cwd()
print(f"📂 Notebook directory: {NOTEBOOK_DIR}")

# Load .env file
_env_file = NOTEBOOK_DIR / ".env"
if _env_file.exists():
    load_dotenv(_env_file, override=True)
    print(f"✅ Loaded .env from {_env_file}")
else:
    print(f"⚠️  No .env file found at {_env_file}")

# ========================================
# PATH RESOLUTION FUNCTION
# ========================================
def _resolve_path(env_var: str, fallback: str) -> str:
    """Resolve environment variables to absolute paths with fallbacks."""
    raw = os.environ.get(env_var, fallback)
    expanded = os.path.expanduser(raw)
    if not os.path.isabs(expanded):
        expanded = str(NOTEBOOK_DIR / expanded)
    return os.path.normpath(expanded)

# ========================================
# GLOBAL CONFIGURATION
# ========================================

# Base directories
BASE_DIR = _resolve_path("DOCKING_BASE_DIR", str(NOTEBOOK_DIR))
RESULTS_BASE_DIR = _resolve_path("RESULTS_BASE_DIR", f"{BASE_DIR}/results/phase3_docking_results")

# GNINA executable
GNINA_BIN = _resolve_path("GNINA_BIN", "/usr/local/bin/gnina")
GNINA_TIMEOUT_SEC = int(os.environ.get("GNINA_TIMEOUT_SEC", "3600"))
GPU_DEVICE = os.environ.get("GNINA_GPU_DEVICE", "0")

# Input files
PROTEIN_PATH = _resolve_path("PROTEIN_PATH", f"{BASE_DIR}/data/receptor.pdb")
REF_LIGAND = _resolve_path("REF_LIGAND", f"{BASE_DIR}/data/ref_ligand.sdf")
LIGAND_SDF_B1 = _resolve_path("LIGAND_SDF_WARHEAD_B1", f"{BASE_DIR}/data/B1_conformers.sdf")
LIGAND_SDF_B2 = _resolve_path("LIGAND_SDF_WARHEAD_B2", f"{BASE_DIR}/data/B2_conformers.sdf")
NUCLEOPHILE_DATA_PATH = _resolve_path("NUCLEOPHILE_DATA_PATH", f"{BASE_DIR}/results/phase2_nucleophile_data.json")

# Docking parameters
FLEX_RESIDUES = os.environ.get("FLEX_RESIDUES", "A:253,A:256,A:257,A:259")
NUM_MODES = int(os.environ.get("NUM_MODES", "10"))
EXHAUSTIVENESS = int(os.environ.get("EXHAUSTIVENESS", "64"))
SEED = os.environ.get("SEED", "42")
AUTOBOX_ADD = float(os.environ.get("AUTOBOX_ADD", "5.0"))  # Cast to float at load
AUTOBOX_EXTEND = int(os.environ.get("AUTOBOX_EXTEND", "1"))  # Cast to int at load

# Covalent docking parameters
WARHEAD_B1_PATTERN = os.environ.get("WARHEAD_B1_SMARTS", "[B;R1;r5](OC)N")
WARHEAD_B2_PATTERN = os.environ.get("WARHEAD_B2_SMARTS", "[B](O)(C)(O)")
COVALENT_BOND_DISTANCE_MIN = float(os.environ.get("COVALENT_BOND_DISTANCE_MIN", "1.5"))
COVALENT_BOND_DISTANCE_MAX = float(os.environ.get("COVALENT_BOND_DISTANCE_MAX", "2.4"))
COVALENT_FORMATION_SCORE_MIN = float(os.environ.get("COVALENT_FORMATION_SCORE_MIN", "0.70"))
AFFINITY_POOR_THRESHOLD = float(os.environ.get("AFFINITY_POOR_THRESHOLD", "-6.5"))

# Covalent reactor atom target (from .env) - DO NOT HARDCODE
COVALENT_REC_ATOM = os.environ.get("COVALENT_REC_ATOM", "A:248:SG")
COVALENT_LIG_ATOM_PATTERN = os.environ.get("COVALENT_LIG_ATOM_PATTERN", "[B;R1;r5](OC)N")
COVALENT_LIG_ATOM_POSITION_STR = os.environ.get("COVALENT_LIG_ATOM_POSITION", "")
COVALENT_LIG_ATOM_POSITION = None
if COVALENT_LIG_ATOM_POSITION_STR.strip():
    try:
        coords = [float(x.strip()) for x in COVALENT_LIG_ATOM_POSITION_STR.split(",")]
        if len(coords) == 3:
            COVALENT_LIG_ATOM_POSITION = tuple(coords)
    except ValueError:
        pass  # Use None if parsing fails
COVALENT_FIX_LIG_ATOM_POSITION = int(os.environ.get("COVALENT_FIX_LIG_ATOM_POSITION", "0")) == 1
COVALENT_BOND_ORDER = int(os.environ.get("COVALENT_BOND_ORDER", "1"))
COVALENT_OPTIMIZE_LIG = True

# Status constants
STATUS_PENDING = "PENDING"
STATUS_RUNNING = "RUNNING"
STATUS_DONE = "DONE"
STATUS_FAILED = "FAILED"
STATUS_DONE_WITH_WARNING = "DONE_WITH_WARNING"

print("\n" + "=" * 80)
print("🧬 COVALENT DOCKING PIPELINE v1.0")
print("=" * 80)
print(f"📂 Base directory:        {BASE_DIR}")
print(f"🔬 GNINA binary:          {GNINA_BIN}")
print(f"🧪 Protein:               {PROTEIN_PATH}")
print(f"📦 B1 Ligands:            {LIGAND_SDF_B1}")
print(f"📦 B2 Ligands:            {LIGAND_SDF_B2}")
print(f"🔬 Nucleophile data:      {NUCLEOPHILE_DATA_PATH}")
print(f"📁 Output directory:      {RESULTS_BASE_DIR}")
print(f"\n⚙️  Parameters:")
print(f"   - Covalent receptor atom: {COVALENT_REC_ATOM}")
print(f"   - Covalent ligand pattern: {COVALENT_LIG_ATOM_PATTERN}")
print(f"   - Covalent bond distance: {COVALENT_BOND_DISTANCE_MIN}-{COVALENT_BOND_DISTANCE_MAX}Å")
print(f"   - Formation score threshold: >{COVALENT_FORMATION_SCORE_MIN:.2f}")
print(f"   - Optimize ligand: {'Yes' if COVALENT_OPTIMIZE_LIG else 'No'}")
print(f"   - GNINA timeout: {GNINA_TIMEOUT_SEC}s")
print(f"   - GPU device: {GPU_DEVICE}")
print(f"   - Poses per run: {NUM_MODES}")
print(f"   - Exhaustiveness: {EXHAUSTIVENESS}")
print("=" * 80 + "\n")

# ========================================
# VALIDATION: Check environment configuration
# ========================================
print("\n🔍 Validating environment configuration...")

# 1. Check numeric parameters have no units
validation_errors = []

try:
    # These should be pure numbers, no units like \"A\" or \"Å\"
    test_float = float(AUTOBOX_ADD)
    test_float = float(AUTOBOX_EXTEND)
    test_float = float(COVALENT_BOND_DISTANCE_MIN)
    test_float = float(COVALENT_BOND_DISTANCE_MAX)
    test_float = float(COVALENT_FORMATION_SCORE_MIN)
    test_float = float(AFFINITY_POOR_THRESHOLD)
    test_int = int(NUM_MODES)
    test_int = int(EXHAUSTIVENESS)
    test_int = int(GNINA_TIMEOUT_SEC)
    print("✅ All numeric parameters valid (no units)")
except ValueError as e:
    validation_errors.append(f"❌ Numeric parameter parsing failed: {e}\\n   Check .env file - parameters should not have units (e.g., use '5.0' not '5.0A')")

# 2. Check COVALENT_REC_ATOM format
if COVALENT_REC_ATOM:
    parts = COVALENT_REC_ATOM.split(':')
    if len(parts) != 3:
        validation_errors.append(f"❌ COVALENT_REC_ATOM format invalid: '{COVALENT_REC_ATOM}'\\n   Expected format: 'CHAIN:RESNUM:ATOM' (e.g., 'A:248:SG')")
    else:
        try:
            int(parts[1])  # resnum should be integer
            print(f"✅ COVALENT_REC_ATOM format valid: {COVALENT_REC_ATOM}")
        except ValueError:
            validation_errors.append(f"❌ COVALENT_REC_ATOM residue number not integer: '{parts[1]}'")

# 3. Check file paths exist
if PROTEIN_PATH and not os.path.exists(PROTEIN_PATH):
    validation_errors.append(f"❌ PROTEIN_PATH not found: {PROTEIN_PATH}")
else:
    print(f"✅ PROTEIN_PATH exists: {PROTEIN_PATH}")

if LIGAND_SDF_B1 and not os.path.exists(LIGAND_SDF_B1):
    validation_errors.append(f"❌ LIGAND_SDF_B1 not found: {LIGAND_SDF_B1}")
else:
    print(f"✅ LIGAND_SDF_B1 exists: {LIGAND_SDF_B1}")

if LIGAND_SDF_B2 and not os.path.exists(LIGAND_SDF_B2):
    validation_errors.append(f"❌ LIGAND_SDF_B2 not found: {LIGAND_SDF_B2}")
else:
    print(f"✅ LIGAND_SDF_B2 exists: {LIGAND_SDF_B2}")

if NUCLEOPHILE_DATA_PATH and not os.path.exists(NUCLEOPHILE_DATA_PATH):
    validation_errors.append(f"❌ NUCLEOPHILE_DATA_PATH not found: {NUCLEOPHILE_DATA_PATH}")
else:
    print(f"✅ NUCLEOPHILE_DATA_PATH exists: {NUCLEOPHILE_DATA_PATH}")

if RESULTS_BASE_DIR and not os.path.exists(RESULTS_BASE_DIR):
    print(f"⚠️  RESULTS_BASE_DIR will be created: {RESULTS_BASE_DIR}")
else:
    print(f"✅ RESULTS_BASE_DIR exists: {RESULTS_BASE_DIR}")

# 4. Check GNINA binary exists
if GNINA_BIN and not os.path.exists(GNINA_BIN):
    validation_errors.append(f"❌ GNINA_BIN not found: {GNINA_BIN}\\n   Please ensure GNINA is installed and path is correct.")
else:
    print(f"✅ GNINA_BIN exists: {GNINA_BIN}")

# 5. Print summary
print("\n" + "=" * 80)
if validation_errors:
    print("❌ VALIDATION FAILED - Configuration errors found:")
    for error in validation_errors:
        print(f"  {error}")
    print("=" * 80)
    raise RuntimeError("Please fix configuration errors before running pipeline")
else:
    print("✅ All validations passed - Configuration is correct!")
    print("=" * 80)

# ========================================
# PROTEIN PATH RESOLUTION (Multi-Target)
# ========================================

def get_protein_path_for_target(target_name: str) -> str:
    """
    Dynamically resolve protein PDB path for a specific target.
    
    Args:
        target_name: Target identifier (e.g., 'PPARA_7BQ2')
    
    Returns:
        Full path to the prepared protein PDB file
    
    Raises:
        ValueError: If target not found in TARGET_MAPPING
    """
    # First, establish the protein data directory
    PROTEIN_DATA_DIR = _resolve_path("PROTEIN_DATA_DIR", f"{BASE_DIR}/data/proteins")
    
    mapping_key = f"TARGET_MAPPING_{target_name}"
    filename = os.environ.get(mapping_key)
    
    if not filename:
        raise ValueError(
            f"Target '{target_name}' not found in environment. "
            f"Please set {mapping_key} in .env file."
        )
    
    full_path = os.path.join(PROTEIN_DATA_DIR, filename)
    
    if not os.path.exists(full_path):
        raise ValueError(
            f"Protein file not found: {full_path}\n"
            f"Check that {mapping_key}={filename} in .env is correct."
        )
    
    return full_path
