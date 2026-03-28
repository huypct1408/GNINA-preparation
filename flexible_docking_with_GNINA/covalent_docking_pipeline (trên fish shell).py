

```
# ============================================================
# 🧬 GNINA Flexible COVALENT Docking Pipeline v2.7.0
# ============================================================
# Thay đổi so với v2.6.2:
#   🔴 COVALENT DOCKING MODE: -cnn_scoring none
#       CNN models chưa calibrated cho covalent docking.
#       Toàn bộ scoring chuyển sang Vinardo (classical).
#
#   🔴 CMD: --cnn_scoring none, --scoring vinardo,
#       --pose_sort_order 0 (sort by affinity)
#
#   🔴 OUTPUT: Xóa tất cả cột CNN (CNNscore, CNN_VS, 
#       CNNaffinity). Thay bằng Vinardo affinity.
#
#   🔴 ARBITRATION: Đơn giản hóa — chỉ dựa trên 
#       thermodynamic sanity (affinity).
#
#   🔴 RANKING: Theo minimizedAffinity (thấp hơn = tốt hơn),
#       KHÔNG còn CNN_VS ranking.
#
#   FAIR compliance:
#     F: Persistent naming (LIG_XXXX), structured dirs
#     A: Open formats (SDF, CSV, XLSX, PDB)
#     I: SMILES + standard prop names in SDF
#     R: Full provenance (command.txt, META.txt, STATUS.txt)
#
#   SMART goals addressed:
#     S: Covalent flexible docking with Vinardo scoring
#     M: Affinity (kcal/mol) as sole quantitative metric
#     A: Timeout-protected, GPU-accelerated minimization
#     R: Scientifically justified (CNN not calibrated)
#     T: Per-ligand timeout + batch progress tracking
# ============================================================

import os
import re
import subprocess
import time
import csv
import shutil
import socket
import traceback
from pathlib import Path
from rdkit import Chem

# =============================================================
# .env loading + LD_LIBRARY_PATH
# =============================================================
NOTEBOOK_DIR = Path.cwd()

_conda_prefix = os.environ.get("CONDA_PREFIX", "")
if _conda_prefix:
    _conda_lib = os.path.join(_conda_prefix, "lib")
    _current_ld = os.environ.get("LD_LIBRARY_PATH", "")
    if _conda_lib not in _current_ld:
        os.environ["LD_LIBRARY_PATH"] = (
            f"{_conda_lib}:{_current_ld}" if _current_ld else _conda_lib
        )
        print(f"🔧 Set LD_LIBRARY_PATH = {os.environ['LD_LIBRARY_PATH']}")

try:
    from dotenv import load_dotenv
    _env_file = NOTEBOOK_DIR / ".env"
    if _env_file.exists():
        load_dotenv(_env_file, override=True)
        print(f"📄 Loaded .env from {_env_file}")
except ImportError:
    print("⚠️ python-dotenv not installed. Using environment variables only.")


def _resolve_path(env_var: str, fallback: str) -> str:
    raw = os.environ.get(env_var, fallback)
    expanded = os.path.expanduser(raw)
    if not os.path.isabs(expanded):
        expanded = str(NOTEBOOK_DIR / expanded)
    return os.path.normpath(expanded)


# =========================
# GLOBAL CONFIG
# =========================
BASE_DIR = _resolve_path("DOCKING_BASE_DIR", str(NOTEBOOK_DIR))
RESULTS_DIR = f"{BASE_DIR}/docking_results"
PROTEIN_PATH = _resolve_path(
    "PROTEIN_PATH",
    f"{BASE_DIR}/data/mmp2_7xjo_ready_for_gnina.pdb",
)
REF_LIGAND = _resolve_path("REF_LIGAND", f"{BASE_DIR}/data/ref_ligand.sdf")
LIGAND_SDF = _resolve_path("LIGAND_SDF", f"{BASE_DIR}/data/ligands_prepared.sdf")
FLEX_RESIDUES = "A:7,A:83,A:130"
SEED = "42"
GPU_DEVICE = os.environ.get("GNINA_GPU_DEVICE", "0")

# ================================================================
# [v2.7.0] Scoring mode — COVALENT
# ================================================================
# CNN models chưa được calibrated cho covalent docking.
# Nguồn: GNINA runtime warning + Ragoza et al. 2022
# → Dùng Vinardo (classical empirical scoring function)
# → Pose sort by affinity (lower = better binding)
SCORING_MODE = "default"         # Classical scoring function
CNN_SCORING = "none"             # Tắt CNN hoàn toàn
POSE_SORT_ORDER = "Energy"           # Sort by affinity column (default Vina/Vinardo)

# ================================================================
# Arbitration thresholds — [v2.7.0] Chỉ affinity-based
# ================================================================
# Không còn CNN_VS, CNNscore → chỉ dùng Vinardo affinity
# để đánh giá thermodynamic sanity.
#
# SMART justification:
#   S: Flag ligands with physically meaningless binding
#   M: Thresholds in kcal/mol — universally measurable
#   A: Conservative thresholds based on literature
#   R: Aligned with Vinardo scoring range
#   T: Applied automatically at summary generation
RED_FLAG_AFFINITY_POOR = -6.5    # kcal/mol — gắn kết quá yếu
RED_FLAG_AFFINITY_MARGINAL = -7.0  # kcal/mol — marginal (info only)

# ================================================================
# [v2.7.0] Vinardo affinity interpretation guide
# ================================================================
# Vinardo affinity (kcal/mol) — lower (more negative) = better
#   < -10.0 : Excellent (very strong predicted binding)
#   -10.0 to -8.0 : Good (strong binding, drug-like)
#   -8.0 to -7.0 : Moderate (worth investigating)
#   -7.0 to -6.5 : Marginal (weak, needs optimization)
#   > -6.5 : Poor (likely not viable)
#   > 0 : Repulsive (physically meaningless)
AFFINITY_EXCELLENT = -10.0
AFFINITY_GOOD = -8.0
AFFINITY_MODERATE = -7.0

# =====================================================
# GNINA timeout
# =====================================================
GNINA_TIMEOUT_SEC = int(os.environ.get("GNINA_TIMEOUT_SEC", "3600"))

# GNINA_BIN — auto-detect
_env_bin = os.environ.get("GNINA_BIN", "")
if _env_bin and os.path.isfile(_env_bin):
    GNINA_BIN = _env_bin
elif shutil.which("gnina"):
    GNINA_BIN = shutil.which("gnina")
else:
    GNINA_BIN = f"{BASE_DIR}/bin/gnina"

# Status constants
STATUS_PENDING = "PENDING"
STATUS_RUNNING = "RUNNING"
STATUS_DONE = "DONE"
STATUS_FAILED = "FAILED"


# =========================
# Subprocess env
# =========================
def _get_subprocess_env() -> dict:
    env = os.environ.copy()
    conda_prefix = env.get("CONDA_PREFIX", "")
    if conda_prefix:
        conda_lib = os.path.join(conda_prefix, "lib")
        current_ld = env.get("LD_LIBRARY_PATH", "")
        if conda_lib not in current_ld:
            env["LD_LIBRARY_PATH"] = (
                f"{conda_lib}:{current_ld}" if current_ld else conda_lib
            )
    return env


# =========================
# STATUS MANAGEMENT — giữ nguyên
# =========================
def write_status(lig_root: str, status: str, **kwargs):
    status_file = os.path.join(lig_root, "STATUS.txt")

    if status == STATUS_RUNNING:
        with open(status_file, "w") as f:
            f.write(f"STATUS={status}\n")
            f.write(f"HOST={socket.gethostname()}\n")
            f.write(f"START_TIME={time.strftime('%Y-%m-%d %H:%M:%S')}\n")
            # [v2.7.0] FAIR — ghi scoring mode vào provenance
            f.write(f"SCORING_MODE={SCORING_MODE}\n")
            f.write(f"CNN_SCORING={CNN_SCORING}\n")
            f.write(f"PIPELINE_VERSION=2.7.0\n")
    else:
        with open(status_file, "a") as f:
            f.write(f"STATUS={status}\n")
            for key, value in kwargs.items():
                f.write(f"{key.upper()}={value}\n")
            f.write(f"END_TIME={time.strftime('%Y-%m-%d %H:%M:%S')}\n")


def read_status(lig_root: str) -> str:
    status_file = os.path.join(lig_root, "STATUS.txt")

    if not os.path.exists(status_file):
        return STATUS_PENDING

    last_status = STATUS_PENDING
    try:
        with open(status_file, "r") as f:
            for line in f:
                if line.startswith("STATUS="):
                    last_status = line.strip().split("=", 1)[1]
    except Exception:
        pass

    return last_status


def get_status_details(lig_root: str) -> dict:
    status_file = os.path.join(lig_root, "STATUS.txt")
    details = {}

    if os.path.exists(status_file):
        with open(status_file, "r") as f:
            for line in f:
                if "=" in line:
                    key, value = line.strip().split("=", 1)
                    details[key] = value

    return details


# =========================
# UTILS — giữ nguyên
# =========================
def sanitize_name(name: str, max_len: int = 80) -> str:
    if not name:
        return "NA"

    name = name.strip().replace(" ", "_")
    name = re.sub(r'[\/:*?"<>|\\]', "", name)
    name = re.sub(r"_+", "_", name)
    name = name.strip("_")

    return name[:max_len] if name else "NA"


# =========================
# PREPARE ROOT FOLDERS — giữ nguyên
# =========================
def prepare_root_folders():
    dirs = [
        f"{RESULTS_DIR}/protein",
        f"{RESULTS_DIR}/reference",
        f"{RESULTS_DIR}/ligands",
        f"{RESULTS_DIR}/summary",
    ]

    for d in dirs:
        os.makedirs(d, exist_ok=True)

    protein_dest = f"{RESULTS_DIR}/protein/receptor.pdb"
    ref_dest = f"{RESULTS_DIR}/reference/ref_ligand.sdf"

    if not os.path.exists(protein_dest):
        shutil.copy(PROTEIN_PATH, protein_dest)
        print(f"✔ Copied receptor to {protein_dest}")

    if not os.path.exists(ref_dest):
        shutil.copy(REF_LIGAND, ref_dest)
        print(f"✔ Copied reference ligand to {ref_dest}")

    # [v2.7.0] FAIR — ghi pipeline config vào summary
    config_file = os.path.join(RESULTS_DIR, "summary", "pipeline_config.txt")
    with open(config_file, "w") as f:
        f.write(f"PIPELINE_VERSION=2.7.0\n")
        f.write(f"DOCKING_TYPE=covalent_flexible\n")
        f.write(f"SCORING_MODE={SCORING_MODE}\n")
        f.write(f"CNN_SCORING={CNN_SCORING}\n")
        f.write(f"POSE_SORT_ORDER=affinity (lower=better)\n")
        f.write(f"FLEX_RESIDUES={FLEX_RESIDUES}\n")
        f.write(f"SEED={SEED}\n")
        f.write(f"EXHAUSTIVENESS=32\n")
        f.write(f"NUM_MODES=10\n")
        f.write(f"TIMEOUT_SEC={GNINA_TIMEOUT_SEC}\n")
        f.write(f"RED_FLAG_AFFINITY_POOR={RED_FLAG_AFFINITY_POOR}\n")
        f.write(f"GNINA_BIN={GNINA_BIN}\n")
        f.write(f"PROTEIN={PROTEIN_PATH}\n")
        f.write(f"REF_LIGAND={REF_LIGAND}\n")
        f.write(f"GENERATED={time.strftime('%Y-%m-%d %H:%M:%S')}\n")
    print(f"✔ Pipeline config saved to {config_file}")


# =========================
# SPLIT LIGANDS — giữ nguyên
# =========================
def split_ligands(input_sdf: str, ligands_root: str) -> list:
    os.makedirs(ligands_root, exist_ok=True)

    cleaned_sdf = os.path.join(ligands_root, "_cleaned_input.sdf")
    with open(input_sdf, "rb") as f_in:
        raw = f_in.read()

    if b"\r" in raw:
        raw = raw.replace(b"\r\n", b"\n").replace(b"\r", b"\n")
        with open(cleaned_sdf, "wb") as f_out:
            f_out.write(raw)
        sdf_to_parse = cleaned_sdf
        print(f"🔧 Found \\r in input — normalized line endings")
    else:
        sdf_to_parse = input_sdf
        print(f"✔ Input SDF has clean line endings — no fix needed")

    suppl = Chem.SDMolSupplier(sdf_to_parse, removeHs=False, sanitize=True)
    ligands = []
    mapping = []
    forced_3d_list = []

    for idx, mol in enumerate(suppl, start=1):
        if mol is None:
            print(f"⚠️ Skipping invalid molecule at index {idx}")
            continue

        lig_id = f"LIG_{idx:04d}"

        conf = mol.GetConformer() if mol.GetNumConformers() > 0 else None
        if conf is not None:
            has_z = any(
                abs(conf.GetAtomPosition(i).z) > 0.001
                for i in range(mol.GetNumAtoms())
            )
            if has_z and not conf.Is3D():
                conf.Set3D(True)
                forced_3d_list.append(lig_id)
                print(f"🔧 {lig_id}: forced 3D tag (had Z coordinates)")

        orig_name = mol.GetProp("_Name") if mol.HasProp("_Name") else "NA"
        safe_name = sanitize_name(orig_name)

        lig_dirname = f"{lig_id}__{safe_name}"
        lig_root = os.path.join(ligands_root, lig_dirname)
        input_dir = os.path.join(lig_root, "input")
        os.makedirs(input_dir, exist_ok=True)

        ligand_sdf = os.path.join(input_dir, "ligand.sdf")
        writer = Chem.SDWriter(ligand_sdf)
        writer.SetForceV3000(False)
        writer.write(mol)
        writer.close()

        smiles = Chem.MolToSmiles(mol)
        with open(os.path.join(lig_root, "META.txt"), "w") as f:
            f.write(f"ID={lig_id}\n")
            f.write(f"DIR_NAME={lig_dirname}\n")
            f.write(f"ORIGINAL_NAME={orig_name}\n")
            f.write(f"SMILES={smiles}\n")
            f.write(f"SDF_INDEX={idx}\n")

        ligands.append(
            {
                "lig_id": lig_id,
                "lig_dirname": lig_dirname,
                "lig_root": lig_root,
                "ligand_sdf": ligand_sdf,
                "orig_name": orig_name,
                "smiles": smiles,
            }
        )
        mapping.append((lig_id, lig_dirname, orig_name, smiles))

    mapping_file = os.path.join(ligands_root, "ligand_mapping.csv")
    with open(mapping_file, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(["ID", "DIR_NAME", "ORIGINAL_NAME", "SMILES"])
        for row in mapping:
            writer.writerow(row)

    if forced_3d_list:
        forced_log = os.path.join(ligands_root, "forced_3d_log.txt")
        with open(forced_log, "w") as f:
            f.write("# Ligands forced from 2D tag to 3D tag\n")
            f.write(f"# Total: {len(forced_3d_list)}\n")
            for lid in forced_3d_list:
                f.write(f"{lid}\n")
        print(f"📝 Forced 3D log: {forced_log} ({len(forced_3d_list)} ligands)")

    count_3d = 0
    count_2d = 0
    for lig in ligands:
        mol = next(Chem.SDMolSupplier(lig["ligand_sdf"], removeHs=False))
        if mol is not None and mol.GetNumConformers() > 0:
            if mol.GetConformer().Is3D():
                count_3d += 1
            else:
                count_2d += 1
    print(f"✔ Validation: {count_3d} molecules 3D, {count_2d} still 2D")
    if count_2d > 0:
        print(f"⚠️ WARNING: {count_2d} molecules still tagged 2D — check input geometry")

    if os.path.exists(cleaned_sdf):
        os.remove(cleaned_sdf)

    print(f"✔ Split {len(ligands)} ligands")
    print(f"✔ Mapping written to {mapping_file}")

    return ligands


# ================================================================
# [v2.7.0] SCORE PARSING — Vinardo/Classical only (NO CNN)
# ================================================================
def parse_top_poses(sdf_path: str, n: int = 3) -> list:
    """
    Extract top N poses from docked SDF.

    [v2.7.0] Covalent docking mode:
      - CNN scoring disabled → CNNscore/CNN_VS/CNNaffinity = None
      - Poses sorted by minimizedAffinity (lower = better)
      - Only classical scoring properties extracted:
          minimizedAffinity, vinardo (if present)
    
    FAIR compliance:
      - Property names match GNINA SDF output exactly
      - Pose rank preserved from file order (no re-sort)
    """
    try:
        suppl = Chem.SDMolSupplier(sdf_path, removeHs=False)
        poses = []

        for pose_idx, mol in enumerate(suppl, start=1):
            if mol is None:
                continue

            pose_data = {"pose_rank": pose_idx}

            # [v2.7.0] Classical scoring properties only
            # CNN properties will be None since --cnn_scoring none
            score_props = [
                "minimizedAffinity",
                "vinardo",
                "CNNscore",       # Will be None — kept for schema consistency
                "CNNaffinity",    # Will be None — kept for schema consistency
                "CNN_VS",         # Will be None — kept for schema consistency
            ]

            for prop in score_props:
                if mol.HasProp(prop):
                    try:
                        pose_data[prop] = float(mol.GetProp(prop))
                    except ValueError:
                        pose_data[prop] = mol.GetProp(prop)
                else:
                    pose_data[prop] = None

            poses.append(pose_data)

            if len(poses) >= n:
                break

        return poses

    except Exception as e:
        print(f"⚠️ Error parsing scores: {e}")
        return []


# ================================================================
# [v2.7.0] classify_red_flag — Affinity-only (no CNN metrics)
# ================================================================
def classify_red_flag(affinity) -> str:
    """
    Thermodynamic Sanity Check — Covalent Docking Mode.

    [v2.7.0] CNN scoring disabled → sanity check dựa hoàn toàn
    trên Vinardo minimizedAffinity.

    Flags:
      🟡 POSITIVE_AFFINITY: affinity > 0 kcal/mol
         → Repulsive interaction — physically meaningless
      🟠 POOR_AFFINITY: affinity ≥ -6.5 kcal/mol
         → Binding too weak for biological significance
      (empty): Pass — thermodynamically reasonable

    SMART:
      S: Identify thermodynamically invalid poses
      M: kcal/mol thresholds
      A: Conservative cutoffs from literature
      R: Consistent with Vinardo scoring range
      T: Applied at summary generation time
    """
    if affinity is None:
        return ""

    if affinity > 0:
        return "🟡 POSITIVE_AFFINITY"

    if affinity >= RED_FLAG_AFFINITY_POOR:
        return "🟠 POOR_AFFINITY"

    return ""


def classify_affinity_tier(affinity) -> str:
    """
    [v2.7.0] Classify affinity into interpretive tiers.
    
    Used in Excel for color-coding and in CSV for quick filtering.
    Tiers based on Vinardo empirical scoring function range.
    """
    if affinity is None:
        return "N/A"
    if affinity > 0:
        return "REPULSIVE"
    if affinity >= RED_FLAG_AFFINITY_POOR:
        return "POOR"
    if affinity >= AFFINITY_MODERATE:
        return "MARGINAL"
    if affinity >= AFFINITY_GOOD:
        return "MODERATE"
    if affinity >= AFFINITY_EXCELLENT:
        return "GOOD"
    return "EXCELLENT"


def get_best_pose_summary(sdf_path: str) -> dict:
    """
    [v2.7.0] Lấy pose tốt nhất (pose 1 = affinity thấp nhất)
    và tính sanity flag.
    
    Với --cnn_scoring none, pose 1 là pose có affinity thấp nhất
    (--pose_sort_order 0 = sort by first scoring column = affinity).
    """
    poses = parse_top_poses(sdf_path, n=1)
    if not poses:
        return {
            "minimizedAffinity": None,
            "vinardo": None,
            "red_flag": "",
            "affinity_tier": "N/A",
        }

    best = poses[0]
    aff = best.get("minimizedAffinity")
    best["red_flag"] = classify_red_flag(aff)
    best["affinity_tier"] = classify_affinity_tier(aff)
    return best


# ================================================================
# parse_best_score — giữ nguyên logic, chỉ dùng minimizedAffinity
# ================================================================
def parse_best_score(sdf_path: str):
    """Extract best minimizedAffinity from docked SDF.
    Returns None nếu không parse được.
    """
    poses = parse_top_poses(sdf_path, n=1)
    if poses and poses[0].get("minimizedAffinity") is not None:
        return poses[0]["minimizedAffinity"]
    return None


# ================================================================
# [v2.7.0] RUN GNINA — Covalent docking command
# ================================================================
def run_gnina(ligand_info: dict, idx: int, total: int) -> bool:
    """Run GNINA covalent flexible docking for a single ligand.

    [v2.7.0] Thay đổi so với v2.6.2:
      - --cnn_scoring none (CNN not calibrated for covalent)
      - --scoring vinardo (classical empirical scoring)
      - --pose_sort_order 0 (sort by affinity, not CNNscore)
      - Removed --cnn_empirical_weight
    """
    lig_id = ligand_info["lig_id"]
    lig_root = ligand_info["lig_root"]
    ligand_sdf = ligand_info["ligand_sdf"]

    out_dir = os.path.join(lig_root, "output")
    log_dir = os.path.join(lig_root, "logs")
    os.makedirs(out_dir, exist_ok=True)
    os.makedirs(log_dir, exist_ok=True)

    out_lig = os.path.join(out_dir, "docked.sdf")
    out_flex = os.path.join(out_dir, "flex_residues.pdb")
    log_file = os.path.join(log_dir, "gnina.log")
    stderr_file = os.path.join(log_dir, "gnina_stderr.log")
    cmd_file = os.path.join(log_dir, "command.txt")

    write_status(lig_root, STATUS_RUNNING)

    print(f"\n🔄 [{idx}/{total}] Docking {lig_id} ...")
    start = time.time()

    # ── [v2.7.0] COVALENT DOCKING COMMAND ──
    # Key changes from v2.6.2:
    #   1. --cnn_scoring none     (was: --cnn_scoring rescore)
    #   2. --scoring vinardo      (explicit classical scoring)
    #   3. --pose_sort_order 0    (was: --pose_sort_order CNNscore)
    #   4. REMOVED: --cnn_empirical_weight 1.0
    cmd = [
        GNINA_BIN,
        "-r",
        f"{RESULTS_DIR}/protein/receptor.pdb",
        "-l",
        ligand_sdf,
        "--autobox_ligand",
        f"{RESULTS_DIR}/reference/ref_ligand.sdf",
        "--autobox_add",
        "5",
        "--autobox_extend",
        "1",
        "--flexres",
        FLEX_RESIDUES,
        "--num_modes",
        "10",
        "--exhaustiveness",
        "64",
        # Thêm vô các function liên quan tới covalent docking
        # ── [THÊM VÀO ĐÂY] KHÓA TỌA ĐỘ KIM LOẠI (COVALENT DOCKING) ──
        "--covalent_rec_atom", "A:301:ZN",
        "--covalent_lig_atom_pattern", "[OX1;$([O]C=O)]",
        "--covalent_lig_atom_position", "6.739,10.721,31.893",
        "--covalent_fix_lig_atom_position",
        "--covalent_optimize_lig",
        # ── [v2.7.0] Classical scoring only ──
        "--scoring",
        SCORING_MODE,              # default hoặc vina
        "--cnn_scoring",
        CNN_SCORING,               # none
        "--pose_sort_order",
        POSE_SORT_ORDER,           # Energy (affinity)
        # ── GPU still used for minimization ──
        "--device",
        GPU_DEVICE,
        "--seed",
        SEED,
        "--atom_term_data",
        "-o",
        out_lig,
        "--out_flex",
        out_flex,
        "--log",
        log_file,
    ]

    # [v2.7.0] FAIR — ghi command đầy đủ cho reproducibility
    with open(cmd_file, "w") as f:
        f.write(f"# GNINA Covalent Docking — Pipeline v2.7.0\n")
        f.write(f"# Generated: {time.strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write(f"# Scoring: {SCORING_MODE} (CNN disabled)\n")
        f.write(f"# Pose sort: affinity (lower = better)\n\n")
        f.write(" \\\n    ".join(cmd))

    try:
        with open(stderr_file, "w") as stderr_f:
            subprocess.run(
                cmd,
                check=True,
                stdout=subprocess.DEVNULL,
                stderr=stderr_f,
                text=True,
                env=_get_subprocess_env(),
                timeout=GNINA_TIMEOUT_SEC,
            )

        if not os.path.exists(out_lig) or os.path.getsize(out_lig) == 0:
            raise RuntimeError("GNINA produced empty output SDF")

        elapsed = (time.time() - start) / 60

        best_score = parse_best_score(out_lig)

        write_status(
            lig_root,
            STATUS_DONE,
            elapsed_min=f"{elapsed:.2f}",
            best_affinity=(
                f"{best_score:.4f}" if best_score is not None else "NA"
            ),
        )

        score_str = f"{best_score:.4f}" if best_score is not None else "N/A"
        tier = classify_affinity_tier(best_score)
        print(
            f"✅ [{idx}/{total}] {lig_id} DONE in {elapsed:.2f} min"
            f" (affinity: {score_str} [{tier}])"
        )
        return True

    except subprocess.TimeoutExpired:
        elapsed = (time.time() - start) / 60
        timeout_min = GNINA_TIMEOUT_SEC / 60
        write_status(
            lig_root,
            STATUS_FAILED,
            elapsed_min=f"{elapsed:.2f}",
            error=f"TIMEOUT after {GNINA_TIMEOUT_SEC}s ({timeout_min:.0f}min)",
        )
        print(
            f"⏰ [{idx}/{total}] {lig_id} TIMEOUT after"
            f" {timeout_min:.0f} min — skipping"
        )
        return False

    except subprocess.CalledProcessError as e:
        elapsed = (time.time() - start) / 60
        write_status(
            lig_root,
            STATUS_FAILED,
            elapsed_min=f"{elapsed:.2f}",
            error=f"GNINA exit code {e.returncode}",
        )
        print(
            f"❌ [{idx}/{total}] {lig_id} FAILED: GNINA exit code {e.returncode}"
        )
        return False

    except Exception as e:
        elapsed = (time.time() - start) / 60
        write_status(
            lig_root,
            STATUS_FAILED,
            elapsed_min=f"{elapsed:.2f}",
            error=str(e),
            traceback=traceback.format_exc().replace("\n", " | "),
        )
        print(f"❌ [{idx}/{total}] {lig_id} FAILED: {e}")
        return False


# ================================================================
# [v2.7.0] EXCEL SUMMARY — Covalent docking (Vinardo only)
# ================================================================
def generate_excel_summary(ligands: list, summary_dir: str):
    try:
        from openpyxl import Workbook
        from openpyxl.styles import Font, PatternFill, Alignment, Border, Side

        USE_OPENPYXL = True
    except ImportError:
        print("⚠️ openpyxl not installed, generating CSV instead")
        USE_OPENPYXL = False

    results = []
    for lig in ligands:
        lig_root = lig["lig_root"]
        out_sdf = os.path.join(lig_root, "output", "docked.sdf")

        status = read_status(lig_root)
        details = get_status_details(lig_root)

        top_poses = []
        best_pose = {}
        if status == STATUS_DONE and os.path.exists(out_sdf):
            top_poses = parse_top_poses(out_sdf, n=3)
            best_pose = get_best_pose_summary(out_sdf)

        results.append(
            {
                "lig_id": lig["lig_id"],
                "lig_dirname": lig["lig_dirname"],
                "orig_name": lig["orig_name"],
                "smiles": lig["smiles"],
                "status": status,
                "elapsed_min": details.get("ELAPSED_MIN", ""),
                "top_poses": top_poses,
                "best_affinity": best_pose.get("minimizedAffinity"),
                "best_vinardo": best_pose.get("vinardo"),
                "affinity_tier": best_pose.get("affinity_tier", "N/A"),
                "red_flag": best_pose.get("red_flag", ""),
            }
        )

    # [v2.7.0] Sort by affinity (lower = better = first)
    results_by_affinity = sorted(
        results,
        key=lambda x: (
            x["best_affinity"] is None,  # None goes last
            x["best_affinity"] if x["best_affinity"] is not None else 0,
        ),
    )

    if USE_OPENPYXL:
        _generate_xlsx_v27(results, results_by_affinity, summary_dir)
    else:
        _generate_csv_fallback_v27(results_by_affinity, summary_dir)


def _generate_xlsx_v27(
    results: list, results_by_affinity: list, summary_dir: str
):
    from openpyxl import Workbook
    from openpyxl.styles import Font, PatternFill, Alignment, Border, Side
    from openpyxl.utils import get_column_letter

    wb = Workbook()

    header_font = Font(bold=True, color="FFFFFF", size=11)
    header_fill_blue = PatternFill(
        start_color="2F5496", end_color="2F5496", fill_type="solid"
    )
    header_fill_green = PatternFill(
        start_color="548235", end_color="548235", fill_type="solid"
    )
    header_fill_red = PatternFill(
        start_color="C00000", end_color="C00000", fill_type="solid"
    )
    header_fill_purple = PatternFill(
        start_color="7030A0", end_color="7030A0", fill_type="solid"
    )
    done_fill = PatternFill(
        start_color="C6EFCE", end_color="C6EFCE", fill_type="solid"
    )
    failed_fill = PatternFill(
        start_color="FFC7CE", end_color="FFC7CE", fill_type="solid"
    )
    excellent_fill = PatternFill(
        start_color="92D050", end_color="92D050", fill_type="solid"
    )
    good_fill = PatternFill(
        start_color="A9D18E", end_color="A9D18E", fill_type="solid"
    )
    moderate_fill = PatternFill(
        start_color="FFEB9C", end_color="FFEB9C", fill_type="solid"
    )
    marginal_fill = PatternFill(
        start_color="FFD93D", end_color="FFD93D", fill_type="solid"
    )
    poor_fill = PatternFill(
        start_color="FFA94D", end_color="FFA94D", fill_type="solid"
    )
    repulsive_fill = PatternFill(
        start_color="FF6B6B", end_color="FF6B6B", fill_type="solid"
    )
    thin_border = Border(
        left=Side(style="thin"),
        right=Side(style="thin"),
        top=Side(style="thin"),
        bottom=Side(style="thin"),
    )
    center_align = Alignment(horizontal="center", vertical="center")

    def _write_header(ws, row, headers, fill):
        for col, header in enumerate(headers, start=1):
            cell = ws.cell(row=row, column=col, value=header)
            cell.font = header_font
            cell.fill = fill
            cell.alignment = center_align
            cell.border = thin_border

    def _set_col_widths(ws, widths: dict):
        for col_letter, width in widths.items():
            ws.column_dimensions[col_letter].width = width

    def _get_tier_fill(tier: str):
        """Return fill color based on affinity tier."""
        tier_fills = {
            "EXCELLENT": excellent_fill,
            "GOOD": good_fill,
            "MODERATE": moderate_fill,
            "MARGINAL": marginal_fill,
            "POOR": poor_fill,
            "REPULSIVE": repulsive_fill,
        }
        return tier_fills.get(tier, None)

    # ── SHEET 1: Intra-Ligand Poses ──
    # [v2.7.0] No CNN columns — only affinity per pose
    ws1 = wb.active
    ws1.title = "Intra-Ligand_Poses"

    ws1.merge_cells("A1:N1")
    note_cell = ws1.cell(
        row=1,
        column=1,
        value=(
            "COVALENT DOCKING — INTRA-LIGAND POSE ANALYSIS: "
            "Poses sorted by Vinardo minimizedAffinity (lower = better). "
            "CNN scoring DISABLED (not calibrated for covalent docking). "
            "Pose 1 = best affinity pose."
        ),
    )
    note_cell.font = Font(bold=True, italic=True, size=10, color="2F5496")

    headers_s1 = [
        "#",                   # A
        "ID",                  # B
        "Original_Name",       # C
        "Status",              # D
        "Time_min",            # E
        "P1_Affinity",         # F
        "P1_Tier",             # G
        "P1_Flag",             # H
        "P2_Affinity",         # I
        "P2_Tier",             # J
        "P3_Affinity",         # K
        "P3_Tier",             # L
        "Dir_Name",            # M
        "SMILES",              # N
    ]

    _write_header(ws1, 2, headers_s1, header_fill_blue)
    ws1.freeze_panes = "A3"

    for row_idx, data in enumerate(results_by_affinity, start=3):
        rank = row_idx - 2
        top_poses = data["top_poses"]

        row_data = [
            rank,
            data["lig_id"],
            data["orig_name"],
            data["status"],
            data["elapsed_min"],
        ]

        # Pose 1
        if len(top_poses) >= 1:
            p = top_poses[0]
            aff = p.get("minimizedAffinity")
            tier = classify_affinity_tier(aff)
            flag = classify_red_flag(aff)
            row_data.extend([aff, tier, flag])
        else:
            row_data.extend(["", "", ""])

        # Pose 2
        if len(top_poses) >= 2:
            p = top_poses[1]
            aff = p.get("minimizedAffinity")
            tier = classify_affinity_tier(aff)
            row_data.extend([aff, tier])
        else:
            row_data.extend(["", ""])

        # Pose 3
        if len(top_poses) >= 3:
            p = top_poses[2]
            aff = p.get("minimizedAffinity")
            tier = classify_affinity_tier(aff)
            row_data.extend([aff, tier])
        else:
            row_data.extend(["", ""])

        row_data.extend([data["lig_dirname"], data["smiles"]])

        for col_idx, value in enumerate(row_data, start=1):
            cell = ws1.cell(row=row_idx, column=col_idx, value=value)
            cell.border = thin_border

            # Status coloring
            if col_idx == 4:
                if value == STATUS_DONE:
                    cell.fill = done_fill
                elif value == STATUS_FAILED:
                    cell.fill = failed_fill

            # Tier coloring (columns G, J, L)
            if col_idx in (7, 10, 12):
                tier_fill = _get_tier_fill(str(value))
                if tier_fill:
                    cell.fill = tier_fill
                    if value in ("POOR", "REPULSIVE"):
                        cell.font = Font(bold=True, color="C00000")

            # Flag coloring (column H)
            if col_idx == 8:
                if "POSITIVE" in str(value):
                    cell.fill = repulsive_fill
                    cell.font = Font(bold=True)
                elif "POOR" in str(value):
                    cell.fill = poor_fill
                    cell.font = Font(bold=True)

    _set_col_widths(ws1, {
        "A": 5, "B": 12, "C": 35, "D": 10, "E": 9,
        "F": 13, "G": 12, "H": 22,
        "I": 13, "J": 12,
        "K": 13, "L": 12,
        "M": 35, "N": 60,
    })

    # ── SHEET 2: Affinity Ranking ──
    # [v2.7.0] Replaces CNN_VS ranking with affinity ranking
    ws2 = wb.create_sheet(title="Affinity_Ranking")

    ws2.merge_cells("A1:I1")
    note2 = ws2.cell(
        row=1,
        column=1,
        value=(
            "COVALENT DOCKING — AFFINITY RANKING: "
            "Ranked by Vinardo minimizedAffinity (lower = stronger binding). "
            "CNN scoring disabled — ranking based purely on classical "
            "force-field scoring. "
            "Tier: EXCELLENT(<-10) | GOOD(-10 to -8) | MODERATE(-8 to -7) "
            "| MARGINAL(-7 to -6.5) | POOR(>-6.5) | REPULSIVE(>0)."
        ),
    )
    note2.font = Font(bold=True, italic=True, size=10, color="548235")

    headers_s2 = [
        "Rank",            # A
        "ID",              # B
        "Original_Name",   # C
        "Affinity_kcal",   # D
        "Tier",            # E
        "Sanity_Flag",     # F
        "Time_min",        # G
        "Dir_Name",        # H
        "SMILES",          # I
    ]

    _write_header(ws2, 2, headers_s2, header_fill_green)
    ws2.freeze_panes = "A3"

    top_done = [r for r in results_by_affinity if r["status"] == STATUS_DONE]

    for row_idx, data in enumerate(top_done, start=3):
        rank = row_idx - 2

        row_data = [
            rank,
            data["lig_id"],
            data["orig_name"],
            data["best_affinity"],
            data["affinity_tier"],
            data["red_flag"],
            data["elapsed_min"],
            data["lig_dirname"],
            data["smiles"],
        ]

        for col_idx, value in enumerate(row_data, start=1):
            cell = ws2.cell(row=row_idx, column=col_idx, value=value)
            cell.border = thin_border

            # Affinity coloring (column D)
            if col_idx == 4 and isinstance(value, (int, float)):
                if value > 0:
                    cell.fill = repulsive_fill
                    cell.font = Font(bold=True, color="C00000")
                elif value >= RED_FLAG_AFFINITY_POOR:
                    cell.fill = poor_fill

            # Tier coloring (column E)
            if col_idx == 5:
                tier_fill = _get_tier_fill(str(value))
                if tier_fill:
                    cell.fill = tier_fill
                    if value in ("POOR", "REPULSIVE"):
                        cell.font = Font(bold=True, color="C00000")

            # Flag coloring (column F)
            if col_idx == 6:
                if "POSITIVE" in str(value):
                    cell.fill = repulsive_fill
                    cell.font = Font(bold=True)
                elif "POOR" in str(value):
                    cell.fill = poor_fill
                    cell.font = Font(bold=True)

    _set_col_widths(ws2, {
        "A": 7, "B": 12, "C": 35, "D": 13, "E": 12,
        "F": 22, "G": 9, "H": 35, "I": 60,
    })

    # ── SHEET 3: Sanity Flags ──
    ws3 = wb.create_sheet(title="Sanity_Flags")

    ws3.merge_cells("A1:I1")
    note3 = ws3.cell(
        row=1,
        column=1,
        value=(
            "COVALENT DOCKING — THERMODYNAMIC SANITY CHECK: "
            "🟡 POSITIVE_AFFINITY = Affinity > 0 kcal/mol (repulsive). "
            f"🟠 POOR_AFFINITY = Affinity ≥ {RED_FLAG_AFFINITY_POOR} "
            "kcal/mol (too weak). "
            "CNN metrics NOT available (disabled for covalent docking). "
            "Manual visual inspection required for flagged ligands."
        ),
    )
    note3.font = Font(bold=True, italic=True, size=10, color="C00000")

    headers_s3 = [
        "#",
        "ID",
        "Original_Name",
        "Flag_Type",
        "Affinity_kcal",
        "Tier",
        "Concern",
        "Dir_Name",
        "SMILES",
    ]

    _write_header(ws3, 2, headers_s3, header_fill_red)
    ws3.freeze_panes = "A3"

    flagged = [r for r in results_by_affinity if r["red_flag"]]

    if flagged:
        for row_idx, data in enumerate(flagged, start=3):
            flag = data["red_flag"]
            aff = data["best_affinity"]

            if "POSITIVE" in flag:
                concern = (
                    f"Affinity={aff:.2f} > 0 (repulsive). "
                    "Ligand repelled from binding pocket — "
                    "physically meaningless pose."
                )
            elif "POOR" in flag:
                concern = (
                    f"Affinity={aff:.2f} ≥ "
                    f"{RED_FLAG_AFFINITY_POOR} kcal/mol. "
                    "Binding too weak for biological significance. "
                    "Check interactions via visualization."
                )
            else:
                concern = ""

            row_data = [
                row_idx - 2,
                data["lig_id"],
                data["orig_name"],
                flag,
                aff,
                data["affinity_tier"],
                concern,
                data["lig_dirname"],
                data["smiles"],
            ]

            for col_idx, value in enumerate(row_data, start=1):
                cell = ws3.cell(row=row_idx, column=col_idx, value=value)
                cell.border = thin_border

                if col_idx == 4:
                    if "POSITIVE" in str(value):
                        cell.fill = repulsive_fill
                        cell.font = Font(bold=True)
                    elif "POOR" in str(value):
                        cell.fill = poor_fill
                        cell.font = Font(bold=True)
    else:
        ws3.cell(
            row=3,
            column=1,
            value="✅ No flagged ligands — all pass thermodynamic sanity check.",
        ).font = Font(italic=True, color="548235", size=12)

    _set_col_widths(ws3, {
        "A": 5, "B": 12, "C": 35, "D": 22, "E": 13,
        "F": 12, "G": 60, "H": 35, "I": 60,
    })

    # ── SHEET 4: Statistics ──
    ws4 = wb.create_sheet(title="Statistics")

    total = len(results)
    done = sum(1 for r in results if r["status"] == STATUS_DONE)
    failed = sum(1 for r in results if r["status"] == STATUS_FAILED)
    pending = sum(1 for r in results if r["status"] == STATUS_PENDING)

    affinities = [
        r["best_affinity"]
        for r in results
        if r["best_affinity"] is not None
    ]

    n_poor = sum(1 for r in results if "POOR" in r.get("red_flag", ""))
    n_positive = sum(
        1 for r in results if "POSITIVE" in r.get("red_flag", "")
    )
    n_flagged = n_poor + n_positive

    # Count tiers
    tier_counts = {}
    for r in results:
        t = r.get("affinity_tier", "N/A")
        tier_counts[t] = tier_counts.get(t, 0) + 1

    n_timeout = sum(
        1
        for r in results
        if r["status"] == STATUS_FAILED
        and "TIMEOUT" in get_status_details(
            os.path.join(
                RESULTS_DIR, "ligands", r["lig_dirname"]
            )
        ).get("ERROR", "")
    )

    stats = [
        ("═══ PIPELINE OVERVIEW ═══", ""),
        ("Pipeline Version", "2.7.0 — Covalent Docking"),
        ("Docking Type", "Flexible Covalent (Vinardo)"),
        ("CNN Scoring", "DISABLED (not calibrated for covalent)"),
        ("Scoring Function", SCORING_MODE),
        ("Pose Sort Order", "minimizedAffinity (lower = better)"),
        ("Total Ligands", total),
        ("Completed", done),
        ("Failed", failed),
        ("  ↳ of which TIMEOUT", n_timeout),
        ("Pending", pending),
        ("GNINA Timeout Setting", f"{GNINA_TIMEOUT_SEC}s ({GNINA_TIMEOUT_SEC/60:.0f}min)"),
        ("", ""),
        ("═══ AFFINITY DISTRIBUTION (Vinardo) ═══", ""),
        (
            "Best Affinity (kcal/mol)",
            f"{min(affinities):.4f}" if affinities else "N/A",
        ),
        (
            "Worst Affinity (kcal/mol)",
            f"{max(affinities):.4f}" if affinities else "N/A",
        ),
        (
            "Mean Affinity (kcal/mol)",
            (
                f"{sum(affinities) / len(affinities):.4f}"
                if affinities
                else "N/A"
            ),
        ),
        (
            "Median Affinity (kcal/mol)",
            (
                f"{sorted(affinities)[len(affinities)//2]:.4f}"
                if affinities
                else "N/A"
            ),
        ),
        ("", ""),
        ("═══ AFFINITY TIERS ═══", ""),
        (
            f"🟢 EXCELLENT (< {AFFINITY_EXCELLENT})",
            tier_counts.get("EXCELLENT", 0),
        ),
        (
            f"🟢 GOOD ({AFFINITY_EXCELLENT} to {AFFINITY_GOOD})",
            tier_counts.get("GOOD", 0),
        ),
        (
            f"🟡 MODERATE ({AFFINITY_GOOD} to {AFFINITY_MODERATE})",
            tier_counts.get("MODERATE", 0),
        ),
        (
            f"🟠 MARGINAL ({AFFINITY_MODERATE} to {RED_FLAG_AFFINITY_POOR})",
            tier_counts.get("MARGINAL", 0),
        ),
        (
            f"🔴 POOR (≥ {RED_FLAG_AFFINITY_POOR})",
            tier_counts.get("POOR", 0),
        ),
        ("🔴 REPULSIVE (> 0)", tier_counts.get("REPULSIVE", 0)),
        ("", ""),
        ("═══ SANITY CHECK ═══", ""),
        ("🟡 Positive Affinity (repulsive)", n_positive),
        (f"🟠 Poor Affinity (≥ {RED_FLAG_AFFINITY_POOR})", n_poor),
        ("✅ Clean (no flags)", done - n_flagged),
        ("", ""),
        ("═══ THRESHOLDS & METHODOLOGY ═══", ""),
        (
            "Poor Affinity threshold",
            f"≥ {RED_FLAG_AFFINITY_POOR} kcal/mol",
        ),
        (
            "CNN metrics",
            "DISABLED — not calibrated for covalent docking",
        ),
        (
            "Ranking method",
            "Vinardo minimizedAffinity (classical force-field)",
        ),
        (
            "Scientific basis",
            "GNINA warning: CNN not calibrated for covalent",
        ),
    ]

    _write_header(ws4, 1, ["Metric", "Value"], header_fill_purple)

    for row_idx, (label, value) in enumerate(stats, start=2):
        label_cell = ws4.cell(row=row_idx, column=1, value=label)
        value_cell = ws4.cell(row=row_idx, column=2, value=value)
        label_cell.border = thin_border
        value_cell.border = thin_border

        if "═══" in str(label):
            label_cell.font = Font(bold=True, size=11, color="7030A0")
            label_cell.fill = PatternFill(
                start_color="E8E0F0", end_color="E8E0F0", fill_type="solid"
            )
            value_cell.fill = PatternFill(
                start_color="E8E0F0", end_color="E8E0F0", fill_type="solid"
            )
        else:
            label_cell.font = Font(bold=True)

    _set_col_widths(ws4, {"A": 42, "B": 52})

    excel_path = os.path.join(summary_dir, "docking_summary.xlsx")
    wb.save(excel_path)
    print(f"✔ Excel summary saved to {excel_path}")

    return excel_path


# ================================================================
# [v2.7.0] CSV Fallback — Covalent docking format
# ================================================================
def _generate_csv_fallback_v27(results_sorted: list, summary_dir: str):
    """
    [v2.7.0] CSV output — no CNN columns, affinity-ranked.
    
    FAIR compliance:
      F: Standardized column names, persistent IDs
      A: Plain CSV — universally readable
      I: SMILES for chemical identity, kcal/mol units
      R: Includes tier, flag, timing for full reproducibility
    """
    csv_path = os.path.join(summary_dir, "docking_summary.csv")

    with open(csv_path, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)

        # [v2.7.0] No CNN columns
        headers = [
            "Rank",
            "ID",
            "Original_Name",
            "Status",
            "Affinity_kcal_mol",
            "Affinity_Tier",
            "Sanity_Flag",
            "Elapsed_Min",
            "Dir_Name",
            "SMILES",
        ]
        writer.writerow(headers)

        for rank, data in enumerate(results_sorted, start=1):
            writer.writerow([
                rank,
                data["lig_id"],
                data["orig_name"],
                data["status"],
                data.get("best_affinity", ""),
                data.get("affinity_tier", ""),
                data.get("red_flag", ""),
                data["elapsed_min"],
                data["lig_dirname"],
                data["smiles"],
            ])

    print(f"✔ CSV summary saved to {csv_path}")
    return csv_path


# =========================
# PROGRESS TRACKING — giữ nguyên
# =========================
def update_progress_csv(ligands: list, summary_dir: str):
    progress_file = os.path.join(summary_dir, "progress.csv")

    with open(progress_file, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(
            [
                "ID",
                "DIR_NAME",
                "STATUS",
                "ELAPSED_MIN",
                "BEST_AFFINITY",
                "START_TIME",
                "END_TIME",
            ]
        )

        for lig in ligands:
            details = get_status_details(lig["lig_root"])
            writer.writerow(
                [
                    lig["lig_id"],
                    lig["lig_dirname"],
                    details.get("STATUS", STATUS_PENDING),
                    details.get("ELAPSED_MIN", ""),
                    details.get("BEST_AFFINITY", ""),
                    details.get("START_TIME", ""),
                    details.get("END_TIME", ""),
                ]
            )


def print_progress_summary(
    finished: list, failed: list, skipped: list, total: int
):
    done = len(finished) + len(skipped)
    pct = (done / total) * 100 if total > 0 else 0
    print(
        f"\n📊 Progress: {done}/{total} ({pct:.1f}%)"
        f" | ✅ {len(finished)} new"
        f" | ⏭️ {len(skipped)} skipped"
        f" | ❌ {len(failed)} failed"
    )


# =========================
# MAIN PIPELINE
# =========================
def main():
    print("=" * 60)
    print("🧬 GNINA Flexible COVALENT Docking Pipeline v2.7.0")
    print("   Scoring: Vinardo (CNN disabled — not calibrated)")
    print("   Ranking: minimizedAffinity (lower = better)")
    print("=" * 60)
    print(f"📂 Base dir:     {BASE_DIR}")
    print(f"🔬 GNINA bin:    {GNINA_BIN}")
    print(f"🧪 Protein:      {PROTEIN_PATH}")
    print(f"📎 Ref ligand:   {REF_LIGAND}")
    print(f"📦 Ligands:      {LIGAND_SDF}")
    print(f"🎯 Flex res:     {FLEX_RESIDUES}")
    print(f"📊 Scoring:      {SCORING_MODE} (--cnn_scoring {CNN_SCORING})")
    print(f"📊 Pose sort:    affinity (--pose_sort_order {POSE_SORT_ORDER})")
    print(f"🚩 Sanity:       Affinity≥{RED_FLAG_AFFINITY_POOR} or >0 → flagged")
    print(f"⏱️  Timeout:      {GNINA_TIMEOUT_SEC}s ({GNINA_TIMEOUT_SEC/60:.0f} min/ligand)")
    print("=" * 60)

    if not os.path.isfile(GNINA_BIN):
        print(f"❌ GNINA binary not found: {GNINA_BIN}")
        return
    if not os.access(GNINA_BIN, os.X_OK):
        print(f"❌ GNINA binary not executable: {GNINA_BIN}")
        return

    prepare_root_folders()

    ligands = split_ligands(
        LIGAND_SDF, ligands_root=f"{RESULTS_DIR}/ligands"
    )

    total = len(ligands)
    if total == 0:
        print("❌ No valid ligands found!")
        return

    finished = []
    failed = []
    skipped = []

    print(f"\n🚀 Starting batch covalent docking: {total} ligands\n")
    start_all = time.time()

    for idx, lig in enumerate(ligands, start=1):
        lig_id = lig["lig_id"]
        lig_root = lig["lig_root"]

        status = read_status(lig_root)

        if status == STATUS_DONE:
            print(f"⏭️ [{idx}/{total}] {lig_id} already DONE — skipping")
            skipped.append(lig_id)
            continue

        if status == STATUS_RUNNING:
            print(
                f"⚠️ [{idx}/{total}] {lig_id} was RUNNING (incomplete)"
                " — retrying"
            )

        success = run_gnina(lig, idx, total)

        if success:
            finished.append(lig_id)
        else:
            failed.append(lig_id)

        if idx % 10 == 0:
            update_progress_csv(ligands, f"{RESULTS_DIR}/summary")
            print_progress_summary(finished, failed, skipped, total)

    elapsed_all = (time.time() - start_all) / 60

    summary_dir = f"{RESULTS_DIR}/summary"

    with open(f"{summary_dir}/finished_ligands.txt", "w") as f:
        f.write("\n".join(finished + skipped))

    with open(f"{summary_dir}/failed_ligands.txt", "w") as f:
        f.write("\n".join(failed))

    update_progress_csv(ligands, summary_dir)

    print("\n📊 Generating Excel summary (Vinardo affinity ranking)...")
    generate_excel_summary(ligands, summary_dir)

    print("\n" + "=" * 60)
    print("📊 FINAL SUMMARY — Covalent Docking v2.7.0")
    print("=" * 60)
    print(f"⏱️  Total time: {elapsed_all:.2f} min")
    print(f"✅ Completed (new): {len(finished)}")
    print(f"⏭️  Skipped (cached): {len(skipped)}")
    print(f"❌ Failed: {len(failed)}")
    print(f"📊 Scoring: {SCORING_MODE} (CNN disabled)")
    print(f"📁 Results: {RESULTS_DIR}")
    print(f"📊 Excel: {summary_dir}/docking_summary.xlsx")
    print(f"📊 CSV:   {summary_dir}/docking_summary.csv")
    print("=" * 60)

    if failed:
        print(
            f"\n⚠️ Failed ligands: {', '.join(failed[:10])}"
            + (
                f"... and {len(failed)-10} more"
                if len(failed) > 10
                else ""
            )
        )


# =============================================================
# Entry point
# =============================================================
if __name__ == "__main__":
    main()
