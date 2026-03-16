# ============================================================
# 🧬 GNINA Flexible Docking Pipeline v2.5 — Linux + Fish Shell
# ============================================================
# Thay đổi so với code gốc Kaggle:
#   1. PATH lấy từ .env thay vì hardcode
#   2. Đảm bảo LD_LIBRARY_PATH cho conda env
#   3. subprocess env đảm bảo gnina tìm được thư viện
#   TOÀN BỘ LOGIC DOCKING + SCORING + EXCEL GIỮ NGUYÊN 100%
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
# [MỚI] .env loading + LD_LIBRARY_PATH đảm bảo
# =============================================================
NOTEBOOK_DIR = Path.cwd()

# Đảm bảo LD_LIBRARY_PATH luôn đúng
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
    """Lấy path từ env var, fallback về giá trị mặc định."""
    raw = os.environ.get(env_var, fallback)
    expanded = os.path.expanduser(raw)
    if not os.path.isabs(expanded):
        expanded = str(NOTEBOOK_DIR / expanded)
    return os.path.normpath(expanded)


# =========================
# GLOBAL CONFIG — [MỚI] path từ .env, logic giữ nguyên
# =========================
BASE_DIR = _resolve_path("DOCKING_BASE_DIR", str(NOTEBOOK_DIR))
RESULTS_DIR = f"{BASE_DIR}/docking_results"
PROTEIN_PATH = _resolve_path(
    "PROTEIN_PATH",
    f"{BASE_DIR}/data/mmp2_7xjo_ready_for_gnina.pdb",
)
REF_LIGAND = _resolve_path("REF_LIGAND", f"{BASE_DIR}/data/ref_ligand.sdf")
LIGAND_SDF = _resolve_path("LIGAND_SDF", f"{BASE_DIR}/data/ligands_prepared.sdf")
FLEX_RESIDUES = "A:182,A:181,A:215,A:262,A:49"
SEED = "42"
GPU_DEVICE = os.environ.get("CUDA_VISIBLE_DEVICES", "0")

# GNINA_BIN — [MỚI] auto-detect
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
# [MỚI] Subprocess env — đảm bảo gnina tìm được thư viện
# =========================
def _get_subprocess_env() -> dict:
    """Tạo env dict cho subprocess, đảm bảo LD_LIBRARY_PATH có conda lib."""
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
# STATUS MANAGEMENT — 100% giữ nguyên code gốc
# =========================
def write_status(lig_root: str, status: str, **kwargs):
    """
    Write status to STATUS.txt.
    For RUNNING: overwrites file (new run)
    For DONE/FAILED: appends to preserve history
    """
    status_file = os.path.join(lig_root, "STATUS.txt")

    if status == STATUS_RUNNING:
        with open(status_file, "w") as f:
            f.write(f"STATUS={status}\n")
            f.write(f"HOST={socket.gethostname()}\n")
            f.write(f"START_TIME={time.strftime('%Y-%m-%d %H:%M:%S')}\n")
    else:
        with open(status_file, "a") as f:
            f.write(f"STATUS={status}\n")
            for key, value in kwargs.items():
                f.write(f"{key.upper()}={value}\n")
            f.write(f"END_TIME={time.strftime('%Y-%m-%d %H:%M:%S')}\n")


def read_status(lig_root: str) -> str:
    """Read the LAST status from STATUS.txt."""
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
    """Parse full status file into dictionary"""
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
# UTILS — 100% giữ nguyên code gốc
# =========================
def sanitize_name(name: str, max_len: int = 80) -> str:
    """Make a filesystem-safe ligand name."""
    if not name:
        return "NA"

    name = name.strip().replace(" ", "_")
    name = re.sub(r'[\/:*?"<>|\\]', "", name)
    name = re.sub(r"_+", "_", name)
    name = name.strip("_")

    return name[:max_len] if name else "NA"


# =========================
# PREPARE ROOT FOLDERS — 100% giữ nguyên code gốc
# =========================
def prepare_root_folders():
    """Create directory structure and copy reference files"""
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


# =========================
# SPLIT LIGANDS — 100% giữ nguyên code gốc
# =========================
def split_ligands(input_sdf: str, ligands_root: str) -> list:
    """Split multi-ligand SDF into per-ligand folders."""
    os.makedirs(ligands_root, exist_ok=True)

    suppl = Chem.SDMolSupplier(input_sdf, removeHs=False)
    ligands = []
    mapping = []

    for idx, mol in enumerate(suppl, start=1):
        if mol is None:
            print(f"⚠️ Skipping invalid molecule at index {idx}")
            continue

        lig_id = f"LIG_{idx:04d}"
        orig_name = mol.GetProp("_Name") if mol.HasProp("_Name") else "NA"
        safe_name = sanitize_name(orig_name)

        lig_dirname = f"{lig_id}__{safe_name}"
        lig_root = os.path.join(ligands_root, lig_dirname)
        input_dir = os.path.join(lig_root, "input")
        os.makedirs(input_dir, exist_ok=True)

        ligand_sdf = os.path.join(input_dir, "ligand.sdf")
        writer = Chem.SDWriter(ligand_sdf)
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

    print(f"✔ Split {len(ligands)} ligands")
    print(f"✔ Mapping written to {mapping_file}")

    return ligands


# =========================
# SCORE PARSING — 100% giữ nguyên code gốc
# =========================
def parse_top_scores(sdf_path: str, n: int = 3) -> list:
    """
    Extract top N poses with all relevant scores from docked SDF.

    Returns: List of dicts with score info, sorted by minimizedAffinity
             (best first = most negative)
    """
    try:
        suppl = Chem.SDMolSupplier(sdf_path, removeHs=False)
        poses = []

        for pose_idx, mol in enumerate(suppl, start=1):
            if mol is None:
                continue

            pose_data = {"pose_rank": pose_idx}

            # Extract all relevant properties
            score_props = [
                "minimizedAffinity",
                "CNNscore",
                "CNNaffinity",
                "CNN_VS",
                "vinardo",
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

        # Sort by minimizedAffinity (lower/more negative = better binding)
        poses.sort(
            key=lambda x: (
                x.get("minimizedAffinity") is None,
                x.get("minimizedAffinity", 0),
            )
        )

        return poses[:n]

    except Exception as e:
        print(f"⚠️ Error parsing scores: {e}")
        return []


def parse_best_score(sdf_path: str) -> float:
    """Extract best minimizedAffinity from docked SDF"""
    top_scores = parse_top_scores(sdf_path, n=1)
    if top_scores and top_scores[0].get("minimizedAffinity") is not None:
        return top_scores[0]["minimizedAffinity"]
    return 0.0


# =========================
# RUN GNINA — 100% giữ nguyên logic, chỉ thêm env= cho subprocess
# =========================
def run_gnina(ligand_info: dict, idx: int, total: int) -> bool:
    """Run GNINA docking for a single ligand."""
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
        "32",
        "--cnn_scoring",
        "rescore",
        "--cnn_empirical_weight",
        "2.0",
        "--pose_sort_order",
        "CNNscore",
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

    with open(cmd_file, "w") as f:
        f.write(" \\\n    ".join(cmd))

    try:
        with open(stderr_file, "w") as stderr_f:
            subprocess.run(
                cmd,
                check=True,
                stdout=subprocess.DEVNULL,
                stderr=stderr_f,
                text=True,
                env=_get_subprocess_env(),  # ← [MỚI] đảm bảo gnina tìm được thư viện
            )

        if not os.path.exists(out_lig) or os.path.getsize(out_lig) == 0:
            raise RuntimeError("GNINA produced empty output SDF")

        elapsed = (time.time() - start) / 60
        best_score = parse_best_score(out_lig)

        write_status(
            lig_root,
            STATUS_DONE,
            elapsed_min=f"{elapsed:.2f}",
            best_affinity=f"{best_score:.4f}" if best_score else "NA",
        )

        print(
            f"✅ [{idx}/{total}] {lig_id} DONE in {elapsed:.2f} min"
            f" (affinity: {best_score:.4f})"
        )
        return True

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


# =========================
# EXCEL SUMMARY GENERATION — 100% giữ nguyên code gốc
# =========================
def generate_excel_summary(ligands: list, summary_dir: str):
    """
    Generate comprehensive Excel summary with:
    - Sheet 1: All ligands with top 3 scores
    - Sheet 2: Top 50 hits ranked by best affinity
    - Sheet 3: Statistics
    """
    try:
        from openpyxl import Workbook
        from openpyxl.styles import Font, PatternFill, Alignment, Border, Side

        USE_OPENPYXL = True
    except ImportError:
        print("⚠️ openpyxl not installed, generating CSV instead")
        USE_OPENPYXL = False

    # Collect all results
    results = []
    for lig in ligands:
        lig_root = lig["lig_root"]
        out_sdf = os.path.join(lig_root, "output", "docked.sdf")

        status = read_status(lig_root)
        details = get_status_details(lig_root)

        top_scores = []
        if status == STATUS_DONE and os.path.exists(out_sdf):
            top_scores = parse_top_scores(out_sdf, n=3)

        results.append(
            {
                "lig_id": lig["lig_id"],
                "lig_dirname": lig["lig_dirname"],
                "orig_name": lig["orig_name"],
                "smiles": lig["smiles"],
                "status": status,
                "elapsed_min": details.get("ELAPSED_MIN", ""),
                "top_scores": top_scores,
                "best_affinity": (
                    top_scores[0].get("minimizedAffinity")
                    if top_scores
                    else None
                ),
            }
        )

    # Sort by best affinity (lower = better)
    results_sorted = sorted(
        results,
        key=lambda x: (
            x["best_affinity"] is None,
            x["best_affinity"] if x["best_affinity"] else 999,
        ),
    )

    if USE_OPENPYXL:
        _generate_xlsx(results, results_sorted, summary_dir)
    else:
        _generate_csv_fallback(results_sorted, summary_dir)


def _generate_xlsx(results: list, results_sorted: list, summary_dir: str):
    """Generate Excel file with openpyxl"""
    from openpyxl import Workbook
    from openpyxl.styles import Font, PatternFill, Alignment, Border, Side

    wb = Workbook()

    # ========== SHEET 1: All Results ==========
    ws1 = wb.active
    ws1.title = "All_Docking_Results"

    headers = [
        "Rank",
        "ID",
        "Original_Name",
        "Status",
        "Elapsed_Min",
        "Score_1",
        "Pose_1",
        "CNNscore_1",
        "CNNaffinity_1",
        "Score_2",
        "Pose_2",
        "CNNscore_2",
        "CNNaffinity_2",
        "Score_3",
        "Pose_3",
        "CNNscore_3",
        "CNNaffinity_3",
        "Best_Affinity",
        "SMILES",
    ]

    # Style settings
    header_font = Font(bold=True, color="FFFFFF", size=11)
    header_fill = PatternFill(
        start_color="2F5496", end_color="2F5496", fill_type="solid"
    )
    good_fill = PatternFill(
        start_color="C6EFCE", end_color="C6EFCE", fill_type="solid"
    )
    bad_fill = PatternFill(
        start_color="FFC7CE", end_color="FFC7CE", fill_type="solid"
    )
    thin_border = Border(
        left=Side(style="thin"),
        right=Side(style="thin"),
        top=Side(style="thin"),
        bottom=Side(style="thin"),
    )

    # Write headers
    for col, header in enumerate(headers, start=1):
        cell = ws1.cell(row=1, column=col, value=header)
        cell.font = header_font
        cell.fill = header_fill
        cell.alignment = Alignment(horizontal="center", vertical="center")
        cell.border = thin_border

    # Freeze header row
    ws1.freeze_panes = "A2"

    # Write data
    for row_idx, data in enumerate(results_sorted, start=2):
        rank = row_idx - 1
        top_scores = data["top_scores"]

        row_data = [
            rank,
            data["lig_id"],
            data["orig_name"],
            data["status"],
            data["elapsed_min"],
        ]

        # Add top 3 scores
        for i in range(3):
            if i < len(top_scores):
                score_data = top_scores[i]
                row_data.extend(
                    [
                        score_data.get("minimizedAffinity", ""),
                        score_data.get("pose_rank", ""),
                        score_data.get("CNNscore", ""),
                        score_data.get("CNNaffinity", ""),
                    ]
                )
            else:
                row_data.extend(["", "", "", ""])

        row_data.append(
            data["best_affinity"] if data["best_affinity"] else ""
        )
        row_data.append(data["smiles"])

        for col_idx, value in enumerate(row_data, start=1):
            cell = ws1.cell(row=row_idx, column=col_idx, value=value)
            cell.border = thin_border

            # Highlight status
            if col_idx == 4:  # Status column
                if value == STATUS_DONE:
                    cell.fill = good_fill
                elif value == STATUS_FAILED:
                    cell.fill = bad_fill

    # Adjust column widths
    column_widths = {
        "A": 6,
        "B": 12,
        "C": 35,
        "D": 10,
        "E": 12,
        "F": 12,
        "G": 8,
        "H": 12,
        "I": 12,
        "J": 12,
        "K": 8,
        "L": 12,
        "M": 12,
        "N": 12,
        "O": 8,
        "P": 12,
        "Q": 12,
        "R": 14,
        "S": 60,
    }
    for col, width in column_widths.items():
        ws1.column_dimensions[col].width = width

    # ========== SHEET 2: Top Hits ==========
    ws2 = wb.create_sheet(title="Top_50_Hits")

    top_hits_headers = [
        "Rank",
        "ID",
        "Original_Name",
        "Best_Affinity",
        "CNNscore",
        "CNNaffinity",
        "SMILES",
        "Notes",
    ]

    for col, header in enumerate(top_hits_headers, start=1):
        cell = ws2.cell(row=1, column=col, value=header)
        cell.font = header_font
        cell.fill = PatternFill(
            start_color="548235", end_color="548235", fill_type="solid"
        )
        cell.alignment = Alignment(horizontal="center", vertical="center")
        cell.border = thin_border

    ws2.freeze_panes = "A2"

    # Write top 50 hits
    top_50 = [r for r in results_sorted if r["status"] == STATUS_DONE][:50]

    for row_idx, data in enumerate(top_50, start=2):
        rank = row_idx - 1
        best_pose = data["top_scores"][0] if data["top_scores"] else {}

        row_data = [
            rank,
            data["lig_id"],
            data["orig_name"],
            best_pose.get("minimizedAffinity", ""),
            best_pose.get("CNNscore", ""),
            best_pose.get("CNNaffinity", ""),
            data["smiles"],
            "",  # Notes column for manual annotation
        ]

        for col_idx, value in enumerate(row_data, start=1):
            cell = ws2.cell(row=row_idx, column=col_idx, value=value)
            cell.border = thin_border

            # Highlight best affinity
            if col_idx == 4 and isinstance(value, (int, float)):
                if value < -8.0:
                    cell.fill = PatternFill(
                        start_color="92D050",
                        end_color="92D050",
                        fill_type="solid",
                    )
                elif value < -6.0:
                    cell.fill = PatternFill(
                        start_color="FFEB9C",
                        end_color="FFEB9C",
                        fill_type="solid",
                    )

    ws2.column_dimensions["A"].width = 6
    ws2.column_dimensions["B"].width = 12
    ws2.column_dimensions["C"].width = 35
    ws2.column_dimensions["D"].width = 14
    ws2.column_dimensions["E"].width = 12
    ws2.column_dimensions["F"].width = 12
    ws2.column_dimensions["G"].width = 60
    ws2.column_dimensions["H"].width = 30

    # ========== SHEET 3: Statistics ==========
    ws3 = wb.create_sheet(title="Statistics")

    total = len(results)
    done = sum(1 for r in results if r["status"] == STATUS_DONE)
    failed = sum(1 for r in results if r["status"] == STATUS_FAILED)
    pending = sum(1 for r in results if r["status"] == STATUS_PENDING)

    done_results = [r for r in results if r["best_affinity"] is not None]
    affinities = [r["best_affinity"] for r in done_results]

    stats = [
        ("Total Ligands", total),
        ("Completed", done),
        ("Failed", failed),
        ("Pending", pending),
        ("", ""),
        ("Best Affinity", min(affinities) if affinities else "N/A"),
        ("Worst Affinity", max(affinities) if affinities else "N/A"),
        (
            "Average Affinity",
            sum(affinities) / len(affinities) if affinities else "N/A",
        ),
        ("", ""),
        (
            "Ligands < -8.0 kcal/mol",
            sum(1 for a in affinities if a < -8.0),
        ),
        (
            "Ligands < -7.0 kcal/mol",
            sum(1 for a in affinities if a < -7.0),
        ),
        (
            "Ligands < -6.0 kcal/mol",
            sum(1 for a in affinities if a < -6.0),
        ),
    ]

    for row_idx, (label, value) in enumerate(stats, start=1):
        ws3.cell(row=row_idx, column=1, value=label).font = Font(bold=True)
        ws3.cell(row=row_idx, column=2, value=value)

    ws3.column_dimensions["A"].width = 25
    ws3.column_dimensions["B"].width = 15

    # Save workbook
    excel_path = os.path.join(summary_dir, "docking_summary.xlsx")
    wb.save(excel_path)
    print(f"✔ Excel summary saved to {excel_path}")

    return excel_path


def _generate_csv_fallback(results_sorted: list, summary_dir: str):
    """Fallback to CSV if openpyxl not available"""
    csv_path = os.path.join(summary_dir, "docking_summary.csv")

    with open(csv_path, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)

        headers = [
            "Rank",
            "ID",
            "Original_Name",
            "Status",
            "Elapsed_Min",
            "Score_1",
            "Pose_1",
            "CNNscore_1",
            "Score_2",
            "Pose_2",
            "CNNscore_2",
            "Score_3",
            "Pose_3",
            "CNNscore_3",
            "Best_Affinity",
            "SMILES",
        ]
        writer.writerow(headers)

        for rank, data in enumerate(results_sorted, start=1):
            top_scores = data["top_scores"]

            row = [
                rank,
                data["lig_id"],
                data["orig_name"],
                data["status"],
                data["elapsed_min"],
            ]

            for i in range(3):
                if i < len(top_scores):
                    row.extend(
                        [
                            top_scores[i].get("minimizedAffinity", ""),
                            top_scores[i].get("pose_rank", ""),
                            top_scores[i].get("CNNscore", ""),
                        ]
                    )
                else:
                    row.extend(["", "", ""])

            row.append(
                data["best_affinity"] if data["best_affinity"] else ""
            )
            row.append(data["smiles"])

            writer.writerow(row)

    print(f"✔ CSV summary saved to {csv_path}")
    return csv_path


# =========================
# PROGRESS TRACKING — 100% giữ nguyên code gốc
# =========================
def update_progress_csv(ligands: list, summary_dir: str):
    """Generate progress.csv with current status of all ligands"""
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
    """Print current progress"""
    done = len(finished) + len(skipped)
    pct = (done / total) * 100 if total > 0 else 0
    print(
        f"\n📊 Progress: {done}/{total} ({pct:.1f}%)"
        f" | ✅ {len(finished)} new"
        f" | ⏭️ {len(skipped)} skipped"
        f" | ❌ {len(failed)} failed"
    )


# =========================
# MAIN PIPELINE — 100% giữ nguyên logic, chỉ thêm info print
# =========================
def main():
    print("=" * 60)
    print("🧬 GNINA Flexible Docking Pipeline v2.5")
    print("=" * 60)
    print(f"📂 Base dir:  {BASE_DIR}")
    print(f"🔬 GNINA bin: {GNINA_BIN}")
    print(f"🧪 Protein:   {PROTEIN_PATH}")
    print(f"📎 Ref lig:   {REF_LIGAND}")
    print(f"📦 Ligands:   {LIGAND_SDF}")
    print(f"🎯 Flex res:  {FLEX_RESIDUES}")
    print("=" * 60)

    # Verify gnina binary exists
    if not os.path.isfile(GNINA_BIN):
        print(f"❌ GNINA binary not found: {GNINA_BIN}")
        return
    if not os.access(GNINA_BIN, os.X_OK):
        print(f"❌ GNINA binary not executable: {GNINA_BIN}")
        return

    # Setup
    prepare_root_folders()

    # Split ligands
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

    print(f"\n🚀 Starting batch docking: {total} ligands\n")
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

    # Final summary
    elapsed_all = (time.time() - start_all) / 60

    summary_dir = f"{RESULTS_DIR}/summary"

    with open(f"{summary_dir}/finished_ligands.txt", "w") as f:
        f.write("\n".join(finished + skipped))

    with open(f"{summary_dir}/failed_ligands.txt", "w") as f:
        f.write("\n".join(failed))

    update_progress_csv(ligands, summary_dir)

    # ========== GENERATE EXCEL SUMMARY ==========
    print("\n📊 Generating Excel summary...")
    generate_excel_summary(ligands, summary_dir)

    print("\n" + "=" * 60)
    print("📊 FINAL SUMMARY")
    print("=" * 60)
    print(f"⏱️  Total time: {elapsed_all:.2f} min")
    print(f"✅ Completed (new): {len(finished)}")
    print(f"⏭️  Skipped (cached): {len(skipped)}")
    print(f"❌ Failed: {len(failed)}")
    print(f"📁 Results: {RESULTS_DIR}")
    print(f"📊 Excel: {summary_dir}/docking_summary.xlsx")
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
