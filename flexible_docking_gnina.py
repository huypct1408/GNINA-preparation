# Code này tối ưu hóa của code trên ở chỗ thêm vào task tracking (biết được khi nào chất gì chạy xong (mất bao lâu) và đang làm tới đâu)
import os
import subprocess
import time
from rdkit import Chem
import socket

# =========================
# GLOBAL CONFIG
# =========================
BASE_DIR = "/kaggle/working"
RESULTS_DIR = f"{BASE_DIR}/docking_results/8skl"

GNINA_BIN = "/kaggle/working/gnina"

PROTEIN_PATH = "/kaggle/input/docking-profile/8skl_receptor_prepared_for_gnina.pdb"
REF_LIGAND   = "/kaggle/input/docking-profile/v2o_ligand_8skl.sdf"
LIGAND_SDF   = "/kaggle/input/docking-profile/ligands_for_8skl_prepared.sdf"

FLEX_RESIDUES = "A:182,A:181,A:215,A:262,A:49"

SEED = "42"
GPU_DEVICE = "0"

# =========================
# UTILS
# =========================
def write_status(lig_root, status):
    with open(os.path.join(lig_root, "STATUS.txt"), "w") as f:
        f.write(status)

def read_status(lig_root):
    status_file = os.path.join(lig_root, "STATUS.txt")
    if os.path.exists(status_file):
        return open(status_file).read().strip()
    return None

# =========================
# PREPARE ROOT FOLDERS
# =========================
def prepare_root_folders():
    os.makedirs(f"{RESULTS_DIR}/protein", exist_ok=True)
    os.makedirs(f"{RESULTS_DIR}/reference", exist_ok=True)
    os.makedirs(f"{RESULTS_DIR}/ligands", exist_ok=True)
    os.makedirs(f"{RESULTS_DIR}/summary", exist_ok=True)

    os.system(f"cp {PROTEIN_PATH} {RESULTS_DIR}/protein/receptor.pdb")
    os.system(f"cp {REF_LIGAND} {RESULTS_DIR}/reference/ref_ligand.sdf")

# =========================
# SPLIT LIGANDS (FILESYSTEM ONLY)
# =========================
def split_ligands(input_sdf, ligands_root):
    suppl = Chem.SDMolSupplier(input_sdf, removeHs=False)
    ligands = []

    for i, mol in enumerate(suppl):
        if mol is None:
            continue

        lig_id = mol.GetProp("_Name") if mol.HasProp("_Name") else f"LIG{i:04d}"
        lig_root = os.path.join(ligands_root, lig_id)
        input_dir = os.path.join(lig_root, "input")

        os.makedirs(input_dir, exist_ok=True)

        ligand_sdf = os.path.join(input_dir, "ligand.sdf")
        writer = Chem.SDWriter(ligand_sdf)
        writer.write(mol)
        writer.close()

        ligands.append((lig_id, ligand_sdf))

    return ligands

# =========================
# RUN GNINA FOR ONE LIGAND
# =========================
def run_gnina(lig_id, ligand_sdf, idx, total):
    lig_root = f"{RESULTS_DIR}/ligands/{lig_id}"

    out_dir = f"{lig_root}/output"
    log_dir = f"{lig_root}/logs"

    os.makedirs(out_dir, exist_ok=True)
    os.makedirs(log_dir, exist_ok=True)

    out_lig  = f"{out_dir}/docked.sdf"
    out_flex = f"{out_dir}/flex_residues.pdb"
    log_file = f"{log_dir}/gnina.log"
    cmd_file = f"{log_dir}/command.txt"

    # --- STATUS: RUNNING ---
    with open(os.path.join(lig_root, "STATUS.txt"), "w") as f:
        f.write("RUNNING\n")
        f.write(f"HOST={socket.gethostname()}\n")
        f.write(f"START_TIME={time.strftime('%Y-%m-%d %H:%M:%S')}\n")

    print(f"\n🔄 [{idx}/{total}] Docking {lig_id} ...")

    start = time.time()

    cmd = [
        GNINA_BIN,
        "-r", f"{RESULTS_DIR}/protein/receptor.pdb",
        "-l", ligand_sdf,
        "--autobox_ligand", f"{RESULTS_DIR}/reference/ref_ligand.sdf",
        "--autobox_add", "10",
        "--flexres", FLEX_RESIDUES,
        "--num_modes", "10",
        "--exhaustiveness", "32",
        "--cnn_scoring", "rescore",
        "--cnn_empirical_weight", "2.0",
        "--pose_sort_order", "CNNscore",
        "--device", GPU_DEVICE,
        "--seed", SEED,
        "--atom_term_data",
        "-o", out_lig,
        "--out_flex", out_flex,
        "--log", log_file,
    ]

    # --- OPTIONAL 1: SAVE COMMAND ---
    with open(cmd_file, "w") as f:
        f.write(" ".join(cmd))

    # --- RUN GNINA ---
    subprocess.run(
        cmd,
        check=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True
    )

    # --- OPTIONAL 3: SANITY CHECK ---
    if not os.path.exists(out_lig) or os.path.getsize(out_lig) == 0:
        raise RuntimeError("GNINA produced empty output SDF")

    elapsed = (time.time() - start) / 60

    # --- STATUS: DONE ---
    with open(os.path.join(lig_root, "STATUS.txt"), "a") as f:
        f.write(f"DONE\n")
        f.write(f"ELAPSED_MIN={elapsed:.2f}\n")
        f.write(f"END_TIME={time.strftime('%Y-%m-%d %H:%M:%S')}\n")

    print(f"✅ [{idx}/{total}] {lig_id} DONE in {elapsed:.2f} min")

# =========================
# MAIN PIPELINE
# =========================
def main():
    prepare_root_folders()

    ligands = split_ligands(
        LIGAND_SDF,
        ligands_root=f"{RESULTS_DIR}/ligands"
    )

    total = len(ligands)
    finished = []
    failed = []

    print(f"\n🚀 GNINA batch docking started: {total} ligands\n")

    start_all = time.time()

    for idx, (lig_id, lig_sdf) in enumerate(ligands, start=1):
        lig_root = f"{RESULTS_DIR}/ligands/{lig_id}"

        # Resume-safe
        status = read_status(lig_root)
        if status and status.startswith("DONE"):
            print(f"⏭️ [{idx}/{total}] {lig_id} already DONE — skipping")
            finished.append(lig_id)
            continue

        try:
            run_gnina(lig_id, lig_sdf, idx, total)
            finished.append(lig_id)

        except Exception:
            write_status(lig_root, "FAILED")
            print(f"❌ [{idx}/{total}] {lig_id} FAILED")
            failed.append(lig_id)

    with open(f"{RESULTS_DIR}/summary/finished_ligands.txt", "w") as f:
        f.write("\n".join(finished))

    with open(f"{RESULTS_DIR}/summary/failed_ligands.txt", "w") as f:
        f.write("\n".join(failed))

    elapsed_all = (time.time() - start_all) / 60

    print("\n==============================")
    print(f"⏱️ Total time: {elapsed_all:.2f} min")
    print(f"✅ Finished: {len(finished)}")
    print(f"❌ Failed:   {len(failed)}")

# =========================
# ENTRY POINT
# =========================
if __name__ == "__main__":
    main()
