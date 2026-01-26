#Code tối ưu hóa version 2.0 (Dưa trên feedback của anh An)
#!/usr/bin/env python3
"""
flexible_docking_execution.py
Optimized GNINA flexible docking pipeline with:
- Task tracking (status, timing, progress)
- Resume capability  
- Robust error handling
"""

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

# =========================
# GLOBAL CONFIG
# =========================
BASE_DIR = "/kaggle/working"
RESULTS_DIR = f"{BASE_DIR}/docking_results/8skl"
GNINA_BIN = "/kaggle/working/gnina"
PROTEIN_PATH = "/kaggle/input/docking-profile/protein_8skl_protonated_chimera.pdb"
REF_LIGAND = "/kaggle/input/docking-profile/v2o_ligand_8skl.sdf"
LIGAND_SDF = "/kaggle/input/docking-profile/ligands_for_8skl_prepared_v2.0.sdf"
FLEX_RESIDUES = "A:182,A:181,A:215,A:262,A:49"
SEED = "42"
GPU_DEVICE = "0"

# Status constants
STATUS_PENDING = "PENDING"
STATUS_RUNNING = "RUNNING"
STATUS_DONE = "DONE"
STATUS_FAILED = "FAILED"

# =========================
# STATUS MANAGEMENT (FIXED)
# =========================
def write_status(lig_root: str, status: str, **kwargs):
    """
    Write status to STATUS.txt.
    For RUNNING: overwrites file (new run)
    For DONE/FAILED: appends to preserve history
    """
    status_file = os.path.join(lig_root, "STATUS.txt")
    
    if status == STATUS_RUNNING:
        # New run - overwrite
        with open(status_file, "w") as f:
            f.write(f"STATUS={status}\n")
            f.write(f"HOST={socket.gethostname()}\n")
            f.write(f"START_TIME={time.strftime('%Y-%m-%d %H:%M:%S')}\n")
    else:
        # Completion - append
        with open(status_file, "a") as f:
            f.write(f"STATUS={status}\n")
            for key, value in kwargs.items():
                f.write(f"{key.upper()}={value}\n")
            f.write(f"END_TIME={time.strftime('%Y-%m-%d %H:%M:%S')}\n")


def read_status(lig_root: str) -> str:
    """
    Read the LAST status from STATUS.txt.
    Returns: PENDING, RUNNING, DONE, or FAILED
    """
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
# UTILS
# =========================
def sanitize_name(name: str, max_len: int = 80) -> str:
    """Make a filesystem-safe ligand name."""
    if not name:
        return "NA"
    
    name = name.strip().replace(" ", "_")
    name = re.sub(r'[\/:*?"<>|\\]', '', name)  # Added backslash
    name = re.sub(r'_+', '_', name)
    name = name.strip("_")
    
    return name[:max_len] if name else "NA"


# =========================
# PREPARE ROOT FOLDERS
# =========================
def prepare_root_folders():
    """Create directory structure and copy reference files"""
    dirs = [
        f"{RESULTS_DIR}/protein",
        f"{RESULTS_DIR}/reference", 
        f"{RESULTS_DIR}/ligands",
        f"{RESULTS_DIR}/summary"
    ]
    
    for d in dirs:
        os.makedirs(d, exist_ok=True)
    
    # Use shutil instead of os.system (safer, cross-platform)
    protein_dest = f"{RESULTS_DIR}/protein/receptor.pdb"
    ref_dest = f"{RESULTS_DIR}/reference/ref_ligand.sdf"
    
    if not os.path.exists(protein_dest):
        shutil.copy(PROTEIN_PATH, protein_dest)
        print(f"✔ Copied receptor to {protein_dest}")
    
    if not os.path.exists(ref_dest):
        shutil.copy(REF_LIGAND, ref_dest)
        print(f"✔ Copied reference ligand to {ref_dest}")


# =========================
# SPLIT LIGANDS
# =========================
def split_ligands(input_sdf: str, ligands_root: str) -> list:
    """
    Split multi-ligand SDF into per-ligand folders.
    
    Returns: List of tuples (lig_id, lig_dirname, lig_root, ligand_sdf_path)
    """
    os.makedirs(ligands_root, exist_ok=True)
    
    suppl = Chem.SDMolSupplier(input_sdf, removeHs=False)
    ligands = []
    mapping = []
    
    for idx, mol in enumerate(suppl, start=1):
        if mol is None:
            print(f"⚠️ Skipping invalid molecule at index {idx}")
            continue
        
        # Canonical ID
        lig_id = f"LIG_{idx:04d}"
        
        # Original name from SDF
        orig_name = mol.GetProp("_Name") if mol.HasProp("_Name") else "NA"
        safe_name = sanitize_name(orig_name)
        
        # Directory structure
        lig_dirname = f"{lig_id}__{safe_name}"
        lig_root = os.path.join(ligands_root, lig_dirname)
        input_dir = os.path.join(lig_root, "input")
        os.makedirs(input_dir, exist_ok=True)
        
        # Write single-ligand SDF
        ligand_sdf = os.path.join(input_dir, "ligand.sdf")
        writer = Chem.SDWriter(ligand_sdf)
        writer.write(mol)
        writer.close()
        
        # Metadata
        smiles = Chem.MolToSmiles(mol)
        with open(os.path.join(lig_root, "META.txt"), "w") as f:
            f.write(f"ID={lig_id}\n")
            f.write(f"DIR_NAME={lig_dirname}\n")
            f.write(f"ORIGINAL_NAME={orig_name}\n")
            f.write(f"SMILES={smiles}\n")
            f.write(f"SDF_INDEX={idx}\n")
        
        ligands.append({
            "lig_id": lig_id,
            "lig_dirname": lig_dirname,
            "lig_root": lig_root,
            "ligand_sdf": ligand_sdf,
            "orig_name": orig_name,
            "smiles": smiles
        })
        
        mapping.append((lig_id, lig_dirname, orig_name, smiles))
    
    # Write mapping file
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
# RUN GNINA FOR ONE LIGAND
# =========================
def run_gnina(ligand_info: dict, idx: int, total: int) -> bool:
    """
    Run GNINA docking for a single ligand.
    
    Returns: True if successful, False otherwise
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
    
    # Mark as RUNNING
    write_status(lig_root, STATUS_RUNNING)
    
    print(f"\n🔄 [{idx}/{total}] Docking {lig_id} ...")
    start = time.time()
    
    cmd = [
        GNINA_BIN,
        "-r", f"{RESULTS_DIR}/protein/receptor.pdb",
        "-l", ligand_sdf,
        "--autobox_ligand", f"{RESULTS_DIR}/reference/ref_ligand.sdf",
        "--autobox_add", "5",
        "--autobox_extend", "1",
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
    
    # Save command for debugging
    with open(cmd_file, "w") as f:
        f.write(" \\\n    ".join(cmd))
    
    # Run GNINA - redirect stderr to file instead of PIPE
    try:
        with open(stderr_file, "w") as stderr_f:
            result = subprocess.run(
                cmd,
                check=True,
                stdout=subprocess.DEVNULL,  # Discard stdout (logs go to file)
                stderr=stderr_f,
                text=True
            )
        
        # Sanity check
        if not os.path.exists(out_lig) or os.path.getsize(out_lig) == 0:
            raise RuntimeError("GNINA produced empty output SDF")
        
        elapsed = (time.time() - start) / 60
        
        # Parse best score from output
        best_score = parse_best_score(out_lig)
        
        # Mark as DONE
        write_status(
            lig_root, 
            STATUS_DONE,
            elapsed_min=f"{elapsed:.2f}",
            best_cnn_score=f"{best_score:.4f}" if best_score else "NA"
        )
        
        print(f"✅ [{idx}/{total}] {lig_id} DONE in {elapsed:.2f} min (score: {best_score:.4f})")
        return True
        
    except subprocess.CalledProcessError as e:
        elapsed = (time.time() - start) / 60
        write_status(
            lig_root, 
            STATUS_FAILED,
            elapsed_min=f"{elapsed:.2f}",
            error=f"GNINA exit code {e.returncode}"
        )
        print(f"❌ [{idx}/{total}] {lig_id} FAILED: GNINA exit code {e.returncode}")
        return False
        
    except Exception as e:
        elapsed = (time.time() - start) / 60
        write_status(
            lig_root, 
            STATUS_FAILED,
            elapsed_min=f"{elapsed:.2f}",
            error=str(e),
            traceback=traceback.format_exc().replace("\n", " | ")
        )
        print(f"❌ [{idx}/{total}] {lig_id} FAILED: {e}")
        return False


def parse_best_score(sdf_path: str) -> float:
    """Extract best minimizedAffinity (or CNN score from docked SDF"""
    try:
        suppl = Chem.SDMolSupplier(sdf_path, removeHs=False)
        best_score = None
        for mol in suppl:
            if mol is None:
                continue
            if mol.HasProp("minimizedAffinity"):
                score = float(mol.GetProp("minimizedAffinity"))
                if best_score is None or score > best_score:
                    best_score = score
        return best_score if best_score else 0.0
    except:
        return 0.0


# =========================
# PROGRESS TRACKING
# =========================
def update_progress_csv(ligands: list, summary_dir: str):
    """Generate progress.csv with current status of all ligands"""
    progress_file = os.path.join(summary_dir, "progress.csv")
    
    with open(progress_file, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow([
            "ID", "DIR_NAME", "STATUS", "ELAPSED_MIN", 
            "BEST_CNN_SCORE", "START_TIME", "END_TIME"
        ])
        
        for lig in ligands:
            details = get_status_details(lig["lig_root"])
            writer.writerow([
                lig["lig_id"],
                lig["lig_dirname"],
                details.get("STATUS", STATUS_PENDING),
                details.get("ELAPSED_MIN", ""),
                details.get("BEST_CNN_SCORE", ""),
                details.get("START_TIME", ""),
                details.get("END_TIME", "")
            ])


def print_progress_summary(finished: list, failed: list, skipped: list, total: int):
    """Print current progress"""
    done = len(finished) + len(skipped)
    pct = (done / total) * 100 if total > 0 else 0
    print(f"\n📊 Progress: {done}/{total} ({pct:.1f}%) | ✅ {len(finished)} new | ⏭️ {len(skipped)} skipped | ❌ {len(failed)} failed")


# =========================
# MAIN PIPELINE
# =========================
def main():
    print("=" * 60)
    print("🧬 GNINA Flexible Docking Pipeline")
    print("=" * 60)
    
    # Setup
    prepare_root_folders()
    
    # Split ligands
    ligands = split_ligands(
        LIGAND_SDF,
        ligands_root=f"{RESULTS_DIR}/ligands"
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
        
        # Check resume status
        status = read_status(lig_root)
        
        if status == STATUS_DONE:
            print(f"⏭️ [{idx}/{total}] {lig_id} already DONE — skipping")
            skipped.append(lig_id)
            continue
        
        if status == STATUS_RUNNING:
            print(f"⚠️ [{idx}/{total}] {lig_id} was RUNNING (incomplete) — retrying")
        
        # Run docking
        success = run_gnina(lig, idx, total)
        
        if success:
            finished.append(lig_id)
        else:
            failed.append(lig_id)
        
        # Update progress every 10 ligands
        if idx % 10 == 0:
            update_progress_csv(ligands, f"{RESULTS_DIR}/summary")
            print_progress_summary(finished, failed, skipped, total)
    
    # Final summary
    elapsed_all = (time.time() - start_all) / 60
    
    # Write final files
    summary_dir = f"{RESULTS_DIR}/summary"
    
    with open(f"{summary_dir}/finished_ligands.txt", "w") as f:
        f.write("\n".join(finished + skipped))
    
    with open(f"{summary_dir}/failed_ligands.txt", "w") as f:
        f.write("\n".join(failed))
    
    update_progress_csv(ligands, summary_dir)
    
    # Print summary
    print("\n" + "=" * 60)
    print("📊 FINAL SUMMARY")
    print("=" * 60)
    print(f"⏱️  Total time: {elapsed_all:.2f} min")
    print(f"✅ Completed (new): {len(finished)}")
    print(f"⏭️  Skipped (cached): {len(skipped)}")
    print(f"❌ Failed: {len(failed)}")
    print(f"📁 Results: {RESULTS_DIR}")
    print("=" * 60)
    
    if failed:
        print(f"\n⚠️ Failed ligands: {', '.join(failed[:10])}" + 
              (f"... and {len(failed)-10} more" if len(failed) > 10 else ""))


# =========================
# ENTRY POINT
# =========================
if __name__ == "__main__":
    main()
