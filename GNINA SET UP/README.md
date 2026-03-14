

# Quy trình GNINA Docking — Phiên bản Linux + Fish Shell (Đã sửa)

---

## Lựa chọn A: Làm mọi thứ trên Linux native (VS Code + Fish Shell)

```
Bước 1:  Kiểm tra máy đã có gì                        ← Lệnh kiểm tra
Bước 2:  Cài Miniconda (nếu chưa có)                  ← Bắt buộc dùng lệnh
Bước 3:  Mở VS Code và tạo project folder (chuột)     ← New Folder trong Explorer
Bước 4:  Tạo file .yml bằng VS Code (chuột)           ← New File → paste → Ctrl+S
Bước 5:  Chạy conda env                                ← Bắt buộc dùng lệnh
Bước 6:  Cài gnina binary                              ← Bắt buộc dùng lệnh
Bước 7:  ⭐ Tạo conda hook cho Fish Shell              ← BƯỚC QUAN TRỌNG NHẤT
Bước 8:  Copy data vào thư mục project (chuột)        ← Kéo thả trong VS Code
Bước 9:  Tạo file .env bằng VS Code (chuột)           ← New File → paste → Ctrl+S
Bước 10: Tạo notebook và chạy                          ← Chọn Kernel → Shift+Enter
```

---

### Bước 1: Kiểm tra máy đã có gì

**⌨️ Mở terminal Fish, paste từng dòng:**
```fish
conda --version
# "command not found" → cần cài (Bước 2)
# Có version rồi → bỏ qua Bước 2

nvidia-smi
# Hiện bảng GPU → có GPU, tốt
# Lỗi → không có GPU, gnina vẫn chạy được với --no_gpu
```

---

### Bước 2: Cài Miniconda (nếu chưa có)

> ⚠️ Bắt buộc dùng lệnh

**⌨️ Paste từng dòng vào Fish terminal:**
```fish
mkdir -p ~/miniconda3

wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh \
  -O ~/miniconda3/miniconda.sh

bash ~/miniconda3/miniconda.sh -b -u -p ~/miniconda3

rm ~/miniconda3/miniconda.sh

~/miniconda3/bin/conda init fish

source ~/.config/fish/config.fish

conda --version
```

---

### Bước 3: Mở VS Code và tạo project folder bằng chuột

**🖱️ GUI — làm hết bằng chuột:**
```
1. Mở VS Code
2. File → Open Folder → điều hướng đến thư mục làm việc của bạn → OK
   Ví dụ: /home/labhhc5/Documents/workspace/D21/Duong Huy/

3. Nhìn cột Explorer bên trái VS Code:
   → Bấm icon "New Folder" → đặt tên: gnina_project → Enter
   → Click vào gnina_project để mở vào trong
   → Bấm "New Folder" thêm 2 lần:
      • Tên thư mục 1: data
      • Tên thư mục 2: bin
```

---

### Bước 4: Tạo file môi trường conda bằng chuột

**🖱️ GUI — tạo file không cần lệnh:**
```
1. Trong cột Explorer VS Code, click vào thư mục gnina_project
2. Bấm icon "New File" → đặt tên: environment_gnina.yml → Enter
3. File mở ra → paste nội dung bên dưới vào → Ctrl+S lưu
```

**Nội dung paste vào file:**
```yaml
# GNINA v1.3.2 Production Environment
# Tested on Linux + Fish Shell + CUDA 12.x

name: gnina_env
channels:
  - pytorch
  - nvidia
  - conda-forge
  - defaults

dependencies:
  # ==================== CORE PYTHON ====================
  - python=3.11
  
  # ==================== DEEP LEARNING ====================
  - pytorch=2.5.1
  - torchvision
  - torchaudio
  - pytorch-cuda=12.1
  
  # ==================== CUDA RUNTIME (cho gnina) ====================
  - cudnn=9
  - cuda-nvtx
  - cudatoolkit
  
  # ==================== CHEMINFORMATICS ====================
  - rdkit=2024.03.1
  - openbabel=3.1.1
  
  # ==================== STRUCTURE REPAIR ====================
  - pdbfixer
  - mdanalysis
  
  # ==================== CRITICAL RUNTIME DEPENDENCIES ====================
  - boost>=1.81.0
  - eigen>=3.4.0
  - glog
  - protobuf>=3.20
  
  # ==================== NUMERICAL COMPUTING ====================
  - numpy>=1.21.0
  - scipy>=1.7.0
  - pandas>=1.3.0
  - scikit-learn>=1.0.0
  
  # ==================== DEVELOPMENT & ANALYSIS ====================
  - matplotlib>=3.4.0
  - jupyter
  - ipykernel
  - ipython
  - tqdm
  - git
  - python-dotenv

  - pip:
    - git+https://github.com/forlilab/molscrub.git
    - py3Dmol>=2.0.0
    - pyquaternion>=0.9.0
    - psutil>=5.8.0
```

> ⚠️ So với bản cũ: **thêm `cudnn=9`, `cuda-nvtx`, `cudatoolkit`, `ipykernel`, `python-dotenv`** — đây là các thư viện gnina cần mà bản cũ thiếu.

---

### Bước 5: Tạo conda env

**🖱️ Mở terminal trong VS Code:**
```
Trong VS Code → chuột phải vào thư mục gnina_project
→ Chọn "Open in Integrated Terminal"
```

**⌨️ Paste lần 1** — tạo env:
```fish
conda env create -f environment_gnina.yml
```
> ⏳ 10–20 phút lần đầu — đừng tắt

**⌨️ Paste lần 2** — sau khi xong:
```fish
conda activate gnina_env

python -m ipykernel install --user --name gnina_env --display-name "GNINA Env"

python -c "from rdkit import Chem; print('✅ RDKit OK')"
```

---

### Bước 6: Cài GNINA binary

> ⚠️ Bắt buộc dùng lệnh

**⌨️ Paste vào terminal (đảm bảo đang ở trong gnina_project):**
```fish
wget -q https://github.com/gnina/gnina/releases/download/v1.3.2/gnina.1.3.2.cuda12.8 \
  -O bin/gnina

chmod +x bin/gnina
```

> ❌ CHƯA TEST `gnina --help` được ở bước này — cần Bước 7 trước!

---

### Bước 7: ⭐ Tạo Conda Hook cho Fish Shell (BƯỚC QUAN TRỌNG NHẤT)

> Đây là bước bản cũ thiếu. Không có bước này, gnina sẽ báo lỗi
> `libcudnn.so.9: cannot open shared object file`
> mỗi khi mở terminal mới.

**⌨️ Paste từng dòng:**

```fish
mkdir -p "$CONDA_PREFIX/etc/conda/activate.d"
```

```fish
mkdir -p "$CONDA_PREFIX/etc/conda/deactivate.d"
```

**Tạo activation hook (`.fish` — KHÔNG phải `.sh`):**
```fish
printf 'set -gx LD_LIBRARY_PATH "$CONDA_PREFIX/lib"\n' > "$CONDA_PREFIX/etc/conda/activate.d/gnina_ld.fish"
```

**Tạo deactivation hook:**
```fish
printf 'set -e LD_LIBRARY_PATH\n' > "$CONDA_PREFIX/etc/conda/deactivate.d/gnina_ld.fish"
```

**Kích hoạt lại để hook có hiệu lực:**
```fish
conda deactivate
conda activate gnina_env
```

**Xác nhận hoạt động:**
```fish
ldd bin/gnina | grep "not found"
```
> Nếu **không hiện gì** = ✅ Thành công

```fish
bin/gnina --help
```
> Hiện danh sách options = ✅ gnina chạy được

---

### ⚠️ Giải thích tại sao cần Bước 7

| Vấn đề | Nguyên nhân | Giải pháp |
|---|---|---|
| gnina cần `libcudnn.so.9`, `libnvToolsExt.so.1` | Các file này nằm trong `$CONDA_PREFIX/lib` nhưng hệ thống không tự biết | Conda hook tự set `LD_LIBRARY_PATH` mỗi khi activate |
| Fish shell dùng dấu **cách** phân cách path, Linux cần dấu **hai chấm** `:` | Hook `.sh` không hoạt động đúng trong Fish | Phải dùng hook `.fish` với `set -gx` |
| Mỗi terminal mới, `LD_LIBRARY_PATH` bị reset | Biến môi trường không lưu giữa các session | Hook tự chạy lại mỗi khi `conda activate` |

---

### Bước 8: Copy data vào thư mục project — kéo thả chuột

**🖱️ GUI — kéo thả hoàn toàn:**
```
1. Mở File Manager (Nautilus/Dolphin)
2. Điều hướng đến thư mục chứa 3 file data của bạn
3. Quay lại VS Code → trong cột Explorer bên trái → click vào thư mục data
4. Kéo 3 file từ File Manager THẢ vào thư mục data trong VS Code:

   • mmp2_7xjo_ready_for_gnina.pdb
   • ref_ligand.sdf
   • ligands_prepared.sdf

5. Thấy 3 file hiện trong data/ → xong!
```

---

### Bước 9: Tạo file .env bằng chuột

**🖱️ GUI:**
```
1. Trong cột Explorer VS Code, click vào thư mục gnina_project (gốc)
2. Bấm icon "New File" → đặt tên: .env → Enter
3. Paste nội dung bên dưới vào → Ctrl+S lưu
```

**Nội dung paste vào file** (thay đường dẫn cho đúng máy bạn):
```env
GNINA_BIN=/home/labhhc5/Documents/workspace/D21/Duong Huy/gnina_project/bin/gnina
DOCKING_BASE_DIR=/home/labhhc5/Documents/workspace/D21/Duong Huy/gnina_project
PROTEIN_PATH=/home/labhhc5/Documents/workspace/D21/Duong Huy/gnina_project/data/mmp2_7xjo_ready_for_gnina.pdb
REF_LIGAND=/home/labhhc5/Documents/workspace/D21/Duong Huy/gnina_project/data/ref_ligand.sdf
LIGAND_SDF=/home/labhhc5/Documents/workspace/D21/Duong Huy/gnina_project/data/ligands_prepared.sdf
CUDA_VISIBLE_DEVICES=0
```

> Không biết đường dẫn chính xác → gõ `pwd` trong terminal khi đang ở gnina_project

---

### Bước 10: Tạo notebook và chạy

**🖱️ GUI:**
```
1. VS Code → chuột phải vào gnina_project → New File
   → đặt tên: docking_pipeline.ipynb → Enter
2. Paste code notebook vào cell đầu tiên (code bên dưới)
3. Góc phải trên cùng → bấm "Select Kernel" → chọn "GNINA Env"
4. Shift+Enter → Chạy!
```

---

## Code Notebook (đã cập nhật cho Linux + Fish)

```python
# ============================================================
# 🧬 GNINA Flexible Docking Pipeline v2.1 — Linux + Fish Shell
# ============================================================
# Thay đổi so với v2.0:
#   1. Thêm đảm bảo LD_LIBRARY_PATH trong subprocess
#   2. Path lấy từ .env (giữ nguyên logic gốc)
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
# [1/3] .env loading + LD_LIBRARY_PATH đảm bảo
# =============================================================
NOTEBOOK_DIR = Path.cwd()

# Đảm bảo LD_LIBRARY_PATH luôn đúng (phòng trường hợp Jupyter
# kernel không kế thừa đầy đủ biến từ conda hook)
_conda_prefix = os.environ.get("CONDA_PREFIX", "")
if _conda_prefix:
    _conda_lib = os.path.join(_conda_prefix, "lib")
    _current_ld = os.environ.get("LD_LIBRARY_PATH", "")
    if _conda_lib not in _current_ld:
        os.environ["LD_LIBRARY_PATH"] = f"{_conda_lib}:{_current_ld}" if _current_ld else _conda_lib
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
# GLOBAL CONFIG
# =========================
BASE_DIR = _resolve_path("DOCKING_BASE_DIR", str(NOTEBOOK_DIR))
RESULTS_DIR = f"{BASE_DIR}/docking_results"
PROTEIN_PATH = _resolve_path("PROTEIN_PATH", f"{BASE_DIR}/data/mmp2_7xjo_ready_for_gnina.pdb")
REF_LIGAND = _resolve_path("REF_LIGAND", f"{BASE_DIR}/data/ref_ligand.sdf")
LIGAND_SDF = _resolve_path("LIGAND_SDF", f"{BASE_DIR}/data/ligands_prepared.sdf")
FLEX_RESIDUES = "A:182,A:181,A:215,A:262,A:49"
SEED = "42"
GPU_DEVICE = os.environ.get("CUDA_VISIBLE_DEVICES", "0")

# GNINA_BIN — resolve
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
# [2/3] Subprocess env — đảm bảo gnina tìm được thư viện
# =========================
def _get_subprocess_env() -> dict:
    """Tạo env dict cho subprocess, đảm bảo LD_LIBRARY_PATH có conda lib."""
    env = os.environ.copy()
    conda_prefix = env.get("CONDA_PREFIX", "")
    if conda_prefix:
        conda_lib = os.path.join(conda_prefix, "lib")
        current_ld = env.get("LD_LIBRARY_PATH", "")
        if conda_lib not in current_ld:
            env["LD_LIBRARY_PATH"] = f"{conda_lib}:{current_ld}" if current_ld else conda_lib
    return env


# =========================
# STATUS MANAGEMENT — 100% giữ nguyên logic gốc
# =========================
def write_status(lig_root: str, status: str, **kwargs):
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
# UTILS — 100% giữ nguyên
# =========================
def sanitize_name(name: str, max_len: int = 80) -> str:
    if not name:
        return "NA"
    name = name.strip().replace(" ", "_")
    name = re.sub(r'[\/:*?"<>|\\]', '', name)
    name = re.sub(r'_+', '_', name)
    name = name.strip("_")
    return name[:max_len] if name else "NA"


# =========================
# PREPARE ROOT FOLDERS — 100% giữ nguyên
# =========================
def prepare_root_folders():
    dirs = [
        f"{RESULTS_DIR}/protein",
        f"{RESULTS_DIR}/reference",
        f"{RESULTS_DIR}/ligands",
        f"{RESULTS_DIR}/summary"
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
# SPLIT LIGANDS — 100% giữ nguyên
# =========================
def split_ligands(input_sdf: str, ligands_root: str) -> list:
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

        ligands.append({
            "lig_id": lig_id,
            "lig_dirname": lig_dirname,
            "lig_root": lig_root,
            "ligand_sdf": ligand_sdf,
            "orig_name": orig_name,
            "smiles": smiles
        })
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
# RUN GNINA — thêm env=_get_subprocess_env() vào subprocess
# =========================
def run_gnina(ligand_info: dict, idx: int, total: int) -> bool:
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

    with open(cmd_file, "w") as f:
        f.write(" \\\n    ".join(cmd))

    try:
        with open(stderr_file, "w") as stderr_f:
            result = subprocess.run(
                cmd,
                check=True,
                stdout=subprocess.DEVNULL,
                stderr=stderr_f,
                text=True,
                env=_get_subprocess_env()  # ← ĐẢM BẢO gnina tìm được thư viện
            )

        if not os.path.exists(out_lig) or os.path.getsize(out_lig) == 0:
            raise RuntimeError("GNINA produced empty output SDF")

        elapsed = (time.time() - start) / 60
        best_score = parse_best_score(out_lig)

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
# PROGRESS TRACKING — 100% giữ nguyên
# =========================
def update_progress_csv(ligands: list, summary_dir: str):
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
    done = len(finished) + len(skipped)
    pct = (done / total) * 100 if total > 0 else 0
    print(f"\n📊 Progress: {done}/{total} ({pct:.1f}%) | ✅ {len(finished)} new | ⏭️ {len(skipped)} skipped | ❌ {len(failed)} failed")


# =========================
# MAIN — 100% giữ nguyên logic
# =========================
def main():
    print("=" * 60)
    print("🧬 GNINA Flexible Docking Pipeline v2.1")
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

        status = read_status(lig_root)

        if status == STATUS_DONE:
            print(f"⏭️ [{idx}/{total}] {lig_id} already DONE — skipping")
            skipped.append(lig_id)
            continue

        if status == STATUS_RUNNING:
            print(f"⚠️ [{idx}/{total}] {lig_id} was RUNNING (incomplete) — retrying")

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


# =============================================================
# [3/3] Chạy trực tiếp
# =============================================================
main()
```

---

## Tổng kết thay đổi so với bản cũ

| Mục | Bản cũ | Bản mới |
|---|---|---|
| Hệ thống | WSL (Windows) | Linux native + Fish shell |
| Shell | Bash | Fish |
| Conda hook | Không có | ⭐ **`.fish` hook** tự set `LD_LIBRARY_PATH` |
| environment.yml | Thiếu `cudnn`, `cuda-nvtx` | Thêm `cudnn=9`, `cuda-nvtx`, `cudatoolkit` |
| Notebook code | `subprocess.run()` không truyền env | Thêm `env=_get_subprocess_env()` + auto-set `LD_LIBRARY_PATH` |
| Bước quan trọng nhất | Không có | **Bước 7: Tạo conda hook** |

---

## Architecture

```
~/Documents/workspace/D21/Duong Huy/gnina_project/
│
├── .env                          ← Tạo bằng chuột (Bước 9)
├── environment_gnina.yml         ← Tạo bằng chuột (Bước 4)
├── docking_pipeline.ipynb        ← ▶️ CHẠY FILE NÀY (Bước 10)
│
├── bin/
│   └── gnina                     ← Cài bằng wget (Bước 6)
│
├── data/
│   ├── mmp2_7xjo_ready_for_gnina.pdb   ← Kéo thả (Bước 8)
│   ├── ref_ligand.sdf                   ← Kéo thả (Bước 8)
│   └── ligands_prepared.sdf             ← Kéo thả (Bước 8)
│
└── docking_results/              ← Tự tạo khi chạy notebook
```

## Quy trình hàng ngày (sau khi setup xong):

```fish
conda activate gnina_env
# Mở VS Code → mở notebook → Shift+Enter
```

**Chỉ 1 lệnh duy nhất trước khi làm việc.** Mọi thứ còn lại tự động.
