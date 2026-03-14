

# 🪟 Setup GNINA trên Windows thuần — Có sẵn gnina.exe

## Giả định

```
Bạn đã có:
  ✅ gnina.exe (Windows build)
  ✅ Anaconda hoặc Miniconda trên Windows
  ✅ VS Code
  ✅ Data files (protein.pdb, ligands.sdf, ref_ligand.sdf)
```

---

## Thứ tự từng bước

### Bước 1: Mở VS Code Terminal (PowerShell)

```
VS Code → Terminal → New Terminal
(hoặc Ctrl + `)
```

Kiểm tra conda:
```powershell
conda --version
```

---

### Bước 2: Tạo project folder

```powershell
# Tạo thư mục project
mkdir D:\gnina_project
mkdir D:\gnina_project\data
mkdir D:\gnina_project\bin

cd D:\gnina_project
```

Copy `gnina.exe` vào `bin\`:
```powershell
# Copy gnina.exe vào project (sửa path nguồn cho đúng)
copy "D:\path\to\gnina.exe" "D:\gnina_project\bin\gnina.exe"

# Verify
D:\gnina_project\bin\gnina.exe --version
```

Copy data:
```powershell
copy "D:\code python\open_protein_ligand_prep_pipeline(v2.0)\output\(READY) mmp2_7xjo_ready_for_gnina.pdb" "D:\gnina_project\data\mmp2_7xjo_ready_for_gnina.pdb"

copy "D:\test_new_venv\data\ref_ligand.sdf" "D:\gnina_project\data\"
copy "D:\test_new_venv\data\ligands_prepared.sdf" "D:\gnina_project\data\"
```

---

### Bước 3: Tạo conda env MỚI

> ⚠️ **Tạo env mới = sandbox riêng, KHÔNG ảnh hưởng base hay env nào khác**

#### Cách A: Tạo bằng lệnh (nhanh, từng bước)

```powershell
# Tạo env mới tên "gnina_env" với Python 3.11
conda create -n gnina_env python=3.11 -y

# Activate
conda activate gnina_env

# Cài packages theo nhóm
conda install -c conda-forge rdkit=2024.03.1 -y
conda install -c conda-forge openbabel=3.1.1 -y
conda install -c conda-forge numpy scipy pandas matplotlib -y
conda install -c conda-forge jupyter ipykernel tqdm -y

# Cài pip packages
pip install python-dotenv py3Dmol psutil pyquaternion
```

#### Cách B: Tạo bằng YAML file (1 lệnh)

Tạo file `D:\gnina_project\environment.yml`:

```powershell
# Tạo file bằng PowerShell
@"
name: gnina_env
channels:
  - conda-forge
  - defaults

dependencies:
  - python=3.11
  - rdkit=2024.03.1
  - openbabel=3.1.1
  - numpy>=1.21.0
  - scipy>=1.7.0
  - pandas>=1.3.0
  - matplotlib>=3.4.0
  - jupyter
  - ipykernel
  - tqdm
  - pip
  - pip:
    - python-dotenv
    - py3Dmol>=2.0.0
    - psutil>=5.8.0
    - pyquaternion>=0.9.0
"@ | Out-File -FilePath "D:\gnina_project\environment.yml" -Encoding UTF8
```

Rồi tạo env:
```powershell
cd D:\gnina_project
conda env create -f environment.yml
```

---

### Bước 4: Đăng ký kernel cho Jupyter

```powershell
conda activate gnina_env

# Đăng ký kernel
python -m ipykernel install --user --name gnina_env --display-name "GNINA Env"

# Verify
python -c "from rdkit import Chem; print('RDKit OK')"
python -c "from dotenv import load_dotenv; print('dotenv OK')"
```

---

### Bước 5: Tạo file .env

Tạo file `D:\gnina_project\.env`:

```powershell
@"
GNINA_BIN=D:/gnina_project/bin/gnina.exe
DOCKING_BASE_DIR=D:/gnina_project
PROTEIN_PATH=D:/gnina_project/data/mmp2_7xjo_ready_for_gnina.pdb
REF_LIGAND=D:/gnina_project/data/ref_ligand.sdf
LIGAND_SDF=D:/gnina_project/data/ligands_prepared.sdf
CUDA_VISIBLE_DEVICES=0
"@ | Out-File -FilePath "D:\gnina_project\.env" -Encoding UTF8 -NoNewline
```

> ⚠️ Dùng **forward slash** `/` trong .env, Python đọc được trên cả Windows

---

### Bước 6: Mở project trong VS Code

```powershell
code D:\gnina_project
```

Hoặc:
```
VS Code → File → Open Folder → chọn D:\gnina_project
```

---

### Bước 7: Tạo notebook

```
File → New File → đặt tên: docking_pipeline.ipynb
Paste code .ipynb (1 cell) vào
```

---

### Bước 8: Chọn Kernel và Chạy

```
╔══════════════════════════════════════════════╗
║  Góc phải trên notebook → "Select Kernel"   ║
║  → Python Environments                       ║
║  → gnina_env (Python 3.11.x)                ║
║                                              ║
║  ⚠️ KHÔNG chọn base hay Python mặc định     ║
╚══════════════════════════════════════════════╝
```

Nhấn **Shift+Enter** hoặc **▶️ Run All**

---

## Cấu trúc cuối cùng

```
D:\gnina_project\                ← MỞ FOLDER NÀY TRONG VS CODE
│
├── .env                         ← Config paths
├── environment.yml              ← Dùng 1 lần tạo env
├── docking_pipeline.ipynb       ← ▶️ CHẠY FILE NÀY
│
├── bin\
│   └── gnina.exe                ← GNINA binary Windows
│
├── data\
│   ├── mmp2_7xjo_ready_for_gnina.pdb
│   ├── ref_ligand.sdf
│   └── ligands_prepared.sdf
│
└── docking_results\             ← Tự tạo khi chạy
    └── 8skl\
```

---

## Tóm tắt lệnh (copy-paste theo thứ tự)

```powershell
# ===== 1. Tạo folders =====
mkdir D:\gnina_project\data
mkdir D:\gnina_project\bin

# ===== 2. Copy files vào (sửa path cho đúng) =====
copy "path\to\gnina.exe"           "D:\gnina_project\bin\"
copy "path\to\protein.pdb"         "D:\gnina_project\data\"
copy "path\to\ref_ligand.sdf"      "D:\gnina_project\data\"
copy "path\to\ligands_prepared.sdf" "D:\gnina_project\data\"

# ===== 3. Tạo conda env =====
conda create -n gnina_env python=3.11 -y
conda activate gnina_env
conda install -c conda-forge rdkit openbabel numpy scipy pandas matplotlib jupyter ipykernel tqdm -y
pip install python-dotenv py3Dmol psutil pyquaternion

# ===== 4. Đăng ký kernel =====
python -m ipykernel install --user --name gnina_env --display-name "GNINA Env"

# ===== 5. Verify =====
python -c "from rdkit import Chem; print('OK')"
D:\gnina_project\bin\gnina.exe --version

# ===== 6. Mở VS Code =====
code D:\gnina_project
```

---

## Troubleshooting Windows

| Lỗi | Fix |
|------|-----|
| `conda không nhận` | Mở Anaconda Prompt thay PowerShell |
| `conda activate không work` | `conda init powershell` rồi restart terminal |
| `gnina.exe không chạy` | Thiếu CUDA / GPU driver → cài CUDA Toolkit |
| `rdkit import lỗi` | Sai kernel → chọn lại "GNINA Env" |
| `Permission denied` | Chạy VS Code as Administrator |
| `.env không load` | Kiểm tra encoding UTF-8 without BOM |

> **Lưu ý quan trọng:** `gnina.exe` cho Windows rất hiếm — bản chính thức chỉ có Linux binary. Nếu bạn không có gnina.exe thật, cần dùng WSL hoặc Colab.
