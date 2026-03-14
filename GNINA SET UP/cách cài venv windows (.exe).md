## Setup GNINA trên Windows thuần — Có sẵn gnina.exe

```
Bước 1:  Tạo folder bằng chuột (Windows Explorer)     ← Không cần mkdir
Bước 2:  Copy file bằng chuột (kéo thả)               ← Không cần lệnh copy
Bước 3:  Tạo conda env                                 ← Bắt buộc dùng lệnh
Bước 4:  Đăng ký kernel                                ← Bắt buộc dùng lệnh
Bước 5:  Tạo file .env bằng VS Code (chuột)           ← New File → paste → Ctrl+S
Bước 6:  Tạo file .yml bằng VS Code (chuột)           ← New File → paste → Ctrl+S
Bước 7:  Mở project trong VS Code (chuột)             ← File → Open Folder
Bước 8:  Tạo notebook và chạy (chuột)                 ← Select Kernel → Shift+Enter
```

### Giả định
```
Bạn đã có:
  ✅ gnina.exe (Windows build)
  ✅ Anaconda hoặc Miniconda trên Windows
  ✅ VS Code
  ✅ Data files (protein.pdb, ligands.sdf, ref_ligand.sdf)
```

---

### Bước 1: Tạo project folder bằng chuột

**🖱️ GUI — làm hết bằng chuột:**
```
1. Nhấn Win+E → mở Windows File Explorer
2. Vào ổ D: → chuột phải vào vùng trống
   → New → Folder → đặt tên: gnina_project → Enter

3. Vào trong gnina_project → tạo tiếp 2 thư mục con:
   → chuột phải → New → Folder → tên: data  → Enter
   → chuột phải → New → Folder → tên: bin   → Enter
```

---

### Bước 2: Copy file vào project bằng chuột

**🖱️ GUI — kéo thả hoàn toàn:**
```
Copy gnina.exe vào bin\:
1. Mở thư mục chứa gnina.exe của bạn
2. Kéo thả gnina.exe vào D:\gnina_project\bin\

Copy 3 file data vào data\:
1. Mở thư mục chứa data (ổ D:\code python\... hoặc D:\test_new_venv\data\)
2. Kéo thả 3 file vào D:\gnina_project\data\:
   • mmp2_7xjo_ready_for_gnina.pdb
   • ref_ligand.sdf
   • ligands_prepared.sdf
```

---

### Bước 3: Mở VS Code vào đúng project

**🖱️ GUI:**
```
1. Mở VS Code
2. File → Open Folder → tìm đến D:\gnina_project → OK
   → Cột Explorer bên trái hiện đúng 3 thứ: bin\, data\, và 2 file vừa kéo vào
```

---

### Bước 4: Tạo file environment.yml bằng chuột

**🖱️ GUI:**
```
1. Trong cột Explorer VS Code → click vào gnina_project (thư mục gốc)
2. Bấm icon "New File" → đặt tên: environment.yml → Enter
3. File mở ra (trắng) → paste nội dung bên dưới → Ctrl+S lưu
```

Nội dung paste vào:
```yaml
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
```

---

### Bước 5: Tạo conda env

> ⚠️ Bắt buộc dùng lệnh — conda cài packages tự động, không có GUI thay thế

**🖱️ GUI mở terminal:**
```
Trong VS Code:
→ Chuột phải vào thư mục gnina_project trong cột Explorer
→ Chọn "Open in Integrated Terminal"
→ Terminal mở ra, tự đứng sẵn ở D:\gnina_project\ luôn
```

**⌨️ Paste lần 1** — tạo env:
```powershell
conda env create -f environment.yml
```
> ⏳ 10–20 phút lần đầu, bình thường — đừng tắt

**⌨️ Paste lần 2** — sau khi xong lần 1:
```powershell
conda activate gnina_env
python -m ipykernel install --user --name gnina_env --display-name "GNINA Env"
python -c "from rdkit import Chem; print('✅ RDKit OK')"
```

> Nếu terminal báo `conda không nhận` → đóng terminal, mở lại bằng **Anaconda Prompt** thay vì PowerShell

---

### Bước 6: Tạo file .env bằng chuột

**🖱️ GUI:**
```
1. Trong cột Explorer VS Code → click vào gnina_project (thư mục gốc)
2. Bấm icon "New File" → đặt tên: .env → Enter
3. File mở ra → paste nội dung bên dưới → Ctrl+S lưu
```

Nội dung paste vào:
```env
GNINA_BIN=D:/gnina_project/bin/gnina.exe
DOCKING_BASE_DIR=D:/gnina_project
PROTEIN_PATH=D:/gnina_project/data/mmp2_7xjo_ready_for_gnina.pdb
REF_LIGAND=D:/gnina_project/data/ref_ligand.sdf
LIGAND_SDF=D:/gnina_project/data/ligands_prepared.sdf
CUDA_VISIBLE_DEVICES=0
```

> ⚠️ Dùng **forward slash `/`** trong .env (không phải `\`) — Python đọc được trên Windows

---

### Bước 7: Tạo notebook và chạy

**🖱️ GUI — làm hết bằng chuột:**
```
1. VS Code → File → New File → đặt tên: docking_pipeline.ipynb → Enter
2. Paste code notebook vào cell đầu tiên
3. Góc phải trên cùng → bấm "Select Kernel"
   → Python Environments
   → gnina_env (Python 3.11.x)
   ⚠️ KHÔNG chọn base hay Python mặc định

4. Shift+Enter → Chạy!
```

---

**Cấu trúc cuối cùng:**
```
D:\gnina_project\                ← MỞ FOLDER NÀY TRONG VS CODE
│
├── .env                         ← Tạo bằng chuột (New File VS Code)
├── environment.yml              ← Tạo bằng chuột (New File VS Code)
├── docking_pipeline.ipynb       ← ▶️ CHẠY FILE NÀY (Kernel: "GNINA Env")
│
├── bin\
│   └── gnina.exe                ← Kéo thả từ Windows Explorer vào
│
├── data\
│   ├── mmp2_7xjo_ready_for_gnina.pdb   ← Kéo thả vào
│   ├── ref_ligand.sdf                   ← Kéo thả vào
│   └── ligands_prepared.sdf             ← Kéo thả vào
│
└── docking_results\             ← Tự tạo khi chạy notebook
```

**Tổng kết: Lệnh bắt buộc chỉ còn 2 chỗ** — tạo conda env và đăng ký kernel. Tất cả còn lại đều dùng chuột.

---

### Troubleshooting Windows

| Lỗi | Fix |
|------|-----|
| `conda không nhận` | Mở **Anaconda Prompt** thay PowerShell |
| `conda activate không work` | Chạy `conda init powershell` → restart terminal |
| `gnina.exe không chạy` | Thiếu CUDA / GPU driver → cài CUDA Toolkit |
| `rdkit import lỗi` | Sai kernel → chọn lại **"GNINA Env"** |
| `Permission denied` | Chuột phải VS Code → **Run as Administrator** |
| `.env không load` | Kiểm tra encoding **UTF-8 without BOM** |

> **Lưu ý quan trọng:** `gnina.exe` cho Windows rất hiếm — bản chính thức chỉ có Linux binary. Nếu không có gnina.exe thật, cần dùng WSL hoặc Colab.
