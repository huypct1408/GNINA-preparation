Đúng rồi! Để tôi viết lại bản hoàn chỉnh, **chuột tối đa — lệnh chỉ khi bắt buộc**:

---

## Lựa chọn A: Làm mọi thứ trong WSL (qua VS Code, vẫn dùng GUI Windows)

```
Bước 1:  VS Code → F1 → "WSL: Connect to WSL"        ← Chuyển sang WSL TRƯỚC
Bước 2:  Tạo thư mục trong VS Code (chuột)            ← New Folder trong Explorer
Bước 3:  Cài conda trong WSL                           ← Bắt buộc dùng lệnh
Bước 4:  Tạo file .yml bằng VS Code (chuột)           ← New File → paste → Ctrl+S
Bước 5:  Chạy conda env                                ← Bắt buộc dùng lệnh
Bước 6:  Cài gnina                                     ← Bắt buộc dùng lệnh
Bước 7:  Copy data bằng kéo thả chuột                 ← Không cần lệnh cp
Bước 8:  Tạo file .env bằng VS Code (chuột)           ← New File → paste → Ctrl+S
Bước 9:  Chạy notebook                                 ← Chọn Kernel → Shift+Enter
```

> Bạn vẫn dùng VS Code trên Windows bình thường, chỉ là "connect" vào WSL

---

### Bước 1: Kiểm tra máy đã có gì

**🖱️ GUI:**
- Nhấn `Win` → gõ **PowerShell** → mở lên
- Nhấn `Win` → gõ **Ubuntu** → mở terminal WSL lên

**⌨️ Lệnh kiểm tra — paste vào, không cần gõ tay:**

Trong **PowerShell**:
```powershell
wsl --list --verbose
conda --version
```

Trong **WSL terminal**:
```bash
conda --version
# "command not found" → cần cài (Bước 2)
# Có version rồi → bỏ qua Bước 2
```

---

### Bước 2: Cài Miniconda trong WSL (nếu chưa có)

> ⚠️ Bước này bắt buộc dùng lệnh — hệ thống cài ngầm bên dưới, không có GUI thay thế được

**⌨️ Paste từng dòng vào WSL terminal:**
```bash
mkdir -p ~/miniconda3
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh \
  -O ~/miniconda3/miniconda.sh
bash ~/miniconda3/miniconda.sh -b -u -p ~/miniconda3
rm ~/miniconda3/miniconda.sh
~/miniconda3/bin/conda init bash
source ~/.bashrc
conda --version
```

> **Conda trong WSL hoàn toàn tách biệt với conda trên Windows.** Cài thêm không ảnh hưởng gì.

---

### Bước 3: Mở VS Code và tạo project folder bằng chuột

**🖱️ GUI — làm hết bằng chuột:**
```
1. Mở VS Code
2. Nhấn F1 → gõ "WSL: Connect to WSL" → Enter
   → Góc trái dưới VS Code hiện: "WSL: Ubuntu"  ← xác nhận đã vào WSL

3. File → Open Folder → gõ: /home/<tên_user_wsl> → OK
   → Bây giờ bạn đang nhìn thấy thư mục home của WSL

4. Nhìn cột Explorer bên trái VS Code:
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
3. File mở ra (trắng trơn) → paste nội dung bên dưới vào → Ctrl+S lưu
```

Nội dung paste vào file:
```yaml
# GNINA v1.3.2 Production Environment (CORRECTED & VALIDATED)
# Fixes applied based on repository evidence verification

name: gnina_env
channels:
  - conda-forge          # Single authoritative channel for consistency
  - defaults             # Fallback only

dependencies:
  # ==================== CORE PYTHON ====================
  - python=3.11
  
  # ==================== DEEP LEARNING ====================
  - pytorch::pytorch=2.8.0
  - pytorch::torchvision           # Fixed: Matches PyTorch 2.8.0, torchvision=0.23.0 
  - pytorch::torchaudio=2.8.0      # Fixed: Matches PyTorch 2.8.0
  - pytorch::pytorch-cuda=12.1     # Required for GPU support
  
  # ==================== CHEMINFORMATICS ====================
  # FIX #3: Pin RDKit more specifically for reproducibility
  - rdkit=2024.03.1      # Specific patch version (tested, stable)
  - openbabel=3.1.1
  
  # ==================== STRUCTURE REPAIR ====================
  - pdbfixer
  - mdanalysis
  
  # ==================== CRITICAL RUNTIME DEPENDENCIES ====================
  - boost>=1.81.0
  - eigen>=3.4.0
  - glog                 # FIX #1: CORRECTED from libgoogle-glog
  - protobuf>=3.20
  
  # ==================== NUMERICAL COMPUTING ====================
  - numpy>=1.21.0
  - scipy>=1.7.0
  - pandas>=1.3.0
  - scikit-learn>=1.0.0
  
  # ==================== DEVELOPMENT & ANALYSIS ====================
  - matplotlib>=3.4.0
  - jupyter
  - ipython
  - tqdm
  - git

pip:
  # Alternative: Install PyTorch via pip for guaranteed ABI compatibility
  # - torch==2.8.0
  # - torchvision==0.17.0
  # - torchaudio==2.8.0
  
  - git+https://github.com/forlilab/molscrub.git
  - py3Dmol>=2.0.0
  - pyquaternion>=0.9.0
  - psutil>=5.8.0
```

---

### Bước 5: Tạo conda env
## Quy trình đúng:

**🖱️ Cách nhanh nhất — không cần gõ đường dẫn:**
```
Trong VS Code, cột Explorer bên trái
→ Chuột phải vào thư mục gnina_project
→ Chọn "Open in Integrated Terminal"
→ Terminal mở ra, tự động đứng sẵn ở ~/gnina_project/ luôn
```

**⌨️ Hoặc cd thủ công trong terminal:**
```bash
cd ~/gnina_project         # ← vào đúng thư mục trước
```

**Sau đó paste bình thường:**
> ⚠️ Bước này bắt buộc dùng lệnh — conda cài packages tự động, không có GUI thay thế

**🖱️ GUI mở terminal:**
```
Trong VS Code → Terminal → New Terminal
(Terminal hiện ra ở dưới cùng — đây đã là WSL rồi)
```

**⌨️ Paste lần 1** — chạy tạo env:
```bash
conda env create -f environment_gnina.yml
```
> ⏳ 10–20 phút lần đầu, bình thường — đừng tắt

**⌨️ Paste lần 2** — sau khi xong lần 1:
```bash
conda activate gnina_env
python -m ipykernel install --user --name gnina_env --display-name "GNINA Env"
python -c "from rdkit import Chem; print('✅ RDKit OK')"
```

---

### Bước 6: Cài GNINA binary

> ⚠️ Bắt buộc dùng lệnh — file Linux binary cần cấp quyền đặc biệt, tải bằng trình duyệt Windows hay bị lỗi phân quyền

**⌨️ Paste vào terminal VS Code:**
```bash
wget -q https://github.com/gnina/gnina/releases/download/v1.3.2/gnina.1.3.2.cuda12.8 \
  -O ~/gnina_project/bin/gnina
chmod +x ~/gnina_project/bin/gnina
~/gnina_project/bin/gnina --version     # ← hiện version = thành công
```

---

### Bước 7: Copy data từ Windows vào WSL — kéo thả chuột

**🖱️ GUI — kéo thả hoàn toàn:**
```
1. Nhấn Win+E → mở Windows File Explorer
2. Điều hướng đến thư mục chứa 3 file data của bạn trên ổ D:
3. Quay lại VS Code → trong cột Explorer bên trái → click vào thư mục data
4. Cầm chuột kéo 3 file từ Windows File Explorer
   THẢ vào thẳng thư mục data trong cột Explorer của VS Code

   3 file cần kéo vào:
   • mmp2_7xjo_ready_for_gnina.pdb
   • ref_ligand.sdf
   • ligands_prepared.sdf

5. Thấy 3 file hiện trong data/ của VS Code → xong!
```

---

### Bước 8: Tạo file .env bằng chuột

**🖱️ GUI — tạo file không cần lệnh:**
```
1. Trong cột Explorer VS Code, click vào thư mục gnina_project (thư mục gốc)
2. Bấm icon "New File" → đặt tên: .env → Enter
3. File mở ra → paste nội dung bên dưới vào
4. Thay <tên_user> bằng tên user WSL thực tế của bạn
   (không biết tên user → gõ lệnh: whoami trong terminal)
5. Ctrl+S lưu lại
```

Nội dung paste vào file:
```env
GNINA_BIN=/home/<tên_user>/gnina_project/bin/gnina
DOCKING_BASE_DIR=/home/<tên_user>/gnina_project
PROTEIN_PATH=/home/<tên_user>/gnina_project/data/mmp2_7xjo_ready_for_gnina.pdb
REF_LIGAND=/home/<tên_user>/gnina_project/data/ref_ligand.sdf
LIGAND_SDF=/home/<tên_user>/gnina_project/data/ligands_prepared.sdf
CUDA_VISIBLE_DEVICES=0
```

---

### Bước 9: Tạo notebook và chạy

**🖱️ GUI — làm hết bằng chuột:**
```
1. VS Code → File → New File → đặt tên: docking_pipeline.ipynb → Enter
2. Paste code notebook vào cell đầu tiên
3. Góc phải trên cùng → bấm "Select Kernel" → chọn "GNINA Env"
4. Shift+Enter → Chạy!
```

---

**Architecture**
```
~/gnina_project/                  ← MỞ FOLDER NÀY TRONG VS CODE (qua WSL: Connect)
│                                    Góc trái dưới VS Code hiện: "WSL: Ubuntu"
│
├── .env                          ← Tạo bằng chuột (New File trong VS Code)
├── environment_gnina.yml         ← Tạo bằng chuột (New File trong VS Code)
├── docking_pipeline.ipynb        ← ▶️ CHẠY FILE NÀY (Kernel: "GNINA Env")
│
├── bin/
│   └── gnina                     ← Cài bằng lệnh wget (bắt buộc)
│
├── data/
│   ├── mmp2_7xjo_ready_for_gnina.pdb   ← Kéo thả từ Windows vào
│   ├── ref_ligand.sdf                   ← Kéo thả từ Windows vào
│   └── ligands_prepared.sdf             ← Kéo thả từ Windows vào
│
└── docking_results/              ← Tự tạo khi chạy notebook
```

**Tổng kết: Lệnh bắt buộc chỉ còn 3 chỗ** — cài conda, chạy conda env, cài gnina binary. Tất cả còn lại đều dùng chuột trong VS Code.
