## Lựa chọn A: Làm mọi thứ trong WSL (qua VS Code, vẫn dùng GUI Windows)

```
Bước 1:  VS Code → F1 → "WSL: Connect to WSL"        ← Chuyển sang WSL TRƯỚC
Bước 2:  Tạo thư mục trong WSL terminal               ← mkdir ~/gnina_project
Bước 3:  Cài conda trong WSL                           ← Miniconda Linux
Bước 4:  Tạo conda env trong WSL                       ← conda env create -f ...
Bước 5:  Cài gnina trong WSL                            ← wget gnina Linux binary
Bước 6:  Copy data từ Windows vào WSL                  ← cp /mnt/d/... ~/gnina_project/
Bước 7:  Chạy notebook trong WSL                       ← Kernel = gnina_env (WSL)

         Tất cả trong 1 thế giới → Mọi thứ thấy nhau → Chạy được
```
> Bạn vẫn dùng VS Code trên Windows bình thường, chỉ là "connect" vào WSL

---

### Bước 1: Kiểm tra máy đã có gì

**🖱️ GUI trước:**
- Nhấn `Win` → gõ **"PowerShell"** → mở lên
- Nhấn `Win` → gõ **"Ubuntu"** (hoặc tên distro WSL của bạn) → mở terminal WSL lên

**⌨️ Rồi gõ lệnh kiểm tra:**

Trong **PowerShell**:
```powershell
wsl --list --verbose       # Xem WSL đang chạy gì
conda --version            # Kiểm tra conda trên WINDOWS
```

Trong **WSL terminal**:
```bash
conda --version
# "command not found" → cần cài (Bước 2)
# Có version rồi → bỏ qua Bước 2
```

---

### Bước 2: Cài Miniconda trong WSL (nếu chưa có)

**🖱️ GUI trước:**
- Mở trình duyệt → vào `https://docs.conda.io/en/latest/miniconda.html`
- Tìm file **Miniconda3 Linux 64-bit** → xem tên file (ví dụ `Miniconda3-latest-Linux-x86_64.sh`)
- *(Không cần tải tay — chỉ để biết tên file đúng)*

**⌨️ Rồi chạy trong WSL terminal:**
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

### Bước 3: Tạo project folder

**🖱️ GUI trước:**
- Mở **File Explorer** (Trên my computer) → điều hướng vào `\\wsl$\Ubuntu\home\<tên_user>\`
- Nhìn xem đã có folder `gnina_project` chưa → nếu chưa có thì dùng lệnh dưới

**⌨️ Tạo bằng lệnh trong WSL terminal:**
```bash
mkdir -p ~/gnina_project/data
mkdir -p ~/gnina_project/bin
cd ~/gnina_project
```

---

### Bước 4: Mở project trong VS Code

**🖱️ GUI — làm hết bằng chuột:**
```
1. Mở VS Code trên Windows bình thường
2. Nhấn F1 (hoặc Ctrl+Shift+P) → gõ "WSL: Connect to WSL" → Enter
3. Góc trái dưới VS Code chuyển thành: "WSL: Ubuntu"  ← xác nhận đã vào WSL
4. File → Open Folder → gõ vào ô đường dẫn: /home/<tên_user>/gnina_project → OK
```

> Từ giờ, **Terminal trong VS Code = WSL terminal**. Bạn gõ lệnh Linux bình thường.

---

### Bước 5: Tạo conda env

**🖱️ GUI trước:**
- Trong VS Code (đang WSL) → nhấn **Ctrl+`` ` ``** để mở Terminal tích hợp
- Terminal hiện ra ở dưới cùng → đây đã là WSL rồi, không cần mở thêm cửa sổ khác

**⌨️ Rồi chạy trong Terminal đó:**
```bash
# Tạo file YAML (copy-paste nguyên khối này)
cat > environment_gnina.yml << 'EOF'
name: gnina_env
channels:
  - conda-forge
  - pytorch
  - defaults

dependencies:
  - python=3.11
  - pytorch=2.2.0
  - pytorch-cuda=12.1
  - rdkit=2024.03.1
  - openbabel=3.1.1
  - pdbfixer
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
EOF

conda env create -f environment_gnina.yml    # ← 10-20 phút lần đầu, bình thường
conda activate gnina_env

python -m ipykernel install --user --name gnina_env --display-name "GNINA Env"

python -c "from rdkit import Chem; print('✅ RDKit OK')"
```

---

### Bước 6: Cài GNINA binary

**🖱️ GUI trước:**
- Mở trình duyệt → vào `https://github.com/gnina/gnina/releases`
- Tìm release **v1.3.2** → xem tên file binary Linux (dạng `gnina.1.3.2`, không có `.exe`)
- *(Không cần tải tay — chỉ để xác nhận link đúng)*

**⌨️ Rồi tải thẳng trong WSL terminal:** 
Dùng `https://github.com/gnina/gnina/releases/download/v1.3.2/gnina.1.3.2` hoặc `https://github.com/gnina/gnina/releases/download/v1.3.2/gnina.1.3.2.cuda12.8`
```bash
wget -q https://github.com/gnina/gnina/releases/download/v1.3.2/gnina.1.3.2 \
  -O ~/gnina_project/bin/gnina
chmod +x ~/gnina_project/bin/gnina
~/gnina_project/bin/gnina --version          # ← hiện version = thành công
```

---

### Bước 7: Copy data từ Windows vào WSL

**🖱️ GUI trước:**
- Mở **File Explorer** → vào thư mục chứa data trên Windows (ví dụ `D:\test_new_venv\data\`)
- Kiểm tra tên file chính xác: `ref_ligand.sdf`, `ligands_prepared.sdf`, `.pdb` — **đừng đánh tay, dễ sai**
- Mở thêm tab File Explorer thứ 2 → vào `\\wsl$\Ubuntu\home\<tên_user>\gnina_project\data\`
- Kéo thả file từ tab Windows sang tab WSL nếu muốn hoàn toàn không dùng lệnh

**⌨️ Hoặc dùng lệnh trong WSL terminal (nhanh hơn):**
```bash
# Trong WSL, ổ D: Windows = /mnt/d/
cp "/mnt/d/code python/open_protein_ligand_prep_pipeline(v2.0)/output/(READY) mmp2_7xjo_ready_for_gnina.pdb" \
   ~/gnina_project/data/mmp2_7xjo_ready_for_gnina.pdb

cp /mnt/d/test_new_venv/data/ref_ligand.sdf \
   ~/gnina_project/data/

cp /mnt/d/test_new_venv/data/ligands_prepared.sdf \
   ~/gnina_project/data/

ls -la ~/gnina_project/data/     # ← xác nhận 3 file đã vào
```

---

### Bước 8: Tạo file .env

**🖱️ GUI trước:**
- Trong VS Code (đang WSL) → **Explorer panel bên trái** (VS Code Explorer (Ctrl+Shift+E)) → chuột phải vào thư mục `gnina_project` → **New File** → đặt tên `.env`
- File `.env` trắng xuất hiện → paste nội dung dưới vào rồi **Ctrl+S** lưu

**⌨️ Hoặc tạo thẳng bằng lệnh:**
```bash
cat > ~/gnina_project/.env << EOF
GNINA_BIN=$HOME/gnina_project/bin/gnina
DOCKING_BASE_DIR=$HOME/gnina_project
PROTEIN_PATH=$HOME/gnina_project/data/mmp2_7xjo_ready_for_gnina.pdb
REF_LIGAND=$HOME/gnina_project/data/ref_ligand.sdf
LIGAND_SDF=$HOME/gnina_project/data/ligands_prepared.sdf
CUDA_VISIBLE_DEVICES=0
EOF

cat .env     # ← xác nhận nội dung đúng
```

---

### Bước 9: Tạo notebook và chạy

**🖱️ GUI — làm hết bằng chuột:**
```
1. VS Code (đang WSL) → File → New File → đặt tên: docking_pipeline.ipynb
2. Paste code notebook vào cell đầu tiên
3. Góc phải trên → "Select Kernel" → chọn "GNINA Env"
4. Shift+Enter → Chạy!
```

---

**Architecture**
```
~/gnina_project/                  ← MỞ FOLDER NÀY TRONG VS CODE (qua WSL: Connect)
│                                    Góc trái dưới VS Code hiện: "WSL: Ubuntu"
│
├── .env                          ← Config paths (dùng Linux paths: /home/...)
├── environment_gnina.yml         ← Dùng 1 lần tạo conda env
├── docking_pipeline.ipynb        ← ▶️ CHẠY FILE NÀY (Kernel: "GNINA Env")
│
├── bin/
│   └── gnina                     ← GNINA binary LINUX (không có .exe)
│
├── data/
│   ├── mmp2_7xjo_ready_for_gnina.pdb
│   ├── ref_ligand.sdf
│   └── ligands_prepared.sdf
│
└── docking_results/              ← Tự tạo khi chạy notebook
    └── .../
```
