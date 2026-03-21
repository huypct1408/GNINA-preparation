quy trinh chuan bi 1HFC 1.56 ANGSTROM STRUCTURE OF MATURE TRUNCATED HUMAN FIBROBLAST COLLAGENASE voi the crystal ligand PLH	(METHYLAMINO-PHENYLALANYL-LEUCYL-HYDROXAMIC ACID)

khi danh gia re-docking 
```
# Check rmsd
!"/home/labhhc5/miniconda3/envs/gnina_env/bin/obrms" \
-f "/home/labhhc5/Documents/workspace/D21/Duong Huy/gnina_project/data/mmp1_1hfc_crystal_ligand.sdf" \
"/home/labhhc5/Documents/workspace/D21/Duong Huy/gnina_project/redocking_results/mmp1_1hfc_redocked.sdf"
```
Thi thu duoc ket qua output:
RMSD mmp1_1hfc_crystal_ligand:mmp1_1hfc_crystal_ligand 2.20929
RMSD mmp1_1hfc_crystal_ligand:mmp1_1hfc_crystal_ligand 2.70654
RMSD mmp1_1hfc_crystal_ligand:mmp1_1hfc_crystal_ligand 2.37947
RMSD mmp1_1hfc_crystal_ligand:mmp1_1hfc_crystal_ligand 6.13337
RMSD mmp1_1hfc_crystal_ligand:mmp1_1hfc_crystal_ligand 7.08122
RMSD mmp1_1hfc_crystal_ligand:mmp1_1hfc_crystal_ligand 6.62858
RMSD mmp1_1hfc_crystal_ligand:mmp1_1hfc_crystal_ligand 6.95338
RMSD mmp1_1hfc_crystal_ligand:mmp1_1hfc_crystal_ligand 3.26056
RMSD mmp1_1hfc_crystal_ligand:mmp1_1hfc_crystal_ligand 6.89451

Chào bạn, 

Dựa trên tài liệu chính thức từ khóa học **IQB-2025 (Python for Cheminformatics-Driven Molecular Docking)**, việc bạn muốn sử dụng trực tiếp tọa độ của ligand từ file tinh thể để chạy re-docking (cognate docking) nhằm đánh giá độ chính xác của thuật toán GNINA là hoàn toàn đúng đắn. 

Tuy nhiên, tài liệu IQB-2025 nhấn mạnh rằng bạn **không thể chỉ cắt trực tiếp ligand từ file PDB và dùng ngay**. Lý do là vì tọa độ ligand trong file PDB tinh thể thường **thiếu nguyên tử Hydro** và **không có thông tin về bậc liên kết (bond orders)**. GNINA cần các thông tin hóa học chuẩn xác này để tính toán năng lượng tương tác.

Để giải quyết vấn đề này, quy trình chuẩn của workshop IQB-2025 (trong Notebook 3) sử dụng **RDKit** để ánh xạ (map) bậc liên kết từ một cấu trúc "chuẩn" (ideal) sang tọa độ 3D của tinh thể, sau đó thêm Hydro. 



Dưới đây là cách thực hiện chi tiết bằng Python theo đúng chuẩn của khóa học:

### Bước 1: Chuẩn bị 2 file đầu vào
1. **File tọa độ gốc (`ligand_fromPDB.pdb`):** Bạn dùng PyMOL hoặc MDAnalysis để trích xuất riêng ligand từ phức hợp protein-ligand (ví dụ 1HFC) và lưu thành file `.pdb`. File này chứa tọa độ 3D thực tế nhưng sai/thiếu liên kết hóa học.
2. **File hóa học chuẩn (`ligand_ideal.sdf`):** Bạn tải file cấu trúc "ideal" của ligand đó từ RCSB PDB. File này có cấu trúc hóa học 2D/3D chính xác nhưng tọa độ không nằm trong túi kết hợp (binding site).

### Bước 2: Chạy mã Python (RDKit) để tạo `crystal_ligand.sdf`
Đoạn mã sau sẽ lấy tọa độ từ file PDB và "ốp" các thuộc tính hóa học từ file Ideal sang, đồng thời thêm nguyên tử Hydro dựa trên tọa độ không gian đó:

```python
from rdkit import Chem
from rdkit.Chem import AllChem

# 1. Tải cấu trúc ideal (để lấy bậc liên kết) và cấu trúc PDB (để lấy tọa độ gốc)
# Tham số removeHs=True để tạm bỏ các Hydro chưa đồng bộ
ideal_mol = Chem.MolFromMolFile("ligand_ideal.sdf", removeHs=True)
pose_mol = Chem.MolFromPDBFile("ligand_fromPDB.pdb", removeHs=True)

# 2. Ánh xạ (Assign) bậc liên kết từ phân tử ideal sang tọa độ của pose_mol
corrected_pose = AllChem.AssignBondOrdersFromTemplate(ideal_mol, pose_mol)

# 3. Thêm nguyên tử Hydro vào cấu trúc đã sửa (dựa trên tọa độ 3D)
corrected_pose_with_H = Chem.AddHs(corrected_pose, addCoords=True)

# 4. (Tùy chọn) Kiểm tra tính toàn vẹn bằng cách so sánh SMILES
assert Chem.MolToSmiles(corrected_pose) == Chem.MolToSmiles(ideal_mol)

# 5. Lưu thành file SDF hoàn chỉnh sẵn sàng cho GNINA
writer = Chem.SDWriter("crystal_ligand_ready_for_redocking.sdf")
writer.write(corrected_pose_with_H)
writer.close()
```

### Bước 3: Đưa vào GNINA để Re-docking
Lúc này, file `crystal_ligand_ready_for_redocking.sdf` đã mang chính xác tư thế gốc của tinh thể. Bạn có thể cấp nó cho GNINA:
* Dùng làm ligand đầu vào: `-l crystal_ligand_ready_for_redocking.sdf`
* Dùng làm tâm tạo hộp lưới (Grid box) tự động: `--autobox_ligand crystal_ligand_ready_for_redocking.sdf`

Cách làm này tuân thủ tuyệt đối tiêu chuẩn học thuật từ **IQB-2025**, đảm bảo bạn đang kiểm tra khả năng hội tụ của GNINA về đúng tư thế tinh thể mà không bị nhiễu do lỗi sai bậc liên kết hay thiếu Hydro.
