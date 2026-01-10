# GNINA-preparation
Protein and Ligand Preparation for GNINA flexible docking
'''
open_prep_pipeline/
├── config.yaml
├── main_pipeline.ipynb
├── input/
│   ├── protein.pdb
│   └── ligands.xlsx
├── output/
│   ├── protein_prepared.pdb
│   └── ligands_prepared.sdf
└── utils/
    ├── __init__.py
    ├── logger.py
    ├── protein_logic.py
    ├── ligand_logic.py
    └── chemistry.py
'''
Tổng quan open_prep_pipeline là một quy trình chuẩn bị có tính tái lập và tuân thủ các nguyên tắc hóa học, được thiết kế cho các nghiên cứu docking protein–ligand. Quy trình này chuẩn hóa quá trình proton hóa và sự toàn vẹn cấu trúc của protein, cũng như trạng thái proton hóa, hình học và tô pô (topology) của ligand trước khi tiến hành docking linh hoạt.

**1. Chuẩn bị protein (utils/protein_logic.py): Tạo ra cấu trúc protein hoàn chỉnh về mặt hóa học và được proton hóa, phù hợp cho docking linh hoạt. Việc chuẩn bị protein được thực hiện bằng PDBFixer và OpenMM:**
   **- Sửa chữa cấu trúc (Structural Repair):**
     + Xác định các gốc axit amin (residue) bị thiếu.
     + Khôi phục các nguyên tử nặng bị thiếu.
     + Hoàn thiện hóa trị.
  **- Proton hóa ở pH sinh lý (Physiological Protonation):**
    + Bổ sung hydro dựa trên pH của dung dịch (bulk solution pH).
    + Trạng thái proton hóa phản ánh các trạng thái chiếm ưu thế tại pH quy định (mặc định là 7.4).
  **- Bảo toàn tô pô (Topology Preservation):**
    + Giữ nguyên mã định danh (ID) gốc và nguyên tử ban đầu.
    + Đầu ra tương thích với các bộ máy docking phía sau.
**-> Đầu ra:** Định dạng: PDB (output/protein_prepared.pdb)

 ** Đặc điểm:**

    + Đầy đủ hóa trị (Valence-complete).

    + Đã hydro hóa (Hydrogenated).

    + Không áp đặt ràng buộc trường lực (force-field constraints).

    + Không áp dụng tối ưu hóa năng lượng (trách nhiệm này thuộc về các bộ máy docking).

  **2. Chuẩn bị ligand (utils/chemistry.py và utils/ligand_logic.py): Tạo ra các cấu trúc ligand đơn trạng thái (single-state), thực tế về mặt vật lý và nhất quán với điều kiện sinh lý.**
 ** - Gán trạng thái Proton hóa (Hiệu chỉnh nhiệt động lực học)**
    + Công cụ: OpenBabel
    + Nguyên tắc: Trạng thái proton hóa được xác định bởi sự chiếm ưu thế pKa–pH, không phải bằng phương pháp liệt kê.
    + Kết quả: Mỗi ligand được đại diện bởi một vi trạng thái (microspecies) chiếm ưu thế tại pH mục tiêu (mặc định 7.4). Điều này giúp tránh: Trung bình hóa điện tích nhân tạo, Docking các trạng thái ion không tồn tại và Sai lệch điểm số do tính chất tĩnh điện không chính xác. 
 ** - Tạo hình học 3D (Hiệu chỉnh hình học):**
    + Công cụ: RDKit (thuật toán nhúng ETKDGv3)
    + Quy trình:

        Thêm hydro sau khi hiệu chỉnh proton hóa.

        Tạo một cấu dạng (conformer) 3D duy nhất.

        Tối ưu hóa năng lượng sử dụng MMFF94 (dự phòng bằng UFF).
    ** - Chú giải và Khả năng truy xuất Mỗi ligand lưu trữ:**

      + Chuỗi SMILES gốc.

      + Chuỗi SMILES đã hiệu chỉnh theo pH.

      + Mã định danh ligand.
    -> Đầu ra: Định dạng SDF (output/ligands_prepared.sdf)

      Đặc điểm: Một ligand → một trạng thái proton hóa → một cấu dạng đã tối ưu hóa năng lượng.

+ Tệp config.yaml tập trung đảm bảo sự nhất quán về pH và tham số.

+ main_pipeline.ipynb là Jupyter Notebook chứa code thực thi việc chuẩn bị protein và ligand

+ File input bao gồm protein.pdb là file protein đã được tách ra khỏi phức hợp tinh thể tải về từ Protein Data Bank bằng PYMOL (đã loại nước, các ion không cần thiết), với ligand.xlsx là file chứa các hợp chất cần chuẩn bị, có 2 cột "SMILES" và "ID" (Trong file xlsx của em dùng cột ID có tên là "Cleaned_Name")
    
