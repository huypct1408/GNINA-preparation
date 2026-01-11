**1. Lời đầu tiên**

Trong các pipeline docking hiện đại (đặc biệt là GNINA – CNN-based docking), chất lượng ligand preparation có ảnh hưởng quyết định đến kết quả docking, đôi khi còn quan trọng hơn cả tham số sampling.

**OpenBabel** thường được sử dụng để:

+ Chuẩn hóa SMILES

+ Thêm hydrogen

+ Proton hóa theo pH

+ Sinh cấu trúc 3D

Tuy nhiên, OpenBabel sử dụng **heuristic-based protonation**, không đảm bảo tạo ra **dominant microspecies** ở pH sinh lý (7.4), đặc biệt với:

+ Heteroaromatic rings

+ Ligand có nhiều tautomer

+ Amino acid–like scaffolds

+ Charged groups (COO⁻, imidazole, pyridine, v.v.)

Để khắc phục điều này, em sử dụng **MolScrub** (theo gợi ý của https://github.com/MolSSI-Education/iqb-2025, workshop education của Protein Data Bank) làm lõi xử lý hóa học (thermodynamic fix), kết hợp với RDKit để sinh hình học 3D (geometric fix).

**2. MolScrub là gì?**

MolScrub là một thư viện xử lý ligand dựa trên **rule-based pKa reactions**, được thiết kế để:

+ Xác định trạng thái **proton hóa chiếm ưu thế (dominant microspecies)** tại một pH xác định

+ Enumerate **tautomers hợp lệ về mặt hóa học**

+ Loại bỏ counter-ions, salts

+ Chuẩn hóa formal charges theo nguyên lý hóa học

| Thuộc tính  | MolScrub                              |
| ----------- | ------------------------------------- |
| Protonation | Rule-based pKa (không heuristic)      |
| pH control  | Chính xác (ph_low = ph_high)          |
| Output      | RDKit Mol object                      |
| Focus       | Thermodynamics (chemical correctness) |

Trong file code ``chemistry 2(molscrub).py``, MolScrub được cấu hình như sau:
```
Scrub(
    ph_low=7.4,
    ph_high=7.4,
    skip_acidbase=False,
    skip_tautomers=False,
    skip_ringfix=False,
    skip_gen3d=True,
    keep_all_frags=False
)
```
Điều này để đảm bảo: 1 ligand → 1 trạng thái proton hóa → đúng tại pH 7.4

**3. Vì sao cần thay thế OpenBabel?**

  **3.1 Các hạn chế của OpenBabel**

    OpenBabel rất mạnh và linh hoạt, nhưng có các hạn chế sau:

    1. Protonation heuristic

    + Không dựa trên pKa thực

    + Có thể tạo nhiều microspecies không chiếm ưu thế
    
    2. Silent permissiveness
    
    + Chấp nhận nhiều cấu trúc hóa học “xấp xỉ đúng”

    + Dễ sinh ra trạng thái valence/charge không tối ưu cho docking

    3. Charge ambiguity

    + Các nhóm acid/base (COOH, imidazole, pyridine) đôi khi không ở trạng thái dominant

Trong flexible docking, những sai lệch nhỏ này bị khuếch đại mạnh do:

  + Protein flexibility

  + Entropy của ligand lớn

  + CNN scoring nhạy với electrostatics

  **3.2 MolScrub Advantages**

  MolScrub khắc phục trực tiếp các điểm yếu trên:

| Tiêu chí                | OpenBabel     | MolScrub       |
| ----------------------- | ------------- | -------------- |
| Protonation logic       | Heuristic     | Rule-based pKa |
| Dominant microspecies   | ❌             | ✅           |
| Chemical rigor          | Trung bình    | Cao            |
| CNN compatibility       | Không đảm bảo | Cao            |
| Docking reproducibility | Trung bình    | Cao            |

**4. Thiết kế**

4.1. Chia làm 2 giai đoạn 

🔹 Giai đoạn 1 — Thermodynamic Fix (MolScrub)
``mol = ph_correct_smiles_molscrub(smiles, ph=7.4)``


+ Xác định đúng charge

+ Chọn tautomer/protonation hợp lệ

+ Lọc bỏ state lỗi (valence, kekulization)

🔹 Giai đoạn 2 — Geometric Fix (RDKit)
``mol_3d = mol_to_3d_mol(mol)``


+ ETKDGv3 embedding

+ MMFF94s / UFF minimization

+ Explicit hydrogens

+ GNINA-compatible SDF

4.2 Chemical Safety Guards

Pipeline chủ động khóa chặt “silent failures” bằng các kiểm tra bắt buộc:

+ ``Chem.SanitizeMol(mol)``

+ Valence check (đặc biệt cho O, N)

+ Xác nhận real 3D geometry bằng Z-span
```
if max(zs) - min(zs) < 0.1:
    FAIL (2D molecule)
```


Code kiểm tra này là 1 cell ở dưới cell thực thi việc chuẩn bị protein và ligand trong ``main_pipeline.ipynb``.

**5. Kết quả kiểm tra**

<img width="600" alt="image" src="https://github.com/user-attachments/assets/53afec1f-3afe-4b93-8d7d-3e9679ab5a48" />

  **5.1 RDKit Sanitation**
    ``Sanitize FAILED    : 0``

      → 100% ligand hợp lệ hóa học theo RDKit

  **5.2 3D Geometry Check**
 ``` 
[LIG 1] Z-span = 7.606 Å
[LIG 5] Z-span = 6.724 Å
[LIG 8] Z-span = 8.117 Å
```


→ Không có ligand 2D
→ Không có “silent fallback” từ MolScrub

**6. Docking Performance: MolScrub vs OpenBabel**

Case flexible docking với GNINA: Cùng ligand (``3-[[2-(2,2-dimethyl-3,4-dihydrochromen-6-yl)acetyl]amino]propanoic acid``) – cùng protein (``PDB id: 8skl``) – cùng parameters – khác preparation

  **6.1. MolScrub-prepared ligand (Top pose)**

| Metric         | Value          |
| -------------- | -------------- |
| CNN pose score | **0.8683**     |
| CNN affinity   | **5.740**      |
| Affinity       | -8.16 kcal/mol |
```
Commandline: /kaggle/working/gnina -r /kaggle/working/docking_results/8skl/protein/receptor.pdb -l /kaggle/working/docking_results/8skl/ligands/LIG_0004__3-[[2-(2,2-dimethyl-3,4-dihydrochromen-6-yl)acetyl]amino]propanoic_acid/input/ligand.sdf --autobox_ligand /kaggle/working/docking_results/8skl/reference/ref_ligand.sdf --autobox_add 10 --flexres A:182,A:181,A:215,A:262,A:49 --num_modes 10 --exhaustiveness 32 --cnn_scoring rescore --cnn_empirical_weight 2.0 --pose_sort_order CNNscore --device 0 --seed 42 --atom_term_data -o /kaggle/working/docking_results/8skl/ligands/LIG_0004__3-[[2-(2,2-dimethyl-3,4-dihydrochromen-6-yl)acetyl]amino]propanoic_acid/output/docked.sdf --out_flex /kaggle/working/docking_results/8skl/ligands/LIG_0004__3-[[2-(2,2-dimethyl-3,4-dihydrochromen-6-yl)acetyl]amino]propanoic_acid/output/flex_residues.pdb --log /kaggle/working/docking_results/8skl/ligands/LIG_0004__3-[[2-(2,2-dimethyl-3,4-dihydrochromen-6-yl)acetyl]amino]propanoic_acid/logs/gnina.log
Flexible residues: A:49 A:181 A:182 A:215 A:262
Using random seed: 42

mode |  affinity  |  intramol  |    CNN     |   CNN
     | (kcal/mol) | (kcal/mol) | pose score | affinity
-----+------------+------------+------------+----------
    1       -8.16       -4.43       0.8683      5.740
    2       -8.80       -5.73       0.8094      5.552
    3       -4.91       -5.50       0.6639      4.561
    4       -4.82       -5.77       0.5908      3.772
    5       -5.16       -6.46       0.5882      4.263
    6       -5.16       -5.80       0.5143      4.248
    7       -6.37       -5.73       0.4571      4.745
    8       -7.77       -5.26       0.3375      4.855
    9       -5.15       -6.16       0.3275      3.943
   10       -5.13       -5.90       0.3224      4.006
```

**6.2. OpenBabel-prepared ligand (Top pose)**

| Metric         | Value          |
| -------------- | -------------- |
| CNN pose score | 0.7180         |
| CNN affinity   | 5.480          |
| Affinity       | -8.09 kcal/mol |

```
Commandline: /kaggle/working/gnina -r /kaggle/working/docking_results/8skl/protein/receptor.pdb -l /kaggle/working/docking_results/8skl/ligands/3-[[2-(2,2-dimethyl-3,4-dihydrochromen-6-yl)acetyl]amino]propanoic acid/input/ligand.sdf --autobox_ligand /kaggle/working/docking_results/8skl/reference/ref_ligand.sdf --autobox_add 10 --flexres A:182,A:181,A:215,A:262,A:49 --num_modes 10 --exhaustiveness 32 --cnn_scoring rescore --cnn_empirical_weight 2.0 --pose_sort_order CNNscore --device 0 --seed 42 --atom_term_data -o /kaggle/working/docking_results/8skl/ligands/3-[[2-(2,2-dimethyl-3,4-dihydrochromen-6-yl)acetyl]amino]propanoic acid/output/docked.sdf --out_flex /kaggle/working/docking_results/8skl/ligands/3-[[2-(2,2-dimethyl-3,4-dihydrochromen-6-yl)acetyl]amino]propanoic acid/output/flex_residues.pdb --log /kaggle/working/docking_results/8skl/ligands/3-[[2-(2,2-dimethyl-3,4-dihydrochromen-6-yl)acetyl]amino]propanoic acid/logs/gnina.log
Flexible residues: A:49 A:181 A:182 A:215 A:262
Using random seed: 42

mode |  affinity  |  intramol  |    CNN     |   CNN
     | (kcal/mol) | (kcal/mol) | pose score | affinity
-----+------------+------------+------------+----------
    1       -8.09       -4.65       0.7180      5.480
    2       -4.93       -6.37       0.5943      4.498
    3       -7.67       -5.51       0.5828      5.158
    4       -7.40       -5.64       0.4930      4.295
    5       -8.34       -4.65       0.4442      5.094
    6       -5.53       -5.90       0.3717      3.704
    7       -4.96       -4.65       0.3433      3.807
    8       -5.90       -5.82       0.3225      4.159
    9       -4.86       -6.23       0.3093      3.834
   10       -5.48       -6.26       0.3064      4.053
```
Theo quan sát của em, em nhận thấy: 
**MolScrub:**

+ CNN score cao hơn đáng kể

+ Pose hội tụ rõ ràng

+ Electrostatics nhất quán với pH 7.4

**OpenBabel:**

+ CNN score thấp hơn

+ Không có basin hội tụ rõ

+ Dễ false positive




  
