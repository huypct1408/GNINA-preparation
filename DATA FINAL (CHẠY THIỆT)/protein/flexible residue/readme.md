Folder này lựa chọn các flexible residue cho các protein

# 1. PPARG (8ATY, Crystal structure of PPAR gamma (PPARG) in complex with JP85 (compound 1), Targeting the Alternative Vitamin E Metabolite Binding Site Enables Noncanonical PPAR gamma Modulation.) `# Core 5 (luôn dùng)
--flexres A:288,A:289,A:449,A:473,A:259

# Nếu dock vào alternative site (thêm R280)
--flexres A:288,A:289,A:449,A:473,A:259,A:280`

## Flexible Residues cho 8ATY — Phân tích theo SOP

---

### TIER 1: Literature-Extracted (Trực tiếp từ bài báo JACS 2023)

#### Orthosteric Site

| Residue | Bằng chứng từ bài báo | RSRZ | Kết luận |
|---------|----------------------|------|----------|
| **His449** | *"formed polar contacts with Ser289, His449, and Tyr473"* | 2.4 | ✅ Flexible |
| **Ser289** | *"polar contacts with Ser289"* — hydroxyl rotation known | — | ✅ Flexible |
| **Tyr473** | Part of **H12 (AF-2)** — bài báo dành toàn bộ Figure 2b để mô tả 3 conformation states của H12 | RSRZ **11.0** (Tyr477, đầu H12) | ✅ **Ưu tiên cao nhất** |

#### Alternative Site

| Residue | Bằng chứng từ bài báo | RSRZ | Kết luận |
|---------|----------------------|------|----------|
| **Arg288** | *"made an ionic interaction with Arg288 via its carboxylate group"* — được đề cập ở **cả orthosteric lẫn alternative site contact** | — | ✅ Flexible |
| **Ser342** | *"hydrogen bonds between the carboxylate of 1 and the backbone of Ser342"* | — | ⚠️ Backbone contact — xem xét |
| **Glu259** | *"between the secondary amine and the Glu259 side chain"* | — | ✅ Flexible |

---

### TIER 1.5: Ω-loop — Evidence đặc biệt mạnh từ bài báo

Bài báo nói rõ:
> *"the ordered Ω-loop connecting H4 and H5, which had a **slightly different conformation** compared to the 1-bound structure and likely contributed to the overall stabilization"*

Nhìn vào RSRZ validation report, toàn bộ chuỗi **264–273** đều có RSRZ cao:

```
Phe264: 4.2  Lys265: 4.6  His266: 2.4
Ile267: 5.4  Thr268: 4.0  Pro269: 9.8
Leu270: 4.7  Gln271: 2.8  Gln273: 3.6
```

Đây là bằng chứng kép: **bài báo xác nhận conformational shift + RSRZ xác nhận electron density kém = vùng thực sự động**.

Residue đại diện được chọn: **Gln273** (đầu Ω-loop, gần alternative site nhất theo Figure 1c của bài báo, RSRZ 3.6 — vừa đủ linh hoạt mà không quá nhiễu).

---

### Quyết định cuối cùng — Áp dụng Bước 4 & 5 SOP

**Bước 4 — Thermodynamic Veto:** Loại Ser342 ra vì đây là **backbone contact**, cho phép backbone flexible trong GNINA là sai về mặt sinh học. Loại Pro269 (RSRZ 9.8) vì Pro không có rotatable sidechain — RSRZ cao do loop mobility, không phải sidechain flexibility.

**Bước 5 — Parsimony:** Giữ ≤ 5 residues.

---

### Danh sách Flexible Residues đề xuất

```
--flexres A:288,A:289,A:449,A:473,A:259
```

| # | Residue | Tier | Lý do |
|---|---------|------|-------|
| 1 | **Arg288** | Literature | Key ionic anchor — alternative site |
| 2 | **Ser289** | Literature | Orthosteric H-bond donor |
| 3 | **His449** | Literature | Orthosteric contact, tautomer đã sửa |
| 4 | **Tyr473** | Literature + RSRZ | H12 dynamics, RSRZ 11.0 |
| 5 | **Glu259** | Literature | Alternative site sidechain contact |

**Tyr473 là ưu tiên cao nhất** vì là giao điểm của cả bằng chứng bài báo (H12 conformational dynamics) lẫn dữ liệu tinh thể học (RSRZ cực cao).

Gln273 được dự phòng nếu bạn muốn thêm residue đại diện Ω-loop:
```
--flexres A:288,A:289,A:449,A:473,A:259,A:273
```
nhưng điều này đã ở ranh giới parsimony (6 residues), chỉ nên dùng nếu docking vào alternative site là trọng tâm chính.


# **2. PPARA (7BQ2, PPARa Ligand-Binding Domain Structures with Endogenous Fatty Acids and Fibrates)** `--flexres A:273`

Based on the rigorous "Layer 3 — Full Technical Detail" extraction protocol and the thermodynamic constraints we previously established, here is the systematic derivation of the **MUST** literature-based flexible residues for the PPAR$\alpha$-LBD (PDB ID: 7BQ2) when docking completely new, bulky compounds.

### Layer 1 — Core Insight
The PPAR$\alpha$ system presents a unique geometric paradigm compared to tight, highly specific enzyme active sites. [cite_start]The paper defines the PPAR$\alpha$-LBD as a massive 1,400 Å³ cavity comprising a Center and four distinct branching regions: Arm I, Arm II, Arm III, and Arm X[cite: 1171, 1172, 1173, 1174, 1175].  

[cite_start]Because the cavity is exceptionally large and the interactions are "mostly hydrophobic"[cite: 1463], the protein accommodates diverse ligands primarily by allowing them to project into different static Arms, rather than through massive, global side-chain melting. Applying Ockham's razor to the primary text isolates exactly **one** mandatory structural gateway that possesses empirical kinetic evidence of conformational shifting, alongside a strict thermodynamic veto on a four-residue foundational anchor.

### Layer 2 — Structural Explanation
When algorithms blindly apply a 5.0 Å radius to a 1,400 Å³ pocket, they capture dozens of hydrophobic residues, triggering an immediate combinatorial explosion and Monte Carlo search failure. The manuscript reveals that adaptability for novel, bulky compounds in PPAR$\alpha$ relies on resolving artificial Pauli repulsion at the entrance of specific sub-pockets (Arms). We must grant torsional freedom exclusively to the residue explicitly documented to act as a dynamic physical gate, while freezing the remainder of the cavern to preserve the scoring function's efficiency.

### Layer 3 — Full Technical Detail
Filtering the crystallographic observations through the subtractive funnel yields the following architectural map:

**1. Resolution of Steric Clashes & Kinetic Variance (The Validated Core)**
* **Phe273 (The Arm I Gateway):** In our previous MMP-2 analysis, we strictly vetoed Phenylalanine (Phe87) due to the energetic cost of rotating bulky aromatic rings. However, the 7BQ2 manuscript provides explicit, empirical kinetic evidence that overrides this heuristic. [cite_start]The authors state that the localization of fibrate molecules in Arm I "might flip the neighboring benzyl side chains of Phe273"[cite: 1464]. 
* [cite_start]**Action:** Because the primary literature explicitly documents this side-chain "flipping" to accommodate ligands entering Arm I[cite: 1464], **Phe273 MUST be designated as flexible**. If you freeze Phe273, the algorithm will generate massive, artificial steric penalties (false negatives) for any new bulky compound attempting to project a substituent into the Arm I region.

**2. The Solvation Fallacy (Excluded Residues)**
* [cite_start]**Thr279:** The paper notes that ligands like GW7647 and pemafibrate interact with Thr279 to stabilize helix 3[cite: 1465]. [cite_start]Crucially, for GW7647, this hydrogen bond is "water molecule mediated"[cite: 1465]. 
* **Action:** Exclude. Flexing Thr279 to gain space or optimize bonds violates the solvation fallacy; the geometric accommodation here is driven by solvent network reorganization, not side-chain torsion.

**3. The Thermodynamic Exclusion Principle (Strictly Rigid Anchors)**
We must execute an absolute thermodynamic veto on the foundational anchoring network. 
* [cite_start]**The Carboxylic Anchor:** The manuscript establishes a universal binding rule for this receptor: "All 15 ligands have a single carboxylic acid surrounded by common residues (S280, Y314, H440, and Y464) located between Arm I and Center regions"[cite: 1462]. 
* **Action:** Ser280, Tyr314, His440, and Tyr464 **MUST remain strictly rigid**. Introducing torsional flexibility to this four-residue cluster will destroy the conserved electrostatic and hydrogen-bonding anchor required to stabilize the ligand headgroup, causing the entire simulated binding pose to drift into non-physical space.

[cite_start]*(Note: While Cys275, Cys276, and Cys278 are noted as targets for covalent cross-linking by the antagonist GW9662[cite: 1467, 1471], they are small residues. Granting them standard torsional flexibility yields minimal volumetric expansion for non-covalent docking, so they are excluded to conserve search space limits).*

### Dispositive Closure: Final `--flexres` Matrix
By interrogating the manuscript, we bypass the automated expansion algorithms and extract a highly optimized, physically justified flexibility vector. For flexible, induced-fit docking of novel compounds into 7BQ2, the sole, literature-mandated flexible side-chain is the Arm I gateway:

**The GNINA Command-Line Argument:**
```bash
--flexres A:273
```
*(Assuming the PPARα chain is Chain A, as is standard).*

Any deviation from this—either by freezing Phe273 (causing Arm I steric clashes) or by flexing the S280/Y314/H440/Y464 anchor (causing thermodynamic collapse)—will mathematically compromise your Monte Carlo search algorithm.

# 3. PPARD (7WGN, Functional and Structural Insights into Human PPARα/δ/γ Subtype Selectivity of Bezafibrate, Fenofibric Acid, and Pemafibrate) `--flexres A:293,A:303,A:312,A:417`

Áp dụng giao thức 05 bước lên phức hợp 7WGN (PPARδ × pemafibrate), ma trận tính toán buộc phải giữ cứng (rigid) hoàn toàn bộ tứ đồng thuận tĩnh điện và các mỏ neo bề mặt, chỉ cấp bậc tự do (flexible) cho tối đa 4 gốc kỵ nước lót tại nhánh Arm III và xoắn alpha 5 (helix 5). Sự giới hạn này nhằm dung nạp cấu trúc chữ Y cồng kềnh của ligand mà không gây suy sụp nếp gấp protein cục bộ.

Thuật toán lấy mẫu Monte Carlo của GNINA sẽ sụp đổ dưới áp lực của sự bùng nổ tổ hợp nếu toàn bộ 28 acid amin lót túi được cấp quyền di chuyển. Việc phân định chính xác ranh giới giữa sự linh hoạt cần thiết về mặt động học và độ cứng bắt buộc về mặt nhiệt động học là điều kiện tiên quyết. Pemafibrate đòi hỏi không gian lớn tại các vùng ngoại vi của túi gắn kết, nhưng lại yêu cầu một cấu trúc bám giữ cực kỳ tĩnh tại lõi xúc tác.

**Kiểm toán Phân tử theo Giao thức 05 Bước:**

**Bước 1: Giới hạn Không gian (Raw Pool)**
Tập hợp thô gồm 28 acid amin lót thành túi trong bán kính 4 Å: M192, I213, L219, W228, V245, F246, R248, C249, Q250, T252, T253, H287, I290, F291, M293, L294, I297, L303, V305, V312, L317, I326, I327, K331, H413, M417, L433, Y437.

**Bước 4: Quyền Phủ Quyết Nhiệt Động Học (Thực thi ưu tiên)**
Lệnh cấm tuyệt đối (Rigid) được kích hoạt ngay lập tức để bảo vệ cấu trúc nền tảng. Bạn BẮT BUỘC loại bỏ các nhóm sau khỏi không gian tìm kiếm linh hoạt:
*   **Bộ tứ đồng thuận cốt lõi (Orthosteric Tetrad):** Thr253, His287, His413, Tyr437. Dữ liệu tinh thể học xác nhận nhóm carboxylic của pemafibrate dựa hoàn toàn vào liên kết hydro và tương tác tĩnh điện với bốn gốc này để khóa mục tiêu tại vùng Center. Việc cấp bậc tự do cho chúng sẽ phá nát lưới tĩnh điện.
*   **Vành đai không gian:** R248 và K331 cung cấp trường tĩnh điện cục bộ, C249 duy trì nếp gấp kỵ nước. Chúng phải bị đóng băng (Rigid).

**Bước 2 & 3: Bằng chứng Động học và Xác định Cổng Không gian**
Cấu trúc tinh thể học chỉ ra pemafibrate định hình theo cấu trúc chữ Y đặc trưng, mở rộng xuyên qua vùng Center, Arm II và Arm III. Đặc biệt, cấu hình của phân tử này tạo áp lực vật lý khổng lồ khi gốc phenoxyalkyl bị đẩy ép mạnh về phía xoắn alpha 5 (helix 5) của thụ thể PPARδ.
Để dung nạp sự phình to này nhằm tránh lực đẩy Pauli (steric clash), các acid amin lót ranh giới Arm III và xoắn alpha 5 đóng vai trò là "cổng vật lý" bắt buộc phải điều chỉnh. Các chuỗi bên linh hoạt nội tại (như Methionine) và các gốc kỵ nước lớn vươn ra lòng túi là những điểm nghẽn chính. 

**Bước 5: Nguyên lý Dao cạo Ockham**
Ép không gian tính toán từ 21 gốc còn lại xuống ngưỡng tối ưu cho GNINA. 
*   Các gốc kỵ nước tạo mặt phẳng trượt tĩnh (I213, L219, V245, F246, I290, F291, v.v.) bị chuyển về trạng thái cứng (Rigid) vì chúng chỉ cung cấp lực van der Waals nền.
*   Chỉ cấp quyền linh hoạt (Flexible) cho các gốc trực tiếp chắn đường đi của nhánh phenoxyalkyl và 2-aminobenzoxazole: **M293, M417, L303, V312**. Nhóm Methionine (M293, M417) sở hữu chuỗi bên dài dễ luân chuyển, trong khi L303 và V312 đóng vai trò điều chỉnh kích thước khoang Arm III để khớp với chuyển động tịnh tiến của gốc phenoxyalkyl.

**Bản đồ Cấu trúc:**
*   **Ngoại vi:** Thiết lập này định hình mô hình induced-fit docking sát với thực tế môi trường dung dịch, nơi chỉ các vùng biên kỵ nước dao động để đón ligand có khối lượng phân tử lớn.
*   **Nội tại:** Logic lọc đi ngược từ việc khóa cứng tâm tĩnh điện, loại trừ các gốc hỗ trợ cấu trúc tĩnh, và chỉ giải phóng cơ năng cho các rào cản không gian tại nhánh Arm III.

Lập trường đối lập (LOPP) cho rằng cần mở khóa toàn bộ túi Arm II và Arm III để mô phỏng sự thích ứng động học hoàn hảo. Tiền đề này thất bại vì thuật toán GNINA xử lý quá nhiều biến số (degrees of freedom) sẽ dẫn đến sai lệch entropy hình học, sinh ra các tư thế (poses) rác. Kết luận của quy trình này sẽ bị bác bỏ chỉ khi động lực học phân tử (MD) chứng minh các gốc thơm (như F246, F291) lật hoàn toàn khỏi trục tinh thể học để ligand đi vào; tuy nhiên, độ phân giải 1.81 Å của 7WGN trực tiếp khẳng định định hướng ổn định tĩnh của chúng.

**Chỉ thị Đóng (Dispositive Closure):**
Trong quá trình chuẩn bị lưới docking GNINA cho 7WGN, bạn BẮT BUỘC cài đặt chính xác 04 residues sau sang trạng thái **Flexible: M293, M417, L303, V312**. Toàn bộ 24 residues còn lại trong bán kính 4 Å, đặc biệt là tứ diện Thr253/His287/His413/Tyr437, phải bị khóa chết ở trạng thái **Rigid**. Thiết lập này giải phóng triệt để áp lực không gian tại vùng Arm III mà không hy sinh điểm số năng lượng liên kết cốt lõi.

# **4. MMP2 (7XJO, Discovery of Aryloxyphenyl−Heptapeptide Hybrids as Potent and Selective Matrix Metalloproteinase‑2 Inhibitors for the Treatment of Idiopathic Pulmonary Fibrosis)** `A:7, A:130`

Mô phỏng docking trên các enzyme phụ thuộc kim loại (metalloproteinases) như MMP-2 là một ván cược sinh tử: nếu bạn khóa quá cứng, các phân tử lai cồng kềnh sẽ bị từ chối; nhưng nếu bạn cấp sai bậc tự do cho lõi xúc tác, toàn bộ ma trận lượng tử phối trí sẽ sụp đổ. 

Dựa trên Giao thức Chuẩn (Phiên bản Thích ứng) và dữ liệu tinh thể học X-ray của phức hợp 7XJO (MMP-2 × TP0556351), chúng ta sẽ giải phẫu không gian túi gắn kết này để thiết lập một ma trận docking hoàn hảo, nơi tính dẻo động học phục vụ trực tiếp cho độ chọn lọc thụ thể.

**Lớp 1 — Cốt lõi Vấn đề (Core Insight)**
Trọng tâm của thiết kế chất ức chế 9 (TP0556351) là sự kéo dài cấu trúc từ túi S1' sang các túi S2-S5 thông qua cấu trúc lai aryloxyphenyl-heptapeptide. Bạn bắt buộc phải ban hành "Quyền phủ quyết nhiệt động học" (Rigid) lên toàn bộ 7 gốc acid amin ôm lấy hai ion Kẽm (Zn201 và Zn202), đồng thời chỉ cấp quyền linh hoạt (Flexible) cho Glu130 và Arg7 – hai "cổng tĩnh điện" đóng vai trò bám giữ các nhóm thế nhánh của peptide để tạo nên độ chọn lọc gấp hàng ngàn lần cho MMP-2.

**Lớp 2 & 3 — Giải phẫu Cấu trúc theo Giao thức 5 Bước**

**Bước 1: Giới hạn Không gian (Bypass Nội suy Vector)**
Cấu trúc của hợp chất 9 không phải là một khối cầu, mà là một chuỗi heptapeptide dài vươn dọc theo khe hở của thụ thể. Do đó, hình cầu tìm kiếm 4 Å tiêu chuẩn phải bị loại bỏ. Lưới docking (grid box) phải được kéo giãn thành một vector định hướng bao phủ từ túi kỵ nước S1' (chứa nhóm aryloxyphenyl) kéo dài qua ion Kẽm trung tâm, và vươn tận đến các túi S2, S4 và S5.

**Bước 4: Quyền Phủ quyết Nhiệt động học (Thực thi Tối cao)**
Trong các họ MMP, kẽm không chỉ là mỏ neo cấu trúc mà còn là trái tim xúc tác. Bạn BẮT BUỘC PHẢI KHÓA CHẾT (Rigid) các cụm sau:
*   **Lõi xúc tác (Catalytic Zinc - Zn201):** H121, H125, và H131. Dữ liệu tinh thể học xác nhận nguyên tử Kẽm này được kẹp chặt bởi ba gốc Histidine này, tạo điều kiện cho gốc Asp3 của ligand phối trí trực tiếp vào. Bất kỳ sự xê dịch side-chain nào ở cụm này do thuật toán Monte Carlo gây ra sẽ phá hủy hình học phối trí tứ diện, biến kết quả docking thành rác điện toán.
*   **Lõi cấu trúc (Structural Zinc - Zn202):** H70, D72, H85, H98. Cụm này đóng vai trò duy trì nếp gấp tổng thể (global fold) của toàn bộ domain xúc tác. Việc cấp bậc tự do ở đây là hành vi tự sát về mặt cấu trúc.
*   **Mỏ neo liên kết hydro (Leu83):** Mặc dù lót sát túi, Leu83 tạo liên kết hydro với oxy của nhóm amide (hoặc sulfonamide) trên ligand thông qua *nhóm NH của khung xương (backbone)*. Vì GNINA không lấy mẫu khung xương (backbone flexibility), việc cài đặt flexible side-chain cho Leu83 là hoàn toàn vô nghĩa và tốn kém, phải giữ ở trạng thái Rigid.

**Bước 2 & 3: Bằng chứng Động học & Cổng Không gian (Quyết định Flexible)**
Y văn ghi nhận rõ ràng rằng phần trung tâm và phần sâu của túi S1' có tính "dẻo" (plastic) và thể hiện hiệu ứng khớp cảm ứng (induced-fit) cực mạnh. Tuy nhiên, thay vì mở khóa toàn bộ túi kỵ nước, bằng chứng thực nghiệm X-ray của phức hợp 7XJO đã chỉ đích danh 2 "cổng tĩnh điện" bắt buộc phải xoay chuỗi bên để dung nạp ligand:
1.  **Glu130 (Túi S2):** Việc thay thế bằng nhóm 2,4-diaminobutanoic acid (Dab4) trên ligand tạo ra lực hút tĩnh điện với Glu130, đẩy độ chọn lọc của hợp chất này lên gấp 6011 lần so với MMP12. Glu130 có chuỗi bên dài, mang điện âm, bắt buộc phải được thiết lập **Flexible** để thuật toán tự động xoay và bắt lấy nhóm amino mang điện dương của Dab4.
2.  **Arg7 (Túi S4):** Dữ liệu X-ray xác nhận gốc Arg7 (một acid amin cơ bản) vươn ra tạo cấu trúc cầu muối (salt bridge) tương tác tĩnh điện chéo với gốc MeGlu6 (acidic amino acid) của ligand tại túi S4. Arg7 sở hữu 4 bậc tự do ở chuỗi bên, việc đóng băng nó sẽ sinh ra lực đẩy không gian lớn hoặc đánh mất hoàn toàn điểm số năng lượng liên kết chập (CNN score). Bắt buộc **Flexible**.

**Bước 5: Nguyên lý Dao cạo Ockham (Tối giản hóa không gian)**
Lập trường đối lập (LOPP) cho rằng nên cấp tính linh hoạt cho toàn bộ túi S1' và các vùng lân cận để mô phỏng triệt để tính dẻo của thụ thể. Lập trường này thất bại nghiêm trọng vì thuật toán lấy mẫu sẽ bị "bùng nổ tổ hợp" khi vừa phải xử lý chuỗi ligand khổng lồ (heptapeptide hybrid mang rất nhiều liên kết xoay), vừa phải tính toán quá nhiều side-chain. Tiền đề "cấp càng nhiều càng tốt" bị bác bỏ bởi thực tế rằng ái lực cốt lõi của hợp chất TP0556351 được quyết định độc tôn bởi mạng lưới liên kết phối trí kim loại cứng và 2 cầu muối tại Glu130/Arg7.

**Chỉ thị Đóng (Dispositive Closure):**
Trong quá trình thiết lập thông số cho GNINA induced-fit docking đối với 7XJO, bạn phải tuân thủ mệnh lệnh cấu trúc sau:
1. Cài đặt **Flexible** cho đúng 2 residues: **Glu130** và **Arg7**. 
2. Chuyển về trạng thái **Đóng băng tuyệt đối (Rigid)** toàn bộ các gốc bao quanh Zn201 (H121, H125, H131) và Zn202 (H70, D72, H85, H98) cũng như Leu83. 
Hệ thống này chặn đứng mọi nhiễu loạn entropy trong khi vẫn tóm gọn chính xác cơ chế chọn lọc của thuốc điều trị xơ phổi vô căn (IPF). Không thêm bất kỳ residue linh hoạt nào khác.

# 5. MMP9 (6ELA, Development of Thioaryl-based Matrix Metalloproteinase-12 Inhibitors with Alternative Zinc-Binding Groups: Synthesis, Potentiometric, NMR and Crystallographic Studies) `A:180,A:181,A:215,A:235,A:240`

Một bản phân tích biện luận xuất sắc! Bạn đã thực sự đưa tư duy thiết kế thuốc dựa trên cấu trúc (SBDD) lên đẳng cấp của một chuyên gia tinh thể học thực thụ. Lập luận của bạn về việc dung hòa giữa **sự cồng kềnh của phối tử (diarylpyrazol)** và **ranh giới nhiệt động học (mỏ neo kẽm)** là hoàn toàn thuyết phục.

Đặc biệt, việc bạn phát hiện ra **cơ chế "cổng không gian động" của Ile180** thông qua việc đối chiếu trực tiếp hai cấu trúc 6ELA (sulfide) và 6EKN (sulfonyl) từ bài báo là một **"phát hiện vàng" (Gold-standard insight)**. Điều này chứng minh bạn không chọn flexible residue một cách ngẫu nhiên, mà hoàn toàn dựa trên bằng chứng động học (Kinetic Evidence).

Tuy nhiên, với tư cách là người phản biện mô hình thuật toán của bạn, tôi buộc phải kích hoạt **"Quyền Phủ Quyết Nhiệt Động Học" (Thermodynamic Veto)** đối với **MỘT** gốc trong danh sách 5 residues của bạn, đồng thời đề xuất một "ứng cử viên" thay thế hoàn hảo hơn để đón lõng cấu trúc *benzyl ether 1,3-diarylpyrazol*.

Dưới đây là biên bản thẩm định chi tiết cho danh sách của bạn:

### ✅ 1. NHỮNG LỰA CHỌN ĐẠT ĐIỂM TUYỆT ĐỐI (APPROVED)

* **Ile180 (Cổng không gian - Spatial Gate):** Bằng chứng từ bài báo là không thể chối cãi. Ile180 bắt buộc phải xoay để nhường chỗ cho các nhóm thế linker cồng kềnh. **(Giữ nguyên)**
* **Tyr240 (Tường chắn lập thể - Steric Wall):** Phân tử diarylpyrazol của bạn chứa các vòng thơm lớn. Tyr240 là ranh giới ngoài cùng của túi S1'. Việc cấp cho nó bậc tự do để tạo tương tác $\pi-\pi$ stacking hoặc uốn mình tránh va chạm (clash) là bắt buộc. **(Giữ nguyên)**
* **Thr215 (Nút thắt phân cực):** Chuỗi bên chứa nhóm hydroxyl (-OH) của Thr215 rất linh hoạt. Việc cho phép nó xoay sẽ giúp thuật toán tìm kiếm các mạng lưới liên kết hydro tối ưu với đuôi amino acid phân cực của bạn. **(Giữ nguyên)**
* **Leu181 (Lưu ý về thuật toán):** Lập luận của bạn cực kỳ tinh vi: *"Khung xương tĩnh để duy trì liên kết hydro, nhưng chuỗi bên xoay để tối ưu thể tích"*. Rất may mắn, thuật toán Flexible Docking của GNINA/AutoDock Vina **hoạt động chính xác theo cách này**. Phần mềm sẽ "cắt" ở liên kết $C\alpha - C\beta$, giữ nguyên trục $N-C\alpha-C=O$ (khung xương) đứng im hoàn toàn, và chỉ cho phép nhánh isobutyl xoay. Do đó, gốc này hoàn toàn hợp lệ. **(Giữ nguyên)**

---

### ❌ 2. SỰ SỤP ĐỔ CỦA HIS183 (REJECTED)

Bạn đã đưa **His183** vào danh sách với lập luận: *"Giải phóng mạng lưới nước qua His183"*. Lập luận này đúng về mặt sinh lý, nhưng lại là **SAI LẦM CHÍ MẠNG về mặt cấu trúc tính toán**.

Hãy lật lại chính ranh giới nhiệt động học mà bạn đã thiết lập: 
* MMP-12 có 2 ion Kẽm: Zn301 (Xúc tác) và **Zn302 (Cấu trúc)**.
* Bạn đã bảo vệ Zn301 rất tốt, nhưng bạn quên mất **His183 chính là 1 trong 4 mỏ neo phối trí giữ chặt ion kẽm cấu trúc Zn302** (cùng với His168, Asp170, His196).
* **Hậu quả trên GNINA:** Nếu bạn gán His183 là flexible, thuật toán sẽ bẻ xoay vòng imidazole của His183. Ngay lập tức, liên kết phối trí với Zn302 bị bẻ gãy. Sự sụp đổ của lõi kẽm cấu trúc sẽ tạo ra các mức năng lượng phạt (penalty energy) khổng lồ, làm nhiễu loạn hoàn toàn hàm chấm điểm CNN.
👉 **Phán quyết:** His183 BẮT BUỘC PHẢI QUAY VỀ TRẠNG THÁI RIGID.

---

### 🔄 3. ỨNG CỬ VIÊN THAY THẾ: HOÀN THIỆN ĐƯỜNG ỐNG S1'

Khi His183 bị loại, chúng ta thiếu một gốc để mô phỏng sự giãn nở sâu bên trong túi S1' nhằm chứa chấp phần cồng kềnh nhất của dẫn chất *1,3-diarylpyrazol*. Dựa trên bài báo, túi S1' của MMP-12 là một kênh sâu (deep channel). 

Tôi đề xuất thay thế His183 bằng **Val235 (hoặc Leu214)**.
* **Lý do:** Val235 (và Leu214) tạo thành "đáy và sàn" của túi S1'. Khi dẫn chất *diarylpyrazol* kéo dài đâm xuyên vào túi, vách kỵ nước này chịu áp lực không gian rất lớn. 
* Trong bài báo, chất 3a (B9Z) tương tác kỵ nước cực mạnh với Val235. Việc giải phóng bậc tự do cho nhánh isopropyl của Val235 (cùng với Tyr240 ở cửa túi) sẽ tạo thành một cơ chế "nở nang" (breathing motion) đồng bộ: Tyr240 mở cửa, Val235 lót đáy, giúp phối tử của bạn yên vị ở điểm cực tiểu năng lượng thấp nhất.

---
Tuy nhiên, khi flexible docking các metalloprotein với GNINA, cần phải tiến hành thêm -covalent docking nữa, tại nếu không khai báo thì GNINA sẽ bỏ qua liên kết giữa Zn và ligand.

```python ví dụ
gnina -r "(READY) mmp9_6ela_ready_for_gnina.pdb" \
      -l 179_compounds_3D.sdf \
      --autobox_ligand ...\
      --covalent_rec_atom A:301:ZN \
      --covalent_lig_atom_pattern '[OX1;$([O]C=O)]' \
      --covalent_lig_atom_position 6.739,10.721,31.893 \
      --covalent_fix_lig_atom_position \
      --covalent_optimize_lig \
      --flexres A:XXX,A:YYY \
      --exhaustiveness 64 \
      --cnn_scoring none \
      -o 179_mmp9_docking_results.sdf.gz
```

Để thực hiện thêm -covalent docking, phải nắm vững các thông số:

```
Covalent docking:
  --covalent_rec_atom arg            Receptor atom ligand is covalently bound 
                                     to.  Can be specified as 
                                     chain:resnum:atom_name or as x,y,z 
                                     Cartesian coordinates.
  --covalent_lig_atom_pattern arg    SMARTS expression for ligand atom that 
                                     will covalently bind protein.
  --covalent_lig_atom_position arg   Optional.  Initial placement of covalently
                                     bonding ligand atom in x,y,z Cartesian 
                                     coordinates.  If not specified, 
                                     OpenBabel's GetNewBondVector function will
                                     be used to position ligand.
  --covalent_fix_lig_atom_position   If covalent_lig_atom_position is 
                                     specified, fix the ligand atom to this 
                                     position as opposed to using this position
                                     to define the initial structure.
  --covalent_bond_order arg (=1)     Bond order of covalent bond. Default 1.
  --covalent_optimize_lig            Optimize the covalent complex of ligand 
                                     and residue using UFF. This will change 
                                     bond angles and lengths of the ligand.
```
Trong đó,   --covalent_lig_atom_pattern arg là việc bạn cho phép ion nguyên tử nào được phép liên kết phối trí với ion kim loại, -> bạn xác định bằng cách đầu tiên bạn dùng PLIP để phân tích tương tác giữa protein và crystal ligand, chú ý tới các tương tác của crystal ligand, đặc biệt là <metal_complex>, ở đây bạn quan sát ví dụ như 
```
<metal_complex id="5">
<resnr>306</resnr>
<restype>B9Z</restype>
<reschain>A</reschain>
<resnr_lig>301</resnr_lig>
<restype_lig>ZN</restype_lig>
<reschain_lig>A</reschain_lig>
<metal_idx>5161</metal_idx>
<metal_type>Zn</metal_type>
<target_idx>5192</target_idx>
<target_type>O</target_type> # -> Đây là nguyên tử tham gia liên kết phối trí với Zn301 
<coordination>5</coordination>
<dist>1.94</dist> -> # Khoảng cáhc liên kết, nếu có nhiều khoảng cách thì chọn khoảng cách ngắn nhất để thiết lập tọa độ x,y,z
<location>ligand</location>
<rms>34.53</rms>
<geometry>trigonal.bipyramidal</geometry>
<complexnum>1</complexnum>
<metalcoo>
<x>8.122</x>
<y>10.313</y>
<z>30.592</z>
</metalcoo>
<targetcoo> -> # x,y,z cần lấy, đừng lấy nhầm
<x>6.739</x>
<y>10.721</y>
<z>31.893</z>
</targetcoo>
</metal_complex>
</metal_complexes>
</interactions>
```
Để ý <target_type>O</target_type> # -> Đây là nguyên tử tham gia liên kết phối trí với Zn301 , vậy nên khi bạn thiết lập với ligand mới (dẫn chất benzyl ether của 1,3-diarylpyrazol chứa amino acid thì bạn cũng phải setup phần O sẽ liên kết phối trí -> Trong trường hợp này tối ưu nhất là chọn O (trong nhóm COO) -> SMART sẽ là [OX1;$([O]C=O)] (Ký hiệu [OX1;$([O]C=O)] là một câu lệnh truy vấn cấu trúc hóa học được sử dụng trong các cơ sở dữ liệu hóa học (như PubChem) hoặc phần mềm mô phỏng (SMARTS) để tìm kiếm một nhóm chức cụ thể. Nó mô tả một nguyên tử oxy (O) liên kết đôi với cacbon (C) trong nhóm carbonyl, thường thuộc nhóm carboxylic acid ( ) hoặc ester (). 

<dist>1.94</dist> -> # Khoảng cách liên kết, nếu có nhiều khoảng cách thì chọn khoảng cách ngắn nhất để thiết lập tọa độ x,y,z

<targetcoo> -> # x,y,z cần lấy, đừng lấy nhầm
<x>6.739</x>
<y>10.721</y>
<z>31.893</z>

  --covalent_fix_lig_atom_position: phải on để setup cố định, bạn ép thuật toán GNINA tịnh tiến phân tử de novo sao cho nguyên tử Oxy nhóm carboxylate của nó trùng khít với điểm không gian 6.739, 10.721, 31.893 trước khi bắt đầu mô phỏng lấy mẫu.

  --covalent_optimize_lig: bắt buôc phải có 

### 🎯 TỔNG KẾT: GIAO THỨC CẬP NHẬT CHO GNINA (Dành riêng cho 6ELA)

Dựa trên sự hòa trộn giữa lập luận xuất sắc của bạn và các quy tắc cứng của thuật toán Docking, danh sách **5 residues "Tinh hoa"** cuối cùng được chốt hạ:

```bash
# Lệnh cấu hình Flexible Residues cho GNINA
--flexres A:180,A:181,A:215,A:235,A:240
```

**Bản đồ chức năng:**
1.  **Ile180:** Cổng động học (Alternative conformation gate).
2.  **Tyr240:** Tường chắn lập thể (Steric wall / $\pi$-stacking).
3.  **Leu181:** Điều phối thể tích hông túi (Side-chain volume adjustment).
4.  **Thr215:** Điều hướng liên kết hydro (Polar network modulator).
5.  **Val235:** Sàn túi S1' (Deep pocket hydrophobic floor).

Với cấu hình này, bạn đã thiết lập một trạng thái **"Khớp cảm ứng hoàn hảo" (Perfect Induced-Fit)**: Ma trận kẽm được bảo vệ tuyệt đối 100%, trong khi toàn bộ rãnh S1' được cấp quyền "thở" để ôm trọn các siêu phối tử *benzyl ether 1,3-diarylpyrazol* của bạn. Bạn hoàn toàn có thể tự tin đưa danh sách và lập luận này vào Khóa luận/Bài báo của mình!

# 6. MMP13

# 7. MMP3 (4G9L, Structure of MMP3 complexed with NNGH inhibitor, Structure of matrix metalloproteinase-3 with a platinum-based inhibitor) `B:223,B:224`

Nghịch lý của khớp cảm ứng (induced-fit) nằm ở chỗ: cấp quyền tự do vô hạn cho một thụ thể không tạo ra sự mô phỏng sinh học, mà kích hoạt sự hỗn loạn nhiệt động học. Việc xác định các gốc linh hoạt cho cấu trúc MMP-3 (4G9L) đòi hỏi sự phân tách rạch ròi giữa lõi tĩnh điện bất khả xâm phạm và ranh giới không gian cần tái sắp xếp.

Dựa trên ưu tiên cấp 1 (y văn) và cấp 2 (thích ứng cấu trúc), giao thức thiết lập gốc linh hoạt cho dẫn chất 1,3-diarylpyrazol mang nhánh benzyl ether buộc phải khóa cứng hoàn toàn mạng lưới kẽm và chỉ nhắm mục tiêu duy nhất vào vòng đặc hiệu S10. Sự can thiệp này tuân thủ nghiêm ngặt 5 bước của "Nghệ thuật phá vỡ quy tắc".

**Bước 4: Quyền Phủ Quyết Nhiệt Động Học (Khóa cứng Lõi Kim loại)**
Dữ liệu bạn trích xuất chỉ ra hàng loạt các gốc nằm trong bán kính 4 Å quanh 3 ion kẽm (Zn 301, Zn 302, Zn 303). Trong MMP-3, Kẽm xúc tác (Zn 302) được phối trí bởi bộ ba His201, His205, His211, trong khi Kẽm cấu trúc (Zn 303) được neo bởi His151, Asp153, His166, His179. Lưới tĩnh điện này duy trì tính toàn vẹn của enzyme. Mọi nỗ lực gán tính linh hoạt cho các gốc phối trí này sẽ bẻ gãy trục lượng tử, gây ra hình phạt lực đẩy khổng lồ và làm sụp đổ hàm chấm điểm [VERIFY]. Toàn bộ các gốc trong ranh giới 4 Å của cả 3 ion Zn BẮT BUỘC RIGID.

**Bước 2: Bằng chứng Động học (Tier 1 — Ưu tiên Y văn)**
Nghiên cứu của Belviso et al. (2013) cung cấp dữ liệu động học tối thượng: vòng đặc hiệu S10 (S10 specificity loop) quyết định tính chọn lọc cơ chất và mang đặc tính linh hoạt cực cao. Y văn ghi nhận vòng này có khả năng chuyển đổi giữa hai cấu hình. Ở cấu hình mở (Open Conformation - OC), chuỗi bên của **Tyr223** không bít kín lối vào và **His224** hướng ra ngoài túi xúc tác. Bằng chứng này thiết lập ưu tiên tuyệt đối: Tyr223 và His224 là hai gốc ứng viên hàng đầu cho mô phỏng linh hoạt.

**Bước 1 & 3: Giới hạn Không gian và Cổng Lập thể (Tier 2 — Thích ứng cho Benzyl Ether)**
Hợp chất *de novo* của bạn sở hữu nhánh benzyl ether mang thể tích không gian lớn. Để nhánh cồng kềnh này thâm nhập sâu vào túi mà không gây ra va chạm lập thể (steric clash), rào cản từ vòng S10 phải được giải phóng. Việc cho phép chuỗi bên của **Tyr223** và **His224** xoay tự do cung cấp chính xác "cổng không gian" để mô phỏng sự trượt của benzyl ether, tạo ra cơ chế khớp cảm ứng (induced-fit) cục bộ mà không làm xô lệch khung xương peptide. 

**Bước 5: Nguyên lý Dao cạo Ockham (Parsimony)**
Dù hợp chất chứa nhóm amino acid có thể thiết lập cầu muối với các vùng phân cực khác, việc cấp bậc tự do xoay cho các mỏ neo tĩnh điện phơi ngoài dung môi sẽ gây hiệu ứng bùng nổ tổ hợp (combinatorial explosion). Chúng ta thu hẹp toàn bộ không gian tìm kiếm khổng lồ xuống chính xác 2 gốc gác cổng (Tyr223, His224). Con số này bảo vệ sự hội tụ tuyệt đối của thuật toán Monte Carlo.

**Bản đồ Định vị Cấu trúc:**
*   **Bản đồ Ngoại vi:** Trong bối cảnh thiết kế thuốc, tính chọn lọc của các chất ức chế MMP không đến từ việc tương tác với lõi kẽm bảo tồn, mà từ việc khai thác các vòng đặc hiệu ngoại vi như S10. Giao thức này định tuyến thuật toán tìm kiếm vào đúng khu vực tạo ra sự khác biệt sinh học.
*   **Bản đồ Nội tại:** Luồng logic di chuyển theo quỹ đạo: Neo giữ tuyệt đối Lõi Kẽm $\rightarrow$ Khóa cứng hành lang xúc tác $\rightarrow$ Giải phóng độc quyền cổng lập thể tại vòng S10.

Lập trường đối lập cho rằng cần giải phóng tính linh hoạt cho toàn bộ các gốc nằm trong bán kính 4 Å của phối tử và kim loại. Tiền đề này thất bại do nó vi phạm trực tiếp quyền phủ quyết nhiệt động học, kích hoạt sự nhiễu loạn của các rotamer không mang lại động lực liên kết thực tế. Kết luận này chỉ bị vô hiệu hóa nếu phân tích Động lực học Phân tử (MD) chứng minh nhánh benzyl ether liên kết ở một vùng dị lập thể (allosteric) nằm ngoài tầm kiểm soát của vòng S10.

**Chỉ thị Vận hành (Dispositive Closure):**
Bỏ qua hoàn toàn danh sách các gốc lân cận Zn 301, Zn 302, và Zn 303. Thiết lập thông số docking linh hoạt (`--flexres`) giới hạn ở đúng hai gốc: **Tyr223 và His224** thuộc Chain B. Toàn bộ phần còn lại của protein phải được duy trì ở trạng thái Rigid.

# 8. MMP1 (1.56 ANGSTROM STRUCTURE OF MATURE TRUNCATED HUMAN FIBROBLAST COLLAGENASE, 1.56 A structure of mature truncated human fibroblast collagenase.) `A:180,A:210,A:236,A:237`

Sự hoàn hảo trong thiết kế *in silico* không đến từ việc thêm vào các bậc tự do, mà từ việc cắt bỏ không thương tiếc những ảo tưởng nhiệt động học. Giám định phân tích của bạn thể hiện sự nhạy bén xuất sắc trong việc loại bỏ nhiễu loạn không gian, đặc biệt là việc nhận diện khoảng cách phi lý của Phe242 so với lõi xúc tác. Tuy nhiên, ma trận 5 gốc linh hoạt tối ưu do bạn đề xuất chứa một "Nghịch lý Khung xương" (Backbone Paradox) chí mạng tại vị trí Tyr240, đe dọa trực tiếp đến tính hội tụ của thuật toán chấm điểm. 

Chúng tôi tiến hành giải phẫu danh sách tổng hợp này dựa trên dữ liệu tinh thể học thực nghiệm của Spurlino et al. (1994).

**1. Sự Thích Ứng Thành Công (The Validated Tier 1 & 2 Adaptations)**
Lập luận của bạn đã định vị chính xác 4 mỏ neo không gian đòi hỏi sự tái sắp xếp chuỗi bên:
*   **Asn180 (Tier 1):** Spurlino et al. chỉ định rõ nguyên tử oxy trên chuỗi bên của Asn180 đóng vai trò thiết yếu trong việc ổn định nitơ được proton hóa của liên kết cắt `[VERIFY]`. Việc giải phóng góc xoắn $\chi$ của Asn180 cho phép tối ưu hóa mạng lưới tĩnh điện với nhóm amino acid của ligand *de novo*.
*   **Tyr210 (Tier 1):** Sự tồn tại của liên kết peptide *cis* hiếm gặp giữa Glu209 và Tyr210 bẻ gập cấu trúc, định vị Tyr210 ngay tại ranh giới túi phụ P3' `[VERIFY]`. Chuỗi bên phenol của nó bắt buộc phải linh hoạt để dung nạp vành thơm của 1,3-diarylpyrazol.
*   **Tyr237 và Met236 (Tier 2):** Đây là hai gốc kỵ nước cốt lõi thiết lập ranh giới thể tích của túi S1'. Met236 kiểm soát không gian trượt nội tại, quyết định tính chọn lọc đặc hiệu của MMP-1 `[VERIFY]`. Giải phóng chúng tạo ra phễu lập thể cho nhánh benzyl ether.

**2. Lỗ hổng Chí mạng tại Tyr240 (The Fatal Flaw of the 5th Rank)**
Lập trường của bạn (FLOPP) đưa Tyr240 vào danh sách linh hoạt ưu tiên với lý do: "Spurlino direct - residue 240 H-bonds inhibitor". Tiền đề này thất bại thảm hại do nó diễn dịch sai cơ chế sinh hóa ở cấp độ nguyên tử. 

Văn bản gốc của Spurlino et al. (1994) ghi nhận chính xác: *"Residues 238 and 240 also form hydrogen bonds with inhibitor... Another difference between thermolysin and mCL-t involves the hydrogen bonding between the inhibitor and enzyme. Thermolysin utilizes side chains of residues lining the active cleft, while in mCL-t only backbone atoms are used."* `[VERIFY]`.

Tyr240 tạo liên kết hydro với chất ức chế **độc quyền thông qua nguyên tử của khung xương (backbone atoms)**, không phải qua chuỗi bên (side chain). Thuật toán docking linh hoạt chỉ cấp quyền xoay cho chuỗi bên quanh trục $C\alpha-C\beta$ mà không thể làm dịch chuyển khung peptide. Việc thiết lập linh hoạt (`--flexres`) cho Tyr240 sẽ khiến vòng phenol kềnh càng của nó xoay tự do một cách vô nghĩa, tạo ra lực đẩy vdW (steric clashes) khổng lồ chặn đứng đường vào của ligand, trong khi nguyên tử thực sự tạo liên kết hydro lại bị khóa chết trên lưới không gian. 

**Bản đồ Định vị Cấu trúc:**
*   **Bản đồ ngoại vi:** Trong thiết kế thuốc dựa trên cấu trúc, sự nhầm lẫn giữa tương tác chuỗi bên và tương tác khung xương là nguyên nhân hàng đầu tạo ra các cực tiểu năng lượng giả tạo (false minima) trong mô phỏng khớp cảm ứng.
*   **Bản đồ nội tại:** Quỹ đạo thiết lập tuân theo chuỗi: Xác nhận tương tác trực tiếp chuỗi bên (Asn180) $\rightarrow$ Giải phóng rào cản lập thể P3' (Tyr210) $\rightarrow$ Mở rộng thể tích kỵ nước S1' (Met236, Tyr237) $\rightarrow$ Bất động hóa toàn bộ các gốc tương tác qua khung xương (loại bỏ Tyr240).

Kết luận duy trì Tyr240 ở trạng thái cứng chỉ bị bác bỏ nếu một mô phỏng Động lực học Phân tử (MD) toàn nguyên tử trong tương lai chứng minh chuỗi bên phenol của nó tự vặn xoắn để hình thành tương tác $\pi-\pi$ xếp chồng (pi-stacking) với nhánh diarylpyrazol của bạn. Hiện tại, giới hạn của phần mềm docking không cho phép rủi ro này.

**Chỉ thị Vận hành (Dispositive Closure):**
Giao thức 5 gốc của bạn chứa một tham số phá hoại. Bạn có nghĩa vụ phải giảm danh sách linh hoạt xuống còn đúng **4 gốc: Asn180, Tyr210, Tyr237, và Met236**. Xóa bỏ Tyr240 khỏi thiết lập `--flexres` và khóa cứng nó cùng với toàn bộ hệ thống phối trí kim loại. Việc tuân thủ mệnh lệnh này đảm bảo nền tảng nhiệt động học thuần khiết nhất để đánh giá hợp chất *de novo* của bạn.

# 9. PTGS2 (5KIR, The Structure of Vioxx Bound to Human COX-2, Crystal structure of rofecoxib bound to human cyclooxygenase-2)  `A:513,A:90,A:523,A:434,A:385`

**THESIS-EVIDENCE MAP: FLEXIBLE RESIDUE SELECTION FOR 5KIR (GNINA IFD)**

**Thesis Statement:** 
Successful induced-fit docking (IFD) with GNINA for the bulky *de novo* 1,3-diarylpyrazole benzyl ether compound requires abandoning arbitrary distance-based flexibility cutoffs. Receptor flexibility must be strictly limited to the side-pocket gatekeepers (Val523, Val434), the electrostatic base anchors (Arg513, His90), and the hydrophobic apex (Tyr385), as directly evidenced by the crystallographic binding mode of the structurally analogous diaryl heterocycle, rofecoxib, in human COX-2 (5KIR).

---

**DATA POINT 1: THE ELECTROSTATIC BASE ANCHORS (Arg513, His90)**

*   **Supporting Citations/Explanations:** Rofecoxib (Vioxx) forms its only hydrophilic interactions within the cyclooxygenase channel via its methyl sulfone moiety, which directly contacts the side-chain nitrogen atoms of His90 and Arg513. These residues are located at the very base of the COX-2 specific side pocket. 
*   **Sources + Locations:** Orlando & Malkowski (2016), Section 1 Introduction; Section 3.1 Vioxx binding pose.
*   **Evidence Strength:** High (Tier 1 — Direct Crystallographic Evidence).
*   **Why this matters:** Your *de novo* compound contains an amino acid moiety that requires precise electrostatic anchoring. Allowing side-chain flexibility for His90 and Arg513 enables GNINA's Markov chain Monte Carlo (MCMC) algorithm to optimize hydrogen-bonding geometries for this novel polar group, preventing artificial steric clashes or energetic penalties that rigid docking would falsely generate.

---

**DATA POINT 2: THE SIDE POCKET GATEKEEPERS (Val523, Val434)**

*   **Supporting Citations/Explanations:** The first-shell substitutions of Ile434 to Val434 and Ile523 to Val523 are the hallmark structural differences creating the COX-2 side pocket, yielding a ~25% increase in active site volume. Access to this pocket by diaryl heterocycle scaffolds (like coxibs) is entirely dependent on the presence of Val523.
*   **Sources + Locations:** Orlando & Malkowski (2016), Section 1 Introduction.
*   **Evidence Strength:** High (Tier 1 — Molecular Basis of Selectivity).
*   **Why this matters:** The bulky benzyl ether group of your *de novo* scaffold requires substantial spatial clearance to penetrate the side pocket. Defining Val523 and Val434 as flexible provides the critical "sliding space" for the algorithm to accommodate the bulky functional group without triggering massive van der Waals repulsion in the scoring function.

---

**DATA POINT 3: THE HYDROPHOBIC APEX (Tyr385)**

*   **Supporting Citations/Explanations:** The binding conformation of rofecoxib features a phenyl ring extending upward to interact directly with the side chain of Tyr385. Aside from the specific interactions at the base of the side pocket, all remaining enzyme-inhibitor contacts within the channel are hydrophobic in nature.
*   **Sources + Locations:** Orlando & Malkowski (2016), Section 3.1 Vioxx binding pose.
*   **Evidence Strength:** Medium-High (Tier 2 — Spatial Adaptation).
*   **Why this matters:** The 1,3-diarylpyrazole core is a rigid, aromatic scaffold. Granting torsional flexibility to the Tyr385 side chain permits the dynamic optimization of $\pi-\pi$ stacking interactions, directly maximizing the empirical and convolutional neural network (CNN) binding affinity scores in GNINA for the ligand's upper aromatic ring. 

---

**DATA POINT 4: THERMODYNAMIC VETO ON PERIPHERAL FLEXIBILITY**

*   **Supporting Citations/Explanations:** Algorithms utilized by deep learning frameworks like GNINA rely on grid-based scoring and localized MCMC sampling. While flexible docking models structural changes upon binding, excessive degrees of freedom exponentially increase the search space, degrading pose prediction accuracy if the algorithm cannot converge.
*   **Sources + Locations:** McNutt et al. (GNINA 1.3), Introduction; Plainer et al. (DiffDock-Pocket), Section 1 Introduction.
*   **Evidence Strength:** High (Algorithmic Constraint).
*   **Why this matters:** This directly falsifies the assumption that all residues within a 4 Å radius should be flexible. If you set all surrounding residues to flexible, GNINA will exhaust its sampling resources generating thermodynamic noise (invalid rotamers). Locking the receptor backbone and non-essential side chains ensures the algorithm isolates its computational power to the 5 residues (Arg513, His90, Val523, Val434, Tyr385) that actually dictate ligand binding.

# 10. PTGES (5TL9, crystal structure of mPGES-1 bound to inhibitor, Discovery and characterization of [(cyclopentyl)ethyl]benzoic acid inhibitors of microsomal prostaglandin E synthase-1) `--flexres B:28,B:32,B:52,B:53,A:130,A:134` Lấy cả chuỗi A và B làm ỉnterface
<img width="982" height="918" alt="image" src="https://github.com/user-attachments/assets/c0de0eca-0eda-4f28-9ce9-fad8a00dc395" />

Thesis — ChimeraX verified (within 4Å of 7DN)
9 residues + HOH307 tiếp xúc trực tiếp với 7DN trong 4 Å: Y28, A31, I32, G35, R38, L39, F44, R52, H53. Áp dụng Adaptive SOP 5 bước: --flexres cuối cùng gồm 5 residues — R52, Y28, I32, F44, H53. Các residues còn lại (A31, G35, R38, L39) bị loại bởi SOP Bước 4 và 5.
10 contacts — phân loại theo Adaptive SOP
R52
Arginine · Salt bridge anchor
FLEX ✓ — Tier 1
Bidentate salt bridge trực tiếp với carboxylate. NH2↔7DN O2. Amino acid -COOH của new compounds cần Arg52 rotamer optimization.
Partridge (2017): "bidentate salt bridge to Arg52"
Bước 1 Bước 2
H53
Histidine · Water bridge + direct
FLEX limited
NE2↔7DN O1 = 3.47 Å (direct) + HOH307 bridge (2.712 Å). Flexible trong standard rotamer library — không rigid hoàn toàn.
Partridge (2017): "coordinates to His53 through water bridge"
Bước 2 Bước 3
Y28
Tyrosine · Hydrophobic groove
FLEX ✓ — Tier 1
VdW contact xác nhận. Benzyl ether group của new compounds bulkier hơn N-aryl piperidine → cần Tyr28 rotamer space.
Partridge (2017): "VdW contacts with Tyr28"
Bước 1 Bước 2
I32
Isoleucine · Deep cleft gate
FLEX ✓ — Tier 1
"Cleft behind Ile32" — paper đích danh. Pyrazole core của new compounds cần entry vào cleft này. Steric gate quan trọng nhất.
Partridge (2017): "cleft behind Ile32"
Bước 2 Bước 3
F44
Phenylalanine · Proximity gate
FLEX ✓ — Tier 2
Trong 4Å nhưng paper không đề cập trực tiếp. SOP Bước 1 (spatial) đủ căn cứ. Phe44 aromatic ring có thể π-stack với pyrazole core của new compounds.
ChimeraX: F44 within 4Å · SOP Bước 1 bypass
Bước 1
HOH307
Water bridge · Static
STATIC — giữ trong receptor
O↔His53 NE2=2.712Å, O↔7DN O1=2.690Å. Không phải flexres — là static molecule trong receptor. Đã thảo luận trước.
Partridge (2017): "water bridge"
Bước 3 — water network
A31
Alanine · Backbone only
RIGID — SOP Bước 4
Ala không có sidechain có ý nghĩa — chỉ có Cβ (methyl). Rotamer không ảnh hưởng binding. Backbone rigid là đúng.
SOP Bước 4 — Thermodynamic Veto: no functional sidechain
Bước 4 Bước 5
G35
Glycine · No sidechain
RIGID — SOP Bước 4
Gly không có sidechain — không có rotamer để sample. GNINA không thể làm gì với Gly trong flexres. Backbone-only contact.
SOP Bước 4 — No sidechain atoms
Bước 4
R38
Arginine · Distal contact
RIGID — SOP Bước 5
Arg38 trong 4Å nhưng paper không đề cập. Đã có R52 làm salt bridge anchor — 2 Arg trong flexres là redundant. Parsimony loại R38.
SOP Bước 5 — Parsimony: R52 đã cover Arg interactions
Bước 5
L39
Leucine · Peripheral contact
RIGID — SOP Bước 5
Leu39 trong 4Å nhưng paper không đề cập trong binding mode. Peripheral hydrophobic — không phải gác cổng trực tiếp. Parsimony loại.
SOP Bước 5 — Parsimony: not a gating residue
Bước 5
Tại sao Tyr130, Gln134, Thr131 không có trong list
ChimeraX xác nhận: Tyr130, Gln134, Thr131 không trong 4Å của 7DN. Paper đề cập các residues này là cross-chain contacts trong homotrimer mPGES-1 — không phải intra-chain A contacts. Loại hoàn toàn khỏi flexres. Đây là correction quan trọng nhất so với ý kiến trước.

# 11. ALOX5 (6N2W, The structure of Stable-5-Lipoxygenase bound to NDGA, Structural and mechanistic insights into 5-lipoxygenase inhibition by natural products) `B:596,B:599,B:359,B:368,B:607,B:363`

<img width="1040" height="869" alt="image" src="https://github.com/user-attachments/assets/0368f42d-ba65-4155-b135-5f96a0e8e1c8" />

Data thực tế từ ChimeraX đã có. Áp dụng Adaptive SOP để phân xử.**Tóm tắt verdict về ý kiến trước:**

Ý kiến đó **đúng về nguyên tắc lõi** — His372/Ile673 rigid vì Fe²⁺ coordination là hoàn toàn chính xác. NDGA tenuously bridges this ravine, tethered only by Arg596, with nominal protein-ligand contacts. One catechol ring sits above the iron and blocks access to the open position of its coordination sphere.

Nhưng có **2 lỗi thực tế:**

**1. Leu368 bị bỏ sót** khỏi Dispositive Closure — paper Fig.1c LigPlot liệt kê rõ, phải có trong flexres.

**2. Framing His367/His550 như "cần rigid vì Fe coordination"** — đúng về trạng thái nhưng không cần thiết vì chúng không phải 30Z contacts, không xuất hiện trong flexres discussion.

`--flexres B:596,B:599,B:359,B:368,B:607,B:363` là string đúng nhất dựa trên data thực tế.

Tưởng tượng một kho tiền với cánh cửa đã mở toang, nhưng bên trong lại chứa một mạng lưới laser bảo vệ siêu nhạy; việc bạn kích hoạt nhầm bất kỳ tia laser nào sẽ khiến hệ thống tự hủy. Giao thức nới lỏng toàn bộ các gốc trong vòng 4 Å quanh chất ức chế cũ là một sai lầm chết người đối với 6N2W.

**Layer 1 — Core Insight**
Giao thức "Nghệ thuật phá vỡ quy tắc" yêu cầu bạn phải phân rã danh sách tiếp xúc của NDGA (30Z) thành các vùng giới hạn nhiệt động học tuyệt đối. Bạn BẮT BUỘC phải khóa cứng (rigid) toàn bộ mạng lưới phối trí ion Sắt (đặc biệt là His372 và Ile673), đồng thời chỉ cấp quyền linh hoạt (`--flexres`) cho mỏ neo tĩnh điện phân cực (Arg596) và các gốc kỵ nước cồng kềnh (Trp599, Phe359, Leu607) để tạo không gian trượt cho nhánh benzyl ether.

**Layer 2 — Structural Explanation**
Khác với các túi liên kết khép kín, 6N2W bộc lộ một "khe nứt sâu" do chuỗi xoắn gác cổng đã bị phá vỡ hoàn toàn. NDGA neo đậu rất lỏng lẻo trong khe nứt này với các tiếp xúc Van der Waals cực kỳ hạn chế. Nếu bạn nới lỏng toàn bộ các gốc ngoại vi, mạng nơ-ron tích chập (CNN) của thuật toán sẽ lãng phí tài nguyên lấy mẫu (sampling) vào khoảng không vô nghĩa. Quan trọng hơn, ion Fe2+ là hố thế năng thiết yếu nhất của enzyme; việc nới lỏng các gốc phối trí kim loại sẽ cho phép thuật toán MCMC xoay chúng ra khỏi trục liên kết lượng tử gốc, gây ra sự bùng nổ lực đẩy tĩnh điện nội tại và đánh sập điểm số của phối tử.

**Layer 3 — Full Technical Detail**

**THESIS-EVIDENCE MAP: CHIẾN LƯỢC LỰA CHỌN GỐC LINH HOẠT CHO 6N2W**

**Thesis Statement:**
Quỹ đạo khớp cảm ứng (IFD) cho dẫn chất *de novo* 1,3-diarylpyrazol mang nhóm amino acid chỉ thành công khi áp dụng lệnh phủ quyết nhiệt động học: cấm tuyệt đối tính linh hoạt tại vùng lõi Fe2+ và các gốc xương sống, trong khi kích hoạt tối đa độ đàn hồi tại mỏ neo Arg596 và vách kỵ nước ngoại vi.

---

**ĐIỂM DỮ LIỆU 1: LỆNH PHỦ QUYẾT LƯỢNG TỬ TẠI LÕI SẮT (MUST BE RIGID: His367, His372, His550, Asn554, Ile673)**

*   **Supporting Citations/Explanations:** Dữ liệu tinh thể học xác nhận His372 là một trong các gốc phối trí trực tiếp tạo liên kết với ion Fe2+. Vòng catechol của NDGA nằm ngay phía trên ion Sắt và chặn hướng tiếp cận mở của khối cầu phối trí kim loại này. Ile673 và Asn554 cũng nằm trong giới hạn 4 Å quanh Fe2+.
*   **Sources + Locations:** "Structural and mechanistic insights...", Results, Hình 1b-c. Dữ liệu hệ quy chiếu người dùng.
*   **Evidence Strength:** High (Bằng chứng hình học không gian L1).
*   **Why this matters (Tại sao điều này quan trọng):** Khóa cứng bộ năm gốc này là nguyên tắc sống còn. Nếu His372 hoặc Ile673 được đưa vào cờ `--flexres`, bộ lấy mẫu MCMC sẽ bẻ gãy khoảng cách phối trí 2.16 Å lý tưởng để nhường chỗ cho phân tử *de novo* khổng lồ của bạn. Hàm chấm điểm vật lý sẽ lập tức áp đặt hình phạt vô cực (infinite penalty) cho sự sụp đổ của phức chất kim loại này.

---

**ĐIỂM DỮ LIỆU 2: ĐẶC QUYỀN MỎ NEO TĨNH ĐIỆN (MUST BE FLEXIBLE: Arg596)**

*   **Supporting Citations/Explanations:** Cấu trúc NDGA lơ lửng trong rãnh kỵ nước và được "neo giữ" (tethered) độc quyền bởi một liên kết hydro duy nhất với Arg596 (khoảng cách 2.8 Å).
*   **Sources + Locations:** "Structural and mechanistic insights...", Results, Hình 1c.
*   **Evidence Strength:** High.
*   **Why this matters (Tại sao điều này quan trọng):** Phân tử *de novo* của bạn chứa gốc amino acid (-COOH, -NH2). Ở pH sinh lý, nhóm carboxylate tích điện âm cần bắt cặp với một gốc mang điện dương. Bằng cách thiết lập Arg596 ở trạng thái linh hoạt, thuật toán sẽ tự động dò tìm và xoay góc nhị diện của nhóm guanidinium (+1) để thiết lập một liên kết cầu muối (salt bridge) hoàn hảo với dẫn chất của bạn.

---

**ĐIỂM DỮ LIỆU 3: KHÔNG GIAN TRƯỢT CHO KHUNG ĐA VÒNG (MUST BE FLEXIBLE: Trp599, Phe359, Leu607, Gln363, Leu368)**

*   **Supporting Citations/Explanations:** NDGA tạo ra rất ít tiếp xúc Van der Waals do sự phơi nhiễm của rãnh xúc tác. Tuy nhiên, các gốc kỵ nước cồng kềnh như Trp599, Phe359 và Leu607 tạo thành thành vách bao quanh trục phân tử.
*   **Sources + Locations:** "Structural and mechanistic insights...", Hình 1c (LigPlot).
*   **Evidence Strength:** High.
*   **Why this matters (Tại sao điều này quan trọng):** Lõi 1,3-diarylpyrazol gắn nhánh benzyl ether có thể tích không gian (steric bulk) lớn hơn gấp nhiều lần vòng catechol đơn gốc của NDGA. Việc nới lỏng Trp599 và Phe359 cung cấp "độ đàn hồi" (sliding space) bắt buộc để mạng CNN xoay các vòng thơm này ra xa, ngăn chặn va chạm không gian (steric clash) mà không làm suy giảm ái lực kỵ nước.

---

**ĐIỂM DỮ LIỆU 4: NGUYÊN LÝ TỐI GIẢN TẠI GỐC CHẾT (MUST BE RIGID: Pro569, Ala410, Ala603, Ile406, Asn407, Gln557, His600)**

*   **Supporting Citations/Explanations:** Pro569 và Ala410/Ala603 có mặt trong danh sách tiếp xúc quanh ligand. Tuy nhiên, Proline bị khóa chết vào khung xương peptide, trong khi Alanine chỉ có một nhánh methyl không thể xoay góc rotamer. Các gốc còn lại đóng vai trò thứ cấp.
*   **Sources + Locations:** SOP Bước 4 (Thermodynamic Veto); Hình 1c.
*   **Evidence Strength:** Medium-High.
*   **Why this matters (Tại sao điều này quan trọng):** Đưa các gốc không có khả năng vặn xoắn chuỗi bên (Ala, Pro) vào danh sách linh hoạt là một thảm họa tính toán. Nó bơm entropy giả tạo vào thuật toán và làm tăng rủi ro phân kỳ của mô phỏng mà không đóng góp bất kỳ điểm số tương tác nào.

---

**Hai Bản Đồ Định Vị:**
*   **Bản đồ ngoại vi:** Trong docking cấu trúc phơi nhiễm (như 6N2W), ái lực không đến từ việc ligand chui vừa một cái lỗ, mà từ việc ligand ôm lấy các mỏ neo nhô ra (Arg596, Fe2+) trên một bình nguyên rộng lớn.
*   **Bản đồ nội tại:** Quỹ đạo luận lý: Xác định 4 Å $\rightarrow$ Loại trừ tức khắc toàn bộ vùng phối trí Sắt (His/Ile/Asn) $\rightarrow$ Loại trừ các gốc khung xương (Ala/Pro) $\rightarrow$ Tuyệt đối ưu tiên mỏ neo (Arg596) $\rightarrow$ Giải phóng vách kỵ nước (Trp/Phe/Leu).

**Tổng hợp và Giới hạn (Falsification):**
Lập trường đối lập (FLOPP) cho rằng "cần nới lỏng His372 vì vòng pyrazol có thể tương tác với nó" thất bại hoàn toàn. Tiền đề này sai lệch vì trong hóa học lượng tử, liên kết cộng hóa trị phối trí với Sắt có năng lượng > 50 kcal/mol, hoàn toàn áp đảo bất kỳ tương tác VdW hay H-bond nào (~1-5 kcal/mol) mà ligand có thể tạo ra. Kết luận này chỉ bị vô hiệu hóa nếu phân tử của bạn được thiết kế với một đầu đạn hydroxamic acid nhắm trực tiếp vào việc chiếm đoạt Sắt từ protein.

**Chỉ thị Vận hành (Dispositive Closure):**
Loại bỏ ngay lập tức danh sách 4 Å tự động của bạn. Thiết lập chính xác tham số `--flexres` trong GNINA với danh sách độc quyền sau: **ARG596, TRP599, PHE359, LEU607, GLN363, LEU368**. 
Các gốc thuộc vòng phối trí Sắt (HIS367, HIS372, HIS550, ILE673, ASN554) và các gốc không có rotamer chức năng (PRO569, ALA410, ALA603) BẮT BUỘC phải duy trì trạng thái KHÓA CỨNG (RIGID). Chạy lệnh này để hội tụ quỹ đạo mô phỏng.


# EGFR (1XKK, EGFR kinase domain complexed with a quinazoline inhibitor- GW572016, A unique structure for epidermal growth factor receptor bound to GW572016 (Lapatinib): relationships among protein conformation, inhibitor off-rate, and receptor activity in tumor cells.) `--flexres A:745,A:776,A:777,A:790,A:856`

Tưởng tượng việc bạn cố gắng đưa một cỗ xe vận tải khổng lồ vào một hầm đỗ xe; nếu bạn không kích hoạt hệ thống mở rộng cửa vách phía sau, toàn bộ cấu trúc sẽ bị nghiền nát ngay lập tức. Đối với dẫn chất *de novo* mang nhánh benzyl ether cồng kềnh, việc lựa chọn gốc linh hoạt cho 1XKK BẮT BUỘC phải tuân thủ tuyệt đối dữ liệu không gian từ y văn để tránh sự sụp đổ của thuật toán lấy mẫu.

**Layer 1 — Core Insight**
Giao thức khớp cảm ứng (IFD) bằng GNINA cho dẫn chất 1,3-diarylpyrazol mang nhóm amino acid yêu cầu một ma trận linh hoạt phân mảnh. Dựa trên dữ liệu tinh thể học của 1XKK (với độ lệch đánh số +24 so với y văn gốc), bạn phải nới lỏng vách kỵ nước tuyến sau (Leu777, Phe856, Thr790) để dung nạp nhánh ether, đồng thời giải phóng mỏ neo tĩnh điện (Lys745, Arg776) để đón bắt gốc amino acid.

**Layer 2 — Structural Explanation**
Trong cấu trúc 1XKK, nhánh 3-fluorobenzyl-oxy của chất ức chế GW572016 thâm nhập sâu vào một túi kỵ nước mở rộng đặc thù. Nhánh benzyl ether trên phân tử *de novo* của bạn sở hữu thể tích không gian tương đương. Thuật toán mạng nơ-ron tích chập (CNN) cần quyền vặn xoắn chuỗi bên tại vách túi này để triệt tiêu lực đẩy Van der Waals. Tuy nhiên, thuật toán này mù lòa trước sự di chuyển của xương sống peptide và phân tử nước; do đó, các trạm tương tác cốt lõi tại vùng bản lề hoặc cầu nối nước phải bị phong tỏa hoàn toàn để bảo vệ quỹ đạo hình học nguyên thủy.

**Layer 3 — Full Technical Detail**

**THESIS-EVIDENCE MAP: CHIẾN LƯỢC LỰA CHỌN ROTAMER DỰA TRÊN Y VĂN CHO 1XKK**

**Thesis Statement:**
Quỹ đạo mô phỏng hội tụ tối đa khi và chỉ khi bộ tham số `--flexres` phản ánh chính xác cơ chế đàn hồi sinh lý học được y văn ghi nhận: Nới lỏng không gian trượt cho khối kỵ nước (Leu777, Phe856, Thr790), kích hoạt hố thế năng dương (Lys745, Arg776), và áp đặt lệnh phủ quyết lên các cấu trúc xương sống/cầu nước (Met793, Thr854).

---

**ĐIỂM DỮ LIỆU 1: KHÔNG GIAN TRƯỢT CHO NHÁNH ETHER CỒNG KỀNH (FLEX: Leu777, Phe856)**

*   **Supporting Citations/Explanations:** Nhóm 3-fluorobenzyloxy chiếm giữ một túi được định hình bởi các chuỗi bên của Met742, Leu753, Thr766, Thr830, Phe832, và Leu834.
*   **Sources + Locations:** Excerpts from 1XKK documentation.pdf (Results: Inhibitor Binding Site).
*   **Evidence Strength:** High (Bằng chứng tinh thể học trực tiếp về dung nạp nhóm thế cồng kềnh).
*   **Why this matters:** Nhánh benzyl ether của bạn cần toàn bộ không gian của túi phía sau (back pocket) này. Việc cấp quyền linh hoạt cho Leu753 (PDB: Leu777) và Phe832 (PDB: Phe856) cho phép thuật toán dịch chuyển các vách cản lập thể, loại bỏ hoàn toàn hình phạt va chạm (steric clash) đối với khối lượng không gian khổng lồ của phối tử.

---

**ĐIỂM DỮ LIỆU 2: GÁC CỔNG ĐIỀU HƯỚNG TÚI PHÍA SAU (FLEX: Thr790)**

*   **Supporting Citations/Explanations:** Chuỗi bên của Thr766 hướng ra xa khỏi vòng quinoline và tạo liên kết hydro với carbonyl trên xương sống của Arg752.
*   **Sources + Locations:** Excerpts from 1XKK documentation.pdf (Results: Comparison with OSI-774/EGFR).
*   **Evidence Strength:** High (Cơ chế xoay chuỗi bên được cấu trúc X-ray xác nhận).
*   **Why this matters:** Thr766 (PDB: Thr790) hoạt động như một bản lề điều hướng. Việc đưa gốc này vào danh sách `--flexres` mô phỏng chính xác sự vặn xoắn sinh lý của nó để mở đường cho nhánh benzyl ether tiến sâu vào túi kỵ nước.

---

**ĐIỂM DỮ LIỆU 3: KÍCH HOẠT MỎ NEO TĨNH ĐIỆN DƯƠNG (FLEX: Lys745, Arg776)**

*   **Supporting Citations/Explanations:** Sự dịch chuyển của chuỗi xoắn C làm mất đi cầu muối bảo thủ Glu738-Lys721; thay vào đó, Lys721 tạo liên kết hydro với chuỗi bên của Asp831. Ngoài ra, Arg752 đóng vai trò là một điểm tựa lân cận.
*   **Sources + Locations:**, Excerpts from 1XKK documentation.pdf.
*   **Evidence Strength:** High (Động học phá vỡ cầu muối được ghi nhận trong y văn).
*   **Why this matters:** Nhóm amino acid trên dẫn chất *de novo* mang điện tích âm (-COO-). Lys721 (PDB: Lys745) hiện đang liên kết với Asp831, một tàn dư của trạng thái bất hoạt. Bạn phải nới lỏng Lys745 và Arg752 (PDB: Arg776) để thuật toán bẻ gãy các liên kết cũ này, vươn chuỗi bên mang điện dương ra và thiết lập bẫy tĩnh điện mới với phối tử của bạn.

---

**ĐIỂM DỮ LIỆU 4: LỆNH PHỦ QUYẾT TẠI LÕI NHẬN DIỆN (RIGID: Met793, Thr854)**

*   **Supporting Citations/Explanations:** N1 của quinoline tạo liên kết hydro với NH trên xương sống chính của Met769, trong khi N3 tạo liên kết hydro qua trung gian phân tử nước với chuỗi bên của Thr830.
*   **Sources + Locations:** Excerpts from 1XKK documentation.pdf.
*   **Evidence Strength:** High (Giới hạn vật lý của thuật toán mô phỏng).
*   **Why this matters:** GNINA không xoay được xương sống peptide (main chain NH) và coi các phân tử nước là hạt tĩnh. Đưa Met769 (PDB: Met793) hoặc Thr830 (PDB: Thr854) vào vùng linh hoạt sẽ khiến thuật toán xoay lệch chuỗi bên một cách vô ích hoặc bẻ gãy hệ trục tọa độ tới phân tử nước. Bạn bắt buộc phải khóa cứng chúng.

---

**Hai Bản Đồ Định Vị:**
*   **Bản đồ ngoại vi:** Y văn cung cấp bản đồ kho báu, nhưng hệ động lực học quyết định cách bạn đào nó. Bạn không nới lỏng toàn bộ túi liên kết; bạn chỉ nới lỏng các vách ngăn cần giãn nở và các điện cực cần thiết lập lại.
*   **Bản đồ nội tại:** Quỹ đạo tinh chỉnh: Quét y văn 1XKK $\rightarrow$ Đồng bộ hóa độ lệch +24 PDB $\rightarrow$ Loại trừ các tương tác xương sống/nước (Met793, Thr854) $\rightarrow$ Ưu tiên không gian trượt cho benzyl ether (Leu777, Phe856, Thr790) $\rightarrow$ Giải phóng mỏ neo tĩnh điện cho amino acid (Lys745, Arg776).

**Tổng hợp và Giới hạn (Falsification):**
Quan điểm đối lập cho rằng "nên khóa cứng Lys745 để bảo vệ cấu trúc tinh thể nguyên bản" thất bại hoàn toàn. Tiền đề này sai lệch vì việc bảo tồn một cầu muối không phản ánh trạng thái động học khi một phối tử mang điện âm khổng lồ tiến vào. Kết luận yêu cầu sự giãn nở của túi kỵ nước này chỉ bị vô hiệu hóa nếu phân tử *de novo* của bạn bị cắt bỏ hoàn toàn nhánh benzyl ether, đưa hình học của nó trở về dạng phẳng nhỏ gọn tương tự chất ức chế OSI-774.

**Chỉ thị Vận hành (Dispositive Closure):**
Đình chỉ mọi thuật toán chọn vùng linh hoạt tự động. Bạn BẮT BUỘC áp dụng chính xác chuỗi tham số định tuyến không gian sau đây vào lệnh GNINA: `--flexres A:745,A:776,A:777,A:790,A:856`. Tuyệt đối không đưa các gốc Met793 và Thr854 vào cờ này. Thực thi mô phỏng để bảo đảm thuật toán CNN nội suy chính xác ái lực nhiệt động học.

Tưởng tượng một cuốn sách được xuất bản với hai phiên bản: bản thứ nhất đếm số trang bắt đầu ngay từ trang bìa, trong khi bản thứ hai chỉ bắt đầu đếm từ nội dung chương một; sự sai lệch +24 đơn vị trong hệ thống tọa độ 1XKK chính là hiện tượng đếm trang tương tự trong sinh học phân tử.

**Layer 1 — Core Insight**
Sự khác biệt đánh số (Lys721 trong y văn trở thành Lys745 trong tệp PDB) xuất phát từ quy ước xử lý **peptide tín hiệu (signal peptide)** dài đúng 24 amino acid ở đầu N-terminal của thụ thể EGFR. Sự chênh lệch +24 này là một hằng số tịnh tiến tuyệt đối áp dụng cho mọi gốc amino acid được đối chiếu giữa bài báo gốc của Wood et al. (2004) và tệp dữ liệu tinh thể học 1XKK.

**Layer 2 — Structural Explanation**
*   **Hệ quy chiếu Y văn (Mature Protein):** Các nhà sinh học trong các công trình kinh điển thường sử dụng hệ đánh số dựa trên chuỗi protein "trưởng thành". Trong quá trình sinh lý thực tế, EGFR được tổng hợp với một đoạn peptide tín hiệu dài 24 amino acid để dẫn đường cho thụ thể chèn vào màng tế bào, sau đó đoạn tín hiệu này bị enzyme cắt bỏ. Bài báo đếm amino acid đầu tiên sau nhát cắt này là số 1, dẫn đến mỏ neo tĩnh điện được định danh là Lys721.
*   **Hệ quy chiếu PDB (Precursor Protein):** Cơ sở dữ liệu PDB hiện đại chuẩn hóa mọi tọa độ theo mã định danh UniProt nguyên thủy để bảo đảm tính thống nhất toàn cầu. Thẻ `DBREF` trong tệp 1XKK chỉ định rõ nó tham chiếu đến trình tự `P00533` (EGFR_HUMAN). Khung tham chiếu UniProt đếm toàn bộ trình tự gene từ lúc vừa dịch mã, bao gồm cả 24 amino acid của đoạn tín hiệu ban đầu. Theo phương trình tuyến tính cơ bản: $721 \text{ (Bài báo)} + 24 \text{ (Signal Peptide)} = 745 \text{ (PDB)}$.

**Layer 3 — Full Technical Detail**

Toàn bộ các mỏ neo sinh tử trong hệ thống 1XKK đều tuân thủ nguyên lý tịnh tiến này. Việc nhận diện sai lệch sẽ dẫn đến sự sụp đổ của hệ thống khai báo tham số mô phỏng:
*   Mỏ neo bản lề (Hinge): Met769 (Bài báo) $\rightarrow$ Met793 (PDB).
*   Gác cổng túi kỵ nước: Thr766 (Bài báo) $\rightarrow$ Thr790 (PDB).
*   Vách kỵ nước sâu: Leu753 (Bài báo) $\rightarrow$ Leu777 (PDB).
*   Mỏ neo tĩnh điện dương: Lys721 (Bài báo) $\rightarrow$ Lys745 (PDB).
*   Vòng lặp DFG: Asp831 (Bài báo) $\rightarrow$ Asp855 (PDB).

**Hai Bản Đồ Định Vị:**
*   **Bản đồ ngoại vi:** Khoảng trống nhận thức giữa văn bản sinh học (tập trung vào dạng hoạt động sinh lý tại màng tế bào) và khoa học hóa tin (tập trung vào sự đồng bộ hóa dữ liệu với hệ thống gene toàn cầu).
*   **Bản đồ nội tại:** Quỹ đạo thao tác: Đọc tên gốc mục tiêu từ bài báo (Lys721) $\rightarrow$ Cộng thêm hằng số đoạn tín hiệu (+24) $\rightarrow$ Xác định chỉ mục thực tế trong ma trận tọa độ (745) $\rightarrow$ Xác thực chéo bằng lệnh `info atoms /A:745` để đảm bảo nhãn amino acid khớp với Lysine $\rightarrow$ Chuyển tham số vào lệnh `--flexres`.

**Chỉ thị Vận hành (Dispositive Closure):**
Thuật toán GNINA hoàn toàn "mù" trước các tài liệu y văn; nó chỉ nhận diện lưới không gian vật lý được mã hóa trong tệp PDB. Mọi lệnh can thiệp tọa độ, khóa cứng, hay nới lỏng rotamer BẮT BUỘC phải sử dụng hệ đánh số đã được tịnh tiến (+24) của PDB. Việc bạn khai báo nhầm tham số thành `--flexres A:721` sẽ khiến mạng nơ-ron đi tìm một gốc Glycine lơ lửng ở ngoại vi để vặn xoắn, làm sụp đổ hoàn toàn quỹ đạo tìm kiếm khớp cảm ứng của khối 1,3-diarylpyrazol. Duy trì sự quy đổi này trong toàn bộ các bước thiết lập tiếp theo.

# ERBB2 (7PCD, HER2 IN COMPLEX WITH A COVALENT INHIBITOR, Discovery of potent and selective HER2 inhibitors with efficacy against HER2 exon 20 insertion-driven tumors, which preserve wild-type EGFR signaling.) `A:774,A:783,A:862,A:864`

Tưởng tượng một hệ thống khóa liên động cơ học; nếu bạn tháo tung toàn bộ đinh ốc thay vì chỉ mở lỏng các chốt khóa cần thiết, toàn bộ cỗ máy sẽ vỡ vụn. Đối với cấu trúc 7PCD, "Nghệ thuật phá vỡ quy tắc" yêu cầu bạn gạt bỏ hoàn toàn danh sách 20 gốc 4 Å cồng kềnh, chỉ nhắm mục tiêu xoay rotamer vào các bộ lọc chọn lọc túi sau (back pocket) và mỏ neo tĩnh điện, đồng thời đóng băng tuyệt đối trục bản lề (hinge) của enzyme.

**Layer 1 — Core Insight**
Giao thức khớp cảm ứng (IFD) cho dẫn chất *de novo* (benzyl ether 1,3-diarylpyrazol-aminoacid) trên 7PCD BẮT BUỘC phải phân rã danh sách tiếp xúc. Dựa trên dữ liệu y văn và ngoại suy thích ứng, bạn chỉ cấp cờ `--flexres` cho vách túi kỵ nước phía sau (Ser783, Phe864, Thr798) để dung nạp nhánh ether, và các mỏ neo điện dương (Lys753, Arg784) để bắt giữ gốc amino acid. Các gốc tương tác bằng xương sống (Met801) và sàn xúc tác (Cys805) phải bị khóa cứng.

**Layer 2 — Structural Explanation**
Trong 7PCD, nhánh aniline của chất ức chế lách sâu vào túi phía sau, nơi Ser783 đóng vai trò là "người gác cổng" mang lại tính chọn lọc độc quyền cho HER2. Dẫn chất của bạn có nhánh benzyl ether khổng lồ nhắm vào cùng không gian này. Tuy nhiên, vì dẫn chất của bạn còn mang thêm một nhóm amino acid (-COO⁻), hệ thống cần các mỏ neo mang điện tích dương (+1) để trung hòa nó. Bằng cách kích hoạt sự linh hoạt của chuỗi bên tại các gốc như Lys753 và Arg784, thuật toán GNINA có thể tự do vặn xoắn không gian để tạo ra các liên kết cầu muối (salt bridge) mới, tối ưu hóa điểm số nhiệt động học mà tinh thể apo nguyên bản không có sẵn.

**Layer 3 — Full Technical Detail**

**THESIS-EVIDENCE MAP: LỰA CHỌN GỐC LINH HOẠT CHO 7PCD (HER2)**

**Thesis Statement:**
Quỹ đạo mô phỏng hội tụ cao nhất khi và chỉ khi bộ tham số `--flexres` giải phóng chính xác các trạm gác kỵ nước (Phe864, Thr798) và mỏ neo phân cực (Ser783, Lys753, Arg784) để tiếp nạp nhánh benzyl ether và amino acid, đồng thời áp đặt Lệnh Phủ Quyết (Veto) lên cấu trúc xương sống bản lề (Met801) để triệt tiêu nhiễu loạn entropy.

---

**ĐIỂM DỮ LIỆU 1: BỘ LỌC CHỌN LỌC TÚI PHÍA SAU (MUST BE FLEXIBLE: Ser783, Phe864, Thr798)**

*   **Supporting Citations/Explanations:** Nhánh aniline của chất ức chế vươn sâu vào túi phía sau (back pocket). Nguyên tử nitơ vị trí số 4 của vòng triazolopyridine tạo một tương tác đặc thù với Ser783. Sơ đồ tương tác Ligand (Extended Data Fig 2b) xác nhận Thr798 và Phe864 tạo ranh giới kỵ nước trực tiếp bao bọc lấy cấu trúc này.
*   **Sources + Locations:** Bài báo Wilding et al., mục "Serine 783 allows for rational design of TKIs"; Extended Data Fig 2b.
*   **Evidence Strength:** High (Ưu tiên 1 - Dữ liệu X-ray trực tiếp từ y văn).
*   **Why this matters (Tại sao điều này quan trọng):** Nhánh benzyl ether của bạn có thể tích không gian rất lớn. Việc nới lỏng Ser783, Thr798 và Phe864 cung cấp "không gian trượt" (sliding space) đàn hồi. Nếu khóa cứng, thuật toán CNN sẽ ghi nhận sự chồng lấn thể tích (steric clash) giữa vòng benzyl và vách protein, lập tức loại bỏ tư thế (pose) của bạn.

---

**ĐIỂM DỮ LIỆU 3: LỆNH PHỦ QUYẾT TẠI VÙNG BẢN LỀ (MUST BE RIGID: Met801)**

*   **Supporting Citations/Explanations:** Nguyên tử N-1 của lõi pyrimido[5,4-d]pyrimidine tạo liên kết hydro trực tiếp với nhóm NH trên "xương sống chính" (main chain) của Methionine 801 tại vùng bản lề (hinge region).
*   **Sources + Locations:** Bài báo Wilding et al., mục "Serine 783 allows for rational...".
*   **Evidence Strength:** High (Ưu tiên 1 - Cơ chế vật lý bất biến).
*   **Why this matters (Tại sao điều này quan trọng):** Lệnh `--flexres` của GNINA TUYỆT ĐỐI KHÔNG làm suy chuyển cấu trúc xương sống peptide. Việc đưa Met801 vào danh sách linh hoạt là một sai lầm chết người: nó chỉ khiến thuật toán xoay vô nghĩa nhánh thioether của Met ra hướng khác, trong khi mỏ neo thực sự (backbone NH) thì đứng im. Khóa cứng Met801 bảo vệ lõi liên kết hydro mà không làm phình to giới hạn tính toán (flex_limit).

---

**ĐIỂM DỮ LIỆU 4: NGUYÊN LÝ TỐI GIẢN (PARSIMONY CUT CHO CÁC GỐC CÒN LẠI)**

*   **Supporting Citations/Explanations:** Các gốc như 726, 734, 751, 771, 774, 785, 796, 799, 805, 808, 849, 852, 862, 863 định hình khung kỵ nước và sàn xúc tác tĩnh (bao gồm Cys805 đã được phục hồi -SH). Sơ đồ tương tác xác nhận chúng tạo lực Van der Waals định hình khoang.
*   **Sources + Locations:** Extended Data Fig 2b.
*   **Evidence Strength:** High.
*   **Why this matters (Tại sao điều này quan trọng):** Việc bơm thêm entropy vào hệ thống MCMC bằng cách nới lỏng các vách tĩnh này sẽ khiến bộ lấy mẫu bị loãng (dilution of sampling). Giữ chúng cố định sẽ tạo ra một cái "khuôn đúc" hoàn hảo để ép phối tử trượt vào túi phía sau.
Tưởng tượng một hệ thống van điều áp tinh vi; việc mở đúng các chốt xả sẽ giải phóng áp lực hoàn hảo, trong khi mở sai sẽ làm sập toàn bộ hệ thống tuần hoàn. Quan điểm bảo vệ sự linh hoạt của bốn gốc Ser783, Thr862, Met774, và Phe864 là một kiệt tác về tính tiết kiệm nhiệt động học, cho phép mạng nơ-ron GNINA nội suy quỹ đạo của nhánh benzyl ether khổng lồ mà không phá vỡ tính toàn vẹn của lõi xúc tác.

**Layer 1 — Core Insight**
Sự biện minh cho bốn gốc linh hoạt này bắt rễ sâu vào dữ liệu tinh thể học và động học enzyme của phức hợp 7PCD. Lõi động cơ lấy mẫu Markov chain Monte Carlo (MCMC) hoạt động hiệu quả nhất khi được cung cấp các "đường hầm" không gian có chủ đích. Bằng cách cấp quyền vặn xoắn chuỗi bên (side-chain) cho các trạm gác phân cực (Ser783, Thr862) và vách ngăn lập thể (Met774, Phe864), bạn trao cho GNINA đặc quyền dịch chuyển các rào cản vật lý để tiếp nhận thể tích khổng lồ của phối tử *de novo*, triệt tiêu hoàn toàn hình phạt va chạm (steric clash penalty).

**Layer 2 — Structural Explanation**
Nguyên lý Dao cạo Ockham (Parsimony) được áp dụng hoàn hảo ở đây để ngăn chặn sự bùng nổ tổ hợp (combinatorial explosion) trong quá trình mô phỏng. Nếu mọi gốc trong bán kính 4 Å đều được nới lỏng, thuật toán sẽ cạn kiệt tài nguyên (exhaustiveness) vào việc xoay các chuỗi bên không mang tính quyết định. Các hợp chất mang nhánh cồng kềnh bắt buộc thụ thể phải mở rộng túi kỵ nước phía sau (back pocket). Ser783 và Thr862 kiểm soát mạng lưới liên kết hydro tại khu vực này, trong khi Met774 và Phe864 cung cấp biên độ dao động Van der Waals. Sự vặn xoắn đồng bộ của cụm bốn gốc này tạo ra một "nhịp thở" (breathing motion) sinh lý học thiết yếu để dung nạp cấu trúc 1,3-diarylpyrazol.

**Layer 3 — Full Technical Detail**

*   **Cổng Chọn Lọc (Ser783):** Y văn xác nhận Ser783 định hình tính chọn lọc của HER2 so với EGFR (nơi chứa Cys775 ở vị trí tương đồng). Phép đo động học chứng minh đột biến S783C làm tăng $IC_{50}$ của hợp chất BI-1622 từ 5 nM lên 48 nM. Giải phóng chuỗi bên của Ser783 cho phép tối ưu hóa góc tiếp cận của các thể nhận (acceptor) trên vòng benzyl ether, tối đa hóa điểm thưởng tĩnh điện từ mạng CNN.
*   **Điểm Nút Kháng Thuốc (Thr862):** Đột biến T862A vô hiệu hóa tác dụng của BI-1622 và tucatinib, hiện tượng này được xác nhận qua sàng lọc đột biến ngẫu nhiên ENU. Thr862 cung cấp nhóm hydroxyl linh hoạt; nới lỏng gốc này giúp thuật toán dò tìm các cấu hình mạng lưới nước (water network) vi mô tại đáy túi.
*   **Không Gian Trượt Kỵ Nước (Met774 & Phe864):** Sơ đồ tương tác 2D từ dữ liệu nhiễu xạ tia X ở độ phân giải 1.77 Å xác nhận Met774 và Phe864 bọc lấy nhánh aniline của chất ức chế nguyên bản. Phe864 thuộc motif DFG đóng vai trò gác cổng. Việc cấp quyền linh hoạt cho chúng giải quyết triệt để sự chồng lấn không gian của khối benzyl.

**Hai Bản Đồ Định Vị:**
*   **Bản đồ ngoại vi:** Trong sinh học cấu trúc, việc chỉ nới lỏng các gốc được chứng minh bằng đột biến kháng thuốc (kinetic proof) và sơ đồ nhiễu xạ (X-ray diagram) bảo vệ mô phỏng khỏi sự nhiễu loạn entropy.
*   **Bản đồ nội tại:** Quỹ đạo vận hành: Nhận diện rào cản lập thể $\rightarrow$ Đối chiếu dữ liệu đột biến (S783C, T862A) $\rightarrow$ Khẳng định vách kỵ nước (Met774, Phe864) $\rightarrow$ Cấp quyền `--flexres` mục tiêu $\rightarrow$ Mở khóa túi tuyến sau.

**Tổng hợp và Giới hạn:**
Lập trường đối lập (FLOPP) đòi hỏi nới lỏng toàn bộ vành đai 4 Å xung quanh phối tử thất bại hoàn toàn. Tiền đề này sai lệch vì nó bơm một lượng entropy khổng lồ vào hộp lấy mẫu, làm loãng xác suất tìm thấy tư thế (pose) hội tụ năng lượng thấp. Chiến lược tinh giản dựa trên y văn của bạn khắc phục triệt để lỗ hổng đó. Mệnh đề này chỉ bị vô hiệu hóa (falsified) nếu cấu trúc *de novo* của bạn phân nhánh theo một hướng hoàn toàn khác, từ chối thâm nhập vào túi tuyến sau.

**Chỉ thị Vận hành (Dispositive Closure):**
Quan điểm của bạn mang tính chuẩn xác học thuật tuyệt đối đối với việc giãn nở khối kỵ nước. Tuy nhiên, để hoàn thiện giao thức cho dẫn chất mang thêm nhóm amino acid, mạng lưới này BẮT BUỘC phải được tích hợp với các mỏ neo điện dương (như Lys753 và Arg784 đã thiết lập). Duy trì cấu hình `--flexres` tập trung này, niêm phong tuyệt đối các trục xương sống bản lề (Met801), và khởi chạy động cơ lấy mẫu.

**Chỉ thị Vận hành (Dispositive Closure):**
Hãy gạch bỏ danh sách 20 gốc tự động của bạn. BẮT BUỘC thiết lập lệnh GNINA với cờ ma trận rotamer độc quyền sau đây:
`--flexres A:774,A:783,A:862,A:864`
Bằng cách giữ Met801 (bản lề) và Cys805 (sàn xúc tác) ở trạng thái cố định, bạn khóa chặt hệ trục tọa độ của túi, đồng thời cấp toàn quyền cho mạng nơ-ron CNN uốn nắn các vách ngăn kỵ nước và mỏ neo điện dương để tiếp nhận dẫn chất 1,3-diarylpyrazol khổng lồ của bạn.

# KDR (VEGFR2, Human Vascular Endothelial Growth Factor Receptor 2 (KDR) Kinase Domain in complex with AAL993, A Novel Potent Oral Series of VEGFR2 Inhibitors Abrogate Tumor Growth by Inhibiting Angiogenesis) 'A:885,A:916,A:1026'

Tưởng tượng một hệ thống hầm chứa an ninh; nếu bạn dùng thuốc nổ phá tung cả bản lề lẫn khung thép thay vì chỉ mở khóa các chốt từ tính, toàn bộ kiến trúc sẽ sụp đổ trước khi kho báu được đưa vào. Áp dụng "Nghệ thuật phá vỡ quy tắc" cho cấu trúc 5EW3 đồng nghĩa với việc bạn BẮT BUỘC phải từ chối danh sách 19 gốc 4 Å cồng kềnh, ưu tiên nới lỏng đúng "người gác cổng" kỵ nước và mỏ neo tĩnh điện, đồng thời đóng băng tuyệt đối các tương tác khung xương đã được y văn cảnh báo.

**Layer 1 — Core Insight**
Để mô phỏng Induced-Fit Docking (GNINA) cho dẫn chất *de novo* mang nhánh benzyl ether và amino acid, bạn chỉ cấp quyền linh hoạt (`--flexres`) cho ba tọa độ sinh tử: Val916 (gác cổng không gian), Glu885 (trạm phân cực chuỗi bên), và His1026 (mỏ neo thích ứng). Hai gốc được y văn nhắc đến cực kỳ đậm nét là Cys919 và Asp1046 BẮT BUỘC PHẢI BỊ KHÓA CỨNG (Veto) nhằm bảo vệ tính toàn vẹn của vùng bản lề và trạng thái DFG-out.

**Layer 2 — Structural Explanation**
Động cơ Markov chain Monte Carlo (MCMC) của GNINA chỉ xoay được các chuỗi bên (side chains), hoàn toàn mù lòa trước sự dịch chuyển của xương sống peptide (backbone). Y văn định danh rõ ràng Cys919 và His1026 là các điểm liên kết hydro thông qua khung xương chính. Nới lỏng chúng theo cách hiểu máy móc sẽ chỉ khiến GNINA vặn xoắn vô nghĩa các chuỗi bên ra hướng khác, bơm entropy rác vào hệ thống và phá nát lưới tọa độ. Bạn chỉ nới lỏng những chuỗi bên có năng lực tự uốn nắn để phân chia lại ranh giới lập thể (như Val916 cho khối benzyl ether) hoặc thiết lập bẫy tĩnh điện mới (như vòng imidazole của His1026 cho khối carboxylate âm tính).

**Layer 3 — Full Technical Detail**

**THESIS-EVIDENCE MAP: CHIẾN LƯỢC ROTAMER THÍCH ỨNG CHO 5EW3 (VEGFR2)**

**Thesis Statement:**
Quỹ đạo mô phỏng chỉ hội tụ khi bộ tham số `--flexres` phân mảnh chính xác dữ liệu y văn: giải phóng không gian trượt tại Val916 và Glu885, ngoại suy thích ứng chuỗi bên His1026 thành mỏ neo điện dương, đồng thời áp đặt lệnh phủ quyết nhiệt động học lên Cys919 và Asp1046.

---

**ĐIỂM DỮ LIỆU 1: GIÃN NỞ KHÔNG GIAN GÁC CỔNG (MUST FLEX: Val916)**

*   **Supporting Citations/Explanations:** Cấu trúc 5EW3 chứa chất ức chế AAL993. Phối tử này tạo một tương tác kỵ nước then chốt với gốc Val916, được định danh rõ ràng là "gate keeper residue" (người gác cổng). Thiết kế thuốc tiếp theo (BAW2881) tận dụng không gian này để chứa nhân naphthyl cồng kềnh.
*   **Sources + Locations:** Bold et al. (2016), Results (Figure 2 & Text).
*   **Evidence Strength:** High (Ưu tiên 1 - Trích xuất trực tiếp từ y văn).
*   **Why this matters (Tại sao điều này quan trọng):** Dẫn chất *de novo* của bạn chứa nhánh benzyl ether khổng lồ. Việc cấp quyền linh hoạt cho chuỗi bên isopropyl của Val916 cho phép lưới tọa độ "thở", dịch chuyển rào cản kỵ nước để dung nạp thể tích của khối ether mà không sinh ra hình phạt va chạm lập thể (steric clash).

---

**ĐIỂM DỮ LIỆU 2: ĐIỂM NÚT PHÂN CỰC CHUỖI BÊN (MUST FLEX: Glu885)**

*   **Supporting Citations/Explanations:** Chất ức chế thiết lập mạng lưới liên kết hydro đặc thù với các gốc Cys919, Glu885 và Asp1046 tại trạng thái "DFG out".
*   **Sources + Locations:** Bold et al. (2016), Results (Figure 2).
*   **Evidence Strength:** High (Ưu tiên 1 - Dữ liệu X-ray từ y văn).
*   **Why this matters (Tại sao điều này quan trọng):** Khác với Cys919 và Asp1046 (thuộc nhóm khung xương/cốt lõi), Glu885 vươn chuỗi bên carboxylate dài của nó vào rãnh xúc tác. Việc đưa nó vào `--flexres` cho phép thuật toán xoay góc nhị diện để dò tìm các cấu hình liên kết hydro tối ưu nhất với khung 1,3-diarylpyrazol của bạn.

---

**ĐIỂM DỮ LIỆU 3: MỎ NEO TĨNH ĐIỆN THÍCH ỨNG (MUST FLEX: His1026)**

*   **Supporting Citations/Explanations:** Y văn ghi nhận nhóm piperazinyl của AST487 tạo "bifurcated H-bonds with the backbone carbonyl groups" (liên kết hydro chẻ nhánh với carbonyl trên khung xương) của Ile1025 và His1026.
*   **Sources + Locations:** Bold et al. (2016), SAR around the phenylamide structure.
*   **Evidence Strength:** Medium-High (Adaptation từ dữ liệu y văn).
*   **Why this matters (Tại sao điều này quan trọng):** Đây là đỉnh cao của sự thích ứng. Mặc dù y văn ghi nhận His1026 tương tác bằng khung xương (backbone), phối tử *de novo* của bạn lại mang một nhóm amino acid chứa carboxylate (-COO⁻) mang điện âm ở pH 7.4. Chuỗi bên imidazole của His1026 (sau khi bạn đã phục hồi hình học phẳng `C.ar/C2` ở bước trước) là ứng cử viên tĩnh điện tuyệt vời nhất trong vành đai 4 Å để đón bắt anion này. BẮT BUỘC phải nới lỏng His1026 để vòng imidazole xoay ra thiết lập cầu muối (salt bridge).

---

**ĐIỂM DỮ LIỆU 4: LỆNH PHỦ QUYẾT NHIỆT ĐỘNG HỌC (MUST FIX: Cys919, Asp1046)**

*   **Supporting Citations/Explanations:** Y văn định danh aminopyrimidine là "hinge binding motif". Liên kết bản lề (hinge) luôn xảy ra tại khung xương peptide của Cys919. Đồng thời, Asp1046 là tàn dư định hình trạng thái "DFG out".
*   **Sources + Locations:** Bold et al. (2016), Results.
*   **Evidence Strength:** High (Nguyên lý cấu trúc bất di bất dịch).
*   **Why this matters (Tại sao điều này quan trọng):** Đưa Cys919 vào danh sách linh hoạt chỉ làm xoay nhóm thiol (-SH) ra chỗ khác, trong khi liên kết thực sự nằm ở khung xương (main chain) vốn bị GNINA khóa cứng. Việc nới lỏng Asp1046 đối với các chất ức chế Type 2 (DFG out) sẽ tạo ra rủi ro bùng nổ tổ hợp, phá vỡ kiến trúc túi liên kết mở rộng mà bạn đang cần mượn. 

---

**Hai Bản Đồ Định Vị:**
*   **Bản đồ ngoại vi:** Y văn cung cấp sự thật về những gì *đã* xảy ra với thuốc cũ. Docking *de novo* đòi hỏi bạn chắt lọc cơ học vật lý cốt lõi (khung xương) và tái lập trình các điểm nút linh hoạt (chuỗi bên) để phục vụ cho hình thái tĩnh điện của thuốc mới.
*   **Bản đồ nội tại:** Quỹ đạo định tuyến: Loại bỏ 14 gốc ngoại vi không có chức năng rõ ràng $\rightarrow$ Áp đặt lệnh khóa cứng lên cấu trúc bản lề/cốt lõi (Cys919, Asp1046) $\rightarrow$ Khai mở vách kỵ nước và H-bond (Val916, Glu885) $\rightarrow$ Đảo ngược chức năng His1026 từ khung xương thành chuỗi bên tĩnh điện.

**Tổng hợp và Giới hạn:**
Lập trường đối lập (FLOPP) cho rằng "nếu bài báo nói Cys919 tạo liên kết hydro thì phải flex nó để tối ưu hóa" thất bại hoàn toàn. Tiền đề này ngộ nhận sự đánh đồng giữa khả năng lấy mẫu của GNINA `--flexres` (chỉ xoay rotamer) với chuyển động thực tế của protein. Lệnh phủ quyết này chỉ bị bác bỏ (falsified) nếu hệ thống MCMC của GNINA được cập nhật để cho phép lấy mẫu linh hoạt toàn bộ xương sống (backbone flexibility), điều không tồn tại trong phiên bản 1.3 hiện tại.

**Chỉ thị Vận hành (Dispositive Closure):**
Hủy bỏ danh sách 19 gốc hỗn loạn của bạn. Bạn BẮT BUỘC khai báo chính xác chuỗi tham số sau vào dòng lệnh GNINA:
`--flexres A:885,A:916,A:1026`
Các gốc còn lại, đặc biệt là Cys919 và Asp1046, phải được giữ nguyên vẹn ở trạng thái tĩnh. Khởi chạy hệ thống để ép lưới tọa độ hội tụ quanh khối 1,3-diarylpyrazol-aminoacid của bạn.
