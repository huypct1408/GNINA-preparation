# Danh sách chất
1. 3-((4H-benzo[d][1,3,2]dioxaborinin-2-yl)methyl)-1-(4-(benzyloxy)-3-methoxybenzyl)-2-methoxy-1,2-azaborolidine

2. 2-((2-methoxy-2,3-dihydrobenzo[d][1,2]oxaborol-3-yl)methyl)-4H-benzo[d][1,3,2]dioxaborinine

# Đích

**1. PPARA (7BQ2)** `A:273, A:275, A:276, A:279`

### 🌟 4 Axit amin linh động BẮT BUỘC (Literature-Verified Tier)
Nhóm này giải quyết trực tiếp vấn đề chênh lệch kích thước (size-mismatch) để 2 chất của bạn có thể lọt được vào túi Arm III:

1. **Phe273 (Phenylalanine 273):** Đóng vai trò là "người gác cổng" (gatekeeper). Nó bắt buộc phải được phép lật/xoay (flip) để tạo khoảng trống cho cái đuôi *benzyloxy* rất dài của chất thứ nhất chui vào sâu trong túi kỵ nước.
2. **Cys275 (Cysteine 275)** và 
3. **Cys276 (Cysteine 276):** Đây là hai bản lề cơ học (mechanical hinges) của túi gắn kết PPARα. Việc cài đặt linh động 2 gốc này cho phép toàn bộ cấu trúc túi (LBD) phình to ra một cách tự nhiên để ôm trọn cấu trúc lõi *dioxaborinine* cồng kềnh của cả 2 chất. Nếu để cứng (rigid), phần mềm sẽ báo lỗi va chạm (steric clash).
4. **Thr279 (Threonine 279):** Đây là mỏ neo (anchor) cực kỳ quan trọng nằm sâu trong nhánh **Arm III**. Cấu trúc của Thr279 rất linh hoạt và sẽ tự động xoay để tạo liên kết hydro (trực tiếp hoặc qua cầu nối nước) với các nguyên tử O, N, B trên lõi phân tử của bạn, từ đó quyết định hoạt lực và độ chọn lọc cực cao.

**2. PPARD (7WGN)** `A:253, A:287, A:413, A:437`

**Layer 1 — Core Insight**
Designing multi-target boronic ligands presents a topological challenge when navigating the promiscuous peroxisome proliferator-activated receptor (PPAR) ligand-binding domains. The 7WGN co-crystal structure demonstrates that pemafibrate binds the Y-shaped Center and Arm II/III regions of PPAR delta. To achieve efficacy, your size-mismatched azaborolidine and dioxaborinine molecules must satisfy the specific consensus hydrogen-bond network of PPAR delta while structurally displacing the spatial volume near helix 5.

**Layer 2 — Structural Explanation**
The literature does not explicitly classify discrete residues as "flexible" within the 7WGN dataset. We execute an adaption-tier protocol to define the obligatory dynamic coordinates for your bulky structures. First, your models must establish local plasticity at the consensus tetrad: Thr253, His287, His413, and Tyr437. These residues act as adaptive anchors that secure the ligand headgroup to stabilize the AF-2 helix 12. Second, structural accommodation in PPAR delta explicitly requires shifting the phenoxyalkyl-equivalent groups toward helix 5. Third, the binding pocket utilizes a water-mediated hydrogen bond or electrostatic interaction to stabilize the 2-aminobenzoxazole group. You must configure your docking grids to allow spatial adjustments at these precise locations to prevent steric clashes with your expansive benzyloxy and oxaborole units.

**Layer 3 — Full Technical Detail**

*External Map:* Fibrate cross-reactivity redefines therapeutic boundaries. Repurposing these molecules bridges dyslipidemia treatments directly to non-alcoholic fatty liver disease interventions. 
*Internal Map:* Receptor activation strictly requires aligning the ligand head with the AF-2 helix 12 consensus tetrad. Concurrently, the hydrophobic tail must occupy the specified Y-shaped Arm II/III geometry.

### Argument-Evidence Map

| Argument | Supporting Citations/Explanation | Source + Location | Strength of Evidence | Assumption | Counter-Evidence |
| :--- | :--- | :--- | :--- | :--- | :--- |
| **Pemafibrate functions as a structural pan-agonist across PPAR subtypes.** | Pemafibrate activates PPAR delta/gamma LBD-mediated transactivation, recruits coactivators, and physically binds PPAR delta/gamma LBDs. |, | High | Biochemical transactivation and co-crystallization successfully translate into physiological pan-agonism. | Pemafibrate requires concentrations three orders of magnitude higher to activate PPAR delta/gamma compared to PPAR alpha. |
| **PPAR delta imposes a unique spatial orientation on specific ligand substituents.** | In the 7WGN structure, the phenoxyalkyl group of pemafibrate is specifically pushed toward helix 5 within the PPAR delta pocket. | | High | Spatial deviation toward helix 5 is mandatory to resolve steric limitations unique to the delta subtype. | The fundamental Y-shaped binding geometry across the Center and Arm II/III regions remains largely identical across all three PPAR subtypes. |
| **AF-2 helix 12 stabilization requires a strict four-residue consensus network.** | Carboxylic groups stabilize AF-2 helix 12 formation through direct hydrogen bonds and electrostatic interactions with Thr253, His287, His413, and Tyr437 in PPAR delta. | | High | Boronic or carboxylic warheads must coordinate this exact tetrad to induce coactivator recruitment. | Fenofibric acid stabilizes PPAR delta thermostability but fails to functionally activate the receptor or produce cocrystals. |
| **Ligand efficacy demands hydrophobic interactions within the core pockets.** | Pemafibrate is further stabilized by extensive hydrophobic interactions within the LBDs of all PPARs. |, | Medium | Hydrophobic packing provides the primary thermodynamic driving force for complex formation. | Water-mediated hydrogen bonds actively participate in anchoring specific functional groups in PPAR delta. |

**Synthesis:** 
Modulating coordinate flexibility at the consensus tetrad dictates the boundary between partial agonism and highly selective receptor activation. The lack of specific structural flexibility forces compounds into inactive conformations. Fenofibric acid increases PPAR delta thermostability but fails to recruit coactivators, proving that simple binding does not guarantee functional activation. 

**Dispositive Closure:** 
When docking your 1,2-azaborolidine and dioxaborinine derivatives against 7WGN, you must mandate side-chain flexibility for Thr253, His287, His413, and Tyr437. You must direct your structural appendages toward the space bounded by helix 5. Failure to model this specific plasticity will yield false-negative collision errors.
