# Danh sách chất
1. 3-((4H-benzo[d][1,3,2]dioxaborinin-2-yl)methyl)-1-(4-(benzyloxy)-3-methoxybenzyl)-2-methoxy-1,2-azaborolidine

2. 2-((2-methoxy-2,3-dihydrobenzo[d][1,2]oxaborol-3-yl)methyl)-4H-benzo[d][1,3,2]dioxaborinine

# Đích

# **1. PPARA (7BQ2)** `A:273, A:275, A:276, A:279`

### 🌟 4 Axit amin linh động BẮT BUỘC (Literature-Verified Tier)
Nhóm này giải quyết trực tiếp vấn đề chênh lệch kích thước (size-mismatch) để 2 chất của bạn có thể lọt được vào túi Arm III:

1. **Phe273 (Phenylalanine 273):** Đóng vai trò là "người gác cổng" (gatekeeper). Nó bắt buộc phải được phép lật/xoay (flip) để tạo khoảng trống cho cái đuôi *benzyloxy* rất dài của chất thứ nhất chui vào sâu trong túi kỵ nước.
2. **Cys275 (Cysteine 275)** và 
3. **Cys276 (Cysteine 276):** Đây là hai bản lề cơ học (mechanical hinges) của túi gắn kết PPARα. Việc cài đặt linh động 2 gốc này cho phép toàn bộ cấu trúc túi (LBD) phình to ra một cách tự nhiên để ôm trọn cấu trúc lõi *dioxaborinine* cồng kềnh của cả 2 chất. Nếu để cứng (rigid), phần mềm sẽ báo lỗi va chạm (steric clash).
4. **Thr279 (Threonine 279):** Đây là mỏ neo (anchor) cực kỳ quan trọng nằm sâu trong nhánh **Arm III**. Cấu trúc của Thr279 rất linh hoạt và sẽ tự động xoay để tạo liên kết hydro (trực tiếp hoặc qua cầu nối nước) với các nguyên tử O, N, B trên lõi phân tử của bạn, từ đó quyết định hoạt lực và độ chọn lọc cực cao.

# **2. PPARD (7WGN)** `A:253, A:287, A:413, A:437`

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

# **3. PPARG (8ATY)** `A:288,A:342,A:259`

Traditional PPARγ therapeutics encounter a biological paradox. Maximal insulin sensitization intrinsically couples with severe pro-adipogenic adverse effects. Resolving this requires shifting focus toward noncanonical allosteric modulation. Topologically, the receptor achieves this via a secondary solvent-exposed cleft bounded by Helix 4, Helix 5, and the Omega loop. 

The conventional premise dictates that PPARγ activation strictly requires direct physical coordination of the Helix 12 AF-2 surface. This position fails because crystallographic evidence (PDB: 8ATZ) demonstrates that alternative site ligands stabilize the active Helix 12 conformation allosterically through the Omega-loop, entirely bypassing direct AF-2 contact. This conclusion would be invalidated if structural variants of your bulky azaborolidines fully occlude the orthosteric pocket while generating classical pro-adipogenic transcriptional profiles.

To properly dock your large 1,2-azaborolidine and dioxaborinine derivatives, you must target the alternative binding site. Static algorithms will fail to capture the required binding volume. You must implement the following literature-verified flexible residues.

### Literature-Tier Flexible Residues

1. **Arg288**: Acts as the primary electrostatic anchor in the alternative site. It dynamically shifts to form direct ionic interactions or water-mediated hydrogen bonds with polar moieties, such as boronic or oxaborole oxygen atoms.
2. **Ser342**: Provides essential backbone and side-chain hydrogen bonding. It shifts to stabilize complex aromatic extensions within the alternative pocket cleft.
3. **Glu259**: Regulates the solvent-exposed periphery. Its side chain establishes critical hydrogen bonds with secondary amines or ether oxygen linkers spanning the alternative pocket.
4. **Ser289, His449, and Tyr473**: The classical orthosteric triad. If your compounds exhibit multi-site occupancy (1:2 stoichiometry like JP85/Compound 1 in PDB 8ATY), these residues must adapt to stabilize the canonical interactions.

### Argument-Evidence Map

| Argument | Supporting Citations/Explanation | Source + Location | Strength of Evidence (Low/Medium/High) | Assumption | Counter-Evidence |
| :--- | :--- | :--- | :--- | :--- | :--- |
| **PPARγ accommodates simultaneous multi-ligand occupancy across distinct sites.** | The LBD permits orthosteric ligands (pioglitazone) and alternative site ligands (Compound 2) to bind concurrently (1:1:1 stoichiometry) without steric competition. | Arifi et al., 2023 | High | Concurrent binding generates synergistic allosteric stabilization rather than independent silent states. | Specific extended ligands (MRL-871) occupy boundary regions overlapping both pockets, blocking secondary ligand entry. |
| **Helix 12 (AF-2) stabilization does not require direct orthosteric contact.** | Compound 2 induces an active Helix 12 state solely through alternative site binding, likely propagated via Omega-loop stabilization. | Arifi et al., 2023 | High | Omega-loop ordering mechanically couples to Helix 12 inward dynamics. | Alternative site ligands like BVT.13 bind the pocket but fail to trigger the active Helix 12 inward swing. |
| **Selective alternative site modulation prevents canonical adipogenesis.** | Targeting the H4/H5 cleft avoids the classical pro-adipogenic gene transcription profiles (e.g., CPT-1, FABP1) triggered by TZD therapies. | Arifi et al., 2023 | High | *In vivo* cellular differentiation pathways directly mirror these *in vitro* mRNAseq hepatocyte profiles. | The alternative site ligand still retains baseline recruitment of select canonical co-activators like PGC-1α. |
| **Alternative site activation uniquely suppresses FOXO signaling.** | The noncanonical binding mode diminishes FOXO3 expression and enhances inactivating Ser253 phosphorylation, driving cells toward a resting state. | Arifi et al., 2023 | Medium | FOXO suppression operates exclusively through genomic PPARγ regulation rather than secondary kinase off-targets. | Orthosteric PPARγ agonists like pioglitazone exert zero regulatory influence on this specific FOXO pathway. |

**Dispositive Closure:**
Configure your docking grid to explicitly encompass the Helix 4/Helix 5 cleft. Mandate side-chain flexibility for Arg288, Glu259, and Ser342. Evaluate your boronic poses based on their ability to stabilize the Omega loop rather than relying exclusively on direct Helix 12 contacts. This topology defines the absolute boundary between classical adipogenesis and selective noncanonical modulation.
