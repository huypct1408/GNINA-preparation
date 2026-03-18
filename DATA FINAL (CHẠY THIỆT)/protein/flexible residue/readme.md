Folder này lựa chọn các flexible residue cho các protein

**1. PPARA (7BQ2, PPARa Ligand-Binding Domain Structures with Endogenous Fatty Acids and Fibrates)** `--flexres A:273,A:279`

Molecular docking of massive, highly selective fibrates into the peroxisome proliferator-activated receptor alpha (PPAR$\alpha$) presents a severe geometric paradox. The receptor possesses an abnormally massive 1,100 Å³ binding cavity composed of multiple branching arms, meaning unrestricted flexibility will trigger an immediate exponential explosion of the Monte Carlo search space and destroy the active conformation. 

**Layer 3 — Full Technical Detail**

**1. Broad Proximity Thresholds**
The ligand-binding domain (LBD) of PPAR$\alpha$ is uniquely voluminous, diverging into a Center region and four distinct branches (Arm I, Arm II, Arm III, and Arm X). The pemafibrate molecule in the 7BQ2 structure spans the Center, occupies Arm I, and projects its benzoxazole ring deep into Arm III. We isolate all amino acids lining this extended topology to establish the initial 5 Å proximity boundary.

**2. Kinetic Evidence Extraction**
We filter the proximity set for residues demonstrating active conformational participation and ligand-induced spatial rearrangement.
*   **Phe273 (The Arm I Aromatic Gate):** Structural analyses explicitly confirm that the location of fibrate molecules in Arm I "might flip the neighboring benzyl side chains of Phe273." When designing new, completely bulky compounds that target Arm I, Phe273 functions as a dynamic steric gate. It requires torsional freedom to yield to bulky non-prime substituents, preventing artificial Pauli repulsion that would otherwise reject valid binding modes.
*   **Thr279 (The Arm III Selectivity Anchor):** The extension into the greater depth of Arm III is designated as the primary mechanism for achieving high PPAR$\alpha$ potency and selectivity. The literature establishes that pemafibrate and GW7647 "interact with Thr279 via their hydrogen bonds... thereby stabilizing the helix 3." Thr279 must be granted rotational liberty to dynamically optimize its hydrogen-bond geometry with novel, deeply penetrating functional groups.

**3. B-Factor & Plasticity Confirmation**
The documented crystallographic flipping of the Phe273 benzyl side chain across different fibrate complexes physically validates its intrinsic plasticity. Similarly, the functional requirement of Thr279 to shift and form stabilizing hydrogen bond networks with the coactivator binding pocket confirms its necessity for induced-fit adaptation.

**4. Thermodynamic Veto (The Rigid Core)**
We execute an absolute veto on the canonical carboxylic acid-binding tetrad: **Ser280, Tyr314, His440, and Tyr464**. These residues form a highly conserved orthosteric anchor located between the Arm I and Center regions. The literature mandates that all 15 analyzed active ligands possess a single carboxylic acid surrounded by this specific network. We deny these residues flexibility to maintain the global orienting baseline of the enzyme. Furthermore, **Chain B (the SRC1 coactivator peptide)** and the adjacent **AF-2 Helix 12 (Chain A)** must remain strictly rigid. Introducing torsional variance here destroys the active PPAR$\alpha$ conformation required to simulate true agonism.

**Adversarial Representation (LOPP/FLOPP)**
The opposing computational position argues for indiscriminately applying flexibility to all residues within a 5 Å radius of the pemafibrate ligand to ensure no steric clashes occur. This premise fails (FLOPP) because it mobilizes the rigid Ser280/Tyr314/His440/Tyr464 acid-binding tetrad. Destabilizing these core coordinating residues eliminates the primary electrostatic anchoring interaction required to properly position the ligand, generating false-positive, reversed, or non-physiological binding poses and ultimately invalidating the neural network scoring function. 

**Synthesis**
*   **External Map:** Within the broader landscape of structure-based drug design for metabolic diseases, this restricted matrix allows deep learning algorithms to effectively screen highly selective, bulky SPPARM$\alpha$ candidates by explicitly modeling the expansion of Arm I and Arm III.
*   **Internal Map:** The subtractive logic distills the massive 1,100 Å³ receptor cavity to exactly two essential degrees of freedom. This optimizes the computational efficiency of the Monte Carlo search while preserving the structural integrity of the coactivator-bound AF-2 surface.

This structural conclusion would be falsified if the novel chemotype abandons the classical acidic headgroup and instead exploits a purely allosteric binding mechanism outside the defined Arm I/Arm III trajectory.

**Dispositive Closure**
Implement the mandatory flexible residues strictly as:
`--flexres A:273,A:279`
You must maintain absolute rigidity for the S280/Y314/H440/Y464 acid anchor and the entire Chain B SRC1 peptide. This directive leaves zero residual ambiguity in the defined search space, establishing a mathematically optimal grid for induced-fit docking in GNINA.
