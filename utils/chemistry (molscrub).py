# Đoạn code này chuẩn bị file ligand (từ SMILES raw), khác với file chemistry.py trước (proton hóa bằng OpenBabel), thì code này sử dụng MolScrub (https://github.com/forlilab/molscrub) hoạt động dựa trên cốt lõi là RDKit
from rdkit import Chem
from rdkit.Chem import AllChem
from molscrub import Scrub
from rdkit.Chem.rdchem import KekulizeException, AtomValenceException

# ------------------------------------------------------------------
# Singleton MolScrub instance (load pKa rules ONCE)
# ------------------------------------------------------------------
_SCRUBBER = None


def _get_scrubber(ph: float) -> Scrub:
    global _SCRUBBER
    if _SCRUBBER is None:
        _SCRUBBER = Scrub(
            ph_low=ph,
            ph_high=ph,
            skip_acidbase=False,
            skip_tautomers=False,
            skip_ringfix=False,
            skip_gen3d=True,       # geometry handled by RDKit
            keep_all_frags=False  # drop salts / counterions
        )
    return _SCRUBBER


def ph_correct_smiles_molscrub(smiles: str, ph: float = 7.4) -> Chem.Mol | None:
    try:
        mol = Chem.MolFromSmiles(smiles, sanitize=False)
        if mol is None:
            return None

        # sanitize input nhẹ nhàng
        Chem.SanitizeMol(
            mol,
            sanitizeOps=Chem.SANITIZE_ALL ^ Chem.SANITIZE_KEKULIZE
        )

        scrubber = _get_scrubber(ph)
        states = scrubber(mol)

        if not states:
            return None

        for state in states:
            try:
                # 🚨 CỰC KỲ QUAN TRỌNG
                Chem.SanitizeMol(
                    state,
                    sanitizeOps=Chem.SANITIZE_ALL ^ Chem.SANITIZE_KEKULIZE
                )

                # Valence sanity check (optional but recommended)
                for atom in state.GetAtoms():
                    if atom.GetAtomicNum() == 8:  # Oxygen
                        if atom.GetTotalValence() > 2 and atom.GetFormalCharge() == 0:
                            raise AtomValenceException()

                return state  # first chemically sane state

            except (KekulizeException, AtomValenceException, ValueError):
                continue

        # No valid state survived
        return None

    except Exception:
        return None

def mol_to_3d_mol(
    mol: Chem.Mol,
    seed: int = 42
) -> Chem.Mol | None:
    """
    Generate GNINA-compatible 3D conformer.
    """

    try:
        # Clone mol to avoid mutating original (good practice)
        mol_3d = Chem.Mol(mol)

        # Ensure explicit hydrogens (safe even if already present)
        mol_3d = Chem.AddHs(mol_3d, addCoords=True)

        # ETKDGv3 embedding
        params = AllChem.ETKDGv3()
        params.randomSeed = seed
        params.useSmallRingTorsions = True

        if AllChem.EmbedMolecule(mol_3d, params) != 0:
            params.useRandomCoords = True
            if AllChem.EmbedMolecule(mol_3d, params) != 0:
                return None

        # Energy minimization
        try:
            if AllChem.MMFFHasAllMoleculeParams(mol_3d):
                AllChem.MMFFOptimizeMolecule(
                    mol_3d,
                    maxIters=200,
                    mmffVariant="MMFF94s"
                )
            else:
                AllChem.UFFOptimizeMolecule(mol_3d, maxIters=200)
        except Exception:
            # Fallback: keep embedded geometry
            pass

        return mol_3d

    except Exception as e:
        print(f"[3D generation failed] {e}")
        return None
