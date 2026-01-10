import subprocess
from rdkit import Chem
from rdkit.Chem import AllChem


def ph_correct_smiles_openbabel(smiles: str, ph: float = 7.4) -> str | None:
    """
    Enforce dominant protonation state at target pH using OpenBabel.

    Chem LibreTexts basis:
    - Protonation state depends on pKa and solution pH
    - Dominant microspecies determines electrostatics

    Returns:
        pH-corrected canonical SMILES or None if failed
    """
    try:
        proc = subprocess.run(
            ["obabel", f"-:{smiles}", "-ismi", "-osmi", "-p", str(ph)],
            capture_output=True,
            text=True,
            check=True
        )
        # OpenBabel may append molecule name → take first token
        return proc.stdout.strip().split()[0]
    except Exception as e:
        print(f"[OpenBabel pH correction failed] {smiles} :: {e}")
        return None


def smiles_to_3d_mol(smiles: str, seed: int = 42) -> Chem.Mol | None:
    """
    Generate a low-energy 3D conformer from a chemically correct SMILES.

    GNINA constraint:
    - Bond lengths/angles are frozen during docking
    - Geometry must be minimized beforehand
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None

    # Explicit H AFTER protonation correction
    mol = Chem.AddHs(mol)

    params = AllChem.ETKDGv3()
    params.randomSeed = seed

    if AllChem.EmbedMolecule(mol, params) != 0:
        return None

    try:
        AllChem.MMFFOptimizeMolecule(mol)
    except ValueError:
        AllChem.UFFOptimizeMolecule(mol)

    return mol
