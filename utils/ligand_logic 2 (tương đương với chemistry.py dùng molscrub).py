import pandas as pd
from rdkit import Chem
from .chemistry import (
    ph_correct_smiles_molscrub,
    mol_to_3d_mol
)


def prepare_ligands(
    excel_file: str,
    smiles_col: str,
    id_col: str,
    output_sdf: str,
    ph: float = 7.4,
    seed: int = 42
):
    """
    Tier-3 ligand preparation with enforced physiological charge.

    Guarantees
    ----------
    - Dominant protonation state at pH 7.4 (MolScrub)
    - Explicit hydrogens
    - GNINA-compatible 3D SDF
    - One ligand → one charged state → one conformer
    """

    df = pd.read_excel(excel_file)
    writer = Chem.SDWriter(output_sdf)

    success, failed = 0, 0

    for _, row in df.iterrows():
        ligand_id = str(row[id_col])
        raw_smiles = str(row[smiles_col]).strip().replace("_x000d_", "")
        raw_smiles = raw_smiles.replace("\n", "").replace("\r", "")

        # --- Step 1: Thermodynamic fix (MolScrub) ---
        mol = ph_correct_smiles_molscrub(raw_smiles, ph)
        if mol is None:
            failed += 1
            continue

        # --- Step 2: Geometric fix (RDKit) ---
        mol = mol_to_3d_mol(mol, seed)
        if mol is None:
            failed += 1
            continue

        # --- Step 3: Annotation ---
        mol.SetProp("_Name", ligand_id)
        mol.SetProp("SMILES_raw", raw_smiles)
        mol.SetProp("Protonation_Method", "MolScrub")
        mol.SetProp("Target_pH", str(ph))

        writer.write(mol)
        success += 1

    writer.close()

    print(f"Prepared {success} ligands")
    print(f"Failed   {failed} ligands")

    return success
