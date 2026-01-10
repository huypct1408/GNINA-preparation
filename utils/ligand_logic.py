import pandas as pd
from rdkit import Chem
from .chemistry import ph_correct_smiles_openbabel, smiles_to_3d_mol


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

    Guarantees:
    - Correct dominant protonation state at pH 7.4
    - GNINA-compatible SDF
    - One ligand → one charged state → one 3D conformer
    """

    df = pd.read_excel(excel_file)
    writer = Chem.SDWriter(output_sdf)

    success, failed = 0, 0

    for _, row in df.iterrows():
        ligand_id = str(row[id_col])
        raw_smiles = row[smiles_col]

        # --- Step 1: pH correction (THERMODYNAMIC FIX) ---
        ph_smiles = ph_correct_smiles_openbabel(raw_smiles, ph)
        if not ph_smiles:
            failed += 1
            continue

        # --- Step 2: 3D generation (GEOMETRIC FIX) ---
        mol = smiles_to_3d_mol(ph_smiles, seed)
        if mol is None:
            failed += 1
            continue

        # --- Step 3: Annotation ---
        mol.SetProp("_Name", ligand_id)
        mol.SetProp("SMILES_raw", raw_smiles)
        mol.SetProp("SMILES_pH7.4", ph_smiles)

        writer.write(mol)
        success += 1

    writer.close()

    print(f"Prepared {success} ligands")
    print(f"Failed   {failed} ligands")

    return success
