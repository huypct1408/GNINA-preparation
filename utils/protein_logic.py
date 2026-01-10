from pdbfixer import PDBFixer
from openmm.app import PDBFile

def prepare_protein(input_pdb, output_pdb, pH=7.4):
    """
    Chem LibreTexts basis:
    - Restores valence completeness
    - Assigns protonation based on bulk pH
    - Prepares topology for flexible docking
    """

    fixer = PDBFixer(filename=input_pdb)

    fixer.findMissingResidues()
    fixer.findNonstandardResidues()
    fixer.replaceNonstandardResidues()
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()

    # Protonation at physiological pH
    fixer.addMissingHydrogens(pH=pH)

    with open(output_pdb, "w") as f:
        PDBFile.writeFile(
            fixer.topology,
            fixer.positions,
            f,
            keepIds=True
        )

    return output_pdb
