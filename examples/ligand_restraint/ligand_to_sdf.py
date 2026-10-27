"""Write one ligand of a PDB file to an SDF, at its crystal coordinates."""

import sys

from rdkit import Chem
from rdkit.Chem import AllChem

pdb_path, residue_name, smiles, sdf_path = sys.argv[1:]

# The HETATM records of the ligand hold its crystal pose.
records = [
    line
    for line in open(pdb_path)
    if line.startswith("HETATM") and line[17:20].strip() == residue_name
]
crystal = Chem.MolFromPDBBlock("".join(records))

# A PDB file has no bond orders: take them from the SMILES.
ligand = AllChem.AssignBondOrdersFromTemplate(Chem.MolFromSmiles(smiles), crystal)

# Add hydrogens, placed from the heavy-atom coordinates.
ligand = Chem.AddHs(ligand, addCoords=True)
Chem.MolToMolFile(ligand, sdf_path)
print(f"{sdf_path}: {ligand.GetNumAtoms()} atoms, formal charge {Chem.GetFormalCharge(ligand)}")
