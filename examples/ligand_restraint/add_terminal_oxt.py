"""Add the terminal OXT atoms of a protein chain, and nothing else."""

import sys

from openmm.app import PDBFile
from pdbfixer import PDBFixer

fixer = PDBFixer(filename=sys.argv[1])
fixer.findMissingResidues()
fixer.missingResidues = {}  # do not model missing residues
fixer.findMissingAtoms()
fixer.missingAtoms = {}  # do not model missing side-chain atoms
print("adding:", {f"{res.name} {res.id}": atoms for res, atoms in fixer.missingTerminals.items()})
fixer.addMissingAtoms()
with open(sys.argv[2], "w") as handle:
    PDBFile.writeFile(fixer.topology, fixer.positions, handle, keepIds=True)
