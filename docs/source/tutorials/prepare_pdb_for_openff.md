# Prepare a PDB for OpenFF and PolyzyMD

In this tutorial you prepare a crystal structure from the Protein Data Bank
(PDB) for a PolyzyMD simulation. You use ubiquitin, PDB entry 1UBQ. At the
end, OpenFF reads the prepared file, and `polyzymd build` makes a solvated
system from it.

You learn these steps:

1. Read what a PDB entry contains before you change it.
2. Remove the crystallographic waters.
3. Add hydrogens with `polyzymd clean-pdb`.
4. Check the file with OpenFF.
5. Build the system with PolyzyMD.

1UBQ is a good first structure. It has one protein chain, no missing residues
and no alternate locations. Many PDB entries need more work: several copies of
the protein, missing residues, ligands or disulfide bonds. The last section of
this tutorial shows two of these steps on PDB entry 181L: a chain that ends
early, and a crystal ligand. For the other cases, see
{doc}`../how_to/troubleshoot_openff_pdb_ingestion` after this tutorial.

:::{admonition} Environment Setup
:class: tip

Run every command of this tutorial in the `build` environment. It holds
PolyzyMD, PDBFixer, OpenMM and OpenFF. From the repository root, activate it
once:

```bash
pixi shell -e build
```
:::

## Step 1: Download the structure

Make a working folder and download 1UBQ:

```bash
mkdir -p ubq_prep/structures
cd ubq_prep
curl -L https://files.rcsb.org/download/1UBQ.pdb -o structures/1UBQ.pdb
```

## Step 2: Read what the file contains

A PDB entry describes a crystal, not a simulation system. Before you change
the file, find out what is in it.

Count the protein atoms and the other records:

```bash
grep -c "^ATOM" structures/1UBQ.pdb
grep "^HETATM" structures/1UBQ.pdb | cut -c18-20 | sort | uniq -c
```

```
602
     58 HOH
```

The file has 602 protein atoms and 58 crystallographic waters (`HOH`). It has
no ligand.

Look for missing residues and missing atoms. The PDB lists them in `REMARK
465` and `REMARK 470` records:

```bash
grep -c "REMARK 465\|REMARK 470" structures/1UBQ.pdb
```

```
0
```

No residue and no heavy atom is missing. Ubiquitin has 76 residues, and the
last residue, Gly 76, has its terminal oxygen `OXT`. So the protein is
complete, and the file has no hydrogens yet.

```{note}
If a structure has `REMARK 465` or `REMARK 470` records, decide whether the
missing parts matter for your question. PolyzyMD does not model missing
residues or heavy atoms. See {doc}`../how_to/troubleshoot_openff_pdb_ingestion`.
```

## Step 3: Remove the waters

PolyzyMD adds its own water box. Remove the crystallographic waters:

```bash
grep -v " HOH " structures/1UBQ.pdb > structures/1UBQ_protein.pdb
```

## Step 4: Add hydrogens

`polyzymd clean-pdb` replaces nonstandard residues with standard residues and
adds the missing hydrogens at the pH that you give. It keeps the chain IDs and
residue numbers.

```bash
polyzymd clean-pdb -i structures/1UBQ_protein.pdb -o structures/ubq_clean.pdb --ph 7.0
```

```
Cleaning PDB: structures/1UBQ_protein.pdb
  pH: 7.0
  Adding missing hydrogens...

Cleaned PDB written to: structures/ubq_clean.pdb
```

Check the result. Every protein atom must be on chain `A`, the PolyzyMD
chain of the protein:

```bash
grep -c "^ATOM" structures/ubq_clean.pdb
grep "^ATOM" structures/ubq_clean.pdb | cut -c22 | sort | uniq -c
```

```
1231
   1231 A
```

The file now has 1231 atoms: the 602 heavy atoms and 629 hydrogens.

## Step 5: Check the file with OpenFF

PolyzyMD reads the protein with `openff.toolkit.Topology.from_pdb()`. OpenFF
matches each residue to a template, with its atoms, bonds and formal charges.
Run the same call directly:

```bash
python -c "
from openff.toolkit import Topology
topology = Topology.from_pdb('structures/ubq_clean.pdb')
charge = sum(atom.formal_charge.m for atom in topology.atoms)
print('OpenFF read', topology.n_molecules, 'molecule,', topology.n_atoms, 'atoms, net charge', charge)
"
```

```
OpenFF read 1 molecule, 1231 atoms, net charge 0
```

OpenFF read the protein as one molecule. The net charge of ubiquitin at pH 7
is 0. If this call fails, the file is not ready. See
{doc}`../how_to/troubleshoot_openff_pdb_ingestion`.

## Step 6: Build the system

Copy the quickstart config and point it at the prepared file. From the
`ubq_prep` folder, with `<repo>` the path of your PolyzyMD clone:

```bash
cp <repo>/examples/quickstart/config.yaml config.yaml
```

Edit these keys of `config.yaml`:

```yaml
name: "ubiquitin_water"
description: "Ubiquitin in water with NaCl"

enzyme:
  name: "ubiquitin"
  pdb_path: "structures/ubq_clean.pdb"
```

Validate the config, then build replicate 1:

```bash
polyzymd validate -c config.yaml
polyzymd build -c config.yaml -r 1
```

The build takes about 30 seconds. It ends with:

```
System built successfully!
Output directory: /home/me/ubq_prep/ubiquitin_300K_run1
```

The build can also print a warning that a few atoms lie between 1 and 2 Å of
a periodic image. Energy minimization removes these contacts before the
simulation starts.

## A chain that ends early, and a crystal ligand

Many PDB entries need two more steps. Entry 181L, T4 lysozyme L99A with a
bound benzene, needs both. Download it into a new folder:

```bash
mkdir -p t4l_prep/structures
cd t4l_prep
curl -L https://files.rcsb.org/download/181L.pdb -o structures/181L.pdb
grep -v "^HETATM" structures/181L.pdb > structures/181L_protein.pdb
```

The `grep` keeps the protein. It removes the waters, the ions, the
crystallization additives and the benzene.

(add-terminal-oxt)=
### Add the terminal oxygen of a chain that ends early

181L has no coordinates for its last two residues, 163 and 164 (`REMARK 465`).
So the chain ends at Lys 162, and Lys 162 has no terminal oxygen `OXT`.
`polyzymd clean-pdb` does not add it, and OpenFF then fails at the last
residue:

```
Input residue A:LYS#0162 contains atoms matching substructures {'PEPTIDE_BOND', 'NO MATCH'}
```

Add the `OXT` atom with PDBFixer before you run `clean-pdb`. Save this script
as `add_terminal_oxt.py`. It adds the terminal atoms and nothing else: no
missing residues and no missing side-chain atoms.

```python
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
```

```bash
python add_terminal_oxt.py structures/181L_protein.pdb structures/181L_oxt.pdb
polyzymd clean-pdb -i structures/181L_oxt.pdb -o structures/t4l_clean.pdb --ph 7.0
```

```
adding: {'LYS 162': ['OXT']}
```

Run the OpenFF check of Step 5 on `structures/t4l_clean.pdb`. It now prints
`OpenFF read 1 molecule, 2603 atoms, net charge 8`.

(ligand-sdf-from-crystal)=
### Make the ligand SDF from the crystal pose

The `substrate:` block of the config needs an SDF file with bond orders and
hydrogens. A PDB file has neither. RDKit, in the `build` environment, makes
the SDF from the `HETATM` records of the ligand and a SMILES string. The SDF
keeps the crystal coordinates, so the ligand stays in its pocket. Save this
script as `ligand_to_sdf.py`:

```python
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
```

Give the PDB file, the residue name of the ligand, its SMILES and the SDF to
write. The residue name of benzene in 181L is `BNZ`:

```bash
python ligand_to_sdf.py structures/181L.pdb BNZ "c1ccccc1" structures/benzene.sdf
```

```
structures/benzene.sdf: 12 atoms, formal charge 0
```

Write the SMILES in the charge state that you want to simulate. RDKit takes
the bond orders and the formal charges from it. For a symmetric molecule,
such as benzene, RDKit also prints `More than one matching pattern found`.
The matches are equivalent, so you can ignore this warning. If an entry has
several copies of the ligand, keep the lines of one copy only, for example
by its chain or residue number.

Point the config at the two files:

```yaml
enzyme:
  name: "t4l_l99a"
  pdb_path: "structures/t4l_clean.pdb"

substrate:
  name: "benzene"
  sdf_path: "structures/benzene.sdf"
  residue_name: "BNZ"
```

## What you did

You read a PDB entry, removed what PolyzyMD does not need, added hydrogens,
and checked the result with OpenFF and with a PolyzyMD build. The file
`structures/ubq_clean.pdb` is ready for a simulation.

## Next steps

- Run the system on this machine with `polyzymd run -c config.yaml -r 1`, as
  in {doc}`../get_started/quickstart`. Make the durations longer for a real
  study.
- Add polymers: {doc}`../how_to/polymers`.
- Prepare a harder structure, with several protein copies, cleaved chains or
  disulfide bonds: {doc}`../how_to/troubleshoot_openff_pdb_ingestion`.
