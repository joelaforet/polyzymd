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
the protein, missing residues, ligands or disulfide bonds. For those cases, see
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
