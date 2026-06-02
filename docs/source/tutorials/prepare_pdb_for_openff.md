# Prepare a PDB for OpenFF and PolyzyMD

Raw Protein Data Bank files often describe crystallographic assemblies, waters,
alternate locations, missing-coordinate records, and historical chain labels that
are not the same thing as a simulation-ready PolyzyMD enzyme input. This tutorial
shows a conservative preparation workflow for the raw 4CHA alpha-chymotrypsin PDB
using PyMOL for biological-system inspection and a small PDBFixer script for
mechanical cleanup.

The workflow deliberately adds **missing hydrogens only**. It does not add
missing amino-acid residues or missing heavy atoms. If your structure is missing
protein residues or heavy atoms that are required for your science, model and
validate them before using PolyzyMD. Other structure-preparation tools can be
appropriate, but this tutorial prescribes PyMOL plus PDBFixer so the exact
selection and cleanup steps are visible.

## What you will produce

Starting from a raw RCSB PDB download:

```text
structures/4CHA.pdb
```

The intended final enzyme file has:

- one alpha-chymotrypsin enzyme copy
- all protein atom records using chain ID `A`
- TER records preserved between the disconnected peptide fragments
- no crystallographic waters
- no second enzyme copy
- explicit hydrogens
- simple structural checks that pass before OpenFF chemistry assignment
- successful `openff.toolkit.Topology.from_pdb()` validation after any required
  upstream residue, atom, protonation, or connectivity curation is complete

```{important}
PolyzyMD-prepared structures must use chain ID `A` for all protein residues.
Substrates use chain `B`, polymers use chain `C`, and solvent chains use `D` and
later. It is acceptable for one enzyme to contain multiple disconnected peptide
fragments after OpenFF ingestion. MD depends on chemical connectivity, not on
whether every protein fragment has a unique chain label.
```

## Prerequisites

- PolyzyMD installed in a pixi environment that includes PDBFixer, OpenMM, and
  OpenFF
- PyMOL for interactive inspection
- the raw 4CHA PDB file downloaded into `structures/4CHA.pdb`
- a scientific decision about whether any missing coordinates must be modeled
  before simulation

Create a clean working directory and download the raw PDB file:

```bash
mkdir -p 4cha_pdb_prep/structures
cd 4cha_pdb_prep
curl -L https://files.rcsb.org/download/4CHA.pdb -o structures/4CHA.pdb
```

## Inspect the raw 4CHA biological contents

Open the raw file in PyMOL:

```text
load structures/4CHA.pdb, raw4cha
hide everything, raw4cha
show cartoon, polymer.protein
show sticks, polymer.protein and resn CYS
show spheres, solvent

select first_copy, raw4cha and polymer.protein and chain A+B+C
select second_copy, raw4cha and polymer.protein and chain E+F+G
select crystal_waters, raw4cha and solvent

color marine, first_copy
color orange, second_copy
color cyan, crystal_waters
zoom first_copy
```

Read the PDB header before deleting anything. The relevant 4CHA records say:

- COMPND assigns one alpha-chymotrypsin enzyme copy across chains `A+B+C`
- COMPND assigns a second enzyme copy across chains `E+F+G`
- mature alpha-chymotrypsin is a cleaved enzyme with multiple peptide fragments
- residues 14-15 and 147-148 are biologically excised activation peptides, so
  TER records between the mature fragments are expected
- REMARK 465 reports missing-coordinate residues in the deposited model, including
  `GLY A 12` and `LEU A 13` in the first copy

Those REMARK 465 records are not a PDBFixer to-do list for this workflow. They
are a signal that the raw crystal model is not automatically simulation-ready. If
those missing coordinates matter for your study, provide a curated PDB where they
have already been modeled and checked. Do not rely on PolyzyMD or the helper
script below to make that biological modeling decision.

## Preview the system selection in PyMOL

Use PyMOL to confirm that the first copy is the enzyme you intend to simulate:

```text
disable second_copy
enable first_copy
show cartoon, first_copy
show sticks, first_copy and resn CYS
distance disulfides, first_copy and name SG, first_copy and name SG, 2.2
zoom first_copy
```

If you want to preview the trimmed system without changing the final preparation
path, create a temporary PyMOL object:

```text
create 4cha_first_copy_preview, raw4cha and polymer.protein and chain A+B+C
disable raw4cha
enable 4cha_first_copy_preview
```

Use this preview for visual inspection only. The final file is written by the
script below so the TER boundaries between fragments are preserved as separate
PDB topology chains while every atom record is relabeled to chain `A` for
PolyzyMD.

## Add hydrogens and write the PolyzyMD PDB

Save this script as `prepare_4cha_for_polyzymd.py` and run it from a PolyzyMD
pixi shell. It removes the second enzyme copy and waters, keeps the first copy's
three protein fragments, relabels all remaining protein chains to `A`, and adds
hydrogens only.

```python
from pathlib import Path

from openmm.app import PDBFile
from pdbfixer import PDBFixer


RAW_PDB = Path("structures/4CHA.pdb")
OUTPUT_PDB = Path("structures/4cha_chymotrypsin_chain_a_openff.pdb")
KEEP_CHAINS = {"A", "B", "C"}
PH = 7.0


def main() -> None:
    """Prepare one 4CHA enzyme copy for PolyzyMD."""
    fixer = PDBFixer(filename=str(RAW_PDB))

    chains_to_remove = [
        chain_index
        for chain_index, chain in enumerate(fixer.topology.chains())
        if chain.id not in KEEP_CHAINS
    ]
    fixer.removeChains(chains_to_remove)

    fixer.removeHeterogens(keepWater=False)
    fixer.addMissingHydrogens(PH)

    for chain in fixer.topology.chains():
        chain.id = "A"

    OUTPUT_PDB.parent.mkdir(parents=True, exist_ok=True)
    with OUTPUT_PDB.open("w") as handle:
        PDBFile.writeFile(fixer.topology, fixer.positions, handle, keepIds=True)

    print(f"Wrote {OUTPUT_PDB}")


if __name__ == "__main__":
    main()
```

Run it:

```bash
python prepare_4cha_for_polyzymd.py
```

```{warning}
This script intentionally does not call `findMissingResidues()`,
`findMissingAtoms()`, `addMissingAtoms()`, or any missing-residue filling method.
If OpenFF later reports missing chemistry, stop and curate the input structure
upstream instead of turning on automatic residue or heavy-atom construction here.
```

## Validate the prepared file with OpenFF

Run this validation snippet before referencing the PDB from `config.yaml`:

```python
from pathlib import Path

from openff.toolkit import Topology


PDB_PATH = Path("structures/4cha_chymotrypsin_chain_a_openff.pdb")


def main() -> None:
    """Validate the prepared PDB for PolyzyMD and OpenFF."""
    lines = PDB_PATH.read_text().splitlines()
    atom_records = [line for line in lines if line.startswith(("ATOM", "HETATM"))]
    ter_count = sum(line.startswith("TER") for line in lines)
    chains = {line[21] for line in atom_records}
    residue_names = {line[17:20].strip() for line in atom_records}
    elements = {line[76:78].strip() for line in atom_records if len(line) >= 78}

    if chains != {"A"}:
        raise SystemExit(f"Expected only chain A, found {sorted(chains)}")
    if ter_count != 3:
        raise SystemExit(f"Expected three TER records, found {ter_count}")
    if "HOH" in residue_names:
        raise SystemExit("Expected no crystallographic waters")
    if "H" not in elements:
        raise SystemExit("Expected explicit hydrogens")

    print("Structural PDB checks passed")
    print(f"  Atom records use chains: {sorted(chains)}")
    print(f"  TER records: {ter_count}")
    print("  Crystallographic waters: absent")
    print("  Explicit hydrogens: present")

    print("Running OpenFF chemistry validation")
    Topology.from_pdb(str(PDB_PATH))
    print("OpenFF chemistry assignment succeeded")


if __name__ == "__main__":
    main()
```

Expected success for a fully curated input is:

```text
Structural PDB checks passed
  Atom records use chains: ['A']
  TER records: 3
  Crystallographic waters: absent
  Explicit hydrogens: present
Running OpenFF chemistry validation
OpenFF chemistry assignment succeeded
```

When this exact script is run on the uncurated raw 4CHA file, the chain, TER,
water, and hydrogen checks pass, but OpenFF can still reject the file. In a
collaborator reproduction, the mechanically prepared PDB had all atom records on
chain `A`, three TER records, waters removed, the second enzyme copy removed, and
hydrogens present, while OpenFF still reported chemistry errors around terminal
CYS#0001 hydrogen and disulfide assignment and SER#0011. That failure is useful:
it separates successful structural preparation from unsuccessful chemistry
assignment and tells you to stop and provide a curated protein model before using
PolyzyMD.

If the chain, TER, water, or hydrogen checks fail, fix the preparation script or
your file paths. If `Topology.from_pdb()` fails after those simple checks pass,
OpenFF is telling you that the file still has unresolved chemistry such as
missing residues, missing heavy atoms, unsupported atom naming, or ambiguous
connectivity. Resolve those issues in a curated input PDB before running
PolyzyMD.

## What curate upstream means

OpenFF failures after the simple structural checks pass are chemistry-modeling
problems in the enzyme input, not PolyzyMD configuration problems. Curating
upstream means producing an externally inspected PDB whose residue templates,
atom names, protonation states, coordinates, and connectivity are chemically
consistent before PolyzyMD loads it.

Common issues to inspect include:

- terminal atom naming and terminal hydrogen counts
- disulfide cysteine protonation and SG-SG connectivity
- missing-coordinate fragments reported in the PDB header
- missing heavy atoms in otherwise retained residues
- biologically excised residues that should remain absent, with TER records kept
  between the mature fragments

Possible curation paths include manual editing and inspection in PyMOL or
ChimeraX, deliberate PDBFixer use with careful review of every modeled atom,
Amber-oriented cleanup with `pdb4amber`, or rebuilding missing structural regions
with tools such as MODELLER, AlphaFold, or SWISS-MODEL when that is scientifically
appropriate. These are external modeling decisions. PolyzyMD should receive the
curated result; do not ask PolyzyMD or this preparation script to infer missing
heavy atoms or missing residues automatically.

After curation, preserve the PolyzyMD conventions from this tutorial: all protein
atom records use chain `A`, and TER records remain between disconnected protein
fragments.

## Use the validated PDB in PolyzyMD

After the OpenFF validation snippet succeeds, point the enzyme block in your
PolyzyMD config at the prepared PDB:

```yaml
enzyme:
  name: "alpha_chymotrypsin_4cha"
  pdb_path: "structures/4cha_chymotrypsin_chain_a_openff.pdb"
```

Then run the normal PolyzyMD config and build checks:

```bash
pixi run -e build polyzymd validate -c config.yaml
pixi run -e build polyzymd build -c config.yaml --dry-run
```

```{warning}
These commands check PolyzyMD configuration and planned build inputs. They do not
prove that OpenFF can assign enzyme chemistry from the PDB. A config validation
and dry-run build can pass even when a real build later fails while loading the
enzyme. Use the OpenFF validation snippet above, or a real PolyzyMD build/load,
to confirm OpenFF chemistry validity.
```

## Why not just run `polyzymd clean-pdb`?

`polyzymd clean-pdb` is a convenience helper for simple PDB cleanup. It is not a
substitute for selecting the biological system, deciding which crystal copy to
simulate, removing unwanted molecules, preserving fragment boundaries, assigning
PolyzyMD chain IDs, or modeling missing residues and heavy atoms. Use it only
when those scientific decisions have already been made and the input is already
close to simulation-ready.
