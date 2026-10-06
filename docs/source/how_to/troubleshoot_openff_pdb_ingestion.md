# Troubleshoot OpenFF PDB ingestion

Use this guide when a PolyzyMD build fails while OpenFF loads an enzyme PDB, or
when direct `openff.toolkit.Topology.from_pdb()` validation fails. For a first
structure that prepares without problems, see
{doc}`../tutorials/prepare_pdb_for_openff`.

```{important}
Do not bypass a charge mismatch or monkeypatch OpenFF to continue. OpenFF is
reporting that the PDB chemistry it inferred does not match a supported residue
graph. Fix or curate the structure first.
```

## Quick triage

1. Reproduce the failure outside PolyzyMD.
2. Run simple structural checks on the PDB.
3. Read the OpenFF error dump for the residue names, atom names, and charges it
   expected versus found.
4. Fix the PDB before you build.
5. If you find an error that this page does not list, open an issue on
   GitHub with the full OpenFF message.

## Check the structure first

Use a small script or text inspection to answer these questions. Some checks are
OpenFF chemistry checks; chain-ID checks are PolyzyMD input conventions that keep
protein, substrate, polymer, and solvent roles unambiguous later in the build.

- Do protein atom records use PolyzyMD chain `A`?
- Are substrate, polymer, and solvent chains kept out of the enzyme PDB or placed
  on the expected chains `B`, `C`, and `D+` when relevant?
- Are TER records present between disconnected protein fragments?
- Are crystallographic waters and unrelated heterogens removed from the enzyme
  PDB unless intentionally retained elsewhere?
- Are hydrogens explicit?
- Do PDB header records report missing residues or missing heavy atoms?
- Do disulfide cysteines have SG-SG connectivity and no SG-bound HG proton?

Example quick check:

```python
import argparse
from pathlib import Path


def summarize_pdb(path: str) -> None:
    lines = Path(path).read_text().splitlines()
    atoms = [line for line in lines if line.startswith(("ATOM", "HETATM"))]
    chains = sorted({line[21] for line in atoms})
    residues = sorted({line[17:20].strip() for line in atoms})
    elements = {line[76:78].strip() for line in atoms if len(line) >= 78}
    print(f"chains: {chains}")
    print(f"TER records: {sum(line.startswith('TER') for line in lines)}")
    print(f"hydrogens present: {'H' in elements}")
    print(f"residue names include CYX: {'CYX' in residues}")
    print(f"SSBOND records: {sum(line.startswith('SSBOND') for line in lines)}")
    print(f"CONECT records: {sum(line.startswith('CONECT') for line in lines)}")


def main() -> None:
    parser = argparse.ArgumentParser(description="Summarize a PDB before OpenFF validation")
    parser.add_argument("pdb", help="Prepared enzyme PDB to inspect")
    args = parser.parse_args()
    summarize_pdb(args.pdb)


if __name__ == "__main__":
    main()
```

## Validate directly with OpenFF

Run direct validation before editing PolyzyMD configuration:

```python
import argparse


def main() -> None:
    parser = argparse.ArgumentParser(description="Validate PDB ingestion with OpenFF")
    parser.add_argument("pdb", help="Prepared enzyme PDB to validate")
    args = parser.parse_args()

    from openff.toolkit import Topology

    topology = Topology.from_pdb(args.pdb)
    print("OpenFF PDB ingestion succeeded")
    print(f"molecules: {topology.n_molecules}")


if __name__ == "__main__":
    main()
```

With pixi:

```bash
pixi run -e build python validate_openff.py prepared_enzyme.pdb
```

If direct validation fails, the problem is in the PDB chemistry OpenFF sees. A
passing `polyzymd validate` or `polyzymd build --dry-run` does not prove OpenFF
can ingest the enzyme PDB.

## Check multi-molecule enzyme inputs

Some valid enzyme PDBs contain more than one disconnected protein molecule.
For example, a homodimer loads as two OpenFF molecules
(`Topology.from_pdb(...).n_molecules` is 2). The validation script above
prints the molecule count.

PolyzyMD keeps every protein molecule, in OpenFF order, and puts all of them
on chain `A`. The substrate is chain `B`, the polymers are chain `C`, and the
solvent starts at chain `D`. The built PDB and GRO files number the protein
residues continuously across the molecules. A homodimer of two monomers
numbered 1-99 becomes chain `A` residues 1-198. See
{doc}`../explanation/residue_assignment`.

## Interpret OpenFF error dumps

OpenFF PDB errors often include a residue-level description of atoms, bonds, and
charges. Focus on:

- the first residue name and number mentioned in the mismatch
- unexpected hydrogens on terminal atoms or cysteine SG atoms
- residues OpenFF matched as a different protonation or connectivity state
- formal-charge totals that differ between the PDB graph and the template
- nearby TER records, SSBOND records, and CONECT records

Do not delete atoms just to make the message disappear. Make a chemically
consistent model and validate it again.

## Common signatures

| Error signature | Likely cause | Diagnostic | Acceptable fix | Caveats |
|---|---|---|---|---|
| `Molecule has more/fewer total formal charges than the matched substructure` | Residue graph, hydrogens, termini, or disulfide state differs from OpenFF's matched template | Inspect the named residue's hydrogens, bonds, and neighboring TER/SSBOND/CONECT records | Correct protonation/connectivity or use a reviewed custom substructure proof of concept | Never ignore the mismatch |
| Failure around `CYS#0001`, terminal `H`, or N-terminal cysteine | N-terminal cysteine/cystine has terminal hydrogens and disulfide state OpenFF does not match cleanly | Check N-terminal atom names, SG-HG absence, and SG-SG bond | Curate the terminal cystine or use a narrow `NCYX` custom substructure proof of concept | Private OpenFF API; not universal |
| Residue names include `CYX`, but OpenFF still fails | `CYX` may be treated as a cysteine alias, not a complete public template solution | Compare residue atoms and SG-SG connectivity | Add/verify disulfide connectivity and hydrogens; consider upstream issue/PR | Do not assume renaming to CYX is sufficient |
| Failure near residues reported in `REMARK 465` | Missing-coordinate residues or missing heavy atoms affect chemistry or termini | Read PDB header and visualize gaps | Model missing regions with an external tool when scientifically appropriate | PolyzyMD should receive a curated result |

## Disulfides

For each disulfide:

1. Confirm the paired cysteine SG atoms are close and intentionally bonded.
2. Remove inappropriate SG-bound `HG` protons from disulfide cysteines.
3. Preserve or add reliable connectivity records. `SSBOND` is useful metadata;
   `CONECT` can make the actual SG-SG bond explicit for parser paths that use it.
4. Validate again with `Topology.from_pdb()`.

## Termini and hydrogen naming

Terminal residues combine residue chemistry with chain-fragment state. A mature
protein can have multiple TER-separated fragments, each with termini. Verify that
terminal hydrogens and atom names match the intended protonation state. This is
especially important for N-terminal cystines.

## Missing residues and heavy atoms

`REMARK 465` and related header records are not instructions to auto-fill a PDB.
They are warnings that the deposited model is incomplete. Decide whether to model
missing residues or heavy atoms with external tools such as PDBFixer, MODELLER,
SWISS-MODEL, AlphaFold-derived models, ChimeraX, or PyMOL workflows, then review
the result before passing it to PolyzyMD.

## Charge mismatch

Treat charge mismatch as a blocker. It means OpenFF's inferred molecule and the
matched substructure disagree. The fix is to make the residue graph, atom names,
bonds, hydrogens, and protonation state consistent, not to suppress the error.

## 4CHA case study

PDB entry 4CHA (alpha-chymotrypsin) shows the problems of a harder
structure. The scripts are in `examples/pdb_preparation/4cha/`.

### Read the header before you delete anything

The 4CHA header records say:

- COMPND puts one enzyme copy on chains `A`, `B` and `C`, and a second copy
  on chains `E`, `F` and `G`.
- Mature alpha-chymotrypsin is cleaved into three peptide fragments.
  Residues 14-15 and 147-148 are excised activation peptides. So TER records
  between the fragments are correct.
- REMARK 465 lists residues without coordinates, such as `GLY A 12` and
  `LEU A 13` of the first copy.

REMARK 465 records are not a list of residues to add. They show that the
crystal model is incomplete. If the missing residues matter for your study,
model them with an external tool and check the result.

### Look at the copies in PyMOL

```text
load structures/4CHA.pdb, raw4cha
hide everything, raw4cha
show cartoon, polymer.protein
show sticks, polymer.protein and resn CYS
show spheres, solvent

select first_copy, raw4cha and polymer.protein and chain A+B+C
select second_copy, raw4cha and polymer.protein and chain E+F+G
color marine, first_copy
color orange, second_copy
distance disulfides, first_copy and name SG, first_copy and name SG, 2.2
zoom first_copy
```

### Prepare one copy

`prepare_4cha.py` does these steps with PDBFixer:

1. It keeps chains `A`, `B` and `C`, the first copy.
2. It removes the waters and the other heterogens.
3. It adds hydrogens at pH 7. It does not add missing residues or heavy
   atoms.
4. It puts every protein atom on chain `A`, and keeps the TER records
   between the three fragments.

```bash
curl -L https://files.rcsb.org/download/4CHA.pdb -o structures/4CHA.pdb
python examples/pdb_preparation/4cha/prepare_4cha.py \
  structures/4CHA.pdb structures/4cha_chain_a.pdb
python examples/pdb_preparation/4cha/validate_openff.py structures/4cha_chain_a.pdb
```

### Result: the structure checks pass, and OpenFF refuses the file

The prepared file has every atom on chain `A`, three TER records, no waters,
no second copy, and hydrogens. `Topology.from_pdb()` still fails. It reports
errors at the terminal hydrogens and the disulfide of the N-terminal
cysteine `CYS#0001`, and at `SER#0011`. The structure is clean, but its
chemistry does not match the OpenFF templates. The file is not ready for
PolyzyMD.

`validate_openff.py --custom-substructures nterminal_cystine_substructure.json`
tests an `NCYX` template for the N-terminal cystine. It uses the private
OpenFF argument `_custom_substructures`. Use it only to test a fix that you
then propose to OpenFF, not for production simulations.

### Curate the structure

To make the file usable, fix its chemistry outside PolyzyMD. Check these
items:

- the atom names and the hydrogen count of each terminus;
- the protonation of each disulfide cysteine (no `HG` on SG), and the SG-SG
  bond;
- the residues without coordinates that the header lists;
- missing heavy atoms in residues that are present;
- excised residues, which stay absent, with TER records between the
  fragments.

Tools for this work are PyMOL or ChimeraX for manual editing, PDBFixer with
a review of each added atom, `pdb4amber`, and MODELLER, AlphaFold or
SWISS-MODEL to rebuild missing regions. After you curate the file, keep the
PolyzyMD conventions: all protein atoms on chain `A`, and TER records between
disconnected fragments. Then validate again with `Topology.from_pdb()`.

`polyzymd clean-pdb` only replaces nonstandard residues and adds hydrogens.
It does not select a copy, remove molecules, set chain IDs or model missing
atoms.

## Report a new error

If you find an OpenFF error that this page and
{doc}`../reference/openff_pdb_ingestion` do not list, open an issue on GitHub.
Include:

- the exact error text, or the shortest unique part of the traceback;
- the PDB entry or the file;
- the steps that you used to prepare it.
