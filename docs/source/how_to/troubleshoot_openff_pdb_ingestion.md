# Troubleshoot OpenFF PDB ingestion

Use this guide when a PolyzyMD build fails while OpenFF loads an enzyme PDB, or
when direct `openff.toolkit.Topology.from_pdb()` validation fails.

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
4. Fix the PDB upstream, or use a documented narrow proof of concept only when
   the user accepts the caveat.
5. Add newly diagnosed error signatures to the durable error catalog.

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

    Topology.from_pdb(args.pdb)
    print("OpenFF PDB ingestion succeeded")


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
| `declared leaving atoms {'HZ3', 'HZ2'} not found in any LYX residue` or similar Pablo crosslink leaving-atom error | A product-state modified PDB was passed to Pablo with reactant-state leaving atom names after PolyzyMD had already removed those atoms | Inspect the emitted product PDB for exact product residue and atom names, for example `LYX:NZ` to `NHX:C047`, and confirm the leaving atoms are absent from the linked residue | Use a product-state `ccd_pablo.crosslinks` entry with exact emitted residue and atom names and `leaving_atoms: [[], []]` when PolyzyMD has already removed the leaving atoms | Do not hardcode canonical lysine names such as `HZ2`/`HZ3`; generated hydrogens may be named `H11`/`H13`, and product PDBs should not ask Pablo to remove them again |
| `Atom A:LYX23.NZ ... has 1 radical electrons, formal charge +1, and 3 bonds` followed by `does not currently support parsing molecules with S- and P-block radicals` | Product-state residue definitions carried protonated lysine `NZ+` charge into an acylated lysine with only `CE`, one `HZ`, and the amide carbon bonded | Inspect the Pablo product definition for `LYX:NZ` charge and bonds; confirm product PDB has removed the extra lysine hydrogens and contains the `LYX:NZ` to `NHX:C047` bond | Generate the product-state `LYX` definition with neutral `NZ` when the crosslinked nitrogen has been deprotonated by leaving-hydrogen removal | Do not treat this as an SBMA sulfonate radical until the atom-level OpenFF diagnostic identifies sulfur or oxygen |

## Product-state modified residues

For already modified conjugation PDBs, PolyzyMD has performed the graph surgery
before Pablo ingestion. The emitted product PDB should contain product residue
names such as `LYX` for the acylated lysine and `NHX` for the reacted
NHS-derived polymer monomer, plus the actual crosslink bond in `CONECT` or
equivalent metadata. Pablo crosslink settings must describe that product state.

Use the exact atom names present in the PDB. For the current NHS-Lys product POC,
the crosslink is `LYX:NZ` to `NHX:C047`. Because the PDB no longer contains the
reactant-side leaving atoms on the linked residues, the Pablo crosslink should use
empty leaving atom groups: `leaving_atoms: [[], []]`.

Do not reuse reactant-state leaving names such as `HZ2`/`HZ3` as a default. The
source protein may have used different hydrogen names, and those atoms should
already be absent from the product PDB.

## Residue-resolved glycan PDB fragments

PolyzyMD's PDB-fragment moiety ingestion is strict about connectivity. A glycan
fragment supplied with `moiety.input_path` must include complete `CONECT` records
for the fragment graph. PolyzyMD does not infer, guess, or repair fragment bonds
from coordinates.

After accepting the explicit graph, PolyzyMD may assign bond orders on that exact
graph from parsed atoms, explicit hydrogens, formal charges, and `CONECT` pairs.
This step must preserve the atom count and undirected bond set exactly; it never
creates or deletes connectivity.

Common fragment-ingestion messages:

- `PDB fragment ingestion requires complete CONECT records; coordinate inference
  is disabled`: add curated `CONECT` records to the source fragment.
- `PDB fragment CONECT references unknown atom serials`: a `CONECT` endpoint is
  not present as an `ATOM` or `HETATM` serial.
- `PDB fragment CONECT contains self bonds`: a source atom is bonded to itself.
- `PDB fragment graph is disconnected`: the fragment has isolated atoms or
  disconnected components after reading only `CONECT` records.
- `PDB fragment CONECT explicit hydrogens must have degree 1`: an explicit
  hydrogen has zero or multiple graph bonds.
- `PDB fragment CONECT graph has atoms above obvious upper valence`: the explicit
  graph is chemically impossible for a common element such as H, C, N, or O.
- `bond orders could not be assigned`: the PDB graph is connected, but PDB atom
  records do not carry enough chemically consistent information to assign formal
  bond orders without radicals or charge changes. Provide an SDF or OpenFF-native
  source for that fragment chemistry.
- `No glycan anomeric motif was found`: N-glycosylation could not find a graph
  motif with an anomeric carbon, an O-H leaving group, and a retained ring oxygen.
- `Ambiguous glycan anomeric motif assignments`: more than one graph motif was
  possible. Use the existing `moiety.link_site` selector to choose the reactive
  atom instead of adding name-based assumptions.

Atom and residue names are used for diagnostics and selectors, not for chemistry
inference. The reducing-end motif can be named `C1/O1/HO1`, use a separate `ROH`
cap, or use source-specific names, provided the `CONECT` graph describes the same
chemistry. For the audited G42666 CONECT fixture, the selected transformation is
anomeric carbon serial 4, leaving oxygen/hydrogen serials 1 and 2, retained ring
oxygen serial 14, removal of one Asn `ND2` hydrogen, and formation of a single
`ND2`-C1 bond.

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

The 4CHA alpha-chymotrypsin preparation can pass structural checks after selecting
one enzyme copy, relabeling protein atoms to chain `A`, preserving TER records,
removing waters, and adding hydrogens. Direct OpenFF validation may still fail
around the N-terminal cystine and nearby residues. That failure separates PDB
cleanup from chemistry assignment.

See `examples/pdb_preparation/4cha/` for proof-of-concept scripts, including an
`NCYX` custom substructure example. The `_custom_substructures` argument is a
private OpenFF API, so this example is evidence for a targeted workaround or
upstream contribution, not a production-ready universal preparation method.

## Keep the error catalog current

When you diagnose a new OpenFF PDB ingestion error, update this page and
{doc}`../reference/openff_pdb_ingestion` before marking the task complete. Add:

- exact error text or shortest unique traceback excerpt
- likely cause
- diagnostic snippet or command
- acceptable fix
- caveats, especially private APIs or structure-specific assumptions

If the user explicitly defers the update, record that deferral in your final
response.
