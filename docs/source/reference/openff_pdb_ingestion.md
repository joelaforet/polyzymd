# OpenFF PDB ingestion reference

This reference separates OpenFF chemistry requirements from PolyzyMD enzyme-input
expectations and lists known error signatures for protein PDB ingestion through
`openff.toolkit.Topology.from_pdb()`.

## OpenFF chemistry requirements

OpenFF does not require PolyzyMD's chain IDs. It requires a PDB whose inferred
chemical graph can be matched to supported residue chemistry.

| Requirement | Expected state | Notes |
|---|---|---|
| Hydrogens | Explicit | OpenFF protein PDB ingestion expects chemically complete hydrogens |
| TER records | Present where fragments are disconnected | Mature cleaved proteins may have multiple fragments |
| Missing residues | Curated intentionally | Header records such as `REMARK 465` require scientific review |
| Disulfides | SG-SG connectivity clear; no SG-HG proton | Verify `SSBOND` and, when needed, `CONECT` records |
| Direct validation | `Topology.from_pdb()` succeeds | Run this before relying on PolyzyMD build steps |

## PolyzyMD enzyme-input expectations

PolyzyMD uses chain IDs to assign biological roles during system building and
analysis. These are project conventions, not OpenFF parser requirements.

| Role | PolyzyMD chain convention | Notes |
|---|---|---|
| Protein/enzyme | `A` | The enzyme PDB passed to OpenFF is usually protein-only on chain `A`; multiple OpenFF enzyme molecules are all retained on chain `A` with continuous output residue numbering |
| Substrate | `B` | Usually kept separate from the enzyme PDB and configured as substrate input |
| Polymer | `C` | Used for conjugates and polymer-specific selections |
| Solvent/ions/other | `D` and later | Usually generated or handled outside the enzyme PDB |

An enzyme PDB can satisfy PolyzyMD chain conventions and still fail OpenFF
ingestion if the residue graph, hydrogens, termini, or disulfide connectivity do
not match supported chemistry.

## Multi-molecule enzyme topology behavior

OpenFF may parse a protein-only enzyme PDB as multiple molecules when the input
contains disconnected protein copies, such as a homodimer. Diagnose this with:

```python
from openff.toolkit import Topology

topology = Topology.from_pdb("enzyme.pdb")
print(topology.n_molecules)
```

PolyzyMD keeps every enzyme molecule, in OpenFF order, before the substrate
and polymer components. Every protein molecule is chain `A`, the substrate is
chain `B`, the polymers are chain `C`, and solvent and ions start at chain
`D`. The built PDB and GRO files number the protein residues continuously
across the molecules. For example, a homodimer whose two monomers both use
residues 1-99 is written as chain `A` residues 1-198.

## OpenFF disulfide behavior

- `CYX` may be accepted as a cysteine-like residue alias during parsing, but it is
  not a stable public OpenFF residue template for every disulfide case.
- Disulfide cysteine SG atoms should be bonded to each other and should not have
  an attached `HG` proton.
- `SSBOND` records identify intended disulfides. `CONECT` records can make the
  SG-SG bond explicit for parser paths that depend on connectivity.
- N-terminal cystines combine terminal hydrogens with disulfide chemistry and can
  expose template/charge mismatches.

## Custom substructures JSON

PolyzyMD's `enzyme.custom_substructures_path` loads JSON and passes it to
`Topology.from_pdb(..., _custom_substructures=...)`.

```{warning}
`_custom_substructures` is a private/experimental OpenFF API. Treat examples as
proofs of concept or upstream-PR candidates, not as stable public OpenFF support.
```

Shape:

```json
{
  "RESNAME": {
    "[SMARTS:1]": ["ATOM1"]
  }
}
```

Each residue name maps to SMARTS patterns, and each SMARTS pattern maps to the
corresponding PDB atom names for that residue.

## Charge diagnostics

Charge mismatch messages are blockers. They usually mean one of these is wrong:

- protonation state or terminal hydrogen count
- disulfide SG-HG or SG-SG bonding
- residue atom naming
- missing heavy atoms or missing residues
- ambiguous TER, SSBOND, or CONECT records
- a custom substructure that does not match the PDB atom graph

Acceptable fixes are chemically explicit: curate the PDB, correct hydrogens and
connectivity, model missing atoms when scientifically justified, or document a
narrow custom-substructure proof of concept. Do not suppress the error.

## Running error catalog

| Exact signature | Likely cause | Diagnostic | Acceptable fix | Caveats |
|---|---|---|---|---|
| `Molecule has more/fewer total formal charges than the matched substructure` | OpenFF matched a residue graph whose formal charge differs from the PDB graph | Inspect the named residue's atom list, bonds, hydrogens, TER records, and disulfide records | Correct residue chemistry or use a reviewed custom substructure proof of concept | Do not ignore; private custom substructures are not stable API |
| Error dump names `CYS#0001`, terminal `H`, or N-terminal cysteine/cystine | N-terminal cysteine has terminal hydrogens plus disulfide chemistry that does not match OpenFF's template | Check SG-HG absence, SG-SG bond, N-terminal hydrogens, and residue naming | Curate the cystine or test a structure-specific `NCYX` custom substructure | Seen in 4CHA proof of concept; not universal |
| Renaming disulfide cysteine to `CYX` does not resolve ingestion | `CYX` aliasing is not equivalent to a complete public template for all contexts | Validate direct OpenFF ingestion and inspect charge mismatch | Fix connectivity/hydrogens or prepare an upstream OpenFF issue/PR | Avoid relying on residue rename alone |
| `Input residue A:LYS#0162 contains atoms matching substructures {'PEPTIDE_BOND', 'NO MATCH'}` at the last residue of a chain | The chain ends early, so its last residue has no terminal `OXT`; `polyzymd clean-pdb` does not add it | Look for `OXT` in the last residue | Add only `OXT` with PDBFixer, then run `clean-pdb` ({ref}`add-terminal-oxt`) | Seen in 181L, which lacks residues 163-164 |
| Failure adjacent to residues listed in `REMARK 465` | Missing-coordinate residues or missing heavy atoms alter termini or local chemistry | Read PDB header and visualize gaps | Model missing regions externally if required for the study | Automatic filling is a modeling decision |

## Report a new error

If you find an OpenFF PDB error that this table does not list, open an issue
on GitHub. Include the exact error text, the likely cause if you know it, and
the steps that you used to prepare the file.
