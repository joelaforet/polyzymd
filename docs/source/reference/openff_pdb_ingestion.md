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
| Protein/enzyme | `A` | The enzyme PDB passed to OpenFF is usually protein-only on chain `A` |
| Substrate | `B` | Usually kept separate from the enzyme PDB and configured as substrate input |
| Polymer | `C` | Used for conjugates and polymer-specific selections |
| Solvent/ions/other | `D` and later | Usually generated or handled outside the enzyme PDB |

An enzyme PDB can satisfy PolyzyMD chain conventions and still fail OpenFF
ingestion if the residue graph, hydrogens, termini, or disulfide connectivity do
not match supported chemistry.

## Source-protein hydrogen canonicalization

Conjugation workflows that ingest a protein through Pablo should canonicalize
source-protein hydrogens before product assembly. PolyzyMD provides
`canonicalize_protein_hydrogens()` for this purpose. The helper runs PDBFixer,
then uses OpenMM `Modeller.delete()` to remove existing hydrogens and
`Modeller.addHydrogens()` to regenerate canonical protein hydrogens at a
user-configurable pH.

The pH defaults to `7.0` and is surfaced in the conjugation POC walkthrough as
`PROTEIN_CANONICALIZATION_PH`. Chain IDs and residue numbers are written with
`keepIds=True` when OpenMM can preserve them.

## Product-state conjugation PDBs with Pablo

For post-crosslink conjugation products, PolyzyMD may have already removed
reactant leaving atoms and written the crosslink connectivity before Pablo reads
the PDB. In that product-state workflow, `ccd_pablo.crosslinks` should describe
the atoms that remain in the emitted PDB, not the reactant atoms that were
removed during assembly.

| Field | Product-state expectation | Example |
|---|---|---|
| Modified lysine residue | Acylated lysine product residue name | `LYX` |
| Reacted NHS-derived monomer | Product monomer residue name, preserving per-monomer residue identity | `NHX` |
| Linking atoms | Exact PDB atom names present in the product file | `LYX:NZ` to `NHX:C047` |
| Leaving atoms | Empty groups when PolyzyMD already removed them | `leaving_atoms: [[], []]` |
| Atom naming | Match emitted PDB names exactly; do not canonicalize unless the residue definition supports it | `C047`, not a guessed `C` |

Do not collapse the whole polymer chain into one residue to make ingestion easier.
Each monomer residue identity should remain explicit so custom residue definitions
and OpenFF parameterization can reason about the actual product graph.

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
| Failure adjacent to residues listed in `REMARK 465` | Missing-coordinate residues or missing heavy atoms alter termini or local chemistry | Read PDB header and visualize gaps | Model missing regions externally if required for the study | Automatic filling is a modeling decision |
| `declared leaving atoms {'HZ3', 'HZ2'} not found in any LYX residue` | Pablo was configured with reactant-state leaving atom names for an already modified product PDB | Inspect the emitted product PDB and the resolved crosslink requirement; confirm the linked residues contain `LYX:NZ` and `NHX:<acyl carbon>` and no leaving atoms | Generate the Pablo crosslink from the resolved product-state plan and use `leaving_atoms: [[], []]` after PolyzyMD has already removed leaving atoms | Do not hardcode `HZ2`/`HZ3`; actual source hydrogens may have different names such as `H11`/`H13`, and product-state ingestion should not remove them again |
| `Atom A:LYX23.NZ ... has 1 radical electrons, formal charge +1, and 3 bonds` | Product-state `LYX` definition kept protonated lysine `NZ+` charge after extra NZ hydrogens were removed and the amide crosslink was added | Inspect the generated `LYX` atom definition and OpenFF atom-level radical diagnostic | Neutralize product-state `LYX:NZ` in the generated Pablo definition when leaving hydrogens are absent | The molecule-level error may mention S/P-block radical support even when the atom-level dump identifies nitrogen |

## Catalog maintenance rule

When a new OpenFF PDB ingestion error is diagnosed, update this table and
{doc}`../how_to/troubleshoot_openff_pdb_ingestion` with the exact error text,
likely cause, diagnostic command, acceptable fix, and caveats before closing the
task unless the user explicitly defers the durable documentation update.
