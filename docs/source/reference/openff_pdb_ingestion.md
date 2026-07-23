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

Fixed PolyzyMD builds retain all enzyme molecules in the original OpenFF order
before substrate and polymer components. Chain assignment still follows the
canonical convention: every protein/enzyme molecule is chain `A`, substrate is
chain `B`, polymers are chain `C`, and solvent/ions start at chain `D`.
Generated PDB/GRO outputs renumber protein residues continuously across retained
enzyme molecules. For example, a homodimer whose two source monomers both use
residues 1-99 is written as chain `A` residues 1-198. Older builds incorrectly
used only `molecule(0)`, which made homodimer outputs look like single-monomer
systems even though OpenFF had loaded multiple protein molecules.

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

For strict GLYCAM N-glycosylation, PolyzyMD may write a Pablo-only product PDB
using attachment-scoped internal residue aliases. This avoids repeated-residue
and branched-residue ambiguity while Pablo parses the full modified topology.
After Pablo succeeds, PolyzyMD restores canonical GLYCAM residue names and atom
names on the topology metadata before native OpenMM GLYCAM parameterization.

## PDB-fragment glycan ingestion contract

Residue-resolved glycan moieties supplied through `moiety.input_path` are treated
as provider-neutral fragments before product-state Pablo/OpenFF ingestion. The
fragment contract is intentionally strict:

| Requirement | Expected state | Failure mode |
|---|---|---|
| Connectivity source | Complete `CONECT` records | Missing `CONECT` raises `coordinate inference is disabled` |
| Coordinate inference | Not performed | Coordinates never create or repair bonds |
| Atom serials | Every bond endpoint must be a known atom serial | Unknown serials fail before motif detection |
| Self bonds | Disallowed | Self-referential `CONECT` records fail explicitly |
| Components | One connected fragment graph | Isolated atoms or disconnected components fail |
| Explicit H degree | Exactly one graph bond per hydrogen | Degree 0 or >1 fails explicitly |
| Common valence | No obvious impossible upper valence | Overbonded H, O, C, N, S, or P fails |
| Bond orders | Assigned only on the accepted graph from atoms, explicit H, charges, and `CONECT` pairs | Assignment that changes atom count, connectivity, radicals, or total formal charge fails |

N-glycosylation motif detection is graph-first. The candidate anomeric carbon
must be bonded to an O-H leaving group and to a retained ring oxygen connected
into the sugar graph. Atom names and residue names are selectors and diagnostics;
they are not used to infer chemistry. A single motif is resolved automatically.
Zero motifs fail with `No glycan anomeric motif was found`. Multiple motifs fail
with `Ambiguous glycan anomeric motif assignments` and should be disambiguated
with the existing `moiety.link_site` reactive atom selector.

For the local G42666 CONECT acceptance fixture, the graph-first detector resolves
anomeric carbon serial 4, leaving oxygen/hydrogen serials 1 and 2, and retained
ring oxygen serial 14. The N-glycosylation transformation removes one selected
Asn `ND2` hydrogen, keeps the glycan reactive carbon, removes exactly the O+H
leaving moiety, and adds one `ND2`-C1 bond. When formal charges are available in
the participating records, the product-state charge bridge is expected to
preserve the net formal charge after leaving-group removal and bond formation.
The same strict graph supplies bond-order assignment for the acetamide carbonyls,
including source serials 6-8 and 33-35 as order 2 bonds.

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
| `Topology.from_pdb()` reports `n_molecules > 1`, but older output contains only one protein copy | Historical PolyzyMD solute combination dropped all enzyme molecules except `molecule(0)` | Compare `Topology.from_pdb(path).n_molecules` and output atom counts/chains | Use a fixed PolyzyMD build; all enzyme molecules are retained before substrate/polymers | Multiple protein molecules intentionally share chain `A`; generated residue IDs are continuous across molecules |
| `declared leaving atoms {'HZ3', 'HZ2'} not found in any LYX residue` | Pablo was configured with reactant-state leaving atom names for an already modified product PDB | Inspect the emitted product PDB and the resolved crosslink requirement; confirm the linked residues contain `LYX:NZ` and `NHX:<acyl carbon>` and no leaving atoms | Generate the Pablo crosslink from the resolved product-state plan and use `leaving_atoms: [[], []]` after PolyzyMD has already removed leaving atoms | Do not hardcode `HZ2`/`HZ3`; actual source hydrogens may have different names such as `H11`/`H13`, and product-state ingestion should not remove them again |
| `Atom A:LYX23.NZ ... has 1 radical electrons, formal charge +1, and 3 bonds` | Product-state `LYX` definition kept protonated lysine `NZ+` charge after extra NZ hydrogens were removed and the amide crosslink was added | Inspect the generated `LYX` atom definition and OpenFF atom-level radical diagnostic | Neutralize product-state `LYX:NZ` in the generated Pablo definition when leaving hydrogens are absent | The molecule-level error may mention S/P-block radical support even when the atom-level dump identifies nitrogen |
| `PDB fragment ingestion requires complete CONECT records; coordinate inference is disabled` | A residue-resolved fragment PDB lacks explicit graph connectivity | Count `CONECT` records and inspect source serial coverage | Add curated complete `CONECT` records upstream | PolyzyMD will not infer bonds from coordinates |
| `PDB fragment CONECT graph connectivity was accepted, but bond orders could not be assigned` | The graph is explicit but the PDB lacks enough chemically consistent information for radical-free bond-order assignment | Inspect explicit hydrogens, formal charges, and ambiguous valences; compare against an SDF/OpenFF source | Provide an SDF/OpenFF-native fragment source or curate the PDB graph and charges | Connectivity is not repaired or altered during this step |
| `Ambiguous glycan anomeric motif assignments` | More than one graph-valid reducing-end motif exists in a glycan fragment | Inspect candidate serials in the error message | Configure `moiety.link_site` to select the intended anomeric carbon | Do not rely on atom names such as `C1` to break ties |
| Product-state failure involving repeated glycan endpoint atoms such as `POU`/`PIN`, or repeated residue names in a branched glycan | Pablo residue matching could not distinguish repeated canonical residue definitions or chain-C endpoint roles | Inspect `assembled_crosslinked.pablo_scoped.pdb`, attachment `product_residue_mappings`, and the canonical-to-scoped alias metadata | Current PolyzyMD writes attachment-scoped aliases and restores canonical names after Pablo. If the failure persists, curate distinct residue/atom identities or reduce ambiguity in the input glycan graph | Do not export or simulate the Pablo-scoped names; they are internal parsing identities only |
| `Product-state Pablo glycan graph for attachment-local residue ... has degree ...` or `... has ... glycan neighbors ... plus a reserved protein crosslink` | The branched glycan exceeds the topology degree Pablo 0.2.2 can represent for a residue definition | Locate the attachment-local residue key and listed atoms in the glycan PDB and product mappings | Move the branch away from the reducing residue, split/curate the glycan, or wait for broader Pablo representability | Strict GLYCAM fails closed; PolyzyMD does not flatten the graph or silently route the glycan through Sage |
| `Product-state Pablo definitions are ambiguous: residue ... has the same non-leaving atom-name selector ... but different chemistry` | Same residue name and retained atom-name selector occurs with different leaving/linking chemistry | Compare generated product-state definitions and repeated residue instances | Use distinct residue names or atom names for chemically different product-state templates, or rely on PolyzyMD scoped aliases when the chemistry is identical | This guards against lower-level Pablo matching failures and prevents accidental residue collapse |
| `Product-state scoped Pablo aliasing could not find emitted residue ...` | Attachment residue mapping points to a residue not present in the emitted product PDB | Inspect `product_residue_mappings`, product PDB residue numbers, insertion codes, and attachment endpoint provenance | Fix the source glycan `CONECT`, site selection, or product assembly so each mapped residue is emitted exactly once | This is a mapping/provenance error, not a force-field fallback opportunity |

## Catalog maintenance rule

When a new OpenFF PDB ingestion error is diagnosed, update this table and
{doc}`../how_to/troubleshoot_openff_pdb_ingestion` with the exact error text,
likely cause, diagnostic command, acceptable fix, and caveats before closing the
task unless the user explicitly defers the durable documentation update.
