# N-glycosylation from SMILES feature plan

This page is the durable implementation plan for the forthcoming Ralph loop that
adds N-glycosylation/conjugation from SMILES to `polyzymd.builders.conjugation`.
It is a handoff document for implementation agents, not a user tutorial.

The public conjugation API is considered stable for this work:

```python
from polyzymd.builders.conjugation import (
    ConjugateBuildRequest,
    ConjugationEngine,
    ConjugationResult,
    build_conjugate,
    build_conjugate_from_config,
)
```

Keep new work behind these entry points. Do not introduce a parallel public API
for glycosylation.

## Goal

Implement a generic moiety-conjugation path that can build an N-glycosylated
protein product from:

- a clean, unmodified protein PDB
- one or more requested attachment sites
- user-provided moiety definitions, usually SMILES plus a residue code
- a reaction mechanism name, initially `n_glycosylation`

The first production target is N-glycosylation of Asn side chains. The current
working mechanism, NHS-Lys conjugation, remains the reference implementation for
placement, product assembly, Pablo/OpenFF ingestion, and final minimization.

## Non-goals

- Do not implement an ingest-existing mode. Users provide an unmodified protein
  and a list of modifications to add.
- Do not require users to provide moiety atom indices. Mechanisms identify
  reactive atoms and leaving atoms from the moiety graph.
- Do not make PolyzyMD own SDF-to-SMILES conversion tooling. PolyzyMD core
  accepts SMILES; documentation may show users how to convert with RDKit.
- Do not support arbitrary multi-residue glycans in the first pass. A SMILES
  moiety is one residue. Polymers remain the only multi-residue modifier type.
- Do not perform multiple minimizations after each attachment unless a later
  validation result shows it is necessary. The target workflow performs one
  final minimization after all attachments are assembled.

## User config UX sketch

Keep the YAML close to the existing conjugation config style. Users identify
protein locations and moiety definitions; the mechanism resolves chemistry.

```yaml
conjugation:
  protein_pdb: src/polyzymd/builders/conjugation/poc/5fyj-monomer.pdb
  output_dir: build/n-glycosylated-5fyj

  force_field:
    protein: amber14-all.xml
    water: amber14/tip3p.xml
    small_molecule: openff_unconstrained-2.2.1.offxml  # Sage is acceptable initially

  modifications:
    - id: glycan_asn_61_chain_a
      mechanism: n_glycosylation
      protein_site:
        chain_id: A
        residue_number: 61
        residue_name: ASN
        atom_name: ND2
      moiety:
        type: smiles
        residue_name: NAG
        smiles: "OC[C@H]1O[C@@H](O)[C@H](O)[C@@H](O)[C@H]1O"

    - id: glycan_asn_122_chain_a
      mechanism: n_glycosylation
      protein_site:
        chain_id: A
        residue_number: 122
        residue_name: ASN
        atom_name: ND2
      moiety:
        type: smiles
        residue_name: NAG
        smiles: "OC[C@H]1O[C@@H](O)[C@H](O)[C@@H](O)[C@H]1O"
```

Python usage should continue through the stable facade:

```python
from polyzymd.builders.conjugation import build_conjugate_from_config

result = build_conjugate_from_config(
    "n_glycosylated_5fyj.yaml",
    output_dir="build/n-glycosylated-5fyj",
)
assert result.success
```

The exact field names may change during implementation, but preserve these UX
rules:

- user sites are protein residue/atom locations
- user moieties are named chemical inputs, not atom-index recipes
- mechanism-specific settings are optional and validated by that mechanism
- multiple attachments are first-class, not a later extension

## Architecture principles

- Keep `build_conjugate()`, `build_conjugate_from_config()`, and
  `ConjugationEngine` as the only user-facing orchestration surface.
- Put N-glycosylation chemistry in a reaction template, for example
  `reactions/n_glycosylation.py`. Generic engine code must not hardcode Asn or
  glycan atom names.
- Reuse the generic Packmol placement strategy developed for NHS-Lys.
- Build all attachments into one product structure, then run conjugate relaxation
  and validation.
- Keep heavy dependencies lazy: RDKit, OpenMM, OpenFF, Pablo, and Packmol-related
  imports belong inside functions or methods.
- Keep product-state residue naming explicit and deterministic.
- Treat SMILES moieties as single-residue fragments. Only polymer inputs may
  expand into multiple residues.

## Proposed data models and config additions

Prefer additive changes to the existing conjugation models.

```python
class ConjugateBuildRequest(BaseModel):
    protein_pdb: Path
    output_dir: Path | None = None
    modifications: list[ConjugationModification]
    # Existing NHS-Lys fields stay supported during migration.


class ConjugationModification(BaseModel):
    id: str
    mechanism: str
    protein_site: ProteinAttachmentSite
    moiety: MoietyDefinition
    settings: dict[str, Any] = Field(default_factory=dict)


class ProteinAttachmentSite(BaseModel):
    chain_id: str
    residue_number: int
    residue_name: str | None = None
    insertion_code: str | None = None
    atom_name: str | None = None


class MoietyDefinition(BaseModel):
    type: Literal["smiles", "polymer"]
    residue_name: str
    smiles: str | None = None
    polymer: PolymerDefinition | None = None
```

Reaction planning should resolve a typed product plan that is shared across
mechanisms:

```python
class ResolvedAttachmentPlan(BaseModel):
    modification_id: str
    mechanism: str
    protein_site: ProteinAtomIdentity
    moiety_residue_name: str
    product_protein_residue_name: str
    new_bonds: list[BondEdit]
    removed_bonds: list[BondEdit]
    removed_atoms: list[AtomIdentity]
    charge_adjustments: list[FormalChargeEdit] = Field(default_factory=list)
    diagnostics: list[str] = Field(default_factory=list)
```

The final `ConjugationResult` should report every attachment, all generated
product files, minimization diagnostics, and any mechanism warnings.

## Reaction template design

Use the existing reaction-template direction from the conjugation refactor plan.
The N-glycosylation implementation should look like another built-in reaction,
not a special mode in the engine.

```python
class NGlycosylationReaction(ReactionTemplate):
    name: ClassVar[str] = "n_glycosylation"
    aliases: ClassVar[tuple[str, ...]] = ("n_linked_glycan", "n_glycan")
    Settings: ClassVar[type[BaseModel]] = NGlycosylationSettings

    def resolve_plan(self, ctx: ReactionPlanningContext) -> ResolvedAttachmentPlan:
        ...

    def build_product_pablo_definitions(self, ctx: ProductPabloContext):
        ...

    def validate_product(self, ctx: ProductValidationContext):
        ...
```

The template owns:

- validation that the target protein residue is Asn-like and has `ND2`
- detection of the glycan reducing-end/anomeric carbon from SMARTS
- detection of the anomeric hydroxyl/leaving atoms
- bond edits that replace reducing-end C-O with C-N to Asn `ND2`
- product protein residue naming, `ASN -> ASX`
- product moiety residue naming from the user-provided `residue_name`

## SMILES and SDF handling

PolyzyMD core should accept SMILES for non-polymer moieties. Users who have SDF
files can convert them outside the core workflow. Documentation may include this
RDKit snippet using the validation asset `test_glycan.sdf`:

```python
from pathlib import Path

from rdkit import Chem

sdf_path = Path("src/polyzymd/builders/conjugation/poc/data/test_glycan.sdf")
supplier = Chem.SDMolSupplier(str(sdf_path), removeHs=False)
mol = supplier[0]
if mol is None:
    raise ValueError(f"Could not read {sdf_path}")

# Preserve stereochemistry because anomeric configuration matters.
smiles = Chem.MolToSmiles(mol, isomericSmiles=True)
print(smiles)
```

Implementation notes:

- RDKit may be used internally to parse and inspect SMILES, but imports must be
  lazy.
- Store the original user SMILES and the canonical/isomeric SMILES used for
  graph matching in diagnostics.
- Fail with an actionable error when stereochemistry is absent or ambiguous at
  the detected anomeric center.
- Do not expose atom-index fields as the primary user path.

## Generic moiety fragment pipeline

The fragment pipeline should be mechanism-neutral:

1. Parse the user moiety definition.
2. For `type: smiles`, create one residue with the user-provided residue code.
3. Generate a 3D conformer if coordinates are not supplied.
4. Assign atom names deterministically and record the atom-name map.
5. Ask the reaction template to identify reactive and leaving atoms.
6. Place each moiety near its target site with the generic Packmol placement
   strategy.
7. Apply all bond/atom edits to produce one product PDB and product-state residue
   definitions.
8. Parameterize the product through Pablo/OpenFF boundaries.
9. Run one final restrained minimization after all attachments.

The fragment pipeline must support multiple independent copies of the same
SMILES. Each copy needs unique chain/residue/atom identities in the assembled
product while preserving the same residue name and chemistry.

## N-glycosylation chemistry assumptions

First implementation assumptions:

- The protein target is an Asn side-chain amide nitrogen, normally atom `ND2`.
- Product protein residue names always change. For N-glycosylation, rename the
  modified protein residue from `ASN` to `ASX`.
- The glycan residue name is user-provided and remains attached to the entire
  SMILES moiety.
- The reaction creates a C-N bond between the reducing-end/anomeric carbon and
  Asn `ND2`.
- The reaction removes one Asn `ND2` hydrogen plus the glycan anomeric
  hydroxyl/leaving atoms.
- The mechanism replaces the reducing-end anomeric C-O connectivity with C-N.
- The reducing-end/anomeric carbon is detected from SMARTS.
- Branching is R-group context. Branch atoms are retained as part of the same
  residue unless the input is explicitly a polymer.
- OpenFF Sage parameters are acceptable for glycan atoms in the first
  implementation.

Do not overfit to one glycan atom ordering. The validation glycan may come from
SDF, but production detection must use graph chemistry rather than source-file
indices.

## Multi-attachment workflow

Multiple attachments are required for this feature. The engine should plan all
attachments before product assembly.

```text
request
→ validate all protein sites and moiety definitions
→ parse/generate all moiety fragments
→ resolve one attachment plan per modification
→ detect conflicts across plans
→ place all moiety fragments
→ assemble one product-state structure
→ build Pablo/OpenFF residue definitions
→ parameterize one product system
→ run one final minimization
→ return one ConjugationResult
```

Conflict checks should include:

- duplicate attachment to the same protein atom
- duplicate residue renaming plans that disagree
- moiety residue/chain identity collisions in the product PDB
- leaving atoms referenced by more than one attachment
- Packmol placement failure or severe steric overlap before minimization

## Implementation milestones

### Phase 0 — freeze API and assets

- Add focused tests that assert the public facade exports
  `build_conjugate`, `build_conjugate_from_config`, `ConjugationEngine`,
  `ConjugationResult`, and `ConjugateBuildRequest`.
- Add tests that the POC assets needed for this loop are present:
  `5fyj-monomer.pdb`, `chainF.pdb`, `chainP.pdb`, `chainX.pdb`, and
  `data/test_glycan.sdf`.
- Do not change chemistry yet.

### Phase 1 — config and model skeleton

- Add `ConjugationModification`, `ProteinAttachmentSite`, and
  `MoietyDefinition` models or equivalent fields.
- Teach config loading to accept multiple modifications while preserving the
  existing NHS-Lys workflow.
- Add validation errors for missing SMILES, missing residue names, and missing
  protein sites.

### Phase 2 — reaction registration

- Add `NGlycosylationReaction` as a built-in reaction template.
- Wire discovery through the existing reaction-library mechanism.
- Keep its implementation mostly validation-only until graph matching lands.

### Phase 3 — SMILES moiety fragments

- Parse SMILES into a single-residue moiety fragment.
- Generate deterministic atom names and atom maps.
- Add a tiny internal representation that can be consumed by Packmol placement,
  product PDB assembly, and Pablo/OpenFF parameterization.

### Phase 4 — SMARTS-based anomeric detection

- Implement reducing-end/anomeric carbon detection with SMARTS.
- Implement anomeric hydroxyl/leaving-atom detection.
- Add clear diagnostics when there are zero or multiple candidate anomeric
  centers.

### Phase 5 — product bond edits and residue naming

- Generate the N-glycosidic C-N bond plan.
- Remove one Asn `ND2` hydrogen and the glycan anomeric hydroxyl/leaving atoms.
- Rename modified protein residues to `ASX`.
- Preserve the user-provided glycan residue name for the whole moiety.

### Phase 6 — multi-attachment assembly

- Resolve and apply all attachment plans in one product assembly.
- Ensure repeated use of the same SMILES creates independent residue instances.
- Add conflict detection before writing product files.

### Phase 7 — parameterization and validation path

- Route the product through Pablo/OpenFF boundaries.
- Use OpenFF Sage for the glycan moiety in the first implementation.
- Run conjugate relaxation and validation after all attachments.
- Record relaxation and parameterization diagnostics in `ConjugationResult`.

### Phase 8 — validation and user-facing docs

- Add a how-to page for N-glycosylation from SMILES after the implementation is
  executable.
- Add a reaction reference entry for `n_glycosylation`.
- Keep this contributor plan as the implementation history and handoff record.

## Testing plan

### Unit tests

- Public API exports and request/result model validation.
- Config parsing for one and multiple `n_glycosylation` modifications.
- Rejection of ingest-existing style inputs.
- Rejection of SMILES moieties without `residue_name`.
- RDKit/SMARTS detection of exactly one anomeric center for the validation
  glycan.
- Ambiguous or missing anomeric-center diagnostics.
- Protein site validation for Asn `ND2`.
- Product residue renaming from `ASN` to `ASX`.
- Bond/atom edit plans: new C-N bond, removed C-O/leaving atoms, removed Asn
  hydrogen.
- Multi-attachment conflict detection.

### Integration tests

Use the real validation assets under
`src/polyzymd/builders/conjugation/poc/`:

- `5fyj-monomer.pdb`
- `chainF.pdb`
- `chainP.pdb`
- `chainX.pdb`
- `data/test_glycan.sdf`

Recommended scenarios:

1. Convert `data/test_glycan.sdf` to isomeric SMILES with RDKit in a test helper
   or fixture.
2. Build one N-glycosylation on `5fyj-monomer.pdb`.
3. Build multiple N-glycosylations on `5fyj-monomer.pdb`.
4. Verify the product PDB has modified Asn residue names `ASX` and user-provided
   glycan residue names.
5. Verify product connectivity contains each expected C-N linkage.
6. Smoke-test Pablo/OpenFF ingestion and final minimization in the pixi build
   environment.

Mark slow or stack-heavy tests with `@pytest.mark.slow` when they require real
OpenMM/OpenFF/RDKit execution.

## Validation commands

Run focused tests while developing:

```bash
pixi run -e build pytest tests/ -v -k "conjugation and glycosylation"
```

Run the existing conjugation tests before each Ralph-loop handoff:

```bash
pixi run -e build pytest tests/test_conjugation_*.py -v
```

Run formatting and lint checks for touched code:

```bash
pixi run -e build ruff check src/polyzymd/builders/conjugation tests/ -v
pixi run -e build black src/polyzymd/builders/conjugation tests/ --check
```

After adding user-facing pages or changing toctrees, the coordinator should run:

```bash
pixi run -e build make -C docs clean html
```

This plan page was validated with a clean docs build before commit.

## Risks and open questions

- SMARTS coverage: one anomeric-detection pattern may not cover all useful
  reducing-end sugars.
- Stereochemistry: missing or inconsistent anomeric stereochemistry should fail
  clearly rather than silently choosing a product.
- Protonation: selecting which Asn `ND2` hydrogen to remove may depend on input
  naming and local preparation.
- Residue naming: `ASX` must not collide with unrelated ambiguous Asn/Asp use in
  downstream tooling.
- Pablo residue definitions: product-state residue templates for `ASX` plus a
  user-named glycan residue may require careful bond-externalization.
- OpenFF Sage coverage: acceptable for a first implementation, but validation
  should record any problematic parameters or charge assignments.
- Packmol placement: closely spaced multiple glycans may need stronger initial
  placement constraints before final minimization.
- Branching: treating all glycan branching as same-residue R-group context may
  limit compatibility with residue-based glycan analysis tools.

## Done criteria

The Ralph loop is complete when:

- `build_conjugate()` and `build_conjugate_from_config()` can build at least one
  N-glycosylated 5FYJ product from a clean protein PDB and SMILES glycan input.
- Multiple N-glycosylation attachments work in one request.
- Users do not provide moiety atom indices.
- The mechanism detects the reducing-end/anomeric carbon and leaving atoms from
  graph chemistry.
- Modified protein residues are renamed `ASX`, and glycan residues use the
  user-provided residue code.
- The product is parameterized with the accepted first-pass OpenFF Sage path.
- One final minimization runs after all attachments.
- Unit and integration tests cover the real 5FYJ plus glycan validation system.
- User-facing documentation explains SMILES input and includes an RDKit SDF to
  isomeric SMILES snippet.
- The docs build is clean after the coordinator runs it.
