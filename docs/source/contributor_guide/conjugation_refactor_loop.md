# Conjugation refactor loop plan

This page is the handoff reference for agents working on the `polyzymd.builders.conjugation`
refactor. The goal is to turn the working NHS-lysine POC into a compact,
general-purpose engine for covalently attaching molecules or polymers to proteins.

## Current objective

Keep the package name `conjugation`, but refactor the current exploratory module
into a succinct, testable workflow:

```text
protein selection + modifier molecule + reaction template
→ resolved reaction plan
→ product-state structure
→ Pablo/OpenFF topology
→ restrained local minimization
```

The working POC already proved the key scientific path for an NHS-lysine polymer
attachment. Future work should preserve that behavior while replacing exploratory
module sprawl with a clean API and reaction-template library.

## Ralph-loop operating protocol

Agents working this loop must repeat these steps until the module is complete:

1. Choose the next smallest phase from the plan below.
2. Implement only that logical unit.
3. Run the phase-specific tests and any affected lint/docs checks.
4. Commit the result before continuing.
5. Use conventional commit style with an explanatory body.

Commit format:

```text
type(scope): short imperative summary

Explain what changed, why it changed, and what validation was run.
Mention known remaining blockers or intentional deferrals.
```

Examples:

```text
refactor(conjugation): introduce engine facade

Adds the public build_conjugate entry point and delegates to the existing
NHS-Lys POC workflow through compatibility shims. Existing import paths are
unchanged.

Validation: pixi run -e build pytest tests/test_conjugation_*.py -v
```

```text
feat(conjugation): add nhs-lys reaction template

Moves NHS-Lys-specific planning out of the generic workflow and into the first
built-in ReactionTemplate implementation. The engine now asks the reaction
template to resolve site/moiety atoms and product-state Pablo definitions.

Validation: focused conjugation tests, ruff, black.
```

Do not batch unrelated cleanups into the same commit. If the working tree
contains unrelated local edits, stage only the intended files and inspect
`git diff --cached --stat` before committing.

## Target package structure

Refactor toward this structure under `src/polyzymd/builders/conjugation/`:

```text
conjugation/
├── __init__.py              # small public facade
├── api.py                   # build_conjugate(), build_conjugate_from_config()
├── engine.py                # orchestrates prepare/place/assemble/pablo/minimize
├── models.py                # public Pydantic models
├── exceptions.py
├── reactions/
│   ├── __init__.py          # get_reaction(), list_reactions()
│   ├── base.py              # ReactionTemplate protocol/base class
│   ├── library.py           # built-in and future external discovery
│   ├── nhs_lys.py           # first working reaction template
│   └── explicit_linkage.py  # simple explicit-link reaction
├── workflow/
│   ├── preparation.py       # protein canonicalization at configurable pH
│   ├── modifier.py          # modifier/polymer loading and generation
│   ├── placement.py         # Packmol/geometric placement
│   ├── assembly.py          # product PDB writing and CONECT/LINK handling
│   ├── pablo.py             # product-state Pablo residue definitions/ingestion
│   ├── parameterization.py  # Interchange construction
│   ├── minimization.py      # restrained OpenMM minimization
│   └── system.py            # config-driven full system workflow
├── io/
│   ├── pdb.py               # PDB atom records, selectors, normalization
│   ├── rdkit.py             # RDKit helpers; lazy imports only
│   ├── openff.py            # OpenFF/Pablo adapter boundaries
│   └── polymerist.py        # Polymerist import/generation adapters
└── _legacy/
    └── ...                  # temporary compatibility wrappers only
```

## Public API target

The stable user-facing API should stay small:

```python
from polyzymd.builders.conjugation import build_conjugate

result = build_conjugate(
    protein_pdb="enzyme.pdb",
    modifier="polymer.sdf",
    reaction="nhs_lys",
    site={"chain_id": "A", "residue_number": 23, "atom_name": "NZ"},
    output_dir="conjugate",
    protein_prep={"ph": 7.0},
)
```

Also support config-driven execution:

```python
from polyzymd.builders.conjugation import build_conjugate_from_config

result = build_conjugate_from_config("config.yaml", output_dir="build/")
```

The public facade should export only stable request/result models, reaction
discovery helpers, and the high-level build functions. Keep exploratory helpers
private or under `_legacy`.

## Core data model

Use one central workflow object rather than many unrelated result types:

```text
ResolvedReactionPlan
```

It should carry:

- selected protein atom identities
- selected modifier atom identities
- leaving atoms
- new/removed/changed bonds
- product residue names
- product-state Pablo definitions or requirements
- product PDB path
- topology/interchange/minimization diagnostics

The generic engine should consume this object; reaction-specific logic should
live inside `ReactionTemplate` implementations.

## Reaction library design

New reaction mechanisms should be contributor-friendly. Define a base class in
`reactions/base.py`:

```python
class ReactionTemplate:
    name: ClassVar[str]
    aliases: ClassVar[tuple[str, ...]] = ()
    Settings: ClassVar[type[BaseModel]]

    def resolve_plan(self, ctx: ReactionPlanningContext) -> ResolvedReactionPlan: ...
    def build_product_pablo_definitions(self, ctx: ProductPabloContext) -> ProductStatePabloLibrary: ...
    def validate_product(self, ctx: ProductValidationContext) -> ConjugationDiagnostics: ...
```

The first built-ins should be:

- `nhs_lys` / `nhs_lys_amide`
- `explicit_linkage`

The generic engine must not hardcode NHS/Lys atom names. NHS-specific atom
selection, leaving-group detection, formal-charge fixes, and geometry checks
belong in `reactions/nhs_lys.py`.

## Migration phases

### Phase 0 — freeze current behavior

Inventory current public imports and lock the working POC behavior with tests.

Gate:

```bash
pixi run -e build pytest tests/test_conjugation_*.py -v
```

Record any expected skips or external-data requirements in the commit body.

### Phase 1 — add the new skeleton and facade

Add `api.py`, `engine.py`, `models.py`, `reactions/`, `workflow/`, and `io/`
without moving major behavior yet. Existing modules should keep working through
wrappers or shims.

Gate:

```bash
pixi run -e build pytest tests/test_conjugation_*.py -v
pixi run -e build ruff check src/polyzymd/builders/conjugation tests/test_conjugation_*.py
pixi run -e build black src/polyzymd/builders/conjugation tests/test_conjugation_*.py --check
```

### Phase 2 — move preparation and minimization

Move current working protein canonicalization and restrained minimization into:

```text
workflow/preparation.py
workflow/minimization.py
```

Keep compatibility imports from `protein_preparation.py` and
`local_minimization.py` for now.

### Phase 3 — move Pablo product-state logic

Move product-state Pablo library generation and ingestion into:

```text
workflow/pablo.py
io/openff.py
```

Compatibility wrappers should remain for `product_pablo.py`, `pablo_adapter.py`,
and related tests until callers migrate.

### Phase 4 — introduce `NhsLysReaction`

Move NHS-Lys-specific logic into:

```text
reactions/nhs_lys.py
```

The engine should call `get_reaction("nhs_lys")` and delegate reaction-specific
planning/product definitions to that template.

Gate: reviewers must confirm that generic workflow modules no longer hardcode
NHS/Lys atom names such as `NZ`, `C047`, `O020`, `H11`, `H13`, `HZ2`, or `HZ3`.

### Phase 5 — consolidate orchestration

Create `ConjugationEngine` in `engine.py`. Existing functions like
`build_conjugated_polymer_system_from_config()` should delegate to the engine or
remain as thin compatibility wrappers.

Disable legacy direct OpenFF fallback by default. Fail loudly when Pablo/OpenFF
cannot ingest the product unless the caller explicitly opts into diagnostics.

### Phase 6 — clean POC artifacts and write docs

Keep `src/polyzymd/builders/conjugation/poc/conjugation_poc_walkthrough.ipynb`
available as the known-working reference POC until the refactored public API can
reproduce the same product-state Pablo/OpenMM minimization result. Do not delete
or overwrite this notebook during early refactor phases.

Once the clean API fully reproduces the notebook workflow, copy or move the
walkthrough into docs/tutorial space and leave a short pointer from the old
location if needed. Remove generated outputs and old exploratory artifacts from
the package tree only after the tutorial copy and API tests are in place.

Add:

- `docs/source/how_to/build_nhs_lys_conjugate.md`
- `docs/source/reference/conjugation_reactions.md`
- `docs/source/contributor_guide/conjugation_reactions.md`

Run a clean docs build after toctree changes:

```bash
pixi run -e build make -C docs clean html
```

## Compatibility rules

Maintain these import paths temporarily with wrappers/deprecation notes:

- `polyzymd.builders.conjugation.system_workflow`
- `polyzymd.builders.conjugation.product_pablo`
- `polyzymd.builders.conjugation.pdb_assembly`
- `polyzymd.builders.conjugation.local_minimization`
- `polyzymd.builders.conjugation.protein_preparation`

Do not delete existing config fields in the first refactor pass. Add new fields
only when the clean engine consumes them.

## Review checklist

Before each commit, verify:

- heavy imports remain lazy (`openmm`, `openff`, `MDAnalysis`, `pdbfixer`, `rdkit`)
- generic engine code does not hardcode NHS-Lys atom names
- residue identity is preserved; no whole-polymer `POLY` shortcut
- `allow_direct_openff_fallback` or equivalent legacy fallback is not silently enabled
- Pablo/OpenFF failures are surfaced with actionable diagnostics
- notebook/tutorial code uses public APIs once the API exists
- tests cover both compatibility imports and the new API path

## Known current scientific boundary

The current working POC reaches local constrained minimization for the NHS-Lys
polymer case. It does not yet claim production-quality charge assignment for all
possible modifiers. Keep the local-minimization smoke path separate from the
production parameterization/charge-template story until validated.
