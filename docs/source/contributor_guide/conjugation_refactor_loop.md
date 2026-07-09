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
→ conjugate relaxation and validation
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
NHS-Lys construction workflow. Internal helper modules remain importable only
from their owning module paths.

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
├── system_workflow.py       # config-driven full-system workflow
├── structure/               # PDB/structure parsing, preparation, and writing
│   ├── __init__.py
│   ├── inspection.py        # pure-Python structure compatibility inspection
│   ├── normalization.py     # clean-PDB chain normalization planning/writing
│   ├── pdb.py               # product PDB writing and CONECT handling
│   └── preparation.py       # protein canonicalization at configurable pH
├── polymer/                 # modifier/polymer recipe and fragment helpers
│   ├── recipe.py
│   ├── fragment.py
│   └── polymerist.py
├── placement.py             # Packmol/geometric placement
├── pablo/                   # Pablo/OpenFF integration boundary
│   ├── ingestion.py         # Pablo ingestion boundary
│   ├── product.py           # product-state Pablo residue definitions
│   ├── parameterization.py  # Interchange construction
│   └── residue_library.py   # Pablo residue-library policy helpers
└── relaxation/              # conjugate relaxation and validation evidence
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

The package facade should export only stable request/result models, the engine,
and the high-level build functions. Reaction discovery and lower-level helpers
are imported from their owning modules, not from the package root.

### Current public API direction

As of Phase 6, the public construction path is:

```python
from polyzymd.builders.conjugation import ConjugationEngine, build_conjugate_from_config
```

`build_conjugate_from_config()` and `ConjugationEngine.build_from_config()` are
the preferred entry points for executable NHS-lysine polymer builds. They
delegate to the known-working config-driven workflow while keeping orchestration
behind the engine boundary. `build_conjugate()` is retained as a higher-level
facade for request/config inputs, but direct molecule/topology construction is
still explicitly pending.

The public API notebook/workflow should use `build_conjugate_from_config()` or
`ConjugationEngine`. Internal helpers such as `system_workflow`, `structure`,
`relaxation`, and the `pablo` package remain direct module imports for maintainers,
but they are not re-exported from `polyzymd.builders.conjugation`.

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

### Phase status snapshot

| Phase | Status | Notes |
|-------|--------|-------|
| 0 | Complete | Current import behavior and POC-adjacent behavior are covered by focused conjugation tests. |
| 1 | Complete | `api.py`, `engine.py`, public models, and the reaction package skeleton exist. |
| 2 | Complete | Pablo/OpenFF boundaries live under the internal `pablo/` package; transitional top-level shims remain removed. |
| 3 | Superseded | Protein preparation and relaxation remain in direct modules; transitional shim packages were removed. |
| 4 | Complete, adapter-first | `NhsLysReaction` owns NHS-Lys defaults and reaction metadata, delegating to the proven implementation. |
| 5 | Complete, conservative | `ConjugationEngine` centralizes public orchestration and delegates executable config builds to the existing system workflow. Direct OpenFF fallback is not silently enabled. |
| 6 | Complete | Preserve old POC assets while keeping the package layout focused on the public API and direct internal modules. |

### Phase 0 — freeze current behavior

Inventory current public imports and lock the working POC behavior with tests.

Gate:

```bash
pixi run -e build pytest tests/test_conjugation_*.py -v
```

Record any expected skips or external-data requirements in the commit body.

### Phase 1 — add the new skeleton and facade

Add `api.py`, `engine.py`, `models.py`, and `reactions/` without moving major
behavior yet. Existing direct modules should keep working while public callers
move to the facade.

Gate:

```bash
pixi run -e build pytest tests/test_conjugation_*.py -v
pixi run -e build ruff check src/polyzymd/builders/conjugation tests/test_conjugation_*.py
pixi run -e build black src/polyzymd/builders/conjugation tests/test_conjugation_*.py --check
```

### Phase 2 — consolidate Pablo/OpenFF boundaries

Product-state Pablo library generation, ingestion, residue-library policy, and
Interchange construction live under `pablo/`. Do not reintroduce top-level shim
wrappers for the old `product_pablo.py`, `pablo_adapter.py`,
`residue_library.py`, or `parameterization.py` paths.

### Phase 3 — group structure helpers

PDB inspection, clean-PDB normalization, product PDB assembly, and protein
canonicalization live under `structure/`. Conjugate relaxation remains under
`relaxation/`. Do not reintroduce transitional shim wrappers for old top-level
structure paths.

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

Create `ConjugationEngine` in `engine.py`. Public callers should enter through
`build_conjugate_from_config()`, `build_conjugate()`, or the engine. Existing
full-system internals remain in `system_workflow.py`.

Remove legacy direct OpenFF fallback. Fail loudly when Pablo/OpenFF cannot ingest
the product.

### Phase 6 — clean POC artifacts and write docs

Keep `src/polyzymd/builders/conjugation/poc/conjugation_poc_walkthrough.ipynb`
available as the known-working reference POC until the refactored public API can
reproduce the same product-state Pablo/OpenMM minimization result. Do not delete
or overwrite this notebook during early refactor phases.

During Phase 6, keep old POC assets for reference while removing transitional
shim packages and fallback experiments from the importable package tree. The
notebook remains the reference for the proven path; new user-facing code should
point toward `build_conjugate_from_config()` or `ConjugationEngine`.

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
- `polyzymd.builders.conjugation._relaxation`

Do not delete existing config fields in the first refactor pass. Add new fields
only when the clean engine consumes them.

## Review checklist

Before each commit, verify:

- heavy imports remain lazy (`openmm`, `openff`, `MDAnalysis`, `pdbfixer`, `rdkit`)
- generic engine code does not hardcode NHS-Lys atom names
- residue identity is preserved; no whole-polymer `POLY` shortcut
- legacy direct OpenFF fallback is absent; Pablo/OpenFF failures surface clearly
- Pablo/OpenFF failures are surfaced with actionable diagnostics
- notebook/tutorial code uses public APIs once the API exists
- tests cover the public API path and direct internal module paths still in use

## Known current scientific boundary

The current production path uses product-state Pablo ingestion, production charge
templates, OpenFF parameterization, and `relax_conjugate()` for restrained
conjugate relaxation. It does not yet claim production-quality charge assignment
for all possible modifiers; new chemistries still need explicit validation of
their product-state templates, charge provenance, and relaxation evidence.
