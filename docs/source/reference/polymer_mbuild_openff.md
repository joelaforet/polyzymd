# mBuild to OpenFF Polymer Adapter Reference

This reference describes the mBuild -> OpenFF conversion boundary used by
PolyzyMD polymer and conjugation workflows. It is intentionally scoped to
current behavior: atomistic molecule conversion, explicit-fragment assembly
support, and adapter boundaries, not a general polymer-design user interface.

## Current status

| Feature | Current support |
|---------|-----------------|
| mBuild to OpenFF molecule conversion | `polyzymd.builders.conjugation.polymer.from_mbuild(compound)` converts an atomistic `mbuild.Compound` into a canonical OpenFF `Molecule`. |
| Existing conjugation adapters | OpenFF molecule and SDF adapters can produce `GeneratedPolymerFragment` objects for existing conjugation paths. |
| Dynamic generation default | Bundled-default dynamic generation uses native linear methacrylate generation and OpenFF charging/cache. Custom `.rxn` workflows retain deprecated Polymerist compatibility behavior. |
| Automatic chemistry suggestions | Not implemented. PolyzyMD does not infer missing chemistry or suggest atoms automatically. |
| Native fragment assembly | `generation_mode: "fragments"` supports linear explicit-fragment assembly from terminal/middle strings with mBuild `Port` and `force_overlap()`. Runtime CGSmiles and a fragment-authoring notebook are not implemented. |

## Public adapter entry points

Use the package-level import for the mBuild conversion boundary. This minimal
example builds valid atomistic H2 with an explicit single bond order:

```python
import mbuild as mb

from polyzymd.builders.conjugation.polymer import from_mbuild

mbuild_compound = mb.Compound(name="hydrogen")
hydrogen1 = mb.Particle(name="H1", element="H", pos=[0.0, 0.0, 0.0])
hydrogen2 = mb.Particle(name="H2", element="H", pos=[0.074, 0.0, 0.0])
mbuild_compound.add([hydrogen1, hydrogen2])
mbuild_compound.add_bond((hydrogen1, hydrogen2))
mbuild_compound.bond_graph[hydrogen1][hydrogen2]["bond_order"] = 1.0

openff_molecule = from_mbuild(mbuild_compound)
```

The OpenFF/SDF adapters are available for existing generated-fragment and
conjugation paths:

```python
from polyzymd.builders.conjugation.polymer import (
    generated_fragment_from_openff_molecule,
    generated_fragment_from_openff_sdf,
)
```

These adapters bridge validated molecules into `GeneratedPolymerFragment`; they
do not replace mechanism-specific conjugation validation, charge auditing, or
the existing build workflow described in {doc}`protein_modification_config`.

## Required molecule information

The mBuild compound must be atomistic. The conversion preserves:

- elements;
- explicit hydrogens;
- coordinates, with mBuild positions interpreted as nanometers;
- bond orders;
- formal charges;
- partial charges when every particle carries a charge;
- atom names;
- residue metadata where mBuild parent compounds provide it.

Explicit bond orders are required by default. If mBuild stores topology-only
bonds as missing or zero-order edges, `from_mbuild()` raises an error rather
than inventing chemistry. A fallback such as `unspecified_bond_order=1.0` must be
an explicit, deliberate caller choice and should be documented in the workflow
that uses it.

The adapter preserves the numeric mBuild coordinate values as nanometer
coordinates on the OpenFF molecule. Downstream OpenFF or RDKit displays may show
the same conformer in angstrom units.

## Validation evidence for the ACB adapter slice

The ACB adapter slice is the fixed SBMA-NHS-EGPMA sequence with labels A-C-B. It
has been checked against Polymerist for:

- formula `C31H42N2O12S`;
- 88 atoms;
- 89 bonds;
- graph isomorphism;
- matching OpenFF parameter assignments;
- matching NAGL charges;
- matching OpenMM energy.

This evidence validates the stated ACB adapter slice. It does not make arbitrary
mBuild compounds, arbitrary reaction chemistry, or incomplete bond-order inputs
production-ready.

## Recipe-defined chemistry

For the ACB recipe, interpret the chemistry as recipe-defined:

- terminal SBMA and EGPMA residues are unsaturated one-port methacrylate caps;
- the middle NHS residue is a saturated two-port residue;
- this is not documented as generic ATRP termination.

Do not generalize this recipe into automatic polymerization behavior or automatic
chemistry suggestions. New recipes need their own molecule construction, atom
identity, bond order, charge, and conjugation validation evidence.

## Relationship to configuration workflows

Current `config.yaml` workflows can use native assembly in two ways:

- `generation_mode: "dynamic"` with bundled `"default"` reactions uses the native
  linear methacrylate backend and OpenFF charging/cache.
- `generation_mode: "fragments"` uses explicit `terminal` and `middle` fragment
  strings keyed by monomer label. Terminal strings contain exactly one `*` atom;
  middle strings contain exactly two `*` atoms, with optional `[*:1]`/`[*:2]`
  direction markers.

Configuration schema validation checks that the `fragments` mapping exists and
has one `terminal`/`middle` pair for each monomer label. Runtime fragment parsing
then enforces dummy atom counts, dummy degrees, directional maps, RDKit parsing,
embedding, sanitization, mBuild port placement, and OpenFF conversion/charging.

For nonlinear, branched, or dendrimer molecules, supply explicit charged SDFs
through `polymers.provided_molecules`. PolyzyMD does not provide runtime
CGSmiles authoring or a fragment-design notebook in this phase.

For the current validated conjugation boundary and build artifacts, see
{doc}`conjugation_support_matrix` and {doc}`../how_to/validate_conjugates`.
