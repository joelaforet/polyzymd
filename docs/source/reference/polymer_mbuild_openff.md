# mBuild to OpenFF Polymer Adapter Reference

This reference describes the approved Phase 12 mBuild -> OpenFF workflow for
polymer conjugation inputs. It is intentionally scoped to what exists now: a
conversion and adapter boundary for atomistic polymer fragments, not a new
polymer-design user interface.

## Current status

| Feature | Current support |
|---------|-----------------|
| mBuild to OpenFF molecule conversion | `polyzymd.builders.conjugation.polymer.from_mbuild(compound)` converts an atomistic `mbuild.Compound` into a canonical OpenFF `Molecule`. |
| Existing conjugation adapters | OpenFF molecule and SDF adapters can produce `GeneratedPolymerFragment` objects for existing conjugation paths. |
| Dynamic generation default | The existing config YAML user experience is unchanged; Polymerist remains the current default backend for dynamic generation in this phase. |
| Automatic chemistry suggestions | Not implemented. PolyzyMD does not infer missing chemistry or suggest atoms automatically. |
| Native polymerization mechanism | Follow-on work. Native polymerization and the fragment-designer notebook are not implemented in this phase. |

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

## Validation evidence for the approved slice

The Phase 12 parity target is the ACB polymer slice: the fixed
SBMA-NHS-EGPMA sequence with labels A-C-B. It has been checked against
Polymerist for:

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

For the approved ACB recipe, interpret the chemistry as recipe-defined:

- terminal SBMA and EGPMA residues are unsaturated one-port methacrylate caps;
- the middle NHS residue is a saturated two-port residue;
- this is not documented as generic ATRP termination.

Do not generalize this recipe into automatic polymerization behavior. New
recipes need their own molecule construction, atom identity, bond order, charge,
and conjugation validation evidence.

## Relationship to configuration workflows

Existing `config.yaml` workflows do not need new keys for this phase. Users who
run dynamic polymer generation through the current YAML interface continue to use
the existing Polymerist-backed path unless a future workflow explicitly documents
a native mBuild/OpenFF backend.

For the current validated conjugation boundary and build artifacts, see
{doc}`conjugation_support_matrix` and {doc}`../how_to/validate_conjugates`.
