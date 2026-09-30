# Native contacts implementation verification

This page records how `polyzymd analyze native_contacts` computes Q and what
its values were checked against.

## What the code does

`polyzymd.analyses.functions.native_contacts`:

1. Finds, once per reference, the pairs of the selected atoms more than
   `min_separation` residues apart (by residue index in the topology) whose
   reference positions are closer than `radius` Å, with MDAnalysis
   `lib.distances.capped_distance` and no periodic images.
2. On each frame, computes the distance of every native pair with MDAnalysis
   `lib.distances.calc_bonds`, using the minimum image of the frame's box.
3. Returns MDAnalysis `analysis.contacts.soft_cut_q` of those distances and
   the reference distances: the mean over pairs of
   `1 / (1 + exp(beta * (r - lambda_constant * r0)))`.

## Agreement with independent implementations

Both checks were run in September 2026 on a 363 K *B. subtilis* lipase A
replicate with a 50:50 SBMA-EGMA polymer.

| Check | Frames | Result |
|---|---|---|
| Heavy atoms at the defaults, against the Best-Hummer-Eaton example of the MDTraj documentation (`mdtraj.compute_distances`, 0.45 nm, residue index difference above 3, β = 50 nm⁻¹, λ = 1.8), same reference frame | 72 | 2,320 native pairs; Q agrees within 1.5 × 10⁻⁷ on every frame |
| Cα atoms of residues 4 to 174 with `radius=8.0` and a prepared crystal structure as `reference_file`, against the continuous Q stored by the analysis scripts of a study of polymer-coated lipases, which use the same switching function | 96 | 457 native pairs in both; Q agrees within 5 × 10⁻⁸ on every frame |

The differences are at the precision of single-precision coordinates.

## What the tests check

`tests/analyses/test_native_contacts.py` checks, on small synthetic systems,
the native pairs at the radius and separation boundaries, Q against values
computed by hand, the region, the minimum image, and the results, settings,
figures and refusals of `polyzymd analyze native_contacts`.

## Scope

The reference must be whole, because native pairs are found without periodic
images. Residue separation counts residues in topology order, so atoms of
different chains are always far enough apart to form a native pair.

## See Also

- [Native Contacts Quick Start](../how_to/analysis_native_contacts_quickstart.md) — commands and settings
- [Shipped analysis functions](../reference/analysis_functions.md) — what each function measures
