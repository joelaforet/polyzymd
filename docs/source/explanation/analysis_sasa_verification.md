# SASA implementation verification

The SASA analysis computes its values with MDTraj. PolyzyMD checked the values
against independent calculations, and works around one defect in MDTraj.

## What the code does

`polyzymd.analyses.functions.sasa` and `residue_sasa` do the following for each
frame:

1. Build, once per universe and context, an MDTraj topology with one atom per
   context atom, in index order, carrying its name, its residue and its element
   from the MDAnalysis universe.
2. Give the frame's context coordinates, converted from Å to nm, to
   `mdtraj.shrake_rupley(mode="atom")` in a call of their own, with the probe
   radius and the number of sphere points. MDTraj takes each atom's radius from
   its element.
3. Sum the per-atom areas of the target atoms, over the whole target for
   `sasa` or per residue for `residue_sasa`, and convert nm² to Å².

`residue_sasa` then averages each residue over the production frames.

`mdtraj.shrake_rupley` puts `n_sphere_points` points on a sphere around each
atom of the **context**. The radius of the sphere is the radius of the atom
from the MDTraj element table plus the probe radius. MDTraj counts the points
that are inside no sphere of another context atom. The SASA of an atom is the
area of its sphere times the fraction of free points. The SASA of the
**target** is the sum over the target atoms. Context atoms that are not in the
target, such as polymer atoms, cover part of the target surface but are not
counted. The calculation does not use periodic images.

The elements come from the loaded universe. PolyzyMD fills them in from the
atom types or names.

## MDTraj gives later frames of one call too much area

In MDTraj 1.11.1, `mdtraj.shrake_rupley` returns a slightly larger area for
every frame after the first one each OpenMP thread computes in a call. With
one thread, given the same coordinates several times in one call, it returns
one value for the first copy and a larger one, about 0.1 percent more, for every
later copy. With several threads, each thread's first frame is right and its
later frames are too large, so the result depends on how many frames are passed
per call and on the number of threads, which it should not.

`scripts/benchmarks/sasa_mdtraj_check.py` shows this and checks which value is
right. It computes each frame three ways:
- with `polyzymd.analyses.functions.sasa`;
- with an independent NumPy/SciPy Shrake-Rupley calculation that uses MDTraj's
  radii and its golden-section spiral of sphere points;
- with MDTraj given the same frame three times in one call.

Run in September 2026 with MDTraj 1.11.1 and one OpenMP thread on a
*B. subtilis* lipase A replicate at 363 K with SBMA-EGMA polymer (2,697 protein
atoms):

| Context | Frame | Independent | `sasa` | MDTraj, 3 copies in one call |
|---|---|---|---|---|
| protein | 249 | 8477.723 Å² | 8477.786 Å² | 8477.801, 8486.039, 8486.048 Å² |
| protein | 250 | 8623.484 Å² | 8623.260 Å² | 8623.275, 8631.693, 8631.702 Å² |
| protein | 251 | 8549.732 Å² | 8549.960 Å² | 8549.966, 8558.327, 8558.333 Å² |
| protein and polymer (10,419 atoms) | 249 | 5103.677 Å² | 5103.740 Å² | 5103.745, 5108.704, 5108.710 Å² |
| protein and polymer | 250 | 5223.710 Å² | 5223.472 Å² | 5223.478, 5228.576, 5228.582 Å² |
| protein and polymer | 251 | 5147.490 Å² | 5147.416 Å² | 5147.423, 5152.463, 5152.469 Å² |

The single-frame value agrees with the independent calculation to within
0.005 percent, and the later copies are 5 to 8 Å² too large. PolyzyMD therefore
gives MDTraj one frame per call. That costs speed, because MDTraj's parallel
loop runs over frames, but it makes every frame's area independent of the
others.

## What the tests check

`tests/analyses/test_sasa.py` checks, on small synthetic systems:
- an isolated atom's SASA is the area of its sphere, from MDTraj's radius plus
  the probe;
- overlapping atoms and context atoms outside the target reduce the target's
  area;
- `residue_sasa` over frames with identical coordinates returns the
  single-frame values exactly, which fails if frames are batched into one
  MDTraj call;
- a target outside its context and an unknown element are refused;
- the command line and `analyze` report the documented results, figures and
  settings.

## Scope

These checks cover the tested selections and trajectories. The areas depend on
MDTraj's radius table, which has no entry for some elements; such atoms are
refused rather than given a guessed radius. Periodic images are not
considered, so atoms occlude each other only within the coordinates as loaded.

## See Also

- [SASA Quick Start](../how_to/analysis_sasa_quickstart.md) — commands and settings
- [Shipped analysis functions](../reference/analysis_functions.md) — what each function measures
