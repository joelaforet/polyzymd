# Why PolyzyMD Guards the Prepared Structure

This page explains three safeguards that sit between system building and
production: the build-time solute/solvent clash assertion, frozen-solute energy
minimization, and the trajectory lineage check. They exist because of a
specific failure that went undetected for a full simulation campaign, and
understanding that failure explains why the defaults are what they are.

## The failure that motivated them

Before PolyzyMD 1.2.1, `solvate_with_packmol()` packed solvent around a copy of
the solute centered in the rectangular PACKMOL brick but assembled the final
topology with the solute at the center of the triclinic cell. For a rhombic
dodecahedron the two centers differ by tens of ångströms, so several hundred
water molecules ended up inside the protein while PACKMOL reported success.

Nothing downstream complained. Energy minimization followed the gradient,
which meant pushing protein atoms away from the trapped waters: the protein
inflated (a lipase's Cα radius of gyration grew from 18.0 to 20.4 Å) before
heating even started, and the heavy-atom position restraints of the heating
stage then held the inflated structure in place. Density, volume, and potential
energy all looked normal. The defect was found only by measuring water-protein
distances in the built PDB months later.

## Safeguard 1: the clash assertion at build time

`solvate_with_packmol()` and `pack_polymers()` now measure the distance from
every packed atom to the nearest solute atom immediately after assembly. If
more than a small number of packed atoms (20, `SOLVATION_CLASH_ATOM_LIMIT`)
lie within half the PACKMOL tolerance of the solute, the build aborts with
`SolvationClashError` before any file is written.

Why a count limit rather than zero tolerance? PACKMOL frequently exits with
code 173 ("imperfect packing") for dense polymer shells and leaves a handful
of atoms slightly inside the tolerance. That is a packing-quality residual that
minimization resolves and PolyzyMD has always accepted. A frame mismatch, by
contrast, produces hundreds to thousands of overlapping atoms. The limit sits
two orders of magnitude from both cases, so it separates them without false
alarms.

## Safeguard 2: frozen-solute minimization

Even a single overlapping water deforms the protein locally during
minimization, because the minimizer moves whatever the gradient points at
regardless of mass. The protocol PolyzyMD implements is therefore: the prepared
structure enters equilibration exactly as built, restraints are applied, and
the configured heating and free-equilibration stages take it from there.

Mechanically, `SimulationRunner.minimize()` runs the minimizer on a temporary
copy of the OpenMM System in which every protein and substrate atom (the
`solute` atom group, hydrogens included) has zero mass. OpenMM never moves a
massless particle, but it also refuses constraints that involve one, so the
copy drops the constraints touching frozen atoms (they remain trivially
satisfied). Solvent and polymers relax against the fixed solute; the relaxed
coordinates are copied back into the real System, whose masses and constraints
are untouched. The runner checks that the solute displacement is zero and
records it in `minimization/phase.json`.

Set `simulation_phases.minimization.freeze_solute: false` to recover the old
behaviour. The setting is runtime-only and does not change the build manifest
hash.

## Safeguard 3: trajectory lineage

A self-resubmitting SLURM chain that is duplicated (two jobs continuing the
same replicate) writes two sets of production segments into one directory.
Concatenating them yields a trajectory that jumps backwards in time. The
analysis loader now reads each segment's raw timestamps and refuses to build a
`ChainReader` unless segment *k* starts exactly one frame interval after
segment *k−1* ends and all segments share the interval. Every segment record
also carries the PolyzyMD version, OpenMM version, and pixi environment that
produced it, so a chain that silently switched software can be identified.

## Related pages

- Configuration keys: {doc}`../reference/configuration`
- Output files and provenance: {doc}`../reference/data_requirements`
- Polymer placement: {doc}`../how_to/polymers`
- Equilibration setup: {doc}`../how_to/equilibration`
