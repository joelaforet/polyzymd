# Why PolyzyMD Guards the Prepared Structure

This page explains four safeguards that sit between system building and
production: the build-time solute/solvent clash assertion, periodic-image
safety with a deterministic box, frozen-solute energy minimization, and the
trajectory lineage check. They exist because of specific failures that went
undetected for a full simulation campaign, and understanding those failures
explains why the defaults are what they are.

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

## Safeguard 2: periodic-image safety and a deterministic box

### The second failure: chains packed outside their own cell

Polymer packing and solvation used to size two different boxes. Chains were
packed into a rectangular box of `solute bbox + 2 * packing.padding`, and only
afterwards was the periodic cell derived from the bounding box of whatever
PACKMOL had produced, plus `solvent.box.padding`. Two things went wrong.

First, **periodic-image overlap**. A rhombic-dodecahedron cell is represented
as a rectangular brick whose `z` height is `sqrt(2)/2` (0.707) times the padded
extent. The rectangular region the chains had been packed into was taller than
that, so chains protruded through the `z` faces of the brick and landed on top
of themselves across the `c` lattice vector. PACKMOL never sees this: it is run
without periodicity and only enforces the tolerance inside its own box. In the
audited production builds this left 11 to 169 atom pairs closer than 1.5 Å to a
periodic image, down to 0.10 Å (one CALB replicate: 17 pairs below 2 Å, 11
below 1.5 Å, 2.1 % of polymer atoms outside the brick). A singular overlap like
that is not something minimization can fix; the runs died with NaN.

Second, **a non-deterministic box**. Deriving the cell from the packed
coordinates makes the cell a function of the PACKMOL seed. One RML replicate
came out with a box 23 % smaller than its siblings and was solvated with 18,865
waters instead of ~25,000, and 56/47 ions instead of 72/63. Replicates of one
condition were therefore not replicates of the same system: "38 pentamers"
meant a different polymer concentration in each one.

### The fix: one box, computed first

The periodic cell is now computed **before** anything is packed, from the
protein and substrate alone:

```
box vectors = shape_matrix @ diag(solute bbox + 2 * (packing.padding
                                                     + solvent.box.padding))
```

Chains are packed inside the rectangular brick of *that* cell, shrunk by the
PACKMOL tolerance, plus an `inside sphere` constraint centred on the solute so
they stay in a shell around the protein rather than in the brick's corners.
Solvation then reuses the same cell and, crucially, does **not** re-centre the
packed topology: a centre-of-geometry shift of an already-framed system is a
rigid translation that pushes atoms back out through the brick faces (measured:
0.02 % to 0.38 % of polymer atoms).

The guarantee is a one-line argument. If two atoms both lie inside the brick
shrunk by `tolerance`, then along the axis of any lattice vector their
separation after a translation of `±L` is at least `L - (L - tolerance)
= tolerance`. Every lattice vector of a reduced-form cell has a non-zero
component along `x`, `y` or `z` that equals a full brick edge, so no image pair
can be closer than the tolerance. Packing in the cell the system will actually
be simulated in is therefore not merely better — it is sufficient.

Because a deterministic cell is a function of the enzyme and substrate only,
replicates of a condition now share their box volume, water count and ion
count. Controls (no polymers) are unaffected: their box is still computed from
the solute at solvation time, so pre-built control bundles remain valid.

### The check: `PeriodicImageClashError`

Arguments are not evidence, so the build measures it. After polymer packing and
again after solvation, every atom is translated by each of the 26 non-zero
lattice vectors and queried against a KD-tree of the untranslated coordinates.
An atom within half the tolerance of an image aborts the build with
`PeriodicImageClashError`, naming the count, the minimum image distance, the
worst pair and the lattice vector; contacts between half the tolerance and the
tolerance only warn. The zero translation is skipped, so an atom is never
compared with itself — only with its images.

## Safeguard 3: frozen-solute minimization

Even a single overlapping water deforms the protein locally during
minimization, because the minimizer moves whatever the gradient points at
regardless of mass. The protocol PolyzyMD implements is therefore: the prepared
structure enters equilibration exactly as built, restraints are applied, and
the configured heating and free-equilibration stages take it from there.

Mechanically, `SimulationRunner.minimize()` runs the minimizer on a temporary
copy of the OpenMM System in which every protein and substrate **heavy** atom
(the `solute_heavy` atom group) has zero mass. OpenMM never moves a massless
particle, so those coordinates come back bit-identical. Solvent and polymers
relax against the fixed solute; the relaxed coordinates are copied back into
the real System, whose masses and constraints are untouched. The runner checks
that the heavy-atom displacement is zero and records it in
`minimization/phase.json`.

### Why the hydrogens stay mobile

Freezing the *whole* solute, hydrogens included, looks like the stricter choice
and is in fact a bug — one this section previously described as the design.
With `constraints=HBonds` the force field emits no harmonic bond term for an
X–H bond; it emits a constraint, and the bond length is whatever that constraint
says. A massless hydrogen keeps the length it had in the input PDB, and those
lengths come from whatever built the structure (a crystallographic hydrogen
placement, PDBFixer, a previous force field) rather than from the force field
about to be simulated. Nothing in a frozen minimization corrects them, because
OpenMM refuses to build a Context for a constraint that involves a massless
particle, so every X–H constraint was simply dropped from the copy.

The result was measurable. In the 2026-09 rebuild, all 1994 protein X–H
constraints of the RML systems were violated in the minimized state — mean
0.113 Å, maximum 0.234 Å, 950 of them by more than 0.1 Å — while the polymer
and water constraints, whose atoms were mobile, were exact to 0.000 Å.
Equilibration then sets those positions on a constrained integrator, assigns
velocities at 60 K, and steps: CCMA has to remove ~2000 violations at once and
occasionally diverges on the very first step, which OpenMM reports as
`Particle coordinate is NaN`. Whether a replicate survived was a matter of the
velocity seed — protein-identical replicates of the same condition failed and
succeeded at random (RML 0:100 replicate 1 ran, replicates 2 and 5 died; the
RML control lost replicate 1 and kept the other four).

Freezing only the heavy atoms fixes this at the source. The hydrogens are the
only solute atoms whose geometry the force field, not the crystal structure,
determines; letting the minimizer place them costs nothing scientifically and
is what the pre-`freeze_solute` behaviour did. Since OpenMM will not accept a
constraint on a massless particle in *any* position — verified empirically:
constructing such a Context raises `A constraint cannot involve a massless
particle` — the copy replaces each frozen-heavy-atom/mobile-hydrogen constraint
with a stiff harmonic bond at the constraint distance
(5 × 10⁵ kJ mol⁻¹ nm⁻²). Dropping the constraint outright is not an option: it
would leave the hydrogen with no bonded term at all and the minimizer would
scatter it. Constraints whose two particles are both frozen are removed, since
neither atom moves; in a real system there are none, because every X–H
constraint has a hydrogen on one end.

### The final check: the constraints of the real System

Because the surrogate bond is an approximation to the constraint, the runner
does not assume the result. After minimization — frozen or not — it measures
the largest `|distance − constraint length|` over every constraint of the real
System and logs it. If that exceeds 0.01 Å it calls
`Context.applyConstraints(1e-6)` to project the coordinates back onto the
constraint manifold and logs the residual. `applyConstraints` distributes a
correction by inverse mass, so a heavy partner can move by roughly 1/12 of it;
the heavy-atom displacement is therefore re-checked afterwards against a 0.02 Å
ceiling instead of the 1e-4 Å ceiling that applies to the minimization itself.
In practice the surrogate leaves a residual well below 0.01 Å, the projection
never runs, and the heavy atoms are bit-identical.

The maximum distance a solute hydrogen travelled is logged and recorded in
`minimization/phase.json` as `hydrogen_max_displacement_angstrom`. A large value
there is not an error — it is the measure of how far the input hydrogens were
from the force field's own bond lengths.

Set `simulation_phases.minimization.freeze_solute: false` to recover the old
behaviour. The setting is runtime-only and does not change the build manifest
hash.

## Safeguard 4: trajectory lineage

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
