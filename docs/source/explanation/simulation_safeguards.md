# Why PolyzyMD guards the prepared structure

PolyzyMD checks the built system before production and checks the
trajectory before analysis. Each check prevents a failure that a real
simulation campaign suffered without any error message. The checks are:

- the solute clash check at build time;
- the periodic-image check, with a box that is fixed before packing;
- energy minimization with the solute heavy atoms frozen;
- the trajectory lineage check;
- the rules that keep a restart chain recoverable and stoppable.

## The failure that the build checks prevent

PACKMOL packs solvent around the solute in a rectangular box. If the final
topology puts the solute at a different position than the one PACKMOL
packed around, the solvent overlaps the protein. For a rhombic dodecahedron,
the two centers can differ by tens of ångströms. In one PolyzyMD build,
several hundred water molecules ended up inside the protein, and PACKMOL
reported success.

Nothing after the build detected the error. Energy minimization pushed
protein atoms away from the trapped waters. The protein inflated before
heating started: a lipase's Cα radius of gyration grew from 18.0 to 20.4 Å.
The heavy-atom restraints of the heating stage then held the inflated
structure in place. Density, volume and potential energy all looked normal.
Only a measurement of water-protein distances in the built PDB found the
defect, months later.

## The clash check at build time

`solvate_with_packmol()` and `pack_polymers()` measure the distance from
each packed atom to the nearest solute atom, directly after assembly. The
build stops with `SolvationClashError` before it writes any file if more than
20 packed atoms (`SOLVATION_CLASH_ATOM_LIMIT`) are closer to the solute than
half the {term}`PACKMOL` tolerance.

The limit is a count, not zero. PACKMOL often exits with code 173
("imperfect packing") for dense polymer shells. It then leaves a few atoms
slightly inside the tolerance, and minimization removes these small overlaps.
A frame mismatch gives hundreds to thousands of overlapping atoms. A limit of
20 is far from both cases, so it separates them.

## The periodic-image check and a fixed box

### The failure: chains packed outside their own cell

A rhombic-dodecahedron cell is stored as a rectangular brick. The `z` height
of the brick is `sqrt(2)/2` (0.707) times the padded extent. If polymer
chains are packed into a rectangular box taller than the brick, they extend
through the `z` faces. Across the `c` lattice vector they then overlap their
own periodic images. PACKMOL cannot see this, because it runs without
periodic boundaries.

In audited production builds, this left 11 to 169 atom pairs closer than
1.5 Å to a periodic image, down to 0.10 Å. One CALB replicate had 17 pairs
below 2 Å, and 2.1 % of its polymer atoms were outside the brick.
Minimization cannot fix an overlap this close. The simulations failed with
NaN coordinates.

A box computed from the packed coordinates has a second fault: its size
depends on the PACKMOL seed. One RML replicate got a box 23 % smaller than
the other replicates. It held 18,865 waters instead of about 25,000, and
56/47 ions instead of 72/63. "38 pentamers" therefore meant a different
polymer concentration in each replicate.

### What PolyzyMD does

PolyzyMD computes the periodic cell before it packs anything. It uses the
protein and the substrate only:

```
edge = solute diameter + 2 * (packing.padding + solvent.box.padding)
box vectors = edge * shape_matrix
```

Every lattice vector is at least one edge long, so the solute starts at least
`2 * padding` from each of its periodic copies. The edge grows further if the
solute's bounding box would come closer than the PACKMOL tolerance to a face of
the brick. PolyzyMD centres the bounding box, not the centre of mass, in the
brick, so opposite faces get the same clearance. A solute that pokes through a
face, or comes within the tolerance of it, overlaps the solvent that PACKMOL
places next to its periodic image. Earlier versions made the rhombic
dodecahedron from `bbox + 2 * padding` per axis. Its brick was too short along
`z` for an elongated protein such as T4 lysozyme (0.12 nm clearance at 1.2 nm
padding), and centring the centre of mass left one face of ubiquitin only
0.06 nm from the solute.

Then it does these steps:

1. It packs the polymer chains inside the rectangular brick of that cell,
   shrunk by the PACKMOL tolerance.
2. It adds an `inside sphere` constraint centered on the solute. The chains
   then stay in a shell around the protein, not in the corners of the brick.
3. It solvates in the same cell, and it does not center the packed system
   again. A second centering moves all atoms rigidly and pushes some back out
   through the brick faces. Measured, this moved 0.02 % to 0.38 % of polymer
   atoms outside.

Packing inside the shrunk brick is enough to prevent image overlaps. Take two
atoms inside the brick shrunk by `tolerance`. Along the axis of a lattice
vector, a translation by `±L` leaves them at least
`L - (L - tolerance) = tolerance` apart. Each lattice vector of a
reduced-form cell has one component along `x`, `y` or `z` equal to a full
brick edge. Therefore no image pair is closer than the tolerance.

The cell depends only on the enzyme and the substrate. The replicates of a
condition therefore have the same box volume, water count and ion count. A
system without polymer gets its box from the solute at solvation time.

### The check: `PeriodicImageClashError`

The build also measures the result. After polymer packing, and again after
solvation, it translates every atom by each of the 26 non-zero lattice
vectors. It compares the translated atoms with a KD-tree of the original
coordinates. An atom closer than half the tolerance to an image stops the
build with `PeriodicImageClashError`. The message gives the count, the
smallest image distance, the worst pair and the lattice vector. A contact
between half the tolerance and the full tolerance gives a warning only.

## Minimization with frozen solute heavy atoms

The minimizer moves any atom that the energy gradient points at. One
overlapping water can therefore deform the protein. PolyzyMD keeps the
prepared protein structure unchanged through minimization. The heating and
equilibration stages then start from the structure as built.

`SimulationRunner.minimize()` does these steps:

1. It makes a temporary copy of the OpenMM System.
2. In the copy, it sets the mass of every protein and substrate heavy atom
   (the `solute_heavy` atom group) to zero. OpenMM never moves a massless
   particle.
3. It minimizes. Solvent, polymer and solute hydrogens relax against the
   fixed heavy atoms.
4. It copies the coordinates back into the real System. The masses and
   constraints of the real System do not change.
5. It checks that the heavy atoms did not move, and records the result in
   `minimization/phase.json`.

To minimize without frozen atoms, set
`simulation_phases.minimization.freeze_solute: false`. The setting applies at
run time and does not change the build manifest hash.

### Why the hydrogens can move

With `constraints=HBonds`, the force field has no harmonic bond term for an
X–H bond. It has a constraint. OpenMM refuses to build a Context in which a
constraint involves a massless particle. If the hydrogens were frozen too,
every X–H constraint would have to be removed from the copy. Each hydrogen
would then keep its bond length from the input PDB, which came from the
crystal structure or from another tool, not from this force field.

This was measured. With frozen hydrogens, all 1994 protein X–H constraints
of the RML systems were violated after minimization: mean 0.113 Å, maximum
0.234 Å, 950 of them by more than 0.1 Å. The polymer and water constraints,
whose atoms could move, were exact. The equilibration integrator then had to
correct about 2000 violations in its first step. The constraint solver
(CCMA) sometimes failed, and OpenMM reported `Particle coordinate is NaN`.
Some replicates survived and others died, at random.

The force field, not the crystal structure, decides the hydrogen positions.
Therefore PolyzyMD lets the minimizer place them. In the copy, it replaces
each constraint between a frozen heavy atom and a mobile hydrogen with a
stiff harmonic bond at the constraint distance (5 × 10⁵ kJ mol⁻¹ nm⁻²).
Without that bond, the hydrogen would have no bonded term, and the minimizer
would move it away. A constraint between two frozen atoms is removed, because
neither atom moves. A real protein has none, because each X–H constraint has
a hydrogen at one end.

### The constraint check after minimization

The harmonic bond only approximates the constraint. Therefore the runner
measures the result. After minimization, frozen or not, it does these steps:

1. It computes the largest `|distance − constraint length|` over all
   constraints of the real System, and logs it.
2. If this value is more than 0.01 Å, it calls
   `Context.applyConstraints(1e-6)` and logs the remaining violation.
3. It checks the heavy-atom displacement again. `applyConstraints` moves
   both atoms of a constraint, in inverse proportion to their mass, so the
   limit at this step is 0.02 Å. During minimization itself the limit is
   1e-4 Å.

In practice the violation stays well below 0.01 Å, step 2 does not run, and
the heavy atoms do not move at all.

The runner records the largest hydrogen displacement as
`hydrogen_max_displacement_angstrom` in `minimization/phase.json`. A large
value is not an error. It shows how far the input hydrogens were from the
force field's bond lengths.

## The trajectory lineage check

A self-resubmitting SLURM chain can be duplicated: two jobs continue the same
replicate. They then write two sets of production segments into one folder.
Joined in order, these segments give a trajectory that jumps back in time.

The analysis loader reads the raw time stamps of each segment. Each segment
must start one frame interval after the previous segment ends, and all
segments must have the same frame interval. The loader repairs two boundary
faults of a restart:

- A segment that starts at the time of the previous segment's last frame. The
  loader leaves out the repeated frame.
- A segment that starts two frame intervals after the previous one ends. The
  loader records the missing frame.

The loader names each repair in the analysis warnings. Any other overlap or
gap, or a change of frame interval, stops the analysis with
`TrajectoryLineageError`.

Each segment record also stores the PolyzyMD version, the OpenMM version and
the pixi environment of the job that wrote it. A chain that changed software
partway is therefore visible.

## A restart chain must be recoverable and stoppable

A production simulation on a preemptable SLURM queue is a chain of jobs.
Each job runs one {term}`segment` and submits the next job. The chain
gives a usable trajectory only if each link can continue after an
interruption. A person must also be able to stop the chain on purpose.

The figure shows the chain of jobs and what each job does.

```{figure} ../_static/diagrams/run_lifecycle.svg
:alt: Each SLURM job checks for STOP, runs the next segment with run-segment, checks progress and submits itself again until production is complete; GROMACS resumes one production run from its checkpoint.
:width: 100%

An OpenMM replicate runs minimization, equilibration and production
segments; every stage and segment writes `progress.json`. Each job checks for
`STOP`, runs `run-segment`, runs `check-progress` and submits itself again
while work remains. A GROMACS job resumes one production run with `-cpi` and
`-maxh`, and centres the trajectory when production is finished.
```

**A chain must not run out of nodes.** A job pinned to a CUDA environment can
land on a node whose driver is too old. The job then submits itself again and
excludes that node. The list of excluded nodes passes to every following
job and grows with each bad node. The retry budget counts one unhealthy
stretch, not the whole simulation. When a job reaches `run-segment`, the node
works, so the next job starts with a full budget and no exclusions. When the
budget runs out, the job names every node it tried. You need this list to
submit again by hand.

**A runtime that does not match is a routing problem.** PolyzyMD pins a
replicate to one pixi environment, OpenMM build, platform and precision,
because a trajectory that mixes them is not defensible. A mismatch is a
property of the node. The job therefore moves to another node, as it does for
an old driver. The chain does not stop.

**Recovery uses portable state first.** An OpenMM `.chk` checkpoint is binary
and reloads reliably only under the build that wrote it. A serialized
`State` XML file reloads under any build. Continuation therefore loads the
newest state XML file first. It uses a checkpoint only when no state XML file
survives. In that case it compares the `openmm_version` and
`pixi_environment` of the previous segment with the running job, and logs the
decision.

**PolyzyMD deletes no partial data.** A hard-killed segment is recognized
only by a stale checkpoint and a missing marker. This is a guess. PolyzyMD
renames such a folder to `production_N.hardkilled-<timestamp>`. The renamed
folder no longer counts as progress, and all its files stay. PolyzyMD does
not make the guess when a `restart_state.xml` shows that the segment can be
recovered.

**A person can stop a chain.** `scancel` sends `SIGTERM`. To a chain built to
survive preemption, `SIGTERM` means "the time limit is near; save and
continue". So `scancel` alone makes the chain submit its next job. The signal
cannot show whether the scheduler or a person sent it. Therefore
`polyzymd cancel` writes a `STOP` marker into the replicate folder, and the
job script checks for the marker before it submits the next job. The marker
is a plain text file that says who wrote it and when. See
{doc}`../how_to/hpc_slurm`.

## Related pages

- {doc}`../reference/configuration` — configuration keys
- {doc}`../reference/data_requirements` — output files and provenance
- {doc}`../how_to/polymers` — polymer placement
- {doc}`../how_to/equilibration` — equilibration stages
