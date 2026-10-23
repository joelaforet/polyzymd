# Why PolyzyMD audits the physics of every build

A simulation is only as good as the system it integrates. A config says which
force field, charge model, water model, salt and restraints the system has.
Before 1.3.0rc5, nothing checked that the built system had them. PolyzyMD now
checks every build against its config, records the result, and refuses to run
files that no passing check vouches for.

This page explains the failures behind the checks, how the checks work, and
what they do not check. {doc}`../reference/physics_checks` lists every check.

## The failures that nothing reported

Each of these happened in a real PolyzyMD campaign. Each built, ran and gave
plausible trajectories, and no step reported an error.

**The protein ran with the wrong force field.** From January 2026 to the
1.3.0 acceptance tests, PolyzyMD loaded ff14SB and Sage into one OpenFF
`ForceField`. In SMIRNOFF the last matching parameter wins, so Sage's generic
parameters replaced ff14SB's on every protein bond, angle, torsion, improper
and Lennard-Jones term. Only the ff14SB charges survived. Every protein in one 18-condition study
was simulated this way. In CALB, for example, the peptide
C–N bond had the Sage length of 1.382 Å instead of ff14SB's 1.335 Å. The
config said ff14SB. Energies, densities and RMSDs looked normal.

**A production run recomputed charges.** A LipA production run charged its
substrate again at run time instead of reading the built system. The
difference was below 1e-7 e, but the run integrated a system no build record
described.

**A restraint held the wrong atom.** A restraint selected its atom with
`pdbindex`, which counts atoms in the file. The documentation said it was the
PDB serial number, which skips one number at each `TER` record. In two
enzymes the restraint held the substrate's carbonyl oxygen, not its carbonyl
carbon.

**The salt changed meaning.** `nacl_concentration` once counted the
neutralizing ions toward the target, and later did not. Builds from before and
after the change hold different numbers of ions for the same config, and no
build record said which meaning it used.

These failures share one property. Each step did what its code said, and the
result was a valid system. Only an independent comparison of the built system
with the config could find them.

## The config states every scientific choice

The audit can only compare the system with what the config says. A default
that PolyzyMD fills in is a choice nobody made. Therefore PolyzyMD has no
scientific defaults:

- every solute names its `force_field` files and its `charges`;
- `charges: {method: nagl}` names its NAGL `model` file;
- `solvent` names its `force_field` and its water `model`;
- each component has a `type`, and the type is never guessed.

The choice of force field is the user's. PolyzyMD recommends none. It checks
that the system carries what the config names, and it records the files, their
hashes and the model versions.

## How the audit checks a build

The audit runs at the end of every build, before any file is committed. It
reads the system that the engine will run, not an intermediate:

- OpenMM: `system.xml`, with the built PDB;
- GROMACS: the `.top`, `.itp` and `.gro` files, read with OpenMM's
  `GromacsTopFile` with the `define`s of the production `.mdp`, and every
  `.mdp` file.

A force that the audit cannot read is an error, not a pass. GROMACS settings
that OpenMM cannot reproduce, such as reaction-field electrostatics, are
refused with the setting named.

### Each component is rebuilt on its own

The builder records each component: its atoms, its force field files and its
charge source. The audit then rebuilds each component's molecules from those
inputs alone, with OpenFF directly. It never calls PolyzyMD's builders, so a
mistake in the builder cannot repeat in the check. It compares the built
system with the rebuild term by term: every bond, angle, torsion, improper,
Lennard-Jones and 1-4 parameter, every constraint, every charge and every
mass. A water box has thousands of copies of one molecule; the audit compares
each copy, so one wrong copy fails the check.

The force-field order bug fails `parameters.bonded` on the protein. The test
suite reproduces that bug and checks that it fails.

A coverage check makes sure that every atom belongs to exactly one component,
and that no bonded term joins two molecules.

### The system is checked against the config

The audit recomputes from the built system what the config asked for:

- the number of chains, ligands and co-solvent molecules;
- the net charge, the salt and the ionic strength, from the realized solvent
  mass, with a formula written separately from the builder's;
- the constraint of every bond to hydrogen;
- the cutoff, switching distance, PME and dispersion correction, against the
  force field file;
- clashes between solvent, polymer and protein under periodic boundaries,
  below half the PACKMOL tolerance the build used;
- the distance from the solute to its own periodic image, and polymer chains
  that touch both the protein and its image;
- each restraint's form, distance and force constant, and the atoms it joins.
  The audit selects each restrained atom again with MDAnalysis, a different
  library than the builder's, inside the component the config names;
- for GROMACS, every `.mdp` file: integrator, thermostat, time step, length,
  temperature or annealing ramp, pressure coupling, constraints, velocity
  generation and the position-restraint switches of each stage.

On the Paper 1 systems (about 90,000 atoms with 78 polymer chains), the whole
audit takes about one minute on a cluster node.

## When a check fails

A check that fails stops the build. For OpenMM, nothing is committed: there is
no `system.xml` to run. The full record is kept as
`build_provenance.failed.json` and `build_provenance.failed.md`. The error
lists each failing check with its component and detail, and names the two ways
forward: correct the config or the inputs, or accept the departure.

Some departures are deliberate. Crystal waters inside a protein fail
`geometry.buried_water` by design. To accept one, add an entry under
`modeling_exceptions` with the check, the component and the reason:

```yaml
modeling_exceptions:
  - check: geometry.buried_water
    component: LipA
    reason: "Crystal waters in the active-site pocket are kept on purpose."
```

The reason is for the reader of the published config. The build record keeps
the accepted finding with its reason and the atoms it covered, so it shows
both what PolyzyMD would have refused and why the user chose it. An entry
that accepts nothing is reported, so stale entries do not stay in a config.

## The record and the gates

The build writes `build_provenance.json` and a readable
`build_provenance.md`. They hold the versions, seeds, config hash, the SHA-256
of every file the engine reads, the component records and the full audit. See
{doc}`../reference/data_requirements`.

The record is what `run` and `submit` check. A replicate that has not started
runs only if its record exists, has a passing audit, is for the same engine,
matches the current config, and lists files whose hashes have not changed.
Each refusal says whether to rebuild or to restore an edited file. A replicate
that has already started continues from its own files, so a running campaign
is not interrupted by a new release.

The checks continue after the build:

- A GROMACS build that passes writes `audited_files.sha256`. Every GROMACS job
  checks it with `sha256sum` before `grompp`, so a file edited after
  submission does not run.
- Each OpenMM production segment compares the System it is about to integrate
  with the audited `system.xml`: every force-field term, the restraints and the
  nonbonded settings. The result is written beside the segment and into
  `progress.json`. A difference is logged as an error and shown by
  `polyzymd status` as `SYSTEM DIVERGED FROM THE AUDITED BUILD`. The run goes
  on, so that a person can decide.

## What the audit does not check

The audit checks that the built system is the one the config describes. It
does not check that the config describes good science:

- whether the force field, charge model or water model suits the system;
- the protonation states and the preparation of the input structure;
- whether the simulation is long enough, or converged;
- how closely GROMACS and OpenMM sample the same ensemble where their
  algorithms differ (see {doc}`../reference/gromacs_openmm`).

A wrong config is built faithfully, and recorded.

## Why this matters when an agent runs the simulations

An agent that sets up and runs simulations can make the same mistakes a
person makes, faster and at a larger scale. The audit gives an agent and the
person who supervises it the same guarantee: what ran is what the config says,
or the run did not start. Each refusal names the check, the component and the
fix, so an agent can act on it. A deliberate departure must be written into
the config with a reason, where a reviewer sees it. `polyzymd status
--format agent` reports a run that diverged from its build.

## Related pages

- {doc}`../reference/physics_checks` — every check, its status and component
- {doc}`../reference/data_requirements` — the build record and run files
- {doc}`../reference/configuration` — solutes, charges and `modeling_exceptions`
- {doc}`../how_to/rewrite_config_rc5` — rewrite a config for 1.3.0rc5
- {doc}`simulation_safeguards` — the build and run safeguards that came before
