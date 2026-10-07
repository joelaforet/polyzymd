# GROMACS and OpenMM: how one config runs on each engine

One `config.yaml` runs on either engine (`engine: openmm` or `engine: gromacs`).
This page gives, for each setting, what each engine runs. A row marked
**1:1** runs the same physics on both engines. A row marked **approximate**
samples the same ensemble with a different algorithm, and PolyzyMD writes a
warning when a GROMACS build uses it. A row marked **differs** does not run
the same way on both engines, and PolyzyMD writes a warning.

## Dynamics

| Setting | OpenMM | GROMACS (`.mdp`) | Match |
|---|---|---|---|
| `thermostat: LangevinMiddle` | `LangevinMiddleIntegrator(T, 1/τ, dt)` | `integrator = sd`, `tau-t = τ`, `ref-t = T`, `tcoupl = no` | 1:1. Both integrate the Langevin equation with friction 1/τ. The splitting schemes differ (BAOAB in OpenMM, leap-frog in GROMACS), at order dt². |
| `thermostat: Langevin` | `LangevinIntegrator(T, 1/τ, dt)` | as above | 1:1, as above |
| `thermostat: NoseHoover` | not supported: runs `LangevinMiddle`, with a warning | `integrator = md`, `tcoupl = nose-hoover` | differs |
| `thermostat: Andersen` | not supported: runs `LangevinMiddle`, with a warning | `integrator = md-vv`, `tcoupl = andersen` | differs |
| `thermostat_timescale` τ (ps) | friction = 1/τ | `tau-t` = τ | 1:1 |
| `time_step` (fs) | integrator step | `dt` (ps) | 1:1 |
| `barostat: MC` | `MonteCarloBarostat(P, T, barostat_frequency)` | `pcoupl = c-rescale`, `pcoupltype = isotropic`, `tau-p = 5 ps`, compressibility 4.5e-5 bar⁻¹ | approximate. GROMACS has no Monte Carlo barostat. Stochastic cell rescaling samples the same NPT ensemble. |
| `barostat: MCA` | not supported: `validate` refuses it | `pcoupl = c-rescale`, `pcoupltype = anisotropic` | differs: runs only on GROMACS |
| `pressure` (atm) | atm | `ref-p` in bar (× 1.01325) | 1:1 |

## Seeds

| What | OpenMM | GROMACS | Match |
|---|---|---|---|
| Starting structure | Packmol and polymer draws seeded with the replicate number | the same, built by PolyzyMD before export | 1:1 |
| Initial velocities | `setVelocitiesToTemperature(T, seed)`, seed from the replicate number | `gen_seed`, seed from the replicate number | 1:1 |
| Thermostat noise | `integrator.setRandomNumberSeed(seed)`, one seed for each stage and segment | `ld_seed`, one seed for each stage | 1:1 |
| Barostat moves | `MonteCarloBarostat.setRandomNumberSeed(seed)`, one seed for each NPT stage and segment | the c-rescale noise uses `ld_seed` | 1:1 |

The replicate number seeds the starting structure, the initial velocities,
the thermostat noise and the barostat moves. Each equilibration stage and production segment
gets its own seed, derived from the replicate number and the stage name
(`polyzymd.simulation.seeds.dynamics_seed`). A single `ld_seed` for all
stages would repeat the same noise in every stage, because GROMACS starts
the step count again at each stage.

The two engines use different random number generators, so one seed does
not give the same trajectory on OpenMM and on GROMACS. A replicate run
again does not give the same trajectory either, even on one engine and
machine. OpenMM's CPU threads, PME and GPUs add up forces in a different
order on each run, so two runs of a replicate differ from the first
minimization on. They agree statistically, not frame by frame. OpenMM has a
`DeterministicForces` platform property, and GROMACS has
`gmx mdrun -reprod`, for bitwise reruns; PolyzyMD sets neither.

## Nonbonded interactions and constraints

PolyzyMD writes the settings of the OpenFF force fields into the `.mdp`, so
that GROMACS computes the same energies as OpenMM through Interchange.

| What | Value on both engines |
|---|---|
| Electrostatics | PME, real-space cutoff 0.9 nm, Ewald tolerance 1e-5, Fourier spacing 0.12 nm |
| Van der Waals | cutoff 0.9 nm, potential switch from 0.8 nm, dispersion correction for energy and pressure |
| Constraints | bonds to hydrogen (`h-bonds`, LINCS in GROMACS); rigid water |

## Equilibration

| What | OpenMM | GROMACS | Match |
|---|---|---|---|
| Minimization with `freeze_solute: true` (the default) | protein and substrate heavy atoms held fixed | `em.mdp` sets `define = -DPOSRES_EM`: position restraints of 1e5 kJ/mol/nm² on the protein and ligand heavy atoms, whatever the equilibration stages restrain | approximate. The same atoms are held; in GROMACS they can move by about 0.01 nm. |
| Position restraints | a harmonic force on the selected atoms | `define = -DPOSRES_...` with restraint `.itp` files | 1:1 |
| Temperature ramp | the target changes every `temperature_interval_steps` steps | `annealing`, a piecewise-linear schedule that changes the target over one step at each boundary | 1:1: the same target at every step |

## Output

| What | OpenMM | GROMACS |
|---|---|---|
| Trajectory | DCD | XTC |
| Analysis topology | `system.prmtop`, written from the built system | the run's `.tpr`, or the `.top` when MDAnalysis cannot read the `.tpr` |

## What `grompp` may warn about

PolyzyMD passes no `-maxwarn` to `gmx grompp`, so every warning stops the
run. A neutral system built by PolyzyMD gives no warning. To accept a warning
that you have read and understood, set `grompp_flags: "-maxwarn 1"` in the
`gromacs:` section of the config.
