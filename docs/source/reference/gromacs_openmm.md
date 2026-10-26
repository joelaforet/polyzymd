# GROMACS and OpenMM: how one config runs on each engine

One `config.yaml` runs on either engine (`engine: openmm` or `engine: gromacs`).
This page gives, for each setting, what each engine runs. A row marked
**1:1** runs the same physics on both engines. A row marked **approximate**
samples the same ensemble with a different algorithm, and PolyzyMD writes a
warning when a GROMACS build uses it. A row marked **differs** does not run
the same way on both engines, and PolyzyMD writes a warning.

## Dynamics

The rows hold for every equilibration stage and for production, including
the later production segments of a chain, which continue with the
thermostat, friction and barostat of the first segment.

| Setting | OpenMM | GROMACS (`.mdp`) | Match |
|---|---|---|---|
| `thermostat: LangevinMiddle` | `LangevinMiddleIntegrator(T, 1/τ, dt)` | `integrator = sd`, `tau-t = τ`, `ref-t = T`, `tcoupl = no` | 1:1. Both integrate the Langevin equation with friction 1/τ. The splitting schemes differ (BAOAB in OpenMM, leap-frog in GROMACS), at order dt². |
| `thermostat: Langevin` | `LangevinIntegrator(T, 1/τ, dt)` | as above | 1:1, as above |
| `thermostat: NoseHoover` | not supported: `validate` refuses it | `integrator = md`, `tcoupl = nose-hoover` | differs: runs only on GROMACS |
| `thermostat: Andersen` | not supported: `validate` refuses it | `integrator = md-vv`, `tcoupl = andersen` | differs: runs only on GROMACS |
| `thermostat_timescale` τ (ps) | friction = 1/τ | `tau-t` = τ | 1:1 |
| `time_step` (fs) | integrator step | `dt` (ps) | 1:1 |
| `barostat: MC` | `MonteCarloBarostat(P, T, barostat_frequency)` | `pcoupl = c-rescale`, `pcoupltype = isotropic`, `tau-p = 5 ps`, compressibility 4.5e-5 bar⁻¹ | approximate. GROMACS has no Monte Carlo barostat. Stochastic cell rescaling samples the same NPT ensemble. |
| `barostat: MCA` | not supported: `validate` refuses it | `pcoupl = c-rescale`, `pcoupltype = anisotropic` | differs: runs only on GROMACS |
| `ensemble: NPT` / `NVT` | a barostat in NPT, none in NVT | `pcoupl` as above in NPT, `pcoupl = no` in NVT | 1:1. An NPT phase without `barostat:` runs `MC`. |
| `barostat_frequency` (steps) | steps between Monte Carlo volume moves | not used: c-rescale keeps GROMACS's coupling interval (`nstpcouple`), with a warning when the value is not 25 | differs |
| `pressure` (atm) | atm | `ref-p` in bar (× 1.01325) | 1:1 |

## Seeds

| What | OpenMM | GROMACS | Match |
|---|---|---|---|
| Starting structure | Packmol and polymer draws seeded with the replicate number | the same, built by PolyzyMD before export | 1:1 |
| Initial velocities | `setVelocitiesToTemperature(T, seed)`, seed from the replicate number | `gen_seed`, seed from the replicate number | 1:1 |
| Thermostat noise | `integrator.setRandomNumberSeed(seed)`, one seed for each stage and segment | `ld_seed`, one seed for each stage | 1:1 |
| Barostat moves | `MonteCarloBarostat.setRandomNumberSeed(seed)`, one seed for each NPT stage and segment | the c-rescale noise uses `ld_seed`, which every stage's `.mdp` sets | 1:1 |

The replicate number seeds the starting structure, the initial velocities,
the thermostat noise and the barostat moves. Each equilibration stage and production segment
gets its own seed, derived from the replicate number and the stage name
(`polyzymd.simulation.seeds.dynamics_seed`). A single `ld_seed` for all
stages would repeat the same noise in every stage, because GROMACS starts
the step count again at each stage.

The two engines use different random number generators, so one seed does
not give the same trajectory on OpenMM and on GROMACS. On OpenMM's CPU,
OpenCL and CUDA platforms, a replicate run again does not give the same
trajectory either, even on one machine: CPU threads, PME and GPUs add up
forces in a different order on each run, so two runs of a replicate differ
from the first minimization on. (The slow Reference platform is
deterministic.) They agree statistically, not frame by frame. OpenMM has a
`DeterministicForces` platform property, and GROMACS has
`gmx mdrun -reprod`, for bitwise reruns; PolyzyMD sets neither.

## Nonbonded interactions and constraints

PolyzyMD writes the settings of the OpenFF force fields into the `.mdp`, so
that GROMACS computes the same energies as OpenMM through Interchange. For a
measured single-point energy of one build on both engines, see
[Engine parity](#engine-parity) below.

| What | OpenMM | GROMACS (`.mdp`) |
|---|---|---|
| Electrostatics | PME, real-space cutoff 0.9 nm, Ewald error tolerance 1e-4 (set by OpenFF Interchange; OpenMM picks the grid from it) | PME, `rcoulomb = 0.9`, `ewald-rtol = 1e-5`, `fourier-spacing = 0.12` |
| Van der Waals | cutoff 0.9 nm, switch from 0.8 nm, dispersion correction | `rvdw = 0.9`, `vdw-modifier = Potential-switch`, `rvdw-switch = 0.8`, `DispCorr = EnerPres` |
| Constraints | bonds to hydrogen; rigid water | `constraints = h-bonds`, LINCS; rigid water |

The cutoff, switch, PME and dispersion correction are the same on both
engines; the physics audit checks them against the force field file
(`nonbonded.settings`). The PME accuracy differs: the GROMACS run files ask
for a tighter reciprocal-space error than OpenMM uses.

## Equilibration

| What | OpenMM | GROMACS | Match |
|---|---|---|---|
| Minimization with `freeze_solute: true` (the default) | protein and substrate heavy atoms held fixed | `em.mdp` sets `define = -DPOSRES_EM`: position restraints of 1e5 kJ/mol/nm² on the protein and ligand heavy atoms, whatever the equilibration stages restrain | approximate. The same atoms are held; in GROMACS they can move by about 0.01 nm. |
| Position restraints | a harmonic force on the selected atoms | `define = -DPOSRES_...` with restraint `.itp` files | 1:1 for the protein, ligand and polymer groups, with one force constant per group and one protein group; `validate` refuses other position restraints on GROMACS |
| Distance restraints (`restraints:`) | applied in every phase | not applied | `validate` refuses enabled distance restraints with `engine: gromacs` |
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

`grompp_flags` reaches only the SLURM job scripts that `polyzymd submit`
writes. A local `polyzymd run` calls `grompp` without it, so the warning still
stops a local run. To run such a build on your own computer, call
`gmx grompp ... -maxwarn 1` on the exported files yourself (see
{doc}`../how_to/run_gromacs`).

(engine-parity)=
## Engine parity: are the two systems the same?

One PolyzyMD build can run on OpenMM or on GROMACS. This section gives the
evidence that both engines receive the same system: what the continuous
integration (CI) tests compare on every pull request, and a single-point
energy of the quickstart system on both engines, with the steps to repeat it.
For the checks that every build passes, see {doc}`physics_checks`.

### What CI checks on every pull request

The workflow
[`.github/workflows/engine-parity.yml`](https://github.com/joelaforet/polyzymd/actions/workflows/engine-parity.yml)
installs GROMACS from Ubuntu and runs the tests below. The audit tests build a small
system that holds every component type: Trp-cage, a substrate, polymer chains
of two monomers, DMSO, water and NaCl (`tests/_support/small_system.py`).

| Test | What it compares |
|---|---|
| `test_both_engines_read_back_the_same_system` in `tests/audit/test_components_and_engines.py` | The OpenMM `System` and the GROMACS `.top`/`.itp` files of one build, atom for atom: every bond, angle, torsion, Lennard-Jones parameter, 1-4 exception, constraint, mass and charge. A single parameter that differs fails the test. |
| `test_both_engines_run_with_the_same_settings`, same file | Constraints (`h-bonds`), nonbonded method (PME), dispersion correction, cutoff and switching distance, read from the OpenMM `System` and from the GROMACS `.mdp`. |
| `tests/audit/test_engine_parity.py` | The whole physics audit passes on the exported GROMACS files. Then each file is edited the way a person might (a bond constant in an `.itp`; the time step, length, temperature, barostat, constraints, integrator, thermostat time constant, annealing or restraint switches of an `.mdp`; a missing stage file), and the right check fails. |
| `tests/exporters/test_gromacs_staged_restraints.py` | Each equilibration stage's `.mdp` is compiled with `gmx grompp` and read back with `gmx dump`, so the test sees the position restraints GROMACS will apply, not the files PolyzyMD meant to write. |
| `tests/config/test_solutes.py` | The config's solutes and solvent, each with its own force field and charges, which both exports start from. |

These tests compare parameters and settings, not energies. CI computes no
energy on either engine.

### Single-point energy of the quickstart system

The shipped quickstart system (`examples/quickstart/config.yaml`: Trp-cage,
TIP3P water, 0.15 M NaCl, 6,138 atoms) was built once and exported both
ways: the OpenMM `System` that `polyzymd run` integrates, and the GROMACS
files that `polyzymd build --format gromacs` writes. Both engines then
evaluated the energy of the same coordinates, read from the exported `.gro`
file.

- OpenMM 8.4, `Reference` platform (double precision), one force group per
  force.
- GROMACS 2026.0 (conda-forge, mixed precision), `gmx grompp` and
  `gmx mdrun` with the `<prefix>_pointenergy.mdp` file of the export.
- PolyzyMD 1.3.0, OpenFF Toolkit 0.18.1. Measured on 2026-10-25.

Energies in kJ/mol. GROMACS reports its own nonbonded parts separately
(LJ-14, Coulomb-14, LJ (SR), dispersion correction, Coulomb (SR) and PME
reciprocal); the table adds them, because OpenMM holds them in one
`NonbondedForce`.

**After a minimization on OpenMM** (2,000 steps, written back to `.gro`):

| Term | OpenMM | GROMACS | Difference | Relative |
|---|---:|---:|---:|---:|
| Bonds | 47.42 | 47.42 | 0.00 | 4.1e-6 |
| Angles | 181.43 | 181.43 | 0.00 | 2.6e-6 |
| Torsions (proper and improper) | 948.04 | 948.04 | 0.00 | 3.8e-7 |
| Nonbonded | -103,000.55 | -102,999.73 | -0.83 | 8.0e-6 |
| Potential energy | -101,823.66 | -101,823.00 | -0.66 | 6.5e-6 |

**At the built coordinates** (before minimization):

| Term | OpenMM | GROMACS | Difference | Relative |
|---|---:|---:|---:|---:|
| Bonds | 43.35 | 43.35 | 0.00 | 1.6e-5 |
| Angles | 163.43 | 163.43 | 0.00 | 2.1e-6 |
| Torsions (proper and improper) | 980.75 | 980.75 | 0.00 | 1.8e-7 |
| Nonbonded | 112,641.94 | 112,627.24 | 14.70 | 1.3e-4 |
| Potential energy | 113,829.47 | 113,814.00 | 15.47 | 1.4e-4 |

How to read these numbers:

- The bonded terms agree to the precision GROMACS prints (six significant
  figures).
- After minimization the potential energies differ by 0.66 kJ/mol in about
  1e5, which is the size of the last digit GROMACS prints for a value of
  that size.
- The built coordinates come straight from Packmol and hold close contacts
  (a repulsive LJ (SR) of 1.1e5 kJ/mol). There the nonbonded sums differ by
  1.3e-4 of the total, from mixed against double precision on very steep
  repulsion and from the two PME implementations. Minimization removes
  these contacts, as it does before every run.
- OpenMM's `NonbondedForce` uses an Ewald error tolerance of 1e-4, which
  OpenFF Interchange sets. The `pointenergy.mdp` also uses `ewald-rtol 1e-4`;
  the run `.mdp` files use 1e-5 (see the tables above).

This comparison is not part of CI. It checks one system with no polymer,
ligand or co-solvent; CI checks the parameters of all of those, term by
term, as the table above says.

#### Repeat the comparison

You need the `build` environment and `gmx` on `PATH`. Copy
`examples/quickstart/config.yaml` and `trpcage.pdb` into an empty folder, save
the script below as `single_point.py` beside them, and run:

```bash
pixi run -e build python single_point.py <folder>
```

It builds the system (about 20 seconds), writes `<folder>/gromacs/`, and
prints both tables. The whole script took about two minutes on a laptop CPU,
most of it in the minimization.

:::{dropdown} `single_point.py`
```python
"""Single-point energy of one PolyzyMD build on OpenMM and on GROMACS.

Builds the shipped quickstart system once, exports it both ways (the OpenMM
System PolyzyMD runs, and the GROMACS files `polyzymd build --format gromacs`
writes), evaluates the energy of the same coordinates (read from the .gro on
both sides) and prints a per-term table in kJ/mol.

Usage: python single_point.py <folder holding config.yaml and trpcage.pdb>
"""

import subprocess
import sys
from pathlib import Path

import openmm
from openmm import app, unit

from polyzymd.builders.system_builder import SystemBuilder
from polyzymd.config.loader import load_config
from polyzymd.exporters.gromacs import GromacsExporter

folder = Path(sys.argv[1]).resolve()
config = load_config(folder / "config.yaml")
builder = SystemBuilder.from_config(config)
builder.build_from_config(config, working_dir=folder / "build", polymer_seed=1)
interchange = builder.create_interchange()
system = interchange.to_openmm(combine_nonbonded_forces=True)
files = GromacsExporter(interchange, config, builder.get_component_info()).export(
    folder / "gromacs"
)
print("atoms", system.getNumParticles())
for k, force in enumerate(system.getForces()):
    force.setForceGroup(k)
names = {k: type(f).__name__ for k, f in enumerate(system.getForces())}
nb = [f for f in system.getForces() if isinstance(f, openmm.NonbondedForce)][0]
print("openmm nonbonded: cutoff", nb.getCutoffDistance(), "switch", nb.getSwitchingDistance(),
      "dispcorr", nb.getUseDispersionCorrection(), "ewald tol", nb.getEwaldErrorTolerance())
gdir = Path(files["top"]).parent
mdp = next(gdir.glob("*_pointenergy.mdp"))
top = Path(files["top"]).name


def run(*a):
    return subprocess.run(a, cwd=gdir, check=True, capture_output=True, text=True)


def evaluate(gro_name, label):
    gro = app.GromacsGroFile(str(gdir / gro_name))
    box = gro.getPeriodicBoxVectors()
    system.setDefaultPeriodicBoxVectors(*box)
    context = openmm.Context(
        system, openmm.VerletIntegrator(0.001), openmm.Platform.getPlatformByName("Reference")
    )
    context.setPeriodicBoxVectors(*box)
    context.setPositions(gro.getPositions())
    omm = {}
    for k, name in names.items():
        e = context.getState(getEnergy=True, groups={k}).getPotentialEnergy()
        omm[name] = omm.get(name, 0.0) + e.value_in_unit(unit.kilojoule_per_mole)
    omm["total"] = context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(
        unit.kilojoule_per_mole
    )
    stem = Path(gro_name).stem + "_sp"
    run("gmx", "grompp", "-f", mdp.name, "-c", gro_name, "-p", top, "-o", stem + ".tpr",
        "-po", stem + "_mdout.mdp")
    run("gmx", "mdrun", "-s", stem + ".tpr", "-deffnm", stem, "-ntmpi", "1", "-ntomp", "4")
    log = (gdir / (stem + ".log")).read_text()
    block = log.split("   Energies (kJ/mol)\n")[1].split("\n\n")[0].splitlines()
    gmx = {}
    for header, values in zip(block[0::2], block[1::2]):
        keys = [header[i : i + 15].strip() for i in range(0, len(header), 15)]
        vals = [float(values[i : i + 15]) for i in range(0, len(values), 15)]
        gmx.update(zip(keys, vals))
    gmx_nb = sum(gmx.get(k, 0.0) for k in (
        "LJ-14", "Coulomb-14", "LJ (SR)", "Disper. corr.", "Coulomb (SR)", "Coul. recip."))
    rows = [
        ("Bonds", omm.get("HarmonicBondForce", 0.0), gmx.get("Bond", 0.0)),
        ("Angles", omm.get("HarmonicAngleForce", 0.0), gmx.get("Angle", 0.0)),
        ("Torsions (proper + improper)", omm.get("PeriodicTorsionForce", 0.0),
         gmx.get("Proper Dih.", 0.0) + gmx.get("Per. Imp. Dih.", 0.0)),
        ("Nonbonded (LJ, Coulomb, 1-4, PME, dispersion)", omm.get("NonbondedForce", 0.0), gmx_nb),
        ("Potential energy", omm["total"], gmx["Potential"]),
    ]
    print(f"\n[{label}]")
    print(f"{'term':48s} {'OpenMM':>14s} {'GROMACS':>14s} {'diff':>10s} {'rel':>9s}")
    for name, a, b in rows:
        rel = abs(a - b) / abs(a) if a else float("nan")
        print(f"{name:48s} {a:14.2f} {b:14.2f} {a - b:10.2f} {rel:9.1e}")
    return context


gro_name = Path(files["gro"]).name
context = evaluate(gro_name, "built coordinates")

# Minimize on OpenMM, write the result as a .gro, and compare again there.
cpu = openmm.Context(system, openmm.VerletIntegrator(0.001), openmm.Platform.getPlatformByName("CPU"))
cpu.setPeriodicBoxVectors(*app.GromacsGroFile(str(gdir / gro_name)).getPeriodicBoxVectors())
cpu.setPositions(context.getState(getPositions=True).getPositions())
openmm.LocalEnergyMinimizer.minimize(cpu, 10.0, 2000)
pos = cpu.getState(getPositions=True).getPositions(asNumpy=True).value_in_unit(unit.nanometer)
lines = (gdir / gro_name).read_text().splitlines()
out = [lines[0], lines[1]]
for line, xyz in zip(lines[2 : 2 + len(pos)], pos):
    out.append(line[:20] + "".join(f"{v:8.3f}" for v in xyz))
out.append(lines[2 + len(pos)])
(gdir / "minimized.gro").write_text("\n".join(out) + "\n")
evaluate("minimized.gro", "after OpenMM minimization")
```
:::

### What is not compared yet

- Sampled observables from runs on both engines, such as density or radius of
  gyration. The two engines integrate with different algorithms (BAOAB
  against leap-frog, Monte Carlo barostat against stochastic cell rescaling)
  and different random number generators, so their trajectories differ frame
  by frame and agree only in distribution. No such comparison is published
  here yet.
- Energies of a system with polymers, a ligand or co-solvents.
