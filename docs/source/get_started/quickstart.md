# Run your first simulation

In this tutorial you simulate the small protein Trp-cage in water, and then
measure its radius of gyration. You start the way every PolyzyMD study
starts: you make a {term}`project`, add a {term}`study` to it, and add a
condition to the study. The condition is the example config in
`examples/quickstart/` of the PolyzyMD repository. The simulation is a few
picoseconds long, so it finishes on a laptop CPU in about two minutes.

You learn these steps:

1. Make a project folder with `polyzymd project init`.
2. Add a condition with `polyzymd study add-condition`.
3. Check the config with `polyzymd validate`.
4. Build and run the simulation with `polyzymd run`.
5. Analyze the project with `polyzymd analyze --project`.

## Before you start

Install PolyzyMD as {doc}`installation` describes. You need the clone of the
repository, because the example files are in it.

:::{admonition} Environment Setup
:class: tip

Run every command of this tutorial in the `build` environment. It holds
OpenMM, PACKMOL and the analysis tools. From the repository root, activate it
once:

```bash
pixi shell -e build
```

The commands then work in any folder of this shell.
:::

## Step 1: Make the project

A project holds the studies of one paper. A study holds the conditions of one
protein that you compare with each other. Make a project with one study,
`trpcage`:

```bash
polyzymd project init ~/pz_quickstart --study trpcage
cd ~/pz_quickstart
```

The output is:

```
created project /home/me/pz_quickstart with studies trpcage
next: fill description:, structures:, regions: and conditions: in each study.yaml (polyzymd study add-condition adds a condition), the analyses and metadata: in project.yaml, then polyzymd project check
```

The project folder holds `project.yaml`, which lists the studies and the
analyses, and the study folder `trpcage/`, with its `study.yaml`. The
command also makes the folder a git repository.

## Step 2: Add the example as a condition

The example folder `examples/quickstart/` holds these files:

| File | What it is |
|---|---|
| `trpcage.pdb` | Trp-cage (PDB 1L2Y, model 1), chain A, with hydrogens, cleaned with `polyzymd clean-pdb` |
| `config.yaml` | The simulation config for OpenMM, on the CPU platform |
| `config_gromacs.yaml` | The same system for GROMACS |
| `README.md` | A short description of the example |

`config.yaml` puts the protein in a rhombic dodecahedron box of TIP3P water
with 0.15 M NaCl, at 300 K and 1 atm. It has one equilibration stage of
0.002 ns (NVT) and a production of 0.004 ns (NPT) that saves 4 frames.

Add it to the study as the condition `Water`. Change `~/polyzymd` to the
folder where you cloned the repository:

```bash
polyzymd study add-condition Water --config ~/polyzymd/examples/quickstart/config.yaml --study trpcage
```

The output is:

```
condition Water: /home/me/pz_quickstart/trpcage/conditions/water/config.yaml, listed in study.yaml
warning: the runs go into runs/trpcage/water unless you set scratch_directory in config.yaml; trajectories can use a lot of disk space, so on a cluster set scratch_directory to scratch storage
next: commit, and run polyzymd study check
```

`add-condition` copies the config and `trpcage.pdb` into
`trpcage/conditions/water/`, and lists the condition in `trpcage/study.yaml`.
The copy writes its simulations into `runs/trpcage/water/` of the project.

:::{warning}
The simulations go into the project folder, in `runs/`, unless you set
`scratch_directory` in the config. Trajectories can use a lot of disk space.
On a cluster, set `scratch_directory` to your scratch storage. Git ignores
`runs/`, so the trajectories are never committed.
:::

## Step 3: Validate the config

```bash
polyzymd validate -c trpcage/conditions/water/config.yaml
```

`validate` reads the config and checks every key and value. It does not build
anything. The output is:

```
Validating configuration: trpcage/conditions/water/config.yaml
Configuration is valid!


Summary:
  Name: trpcage_water
  Engine: openmm
  Enzyme: trpcage
  Substrate: None (apo simulation)
  Polymers: Disabled
  Co-solvents: none
  Temperature: 300.0 K
  Pressure: 1.0 atm

Simulation phases:
  Equilibration: 0.002000 ns across 1 stage(s)
    - equil: 0.002 ns (NVT)
  Production: 0.004 ns (NPT)
```

## Step 4: Build and run the simulation

```bash
polyzymd run -c trpcage/conditions/water/config.yaml -r 1
```

`-r 1` selects replicate 1. The replicate number seeds the starting structure.
`polyzymd run` does these steps on this machine:

1. It builds the system: it solvates the protein with PACKMOL, adds the ions
   and assigns force-field parameters with OpenFF.
2. It minimizes the energy. The protein heavy atoms stay fixed.
3. It runs the equilibration stage.
4. It runs production.

The command prints a log of each step. It takes about one to two minutes. The
last lines are:

```
OpenMM simulation completed successfully.
Output directory: /home/me/pz_quickstart/runs/trpcage/water/trpcage_300K_run1
```

The {term}`replicate folder` `trpcage_300K_run1/` now holds the simulation:

```text
runs/trpcage/water/trpcage_300K_run1/
├── build_manifest.json
├── progress.json                # the stages and segments that have run
├── solvated_system.pdb          # the built system, for viewers
├── system.prmtop                # the topology that the analyses read
├── system.xml                   # the OpenMM System
├── minimization/
├── equilibration_0_equil/
└── production_0/
    ├── production_0_trajectory.dcd
    └── ...
```

## Step 5: Measure the radius of gyration

List the analysis in `project.yaml`. Open the file and replace the line
`analyses: {}` with these two lines:

```yaml
analyses:
  rg: {}
```

`rg: {}` runs the radius of gyration with its default settings in every study
of the project. Commit the inputs, so the report can name the commit that made
it:

```bash
git add -A
git commit -m "Add the Water condition and the rg analysis"
```

Then run every analysis of the project:

```bash
polyzymd analyze --project .
```

`analyze` measures the radius of gyration of the protein on each production
frame, and takes the mean of each replicate. The `trpcage` study uses every
production frame, because its `study.yaml` sets `equilibration: 0ns`. A real
study removes the first part of production. The output is:

```
log: /home/me/pz_quickstart/logs/polyzymd-analyze-20261006-171142-pid18811.log
== study trpcage
== rg
# polyzymd analyze rg  metric mean_rg  unit A  eq 0ns  conditions 1  replicates 1  protocol rg/2
Water  n 1  mean 7.138  sem na  ci95 na  values 7.138  replicates 1  g 1  n_eff 4  eq_detected 0.001 ns
warning: condition Water: replicates 1 have fewer than 20 effective samples, so the start of an equilibrated region cannot be detected reliably; values and statistics are unaffected
warning: condition Water has one replicate, so it has no interval
verdict: Water mean_rg 7.138 A (no interval, n 1)
```

Read the lines in this order:

1. `verdict:` gives the result: a mean radius of gyration of about 7 Å.
   Your value can differ a little, because each run of this short test
   gives a slightly different value.
2. The `Water` line gives `n 1`, one replicate. With one replicate there is
   no confidence interval, so `sem` and `ci95` are `na`.
3. The two `warning:` lines come from the short test: one replicate of 4
   frames. A real study has several replicates and many frames.

For the meaning of each field, see {ref}`polyzymd analyze <cli-analyze>`.

PolyzyMD stored the result in `trpcage/results/rg/`:

| Path | Holds |
|---|---|
| `polyzymd_results/` | The value of each frame and a record of what produced it |
| `report.json` | The full report |
| `figures/` | The figures of the analysis |

If you run the command again, PolyzyMD reads the stored values and does not
read the trajectory.

## What you did

You made a project, added a condition to its study, validated the config,
simulated a protein in water and measured one quantity with its provenance.
The steps for a real simulation are the same. Only the protein, the durations
and the hardware change.

## Next steps

- **Your own protein.** Follow {doc}`../tutorials/own_system`. It starts
  from `polyzymd study add-condition --new`, which writes a template config.
- **GROMACS.** If `gmx` is installed, add `config_gromacs.yaml` as a
  condition and run it. See {doc}`../how_to/gromacs_export`.
- **A cluster.** Run long simulations on GPUs with `polyzymd submit`. See
  {doc}`../how_to/hpc_slurm` and {doc}`../how_to/monitor_simulations`.
- **More conditions.** Copy a condition with
  `polyzymd study add-condition NAME --from Water`, and change what differs,
  such as a polymer or a co-solvent. See {doc}`../how_to/study_folder` and
  {doc}`../how_to/polymers`.
- **More studies.** Add a study for another protein with
  `polyzymd project add-study LABEL`. See {doc}`../how_to/project`.
- **More analyses.** See {doc}`../how_to/analysis_chooser` and
  {doc}`../tutorials/first_analysis`.
