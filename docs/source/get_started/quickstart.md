# Run your first simulation

In this tutorial you simulate the small protein Trp-cage in water, and then
measure its radius of gyration. You use the input files in
`examples/quickstart/` of the PolyzyMD repository. The simulation is a few
picoseconds long, so it finishes on a laptop CPU in about two minutes.

You learn these steps:

1. Check a simulation config with `polyzymd validate`.
2. Build and run a simulation with `polyzymd run`.
3. Make a {term}`study` folder with `polyzymd study init`.
4. Analyze the study with `polyzymd analyze`.

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

## Step 1: Copy the example

Copy the example folder out of the repository, and go into the copy:

```bash
cp -r examples/quickstart ~/pz_quickstart
cd ~/pz_quickstart
```

The folder holds these files:

| File | What it is |
|---|---|
| `trpcage.pdb` | Trp-cage (PDB 1L2Y, model 1), chain A, with hydrogens, cleaned with `polyzymd clean-pdb` |
| `config.yaml` | The simulation config for OpenMM, on the CPU platform |
| `config_gromacs.yaml` | The same system for GROMACS |
| `README.md` | A short description of the example |

`config.yaml` puts the protein in a rhombic dodecahedron box of TIP3P water
with 0.15 M NaCl, at 300 K and 1 atm. It has one equilibration stage of
0.002 ns (NVT) and a production of 0.004 ns (NPT) that saves 4 frames.

## Step 2: Validate the config

```bash
polyzymd validate -c config.yaml
```

`validate` reads the config and checks every key and value. It does not build
anything. The output is:

```
Validating configuration: config.yaml
Configuration is valid!


Summary:
  Name: trpcage_water
  Enzyme: trpcage
  Substrate: None (apo simulation)
  Polymers: Disabled
  Temperature: 300.0 K
  Pressure: 1.0 atm

Simulation phases:
  Equilibration: 0.002000 ns across 1 stage(s)
    - equil: 0.002 ns (NVT)
  Production: 0.004 ns (NPT)
```

## Step 3: Build and run the simulation

```bash
polyzymd run -c config.yaml -r 1
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
Output directory: /home/me/pz_quickstart/trpcage_300K_run1
```

The {term}`replicate folder` `trpcage_300K_run1/` now holds the simulation:

```text
trpcage_300K_run1/
├── build_manifest.json
├── solvated_system.pdb          # the built system, for viewers
├── system.prmtop                # the topology that the analyses read
├── system.xml                   # the OpenMM System
├── minimization/
├── equilibration_0_equil/
└── production_0/
    ├── production_0_trajectory.dcd
    └── ...
```

## Step 4: Make a study folder

A study holds conditions that you compare with each other, the analysis
protocol and the stored results. This study has one condition, `Water`:

```bash
polyzymd study init study --condition "Water=config.yaml" --equilibration 0ns
```

`--equilibration 0ns` tells the analyses to use every production frame. A
real study removes the first part of production. The output is:

```
created /home/me/pz_quickstart/study
condition Water: conditions/water/config.yaml, with 1 input files copied to structures/
git: committed 2fa172f3f01e
next: polyzymd study check /home/me/pz_quickstart/study
```

`study init` copies the config and `trpcage.pdb` into
`study/conditions/water/`. It records where the simulation is in
`study/data.local.yaml`. It also makes the folder a git repository. The
commit hash in your output is different.

## Step 5: Measure the radius of gyration

```bash
polyzymd analyze rg --study study
```

`analyze rg` measures the radius of gyration of the protein on each
production frame, and takes the mean of each replicate. The output is:

```
log: /home/me/pz_quickstart/study/logs/polyzymd-analyze-20261006-124657-pid19176.log
note: /home/me/pz_quickstart/study/study.yaml does not list rg; running it with its defaults.
# polyzymd analyze rg  metric mean_rg  unit A  eq 0ns  conditions 1  replicates 1  protocol rg/2
Water  n 1  mean 7.302  sem na  ci95 na  values 7.302  replicates 1  g 1  n_eff 4  eq_detected 0.001 ns
warning: condition Water: replicates 1 have fewer than 20 effective samples, so the start of an equilibrated region cannot be detected reliably; values and statistics are unaffected
warning: condition Water has one replicate, so it has no interval
verdict: Water mean_rg 7.302 A (no interval, n 1)
```

Read the lines in this order:

1. `verdict:` gives the result: a mean radius of gyration of about 7.3 Å.
   The replicate number seeds the run, so on the same platform and software
   versions you get the same value. Other versions or hardware give a
   slightly different value.
2. The `Water` line gives `n 1`, one replicate. With one replicate there is
   no confidence interval, so `sem` and `ci95` are `na`.
3. The two `warning:` lines come from the short test: one replicate of 4
   frames. A real study has several replicates and many frames.

For the meaning of each field, see {ref}`polyzymd analyze <cli-analyze>`.

PolyzyMD stored the result in `study/results/rg/`:

| Path | Holds |
|---|---|
| `polyzymd_results/` | The value of each frame and a record of what produced it |
| `report.json` | The full report |
| `figures/` | The figures of the analysis |

If you run the command again, PolyzyMD reads the stored values and does not
read the trajectory.

## What you did

You validated a config, built and simulated a protein in water, made a study
folder and measured one quantity with its provenance. The steps for a real
simulation are the same. Only the protein, the durations and the hardware
change.

## Next steps

- **Your own protein.** Clean the structure with
  {doc}`../tutorials/prepare_pdb_for_openff`. Then copy `config.yaml`, set
  `enzyme.pdb_path`, and make the durations longer. For each key, see
  {doc}`../reference/configuration`. `polyzymd init --name my_simulation`
  makes a simulation folder with a template config.
- **GROMACS.** If `gmx` is installed, run
  `polyzymd run -c config_gromacs.yaml -r 1`. See
  {doc}`../how_to/gromacs_export`.
- **A cluster.** Run long simulations on GPUs with `polyzymd submit`. See
  {doc}`../how_to/hpc_slurm` and {doc}`../how_to/monitor_simulations`.
- **More conditions.** Add a condition, such as a polymer or a co-solvent,
  with `polyzymd study add-condition`. See {doc}`../how_to/study_folder` and
  {doc}`../how_to/polymers`.
- **More analyses.** See {doc}`../how_to/analysis_chooser` and
  {doc}`../tutorials/first_analysis`.
