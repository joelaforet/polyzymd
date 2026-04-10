# Run Your First PolyzyMD Simulation

This tutorial walks through one complete first run: create a project scaffold,
add an enzyme structure, write a minimal configuration, validate it, and make
sure PolyzyMD can build the system.

## What you need

- PolyzyMD installed with `pixi` as described in {doc}`installation`
- a simulation-ready enzyme PDB file
- the repository available locally so `pixi run -e ...` can find `pixi.toml`

By the end, you will have a working project directory and a validated
`config.yaml` that is ready for local building or HPC submission.

## Step 1: Create a new project scaffold

From the PolyzyMD repository root, run:

```bash
pixi run -e build polyzymd init --name my_first_simulation
```

This creates a directory like:

```text
my_first_simulation/
|- config.yaml
|- structures/
|- job_scripts/
`- slurm_logs/
```

Move into the new project:

```bash
cd my_first_simulation
```

## Step 2: Add your enzyme structure

Copy your prepared enzyme structure into `structures/`:

```bash
cp /path/to/enzyme.pdb structures/enzyme.pdb
rm structures/*.placeholder.txt
```

Your PDB should already be ready for simulation:

- hydrogens added
- missing atoms or residues fixed
- intended protonation state chosen
- alternate conformers removed

If you are also using a substrate, place its SDF file in `structures/` now, but
this tutorial keeps the first run to an enzyme-only example.

## Step 3: Replace the template with a minimal config

Open `config.yaml` and replace the template contents with:

```yaml
name: "my_first_simulation"
description: "First PolyzyMD tutorial run"

enzyme:
  name: "MyEnzyme"
  pdb_path: "structures/enzyme.pdb"

solvent:
  primary:
    type: "water"
    model: "tip3p"
  co_solvents: []
  ions:
    neutralize: true
    nacl_concentration: 0.15
  box:
    padding: 1.2
    shape: "rhombic_dodecahedron"
    target_density: 1.0
    tolerance: 2.0

thermodynamics:
  temperature: 300.0
  pressure: 1.0

simulation_phases:
  equilibration_stages:
    - name: "heating"
      duration: 0.2
      samples: 20
      ensemble: "NVT"
      temperature_start: 60.0
      temperature_end: 300.0
      temperature_increment: 1.0
      temperature_interval: 1200.0
      position_restraints:
        - group: "protein_heavy"
          force_constant: 4184.0
    - name: "free_equilibration"
      duration: 0.8
      samples: 80
      ensemble: "NPT"
      temperature: 300.0
  production:
    ensemble: "NPT"
    duration: 10.0
    samples: 250
    time_step: 2.0
    thermostat: "LangevinMiddle"
    thermostat_timescale: 1.0
    barostat: "MC"
    barostat_frequency: 25

output:
  projects_directory: "."
  scratch_directory: null
  job_scripts_subdir: "job_scripts"
  slurm_logs_subdir: "slurm_logs"
  naming_template: "{enzyme}_{temperature}K_run{replicate}"
  save_checkpoint: true
  save_state_data: true
  trajectory_format: "dcd"
```

This is intentionally small, but it still uses the staged equilibration model
that PolyzyMD now requires for all simulations.

## Step 4: Validate the config

Run:

```bash
pixi run -e build polyzymd validate -c config.yaml
```

Success looks like a clean validation report with your enzyme name,
temperature, and equilibration/production phases listed.

If validation fails, fix the reported paths or field values before continuing.

## Step 5: Check that the system can build

Start with a dry run to see the full validation report:

```bash
pixi run -e build polyzymd build -c config.yaml --dry-run
```

If that succeeds, run the actual build. The commands differ depending on
which simulation engine you plan to use:

`````{tab-set}
````{tab-item} OpenMM (default)
Build the system for OpenMM simulation:

```bash
pixi run -e build polyzymd build -c config.yaml
```

This prepares the solvated system (PDB + OpenMM XML) for the configured
equilibration and production phases.
````

````{tab-item} GROMACS
Build the system and export GROMACS input files:

```bash
pixi run -e build polyzymd build -c config.yaml --format gromacs
```

This builds the solvated system and exports `.gro`, `.top`, `.itp`, `.mdp`, and
a run script to `replicate_1/gromacs/`.

```{tip}
For the full GROMACS workflow, see {doc}`../how_to/gromacs_export`.
`polyzymd submit` is OpenMM-only for now.
Use `polyzymd build --format gromacs` when you only want input files,
or `polyzymd run --engine gromacs` when you want PolyzyMD to build and run locally.
Integrated GROMACS SLURM submission is planned for v1.4.0.
```

For multiple replicates (each with an independently built system):

```bash
pixi run -e build polyzymd build -c config.yaml --format gromacs --replicates 1-3
```
````
`````

## Step 6: Decide how you want to run production

`````{tab-set}
````{tab-item} OpenMM — local
For a local smoke test on a GPU-enabled workstation, you can run one segment:

```bash
pixi run -e cuda-12-4 polyzymd run-segment -c config.yaml -r 1
```
````

````{tab-item} OpenMM — HPC
For an HPC workflow, generate job scripts first:

```bash
pixi run -e cuda-12-4 polyzymd submit -c config.yaml --preset aa100 --replicates 1 --dry-run
```

If the dry run looks right, submit for real with the CUDA environment that
matches your cluster.
````

````{tab-item} GROMACS — local
Run the full GROMACS workflow locally (build + EM + equilibration + production):

```bash
pixi run -e build polyzymd run -c config.yaml --engine gromacs --replicates 1
```

If you only need the input files without running GROMACS, use
`build --format gromacs` instead (see Step 5 above).
````

````{tab-item} GROMACS — HPC
Export GROMACS files and submit manually with your cluster's SLURM scripts:

```bash
pixi run -e build polyzymd build -c config.yaml --format gromacs --replicates 1-3
```

```{note}
`polyzymd submit` does not yet support GROMACS as a simulation engine.
Export files with `build --format gromacs` and submit manually via SLURM.
Integrated GROMACS submission support is planned for v1.4.0.
```

See {doc}`../how_to/gromacs_export` for the full GROMACS HPC workflow.
````
`````

## What you should have now

At this point you should have:

- a project directory created by `polyzymd init`
- a minimal, validated `config.yaml`
- a successful build dry run or real build
- a clear next step for local execution or SLURM submission

## Where to go next

- Add a substrate or polymers: {doc}`../reference/configuration`
- Export to GROMACS: {doc}`../how_to/gromacs_export`
- Add distance restraints: {doc}`../how_to/restraints`
- Tune staged equilibration: {doc}`../how_to/equilibration`
- Run on a cluster: {doc}`../how_to/hpc_slurm`

<!-- IMAGE OPPORTUNITY: Add a screenshot of the generated project scaffold just
after `polyzymd init`, with `config.yaml`, `structures/`, `job_scripts/`, and
`slurm_logs/` annotated. -->
