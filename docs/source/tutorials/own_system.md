# Set up your own system

In this tutorial you set up simulations of your own protein. You make a
{term}`project`, write a template config for a first condition, fill it in
and check it. Then you copy the condition to make a second one, run both and
analyze them.

The example protein is Trp-cage, from `examples/quickstart/trpcage.pdb` of
the PolyzyMD repository. Use your own cleaned PDB file instead. The two
conditions are Trp-cage in water and in 2 M urea.

## Before you start

- Do {doc}`../get_started/quickstart` first.
- Clean your structure for OpenFF, as {doc}`prepare_pdb_for_openff` shows.

:::{admonition} Environment Setup
:class: tip

Run every command of this tutorial in the `build` environment. From the
repository root, activate it once:

```bash
pixi shell -e build
```

The commands then work in any folder of this shell.
:::

## Step 1: Make the project

Make a project with one {term}`study`. A study holds the conditions of one
protein, so name it after the protein:

```bash
polyzymd project init ~/my_paper --study trpcage
cd ~/my_paper
```

```
created project /home/me/my_paper with studies trpcage
next: fill description:, structures:, regions: and conditions: in each study.yaml (polyzymd study add-condition adds a condition), the analyses and metadata: in project.yaml, then polyzymd project check
```

A study label is also its folder name. Use lower case letters, digits and
`_`. For a second protein, add a study later with
`polyzymd project add-study LABEL`.

## Step 2: Add the first condition

The first condition of a study is its control. Every other condition is
compared with it. Write a template config for it:

```bash
polyzymd study add-condition Water --new --study trpcage
```

```
condition Water: /home/me/my_paper/trpcage/conditions/water/config.yaml, listed in study.yaml
warning: the runs go into runs/trpcage/water unless you set scratch_directory in config.yaml; trajectories can use a lot of disk space, so on a cluster set scratch_directory to scratch storage
next: fill in the config, check it with polyzymd validate, commit, and run polyzymd study check
```

The command writes `trpcage/conditions/water/config.yaml` and the folder
`trpcage/conditions/water/structures/`. It lists the condition in
`trpcage/study.yaml`.

:::{warning}
The simulations go into the project folder, in `runs/`, unless you set
`scratch_directory` in the config. Trajectories can use a lot of disk space.
On a cluster, set `scratch_directory` to your scratch storage. Git ignores
`runs/`.
:::

Copy your PDB file into the `structures/` folder of the condition, and
delete the placeholder files there:

```bash
cp ~/polyzymd/examples/quickstart/trpcage.pdb trpcage/conditions/water/structures/
rm trpcage/conditions/water/structures/place_*_here.placeholder.txt
```

## Step 3: Fill in the config

Open `trpcage/conditions/water/config.yaml`. Each section has comments that
explain its keys. Change these fields:

| Field | What to write |
|---|---|
| `name` | A short name of the simulation, such as `water` |
| `description` | One sentence about the system |
| `enzyme.name` | The name of the protein. It is part of the run folder names |
| `enzyme.pdb_path` | The PDB file, relative to the config: `structures/trpcage.pdb` |
| `thermodynamics.temperature` | The temperature in K. Also change the temperatures of the equilibration stages |
| `simulation_phases.production.duration` | The production length in ns. The template has 100 ns |
| `output.scratch_directory` | On a cluster, your scratch storage. Leave it `null` to write the runs into `runs/` |
| `openmm.platform` | `CUDA` for a GPU, as `polyzymd submit` needs. `CPU` to run with `polyzymd run` on a laptop |

The optional sections, such as `substrate:`, `polymers:` and `restraints:`,
are commented out. Remove the leading `# ` of a section to use it. For every
key, see {doc}`../reference/configuration`.

## Step 4: Check the config

`validate` checks every key and value, and warns about each file that does
not exist:

```bash
polyzymd validate -c trpcage/conditions/water/config.yaml
```

If `enzyme.pdb_path` still names the template file, the output has a
warning:

```
Validating configuration: trpcage/conditions/water/config.yaml
Configuration is valid!

Referenced file warnings:
  Warning: Missing enzyme PDB: /home/me/my_paper/trpcage/conditions/water/structures/protein_X.pdb
```

Fix each warning, and run `validate` again until it prints no warning.

Then let `build --dry-run` print what it would build, without building it:

```bash
polyzymd build -c trpcage/conditions/water/config.yaml -r 1 --dry-run
```

Check the components, the phases and the folders in its report. These lines
say where the runs go:

```
Directories:
  Projects: /home/me/my_paper/runs/trpcage/water
  Scratch: /home/me/my_paper/runs/trpcage/water

Per-Replicate Output:
  Replicate 1:
    Working dir: /home/me/my_paper/runs/trpcage/water/trpcage_apo_none_100ns_300K_run1
```

The last line of the report is `Validation passed. Ready to build.`

## Step 5: Add a second condition

Copy the first condition, and change only what differs:

```bash
polyzymd study add-condition "Urea 2 M" --from Water --study trpcage
```

```
condition Urea 2 M: /home/me/my_paper/trpcage/conditions/urea_2_m/config.yaml, listed in study.yaml
warning: the runs go into runs/trpcage/urea_2_m unless you set scratch_directory in config.yaml; trajectories can use a lot of disk space, so on a cluster set scratch_directory to scratch storage
next: commit, and run polyzymd study check
```

`--from` copies the config of `Water` and the input files it names into
`trpcage/conditions/urea_2_m/`. The copy writes its runs into its own folder,
`runs/trpcage/urea_2_m/`. If you set `scratch_directory` in the first
config, set it again in the copy.

In the copy, set `name: "urea_2_m"` and replace `co_solvents: []` with the
urea:

```yaml
  co_solvents:
    - name: "urea"
      concentration: 2.0                # molar
```

Check the copy with `validate` and `build --dry-run`, as in step 4. The
summary of `validate` now gives `Co-solvents: urea`.

Commit the inputs:

```bash
git add -A
git commit -m "Add the Water and Urea 2 M conditions"
```

## Step 6: Run the simulations

Run three replicates of each condition. The replicate number is the random
seed.

On a cluster with GPUs, submit them to SLURM. Choose the preset of your
cluster:

```bash
polyzymd submit -c trpcage/conditions/water/config.yaml -r 1-3 --preset aa100
polyzymd submit -c trpcage/conditions/urea_2_m/config.yaml -r 1-3 --preset aa100
```

See {doc}`../how_to/hpc_slurm` and {doc}`../how_to/monitor_simulations`.

On a workstation, run one replicate after the other with `polyzymd run`:

```bash
polyzymd run -c trpcage/conditions/water/config.yaml -r 1
```

When the runs are done, check that PolyzyMD finds them:

```bash
polyzymd project check .
```

Each condition line names its runs, such as
`control Water: runs [1, 2, 3] under /home/me/my_paper/runs/trpcage/water (from config)`.

## Step 7: Analyze

In `trpcage/study.yaml`, set `equilibration:` to the part of production to
discard, such as `20ns`. To choose it, see {doc}`../how_to/equilibration`.

In `project.yaml`, list the analyses, such as the radius of gyration and the
RMSF of the C-alpha atoms:

```yaml
analyses:
  rg: {}
  rmsf:
    selection: protein and name CA
```

Commit, and run every analysis of the project:

```bash
git add -A
git commit -m "Set the equilibration window and the analyses"
polyzymd analyze --project .
```

The report of each analysis compares `Urea 2 M` with the control, `Water`.
To read the reports, see {doc}`first_analysis`.

## What you did

You made a project for your own protein, filled in a template config, checked
it, copied it to make a second condition, and ran and analyzed both. To add a
polymer, see {doc}`../how_to/polymers`. To publish the project, see
{doc}`../how_to/project`.
