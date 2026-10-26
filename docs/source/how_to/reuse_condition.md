# Reuse runs by sharing a condition between studies

```{index} !reuse trajectories, !reuse runs, !share a condition between studies, !use the same simulations in two studies, !simulations shared by studies
```

Analyze simulations you have already run in a second study, without running
them again. One condition, with its config and its replicates, then belongs
to two studies. Each study analyzes the same trajectories with its own
equilibration window, regions, control and factors.

Reuse a condition when the same simulations answer more than one question.
For example, a project has a study `solvents` that compares a protein with
SBMA chains in water and in 10 % DMSO, at 300 K. A second study,
`temperatures`, compares the same system in water at 300 K, 330 K and 360 K.
The 300 K water condition belongs to both studies, and it is simulated once.

## What the project keeps, and where the trajectories are

A project folder holds what describes the simulations and the analyses: each
condition's `config.yaml` and its input files, `study.yaml`, `project.yaml`,
the analysis code and the stored results. Git tracks these, and freeze
publishes them.

The trajectories are not part of the project. Each config says where its runs
are written: `runs/<study>/<condition>/` in the project folder by default, which
git ignores, or its `scratch_directory`. When you move the runs, for example to
long-term storage, or download them, each study's `data.local.yaml` says where
they are on that machine (see {doc}`study_folder`).

So reusing a condition copies no trajectory. The second study lists the first
study's config, and it finds the runs where that config, or `data.local.yaml`,
says they are.

## Before you start

- Both studies are in one project. To bring existing studies into a project,
  see {doc}`move_studies_into_project`.
- The condition to reuse is listed in the first study, here `Water` in
  `solvents`.

## Step 1: Reuse the condition

From the project folder, give the new study's label for the condition, and
the study and condition to reuse:

```bash
polyzymd study add-condition "Water 300 K" --reuse "solvents:Water" --study temperatures
```

```
condition Water 300 K: reuses solvents/conditions/water/config.yaml, listed in study.yaml; nothing was copied
runs: shared with the other study. Both analyze the same replicates, each with its own equilibration window and control; a replicate run or submitted from either goes into the same folder
next: add factors: to the entry in study.yaml if the study compares by them, commit, and run polyzymd study check
```

`temperatures/study.yaml` now lists the config of the other study:

```yaml
conditions:
  "Water 300 K": "../solvents/conditions/water/config.yaml"
```

When the label is the same in both studies, `--reuse solvents` is enough.

## Step 2: Add factors, and the other conditions

If the study compares its conditions by a factor, write the entry with its
factors:

```yaml
conditions:
  "Water 300 K": {config: ../solvents/conditions/water/config.yaml, factors: {temperature_K: 300}}
```

Add the conditions that are new to this study as usual. To start from the
reused condition, copy it with `--from`. The copy is a new condition of
`temperatures`, with its own runs:

```bash
polyzymd study add-condition "Water 330 K" --from "Water 300 K" --study temperatures
```

Then change the copy's temperatures: `thermodynamics.temperature` and the
temperatures of the equilibration stages.

## Step 3: Check the study

```bash
git add -A
git commit -m "Reuse the 300 K water condition in the temperature study"
polyzymd study check temperatures
```

`study check` lists the reused condition's replicates, in the folder of the
`solvents` study's runs.

## What is shared, and what is not

| Shared by both studies | Each study's own |
|---|---|
| The condition's `config.yaml` and input files | The equilibration window (`equilibration:` in `study.yaml`) |
| The replicates and their trajectories | The control, the factors and the comparisons |
| The run folder: a replicate run or submitted with this config goes there | The structures and the regions |
| | The stored results, in each study's `results/` |

So each study can discard a different equilibration window from the same
trajectories, or compare the condition with a different control.

The config is one file. Change it only for a reason that holds in both
studies: a replicate run after a change of its physics is not comparable with
the earlier ones, in either study.

## Run more replicates

Run or submit more replicates with the config, from either study. They go into
the shared run folder, and both studies analyze them:

```bash
polyzymd run -c solvents/conditions/water/config.yaml -r 4
```

## After you move the runs

`--reuse` copies the other study's `data.local.yaml` entry for the condition,
if it has one. After you move the runs again, run `polyzymd study locate` in
each study that lists the condition, or edit each `data.local.yaml`.

## Publish

`polyzymd project freeze` records the reused condition in both studies. Its
trajectories are listed once in `deposit/trajectories.csv` and in the upload
table of `deposit/UPLOAD.md`, named for both conditions, for example
`solvents / Water; temperatures / Water 300 K`. Upload them once. See
{doc}`study_freeze`.

## Limits

- A condition is reused within one project. A study in another project, or in
  no project, copies the config with `--config` instead: that makes a new
  condition, whose new replicates go into its own runs folder.
- A label is listed once per study, so a study cannot reuse its own
  condition. To copy one of its conditions and change it, use `--from`.
- `add-condition --config` with the config of another study's condition makes
  a copy, not a reuse, and says so: it prints the `--reuse` command to use
  instead.
