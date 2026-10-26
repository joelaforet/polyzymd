# Reuse runs by sharing a condition between studies

```{index} !reuse trajectories, !reuse runs, !share a condition between studies, !use the same simulations in two studies, !simulations shared by studies
```

Analyze simulations you have already run in a second study, without running
them again: copy the condition's config into the second study. Copying a
config is cheap; the copy keeps reading the existing trajectories. Each study
then analyzes the same trajectories with its own equilibration window,
regions, control and factors.

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

So reusing a condition copies a config, never a trajectory.

## Step 1: Copy the condition's config

From the project folder, add the condition to the second study with the
config of the first:

```bash
polyzymd study add-condition "Water 300 K" --config solvents/conditions/water/config.yaml --study temperatures
```

```
condition Water 300 K: temperatures/conditions/water_300_k/config.yaml, listed in study.yaml
data Water 300 K: /home/me/Paper_1/runs/solvents/water (from the config's scratch_directory)
runs: the copy reuses the runs in /home/me/Paper_1/runs/solvents/water; nothing is run again, and new replicates of either config go there too
```

The copy's `scratch_directory` points at the existing runs. So analysis,
`polyzymd status`, `run` and `submit` find the same replicates from either
study, and a replicate run with either config goes into the same folder.

## Step 2: Add the other conditions

Add the conditions that are new to this study as usual. To start from the
reused condition, copy it with `--from`. A copy made with `--from` is a new
condition, with runs of its own:

```bash
polyzymd study add-condition "Water 330 K" --from "Water 300 K" --study temperatures
```

Then change the copy's temperatures: `thermodynamics.temperature` and the
temperatures of the equilibration stages. Give the conditions their factors in
`study.yaml` if the study compares by them.

## Step 3: Check the study

```bash
git add -A
git commit -m "Reuse the 300 K water condition in the temperature study"
polyzymd study check temperatures
```

`study check` lists the reused condition's replicates in the folder of the
`solvents` study's runs.

## What is shared, and what is not

| Shared by both studies | Each study's own |
|---|---|
| The replicates and their trajectories | The equilibration window (`equilibration:` in `study.yaml`) |
| The run folder: a replicate run or submitted with either config goes there | The control, the factors and the comparisons |
| | The structures and the regions |
| | The stored results, in each study's `results/` |
| | Its copy of the config |

So each study can discard a different equilibration window from the same
trajectories, or compare the condition with a different control.

The two copies of the config must describe the same simulation. If you edit
one, edit the other alike: a replicate run with one is not a run of the other.
`polyzymd study check` warns when two conditions read the same runs but their
configs differ:

```
warning: control Water 300 K: its runs are also those of solvents / Water, whose config differs, so replicates run with one config are not runs of the other; copy that config again, or edit the two alike
```

## After you move the runs

Run `polyzymd study locate` in each study that reads them, or edit each
`data.local.yaml`.

## Publish

`polyzymd project freeze` records the condition in both studies. Its
trajectories are listed once in `deposit/trajectories.csv` and in the upload
table of `deposit/UPLOAD.md`, named for both conditions, for example
`solvents / Water; temperatures / Water 300 K`. Upload them once. See
{doc}`study_freeze`.
