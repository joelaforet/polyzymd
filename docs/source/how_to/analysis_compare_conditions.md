# Compare simulation conditions

Find out whether the conditions of a study differ, for example an enzyme with
and without a polymer. You need finished simulations of each condition.

The steps are:

1. Give the `config.yaml` of each condition, the control first.
2. Run `polyzymd analyze NAME -c ... -c ...` for each analysis.
3. Read the value of each condition and its comparison with the control.
4. Find the stored results and the figures.

For a guided lesson, see {doc}`../tutorials/analysis_complete_workflow`. For
the settings of each analysis, see its how-to under {doc}`index`.

:::{admonition} Environment Setup
:class: tip

All analysis commands below assume you have activated the PolyzyMD analysis
pixi environment:

```bash
pixi shell -e analysis
```

Alternatively, prefix each command with `pixi run -e analysis`.
:::

:::{admonition} Resource requirements
:class: important

`polyzymd analyze` loads trajectories. It can use much memory, CPU time and
disk I/O. On a shared cluster, run it in a batch job or an interactive job,
not on a login node. For a batch script, see {doc}`hpc_execution`.
:::

## Study or configs

This guide gives the configs with `-c`. In a {term}`study`, give the study
instead:

```bash
polyzymd analyze hydrogen_bonds --study my_study
```

The study names the conditions, the control, the replicates and the
equilibration window. The results go to `<study>/results/<name>/`, and
`polyzymd project freeze` can publish them. Use `--study` or `--project` for
the results that you keep.

Use `-c` for a quick look without a study: to check a simulation while it
runs, to try a setting, or to compare configs that belong to no study. The
results then go to the current folder, or to `--output-dir`. For a lesson
with both, see {doc}`../tutorials/first_analysis`.

## Before you start

Make sure that each condition has:

- a simulation `config.yaml`;
- finished trajectories for the replicates that you want to compare.

`polyzymd analyze` stores the result of each replicate under
`polyzymd_results/`. A later command reuses it if the function, the settings,
the config, the input files, the equilibration window and the frames are
unchanged. If the trajectory of a replicate grew, PolyzyMD measures the
replicate again. To measure every replicate again, add `--recompute`.

## Analyze simulations that are still running

An analysis leaves out each OpenMM production segment that `progress.json`
records as running or failed. The trajectory file of such a segment ends at
the last flush. A warning names the segments that PolyzyMD left out.

- **The window ends before the left-out segments.** The results describe the
  completed part of the simulation. When the simulation finishes, run the
  command again with `--recompute` to use the full window.
- **The left-out segment makes a gap.** A left-out segment lies between two
  kept segments, so the lineage check refuses the replicate. Wait for the
  simulation to finish. Or load the replicate with `require_complete=False`:

  ```python
  from polyzymd.analyses.shared.loader import TrajectoryLoader
  u = TrajectoryLoader(config).load_universe(replicate=1, require_complete=False)
  ```

  PolyzyMD then reads the incomplete segments as they are. It lists them in
  `incomplete_segments` on the layout.

GROMACS records no status for each segment. PolyzyMD reads a production XTC
that GROMACS is still writing as it is. Make sure that the job has finished
before you analyze a GROMACS simulation.

## Step 1: Run one comparison

Give one `-c` for each condition, the control first. Name the conditions with
`--label`, in the same order:

```bash
polyzymd analyze hydrogen_bonds \
  -c ../SBMA_100_enzyme_DMSO/config.yaml \
  -c ../EGMA_100_enzyme_DMSO/config.yaml \
  -c ../SBMA_50_enzyme_DMSO/config.yaml \
  --label "100% SBMA" --label "100% EGMA" --label "50% SBMA" \
  --replicates 1-3 --eq 10ns
```

The command does these steps:

1. It reads each condition from its `config.yaml`. It uses the replicates of
   `--replicates`, or without it, every replicate it finds.
2. It discards the first 10 ns of each replicate.
3. On each production frame, it counts the hydrogen bonds between the protein
   (chain A) and the polymer (chain C). These are the default groups.
4. It prints the mean of each condition with its 95 % confidence interval.
5. It compares each condition with the control by Welch's t test. It corrects
   the p values with the {term}`Benjamini-Hochberg` method.

A condition without polymer can be the control. Its replicates report 0
hydrogen bonds and stay in the statistics.

Without `--label`, each condition takes the name of the folder that holds its
config. For the meaning of each output line, see {doc}`analysis_agent_protocol`.
For the other settings, such as named groups and summaries, see
{doc}`hydrogen_bonds`.

## Step 2: Report another result

An analysis with several results reports one at a time. The JSON report lists
the others in `all_runs`. To report a different result, use `--run`:

```bash
polyzymd analyze hydrogen_bonds -c A/config.yaml -c B/config.yaml --eq 10ns \
  --run protein_polymer_residues
```

PolyzyMD reads back the replicate results that the first command stored. Only
a result that is not measured yet loads the trajectories.

## Step 3: Save the full report

`--format json` prints the whole `ProtocolReport`. `-o` also writes it to a
file:

```bash
polyzymd analyze hydrogen_bonds -c A/config.yaml -c B/config.yaml --eq 10ns \
  --format json -o reports/hydrogen_bonds.json
```

## Step 4: Run the other analyses

Every analysis takes the same `-c`, `--label`, `--replicates` and `--eq`
options. Run each analysis on the same conditions:

```bash
polyzymd analyze rmsf -c A/config.yaml -c B/config.yaml --eq 10ns
polyzymd analyze contacts -c A/config.yaml -c B/config.yaml --eq 10ns \
  --set polymer_selection='chainid C' --set protein_selection='chainid A'
polyzymd analyze distances -c A/config.yaml -c B/config.yaml --eq 10ns \
  --set pairs=pairs.yaml
```

For the distances of a catalytic triad, write the pairs to `pairs.yaml`. For
its hydrogen bonds, see {doc}`analysis_triad_quickstart`.

## Step 5: Check the outputs

After you run the commands in `polymer_stability_study/`, the folder holds
files such as these:

```text
polymer_stability_study/
├── polyzymd_results/
│   └── hydrogen_bonds_protein_polymer/
│       ├── 100_SBMA/
│       │   ├── replicate_1/
│       │   │   ├── record.json
│       │   │   └── values.npz
│       │   └── ...
│       └── 100_EGMA/
│           └── ...
└── figures/
    └── hydrogen_bonds/
        └── hbonds_protein_polymer_mean_hbonds_comparison.png
```

Each folder name is the result or condition label. PolyzyMD replaces each
series of characters other than letters, digits, `.`, `+` and `-` with `_`.

- `record.json` records what was measured and from which inputs.
- `values.npz` holds the values of the replicate.

`--output-dir` puts `polyzymd_results/` and `figures/` in another folder.
`--no-plots` draws no figures.

## From Python

In Python, `analyze` returns the same report:

```python
from polyzymd.analyses import analyze

report = analyze(
    "hydrogen_bonds",
    ["A/config.yaml", "B/config.yaml"],
    labels=["100% SBMA", "100% EGMA"],
    equilibration="10ns",
)
print(report.to_agent_text())
```

To measure your own function on the same conditions, see {doc}`study_api`.

## Troubleshooting

### `config` path not found

PolyzyMD resolves a relative `-c` path from the current folder of the shell.

### `replicates ... have no run directory`

The replicates of `--replicates` have no replicate folder under the scratch
directory of the config. Do one of these:

- Leave out `--replicates` to use every replicate that PolyzyMD finds.
- Run the simulations first.
- If the replicate folders are in another place, name the place in
  `data.local.yaml`, with `polyzymd study locate DIR`, or with `--data`.

### `the control ... has no replicate where every selection matches atoms`

A selection of the analysis matched no atoms in any replicate of the control,
so PolyzyMD gives a summary of the other conditions and does not compare them.
An empty polymer selection does not cause this warning: a replicate without
polymer reports 0. Check the other selections of the control, such as
`protein_selection` or the groups of a hydrogen-bond summary. Or give a
condition in which every selection matches atoms first.

## See also

- {doc}`../tutorials/analysis_complete_workflow`
- {doc}`analysis_agent_protocol`
- {doc}`../explanation/analysis_statistics_best_practices`
