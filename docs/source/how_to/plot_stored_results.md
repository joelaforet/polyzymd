# Plot stored results with NumPy and matplotlib

```{index} !plot results with matplotlib, !results as numpy array, !own plots, !matplotlib, !numpy, !read stored results
```

Make your own figures from the values that `polyzymd analyze` stored, with
NumPy and matplotlib or any other plotting library. PolyzyMD is the data
source: you read its results as arrays and draw them as you like. No
trajectory is loaded, so this works on a copy of a study, or a frozen deposit,
that has no trajectories.

To change PolyzyMD's own figures instead, see {doc}`publication_plots`. To
measure a new quantity and plot it, see {doc}`custom_artifact_plotting`.

The examples use the project of {doc}`../tutorials/analysis_complete_workflow`:
the study `trpcage` with the conditions `Water` and `SBMA`, three replicates
each, after `polyzymd analyze --project .`. Run them from the project folder,
in the `build` or `analysis` environment, for example in a notebook started
with `pixi run -e analysis jupyter lab`.

```{figure} ../_static/figures/plot_stored_results.webp
:figclass: pz-render
:class: dark-light
:alt: Three plots. Left, the radius of gyration of each replicate against time, blue for Water and orange for SBMA. Middle, the RMSF of each residue, a line per condition with a shaded 95 percent interval. Right, the mean radius of gyration of each condition as a bar with its interval and the three replicate values as dots.

The three figures this page draws, from the tutorial's short runs: a time
series per replicate, a per-residue profile with its 95% interval, and one
value per replicate with the mean and its interval.
```

## Step 1: Read the stored results

```python
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats

import polyzymd as pz

study = pz.Study("trpcage")
rg = study.results("rg")
rg.table.head()
```

```
  name condition  replicate  part label  frame  time_ns     value unit
0   rg      SBMA          1  None  None      0    0.001  7.291977    A
1   rg      SBMA          1  None  None      1    0.002  7.221196    A
2   rg      SBMA          1  None  None      2    0.003  7.202337    A
3   rg      SBMA          1  None  None      3    0.004  7.255366    A
4   rg      SBMA          2  None  None      0    0.001  7.322464    A
```

`study.results(run)` reads `trpcage/results/<run>/`, where `polyzymd analyze`
stored the run. Its `table` is a pandas DataFrame with one row per stored
value:

| Column | Holds |
|---|---|
| `name` | The measured quantity. `rg.names` lists them |
| `condition`, `replicate` | Where the value comes from |
| `frame`, `time_ns` | For a per-frame series: the frame and its time |
| `label` | For a per-residue (or other labelled) result: the residue |
| `part` | For a result with several parts, such as RMSF's `rmsf`, `offset` and `rmsd_per_residue` |
| `value`, `unit` | The value and its unit |

`study.labels` gives the conditions in the study's order, control first.

## Step 2: Plot a time series per replicate

Select rows with pandas, and take NumPy arrays with `.to_numpy()`:

```python
colors = dict(zip(study.labels, plt.rcParams["axes.prop_cycle"].by_key()["color"]))
fig, ax = plt.subplots(figsize=(6, 3.5))
for (condition, replicate), rows in rg.table.groupby(["condition", "replicate"]):
    time = rows["time_ns"].to_numpy()
    value = rows["value"].to_numpy()
    ax.plot(time, value, marker="o", color=colors[condition],
            label=condition if replicate == 1 else None)
ax.set_xlabel("Time (ns)")
ax.set_ylabel("Radius of gyration (Å)")
ax.legend(frameon=False)
```

The table holds the frames that the analysis measured.

## Step 3: Plot a per-residue profile

A per-residue result has one row per residue. Pivot it into an array with one
row per replicate and one column per residue, then compute what you want to
draw, here the mean and its Student t 95% interval across the replicates:

```python
rmsf = study.results("rmsf").table
rmsf = rmsf[rmsf["part"] == "rmsf"]
fig, ax = plt.subplots(figsize=(6, 3.5))
for condition in study.labels:
    rows = rmsf[rmsf["condition"] == condition]
    grid = rows.pivot(index="replicate", columns="label", values="value")
    residues = grid.columns.to_numpy()
    values = grid.to_numpy()  # one row per replicate, one column per residue
    mean = values.mean(axis=0)
    half = stats.t.ppf(0.975, len(values) - 1) * stats.sem(values, axis=0)
    ax.plot(residues, mean, color=colors[condition], label=condition)
    ax.fill_between(residues, mean - half, mean + half, color=colors[condition], alpha=0.25)
ax.set_xlabel("Residue")
ax.set_ylabel("RMSF (Å)")
ax.legend(frameon=False)
```

The replicate is the sampling unit: compute the interval across replicates,
not across frames.

## Step 4: Plot one value per replicate, with the reported statistics

`results.report` is the report that `polyzymd analyze` printed. Each of its
conditions has the mean, its 95% interval and the replicate values, so a
figure can show exactly the numbers of the report:

```python
report = rg.report
fig, ax = plt.subplots(figsize=(3.5, 3.5))
for i, row in enumerate(report.conditions):
    errors = None if row.ci95 is None else [[row.mean - row.ci95[0]], [row.ci95[1] - row.mean]]
    ax.bar(i, row.mean, yerr=errors, capsize=4, color=colors[row.label], alpha=0.7)
    ax.scatter(np.full(len(row.replicate_values), i), row.replicate_values,
               color="black", s=12, zorder=3)
ax.set_xticks(range(len(report.conditions)), [row.label for row in report.conditions])
ax.set_ylabel(f"Mean radius of gyration ({report.unit})")
```

`row.ci95` is `None` when no interval can be estimated, as with one replicate.
For the per-replicate values of every quantity as a table, with a column for
each factor of the study, use `study.replicate_table("rg")`.

## Step 5: Save the figures

Save your figures outside `results/`, so the stored values and their records
stay as PolyzyMD wrote them. A study's `figures/` folder is deposited by
freeze:

```python
fig.savefig("trpcage/figures/rg_means.pdf", bbox_inches="tight")
```

## Read the files without PolyzyMD

Each replicate's values are plain NumPy files, so any tool can read them:

```
trpcage/results/rg/polyzymd_results/<name>/<condition>/replicate_<n>/
    record.json    what was measured, how, and on which files
    series.npz     a per-frame series: arrays values, frames and times (ns)
    values.npz     a per-replicate value or profile: array values
    labels.json    the labels, when the function returned its own
```

```python
with np.load("trpcage/results/rg/polyzymd_results/rg/<condition>/replicate_1/series.npz") as data:
    times, values = data["times"], data["values"]
```

`study.results` does this for every replicate and checks the values against
the report, so prefer it when PolyzyMD is installed.

## See also

- {doc}`study_yaml`: the Study object at a glance, and results of several
  studies with `pz.Project`.
- {doc}`custom_artifact_plotting`: measure a new quantity, then plot it.
- {doc}`publication_plots`: change PolyzyMD's own figures.
