# Create Custom Plots from Analysis Artifacts

Want to use the PolyzyMD artifacts to make your own plots? Use this guide when
you already have cached analysis artifacts and sidecars and want to make your
own matplotlib plots from those existing results. The example loads cached
hydrogen-bond aggregate artifacts and combines the `ser_his` and `asp_his`
observables on one graph without rerunning the analysis.

This workflow is intended for JupyterLab, Jupyter Notebook, VS Code notebooks,
or an IPython session.

::::{tip}
Launch your notebook server from the PolyzyMD analysis pixi environment, or select a
kernel created from that environment:

```bash
pixi run -e analysis jupyter lab
```
::::

```{important}
This is a post-processing workflow for existing artifacts and sidecars. It does
not customize a plugin's `plot()` method, load trajectories, or rerun
MDAnalysis. To adjust PolyzyMD's standard analysis plots, start with
{doc}`publication_plots` instead.
```

The code below reads small JSON artifacts from
`analysis/<condition-directory>/hydrogen_bonds/aggregated/result.json`.

## Prepare a notebook context cell

Add a short Markdown cell at the top of the notebook so exported notebooks keep
the intent of the figure clear.

````markdown
## Custom hydrogen-bond occupancy plot

This notebook loads existing PolyzyMD hydrogen-bond aggregate artifacts and
plots the `ser_his` and `asp_his` observables together. It does not rerun the
hydrogen-bond analysis.
````

## Import notebook dependencies

```python
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from polyzymd.analyses.mda import ArtifactStore
```

## Set the project paths and observables

Edit `project_dir`, `conditions`, and `condition_dirs` to match the directory
that contains your `comparison.yaml` and `analysis/` tree. The condition labels
should match the labels in `comparison.yaml`. The directory values should match
the corresponding directories under `analysis/`.

```python
project_dir = Path("/path/to/polyzymd/project").expanduser().resolve()

conditions = [
    "No Polymer",
    "100% SBMA",
    "100% EGMA",
    "1% EGPMA",
    "2% EGPMA",
    "5% EGPMA",
    "10% EGPMA",
]

condition_dirs = {
    "No Polymer": "no_polymer",
    "100% SBMA": "100_sbma",
    "100% EGMA": "100_egma",
    "1% EGPMA": "1_egpma",
    "2% EGPMA": "2_egpma",
    "5% EGPMA": "5_egpma",
    "10% EGPMA": "10_egpma",
    # Edit these to match directories under analysis/
}

observable_names = ["hbonds_ser_his", "hbonds_asp_his"]
```

If you are unsure how labels map to directories, list the available analysis
directories and update `condition_dirs` to match them.

```python
for path in sorted((project_dir / "analysis").iterdir()):
    if path.is_dir():
        print(path.name)
```

This example is useful when a `within: catalytic_triad` summary includes an
unwanted Ser-Asp component. Loading the named observables directly lets you plot
only the Ser-His and Asp-His counts.

## Validate that expected artifacts exist

```python
artifact_paths = {
    condition: project_dir
    / "analysis"
    / condition_dirs[condition]
    / "hydrogen_bonds"
    / "aggregated"
    / "result.json"
    for condition in conditions
}

missing = [path for path in artifact_paths.values() if not path.exists()]
if missing:
    raise FileNotFoundError(
        "Missing hydrogen-bond aggregate artifacts:\n"
        + "\n".join(str(path) for path in missing)
    )
```

## Load and validate hydrogen-bond artifacts

```python
condition_artifacts = {}

for condition, artifact_path in artifact_paths.items():
    artifact_dir = artifact_path.parent
    artifact = ArtifactStore(artifact_dir).read_condition_result("result.json")

    if artifact.analysis_name != "hydrogen_bonds":
        raise ValueError(
            f"Expected a hydrogen_bonds artifact for {condition!r}, "
            f"got {artifact.analysis_name!r}"
        )

    observables = artifact.payload.get("observables")
    if not isinstance(observables, list):
        raise TypeError(
            f"Expected artifact.payload['observables'] to be a list for {condition!r}"
        )

    condition_artifacts[condition] = artifact
```

## Extract a tidy DataFrame

This cell extracts one row per condition and observable. A condition aggregate
carries, for every observable:

- `mean`: mean across replicates of that replicate's mean bonds per frame.
- `sem`: standard error across replicates.
- `replicate_values`: one value per replicate, the value the tests use.

For this example, `hbonds_ser_his` is the number of Ser-His hydrogen bonds per
frame and `hbonds_asp_his` the number of Asp-His bonds per frame.

```{note}
Every contract plugin writes this same payload shape, so the cell below works
unchanged for `sasa`, `rmsf` and the others. Profiles carry `profile_mean` and
`profile_sem` instead of `mean` and `sem`; skip them by testing `kind`.
```

```python
rows = []

for condition, artifact in condition_artifacts.items():
    by_name = {
        observable["name"]: observable
        for observable in artifact.payload["observables"]
        if isinstance(observable, dict) and "name" in observable
    }

    missing = [name for name in observable_names if name not in by_name]
    if missing:
        raise KeyError(
            f"Missing observables for {condition!r}: {missing}. "
            f"Available observables: {sorted(by_name)}"
        )

    for name in observable_names:
        observable = by_name[name]
        rows.append(
            {
                "condition": condition,
                "observable": name,
                "mean": observable["mean"],
                "sem": observable["sem"],
                "replicate_values": observable["replicate_values"],
            }
        )

df = pd.DataFrame(rows)
df
```

## Plot grouped bars with SEM and replicate dots

The bars show the mean across replicates, error bars show the SEM across
replicates, and black dots show the per-replicate values.

```python
fig, ax = plt.subplots(figsize=(10, 5))

x = np.arange(len(conditions))
width = 0.8 / len(observable_names)
colors = dict(zip(observable_names, plt.get_cmap("tab10").colors))

for idx, name in enumerate(observable_names):
    offset = (idx - (len(observable_names) - 1) / 2) * width
    bar_x = x + offset

    means = []
    sems = []
    replicate_series = []

    for condition in conditions:
        row = df[(df["condition"] == condition) & (df["observable"] == name)].iloc[0]
        means.append(row["mean"])
        sems.append(row["sem"])
        replicate_series.append(row["replicate_values"])

    ax.bar(
        bar_x,
        means,
        width=width,
        yerr=sems,
        capsize=4,
        label=name,
        color=colors[name],
        edgecolor="black",
        linewidth=0.6,
        alpha=0.85,
    )

    for xpos, values in zip(bar_x, replicate_series):
        values = np.asarray(values, dtype=float)
        jitter = np.linspace(-width * 0.25, width * 0.25, num=len(values))
        ax.scatter(
            np.full_like(values, xpos) + jitter,
            values,
            color="black",
            s=24,
            zorder=3,
            alpha=0.8,
        )

ax.set_xticks(x)
ax.set_xticklabels(conditions, rotation=35, ha="right")
ax.set_ylabel("Hydrogen bonds per frame")
ax.set_title("Catalytic triad hydrogen bonds")
ax.set_ylim(bottom=0)
ax.legend(title="Observable")
ax.grid(axis="y", alpha=0.25)
fig.tight_layout()
```

## Save outside the canonical analysis tree

Save notebook-generated figures under a separate directory such as
`figures/custom/`. Avoid writing custom outputs inside the canonical `analysis/`
tree, which PolyzyMD owns.

```python
output_dir = project_dir / "figures" / "custom"
output_dir.mkdir(parents=True, exist_ok=True)

figure_path = output_dir / "hbond_ser_his_asp_his.png"
fig.savefig(figure_path, dpi=300, bbox_inches="tight")
figure_path
```

## Adapt this pattern

- Change condition labels and order by editing the `conditions` list. Use the
  display labels from `comparison.yaml`.
- Change artifact directory names by editing `condition_dirs` to match the
  directories under `analysis/`.
- Change which hydrogen-bond summaries appear by editing `observable_names`.
  An observable is named `hbonds_<summary>`, so the summary part must match a
  summary configured in `comparison.yaml`.
- If artifacts live on an HPC filesystem, copy the small JSON artifact files to
  local storage before opening the notebook to improve responsiveness.
- Use the same `ArtifactStore(...).read_condition_result("result.json")` pattern
  for other condition-level artifacts, then inspect `artifact.payload` for the
  fields you want to plot.

## Troubleshoot common notebook issues

### `ImportError: No module named polyzymd`

Start the notebook server from the analysis pixi environment:

```bash
pixi run -e analysis jupyter lab
```

If you use VS Code or an existing Jupyter server, select the kernel associated
with the PolyzyMD `analysis` environment.

### `FileNotFoundError` for `result.json`

Check that the hydrogen-bond analysis has been run and that `project_dir` points
to the directory containing `analysis/`:

```bash
pixi run -e analysis polyzymd compare run hydrogen_bonds -f comparison.yaml
```

Also list the directories under `analysis/` and compare them with
`condition_dirs`:

```python
for path in sorted((project_dir / "analysis").iterdir()):
    if path.is_dir():
        print(path.name)
```

### `KeyError` for an observable name

Print the available observable names from each artifact and compare them with
the `summaries` section of `comparison.yaml`.

```python
for condition, artifact in condition_artifacts.items():
    available = [item.get("name") for item in artifact.payload["observables"]]
    print(condition, available)
```

## See also

- {doc}`publication_plots` for settings that control PolyzyMD's standard
  analysis plots.
- {doc}`hydrogen_bonds` for configuring hydrogen-bond summaries.
- {doc}`../reference/analysis_hydrogen_bonds_reference` for hydrogen-bond
  settings and the generated output files and plots.
- {doc}`../reference/comparison_yaml` for `comparison.yaml` schema details.
- {doc}`../reference/analysis_comparison_reference` for comparison output paths
  and plotting behavior.
