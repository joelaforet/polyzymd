# Create Custom Plots from Analysis Artifacts

Want to use the PolyzyMD artifacts to make your own plots? Use this guide when
you already have cached analysis artifacts and sidecars and want to make your
own matplotlib plots from those existing results. The example loads cached
hydrogen-bond aggregate artifacts and combines the `ser_his` and `asp_his`
summaries on one graph without rerunning the analysis.

This workflow is intended for JupyterLab, Jupyter Notebook, VS Code notebooks,
or an IPython session.

::::{tip}
Launch your notebook server from the PolyzyMD pixi environment, or select a
kernel created from that environment:

```bash
pixi run -e build jupyter lab
```
::::

```{important}
This is a post-processing workflow for existing artifacts and sidecars. It does
not customize a plugin's `plot()` method, load trajectories, or rerun
MDAnalysis. To customize PolyzyMD's regenerated built-in plots, start with
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
plots the `ser_his` and `asp_his` summaries together. It does not rerun the
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

## Set the project paths and summaries

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

summary_names = ["ser_his", "asp_his"]
```

If you are unsure how labels map to directories, list the available analysis
directories and update `condition_dirs` to match them.

```python
for path in sorted((project_dir / "analysis").iterdir()):
    if path.is_dir():
        print(path.name)
```

This example is useful when the built-in hydrogen-bond plot for
`within: catalytic_triad` includes an unwanted Ser-Asp component. Loading the
named summaries directly lets you plot only Ser-His and Asp-His occupancy.

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

    summaries = artifact.payload.get("summaries")
    if not isinstance(summaries, list):
        raise TypeError(
            f"Expected artifact.payload['summaries'] to be a list for {condition!r}"
        )

    condition_artifacts[condition] = artifact
```

## Extract a tidy DataFrame

This cell extracts one row per condition and summary. The occupancy fields mean:

- `mean_fraction_with_any`: mean across replicates of the fraction of analyzed
  frames with at least one H-bond matching that summary.
- `sem_fraction_with_any`: SEM across replicates.
- `per_replicate_fraction_with_any`: one occupancy value per replicate.

For this example, `ser_his` is the fraction of frames with at least one Ser-His
H-bond, and `asp_his` is the fraction of frames with at least one Asp-His
H-bond.

```{note}
Artifact payload shapes are plugin- and version-specific. This example uses the
documented hydrogen-bond summary fields from the current hydrogen-bond artifact
output.
```

```python
rows = []

for condition, artifact in condition_artifacts.items():
    summaries_by_name = {
        summary["name"]: summary
        for summary in artifact.payload["summaries"]
        if isinstance(summary, dict) and "name" in summary
    }

    missing_summaries = [
        summary_name
        for summary_name in summary_names
        if summary_name not in summaries_by_name
    ]
    if missing_summaries:
        available = sorted(summaries_by_name)
        raise KeyError(
            f"Missing summaries for {condition!r}: {missing_summaries}. "
            f"Available summaries: {available}"
        )

    for summary_name in summary_names:
        summary = summaries_by_name[summary_name]
        replicate_values = summary["per_replicate_fraction_with_any"]

        rows.append(
            {
                "condition": condition,
                "summary": summary_name,
                "mean_fraction_with_any": summary["mean_fraction_with_any"],
                "sem_fraction_with_any": summary["sem_fraction_with_any"],
                "per_replicate_fraction_with_any": replicate_values,
            }
        )

df = pd.DataFrame(rows)
df
```

## Plot grouped bars with SEM and replicate dots

The bars show mean occupancy, error bars show SEM across replicates, and black
dots show per-replicate occupancy values.

```python
fig, ax = plt.subplots(figsize=(10, 5))

x = np.arange(len(conditions))
width = 0.8 / len(summary_names)
colors = dict(zip(summary_names, plt.get_cmap("tab10").colors))

for idx, summary_name in enumerate(summary_names):
    offset = (idx - (len(summary_names) - 1) / 2) * width
    bar_x = x + offset

    means = []
    sems = []
    replicate_series = []

    for condition in conditions:
        row = df[(df["condition"] == condition) & (df["summary"] == summary_name)].iloc[0]
        means.append(row["mean_fraction_with_any"])
        sems.append(row["sem_fraction_with_any"])
        replicate_series.append(row["per_replicate_fraction_with_any"])

    ax.bar(
        bar_x,
        means,
        width=width,
        yerr=sems,
        capsize=4,
        label=summary_name,
        color=colors[summary_name],
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
ax.set_ylabel("Fraction of analyzed frames with ≥1 H-bond")
ax.set_title("Catalytic triad hydrogen-bond occupancy")
ax.set_ylim(bottom=0)
ax.legend(title="Summary")
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

figure_path = output_dir / "hbond_ser_his_asp_his_occupancy.png"
fig.savefig(figure_path, dpi=300, bbox_inches="tight")
figure_path
```

## Adapt this pattern

- Change condition labels and order by editing the `conditions` list. Use the
  display labels from `comparison.yaml`.
- Change artifact directory names by editing `condition_dirs` to match the
  directories under `analysis/`.
- Change which hydrogen-bond summaries appear by editing `summary_names`.
- Summary names must match the names configured in `comparison.yaml`.
- If artifacts live on an HPC filesystem, copy the small JSON artifact files to
  local storage before opening the notebook to improve responsiveness.
- Use the same `ArtifactStore(...).read_condition_result("result.json")` pattern
  for other condition-level artifacts, then inspect `artifact.payload` for the
  fields you want to plot.

## Troubleshoot common notebook issues

### `ImportError: No module named polyzymd`

Start the notebook server from the pixi environment:

```bash
pixi run -e build jupyter lab
```

If you use VS Code or an existing Jupyter server, select the kernel associated
with the PolyzyMD `build` environment.

### `FileNotFoundError` for `result.json`

Check that the hydrogen-bond analysis has been run and that `project_dir` points
to the directory containing `analysis/`:

```bash
pixi run -e build polyzymd compare run hydrogen_bonds -f comparison.yaml
```

Also list the directories under `analysis/` and compare them with
`condition_dirs`:

```python
for path in sorted((project_dir / "analysis").iterdir()):
    if path.is_dir():
        print(path.name)
```

### `KeyError` for a summary name

Print the available summary names from each artifact and compare them with the
`summaries` section of `comparison.yaml`.

```python
for condition, artifact in condition_artifacts.items():
    available = [summary.get("name") for summary in artifact.payload["summaries"]]
    print(condition, available)
```

## See also

- {doc}`publication_plots` for built-in plot settings and regenerated plots.
- {doc}`hydrogen_bonds` for configuring hydrogen-bond summaries.
- {doc}`../reference/analysis_hydrogen_bonds_reference` for hydrogen-bond
  settings and documented output fields.
- {doc}`../reference/comparison_yaml` for `comparison.yaml` schema details.
- {doc}`../reference/analysis_comparison_reference` for comparison output paths
  and plotting behavior.
