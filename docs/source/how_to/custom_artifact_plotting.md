# Create Custom Plots from Study Results

Use this guide when you want a figure that PolyzyMD does not draw, from values
measured with the study API. The example measures, for each protein residue,
the fraction of production frames in which it has a hydrogen bond to the
polymer, and draws chosen residues side by side for every condition, with each
condition's 95% interval and every replicate value.

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
The first run of the measurement loads the trajectories. Each replicate's
values and their record are stored under `polyzymd_results/` in the current
directory, and a later call with the same function, selections and settings
reads them back instead of measuring again. To change PolyzyMD's standard
figures instead, start with {doc}`publication_plots`.
```

## Import notebook dependencies

```python
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import polyzymd as pz
from polyzymd.analyses.functions import residue_hbond_occupancy
from polyzymd.analyses.shared.plotting import add_uncertainty_footnote
```

## Measure per-residue occupancy

Give each condition its simulation `config.yaml`, control first.
`residue_hbond_occupancy` returns one value per residue of its first group, so
`labels` names each entry by its residue ID:

```python
study = pz.Study.from_configs(
    {"No polymer": "noPoly/config.yaml", "SBMA 50%": "SBMA50/config.yaml"},
    equilibration="200ns",
)
occupancy = study.per_replicate(
    residue_hbond_occupancy,
    pz.select("chainid A"),
    pz.select("chainid C"),
    unit=None,
    labels=lambda u: [int(resid) for resid in u.select_atoms("chainid A").residues.resids],
    bounds=(0.0, 1.0),
)
```

`occupancy.values` holds each condition's replicate arrays, in the order of
`occupancy.labels`. `occupancy.compare()` tests every residue of every
condition against the control, with the Benjamini-Hochberg correction over
all of them; print `occupancy.compare().to_agent_text()` before choosing the
residues to plot.

## Extract a tidy DataFrame

`occupancy.summary()` gives one row per condition and residue, with the mean,
its Student t 95% interval and every replicate value. Keep the residues to
plot:

```python
residues = ["77", "133", "156"]
rows = [row for row in occupancy.summary().conditions if row.entry in residues]
frame = pd.DataFrame(
    {
        "condition": row.label,
        "residue": row.entry,
        "mean": row.mean,
        "ci_low": row.ci95[0] if row.ci95 else np.nan,
        "ci_high": row.ci95[1] if row.ci95 else np.nan,
        "values": row.replicate_values,
        "n": row.n_replicates,
    }
    for row in rows
)
frame
```

`entry` holds the label as text. When no interval can be estimated, for a
condition with one replicate or with the same value in every replicate,
`ci95` is `None`, so `ci_low` and `ci_high` are `NaN` and no error bar is
drawn. An interval can also reach below 0 or above 1; the `warning:` lines of
`occupancy.compare().to_agent_text()` name those residues.

## Plot grouped bars with 95% intervals and replicate dots

```python
conditions = study.labels
x = np.arange(len(residues))
width = 0.8 / len(conditions)

fig, ax = plt.subplots(figsize=(7, 4))
for i, condition in enumerate(conditions):
    part = frame[frame["condition"] == condition].set_index("residue").loc[residues]
    positions = x + (i - (len(conditions) - 1) / 2) * width
    errors = np.array([part["mean"] - part["ci_low"], part["ci_high"] - part["mean"]])
    ax.bar(positions, part["mean"], width, yerr=errors, capsize=3, label=condition, alpha=0.8)
    for position, values in zip(positions, part["values"]):
        ax.scatter(np.full(len(values), position), values, color="black", s=12, zorder=3)

ax.set_xticks(x, [f"Residue {residue}" for residue in residues])
ax.set_ylabel("Fraction of frames with a hydrogen bond to the polymer")
ax.set_ylim(0, 1.05)
ax.legend(frameon=False)
counts = frame["n"].unique()
add_uncertainty_footnote(
    fig,
    n_replicates=int(counts[0]) if len(counts) == 1 else None,
    equilibration=study[conditions[0]].equilibration,
)
fig.tight_layout()
```

`add_uncertainty_footnote` writes a sentence under the axes saying that the
error bars are the 95% Student t confidence interval of the mean across the
replicates and that the points are the per-replicate values, as on every
PolyzyMD figure; with `n_replicates=None` it says that n is per condition.
Its `drawn=` and `of=` keywords name the mark and the
quantity, for example `drawn="Band", of="the condition mean at each residue"`.

## Save outside the results folder

Keep custom figures apart from `polyzymd_results/`, so the stored values and
their records stay as PolyzyMD wrote them:

```python
from pathlib import Path

output_dir = Path("custom_figures")
output_dir.mkdir(exist_ok=True)
fig.savefig(output_dir / "hbond_occupancy_selected_residues.png", dpi=300, bbox_inches="tight")
```

## Adapt this pattern

- Any function measured with `study.per_replicate` or `study.timeseries`
  gives values of the same kind; `series.reduce("mean")` turns a time series
  into one value per replicate before `summary()`.
- For one number per replicate, such as the `mean_hbonds` row of
  `functions.hydrogen_bonds` in {doc}`hydrogen_bonds`, each summary row has
  `entry` set to `None`: plot one bar per condition.
- `occupancy.over_labels("mean", labels=[76, 132, 155])` averages the chosen
  residues into one value per replicate, which `compare()` then tests as one
  quantity.

## Troubleshoot common notebook issues

### `ImportError: No module named polyzymd`

Start the notebook from the analysis environment, or select a kernel created
from it:

```bash
pixi run -e analysis jupyter lab
```

### `KeyError` for a residue

The residues in `residues` must be labels of `occupancy`, written as text.
Print `occupancy.labels` to see them, and check that the residue is in the
first selection passed to `study.per_replicate`.

## See also

- {doc}`hydrogen_bonds` for the hydrogen-bond functions and the
  `polyzymd analyze hydrogen_bonds` command.
- {doc}`publication_plots` for settings that control PolyzyMD's standard
  analysis plots.
- {doc}`../reference/analysis_functions` for every shipped function and what
  it returns.
