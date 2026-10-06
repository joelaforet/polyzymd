# Run Rg analysis

Measure the radius of gyration (Rg) of a selection on each production frame of
each replicate. Then compare the conditions, with one value per replicate.

Rg measures how compact a selection is. It does not change when the selection
moves or rotates, so it needs no alignment and no reference structure.

:::{admonition} Environment Setup
:class: tip

All analysis commands below assume you have activated the PolyzyMD analysis
pixi environment:

```bash
pixi shell -e analysis
```

Alternatively, prefix each command with `pixi run -e analysis`.
:::

## From the command line

```bash
polyzymd analyze rg -c noPoly/config.yaml -c SBMA50/config.yaml \
  --label "No polymer" --label "SBMA 50%" --eq 200ns
```

The first `-c` is the control. The command does these steps:

1. It measures the mass-weighted Rg of the `protein` atoms on each frame after
   the equilibration window.
2. It takes the mean over frames of each replicate.
3. It compares each condition with the control by Welch's t test. It corrects
   the p values with the {term}`Benjamini-Hochberg` method.

To measure a different selection, add `--set selection='protein and name CA'`.

Useful options:

- `--format json` prints the full report.
- `--replicates 1-3` uses only some replicates.
- `--recompute` ignores stored results.

For the line format and the verdict words, see
{ref}`polyzymd analyze <cli-analyze>`.

## From Python

```python
import polyzymd as pz
from polyzymd.analyses.functions import radius_of_gyration

study = pz.Study.from_configs(
    {"No polymer": "noPoly/config.yaml", "SBMA 50%": "SBMA50/config.yaml"},
    equilibration="200ns",
)
rg = study.timeseries(radius_of_gyration, pz.select("protein"), unit="Å")
print(rg.reduce("mean").compare(control="No polymer").to_agent_text())
```

`radius_of_gyration(atoms)` calls MDAnalysis `AtomGroup.radius_of_gyration()`
on the current frame. You can use any function of an `AtomGroup` that returns
one number in its place, for example the Rg of one polymer chain. See
{doc}`study_api`.

PolyzyMD stores the per-frame values of each replicate under
`polyzymd_results/<name>/<condition>/replicate_<n>/`. A record beside the
values names the function, the selection, the input files and the frames
used. The next command reuses the values only if all of these match.

## Before you interpret the numbers

`radius_of_gyration` does not unwrap molecules that cross the periodic
boundary. Make sure that the selected atoms are whole in the trajectory. Plot
the time series of a few replicates before you choose the equilibration
window. For what Rg shows and does not show, see
{doc}`../explanation/analysis_rg_best_practices`.
