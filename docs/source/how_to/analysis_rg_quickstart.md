# Rg analysis: quick start

Measure the radius of gyration (Rg) of a selection on every production frame of
every replicate, and compare conditions with the replicate as the sampling
unit.

```{versionadded} 1.3.0
Rg analysis was added in PolyzyMD 1.3.0.
```

```{note}
This page focuses on getting results quickly. For what each shipped function
measures, see {doc}`../reference/analysis_functions`; for the study API behind
it, see {doc}`../explanation/analysis_api`.
```

:::{admonition} Environment Setup
:class: tip

All analysis commands below assume you have activated the PolyzyMD analysis
pixi environment:

```bash
pixi shell -e analysis
```

Alternatively, prefix each command with `pixi run -e analysis`.
:::

```{tip}
Rg complements RMSD and RMSF:

- **Rg** answers compactness questions
- **RMSD** answers reference-deviation questions
- **RMSF** answers per-residue flexibility questions

Rg is translation and rotation invariant, so it does not require alignment or
reference structures.
```

## From the command line

```bash
polyzymd analyze rg -c noPoly/config.yaml -c SBMA50/config.yaml \
  --label "No polymer" --label "SBMA 50%" --eq 200ns
```

The first `-c` is the control. For each replicate, the mass-weighted radius of
gyration of the `protein` atoms is measured on every frame after the
equilibration window and averaged over those frames. The replicate means are
then summarised per condition, and every other condition is compared with the
control by Welch's t test with the Benjamini-Hochberg correction. Measure a
different selection with `--set selection='protein and name CA'`.

Add `--format json` for the full report, `--replicates 1-3` to use only some
replicates, and `--recompute` to ignore stored results. The line format and the
verdict words are described under {ref}`polyzymd analyze <cli-analyze>`.

## From Python

The same analysis written with the {doc}`study API <../explanation/analysis_api>`:

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
on the current frame. Any function of an `AtomGroup` that returns one number
works in its place, for example the radius of gyration of one polymer chain.

Each replicate's per-frame values are stored under
`polyzymd_results/<name>/<condition>/replicate_<n>/` with a record of the
function, the selection, the input files and the frames used, and they are
reused on the next run only when all of those match.

## Before interpreting the numbers

Rg does not unwrap molecules split across periodic boundaries, so check that
the selected atoms are whole in the trajectory. Plot a few replicates' time
series before choosing the equilibration window.
{doc}`../explanation/analysis_rg_best_practices` covers what Rg does and does
not show, and how to read its time series.
