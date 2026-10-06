# Measure your own quantity on every replicate

Write a Python function that measures one frame or one replicate. PolyzyMD
runs it on the production frames of every replicate of every condition,
stores each result with a record of its inputs, and compares the conditions
with one value per replicate.

For the full signatures, see {doc}`../reference/study_api`. For how the
statistics work and why, see {doc}`../explanation/analysis_api`. To run your
function from `study.yaml` instead of Python, see {doc}`study_yaml`.

:::{admonition} Environment Setup
:class: tip

The examples on this page assume you have activated the PolyzyMD analysis
pixi environment, which provides MDAnalysis and the statistics dependencies:

```bash
pixi shell -e analysis
```

Alternatively, prefix each command with `pixi run -e analysis`, for example
`pixi run -e analysis python my_analysis.py`.
:::

## 1. Load the study

Give the config of each condition, the control first, and the equilibration
window:

```python
import polyzymd as pz

study = pz.Study.from_configs(
    {"No polymer": "noPoly/config.yaml", "SBMA 50%": "SBMA50/config.yaml"},
    equilibration="100ns",
)
```

PolyzyMD finds the replicates of each condition from its config. It joins the
production segments of each replicate in order, and removes the first 100 ns.

To load a study folder, use `pz.Study("path/to/study.yaml")`. The
`study.yaml` gives the conditions and the window.

To measure every fifth production frame, for an expensive measurement, add
`stride=5`:

```python
study = pz.Study.from_configs(configs, equilibration="100ns", stride=5)
```

Each replicate gives its `Universe` and the frames that PolyzyMD analyzes:

```python
for replicate in study["SBMA 50%"].replicates:
    u = replicate.universe()    # the production segments, in order
    frames = replicate.frames   # trajectory frame indices after the window
    times = replicate.times     # their times, in ns
```

## 2. Measure one value per frame

Write a function of one or more `AtomGroup`s that returns one number for the
current frame. Pass it to `study.timeseries`:

```python
def radius_of_gyration(atoms):
    return float(atoms.radius_of_gyration())

rg = study.timeseries(radius_of_gyration, pz.select("protein"), unit="Å", bounds=(0.0, None))
```

- `pz.select("protein")` stands for the `AtomGroup` of that selection in each
  replicate. Several `pz.select` arguments reach the function in the order
  that you give them.
- `pz.universe()` stands for the `Universe` of each replicate.
- Other arguments and keyword arguments go to the function unchanged.
  PolyzyMD records them.
- `unit` is the unit of the value. Give `unit=None` for a dimensionless
  quantity.
- `bounds` gives the lowest and highest possible value. Use `None` for no
  limit.

The function must return one number per frame. For one value per residue,
use `study.per_replicate` (step 5).

To measure distances between two atoms, use the shipped `pair_distance`:

```python
from polyzymd.analyses.functions import pair_distance

distance = study.timeseries(
    pair_distance,
    pz.select("protein and resid 77 and name OG"),
    pz.select("protein and resid 156 and name NE2"),
    unit="Å",
    bounds=(0.0, None),
)
```

## 3. Measure against a reference structure

Pass `pz.reference(mode, selection, ...)` where the function takes the
reference atoms:

```python
from polyzymd.analyses.functions import rmsd

ca = "protein and name CA"
deviation = study.timeseries(
    rmsd, pz.select(ca), pz.reference("external", ca, file="structures/1ISP.pdb"),
    unit="Å", bounds=(0.0, None),
)
```

| Mode | The reference is |
|---|---|
| `external` | The `selection` atoms of `file` |
| `frame` | Production frame `frame`, counted from 1 after the window |
| `average` | The mean positions, after superposition on the `alignment` atoms |
| `centroid` | The production frame closest to the average structure of the `alignment` atoms |

PolyzyMD builds the reference once per replicate, in a separate universe. The
trajectory does not change.

## 4. Turn each replicate into one value

Reduce each time series to one value per replicate:

```python
mean_rg = rg.reduce("mean")
```

`reduce` takes `"mean"`, `"fraction"` (for a series of 0 and 1) or `"std"`.
It also takes a function of the values and times (in ns) of one replicate
that returns one number:

```python
def final_quarter_mean(values, times):
    return float(values[len(values) * 3 // 4:].mean())

late_rg = rg.reduce(final_quarter_mean, unit="Å")
```

To compute a new series from stored series, use `transform`. It reads no
trajectory:

```python
below = distance.transform(lambda d, cutoff: d < cutoff, cutoff=4.0, unit=None)
fraction_below = below.reduce("fraction")
```

Give values such as the cutoff as keyword arguments, as above. PolyzyMD then
records them. A value that the function reads from an enclosing scope is not
recorded.

## 5. Measure one value per replicate

When the function needs all frames at once, use `study.per_replicate`.
PolyzyMD calls the function once per replicate, with the keyword argument
`frames`, the production frame indices:

```python
from MDAnalysis.analysis.hydrogenbonds import HydrogenBondAnalysis


def mean_hbonds(u, frames):
    """Mean number of protein-polymer hydrogen bonds per frame."""
    groups = "protein or resname SBM"
    analysis = HydrogenBondAnalysis(u, between=["protein", "resname SBM"])
    analysis.hydrogens_sel = analysis.guess_hydrogens(groups)
    analysis.acceptors_sel = analysis.guess_acceptors(groups)
    analysis.run(frames=frames)
    return len(analysis.results.hbonds) / len(frames)


hbonds = study.per_replicate(mean_hbonds, pz.universe(), unit=None)
```

This form runs any MDAnalysis analysis class on the production frames. For
protein-polymer hydrogen bonds, the shipped `polyzymd analyze hydrogen_bonds`
selects the donors and acceptors for you. See {doc}`hydrogen_bonds`.

A function can return one value per label, such as one value per residue.
Name the labels with `labels=`, a list or a function of the `Universe`:

```python
from polyzymd.analyses.functions import MS_PARTS, RMS_PARTS, rms_decomposition

ca = "protein and name CA"
rows = study.per_replicate(
    rms_decomposition, pz.select(ca), pz.select(ca), pz.reference("average", ca),
    unit="Å", labels=lambda u: u.select_atoms(ca).residues.resids, parts=RMS_PARTS + MS_PARTS,
)
rmsf = rows["rmsf"]                                                 # one value per residue
mean_rmsf = rmsf.over_labels("mean")                                # one value per replicate
loop_rmsf = rmsf.over_labels("mean", "loop_rmsf", labels=range(40, 51))
```

- PolyzyMD lines up the replicates by label, never by position.
- A label that one replicate lacks is an error. To fill it, give `missing=`
  with a value.
- `parts=` names the rows of a function that returns several quantities in
  one pass. `rows[part]` is the result of each row.
- `over_labels(how, metric, labels)` turns each profile into one number,
  over every label or over the labels that you give.

## 6. Compare the conditions

```python
print(mean_rg.summary().to_agent_text())                    # each condition
print(mean_rg.compare(control="No polymer").to_agent_text())  # each condition against the control
```

- `summary()` gives, for each condition, n, the mean, the standard error, the
  95 % confidence interval and every replicate value.
- `compare()` gives, for each condition against the control, the difference
  of the means with its 95 % interval, the p value, the
  {term}`Benjamini-Hochberg` adjusted p value, and Cohen's d and Hedges' g.
- `compare(test="student")` uses Student's t test instead of Welch's t test.
- `conditions=["SBMA 50%", "EGMA 50%"]` limits a summary or a comparison to
  some conditions.
- For a labelled result, each label is one test. `untested=["coil"]` leaves a
  label out of the tests, for example the last of a set of fractions that sum
  to 1.

To report several results, print each report:

```python
for result in (mean_rg, fraction_below, mean_rmsf):
    print(result.compare(control="No polymer").to_agent_text())
```

## 7. Draw figures

```python
rg.plot()                            # the series of each replicate against time
rg.plot_distribution(threshold=None) # the distribution of frame values of each condition
mean_rg.plot()                       # the mean of each condition, with every replicate value
pz.plot_values([fraction_below, hbonds], labels=["Ser-His below 4 Å", "Hydrogen bonds"])
```

The figures read the stored results. They go to `figures/` beside
`polyzymd_results/`, or to `output_dir=`. `pz.plot_values` and
`pz.plot_distributions` put several results that share a unit in one figure.
Results with different units are refused.

## 8. Find the stored results

PolyzyMD stores each replicate in
`polyzymd_results/<name>/<condition>/replicate_<n>/`, in the current folder
or in `output_dir=`. The name is the name of the function, unless you give
`name=`. The next call with the same function, arguments and input files
reads the stored values and does not read the trajectory. `recompute=True`
measures every replicate again.

A function defined in a notebook cell, or a lambda, has no stable source
file. PolyzyMD then hashes its bytecode and prints a warning. Put the
functions of a study in a `.py` file.

## Your own statistics

The tests of `compare()` use one value per replicate. For a test that
`compare()` does not make, read that table and use your own code:

```python
table = pz.Study("my_study/study.yaml").replicate_table("lid_opening")
```

`replicate_table` reads stored results of a study folder. See
{doc}`study_yaml` and {doc}`project`.
