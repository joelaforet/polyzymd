# Analysing a set of simulations

PolyzyMD runs your analysis code on every replicate of every condition in a
study and turns the results into comparisons you can report. You write the
science, as a function of an MDAnalysis `AtomGroup` or `Universe`. PolyzyMD
finds the replicates from each `config.yaml`, loads each one as a `Universe`,
removes the equilibration window, runs your function over the production
frames, stores every per-replicate result with a record of how it was made,
and computes intervals and tests with the replicate as the sampling unit.

That interface is described below. The shipped analyses, such as radius of
gyration or hydrogen-bond occupancy, are ordinary functions written against the
same interface, so anything they do, your own code can do too.

## What you write and what PolyzyMD does

| You write | PolyzyMD does |
|---|---|
| A function that measures one frame, or an MDAnalysis analysis class you already use | Loads each replicate as a `Universe` with its production segments in order and its equilibration window removed |
| The atom selections your function needs, as selection strings | Builds the `AtomGroup` for each replicate's `Universe` and runs your function on the production frames |
| Optionally, how one replicate's frames become one value | Stores every per-frame and per-replicate result, keyed by the code, arguments and input files that produced it, and reuses it only when all of those match |
| Optionally, your own aggregation or comparison | Computes the mean, standard error and 95 percent Student t interval across replicates, runs the tests between conditions, and corrects for multiple comparisons |

## Load the replicates

```python
import polyzymd as pz

study = pz.Study.from_configs(
    {"No polymer": "noPoly/config.yaml", "SBMA 50%": "SBMA50/config.yaml"},
    equilibration="100ns",
)
```

A study folder marked by `study.yaml` loads the same way with
`pz.Study("path/to/study")`. The first condition is the control unless you name
another one when you compare. The equilibration window is measured in
simulation time from the start of each replicate's production trajectory, with
its segments joined in order.

Every replicate gives you its `Universe` and the frames PolyzyMD will analyse:

```python
for replicate in study["SBMA 50%"].replicates:
    u = replicate.universe()      # production segments in order
    frames = replicate.frames     # frame indices after the equilibration window
```

You can use these universes for anything, including analyses PolyzyMD does not
support, such as a principal component basis fitted on all replicates at once.

## Measure something on every frame

This follows MDAnalysis
[`AnalysisFromFunction`](https://docs.mdanalysis.org/stable/documentation_pages/analysis/base.html):
your function receives `AtomGroup` arguments positioned at the current frame
and returns a number or a one-dimensional array.

```python
def radius_of_gyration(atoms):
    return atoms.radius_of_gyration()

rg = study.timeseries(radius_of_gyration, pz.select("protein"), unit="Å")
```

`pz.select("protein")` stands for an `AtomGroup`. The string is an ordinary
MDAnalysis selection string. An `AtomGroup` belongs to one `Universe`, so
PolyzyMD builds it again for each replicate from the selection string, and
records the string. Several `pz.select` arguments reach your function in the
order you give them. Use `pz.universe()` where your function takes
the `Universe` itself. Every other argument is passed to your function
unchanged and recorded.

PolyzyMD runs `AnalysisFromFunction` on each replicate with that replicate's
production frames, so your function never sees an equilibration frame. Pass
`step=5` to measure every fifth production frame, for an expensive measurement
such as solvent-accessible surface area; the step is recorded with the result.
The result holds, for each replicate, the per-frame values, the frame indices
and the simulation times. Give `unit=None` for a dimensionless quantity.

A function that returns an array gives one value per frame for each entry of
the array. Name the entries with `labels`, either a list or a function of the
`Universe`:

```python
def per_residue_sasa(u): ...

sasa = study.timeseries(
    per_residue_sasa, pz.universe(),
    unit="Å^2", labels=lambda u: u.select_atoms("protein").residues.resids,
)
```

Labels are how replicates are lined up. Values are aggregated per label, never
per position, so residue 57 in one replicate is always compared with residue
57 in another.

## Use an MDAnalysis analysis you already run

Many analyses already exist as MDAnalysis `AnalysisBase` classes, in
MDAnalysis itself or in an MDAKit. Pass the class, its arguments, and a
function that reads the result you want from the finished analysis:

```python
from MDAnalysis.analysis.hydrogenbonds import HydrogenBondAnalysis

hbonds = study.run(
    HydrogenBondAnalysis, pz.universe(),
    donors_sel="protein", acceptors_sel="resname SBM",
    timeseries=lambda analysis: analysis.count_by_time(),
)
```

PolyzyMD constructs the class for each replicate, calls
`run(frames=replicate.frames)`, and passes the finished object to your
function. Use `timeseries=` when the function returns one value per frame, and
`value=` when it returns one value, or one labelled array, for the whole
replicate.

## Turn each replicate into one value

Comparisons between conditions need one value per replicate, or one value per
label. A time series becomes one with `reduce`:

```python
mean_rg = rg.reduce("mean")
```

The named reductions are `"mean"`, `"fraction"` (the mean of a series of 0 and
1, rejected for any other values) and `"std"` (the sample standard deviation
over frames). A stored time series can be transformed frame by frame before it
is reduced, without measuring the trajectory again. For example, the fraction
of frames in which a distance is below 4 Å comes from the stored distances:

```python
distance = study.timeseries(atom_distance, pz.select("resid 77 and name OG"),
                            pz.select("resid 156 and name NE2"), unit="Å")
mean_distance = distance.reduce("mean")
below_4 = distance.transform(lambda d: d < 4.0, unit=None).reduce("fraction")
```

The transform is recorded with the result like any other function.

Besides the named reductions, any function that takes the per-frame values and
times of one replicate and returns a number or a labelled array works as a
reduction, for example a residence time or the slope of a mean squared
displacement:

```python
def mean_lifetime(values, times): ...

lifetime = contact.reduce(mean_lifetime, unit="ns", bounds=(0.0, None))
```

`reduce`, `per_replicate` and the `value=` form of `study.run` all take
`bounds=(low, high)` for a quantity with a physical limit; use `None` for a
side with no limit. A reduction named `"fraction"` has the bounds 0 and 1.

When you want to compute the per-replicate value yourself from the `Universe`,
use `study.per_replicate`. Your function receives the `Universe` and the
production frames and returns a number or a labelled array:

```python
def my_quantity(u, frames): ...

values = study.per_replicate(my_quantity, pz.universe(), unit="kcal/mol")
```

PolyzyMD cannot check how such a value was computed. It records the function
and its arguments like any other result.

## Compare conditions

```python
summary = mean_rg.summary()                     # one row per condition
report = mean_rg.compare(control="No polymer")  # every condition against the control
print(report.to_agent_text())

polymers = contacts.compare(control="SBMA 50%", conditions=["SBMA 50%", "EGMA 50%"])
print(pz.report(mean_distance, below_4, polymers).to_agent_text())
```

`conditions=` limits a summary or a comparison to some of the study's
conditions, for example when the control has no polymer and the quantity
concerns the polymer. `summary()` and `compare()` results both print with
`to_agent_text()`, and `pz.report` joins several of them into one text block
and one JSON document, each result under its own name. A result is named after
its function unless you pass `name=`.

`summary()` gives, for each condition, the number of replicates, the mean, the
standard error and the 95 percent interval, along with every replicate value.
`compare()` gives, for each condition against the control, the difference in
means with its 95 percent interval, the p value, the p value adjusted for
multiple comparisons, and Cohen's d. For a labelled result there is one row per
label, and every label is one test in the correction. Leave a label out with
`untested=["coil"]` when its value is fixed by the others, for example the last
of a set of fractions that sum to one. It is still summarised with its interval
but takes no part in the tests or the correction. A label with the same value
in every replicate of both conditions, such as a residue that is never in
contact, has no variance to test; it is reported as not testable and is left
out of the correction automatically.

The statistics follow Grossfield et al. (2018), the LiveCoMS best practices for
quantifying uncertainty in molecular simulations:

- The replicate is the sampling unit. Frames inside one replicate are
  correlated, so their number never shrinks an interval. The statistical
  inefficiency and effective sample size of each replicate's time series,
  computed with `pymbar.timeseries`, are reported as diagnostics.
- The equilibration window is the one you set, applied to every replicate of
  every condition. For each time series, the start of the equilibrated region
  found by `pymbar.timeseries.detect_equilibration` (Chodera 2016) is also
  reported. A warning is added when that point falls after your window, which
  suggests the window is too short for that replicate. The window you set is
  never changed automatically.
- The 95 percent interval is the mean plus or minus the Student t coverage
  factor for n replicates times the standard error.
- Tests are Welch's t test by default, or Student's t test, with the
  Benjamini-Hochberg correction across every test in one comparison.
- Grossfield et al. recommend plotting every point when there are fewer than
  10 independent measurements, so every figure shows each replicate value next to
  the mean and interval.
- Grossfield et al. point out that a quantity with a strict upper or lower
  limit is not Gaussian, and that a t interval can then extend past the limit.
  They recommend bootstrapping for such quantities, but also caution that
  bootstrap intervals are unreliable for small samples, which is the usual
  case of three to five replicates. PolyzyMD therefore keeps the t interval
  and adds a warning when
  the interval extends past a limit, using the `bounds` of the result. When
  every replicate has
  the same value, the interval is reported as not estimable rather than as a
  zero-width interval.

## Write your own aggregation or comparison

A custom aggregation or comparison can work from the stored results, which is
fast because nothing is recomputed, or from the universes, when it needs the
structures:

```python
def my_comparison(results): ...          # stored results, keyed by condition

report = mean_rg.compare(control="No polymer", method=my_comparison)

def my_structural_comparison(study): ...  # replicates and their universes

report = study.compare(my_structural_comparison)
```

Either way the function and its arguments are recorded with the report.

## What is recorded

Every per-replicate result is stored under the study's results folder with:

- the name, module and source-code hash of the function, and of the reduction
  if there is one;
- every argument, including the selection strings;
- the simulation config hash, and the size and hash of every topology and
  trajectory file read;
- the equilibration window, the frame indices and the simulation times used;
- the PolyzyMD, MDAnalysis, NumPy and Python versions;
- the unit and the labels.

A stored result is reused only when every one of these matches the new
request. Changing the function, an argument or an input file, or extending a
trajectory, recomputes the affected replicates. `recompute=True` forces a
recomputation. A function defined in a notebook cell or as a lambda has no
stable source file, so its code is hashed from its bytecode and a warning says
that the record identifies it less reliably.

## What is out of scope

- An analysis that needs all replicates at once inside the measurement, such
  as a principal component basis, clustering or a Markov state model. Load the
  universes with `replicate.universe()` and fit such models yourself.
- A quantity defined between two conditions rather than for one replicate, such
  as a divergence between two distributions. Write it as a custom comparison.
- Intervals that are exact for a non-Gaussian quantity from three to five
  replicates. No method gives these, so an interval that extends past a limit
  is flagged with a warning instead.

## Shipped analyses

The analyses PolyzyMD ships are functions in `polyzymd.analyses` written
against this interface, and `polyzymd analyze NAME` runs them from the command
line. Each one calls MDAnalysis or MDTraj for the per-frame measurement. The
list, with the measurement each one makes and the reduction it uses, is in
{doc}`../reference/analysis_functions`.

## References

Grossfield, A.; Patrone, P. N.; Roe, D. R.; Schultz, A. J.; Siderius, D. W.;
Zuckerman, D. M. Best Practices for Quantifying Uncertainty and Sampling
Quality in Molecular Simulations. *Living J. Comput. Mol. Sci.* **2018**, 1 (1),
5067. doi:10.33011/livecoms.1.1.5067

Chodera, J. D. A Simple Method for Automated Equilibration Detection in
Molecular Simulations. *J. Chem. Theory Comput.* **2016**, 12 (4), 1799-1805.
doi:10.1021/acs.jctc.5b00784
