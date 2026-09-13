---
name: polyzymd-extend
description: Add a new analysis to PolyzyMD in about 40 lines using the observable contract. Use when an analysis does not exist yet and you need to measure something new from trajectories, when asked to "add an analysis", "extend PolyzyMD", or "write a plugin". Do not copy an existing plugin package, and do not write persistence, statistics or plotting code.
---

# Add an analysis to PolyzyMD

## 1. Scaffold it

```bash
pixi run -e analysis polyzymd new-analysis <name> --style contract
```

This writes `src/polyzymd/analyses/<name>.py` and
`tests/analyses/plugins/test_<name>.py`. Edit those two files and nothing else.

## 2. The contract

```python
class MySettings(BaseModel):
    selection: str = "protein"

class My:
    name: ClassVar[str] = "my"
    Settings: ClassVar[type[BaseModel]] = MySettings
    references: ClassVar[tuple[str, ...]] = ("Author 2020, Journal 1:1, doi:...",)

    def compute(self, universe, frames, settings) -> Sequence[Observable]:
        group = universe.select_atoms(settings.selection)
        values = [measure(group) for _ in iter_frames(universe, frames)]
        return [Observable(name="my", kind="mean_of_timeseries", unit="A", values=values)]

MyAnalysis = contract_analysis(My)   # keep this line; discovery needs it
```

`compute` runs once per replicate and returns raw per-frame numbers. Never
average across replicates, test anything, write a file, or import matplotlib.
The framework does all of that from `kind`.

## 3. Pick the kind

| kind | values are | replicate value | reported as |
|---|---|---|---|
| `mean_of_timeseries` | one number per frame | mean | mean, SEM, 95 percent CI over replicates |
| `fluctuation` | one number per frame | sample standard deviation | same, on the fluctuation |
| `fraction` | 0 or 1 per frame | occupancy | same, on the fraction |
| `profile` | one number per index | the vector | per-index mean and SEM, no pairwise test yet |

There is no distribution kind. Express a shape as a `profile` over histogram
bins; the raw per-frame series is kept in the NPZ sidecar either way.

The kind also picks the figures, so plots come free and you never write a
plotter: the time-series kinds get a comparison bar chart and a per-frame panel,
a `fraction` gets bars on a `[0, 1]` axis, and a `profile` gets a line per
condition with a 95 percent band, or grouped bars when its index names at most
30 categories. To change a figure, subclass `ContractPlotSettings` and attach it
to the plugin as `PlotSettings: ClassVar[type[BasePlotSettings]] = MyPlotSettings`;
leave it off and the plugin uses the default, whose fields are `error_bar`,
`figsize`, `show_replicates` and `max_categories_for_bars`.

If one observable is a function of others you already report, such as the last
class of a set of fractions that sums to one, declare it `tested=False`. It is
still aggregated and reported with its uncertainty, but it stays out of the
pairwise tests and out of the Benjamini-Hochberg family, so it cannot weaken the
adjusted p-values of the quantities that carry independent information.
A profile is not tested pairwise, so give it a comparable scalar with
`reduce="mean_over_index"` or `"sum_over_index"`, reported as `<name>_mean` or
`<name>_total`. Add `reduced_kind="fluctuation"` or `"fraction"` when that
scalar is one, and set `n_frames` on the profile so the scalar counts frames
rather than indices. If your answer depends on a file the framework does not
load, such as a reference structure your settings name, add
`identity_files(settings) -> Sequence[Path]` so replacing that file recomputes
the replicate.

Every observable states a `unit`, and the scaffold placeholder `"TODO"` is
rejected, so the generated tests fail until you replace it. A `profile` also
states an `index`, one entry per value (residue IDs, bin centres), and may set
`index_label` to name the x axis of its figure.

## 4. Citations

Put the method paper in the module docstring under a NumPy `References`
heading and repeat it in `references`; the comparison artifact carries it.

## 5. The one test to write

Two fixtures in `tests/analyses/conftest.py` do the setup: `synthetic_universe`
is four unit-mass atoms on a cross, whose radius of gyration is exactly 1.0, and
`run_contract_analysis(AnalysisClass, settings, universe)` runs the real
lifecycle over three replicates and returns the `ConditionArtifact`.

```python
def test_aggregates(synthetic_universe, run_contract_analysis):
    artifact = run_contract_analysis(MyAnalysis, MySettings(), synthetic_universe)
    aggregate = ObservableAggregate.model_validate(artifact.payload["observables"][0])
    assert aggregate.n_replicates == 3
    assert aggregate.replicate_values == pytest.approx([1.0, 1.0, 1.0])
    assert aggregate.mean == pytest.approx(1.0) and aggregate.sem == pytest.approx(0.0)
```

Assert on the aggregate fields, not on the raw observable: `mean`, `sem`,
`ci95_low`, `ci95_high`, `n_replicates`, `replicate_values`, and for a profile
`profile_mean`, `profile_sem` and `index`. A comparison entry carries `delta`,
`p_value`, `p_adjusted`, `significant`, `testable` and `note`. The scaffolded
test has both tests already; replace the numbers.

```bash
PYTHONPATH=$PWD/src pixi run -e test pytest tests/analyses/plugins/test_<name>.py -q
pixi run -e analysis polyzymd compare run <name> -f comparison.yaml
```

## 6. Check the science

Run the `livecoms-check` skill before opening a pull request. It fails a number
reported without a unit, an uncertainty that does not say what it is, a test on
frames instead of replicates, and an uncited method. The contract handles the
first three as long as you do not aggregate inside `compute`.
