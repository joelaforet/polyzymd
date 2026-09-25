# Write an analysis plugin

An analysis in PolyzyMD is one module holding a settings model and a
`compute()` that returns observables. The framework owns everything else:
loading the universe, choosing the production frames, caching the replicate,
aggregating across replicates, testing across conditions, drawing the figures
and formatting the report. A new analysis is about forty lines.

This page is the whole contributor path. The
[checklist](checklist.md) is what to run before opening the pull request.

## Where it lives

An analysis written for a study lives in the study's `analyses/` folder; see
{doc}`../../how_to/study_layout`. An analysis contributed to PolyzyMD itself
lives in `src/polyzymd/analyses/`. The contract is the same in both places.

## Scaffold it

Inside a study:

```bash
polyzymd new-analysis solvent_shell
```

That writes `analyses/solvent_shell.py` and `analyses/test_solvent_shell.py`.
In a PolyzyMD checkout, or with `--builtin`, it writes
`src/polyzymd/analyses/solvent_shell.py` and
`tests/analyses/plugins/test_solvent_shell.py` instead. Edit those two files
and nothing else. There is no registry entry and no import to add; discovery
finds the module.

## The contract

```python
from typing import Any, ClassVar, Sequence

from pydantic import BaseModel

from polyzymd.analyses import Observable, contract_analysis, iter_frames


class SolventShellSettings(BaseModel):
    selection: str = "protein"


class SolventShell:
    name: ClassVar[str] = "solvent_shell"
    Settings: ClassVar[type[BaseModel]] = SolventShellSettings
    references: ClassVar[tuple[str, ...]] = ("Author 2020, Journal 1:1, doi:10.0/x",)

    def compute(
        self, universe: Any, frames: Any, settings: SolventShellSettings
    ) -> Sequence[Observable]:
        group = universe.select_atoms(settings.selection)
        values = [measure(group) for _ in iter_frames(universe, frames)]
        return [
            Observable(
                name="solvent_shell",
                kind="mean_of_timeseries",
                unit="A",
                values=values,
            )
        ]


SolventShellAnalysis = contract_analysis(SolventShell)
```

The last line is optional. Discovery registers a module-level class that has
`name`, `Settings` and `compute()` whether or not it is wrapped, and wrapping it
yourself gives a name to import in a script.

`compute()` runs once per replicate and returns raw per-frame numbers. Never
average across replicates inside it, never run a statistical test, never write
a file and never import matplotlib. The framework does all of that, and it does
it from the `kind` of each observable.

## Pick the kind

| kind | values are | replicate value | reported as |
| --- | --- | --- | --- |
| `mean_of_timeseries` | one number per frame | mean | mean, SEM and 95 percent interval across replicates |
| `fluctuation` | one number per frame | sample standard deviation | the same, on the fluctuation |
| `fraction` | 0 or 1 per frame | occupancy | the same, on the fraction |
| `profile` | one number per index | the whole vector | per-index mean and SEM |

There is no distribution kind. Express a shape as a `profile` over histogram
bins. The raw per-frame series is kept in an NPZ sidecar either way.

The kind also picks the figures, so a plugin never writes a plotter. The
time-series kinds get a comparison bar chart and a per-replicate time-series
panel, a `fraction` gets bars on a `[0, 1]` axis, and a `profile` gets a line
per condition with a 95 percent band, or grouped bars when its index names at
most thirty categories. To change a figure, subclass `ContractPlotSettings` and
attach it to the plugin as
`PlotSettings: ClassVar[type[BasePlotSettings]] = MyPlotSettings`.

A `profile` carries no single comparable number, so give it one with
`reduce="mean_over_index"` or `"sum_over_index"`, reported as `<name>_mean` or
`<name>_total`. Add `reduced_kind="fluctuation"` or `"fraction"` when that
scalar is one of those, and set `n_frames` on the profile so the scalar counts
frames rather than indices.

When one observable is a function of others the plugin already reports, such as
the last class of a set of fractions that sums to one, declare it
`tested=False`. It is still aggregated and reported with its uncertainty, and it
stays out of the pairwise tests and out of the Benjamini-Hochberg family, so it
cannot weaken the adjusted p-values of the quantities that carry independent
information.

Every observable states a `unit`. The scaffold placeholder `"TODO"` is
rejected, so the generated tests fail until it is replaced. A `profile` also
states an `index`, one entry per value, and may set `index_label` to name the x
axis of its figure.

If the answer depends on a file the framework does not load, such as a
reference structure the settings name, define
`identity_files(settings) -> Sequence[Path]`. Those files join the replicate
identity block, so replacing one recomputes the replicate.

## Cite the method

Put the method paper in the module docstring under a NumPy `References`
heading and repeat it in `references`. The comparison artifact carries it.

## Write the tests

`polyzymd.analyses.testing` does the setup without trajectories or disk.
`synthetic_universe()` is four unit-mass atoms on a cross, whose radius of
gyration is exactly 1.0. `run_in_memory(plugin, settings, universes)` runs
`compute()` on three replicates, or on one universe per replicate if you pass a
list, and returns the aggregates the framework would report.

```python
from polyzymd.analyses.testing import run_in_memory, synthetic_universe


def test_aggregates():
    settings = SolventShellSettings()

    aggregates = run_in_memory(SolventShell, settings, synthetic_universe())

    assert aggregates[0].n_replicates == 3
    assert aggregates[0].mean == pytest.approx(1.0)
    assert aggregates[0].sem == pytest.approx(0.0)
```

Assert on the aggregate fields rather than on the raw observable: `mean`,
`sem`, `ci95_low`, `ci95_high`, `n_replicates`, `replicate_values`, and for a
profile `profile_mean`, `profile_sem` and `index`. A comparison entry carries
`delta`, `p_value`, `p_adjusted`, `significant`, `testable` and `note`. The
scaffolded test file already holds both tests; replace the numbers.

```bash
pytest analyses/test_solvent_shell.py -q                       # in a study
PYTHONPATH=$PWD/src pixi run -e test pytest tests/analyses/plugins/test_solvent_shell.py -q
pixi run -e analysis polyzymd compare run solvent_shell -f comparison.yaml
```

A plugin defined in a script or notebook can skip the file entirely:
`analyze(SolventShell, ["A/config.yaml", "B/config.yaml"], equilibration="10ns")`
registers it and runs it like any other analysis.

## What to import

A plugin imports from three places and nowhere else.

- `polyzymd.analyses` for `Observable`, `iter_frames` and `contract_analysis`.
- `polyzymd.analyses.shared` for a helper that already exists, such as
  alignment, topology checks or amino acid classification.
- `polyzymd.analyses.testing` in its tests.

`polyzymd.analyses._framework` is internal. A plugin that imports from it is
reaching past the contract, and the next change to the framework will break it.

```{toctree}
:hidden:
:maxdepth: 1

Analysis plugin contribution checklist <checklist>
```
