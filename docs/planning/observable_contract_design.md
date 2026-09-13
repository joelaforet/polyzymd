# The observable contract

A new PolyzyMD analysis costs a 337-line scaffold, nine hook overrides and two
vocabularies (`Analysis` in `base.py` and the `mda/` artifact layer). This
document defines the contract that replaces that with a settings model and one
function, and how the nine shipped plugins move onto it.

## What an Observable is

```python
class Observable(BaseModel):
    name: str                  # unique within the plugin
    kind: ObservableKind       # one of five, fixes the statistics
    values: list[float]        # per frame, or per index for a profile
    unit: str | None
    index: list[float] | None  # residue IDs or bin centres, profile only
    higher_is_better: bool | None
```

A plugin returns one `Observable` per reported quantity per replicate, holding
raw per-frame numbers. It never averages across replicates, runs a test or
writes a file. The replicate stays the sampling unit because the plugin has no
way to cross that boundary.

## The five kinds

| kind | replicate value | condition uncertainty | comparison |
|---|---|---|---|
| `mean_of_timeseries` | mean of the series | mean, SEM and Student t 95 percent interval over replicates | t-test on replicate means |
| `fluctuation` | sample standard deviation | same, on the fluctuation | same |
| `fraction` | mean of the 0/1 series | same, on the fraction | same |
| `profile` | the per-index vector | per-index mean and SEM across replicates | not tested pairwise yet |

A fifth kind, `distribution`, was dropped before merge: it reduced exactly like
`mean_of_timeseries` while the documentation promised a shape test that did not
exist. A distribution is expressed today as a `profile` over histogram bins. A
dedicated kind returns when a plugin needs a real shape statistic, with the test
that goes with it.

A `profile` may declare `reduce` as `"mean_over_index"` or `"sum_over_index"`,
which reports a second observable named for the operation, `<name>_mean` or
`<name>_total`, holding the profile reduced over its index. Its kind is the
observable's `reduced_kind`, `mean_of_timeseries` unless the plugin says
otherwise, so a profile of fluctuations declares `reduced_kind="fluctuation"`
and a profile of occupancies declares `"fraction"`. That scalar is what the
pairwise test uses, because profiles are not tested pairwise.

A fluctuation over one frame has no estimate, and aggregation raises
`PluginContractError` when fewer than two replicates remain estimable. A pair
where either side has one replicate carries `testable=False` and the note
"single replicate", as the framework's `PairwiseResult` does, and prints as "not
testable". A `control_label` naming no compared condition raises and lists the
labels that exist.

Correlation inside a replicate never shrinks an error bar. The shared
`statistical_inefficiency` gives g and N_eff per replicate as a diagnostic
(`n_eff_min`), so a reader can see a barely decorrelated replicate.

Every test in one run forms a single Benjamini-Hochberg family, across all
observables and all pairs, using `benjamini_hochberg` from
`shared/inferential_statistics.py`. With `posthoc_method: tukey_hsd` and three
or more conditions, Tukey's test adjusts family-wise per observable instead. The
t-test method (`student` or `welch`) and alpha come from the comparison config,
so whether a result is corrected no longer depends on which plugin ran.

`aggregate_observables` uses `compute_sem` from `shared/statistics.py`. The
Student t half width sits in one private helper carrying a TODO to call
`mean_sem_ci` once the confidence-intervals branch lands; there is no second
public interval estimator.

## What the framework writes to disk

Nothing new. Replicates, conditions and comparisons stay `ReplicateArtifact`,
`ConditionArtifact` and `ComparisonArtifact` from `mda/artifacts.py`, written
through `mda/store.py`. Each payload is a list of observable records instead of
a per-plugin result model.

The runner writes a `provenance.identity` block the plugin cannot omit or get
wrong: `polyzymd_version`, `plugin`, `plugin_code_hash` (SHA-256 of the plugin
module source, so a fix in a module-level helper invalidates the cache too),
`settings_fingerprint`, `config_hash`, `equilibration`, `inputs`, the topology
and trajectory `FileIdentity` records the universe provider already computes,
and `settings_files`, the `FileIdentity` of every path a plugin declares
through the optional `identity_files(settings)` hook. That hook exists because
a plugin can depend on a file the framework never loads, such as an external
reference structure, and replacing its contents must recompute the replicate. The per-frame series goes to an `observables.npz` sidecar beside
`result.json`, hashed and validated by the store, which is what later makes
generic time-series and shape plots possible without re-running the trajectory.

A replicate is reused only when every one of those fields matches, `inputs`
included, so extending a trajectory and rerunning recomputes rather than
reporting a stale number. Reading the current file identity costs one provider
call and does not load the trajectory. This closes the M4 gap, where
`compare run` recomputed every replicate every time.

## How a contract plugin is discovered

A plugin is a plain class with `name`, `Settings`, `compute` and optional
`references`, wrapped in one line:

```python
MyAnalysis = contract_analysis(My)
```

`contract_analysis` generates the `Analysis` subclass that the existing
discovery scan finds, so nothing in `discovery.py` changes yet. It checks the
plugin against the runtime-checkable `AnalysisProtocol` first and names the
attributes that are missing, so the protocol is enforced rather than described.
`Observable` rejects the scaffold's placeholder unit, so an analysis cannot
reach disk without stating what its numbers mean.

A module-level `compute` function plus a `Settings` class would save the class
statement and the `self` argument, about three lines. A class was chosen anyway
because it keeps `name`, `Settings`, `references` and `compute` in one place an
agent can read and edit without scanning a module, because a second analysis can
live in the same file, and because the wrapper holds one object rather than
three module attributes found by name. The scaffold writes the class header, so
the agent pays nothing for it.

Out of tree, the same file installed from any package registers through the
entry-point group `polyzymd.analyses`:

```toml
[project.entry-points."polyzymd.analyses"]
my = "my_package:MyAnalysis"
```

`discovery._discover_plugins` then merges the entry points of that group with
the in-tree scan, about ten lines, kept out of this change because seven
branches are in flight over that file.

## What an agent writes

```python
"""Radius of gyration, written against the observable contract.

References
----------
Michaud-Agrawal, N. et al. (2011). MDAnalysis. *J Comput Chem*, 32(10),
2319-2327. doi:10.1002/jcc.21787
"""

from typing import Any, ClassVar, Sequence

from pydantic import BaseModel, Field

from polyzymd.analyses.contract import Observable, contract_analysis, iter_frames
from polyzymd.analyses.exceptions import ReplicateError


class RgRun(BaseModel):
    label: str = Field(min_length=1)
    selection: str = Field(min_length=1)


class RgSettings(BaseModel):
    runs: list[RgRun] = Field(min_length=1)


class RgContract:
    name: ClassVar[str] = "rg2"
    Settings: ClassVar[type[BaseModel]] = RgSettings
    references: ClassVar[tuple[str, ...]] = ("Michaud-Agrawal 2011, doi:10.1002/jcc.21787",)

    def compute(self, universe: Any, frames: Any, settings: RgSettings) -> Sequence[Observable]:
        groups = {run.label: universe.select_atoms(run.selection) for run in settings.runs}
        empty = sorted(label for label, group in groups.items() if len(group) == 0)
        if empty:
            raise ReplicateError(f"rg2: selections {empty} matched no atoms")
        series: dict[str, list[float]] = {label: [] for label in groups}
        for _ in iter_frames(universe, frames):
            for label, group in groups.items():
                series[label].append(float(group.radius_of_gyration()))
        return [
            Observable(name=label, kind="mean_of_timeseries", unit="A", values=values)
            for label, values in series.items()
        ]


Rg2Analysis = contract_analysis(RgContract)
```

That is the whole of `rg_contract/`: 84 lines with its docstrings, 26 lines of
code, against 3,881 lines in `rg/`.

## What an agent runs

```bash
pixi run -e analysis polyzymd new-analysis my_metric --style contract
pixi run -e analysis polyzymd compare run my_metric -f comparison.yaml
```

The scaffold emits a 59-line plugin (25 lines of code) and a 50-line test with
two known answers, both failing until the author replaces the placeholder unit.
The test uses two fixtures from `tests/analyses/conftest.py`,
`synthetic_universe` and `run_contract_analysis`, so the end-to-end check over
three replicates is five lines and needs no trajectory. Reading the 74-line `polyzymd-extend`
skill costs about 800 tokens and the two generated files about 800 more, so a
new analysis costs an agent roughly 1,600 tokens of reading and a few hundred of
editing. The current scaffold is 337 plus 277 lines, about 4,900 tokens to read
before a line is changed.

## What rg2 does not cover yet

Fragment mode (`calculation_mode: fragments`, mass or equal weighting, the
single-fragment and missing-bond fallbacks) is about 220 lines in `rg/_mda.py`
and returns as one `mean_of_timeseries` for the reduction plus a `profile` over
histogram bins for the spread, roughly 40 lines. Fragment Rg histograms (about
180 lines of aggregation and sidecar code) become framework profile aggregation,
no plugin code. The three plot families (`plot_rg_timeseries`,
`plot_rg_comparison_bars`, `plot_rg_distributions`, 790 lines) become the
generic plotters keyed on kind, phase 2 of the review, about 300 lines shared by
every plugin. An empty selection is one raise here instead of a
`RgSkippedRunPayload` ladder.

## Porting, one plugin per pull request

Order: rmsf, rg, rmsd, sasa, distances, secondary_structure, catalytic_triad,
hydrogen_bonds, and contacts last because residence times need an optional
`Aggregator` protocol. Each pull request rewrites one plugin against the contract, keeps its name, and
deletes that plugin's `compare()`, `aggregate()`, `extract_metrics()`,
`_comparison_results.py`, `_formatters.py`, `_plot_settings.py` and its copies of
the five helpers the review lists as duplicated. Its tests move to known answers
on the observable.

Measured deletion targets: `_comparison_results.py` across five plugins 1,427
lines, `_formatters.py` 1,875, `_plotters.py` with `_plot_settings.py` about
3,900, per-plugin `compare()` 832. The framework then loses
`_framework/compare.py` (152), most of `comparison_models.py` (252) and
`mda/plugin.py` (322), `_framework/contract.py` (111) and the `__module__`
rewriting, with `mda/comparison.py` (733) and much of `stats.py` (1,291)
replaced by `compare_observables`. Between 10,000 and 12,000 lines out against the
1,351 added here, 674 of which are code and the rest docstrings.

## Why there are three models, not two

`ObservableEstimate` was kept rather than folded into its neighbours. It carries
`value`, `n_frames`, the statistical inefficiency and `n_eff`, none of which
exist on `Observable` (raw per-frame values) or `ObservableAggregate` (across
replicates). Merging it would give one model that is sometimes raw and sometimes
reduced, with half its fields empty in each state, and the persisted replicate
payload would lose its type.

## One hook the base class lacks

Replicate cache reuse is implemented by overriding the private
`Analysis._run_compute_stage`. The decision belongs to the lifecycle, not the
plugin: `AnalysisLifecycle.run_replicate_once` should ask the framework whether
a stored artifact's identity block still matches before calling the compute
stage at all, which also gives `--recompute` and the HPC worker path one owner.
That change touches `_framework/lifecycle.py` and waits for the cache-freshness
branch.
