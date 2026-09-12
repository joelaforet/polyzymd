# The observable contract

Today a new PolyzyMD analysis costs a 337-line scaffold, nine hook overrides and
two vocabularies (`Analysis` in `base.py` and the `mda/` artifact layer). This
document defines the contract that replaces all of it with a settings model and
one function, and sets out how the nine shipped plugins move onto it.

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
raw per-frame numbers. It never averages across replicates, never runs a test,
never writes a file. The replicate stays the sampling unit because the plugin
has no way to cross that boundary.

## The five kinds

| kind | replicate value | condition uncertainty | comparison |
|---|---|---|---|
| `mean_of_timeseries` | mean of the series | mean, SEM and Student t 95 percent interval over replicates | t-test on replicate means |
| `fluctuation` | sample standard deviation | same, on the fluctuation | same |
| `fraction` | mean of the 0/1 series | same, on the fraction | same |
| `distribution` | mean, with the full series kept in the NPZ sidecar | same | same, shape comparison pending |
| `profile` | the per-index vector | per-index mean and SEM across replicates | not tested pairwise yet |

Correlation inside a replicate never shrinks an error bar. The shared
`statistical_inefficiency` gives g and N_eff per replicate, as a diagnostic
(`n_eff_min`) so a reader can see when one replicate is barely decorrelated.

Every test in one run forms a single Benjamini-Hochberg family, across all
observables and all pairs, using `benjamini_hochberg` from
`shared/inferential_statistics.py`. With `posthoc_method: tukey_hsd` and three
or more conditions, Tukey's test adjusts family-wise per observable instead. The
t-test method (`student` or `welch`) and alpha come from the comparison config,
so whether a result is corrected no longer depends on which plugin ran.

`aggregate_observables` uses `compute_sem` from `shared/statistics.py`. The
Student t half width is computed in one private helper carrying a TODO to call
`mean_sem_ci` once the confidence-intervals branch lands; there is deliberately
no second public interval estimator.

## What the framework writes to disk

Nothing new. Replicates are `ReplicateArtifact`, conditions are
`ConditionArtifact`, comparisons are `ComparisonArtifact`, all from
`mda/artifacts.py`, all written through `mda/store.py`. The payload of each is
a list of observable records instead of a per-plugin result model.

The runner writes a `provenance.identity` block the plugin cannot omit or get
wrong: `polyzymd_version`, `plugin`, `plugin_code_hash` (SHA-256 of the plugin
class source, so a fixed bug invalidates the cache), `settings_fingerprint`,
`config_hash` and `equilibration`, plus `inputs`, the topology and trajectory
`FileIdentity` records the universe provider already computes. The per-frame
series goes to an `observables.npz` sidecar beside `result.json`, hashed and
validated by the store, which is what later makes generic time-series and
distribution plots possible without re-running the trajectory.

A replicate is reused when its stored identity matches on version, code hash,
settings fingerprint, config hash and equilibration. This closes the gap the
review names in M4, where `compare run` recomputed every replicate every time.

## How a contract plugin is discovered

A plugin is a plain class with `name`, `Settings`, `compute` and optional
`references`, wrapped in one line:

```python
MyAnalysis = contract_analysis(My)
```

`contract_analysis` generates the `Analysis` subclass that the existing
discovery scan finds, so nothing in `discovery.py` changes yet.

A module-level `compute` function plus a `Settings` class would save the class
statement and the `self` argument, about three lines. A class was chosen anyway
because it keeps `name`, `Settings`, `references` and `compute` in one place an
agent can read and edit without scanning a module, because a second analysis can
live in the same file, and because the wrapper needs one object to hold rather
than three module attributes to find by name. The token cost is the same for
the agent: the scaffold writes the class header.

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

from polyzymd.analyses.contract import Observable, iter_frames
from polyzymd.analyses.contract_runner import contract_analysis
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

The scaffold emits a 53-line plugin (20 lines of code, 1.6 kB) and a 44-line
test with two known answers (1.7 kB). Reading the 74-line `polyzymd-extend`
skill costs about 800 tokens and the two generated files about 800 more, so a
new analysis costs an agent roughly 1,600 tokens of reading and a few hundred of
editing. The current scaffold is 337 plus 277 lines, about 4,900 tokens to read
before a line is changed.

## What rg2 does not cover yet

Fragment mode (`calculation_mode: fragments`, mass or equal weighting, the
single-fragment and missing-bond fallbacks) is about 220 lines in `rg/_mda.py`
and returns as a `distribution` observable per fragment plus one
`mean_of_timeseries` for the reduction, roughly 40 lines. Fragment Rg histograms
(about 180 lines of aggregation and sidecar code) become a framework
`distribution` aggregation, no plugin code. The three plot families
(`plot_rg_timeseries`, `plot_rg_comparison_bars`, `plot_rg_distributions`, 790
lines) become the generic plotters keyed on kind, phase 2 of the review, about
300 lines shared by every plugin. Per-selection skip handling for an empty
selection is one raise here instead of a `RgSkippedRunPayload` ladder.

## Porting, one plugin per pull request

Order: rmsf, rg, rmsd, sasa, then distances, then secondary_structure,
catalytic_triad, hydrogen_bonds, and contacts last because residence times need
an optional `Aggregator` protocol. Each pull request rewrites one plugin against
the contract, keeps its name, and deletes the plugin's `compare()`, `aggregate()`,
`extract_metrics()`, `_comparison_results.py`, `_formatters.py`,
`_plot_settings.py` and its copy of `_apply_fdr_correction`,
`_coerce_and_validate_aggregated_result`, `mdanalysis_version`,
`_combined_warnings` and `_effective_timestep_ps`. Its tests move to known
answers on the observable, not on the plugin's own result model.

Measured targets for the deletions: `_comparison_results.py` across five plugins
1,427 lines, `_formatters.py` 1,875, `_plotters.py` and `_plot_settings.py`
about 3,900, per-plugin `compare()` methods 832 (rmsd 174, rg 185, sasa 222,
distances 251). After the ports, the framework layer loses
`_framework/compare.py` (152), most of `_framework/comparison_models.py` (252),
most of `mda/plugin.py` (322), `_framework/contract.py` (111) and the
`__module__` rewriting in `base.py`, with `mda/comparison.py` (733) and much of
`stats.py` (1,291) replaced by `compare_observables`. That is between 10,000 and
12,000 lines out against the 1,235 added here.

## One hook the base class lacks

Replicate cache reuse is implemented by overriding the private
`Analysis._run_compute_stage`. The decision belongs to the lifecycle, not the
plugin: `AnalysisLifecycle.run_replicate_once` should ask the framework whether
a stored artifact's identity block still matches before calling the compute
stage at all, which also gives `--recompute` and the HPC worker path one owner.
That change touches `_framework/lifecycle.py` and waits for the cache-freshness
branch.
