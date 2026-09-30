# Analysis module rules

## The real tree

Verified against `src/polyzymd/analyses/` on 2026-09-29, after slice 7 of the
v1.3 analysis refactor. Every file named here exists. If you add or delete a
module, update this list in the same commit.

```
src/polyzymd/analyses/
├── study.py             # Study, Condition, Replicate: replicates as MDAnalysis universes
├── timeseries.py        # Study.timeseries, Study.per_replicate, Timeseries, ReplicateValues
├── functions.py         # Shipped measurements: radius_of_gyration, rmsd, pair_distance,
│                        # all_below, rmsf, rms_deviation, rms_decomposition, sasa,
│                        # residue_sasa, dssp_occupancy
├── reference.py         # pz.reference: external, frame, average and centroid references
├── figures.py           # Figures drawn from stored study results
├── protocols.py         # polyzymd analyze: FUNCTION_ANALYSES, ProtocolReport, plugin path
├── base.py              # Public import surface for plugin authors
├── discovery.py         # pkgutil auto-discovery of plugins
├── orchestrator.py      # Plugin engine: compute, aggregate, compare, plot
├── stats.py             # default_scalar_comparison, format_scalar_comparison
├── exceptions.py        # Typed analysis errors
├── _framework/          # aggregate_validation, cache_identity, compare,
│                        # comparison_models, contexts, contract, io,
│                        # lifecycle, results_base
├── mda/                 # aggregation, artifacts, base, comparison,
│                        # frame_selection, job, lifecycle, plugin, store, universe
├── shared/              # aa_classification, autocorrelation, centroid, diagnostics,
│                        # inferential_statistics, loader, multi_run_formatting,
│                        # paths, plotting, selections, statistics, topology,
│                        # window, groupings/, selectors/
├── contacts/            # __init__, _aggregator, _comparison,
│                        # _comparison_results, _events, _filters, _formatters,
│                        # _identity, _lifecycle, _mda, _plot_settings, _plotters
└── hydrogen_bonds/      # __init__, _mda, _models, _plotters
```

rg, rmsd, rmsf, rms_deviation, distances, catalytic_triad, sasa and
secondary_structure are no longer plugins. They are functions in `functions.py`
run through the study API, and `polyzymd analyze <name>` runs them through
`protocols._analyze_function`. `polyzymd analyze contacts` also runs through
the study API (`functions.residue_contacts`, `protocols._analyze_contacts`) for
contact fractions and coverage; the contacts plugin stays only for residence
times through `polyzymd compare run contacts`, until those are ported. Only
contacts and hydrogen_bonds are still plugins; the refactor ports them next and
then removes the plugin framework.

There is no `_results.py`, `_cache.py`, `_paths.py` or `_plotting.py` in any
plugin package. Result models live in `_models.py` (hydrogen bonds), in
`_comparison_results.py`, or in the plugin `__init__.py`. Cache and path
handling belong to the framework artifact layer, not to plugins.

## Adding a measurement

Write a new measurement as a function, not a plugin, and run it through the
study API; see `docs/source/explanation/analysis_api.md`.

- A per-frame function takes MDAnalysis `AtomGroup`s at one frame and returns
  one number. Run it with `study.timeseries(fn, pz.select(...), ...)`.
- A per-replicate function also takes `frames` and returns a number or one value
  per label, such as per residue. Run it with `study.per_replicate(fn, ...,
  labels=..., parts=...)`.
- Call MDAnalysis or MDTraj for the measurement; do not reimplement them.
- To ship it in `polyzymd analyze`, add it to `FUNCTION_ANALYSES` with its
  settings and a dispatcher in `protocols.py`, as the existing analyses do, with
  its figures, a quick-start page and a real-data parity check.

## Loading trajectories

`polyzymd.analyses.shared.loader.TrajectoryLoader` is the canonical universe
loader. It resolves topology and trajectory files for a replicate, checks
segment lineage, and builds the MDAnalysis universe.

`polyzymd.analyses.mda.universe.UniverseProvider` wraps `TrajectoryLoader`. It
takes a `SimulationConfig`, instantiates the loader lazily, and adds input
provenance (`UniverseProvenance`) to each load. Plugins and framework code use
`UniverseProvider`. Nothing else should build a `Universe` directly, and no
plugin should construct file paths by hand.

## Public import surface

Import from `polyzymd.analyses.base`. It re-exports `Analysis`, the four
lifecycle contexts, `MetricValue`, the comparison models, `PluginContractError`
and `SlurmResourceHint`. Do not import `_framework/` modules from a plugin.

## Adding a plugin (contacts and hydrogen_bonds only)

1. Run `polyzymd new-analysis <name>` to scaffold the package and its tests, or
   write the package by hand under `src/polyzymd/analyses/`.
2. Define a `Settings` class as a Pydantic v2 `BaseModel`.
3. Subclass `Analysis` and set `name` and `Settings` as `ClassVar`s.
4. With `has_compute_stage=True`, implement `build_mda_jobs()` and, when the
   plugin needs one, `build_mda_collector()`. Put `AnalysisBase` subclasses in
   `_mda.py`.
5. Set `has_compute_stage=False` for a compare-only plugin. Setting
   `has_aggregate_stage=True` with `has_compute_stage=False` raises
   `PluginContractError`.
6. Implement `aggregate()` only when `has_aggregate_stage=True`.
7. Discovery is automatic through `pkgutil`. There is no registry to edit.

## Lifecycle hooks and contexts

| Hook | When | Context | Returns |
|------|------|---------|---------|
| `build_mda_jobs()` plus `build_mda_collector()` | `has_compute_stage=True` | `ReplicateContext` | `ReplicateArtifact` through the collector |
| `aggregate()` | `has_aggregate_stage=True` | `AggregateContext` | Pydantic model or dict |
| `compare()` | Once per analysis | `ComparisonContext` | Pydantic model, or `None` |
| `plot()` | Once per analysis | `PlotContext` | `list[Path]` |
| `extract_metrics()` | Default compare path | `ComparisonContext` | `dict[str, MetricValue]` |
| `filter_conditions()` | Optional | conditions | filtered conditions |
| `format()` | Optional | comparison result | CLI text |

Contexts carry what a plugin needs. Never load a config inside a plugin.
`PlotContext.plot_settings` is always a valid `PlotSettings`, so do not guard
against `None`. A hook that returns a type outside the contract raises
`PluginContractError`.

## Two comparison paths

The simple path implements `extract_metrics()` and lets `stats.py` run the
t-tests, the ANOVA, the Benjamini-Hochberg correction and the ranking. The
custom path overrides `compare()` and returns its own saveable model. Both
remaining plugins, contacts and hydrogen_bonds, take the custom path.

## Results and plotting

Persist replicate, condition and comparison artifacts through `ArtifactStore`.
Large arrays and event tables go in validated sidecars that the artifact
payload refers to. Do not invent a plugin-specific cache filename scheme.

`plot()` reads cached artifacts and sidecars. It must not reload a trajectory
or rerun an analysis.

## Statistical contract

The replicate is the sampling unit for every cross-condition test and every
condition-level uncertainty. Equilibration is one global value applied
uniformly, and no diagnostic is allowed to select data. Every metric carries a
unit and a stated uncertainty. Invoke the `livecoms-check` skill in
`.claude/skills/` before committing analysis code.
