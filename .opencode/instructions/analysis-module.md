# Analysis module rules

## The real tree

Verified against `src/polyzymd/analyses/` on 2026-09-30, after the
hydrogen-bonds slice of the v1.3 analysis refactor. Every file named here exists. If you add or delete a
module, update this list in the same commit.

```
src/polyzymd/analyses/
├── study.py             # Study, Condition, Replicate: replicates as MDAnalysis universes
├── timeseries.py        # Study.timeseries, Study.per_replicate, Timeseries, ReplicateValues
├── functions.py         # Shipped measurements: radius_of_gyration, rmsd, pair_distance,
│                        # all_below, native_contacts, rmsf, rmsd_per_residue,
│                        # rms_decomposition, sasa, residue_sasa, dssp_occupancy,
│                        # residue_contacts, residue_occlusion, contact_events,
│                        # restricted_mean_lifetime, contact_lifetimes, event_lifetimes,
│                        # hbond_atoms, hydrogen_bonds, hbond_lifetimes,
│                        # residue_hbond_occupancy, residue_pair_hbond_occupancy, hbond_count
├── reference.py         # pz.reference: external, frame, average and centroid references
├── figures.py           # Figures drawn from stored study results
├── protocols.py         # polyzymd analyze: FUNCTION_ANALYSES, ProtocolReport, plugin path
├── base.py              # Plugin framework: public import surface for plugin authors
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
```

There is no plugin package under `analyses/` any more.

rg, rmsd, rmsf, rmsd_per_residue, distances, sasa, secondary_structure,
native_contacts, contacts and hydrogen_bonds are not plugins. They are functions in
`functions.py` run through the study API, listed with their settings in
`protocols.FUNCTION_ANALYSES`, and `polyzymd analyze <name>` runs them through
`protocols._analyze_function`. `polyzymd analyze contacts` is dispatched by
`protocols._analyze_contacts`: `functions.residue_occlusion` for
`method=occlusion`, the default, and `functions.residue_contacts` for
`method=distance` give coverage, contact fractions and occluded area, and
`functions.contact_lifetimes` gives `--run mean_lifetime`, `lifetime_events`
and `censored_fraction`. `polyzymd analyze hydrogen_bonds` is dispatched by
`protocols._analyze_hydrogen_bonds`: `functions.hydrogen_bonds` gives each
summary's counts, `functions.hbond_lifetimes` its lifetimes,
`functions.residue_hbond_occupancy` its `_residues` and
`functions.residue_pair_hbond_occupancy` (with `labels="returned"`) its `_pairs`.

Retired names: a `plugins.<name>` or `plot_settings.<name>` block for any of
them in `comparison.yaml` is ignored with a warning
(`config.comparison.RETIRED_PLUGINS`, `RETIRED_PLUGIN_WARNING`) that names the
command, the Python function, the analyze protocol page and the agent skill,
and `polyzymd compare run <name>` exits with the `polyzymd analyze` command
(`cli._compare_utils.analyze_command_for`). `catalytic_triad` is a routine on
the study API (`docs/source/how_to/analysis_triad_quickstart.md`), not an
analysis: `polyzymd analyze catalytic_triad` raises `ProtocolError`
(`protocols._refuse_retired`) pointing to that page and to
`polyzymd analyze distances --set pairs=...`. `polyzymd analyze -f
comparison.yaml` raises `ProtocolError` with the equivalent `-c` command
(`cli.analyze._refuse_comparison_file`). Warnings for retired things name the
replacement, the docs page to read and what to point an agent at.

No shipped analysis is a plugin any more. The plugin framework (`base.py`,
`discovery.py`, `orchestrator.py`, `stats.py`, `mda/`, `_framework/`, the
`polyzymd compare` commands, `polyzymd new-analysis` and
`workflow/analysis_slurm.py`) stays until a later slice removes it; its tests
register the toy plugin of `tests/analyses/test_protocols.py` (`_install_toy`).
Do not add plugins.

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
code should construct file paths by hand.

## Plugin framework (being removed)

The rest of this file describes the plugin framework for the code that still
uses it. Import from `polyzymd.analyses.base`; it re-exports `Analysis`, the
four lifecycle contexts, `MetricValue`, the comparison models,
`PluginContractError` and `SlurmResourceHint`. Discovery is automatic through
`pkgutil`.

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
custom path overrides `compare()` and returns its own saveable model.

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
