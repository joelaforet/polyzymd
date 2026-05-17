# Analyses Plugin System API

This page documents the `polyzymd.analyses` package — the plugin system for
adding new analysis types to PolyzyMD.

PolyzyMD analyses are MDAnalysis analyses lifted from one trajectory to
condition/replicate ensembles. New contributor analyses usually use the
MDAnalysis job/artifact API: a plugin subclasses `Analysis`, builds
`MDAAnalysisJob` objects, maps completed MDAnalysis results to canonical
artifacts, and lets the framework handle cache identity, aggregation, scalar
comparison, CLI wiring, and plotting orchestration.

The public MDAnalysis extension-layer reference is on {doc}`analyses_mda`.
The old scalar Measurement API was removed in v1.3 after the MDAnalysis
job/artifact lifecycle became the only supported contributor extension path.

## Public API

```{eval-rst}
.. currentmodule:: polyzymd.analyses
```

### Discovery

- `list_analyses()` — return dict of `{name: class}` for all discovered plugins
- `list_all_names()` — return list of all names including aliases
- `get_analysis(name)` — get a plugin class by name or alias
- `clear_cache()` — reset the discovery cache

### Orchestration

- `run_analysis(analysis: Analysis, condition: Condition, settings: Any, equilibration: str = "0ns", output_dir: Path | None = None, recompute: bool = False)` — run compute + aggregate for one condition
- `run_comparison(analysis: Analysis, config: ComparisonConfig, recompute: bool = False, equilibration: str | None = None)` — run full lifecycle (compute + aggregate + compare + plot)
- `run_all_comparisons(config: ComparisonConfig, analysis_names: list[str] | None = None, recompute: bool = False, equilibration: str | None = None)` — run multiple analyses from one comparison config

`ComparisonConfig` is defined in `polyzymd.config.comparison`.

### Base Class

- `Analysis` — abstract base class all plugins inherit from
- `build_mda_jobs(ctx)` — hook implemented by compute-stage plugins to construct `MDAAnalysisJob` objects for one replicate
- `build_mda_collector(ctx)` — hook that maps completed jobs to a `ReplicateArtifact`
- `extract_metrics(summary)` — optional hook for default scalar comparison of condition summaries

`polyzymd.analyses.base` is the stable public facade for contributor imports.
It re-exports the base class, lifecycle contexts, metric descriptors, and
comparison result models while delegating implementation to private framework
modules. Contributor plugins should import from
`polyzymd.analyses.base`, not from modules named `_analysis_*`, `_contexts`,
or `_comparison_models`.

### Discovery

Discovery supports both contributor-friendly single-file plugins such as
`polyzymd.analyses.my_metric` and package plugins such as
`polyzymd.analyses.contacts`. Use a package when the plugin needs private helper
modules for MDAnalysis jobs, results, plotting, or formatting; use a single file
for the default function-adapter MDA job pattern.

### MDAnalysis extension layer

The `polyzymd.analyses.mda` namespace provides the public job/artifact layer:

- `MDAAnalysisJob` wraps one `AnalysisBase`-compatible object or simple function adapter
- `FrameSelection` maps PolyzyMD windows to MDAnalysis `run()` kwargs
- `MDAFunctionAdapter` powers the default scaffold function path
- `ReplicateArtifact`, `ConditionArtifact`, and `ComparisonArtifact` define cache envelopes
- `ArtifactStore` persists `result.json`, manifests, and sidecars
- `aggregate_replicate_artifacts()` provides default aggregation from `payload["metrics"]`
- `compare_condition_artifacts()` compares condition artifacts with replicate-level statistics

See {doc}`analyses_mda` for autodoc reference.

### Context Objects

- `ReplicateContext` — passed to per-replicate compute hooks
- `AggregateContext` — passed to `aggregate()`
- `ComparisonContext` — passed to `compare()`
- `PlotContext` — passed to `plot()`
- `Condition` — represents one simulation condition

### Result Models

- `ComparisonResult` — universal comparison result (Pydantic model with `.save()/.load()`)
- `ConditionSummary` — per-condition statistics
- `PairwiseResult` — pairwise t-test result
- `ANOVAResult` — ANOVA result
- `MetricValue` — scalar metric descriptor for default comparison
- `BaseComparisonResult` — abstract base for custom plugin comparison results
- `BaseConditionSummary` — abstract base for per-condition summaries

## Available Plugins

| Plugin | Module | Comparison Style |
|--------|--------|-----------------|
| `rmsd` | `analyses.rmsd` | MDAnalysis job artifacts with custom per-run comparison |
| `rg` | `analyses.rg` | Custom `AnalysisBase` artifacts with custom per-run comparison |
| `rmsf` | `analyses.rmsf` | Profile artifacts with default scalar comparison |
| `catalytic_triad` | `analyses.catalytic_triad` | Pair-distance artifact reducer with default scalar comparison |
| `secondary_structure` | `analyses.secondary_structure` | Categorical matrix artifacts with default scalar comparison |
| `sasa` | `analyses.sasa` | AnalysisBase-compatible SASA wrapper with custom per-run comparison |
| `distances` | `analyses.distances` | Pair-distance artifacts with custom comparison |
| `contacts` | `analyses.contacts` | Sparse contact-event artifacts with custom comparison |
| `hydrogen_bonds` | `analyses.hydrogen_bonds` | MDAnalysis hydrogen-bond event artifacts with custom summaries |

### Contacts package facade

`polyzymd.analyses.contacts` exposes the public `ContactsAnalysis` plugin and
its supported settings/result classes. The package is intentionally split into
private helper modules for contact-event detection, artifact aggregation,
condition filtering, custom comparison, plotting, result models, MDAnalysis job
helpers, and sidecar loading. These `contacts/_*.py` modules are internal
implementation details and are not separate contributor API entry points.

## Shared Utilities

The `analyses/shared/` package provides reusable infrastructure used across
plugins.

### Convergence Diagnostics

The `convergence` module provides a sliding-window slope heuristic for
detecting sustained convergence in timeseries data (e.g., RMSD traces).

```{eval-rst}
.. automodule:: polyzymd.analyses.shared.convergence
   :members:
   :undoc-members:
   :show-inheritance:
   :no-index:
```

### Multi-Run Comparison Orchestration

The `multi_run_comparison` module provides helpers for plugins that compare
multiple named runs (e.g., per-chain RMSD, per-domain Rg) across conditions.
Used by the RMSD, Rg, and SASA plugins.

```{eval-rst}
.. automodule:: polyzymd.analyses.shared.multi_run_comparison
   :members:
   :undoc-members:
   :show-inheritance:
```

### Multi-Run Formatting

The `multi_run_formatting` module provides text and markdown formatting helpers
for multi-run analysis CLI output — ranked tables, pairwise lines, and ANOVA
summaries.

```{eval-rst}
.. automodule:: polyzymd.analyses.shared.multi_run_formatting
   :members:
   :undoc-members:
   :show-inheritance:
```

## Related Documentation

- [Extending PolyzyMD with MDAnalysis-native analyses](../contributor_guide/extending_analyses.md) — contributor guide
- [MDAnalysis Extension-Layer API](analyses_mda.md) — API reference
- [API Overview](overview.rst)
