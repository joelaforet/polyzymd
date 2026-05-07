# Analyses Plugin System API

This page documents the `polyzymd.analyses` package — the plugin system for
adding new analysis types to PolyzyMD.

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

`polyzymd.analyses.base` is the stable public facade for contributor imports.
It re-exports the base class, lifecycle contexts, metric descriptors, and
comparison result models while delegating implementation to private framework
modules. Contributor plugins should import from `polyzymd.analyses.base`, not
from modules named `_analysis_*`, `_contexts`, or `_comparison_models`.

### Context Objects

- `ReplicateContext` — passed to runner-backed replicate hooks
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
| `rmsd` | `analyses.rmsd` | Custom (per-run) |
| `rg` | `analyses.rg` | Custom (per-run) |
| `rmsf` | `analyses.rmsf` | Default (scalar) |
| `catalytic_triad` | `analyses.catalytic_triad` | Default (scalar) |
| `secondary_structure` | `analyses.secondary_structure` | Default (scalar) |
| `sasa` | `analyses.sasa` | Custom (per-run) |
| `distances` | `analyses.distances` | Custom |
| `contacts` | `analyses.contacts` | Custom |
| `hydrogen_bonds` | `analyses.hydrogen_bonds` | Default (scalar) |
| `binding_free_energy` | `analyses.binding_free_energy` | Custom (experimental) |
| `polymer_affinity` | `analyses.polymer_affinity` | Custom (experimental) |
| `polymer_bridging` | `analyses.polymer_bridging` | Custom (experimental) |

### Contacts package facade

`polyzymd.analyses.contacts` exposes the public `ContactsAnalysis` plugin and
its supported settings/result classes. The package is intentionally split into
private helper modules for cache handling, lifecycle dispatch, condition
filtering, custom comparison, binding-preference support, plotting, result
models, and trajectory runners. These `contacts/_*.py` modules are internal
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

- [Extending the Analysis Framework](../contributor_guide/extending_analyses.md) — tutorial
- [Analysis API Overview](analysis.md) — entry point for the analysis API
- [API Overview](overview.rst)
