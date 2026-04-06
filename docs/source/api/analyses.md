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

- `run_analysis(name, config_path)` — run a single analysis (compute + aggregate)
- `run_comparison(name, config_path)` — run full lifecycle (compute + aggregate + compare + plot)
- `run_all_comparisons(config_path)` — run all configured analyses

### Base Class

- `Analysis` — abstract base class all plugins inherit from

### Context Objects

- `ReplicateContext` — passed to `compute_replicate()`
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
| `exposure` | `analyses.exposure` | Custom |
| `binding_free_energy` | `analyses.binding_free_energy` | Custom |
| `polymer_affinity` | `analyses.polymer_affinity` | Custom |
| `polymer_bridging` | `analyses.polymer_bridging` | Custom (experimental) |

## Related Documentation

- [Extending the Analysis Framework](../tutorials/extending_analyses.md) — tutorial
- [Analysis Calculator API](analysis.md) — underlying compute layer
- [Compare API](compare.md) — statistics and formatting
- [API Overview](overview.rst)
