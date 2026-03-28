# Compare Module

This page maps the current `polyzymd.compare` package.

The comparison framework is now driven by plugins in `polyzymd.analyses`.
The `polyzymd.compare` package provides shared building blocks that those
plugins can reuse for configuration, statistics, result models, formatting, and
plotting.

## Where to Start

- New plugin author: [Extending the Analysis Framework](../tutorials/extending_analyses.md)
- Plugin system API: [Analyses Plugin System API](analyses.md)
- User workflow: [How to Compare Simulation Conditions](../tutorials/analysis_compare_conditions.md)

## Package Map

### Configuration and CLI

- `polyzymd.compare.config` - `ComparisonConfig`, condition models, plot settings
- `polyzymd.compare.settings` - shared settings models used by comparison features
- `polyzymd.compare.cli` - `polyzymd compare` subcommands
- `polyzymd.compare.cli_utils` - config loading and validation helpers

### Statistics and IO

- `polyzymd.compare.statistics` - t-tests, ANOVA, effect sizes, ranking helpers
- `polyzymd.compare.io.results` - canonical comparison result discovery and loading
- `polyzymd.compare.io.paths` - comparison and figure path helpers

### Result Models

- `polyzymd.analyses.rmsf._comparison_results`
- `polyzymd.analyses.catalytic_triad._comparison_results`
- `polyzymd.analyses.contacts._comparison_results`
- `polyzymd.analyses.distances._comparison_results`
- `polyzymd.analyses.secondary_structure._comparison_results`
- `polyzymd.analyses.binding_free_energy._comparison_results`
- `polyzymd.analyses.exposure._comparison_results`
- `polyzymd.analyses.polymer_affinity._comparison_results`

### Formatting Helpers

- `polyzymd.analyses.contacts._formatters`
- `polyzymd.analyses.distances._formatters`
- `polyzymd.analyses.exposure._formatters`
- `polyzymd.analyses.binding_free_energy._formatters`
- `polyzymd.analyses.polymer_affinity._formatters`

Catalytic triad and RMSF formatting is inline in the plugin's `format()` method;
secondary structure formatting is likewise inline.

### Plotting Helpers

- `polyzymd.analyses.shared.plotting` - shared figure utilities used by analysis plugins

Plotting logic lives inside each analysis plugin's `plot()` method. Shared
helpers for axis styling, legends, color palettes, and figure saving are in
`polyzymd.analyses.shared.plotting`.

## How This Fits Together

- A plugin in `polyzymd.analyses` owns the user-facing workflow for one analysis.
- The orchestrator resolves `plugins:` settings from `comparison.yaml`.
- Shared compare-layer models and utilities live in `polyzymd.compare`.
- Canonical comparison caches are written under `comparison/<analysis>/result.json`.
- Plotting helpers rediscover cached comparison data from disk.
- Each plugin's `plot()` method produces publication-quality figures.

## Related Pages

- [API Overview](overview.rst)
- [Analyses Plugin System API](analyses.md)
- [Analysis Calculator API](analysis.md)
