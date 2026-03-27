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

- `polyzymd.compare.results.rmsf`
- `polyzymd.compare.results.triad`
- `polyzymd.compare.results.contacts`
- `polyzymd.compare.results.distances`
- `polyzymd.compare.results.secondary_structure`
- `polyzymd.compare.results.binding_free_energy`
- `polyzymd.compare.results.exposure`
- `polyzymd.compare.results.polymer_affinity`

### Formatting Helpers

- `polyzymd.compare.formatters`
- `polyzymd.compare.contacts_formatters`
- `polyzymd.compare.distances_formatters`
- `polyzymd.compare.exposure_formatters`
- `polyzymd.compare.triad_formatters`
- `polyzymd.compare.binding_free_energy_formatters`
- `polyzymd.compare.polymer_affinity_formatters`

### Plotting Helpers

- `polyzymd.compare.plotter` - shared plotter base classes and registry
- `polyzymd.compare.plotters.rmsf`
- `polyzymd.compare.plotters.triad`
- `polyzymd.compare.plotters.contacts`
- `polyzymd.compare.plotters.contacts_binding_preference`
- `polyzymd.compare.plotters.contacts_coverage`
- `polyzymd.compare.plotters.contacts_grouped_bars`
- `polyzymd.compare.plotters.contacts_profiles`
- `polyzymd.compare.plotters.distances`
- `polyzymd.compare.plotters.secondary_structure`
- `polyzymd.compare.plotters.exposure`
- `polyzymd.compare.plotters.binding_free_energy`
- `polyzymd.compare.plotters.polymer_affinity`

## How This Fits Together

- A plugin in `polyzymd.analyses` owns the user-facing workflow for one analysis.
- The orchestrator resolves `plugins:` settings from `comparison.yaml`.
- Shared compare-layer models and utilities live in `polyzymd.compare`.
- Canonical comparison caches are written under `comparison/<analysis>/result.json`.
- Plotting helpers rediscover cached comparison data from disk.

## Related Pages

- [API Overview](overview.rst)
- [Analyses Plugin System API](analyses.md)
- [Analysis Calculator API](analysis.md)
