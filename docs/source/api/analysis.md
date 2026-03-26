# Analysis API

This page is the entry point for the analysis-side API.

PolyzyMD analysis code is organized into two packages:

| Package | Role |
|---------|------|
| `polyzymd.analysis` | Per-condition **calculators** (RMSF, distances, contacts, etc.) |
| `polyzymd.analyses` | **Plugin system** — unified analysis lifecycle (see [Analyses Plugin API](analyses.md)) |

## Core Analysis Interfaces

- `polyzymd.analysis`
- `polyzymd.analysis.config`
- `polyzymd.analysis.core.registry`
- `polyzymd.analysis.results.base`

These modules define the reusable contracts for analysis settings, analyzers,
and serialized result objects.

## Stable Analysis Implementations

- `polyzymd.analysis.rmsf`
- `polyzymd.analysis.distances`
- `polyzymd.analysis.contacts`
- `polyzymd.analysis.triad`
- `polyzymd.analysis.secondary_structure`

## Plugin System

For the unified analysis lifecycle (compute → aggregate → compare → plot →
format), see the [Analyses Plugin API](analyses.md).

## Related Documentation

- [Analyses Plugin API](analyses.md)
- [Configuration API](config.md)
- [Compare API](compare.md)
- [API Overview](overview.rst)
- [Extending the Analysis Framework](../tutorials/extending_analyses.md)

## Notes

- The analysis API is still evolving rapidly as the post-simulation analysis
  stack stabilizes.
- Stable release workflows are RMSF, contacts, distances, catalytic triad, and
  secondary structure.
- Experimental workflows remain documented separately in the user docs and are
  clearly labeled there.

<!-- IMAGE OPPORTUNITY: Add a compact analysis architecture diagram showing
analysis core -> analyses plugin -> comparison -> plot generation. -->
