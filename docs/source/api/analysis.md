# Analysis API

This page is the entry point for the analysis-side API.

PolyzyMD analysis code is organized around a small shared core plus
analysis-specific implementations for stable and experimental workflows.

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

## Related Documentation

- [Configuration API](config.md)
- [Compare API](compare.md)
- [API Overview](overview.rst)
- [Analyze a Multi-Condition Study](../tutorials/analysis_complete_workflow.md)

## Notes

- The analysis API is still evolving rapidly as the post-simulation analysis
  stack stabilizes.
- Stable release workflows are RMSF, contacts, distances, catalytic triad, and
  secondary structure.
- Experimental workflows remain documented separately in the user docs and are
  clearly labeled there.

<!-- IMAGE OPPORTUNITY: Add a compact analysis architecture diagram showing
analysis core -> per-condition analysis -> comparison -> plot generation. -->
