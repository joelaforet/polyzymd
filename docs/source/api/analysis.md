# Analysis API

This page is the entry point for the analysis-side API.

All analysis code lives in `polyzymd.analyses`:

| Package | Role |
|---------|------|
| `polyzymd.analyses` | **Plugin system** — unified analysis lifecycle (compute → aggregate → compare → plot → format) |
| `polyzymd.analyses.shared` | Reusable utilities (trajectory loading, alignment, autocorrelation, statistics) |
| `polyzymd.analyses.<name>._results` | Pydantic result models for serialization (existing plugins) |

## Plugin System

For the unified analysis lifecycle, see the [Analyses Plugin API](analyses.md).

## Related Documentation

- [Analyses Plugin API](analyses.md)
- [Configuration API](config.md)
- [API Overview](overview.rst)
- [Extending the Analysis Framework](../contributor_guide/extending_analyses.md)

## Notes

- **Stable plugins** (9): RMSD, Rg, RMSF, catalytic triad, secondary structure,
  SASA, distances, contacts, and hydrogen bonds.

<!-- IMAGE OPPORTUNITY: Add a compact analysis architecture diagram showing
analyses plugin -> comparison -> plot generation. -->
