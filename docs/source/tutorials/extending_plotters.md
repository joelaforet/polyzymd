# Plotting in an Analysis Plugin

PolyzyMD uses one extension workflow for new contributions: add a single
`Analysis` subclass in `src/polyzymd/analyses/`.

If your contribution needs figures, put that logic in the plugin's `plot()`
method.

Use {doc}`extending_analyses` as the main contributor guide.

## What to Implement

- `plot(ctx)` to generate figures for your analysis
- optional helper functions — these can live in the same file for simple
  plugins, or be extracted to a `_plotters.py` module within the package
- optional `format()` for matching CLI output and figures

## Recommended Pattern

- load any cached comparison data from the plugin context
- write figures into the analysis-specific output directory provided by `ctx`
- keep compute, compare, and plot behavior close together so contributors only
  need to reason about one package

## When to Extract `_plotters.py`

As a plugin grows, extracting plotting functions into a dedicated
`_plotters.py` module within the package keeps `__init__.py` focused on the
Analysis lifecycle. This is the established pattern for rmsf, distances,
secondary_structure, and contacts. Consider extracting when your plugin has
3+ plot functions or `__init__.py` exceeds ~500 lines.

See {ref}`plotters-extraction` in the extending analyses tutorial for details.

## Next Step

- Read {doc}`extending_analyses`
