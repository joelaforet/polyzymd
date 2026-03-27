# Plotting in an Analysis Plugin

PolyzyMD uses one extension workflow for new contributions: add a single
`Analysis` subclass in `src/polyzymd/analyses/`.

If your contribution needs figures, put that logic in the plugin's `plot()`
method.

Use {doc}`extending_analyses` as the main contributor guide.

## What to Implement

- `plot(ctx)` to generate figures for your analysis
- optional helper functions in the same file when that keeps the workflow clear
- optional `format()` for matching CLI output and figures

## Recommended Pattern

- load any cached comparison data from the plugin context
- write figures into the analysis-specific output directory provided by `ctx`
- keep compute, compare, and plot behavior close together so contributors only
  need to reason about one file

## Next Step

- Read {doc}`extending_analyses`
