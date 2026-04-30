# Comparison Logic in an Analysis Plugin

PolyzyMD uses one extension workflow for new contributions: add a single
`Analysis` subclass in `src/polyzymd/analyses/`.

If your contribution needs custom cross-condition statistics, put that logic in
the plugin's `compare()` method.

Use {doc}`extending_analyses` as the main contributor guide.

## What to Implement

- `build_runner()` and `summarize_replicate()` for per-replicate work
- `aggregate()` for per-condition summaries
- `compare()` when the default scalar comparison path is not enough
- `format()` if you want custom CLI output

## Keep It Simple

For a new community contribution, the preferred shape is:

1. create one package in `src/polyzymd/analyses/<name>/`
2. keep settings in the plugin's inner `Settings` model
3. use `plugins:` in `comparison.yaml`
4. let the orchestrator manage cache files under `comparison/<analysis>/result.json`

That keeps the full user-facing workflow in one place and makes reviews easier.

## Next Step

- Read {doc}`extending_analyses`
