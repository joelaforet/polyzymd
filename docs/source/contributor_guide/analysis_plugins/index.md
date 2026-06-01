# New Analysis Contributor Path

This landing page is for contributors who want to add a new PolyzyMD analysis
plugin or make a safe change to an existing plugin. It is a map of the
contributor path, not the tutorial itself.

Use it if you know Python and molecular simulation concepts, but are still
learning how PolyzyMD connects MDAnalysis trajectory work to replicate
artifacts, condition aggregation, cross-condition comparison, and plotting.

## Start

Begin here when you are new to contributing analysis code.

- [Set Up a Contributor Environment](../setup.md) — install the pixi-managed
  environment used by docs, tests, and local development.
- [Contributing to PolyzyMD](../contributing.md) — learn the branch, review, and
  verification expectations for small, reviewable changes.
- [Architecture](../../explanation/architecture.md) — orient to the package
  layout before narrowing in on analysis plugins.
- [How PolyzyMD analysis plugins work](architecture.md) — understand the plugin
  lifecycle, public facades, artifacts, sidecars, and artifact-only plotting.

## Build

Use these pages when you are creating your first plugin or replacing scaffolded
placeholder logic with real analysis code.

- [Scaffold your first analysis plugin](first_scaffold.md) — generate a
  single-file plugin and matching tests with a throwaway name such as
  `solvent_shell`.
- [Build a simple scalar analysis plugin](simple_scalar_plugin.md) — connect one
  function job, collector, default aggregation path, and scalar comparison.
- [Extend PolyzyMD with MDAnalysis-native analyses](../extending_analyses.md) —
  use the full implementation guide for current scaffold patterns, public
  imports, `MDAAnalysisJob`, collectors, artifacts, default comparison, and
  artifact-only plotting.

## Scale

Use these pages when a working plugin needs richer stored data or a larger module
shape.

- [Store large analysis outputs with artifact sidecars](sidecars.md) — persist
  NPZ arrays, CSV tables, and other bulky outputs beside artifacts without
  bloating JSON.
- [Convert a plugin into an advanced package](advanced_package.md) — split a
  single-file plugin into `_mda.py`, `_plotters.py`, `_models.py`, or
  `_formatters.py` helpers when reviewability improves.

## Verify

Use these pages before opening a pull request or reviewing one.

- [Test your analysis plugin contribution](testing.md) — cover discovery,
  settings validation, MDAnalysis jobs, collectors, artifacts, sidecars,
  aggregation, comparison metrics, and artifact-only plotting.
- [Analysis plugin contribution checklist](checklist.md) — confirm the plugin
  contract, import boundaries, tests, docs updates, and PR-ready commands.

## Reference

Use these lookup pages when you need exact import paths, commands, or settings.

- [API map for analysis plugin contributors](api_map.md) — find stable public
  imports and the nearest API reference page without browsing internals.
- [Analysis plugin API](../../api/analyses.md),
  [analysis base classes](../../api/analyses_base.md),
  [MDAnalysis integration](../../api/analyses_mda.md), and
  [shared analysis utilities](../../api/analyses_shared.md) — inspect the API
  reference directly.
- [CLI reference](../../reference/cli_reference.md) and
  [analysis plugin settings reference](../../reference/analysis_plugin_settings.md)
  — look up commands and `comparison.yaml` plugin options.

## Public and private import guardrails

Contributor-facing plugin code should use public facades:

- `polyzymd.analyses.base` for `Analysis`, lifecycle contexts,
  `MetricValue`, and comparison result models.
- `polyzymd.analyses.mda` for MDAnalysis jobs, frame selection, artifacts,
  artifact stores, aggregation helpers, and comparison helpers.
- documented utilities from `polyzymd.analyses.shared` when an existing shared
  helper fits the task.

Do not import from `polyzymd.analyses._framework` in contributor plugins. That
package is an internal implementation detail behind the public facades.

```{toctree}
:hidden:
:maxdepth: 1

How PolyzyMD analysis plugins work <architecture>
Scaffold your first analysis plugin <first_scaffold>
Build a simple scalar analysis plugin <simple_scalar_plugin>
Store large analysis outputs with artifact sidecars <sidecars>
Convert a plugin into an advanced package <advanced_package>
Test your analysis plugin contribution <testing>
API map for analysis plugin contributors <api_map>
Analysis plugin contribution checklist <checklist>
```
