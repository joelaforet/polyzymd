# New Analysis Contributor Path

This landing page is for contributors who want to add a new PolyzyMD analysis
plugin or make a safe change to an existing plugin. It is a map of the
contributor path, not the tutorial itself.

Use it if you know Python and molecular simulation concepts, but are still
learning how PolyzyMD connects MDAnalysis trajectory work to replicate
artifacts, condition aggregation, cross-condition comparison, and plotting.

## Recommended order

1. **Set up the contributor environment.** Start with
   [Set Up a Contributor Environment](../setup.md) so commands run inside the
   managed `pixi` environment.
2. **Read the project contribution basics.** Use
   [Contributing to PolyzyMD](../contributing.md) for project workflow and review
   expectations.
3. **Orient to the package architecture.** Read
   [Architecture](../../explanation/architecture.md) for the high-level package
   layout before focusing on analyses.
4. **Understand how analysis plugins work.** Read
   [How PolyzyMD analysis plugins work](architecture.md) to learn why plugins
   use `MDAAnalysisJob`, collectors, artifacts, sidecars, and artifact-only
   plotting before you start editing code.
5. **Scaffold your first plugin.** Follow
   [Scaffold your first analysis plugin](first_scaffold.md) to generate one
   default single-file plugin and matching tests with the throwaway name
   `solvent_shell`.
6. **Build one small scalar plugin.** Follow
   [Build a simple scalar analysis plugin](simple_scalar_plugin.md) to learn how
   a function job, collector, default aggregation, and `MetricValue` fit together.
7. **Store arrays and tables as sidecars.** Use
   [Store arrays and tables with sidecars](sidecars.md) when your collector needs
   to persist NPZ arrays, CSV tables, or other bulky data beside an artifact.
8. **Split larger plugins into packages when needed.** Use
   [Convert a plugin into an advanced package](advanced_package.md) when a
   working single-file plugin has grown enough to need `_mda.py`, `_plotters.py`,
   `_results.py`, or `_formatters.py` helpers.
9. **Use the current full analysis extension guide.** The supported detailed
   guide is [Extend PolyzyMD with MDAnalysis-native analyses](../extending_analyses.md).
   It covers the current scaffold, public imports, `MDAAnalysisJob`, collectors,
   artifacts, default comparison, and artifact-only plotting.
10. **Look up stable API surfaces when needed.** Use the API pages for
    [analysis plugins](../../api/analyses.md),
    [analysis base classes](../../api/analyses_base.md),
    [MDAnalysis integration](../../api/analyses_mda.md), and
    [shared analysis utilities](../../api/analyses_shared.md).
11. **Check CLI and configuration reference material.** Use the
     [CLI reference](../../reference/cli_reference.md) and
     [analysis plugin settings reference](../../reference/analysis_plugin_settings.md)
     when you need command or configuration lookup.

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

## Remaining planned path pages

The New Analysis Contributor Path will continue to grow into smaller pages with
one purpose each. These pages are planned and are intentionally not linked until
they exist:

- Testing analysis plugin contributions — how-to guide.
- Pre-PR plugin checklist — reference.
- Stable API map for analysis contributors — reference.

Until those pages are written, use
[Extend PolyzyMD with MDAnalysis-native analyses](../extending_analyses.md) as
the complete implementation guide.

```{toctree}
:hidden:
:maxdepth: 1

How PolyzyMD analysis plugins work <architecture>
Scaffold your first analysis plugin <first_scaffold>
Build a simple scalar analysis plugin <simple_scalar_plugin>
Store arrays and tables with sidecars <sidecars>
Convert a plugin into an advanced package <advanced_package>
```
