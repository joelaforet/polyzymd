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
4. **Use the current full analysis extension guide.** The supported detailed
   guide is [Extend PolyzyMD with MDAnalysis-native analyses](../extending_analyses.md).
   It covers the current scaffold, public imports, `MDAAnalysisJob`, collectors,
   artifacts, default comparison, and artifact-only plotting.
5. **Look up stable API surfaces when needed.** Use the API pages for
   [analysis plugins](../../api/analyses.md),
   [analysis base classes](../../api/analyses_base.md),
   [MDAnalysis integration](../../api/analyses_mda.md), and
   [shared analysis utilities](../../api/analyses_shared.md).
6. **Check CLI and configuration reference material.** Use the
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

## Planned path pages

The New Analysis Contributor Path will grow into smaller pages with one purpose
each. These pages are planned and are intentionally not linked until they exist:

- Conceptual overview of the analysis plugin lifecycle — explanation.
- First scaffold walkthrough — tutorial.
- Simple scalar analysis plugin — tutorial.
- Sidecars for arrays, tables, and event streams — how-to guide.
- Advanced package layout for larger plugins — how-to guide.
- Testing analysis plugin contributions — how-to guide.
- Pre-PR plugin checklist — reference.
- Stable API map for analysis contributors — reference.

Until those pages are written, use
[Extend PolyzyMD with MDAnalysis-native analyses](../extending_analyses.md) as
the complete implementation guide.
