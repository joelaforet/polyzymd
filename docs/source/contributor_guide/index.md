# Contributor Guide

This section is for developers and scientific contributors who want to extend
PolyzyMD, understand the codebase, or add new capabilities.

## Start Here

- [Set Up a Contributor Environment](setup.md)
- [Contributing to PolyzyMD](contributing.md)
- [Packaging and Distribution Notes](packaging.md)
- [Architecture](../explanation/architecture.md)

## New Analysis Contributor Path

No shipped analysis uses the plugin framework any more, and the framework is
being removed: a new analysis is a plain function of MDAnalysis atom groups or a
`Universe` run through `Study.timeseries` or `Study.per_replicate`; see
[Analysing a set of simulations](../explanation/analysis_api.md).

Use the [New Analysis Contributor Path](analysis_plugins/index.md) if you
maintain the plugin framework or a plugin of your own. It gives new analysis
contributors a recommended reading order, points to the current full extension
guide, and highlights the public APIs to use before you start editing plugin
code.

Start there when you need to understand how a PolyzyMD analysis moves from
MDAnalysis trajectory work to replicate artifacts, condition aggregation,
comparison, and artifact-only plotting.

## Extension Workflows

- **[Extend the Analysis Framework](extending_analyses.md)** —
  how to add an MDAnalysis-native analysis plugin with `MDAAnalysisJob`,
  `ReplicateArtifact`, default artifact aggregation, comparison, formatting,
  and artifact-only plotting
- **[Store large analysis outputs with artifact sidecars](analysis_plugins/sidecars.md)** —
  how to persist arrays, tables, and plotting inputs through registered
  `ArtifactStore` sidecars instead of plugin-specific cache files

## Contributor Mindset

Use this section when you need to understand internal design, extension points,
or project maintenance patterns. For command lookup, switch to
[Reference](../reference/index.md).

<!-- IMAGE OPPORTUNITY: Add a high-level package architecture diagram showing
config -> builders -> simulation -> workflow -> analyses. -->

```{toctree}
:hidden:
:maxdepth: 1

Contributing to PolyzyMD <contributing>
Set Up a Contributor Environment <setup>
Packaging and Distribution Notes <packaging>
New Analysis Contributor Path <analysis_plugins/index>
Extend the Analysis Framework <extending_analyses>
```
