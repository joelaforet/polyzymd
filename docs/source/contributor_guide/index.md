# Contributor Guide

This section is for developers and scientific contributors who want to extend
PolyzyMD, understand the codebase, or add new capabilities.

## Start Here

- [Set Up a Contributor Environment](setup.md)
- [Contributing to PolyzyMD](contributing.md)
- [Packaging and Distribution Notes](packaging.md)
- [Architecture](../explanation/architecture.md)

## New Analysis Contributor Path

Use the [New Analysis Contributor Path](analysis_plugins/index.md) if you want
to add a new analysis plugin or modify an existing one. It gives new analysis
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
