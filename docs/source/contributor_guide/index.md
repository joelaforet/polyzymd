# Contributor Guide

This section is for developers and scientific contributors who want to extend
PolyzyMD, understand the codebase, or add new capabilities.

## Start Here

- [Set Up a Contributor Environment](setup.md)
- [Contributing to PolyzyMD](contributing.md)
- [Packaging and Distribution Notes](packaging.md)
- [Architecture](../explanation/architecture.md)

## Add an analysis

An analysis is a Python function of an MDAnalysis `AtomGroup` or `Universe`,
run on every replicate by `Study.timeseries` or `Study.per_replicate`.
[Add an analysis](adding_an_analysis.md) shows how to write one, ship it in
`polyzymd.analyses.functions`, expose it through `polyzymd analyze` and test
it; [Study API](../reference/study_api.md) describes
the study API in full.

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
Add an Analysis <adding_an_analysis>
```
