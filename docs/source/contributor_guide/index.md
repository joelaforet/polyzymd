# Contributor Guide

This section is for developers and scientific contributors who want to extend
PolyzyMD, understand the codebase, or add new capabilities.

## Start Here

- [Set Up a Contributor Environment](setup.md)
- [Contributing to PolyzyMD](contributing.md)
- [Packaging and Distribution Notes](packaging.md)
- [Architecture](../explanation/architecture.md)

## Add an analysis

[Write an analysis plugin](analysis_plugins/index.md) is the whole path. An
analysis is one module holding a settings model and a `compute()` that returns
observables; the framework owns caching, aggregation, testing, plotting and
formatting. The [checklist](analysis_plugins/checklist.md) is what to run before
opening the pull request.

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
Write an analysis plugin <analysis_plugins/index>
```
