# Contributor Guide

This section is for developers and scientific contributors who want to extend
PolyzyMD, understand the codebase, or add new capabilities.

## Start Here

- [Set Up a Contributor Environment](setup.md)
- [Contributing to PolyzyMD](../tutorials/contributing.md)
- [Packaging and Distribution Notes](../tutorials/packaging.md)
- [Architecture](../tutorials/architecture.md)

## Extension Workflows

- **[Extend the Analysis Framework](../tutorials/extending_analyses.md)** —
  how to add a new analysis type in one file, including compute, comparison,
  formatting, and plotting

## Contributor Mindset

Use this section when you need to understand internal design, extension points,
or project maintenance patterns. For command lookup, switch to
[Reference](../reference/index.md).

<!-- IMAGE OPPORTUNITY: Add a high-level package architecture diagram showing
config -> builders -> simulation -> workflow -> analysis -> analyses -> compare. -->

```{toctree}
:hidden:
:maxdepth: 1

Contributing to PolyzyMD <../tutorials/contributing>
Set Up a Contributor Environment <setup>
Packaging and Distribution Notes <../tutorials/packaging>
Extend the Analysis Framework <../tutorials/extending_analyses>
Comparison Logic in an Analysis Plugin <../tutorials/extending_comparators>
Plotting in an Analysis Plugin <../tutorials/extending_plotters>
```
