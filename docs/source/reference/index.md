# Reference

Reference pages are for lookup. Use this section when you need commands,
configuration fields, analysis settings, API signatures, or benchmark data.

```{tip}
For a guided walkthrough, go to {doc}`../tutorials/index`.
For a task-oriented solution, go to {doc}`../how_to/index`.
```

::::{grid} 2
:gutter: 3

:::{grid-item-card} CLI Commands
:link: cli_reference
:link-type: doc

All `polyzymd` commands, flags, and options.
:::

:::{grid-item-card} Configuration YAML
:link: configuration
:link-type: doc

Every key in `config.yaml` — types, defaults, and constraints.
:::

:::{grid-item-card} PDB Input Requirements
:link: openff_pdb_ingestion
:link-type: doc

OpenFF chemistry requirements and PolyzyMD chain conventions for enzyme PDBs.
:::

:::{grid-item-card} Analysis Reference
:link: analysis_functions
:link-type: doc

Shipped analysis functions, the report schema, the comparison tests, and the archive of experimental analyses.
:::

:::{grid-item-card} API Documentation
:link: ../api/index
:link-type: doc

Module-level Python API for config, builders, simulation, workflow, and
analysis.
:::

::::

## CLI & Configuration

```{toctree}
:maxdepth: 1

CLI Reference <cli_reference>
Configuration Reference <configuration>
```

## Input Data & PDB Requirements

```{toctree}
:maxdepth: 1

Data Requirements & Directory Layout <data_requirements>
OpenFF PDB Ingestion Reference <openff_pdb_ingestion>
Benchmarks <benchmarks>
```

## Analysis Reference

```{toctree}
:maxdepth: 1

ProtocolReport Schema <analysis_protocol_report>
Shipped analysis functions <analysis_functions>
Comparison Tests Reference <posthoc_testing>
Experimental Analyses Archive <experimental_analyses_archive>
```

## API Reference

Full Python API documentation.

```{toctree}
:maxdepth: 2

API Reference <../api/index>
Package and workflow APIs <../api/package_api>
Analysis APIs <../api/analysis_api>
```
