# Reference

Reference pages are for lookup. Use this section when you need commands,
configuration fields, plugin settings, API signatures, or benchmark data.

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

:::{grid-item-card} Protein Modification Config
:link: protein_modification_config
:link-type: doc

Lookup fields for covalent modifications and polymer conjugation.
:::

:::{grid-item-card} Analysis Plugin Reference
:link: analysis_plugin_settings
:link-type: doc

Per-plugin settings, comparison YAML schema, and post-hoc testing options.
:::

:::{grid-item-card} API Documentation
:link: ../api/index
:link-type: doc

Module-level Python API for config, builders, simulation, workflow, and
analysis.
:::

::::

## CLI, Configuration & Data

```{toctree}
:maxdepth: 1

CLI Reference <cli_reference>
Configuration Reference <configuration>
Protein Modification Config Reference <protein_modification_config>
Data Requirements & Directory Layout <data_requirements>
Benchmarks <benchmarks>
```

## Analysis Plugin Reference

```{toctree}
:maxdepth: 1

Comparison YAML Schema <comparison_yaml>
Analysis Plugin Settings Reference <analysis_plugin_settings>
Comparison and Plotting Reference <analysis_comparison_reference>
RMSD Plugin Reference <analysis_rmsd_reference>
Rg Plugin Reference <analysis_rg_reference>
RMSF Plugin Reference <analysis_rmsf_reference>
Catalytic Triad Plugin Reference <analysis_triad_reference>
Distances Plugin Reference <analysis_distances_reference>
Contacts Plugin Reference <analysis_contacts_reference>
Post-Hoc Testing Reference <posthoc_testing>
```

## API Reference

Full Python API documentation.

```{toctree}
:maxdepth: 1

API Reference <../api/index>
```
