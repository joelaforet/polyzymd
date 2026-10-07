# Reference

Reference pages list commands, configuration keys, settings, file layouts and
API signatures. Use them to look up a fact.

::::{grid} 2
:gutter: 3

:::{grid-item-card} CLI commands
:link: cli_reference
:link-type: doc

Every `polyzymd` command and option.
:::

:::{grid-item-card} Configuration
:link: configuration
:link-type: doc

Every key of `config.yaml`, with its type, default and unit.
:::

:::{grid-item-card} Analyses
:link: analysis_functions
:link-type: doc

The shipped analysis functions, the study API, the report schema and the
comparison tests.
:::

:::{grid-item-card} Glossary
:link: glossary
:link-type: doc

The terms these docs use, such as study, condition, replicate, NAGL and
n_eff.
:::

::::

## Commands and configuration

```{toctree}
:maxdepth: 1

cli_reference
configuration
gromacs_openmm
```

## Input data

```{toctree}
:maxdepth: 1

data_requirements
openff_pdb_ingestion
benchmarks
```

## Analysis

```{toctree}
:maxdepth: 1

study_api
analysis_functions
analysis_protocol_report
posthoc_testing
```

## Terms and sources

```{toctree}
:maxdepth: 1

glossary
references
```

## Python API

```{toctree}
:maxdepth: 2

../api/index
../api/package_api
../api/analysis_api
```
