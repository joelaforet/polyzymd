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
:caption: Commands and configuration
:maxdepth: 1

cli_reference
configuration
physics_checks
gromacs_openmm
```

## Input data

```{toctree}
:caption: Input data
:maxdepth: 1

data_requirements
openff_pdb_ingestion
benchmarks
```

## Analysis

```{toctree}
:caption: Analysis
:maxdepth: 1

study_api
analysis_functions
analysis_protocol_report
posthoc_testing
```

(what-is-validated)=
## What is validated, and what is not

```{index} !validation, !validated, !verification
```

PolyzyMD checks that it does what you ask:

- **Each build.** The physics audit compares the built system, term by term,
  with references rebuilt from the force fields, charges, water model and ions
  that the config names, and refuses a build that differs
  ({doc}`physics_checks`).
- **Each analysis.** Each was compared, on test systems, with an independent
  calculation, such as `gmx rmsf` for RMSF and MDTraj for SASA and native
  contacts. The records below give the comparisons and their agreement.
- **The two engines.** The same system on OpenMM and on GROMACS is compared in
  {ref}`engine-parity`.
- **Every change.** An automated test suite runs on each change to the code.

PolyzyMD does not validate the science of a simulation. The force fields,
charge models and water model are your choice, and whether a result agrees
with experiment is yours to establish. A built polymer chain or a packed box is
a starting structure, not an equilibrated one: equilibrate before you analyze,
and check it ({doc}`../explanation/convergence_detection`).

## Verification records

Each analysis compared, on test systems, with an independent calculation.

```{toctree}
:caption: Verification records
:maxdepth: 1

analysis_rmsf_verification
analysis_sasa_verification
analysis_contacts_verification
analysis_native_contacts_verification
analysis_hydrogen_bonds_verification
```

## Terms and sources

```{toctree}
:caption: Terms and sources
:maxdepth: 1

glossary
references
```

## Python API

```{toctree}
:caption: Python API
:maxdepth: 2

../api/index
../api/package_api
../api/analysis_api
```
