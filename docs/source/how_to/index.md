# How-to guides

Each guide solves one task. Use a guide when you know what you want to do.
If you are new to PolyzyMD, do the {doc}`../get_started/quickstart` first.

::::{grid} 2
:gutter: 3

:::{grid-item-card} Set up a simulation
:link: polymers
:link-type: doc

Add polymers, restraints and equilibration stages. Fix a PDB that OpenFF
refuses.
:::

:::{grid-item-card} Run on a cluster
:link: hpc_slurm
:link-type: doc

Submit simulations and analysis jobs with SLURM, and monitor them.
:::

:::{grid-item-card} Organize and publish studies
:link: study_folder
:link-type: doc

Make a study folder for each set of compared conditions, group studies into
a project, and freeze them for publication.
:::

:::{grid-item-card} Analyze trajectories
:link: analysis_chooser
:link-type: doc

Choose an analysis, run it on a study, and measure your own quantity.
:::

::::

## Set up a simulation

```{toctree}
:maxdepth: 1

polymers
dynamic_polymers
restraints
equilibration
troubleshoot_openff_pdb_ingestion
```

## Run simulations

```{toctree}
:maxdepth: 1

hpc_slurm
monitor_simulations
hardware_platforms
run_gromacs
site_cu_boulder
```

## Organize and publish studies

```{toctree}
:maxdepth: 1

study_folder
study_yaml
project
move_studies_into_project
study_freeze
```

## Analyze trajectories

```{toctree}
:maxdepth: 1

analysis_chooser
analysis_agent_protocol
analysis_compare_conditions
study_api
hpc_execution
```

## Analyses

```{toctree}
:maxdepth: 1

analysis_rmsd_quickstart
analysis_rg_quickstart
analysis_rmsf_quickstart
analysis_distances_quickstart
analysis_contacts_quickstart
analysis_sasa_quickstart
analysis_secondary_structure_quickstart
analysis_native_contacts_quickstart
hydrogen_bonds
analysis_triad_quickstart
```

## Figures

```{toctree}
:maxdepth: 1

publication_plots
custom_artifact_plotting
```

## Fix problems

```{toctree}
:maxdepth: 1

troubleshooting
```
