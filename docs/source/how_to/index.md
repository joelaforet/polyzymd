# How-To Guides

How-to guides help you accomplish a specific task. Use this section when you
already know your goal and need the shortest practical route to it.

```{tip}
If you're new to PolyzyMD, start with {doc}`../tutorials/index` instead.
These guides assume you already have a working installation and know the
basics.
```

::::{grid} 2
:gutter: 3

:::{grid-item-card} Simulation Setup
:link: polymers
:link-type: doc

Add polymers, prepare inputs, configure restraints, or export to GROMACS.
:::

:::{grid-item-card} Build Input Troubleshooting
:link: troubleshoot_openff_pdb_ingestion
:link-type: doc

Diagnose OpenFF PDB ingestion failures before running a PolyzyMD build.
:::

:::{grid-item-card} Analysis Workflows
:link: analysis_chooser
:link-type: doc

Pick the right analysis plugin and get results in a few commands.
:::

:::{grid-item-card} HPC & Job Submission
:link: hpc_slurm
:link-type: doc

Submit simulations or analysis jobs via SLURM.
:::

:::{grid-item-card} Plots & Troubleshooting
:link: publication_plots
:link-type: doc

Adjust themes, colors, and figure sizes for publication-quality output.
:::

::::

## Simulation Setup & Input Preparation

Prepare build inputs and configure your enzyme-polymer conjugate systems.

```{toctree}
:maxdepth: 1

Add Polymers to a Simulation <polymers>
Generate Polymers from SMILES <dynamic_polymers>
Troubleshoot OpenFF PDB Ingestion <troubleshoot_openff_pdb_ingestion>
Add Distance Restraints <restraints>
Set Up Equilibration Stages <equilibration>
Run GROMACS Simulations on HPC Clusters <gromacs_export>
```

## HPC & Job Submission

Run on SLURM clusters — from individual jobs to daisy-chained workflows.

```{toctree}
:maxdepth: 1

Run Simulations on SLURM Clusters <hpc_slurm>
Run OpenMM on Other Hardware <hardware_platforms>
Submit Analysis Jobs to SLURM <hpc_execution>
```

## Stable Analysis Workflows

These analysis plugins are well-tested and recommended for production use.
See also {doc}`../tutorials/sasa_analysis` for a guided SASA walkthrough.

```{toctree}
:maxdepth: 1

Which Analysis Should I Run? <analysis_chooser>
Compare Simulation Conditions <analysis_compare_conditions>
Run RMSD Analysis <analysis_rmsd_quickstart>
Run Rg Analysis <analysis_rg_quickstart>
Run RMSF Analysis <analysis_rmsf_quickstart>
Run Distance Analysis <analysis_distances_quickstart>
Run Contacts Analysis <analysis_contacts_quickstart>
Run SASA Analysis <analysis_sasa_quickstart>
Analyze Hydrogen Bonds <hydrogen_bonds>
Run Catalytic Triad Analysis <analysis_triad_quickstart>
```

## Plots & Troubleshooting

Customize publication-quality figures and debug common issues.

```{tip}
Want to use the PolyzyMD artifacts to make your own plots? See
{doc}`custom_artifact_plotting` for a post-processing workflow that reads cached
analysis artifacts and sidecars. Start with {doc}`publication_plots` if you want
to customize PolyzyMD's standard analysis plots.
```

```{toctree}
:maxdepth: 1

Customizing Plots for Publication <publication_plots>
Create Custom Plots from Analysis Artifacts <custom_artifact_plotting>
Broken Molecule Debugging <broken_molecules_debugging>
Troubleshoot Common Problems <troubleshooting>
```
