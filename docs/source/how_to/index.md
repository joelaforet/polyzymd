# How-To Guides

How-to guides help you accomplish a specific task. Use this section when you
already know your goal and need the shortest practical route to it.

## Simulation Setup

- [Add Polymers to a Simulation](polymers.md)
- [Generate Polymers from SMILES](dynamic_polymers.md)
- [Add Distance Restraints](restraints.md)
- [Set Up Equilibration Stages](equilibration.md)
- [Export a System to GROMACS](gromacs_export.md)

## HPC and Job Submission

- [Run Simulations on SLURM Clusters](hpc_slurm.md)
- [Submit Analysis Jobs to SLURM](hpc_execution.md)

## Stable Analysis Workflows

These analyses are production-ready with validated statistics:

- [Compare Simulation Conditions](analysis_compare_conditions.md) -- set up
  `comparison.yaml` and run cross-condition analysis
- [Run RMSD Analysis](analysis_rmsd_quickstart.md)
- [Run Rg Analysis](analysis_rg_quickstart.md)
- [Run RMSF Analysis](analysis_rmsf_quickstart.md)
- [Run Distance Analysis](analysis_distances_quickstart.md)
- [Run Contacts Analysis](analysis_contacts_quickstart.md)
- [Run SASA Analysis](../tutorials/sasa_analysis.md)
- [Analyze Hydrogen Bonds](../tutorials/hydrogen_bonds.md)
- [Run Catalytic Triad Analysis](analysis_triad_quickstart.md)

## Experimental Analysis Workflows

```{warning}
These analyses are functional but their interpretation guidance is still
evolving. Results should be treated as exploratory.
```

- [Analyze Binding Preference](analysis_binding_preference.md)
- [Analyze Binding Free Energy](analysis_binding_free_energy.md)
- [Analyze Polymer Affinity](analysis_polymer_affinity.md)
- [Analyze Polymer Bridging](analysis_polymer_bridging.md)
- [Analyze Exposure Dynamics](analysis_exposure_dynamics.md)

## Recipes and Troubleshooting

- [Contacts Cookbook](analysis_contacts_cookbook.md) -- advanced selection
  patterns, custom criteria, and Python API usage
- [Broken Molecule Debugging](broken_molecules_debugging.md)
- [Troubleshoot Common Problems](troubleshooting.md)

```{toctree}
:hidden:
:maxdepth: 1

Add Polymers to a Simulation <polymers>
Generate Polymers from SMILES <dynamic_polymers>
Add Distance Restraints <restraints>
Set Up Equilibration Stages <equilibration>
Export a System to GROMACS <gromacs_export>
Run Simulations on SLURM Clusters <hpc_slurm>
Submit Analysis Jobs to SLURM <hpc_execution>
Compare Simulation Conditions <analysis_compare_conditions>
Run RMSD Analysis <analysis_rmsd_quickstart>
Run Rg Analysis <analysis_rg_quickstart>
Run RMSF Analysis <analysis_rmsf_quickstart>
Run Distance Analysis <analysis_distances_quickstart>
Run Contacts Analysis <analysis_contacts_quickstart>
Run SASA Analysis <../tutorials/sasa_analysis>
Analyze Hydrogen Bonds <../tutorials/hydrogen_bonds>
Run Catalytic Triad Analysis <analysis_triad_quickstart>
Analyze Binding Preference <analysis_binding_preference>
Analyze Binding Free Energy <analysis_binding_free_energy>
Analyze Polymer Affinity <analysis_polymer_affinity>
Analyze Polymer Bridging <analysis_polymer_bridging>
Analyze Exposure Dynamics <analysis_exposure_dynamics>
Contacts Cookbook <analysis_contacts_cookbook>
Broken Molecule Debugging <broken_molecules_debugging>
Troubleshoot Common Problems <troubleshooting>
```
