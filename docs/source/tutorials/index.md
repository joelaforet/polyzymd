# Tutorials

Each tutorial is one lesson with a known result. Do them in this order.

1. {doc}`../get_started/quickstart`: simulate Trp-cage in water and measure its
   radius of gyration. The input files ship with PolyzyMD.
2. {doc}`prepare_pdb_for_openff`: download a structure from the PDB and clean
   it for OpenFF.
3. {doc}`own_system`: make a project for your own protein, fill in a config,
   check it, add a second condition, and run and analyze both.
4. {doc}`first_analysis`: run RMSF on your own finished simulations and read
   the stored result.
5. {doc}`analysis_complete_workflow`: compare conditions of a study and make
   the figures.
6. {doc}`sasa_analysis`: compare the solvent-accessible surface of the protein
   with and without polymer.

Tutorials 4 to 6 need finished production simulations of your own.

```{toctree}
:hidden:
:maxdepth: 1

prepare_pdb_for_openff
own_system
first_analysis
analysis_complete_workflow
sasa_analysis
```
