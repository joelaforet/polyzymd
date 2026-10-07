# Tutorials

Each tutorial is one lesson with a known result. Do them in this order.

1. {doc}`../get_started/quickstart`: simulate Trp-cage in water and measure its
   radius of gyration. The input files ship with PolyzyMD.
2. {doc}`first_analysis`: run two more replicates of the quickstart, and
   measure their RMSF on the study.
3. {doc}`analysis_complete_workflow`: add a condition with two short SBMA
   chains, and compare it with water.
4. {doc}`sasa_analysis`: measure how much of the protein surface the polymer
   covers.
5. {doc}`prepare_pdb_for_openff`: download a structure from the PDB and clean
   it for OpenFF.
6. {doc}`own_system`: make a project for your own protein, fill in a config,
   check it, add a second condition, and run and analyze both.

Tutorials 2 to 4 continue in the project folder of the quickstart. Their
simulations run on a laptop CPU in minutes.

```{toctree}
:hidden:
:maxdepth: 1

Run your first simulation <../get_started/quickstart>
first_analysis
analysis_complete_workflow
sasa_analysis
prepare_pdb_for_openff
own_system
```
