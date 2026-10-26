# Tutorials

Each tutorial is one lesson with a known result. Do them in this order. The
times are for a laptop CPU and include the simulations.

## Pick your track

::::{grid} 1 1 2 2
:gutter: 2

:::{grid-item-card} New to PolyzyMD
:link: ../get_started/quickstart
:link-type: doc

Install, then do tutorials 1 to 4 in order. They run on a laptop CPU.
:::

:::{grid-item-card} GROMACS user
:link: ../get_started/quickstart
:link-type: doc

Do tutorial 1 with `config_gromacs.yaml`, as its GROMACS tip says. Then read
Run GROMACS simulations and the GROMACS and OpenMM reference.
:::

:::{grid-item-card} Cluster user
:link: ../get_started/installation
:link-type: doc

Read the HPC part of the installation, do tutorial 1 on a laptop or in a
compute job, then read Run on a SLURM cluster.
:::

:::{grid-item-card} Extending PolyzyMD
:link: ../contributor_guide/index
:link-type: doc

Do tutorials 1 to 3, then read the architecture and add your own analysis,
co-solvent or check, as the contributor guide shows.
:::

::::

The GROMACS track continues in {doc}`../how_to/run_gromacs` and
{doc}`../reference/gromacs_openmm`; the cluster track in
{doc}`../how_to/hpc_slurm`.

## The learning path

1. **{doc}`../get_started/quickstart`**
   - *You learn:* the project → study → condition structure, and the
     validate → run → analyze loop.
   - *You run:* Trp-cage in water, one replicate, and its radius of gyration.
   - *Time:* about 10 minutes, 2 of them simulating.
   - *Engines:* OpenMM on the CPU as shipped. For GROMACS, add
     `config_gromacs.yaml` instead. For an NVIDIA GPU, see the environment
     tip of the tutorial.
2. **{doc}`first_analysis`**
   - *You learn:* replicates, and how a study-level analysis reports a mean
     with its interval and stores its results.
   - *You run:* replicates 2 and 3 of the quickstart, and RMSF on the study.
   - *Time:* about 10 minutes, 3 of them simulating.
   - *Engines:* the engine of your quickstart condition.
3. **{doc}`analysis_complete_workflow`**
   - *You learn:* how a second condition is compared with the control, and
     how to read the comparison and its warnings.
   - *You run:* Trp-cage with two short SBMA chains, three replicates, and
     four analyses with one command.
   - *Time:* about 20 minutes, 5 to 7 of them building and simulating.
   - *Engines:* OpenMM or GROMACS. The copied condition keeps the engine of
     `Water`.
4. **{doc}`sasa_analysis`**
   - *You learn:* how the context of a SASA calculation shows polymer
     shielding, and how to compute a new quantity in Python.
   - *You run:* two SASA analyses and a short Python script. No new
     simulation.
   - *Time:* about 10 minutes.
   - *Engines:* either; it reads the runs of tutorial 3.
5. **{doc}`prepare_pdb_for_openff`**
   - *You learn:* what OpenFF needs from a PDB file, and how to get there
     from a crystal structure.
   - *You run:* the cleanup of ubiquitin (1UBQ) and a test build.
   - *Time:* about 20 minutes.
   - *Engines:* either; the cleaned file serves both.
6. **{doc}`own_system`**
   - *You learn:* how to start a project for your own protein from a
     template config, and add a co-solvent condition.
   - *You run:* a water and a 2 M urea condition, with dry runs first.
   - *Time:* about 30 minutes to set up. The simulations take as long as your
     protein and durations need.
   - *Engines:* the template is for OpenMM. Set `engine: gromacs` for
     GROMACS.
7. **{doc}`star_sited_monomers`**
   - *You learn:* how to write the star-sited SMILES that dynamic polymers
     are built from.
   - *You run:* a notebook that draws your monomer and writes its SMILES.
   - *Time:* about 15 minutes.
   - *Engines:* none; no simulation.

Tutorials 2 to 4 continue in the project folder of the quickstart.

For one task at a time, such as {doc}`co-solvents <../how_to/cosolvents>`,
{doc}`restraints <../how_to/restraints>` or
{doc}`polymers from SMILES <../how_to/dynamic_polymers>`, see the
{doc}`../how_to/index`. For the shipped configs of each kind of system, see
{doc}`../get_started/examples`.

```{toctree}
:hidden:
:maxdepth: 1

Run your first simulation <../get_started/quickstart>
first_analysis
analysis_complete_workflow
sasa_analysis
prepare_pdb_for_openff
own_system
star_sited_monomers
```
