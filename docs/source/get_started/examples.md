# Run the examples

The `examples/` folder of the PolyzyMD repository holds one small system for
each thing PolyzyMD builds: a protein in water, a protein with polymer chains,
a protein in a co-solvent, and a protein with a ligand and a distance
restraint. Each example runs end to end on a laptop in a few minutes, with
the same five commands as {doc}`quickstart`. Use them to see what PolyzyMD
does before you set up your own system, or as a starting config.

## The examples

| Example | What it shows | Engine | Time on a laptop | Explained in |
|---|---|---|---|---|
| [`quickstart/`](https://github.com/joelaforet/polyzymd/tree/main/examples/quickstart) | Trp-cage in water with NaCl | OpenMM, GROMACS | about 2 min (CPU) | {doc}`quickstart`, {doc}`../how_to/run_gromacs` |
| [`polymer/`](https://github.com/joelaforet/polyzymd/tree/main/examples/polymer) | Trp-cage with two SBMA chains built from SMILES | OpenMM, GROMACS | about 1.5 min (OpenMM, CPU); under 1 min (GROMACS, CPU) | {doc}`../tutorials/analysis_complete_workflow`, {doc}`../how_to/dynamic_polymers` |
| [`cosolvent/`](https://github.com/joelaforet/polyzymd/tree/main/examples/cosolvent) | Trp-cage in 2 M urea | OpenMM | about 1.5 min (CPU) | {doc}`../tutorials/own_system`, {doc}`../reference/configuration` |
| [`ligand_restraint/`](https://github.com/joelaforet/polyzymd/tree/main/examples/ligand_restraint) | T4 lysozyme L99A (PDB 181L) with benzene and a flat-bottom distance restraint | OpenMM | about 3.5 min (CPU); about 1.5 min (GPU) | {doc}`../tutorials/prepare_pdb_for_openff`, {doc}`../how_to/restraints` |

The times are for the shipped configs, build included, measured on a laptop
with a 24-thread CPU and, where given, an RTX 5070 Ti laptop GPU. Each run is
a few picoseconds: long enough to check that the system builds, passes the
physics audit and runs, not long enough for a result.

## Run an example

Each folder has a `README.md` with its exact commands. They follow the
quickstart: make a project, add the example config as a condition, validate,
run and analyze. For the polymer example, from the `examples/` folder of your
clone, in the `build` environment:

```bash
polyzymd project init paper --study trpcage
polyzymd study add-condition SBMA --config polymer/config.yaml --study paper/trpcage
polyzymd validate -c paper/trpcage/conditions/sbma/config.yaml
polyzymd run -c paper/trpcage/conditions/sbma/config.yaml -r 1
```

`add-condition` copies the config and the input files it names into the
condition folder, so the project does not depend on the repository. The run
goes into `paper/runs/trpcage/sbma/`. To analyze it, list the analyses in
`paper/project.yaml` (the README names them) and run
`polyzymd analyze --project paper`.

:::{admonition} Environment setup
:class: tip

Run the commands in the `build` environment (`pixi shell -e build` in the
clone). The shipped configs use the OpenMM CPU platform, so they need no GPU.
The GROMACS configs need `gmx` on `PATH`.
:::

## What to check in each run

Every build writes `build_provenance.md` into the replicate folder. Read it
after the run:

- **Components** lists each solute, the solvent and the ions, with the force
  field and the charge method of each.
- **Physics audit** must say `PASSED`. Its rows check the charges, masses and
  parameters of every component against their sources, the ion counts, the
  overlaps and the box. Each example README names the rows that matter for
  that system. See {doc}`../explanation/physics_audit`.
- **Methods to cite** lists the software, force fields and models the build
  used.
- **Inputs and seeds** gives the seeds and, for polymers, the mBuild version.

## Make an example a real run

The shipped durations are a few picoseconds. For a real study:

1. Lengthen the production `duration` (ns) and the equilibration stages. See
   {doc}`../how_to/run_longer` and {doc}`../how_to/equilibration`.
2. Run three or more replicates (`-r 1-3`).
3. Run on a GPU: set `openmm.platform: CUDA` and follow
   {doc}`../how_to/hardware_platforms`, or submit to a cluster with
   {doc}`../how_to/hpc_slurm`.

The test `tests/test_examples.py` checks that every example config validates
and names files that exist. With `POLYZYMD_RUN_EXAMPLES=1` it also runs each
example and checks that its physics audit passes.
