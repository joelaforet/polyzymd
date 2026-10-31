# Get started

Install PolyzyMD, then run the first simulation tutorial. Together they take
about 15 minutes.

## Quick install

PolyzyMD uses [pixi](https://pixi.sh) to manage its environments. Install
pixi, then clone PolyzyMD and install the `build` environment:

```bash
# Install pixi
curl -fsSL https://pixi.sh/install.sh | sh
source ~/.bashrc

# Clone and install PolyzyMD
git clone https://github.com/joelaforet/polyzymd.git
cd polyzymd
pixi install -e build
pixi shell -e build
```

Check the install:

```bash
polyzymd --help
polyzymd info
```

If both commands print output and no error, the install works. For GPU
clusters and installation problems, see {doc}`installation`.

## Which pixi environment to use

PolyzyMD has one environment for each kind of work:

| Task | Environment | Example command |
|------|-------------|-----------------|
| Validate configs, build systems, run a short simulation on this machine | `build` | `pixi run -e build polyzymd run -c config.yaml -r 1` |
| Submit OpenMM simulations to SLURM | `build` | Build first, in an interactive job: `pixi run -e build polyzymd build -c config.yaml -r 1-3`; then `pixi run -e build polyzymd submit -c config.yaml -r 1-3 --preset aa100 --pixi-env auto` |
| Analyze trajectories and make plots | `build`, or the smaller `analysis` | `pixi run -e analysis polyzymd analyze rmsf --study my_study` |
| Run OpenMM inside a SLURM job on an NVIDIA GPU | a `sim-cuda-*` environment | The job scripts of `polyzymd submit` select it |
| Run OpenMM on the NVIDIA GPU of this machine, after `polyzymd build` in `build` | a `sim-cuda-*` environment | `pixi run -e sim-cuda-12-6 polyzymd run -c config.yaml -r 1` |

The `build` environment includes the analysis tools, so `polyzymd analyze`
also runs there. To activate an environment once, instead of a prefix on each
command, use `pixi shell -e <env>`.

The example commands work inside the clone of the repository. `pixi run` and `pixi shell` look for `pixi.toml` in the current folder and
its parent folders. In a project folder outside the clone, they stop with
`could not find pixi.toml`. Give them the clone's `pixi.toml` with
`--manifest-path`:

```bash
pixi run --manifest-path ~/polyzymd/pixi.toml -e build polyzymd validate -c config.yaml
```

Or start `pixi shell -e build` once in the clone and then `cd` to your
project. The shell keeps its environment in every folder. To build in `build`
and run in a `sim-cuda-*` environment on the same machine, open one shell for
each environment, or use `--manifest-path`.

The SLURM job scripts for OpenMM need an NVIDIA GPU. For a local NVIDIA GPU,
a CPU, an AMD GPU, another cluster or a new driver, see
{doc}`../how_to/hardware_platforms`.

(scope-and-limitations)=
## Scope and limitations

PolyzyMD 1.3 builds explicit-solvent systems around one protein and runs plain
molecular dynamics of them on OpenMM or GROMACS. Check that your study fits
before you start.

**What it builds.** One protein from a PDB file, with at most one small
molecule (a ligand), one polymer entry (any number of chains of one or more
monomers) and any number of co-solvents. The solvent is TIP3P water with NaCl
ions, in a cubic or rhombic dodecahedral box. See
{doc}`../reference/configuration`.

**Force fields and charges.** Only OpenFF SMIRNOFF force fields (`.offxml`
files), applied through OpenFF Interchange, such as the OpenFF port of Amber
ff14SB for the protein and Sage for everything else. Each solute names one
file, and every molecule except the protein uses the file of the solvent.
Charges come from the force field's library (protein), files PolyzyMD ships
(`preset`), NAGL, AM1-BCC or Espaloma. PolyzyMD reads no AMBER, CHARMM or
GROMACS topologies as input.

**Dynamics.** NVT and NPT stages with Langevin dynamics on both engines.
OpenMM uses a Monte Carlo barostat; GROMACS, which has none, runs stochastic
cell rescaling in its place. Nosé-Hoover, Andersen and anisotropic pressure
coupling run only on GROMACS. Position restraints work on both engines, with
fewer groups on GROMACS; distance restraints only on OpenMM. No NVE. Each build is audited against its config
(see {doc}`../explanation/physics_audit`).

**Engines and hardware.** OpenMM 8 on a CPU or an NVIDIA GPU (CUDA), on Linux.
GROMACS, which you install yourself. The generated OpenMM SLURM jobs need an
NVIDIA GPU. macOS works for setup and analysis, without CUDA. See
{doc}`../reference/gromacs_openmm` and {doc}`../how_to/hardware_platforms`.

**Statistics.** One value per replicate, a Student-t interval, Welch's test
and the Benjamini-Hochberg correction. Every replicate starts from the same
protein structure, so the spread measures sampling from that start. There is
no block averaging, bootstrap or jackknife. See
{doc}`../explanation/analysis_statistics_best_practices`.

**Out of scope.** Enhanced sampling, free-energy calculations, membranes,
non-OpenFF force fields, water models other than TIP3P, salts other than NaCl,
and bitwise-reproducible reruns.

(shipped-examples)=
## Shipped examples

```{index} !examples, !example configs
```

The `examples/` folder of the PolyzyMD repository holds one small system for
each thing PolyzyMD builds: a protein in water, a protein with polymer chains,
a protein in a co-solvent, and a protein with a ligand and a distance
restraint. Each example runs end to end on a laptop in a few minutes, with
the same five commands as {doc}`quickstart`. Use them to see what PolyzyMD
does before you set up your own system, or as a starting config.

| Example | What it shows | Engine | Time on a laptop | Explained in |
|---|---|---|---|---|
| [`quickstart/`](https://github.com/joelaforet/polyzymd/tree/main/examples/quickstart) | Trp-cage in water with NaCl | OpenMM, GROMACS | about 2 min (CPU) | {doc}`quickstart`, {doc}`../how_to/run_gromacs` |
| [`polymer/`](https://github.com/joelaforet/polyzymd/tree/main/examples/polymer) | Trp-cage with two SBMA chains built from SMILES | OpenMM, GROMACS | about 1.5 min (OpenMM, CPU); under 1 min (GROMACS, CPU) | {doc}`../tutorials/analysis_complete_workflow`, {doc}`../how_to/dynamic_polymers` |
| [`cosolvent/`](https://github.com/joelaforet/polyzymd/tree/main/examples/cosolvent) | Trp-cage in 2 M urea | OpenMM | about 1.5 min (CPU) | {doc}`../tutorials/own_system`, {doc}`../reference/configuration` |
| [`ligand_restraint/`](https://github.com/joelaforet/polyzymd/tree/main/examples/ligand_restraint) | T4 lysozyme L99A (PDB 181L) with benzene and a flat-bottom distance restraint | OpenMM | about 3.5 min (CPU); about 1.5 min (GPU) | {doc}`../tutorials/prepare_pdb_for_openff`, {doc}`../how_to/restraints` |
| [`modified_proteins/`](https://github.com/joelaforet/polyzymd/tree/main/examples/modified_proteins) | Ubiquitin with a GLYCAM06j glycan on Asn60, and with an octanoyl amide on Lys63 | OpenMM | about 2 min (GPU), build about 1 min | {doc}`../tutorials/glycosylated_protein`, {doc}`../how_to/attach_groups` |

The times are for the shipped configs, build included, measured on a laptop
with a 24-thread CPU and, where given, an RTX 5070 Ti laptop GPU. Each run is
a few picoseconds: long enough to check that the system builds, passes the
physics audit and runs, not long enough for a result.

### Run an example

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

### What to check in each run

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

### Make an example a real run

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

## Where to go next

::::{grid} 2
:gutter: 3

:::{grid-item-card} Run a first simulation
:link: quickstart
:link-type: doc

Simulate Trp-cage in water with the shipped example, and measure its radius
of gyration.
:::

:::{grid-item-card} Analyze a study
:link: ../tutorials/first_analysis
:link-type: doc

Run more replicates of the quickstart and measure their RMSF.
:::

:::{grid-item-card} Compare conditions
:link: ../tutorials/analysis_complete_workflow
:link-type: doc

Add a condition with a polymer and compare it with water.
:::

:::{grid-item-card} Your own analyses and molecules

Measure your own quantity on every replicate ({doc}`../how_to/study_api`,
{doc}`../how_to/own_analysis`), build your own polymers and co-solvents
({doc}`../how_to/dynamic_polymers`, {doc}`../how_to/cosolvents`), or run a
build with your own protocol ({doc}`../how_to/own_protocol`).
:::

:::{grid-item-card} Start from an example
:link: shipped-examples
:link-type: ref

The shipped configs: protein in water, with a polymer, with a co-solvent,
and with a ligand and a restraint.
:::

:::{grid-item-card} Do a specific task
:link: ../how_to/index
:link-type: doc

Find how-to guides for polymers, restraints, GROMACS, SLURM and more.
:::

::::

```{toctree}
:caption: Where to go next
:hidden:
:maxdepth: 1

Install PolyzyMD with pixi <installation>
Run your first simulation <quickstart>
```
