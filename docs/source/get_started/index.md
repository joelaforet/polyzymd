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
| Submit OpenMM simulations to SLURM | `build` | `pixi run -e build polyzymd submit -c config.yaml --preset aa100 --pixi-env auto` |
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

:::{grid-item-card} Evaluate PolyzyMD
:link: evaluate
:link-type: doc

For reviewers: one run, the build record, the audit and the statistics, in
order.
:::

:::{grid-item-card} Start from an example
:link: examples
:link-type: doc

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
examples
evaluate
```
