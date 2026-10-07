# Get Started

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
| Analyze trajectories and make plots | `analysis` | `pixi run -e analysis polyzymd analyze rmsf -c config.yaml --eq 10ns` |
| Run OpenMM inside a SLURM job on an NVIDIA GPU | a `sim-cuda-*` environment | The job scripts of `polyzymd submit` select it |

To activate an environment once, instead of a prefix on each command, use
`pixi shell -e <env>`.

The SLURM job scripts for OpenMM need an NVIDIA GPU. For a CPU, an AMD GPU,
another cluster or a new driver, see {doc}`../how_to/hardware_platforms`.

## Where to go next

::::{grid} 2
:gutter: 3

:::{grid-item-card} Run a first simulation
:link: quickstart
:link-type: doc

Simulate Trp-cage in water with the shipped example, and measure its radius
of gyration.
:::

:::{grid-item-card} Analyze finished simulations
:link: ../tutorials/first_analysis
:link-type: doc

Run your first analysis in three steps.
:::

:::{grid-item-card} Compare conditions
:link: ../tutorials/analysis_complete_workflow
:link-type: doc

Analyze a study with several conditions.
:::

:::{grid-item-card} Do a specific task
:link: ../how_to/index
:link-type: doc

Find how-to guides for polymers, restraints, GROMACS, SLURM and more.
:::

::::

```{toctree}
:hidden:
:maxdepth: 1

Install PolyzyMD with pixi <installation>
Run your first simulation <quickstart>
```
