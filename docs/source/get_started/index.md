# Get Started

Welcome to PolyzyMD! This page gets you from zero to a working setup in about
10 minutes.

## Quick Install

PolyzyMD uses [pixi](https://pixi.sh) for environment management. Install pixi
(if you don't already have it), then clone and set up the project:

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

Verify the install:

```bash
polyzymd --help
polyzymd info
```

If both commands print output without errors, you are ready to go.

For the full installation guide (including GPU cluster setup and
troubleshooting), see {doc}`../tutorials/installation`.

## Create Your First Project

PolyzyMD organizes simulations into project directories. Create one with:

```bash
polyzymd init --name my_first_project
cd my_first_project
```

This creates a directory with a template `config.yaml`, a `structures/` folder
for your input PDB files, and `job_scripts/` and `slurm_logs/` directories for
HPC workflows.

## What's Inside `config.yaml`

The `config.yaml` file controls everything about your simulation: which enzyme
to simulate, solvent conditions, thermodynamic parameters, equilibration stages,
and production settings. You edit this file to describe your system, then
PolyzyMD handles the rest.

For a detailed walkthrough of writing your first config and running a build, see
{doc}`../tutorials/quickstart`.

## Where to Go Next

Pick the path that matches your situation:

::::{grid} 2
:gutter: 3

:::{grid-item-card} I want to set up and run a simulation
:link: ../tutorials/quickstart
:link-type: doc

Follow the first simulation tutorial step by step.
:::

:::{grid-item-card} I have trajectories and want to analyze them
:link: ../tutorials/first_analysis
:link-type: doc

Run your first analysis in five steps.
:::

:::{grid-item-card} I want to compare multiple conditions
:link: ../tutorials/analysis_complete_workflow
:link-type: doc

Set up a multi-condition comparison study.
:::

:::{grid-item-card} I need a specific task done
:link: ../how_to/index
:link-type: doc

Jump to how-to guides for polymers, restraints, GROMACS export, SLURM, and more.
:::

::::

```{toctree}
:hidden:
:maxdepth: 1

Install PolyzyMD with pixi <../tutorials/installation>
Run Your First PolyzyMD Simulation <../tutorials/quickstart>
```
