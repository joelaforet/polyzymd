# Install PolyzyMD with pixi

Use this page when you want a working PolyzyMD install quickly.

By the end of this guide, you will have a pixi-managed environment where
`polyzymd --help` and `polyzymd info` run successfully.

```{note}
PolyzyMD's simulation stack depends on packages distributed through
conda-forge. [pixi](https://pixi.sh) handles that solver stack for you, so you
do not need conda or mamba installed separately.
```

## Before You Start

- Linux x86-64 for CUDA simulation environments
- Git
- a shell where you can install `pixi`

```{note}
macOS is supported for non-CUDA setup, development, and trajectory analysis
where the required dependencies solve. The CUDA simulation environments
(`sim-cuda-12-4` and `sim-cuda-12-6`) are Linux-only and intended for GPU
clusters.
```

```{note}
AmberTools/AM1-BCC support in PolyzyMD's pixi workflow is Linux-focused.
AmberTools is not part of the default v1.3 NumPy 2 build/analysis/test solve,
and it is not expected to be available in all macOS solves. Use the recommended
NAGL/OpenFF charging path or provide pre-charged molecules. If you need
AmberTools/AM1-BCC, please open an issue so a dedicated environment can be
designed for that platform.
```

If you are setting up PolyzyMD on a GPU cluster, read this page first and then
use the short HPC section below.

## 1. Install pixi

If `pixi` is not already available on your system:

```bash
curl -fsSL https://pixi.sh/install.sh | sh
source ~/.bashrc
```

Check that it is available:

```bash
pixi --version
```

## 2. Install on a Local Machine

The default local setup uses the `build` environment. This is the right choice
for project scaffolding, validation, system building, documentation work, and
most local development tasks.

```bash
git clone https://github.com/joelaforet/polyzymd.git
cd polyzymd

pixi install -e build
pixi shell -e build

polyzymd --help
polyzymd info
```

If those last two commands work, your install is ready.

<!-- IMAGE OPPORTUNITY: Add a terminal screenshot showing `pixi install -e
build`, `pixi shell -e build`, and a successful `polyzymd info` output. -->

## 3. If You Are Installing on HPC

PolyzyMD v1.3 uses a split environment workflow on clusters:

- Use `build` on a login/build node to prepare and validate systems.
- Use `sim-cuda-12-4` or `sim-cuda-12-6` on GPU nodes for OpenMM simulation
  execution only.
- Use `analysis` after trajectories are produced to run comparisons and plots.

This is different from the older `cuda-*` workflow where one CUDA environment
was used for almost everything. The split is intentional: CUDA 12.4 OpenMM
builds require a NumPy 1.x runtime, while package development and analysis use
Python 3.12 with NumPy 2. Your project directory, prepared system files,
checkpoints, and trajectories are shared across environments.

Install the system-preparation environment first:

```bash
git clone https://github.com/joelaforet/polyzymd.git
cd polyzymd

pixi install -e build
pixi run -e build polyzymd validate -c config.yaml
pixi run -e build polyzymd build -c config.yaml
```

Then install the simulation runtime that matches your cluster:

```bash
# Pick the environment that matches your cluster
pixi install -e sim-cuda-12-6
# or
pixi install -e sim-cuda-12-4

pixi shell -e sim-cuda-12-6
polyzymd info
```

Built-in presets include configurations for CU Boulder Alpine/Blanca and PSC Bridges2.
For submission details and cluster-specific notes, see
{doc}`../how_to/hpc_slurm`.

After simulations finish, switch to the analysis environment:

```bash
pixi install -e analysis
pixi run -e analysis polyzymd compare validate -f comparison.yaml
pixi run -e analysis polyzymd compare run rmsf -f comparison.yaml
pixi run -e analysis polyzymd compare plot-all -f comparison.yaml
```

### Choose the Right CUDA Environment

On a GPU node, check the driver version:

```bash
nvidia-smi | head -1
```

Use this mapping:

| CUDA driver version | Environment | Example clusters |
|---------------------|-------------|------------------|
| 12.6 | `sim-cuda-12-6` | PSC Bridges2 |
| 12.4 | `sim-cuda-12-4` | CU Boulder Blanca |

## 4. Verify the Commands You Can Use

In the `build` environment, these commands should work directly:

| Command | Works in `build`? | Notes |
|---------|-------------------|-------|
| `polyzymd init` | Yes | Project scaffolding |
| `polyzymd build` | Yes | System building |
| `polyzymd validate` | Yes | Config validation |
| `polyzymd clean-pdb` | Yes | Convenience PDB cleanup, not biological-system selection |
| `polyzymd info` | Yes | Version/dependency summary |
| `polyzymd submit --engine openmm` | No | OpenMM submission targets a `sim-cuda-*` runtime |
| `polyzymd submit --engine gromacs` | Yes | GROMACS submission can run from `build`; SLURM runs GROMACS in the external cluster environment |
| `polyzymd run-segment` | No | OpenMM execution requires a `sim-cuda-*` runtime |

The `analysis` environment is the supported environment for `polyzymd compare`
commands. It contains MDAnalysis, MDTraj, pandas, SciPy, scikit-learn,
matplotlib, seaborn, Python 3.12, and NumPy 2 without the CUDA runtime pins.

## Common Installation Checks

If something looks wrong, try these first:

### `pixi: command not found`

```bash
curl -fsSL https://pixi.sh/install.sh | sh
source ~/.bashrc
which pixi
```

### `polyzymd: command not found`

```bash
pixi shell -e build
which polyzymd
```

### OpenMM cannot see CUDA

Make sure you activated the CUDA environment that matches your cluster:

```bash
pixi shell -e sim-cuda-12-6
python -c "import openmm; print(openmm.Platform.getPlatformByName('CUDA').getName())"
```

For broader installation and workflow issues, see {doc}`../how_to/troubleshooting`.

## Contributor Setup

If you are modifying PolyzyMD itself, use the contributor environment guide:

- {doc}`../contributor_guide/setup`

## Next Step

Continue to {doc}`quickstart` to create a project scaffold and run your first
simulation.
