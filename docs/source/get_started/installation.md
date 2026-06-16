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

- Linux x86-64
- Git
- a shell where you can install `pixi`

```{note}
Linux pixi environments include AmberTools 24.8 for legacy AM1-BCC/OpenFF
backend compatibility. macOS pixi environments omit AmberTools by default
because current AmberTools builds conflict with the Python 3.12 / NumPy 2 full
stack in the cross-platform lock. On macOS, use the recommended NAGL/OpenFF
charging path or provide pre-charged molecules. If you need AmberTools/AM1-BCC
on macOS, please open an issue so a dedicated environment can be designed.
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

Use one of the CUDA environments instead of `build`.

```bash
git clone https://github.com/joelaforet/polyzymd.git
cd polyzymd

# Pick the environment that matches your cluster
pixi install -e cuda-12-6
# or
pixi install -e cuda-12-4

pixi shell -e cuda-12-6
polyzymd info
```

Built-in presets include configurations for CU Boulder Alpine/Blanca and PSC Bridges2.
For submission details and cluster-specific notes, see
{doc}`../how_to/hpc_slurm`.

### Choose the Right CUDA Environment

On a GPU node, check the driver version:

```bash
nvidia-smi | head -1
```

Use this mapping:

| CUDA driver version | Environment | Example clusters |
|---------------------|-------------|------------------|
| 12.6 | `cuda-12-6` | PSC Bridges2 |
| 12.4 | `cuda-12-4` | CU Boulder Blanca |

## 4. Verify the Commands You Can Use

In the `build` environment, these commands should work directly:

| Command | Works in `build`? | Notes |
|---------|-------------------|-------|
| `polyzymd init` | Yes | Project scaffolding |
| `polyzymd build` | Yes | System building |
| `polyzymd validate` | Yes | Config validation |
| `polyzymd clean-pdb` | Yes | Convenience PDB cleanup, not biological-system selection |
| `polyzymd info` | Yes | Version/dependency summary |
| `polyzymd submit` | No | Requires a CUDA environment |
| `polyzymd run-segment` | No | Requires a CUDA environment |

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
pixi shell -e cuda-12-6
python -c "import openmm; print(openmm.Platform.getPlatformByName('CUDA').getName())"
```

For broader installation and workflow issues, see {doc}`../how_to/troubleshooting`.

## Contributor Setup

If you are modifying PolyzyMD itself, use the contributor environment guide:

- {doc}`../contributor_guide/setup`

## Next Step

Continue to {doc}`quickstart` to create a project scaffold and run your first
simulation.
