# Installation

This guide covers installing PolyzyMD and its dependencies using
[pixi](https://pixi.sh), a fast, reproducible package manager.

## Prerequisites

- **Linux x86-64** (the only platform with tested environments)
- **pixi** (installed below)

```{note}
PolyzyMD's simulation stack (OpenMM, OpenFF Toolkit, etc.) is only available
via conda-forge channels. Pixi handles this automatically - you do not need
conda or mamba installed.
```

## Install pixi

If you don't already have pixi:

```bash
curl -fsSL https://pixi.sh/install.sh | sh
```

Then restart your shell or run `source ~/.bashrc` (or `~/.zshrc`).

## User Installation (local / no GPU)

Use the **build** environment for system building, validation, and PDB
preparation on a local machine without a GPU.

```bash
# Clone the repository
git clone https://github.com/joelaforet/polyzymd.git
cd polyzymd

# Install the build environment
pixi install -e build

# Activate the environment
pixi shell -e build

# Verify installation
polyzymd --help
polyzymd info
```

Inside the pixi shell, the `polyzymd` CLI is on PATH and works directly.

### What works in the build environment

| Command | Works? | Notes |
|---------|--------|-------|
| `polyzymd build` | Yes | System building (no GPU needed) |
| `polyzymd validate` | Yes | Config validation |
| `polyzymd init` | Yes | Project scaffolding |
| `polyzymd info` | Yes | Version/dependency info |
| `polyzymd clean-pdb` | Yes | PDB cleanup |
| `polyzymd submit` | No | Requires a CUDA environment |
| `polyzymd run-segment` | No | Requires a CUDA environment |

## HPC Installation

```{warning}
The SLURM presets in PolyzyMD are designed for **CU Boulder Alpine/Blanca**
and **PSC Bridges2**. Other clusters will need custom presets — see
{doc}`hpc_slurm` for details.
```

On HPC clusters with GPUs, use one of the **CUDA environments**. Pick the
environment that matches your cluster's CUDA driver version.

### How to find your CUDA version

On a GPU node, run:

```bash
nvidia-smi | head -1
```

The driver version in the top-right maps to a maximum CUDA version:

| Cluster | CUDA driver | Environment |
|---------|-------------|-------------|
| PSC Bridges2 | 12.6 | `cuda-12-6` |
| CU Boulder Blanca | 12.4 | `cuda-12-4` |

### Setup steps

```bash
# 1. Install pixi (if not already installed)
curl -fsSL https://pixi.sh/install.sh | sh
source ~/.bashrc

# 2. Clone the repo to your projects directory
cd /projects/$USER
git clone https://github.com/joelaforet/polyzymd.git
cd polyzymd

# 3. Install the environment matching your cluster
pixi install -e cuda-12-6    # Bridges2
# or
pixi install -e cuda-12-4    # CU Boulder Blanca

# 4. Verify installation (on a GPU node)
pixi shell -e cuda-12-6
polyzymd info
```

```{tip}
You only need to run `pixi install` once. After that, `pixi shell -e <env>`
is all you need to activate the environment in each new session.
```

### Submitting jobs

Inside the pixi shell, submission works as usual:

```bash
pixi shell -e cuda-12-6

polyzymd submit -c config.yaml \
    --replicates 1-3 \
    --preset aa100 \
    --dry-run
```

The generated SLURM scripts automatically activate the correct pixi
environment via `pixi shell-hook` — no `module load` or `conda activate`
needed in job scripts.

## Developer Installation

For contributing to PolyzyMD or modifying the source code:

```bash
# Clone the repository
git clone https://github.com/joelaforet/polyzymd.git
cd polyzymd

# Install the build environment (editable mode is automatic via pixi.toml)
pixi install -e build

# Activate
pixi shell -e build

# Run tests
pytest tests/ -x -q --tb=short --ignore=tests/test_packmol.py

# Run linting
ruff check src/
ruff format --check src/
```

PolyzyMD is installed in editable mode by default (configured in `pixi.toml`
via `polyzymd = { path = ".", editable = true }`). Any source changes are
immediately reflected without reinstalling.

## Optional Dependencies

### AmberTools

AmberTools is included in all pixi environments by default (pinned to
version 23.3 for compatibility with the OpenFF stack). No additional
installation is needed.

To use AM1-BCC charges via AmberTools in your configuration:

```yaml
substrate:
  charge_method: "am1bcc"
```

### OpenEye Toolkit (commercial, faster AM1-BCC)

If you have an OpenEye license, install the toolkit in your pixi environment:

```bash
pixi shell -e build
pip install openeye-toolkits
```

## Common Installation Issues

### "pixi: command not found"

Make sure pixi is installed and on PATH:

```bash
curl -fsSL https://pixi.sh/install.sh | sh
source ~/.bashrc
which pixi
```

### "polyzymd: command not found"

Make sure you've activated the pixi environment:

```bash
pixi shell -e build        # local use
pixi shell -e cuda-12-6    # HPC use
which polyzymd             # should show path inside .pixi/
```

### pixi install fails on HPC login node

Some HPC systems restrict network access on login nodes. Try:

1. Use a compile/data-transfer node if available
2. Or pre-download on a machine with internet access and `rsync` the
   repository (including `.pixi/`) to the cluster

### OpenMM can't find CUDA

Make sure you're using the correct CUDA environment for your cluster:

```bash
# Check your GPU's CUDA version
nvidia-smi | head -1

# Use the matching environment
pixi shell -e cuda-12-6    # NOT cuda-12-4 on a 12.6 system
python -c "import openmm; print(openmm.Platform.getPlatformByName('CUDA').getName())"
```

### Adding support for a new cluster

If your cluster has a different CUDA driver version, you can add a new
environment to `pixi.toml` following the pattern of the existing CUDA
features. See the comments in `pixi.toml` for instructions, or file a PR.

## Updating PolyzyMD

To update to the latest version:

```bash
cd /path/to/polyzymd
git pull
pixi install -e <your-environment>
```

Since PolyzyMD is installed in editable mode from the local clone,
`git pull` picks up all code changes. Running `pixi install` again
updates any changed conda dependencies.

## Uninstalling

```bash
# Remove the pixi environments
cd /path/to/polyzymd
rm -rf .pixi

# Or remove the entire clone
rm -rf /path/to/polyzymd
```

## Next Steps

Once installed, continue to the {doc}`quickstart` guide to run your first simulation.
