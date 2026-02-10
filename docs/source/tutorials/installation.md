# Installation

This guide covers installing PolyzyMD and its dependencies.

## Prerequisites

- **Python 3.10+**
- **conda or mamba** (required for OpenMM/OpenFF dependencies)

```{note}
OpenMM, OpenFF Toolkit, and other core simulation dependencies are **only available via conda/mamba**, not PyPI. A conda-based environment is required.
```

## User Installation

### Step 1: Create Conda Environment

#### Option A: Using mamba (Recommended)

[Mamba](https://mamba.readthedocs.io/) is a faster drop-in replacement for conda:

```bash
mamba create -n polyzymd-env python=3.11 openmm openff-toolkit openff-interchange \
    openff-nagl openff-nagl-models openff-forcefields openff-units packmol \
    mbuild mdtraj numpy scipy pandas pydantic pyyaml click tqdm -c conda-forge
mamba activate polyzymd-env
```

#### Option B: Using conda

```bash
conda create -n polyzymd-env python=3.11 openmm openff-toolkit openff-interchange \
    openff-nagl openff-nagl-models openff-forcefields openff-units packmol \
    mbuild mdtraj numpy scipy pandas pydantic pyyaml click tqdm -c conda-forge
conda activate polyzymd-env
```

#### Option C: Using our environment file

For a fully reproducible environment matching our CI:

```bash
# Using mamba (recommended)
mamba env create -f https://raw.githubusercontent.com/joelaforet/polyzymd/main/devtools/conda-envs/polyzymd-env.yml
mamba activate polyzymd-env

# Or using conda
conda env create -f https://raw.githubusercontent.com/joelaforet/polyzymd/main/devtools/conda-envs/polyzymd-env.yml
conda activate polyzymd-env
```

### Step 2: Install PolyzyMD

With your environment activated:

```bash
pip install polyzymd
```

To install a specific version:

```bash
pip install polyzymd==1.0.0
```

### Verify Installation

```bash
# Check CLI is available
polyzymd --help

# Check version and dependencies
polyzymd info
```

Expected output:

```
PolyzyMD - Molecular Dynamics for Enzyme-Polymer Systems
Version: 1.0.0

Dependencies:
  OpenMM: 8.1.1
  OpenFF Toolkit: 0.16.x
  OpenFF Interchange: 0.4.x
  Pydantic: 2.x.x

Example configs: polyzymd/configs/examples/
```

## Developer Installation

For contributing to PolyzyMD or modifying the source code:

```bash
# Clone the repository
git clone https://github.com/joelaforet/polyzymd.git
cd polyzymd

# Create environment from our environment file
mamba env create -f devtools/conda-envs/polyzymd-env.yml
mamba activate polyzymd-env

# Install in editable mode with dev dependencies
pip install -e ".[dev]"

# (Optional) Install pre-commit hooks
pre-commit install
```

## HPC Installation

```{warning}
The HPC/SLURM functionality in PolyzyMD was designed with **CU Boulder's Alpine and Blanca** clusters in mind. The SLURM presets and job submission scripts **will not work out of the box** for other systems like Bridges2 or XSEDE resources. Support for additional HPC systems is planned for a future release.
```

On HPC systems, you typically need to:

### 1. Load the anaconda/mamba module

```bash
module load anaconda
# or: module load mambaforge
# or: module load miniforge
```

### 2. Create the environment

```bash
# Navigate to your projects directory (long-term storage)
cd /projects/$USER

# Create the environment
mamba env create -f https://raw.githubusercontent.com/joelaforet/polyzymd/main/devtools/conda-envs/polyzymd-env.yml
```

### 3. Install PolyzyMD

```bash
mamba activate polyzymd-env
pip install polyzymd
```

### 4. Verify installation

```bash
polyzymd info
```

### CU Boulder Specific Notes

On Alpine/Blanca, you may want to install to a shared project space:

```bash
# Load modules
module load mambaforge

# Create environment in project space
mamba env create -f https://raw.githubusercontent.com/joelaforet/polyzymd/main/devtools/conda-envs/polyzymd-env.yml \
    -p /projects/$USER/envs/polyzymd-env

# Activate with full path
mamba activate /projects/$USER/envs/polyzymd-env

# Install polyzymd
pip install polyzymd
```

## Note for uv Users

[uv](https://github.com/astral-sh/uv) is a fast Python package manager, but it only works with packages available on PyPI. Unfortunately, **OpenMM, OpenFF Toolkit, and other core simulation dependencies are not available on PyPI** — they must be installed via conda or mamba.

If you prefer uv for package management, you can:

1. Create a conda/mamba environment with the simulation stack (OpenMM, OpenFF, etc.)
2. Activate that conda environment
3. Use uv for additional pure-Python packages within that environment

However, `pip install polyzymd` within an activated conda environment is the simplest and recommended approach.

## Optional Dependencies

### AmberTools (for AM1-BCC charging)

PolyzyMD defaults to NAGL for partial charge assignment, which is fast and doesn't require additional dependencies. If you need AM1-BCC charges via AmberTools:

```bash
mamba install -c conda-forge ambertools
```

```{note}
AmberTools has limited availability on macOS and can cause dependency conflicts. We recommend using a separate environment if needed.
```

Then use in your configuration:

```yaml
substrate:
  charge_method: "am1bcc"
```

Or in code:

```python
from polyzymd.utils.charging import get_charger

charger = get_charger("am1bcc", toolkit="ambertools")
charged_mol = charger.charge_molecule(molecule)
```

### OpenEye Toolkit (commercial, faster AM1-BCC)

If you have an OpenEye license:

```bash
mamba install -c openeye openeye-toolkits
```

### MDAnalysis (for advanced trajectory analysis)

```bash
mamba install -c conda-forge mdanalysis
```

## Common Installation Issues

### "Module not found: openmm"

OpenMM must be installed via conda, not pip:

```bash
mamba install -c conda-forge openmm
```

### "polyzymd: command not found"

Make sure you've activated the correct environment:

```bash
mamba activate polyzymd-env
which polyzymd  # Should show path in your env
```

### Environment solver takes forever

Use mamba instead of conda for faster environment solving:

```bash
# Install mamba if you don't have it
conda install -n base -c conda-forge mamba

# Then use mamba for environment creation
mamba env create -f polyzymd-env.yml
```

### Conflicts with existing environment

Create a fresh environment rather than installing into an existing one:

```bash
mamba create -n polyzymd-env --clone base  # Don't do this
mamba env create -f polyzymd-env.yml       # Do this instead
```

## Updating PolyzyMD

To update to the latest version:

```bash
mamba activate polyzymd-env
pip install --upgrade polyzymd
```

## Uninstalling

```bash
# Remove the package
pip uninstall polyzymd

# Remove the entire environment (optional)
mamba deactivate
mamba env remove -n polyzymd-env
```

## Next Steps

Once installed, continue to the {doc}`quickstart` guide to run your first simulation.
