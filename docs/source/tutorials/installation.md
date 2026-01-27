# Installation

This guide covers installing PolyzyMD and its dependencies.

## Prerequisites

PolyzyMD requires:

- **Python 3.10+**
- **OpenMM 8.0+** (for MD simulations)
- **OpenFF Toolkit 0.14+** (for force field assignment)
- **OpenFF Interchange 0.3+** (for system building)

These packages have complex dependencies and are best installed via conda.

## Installation Methods

### Method 1: Conda Environment (Recommended)

We recommend creating a dedicated conda environment with all dependencies:

```bash
# Create environment with OpenMM and OpenFF
conda create -n polyzymd-env python=3.11 \
    openmm openff-toolkit openff-interchange \
    -c conda-forge

# Activate the environment
conda activate polyzymd-env

# Clone and install PolyzyMD
git clone https://github.com/joelaforet/polyzymd.git
cd polyzymd
pip install -e .
```

### Method 2: Existing Environment

If you already have a conda environment with OpenMM/OpenFF (e.g., `polymerist-env`):

```bash
conda activate polymerist-env
git clone https://github.com/joelaforet/polyzymd.git
cd polyzymd
pip install -e .
```

### Method 3: pip Only (Advanced)

```{warning}
Installing OpenMM via pip can be problematic. Conda is strongly recommended.
```

```bash
pip install openmm openff-toolkit openff-interchange
pip install git+https://github.com/joelaforet/polyzymd.git
```

## HPC Installation

On HPC systems, you typically need to:

1. **Load the anaconda module**:
   ```bash
   module load anaconda
   # or: module load miniconda
   ```

2. **Clone to your projects directory** (long-term storage):
   ```bash
   cd /projects/$USER
   git clone https://github.com/joelaforet/polyzymd.git
   cd polyzymd
   ```

3. **Install in your environment**:
   ```bash
   conda activate polymerist-env
   pip install -e .
   ```

4. **Verify installation**:
   ```bash
   polyzymd info
   ```

## Verifying Installation

After installation, verify everything is working:

```bash
# Check CLI is available
polyzymd --help

# Check version and dependencies
polyzymd info
```

Expected output:
```
PolyzyMD - Molecular Dynamics for Enzyme-Polymer Systems
Version: 0.1.0

Dependencies:
  OpenMM: 8.1.1
  OpenFF Toolkit: 0.14.4
  OpenFF Interchange: 0.3.18
  Pydantic: 2.5.0

Example configs: polyzymd/configs/examples/
```

## Common Installation Issues

### "README.md not found"

If you see this error during `pip install -e .`:

```
OSError: Readme file does not exist: README.md
```

Make sure you have the latest version of the repository:

```bash
git pull origin main
```

### Missing OpenMM

If `polyzymd info` shows "OpenMM: NOT INSTALLED":

```bash
conda install -c conda-forge openmm
```

### Module Import Errors

If you get import errors when running `polyzymd`:

```bash
# Make sure you're in the right environment
conda activate polyzymd-env

# Reinstall in development mode
pip install -e .
```

## Development Installation

For contributing to PolyzyMD:

```bash
# Clone the repository
git clone https://github.com/joelaforet/polyzymd.git
cd polyzymd

# Install with development dependencies
pip install -e ".[dev]"

# Install pre-commit hooks (optional)
pre-commit install
```

## Next Steps

Once installed, continue to the {doc}`quickstart` guide to run your first simulation.
