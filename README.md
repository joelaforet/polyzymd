<p align="center">
  <img src="docs/source/_static/logo.png" alt="PolyzyMD Logo" width="400">
</p>

<h1 align="center">PolyzyMD</h1>

<p align="center">
  <a href="https://pypi.org/project/polyzymd/"><img src="https://img.shields.io/pypi/v/polyzymd.svg" alt="PyPI"></a>
  <a href="https://github.com/joelaforet/polyzymd/actions/workflows/ci.yml"><img src="https://github.com/joelaforet/polyzymd/actions/workflows/ci.yml/badge.svg" alt="CI"></a>
  <a href="https://www.python.org/downloads/"><img src="https://img.shields.io/badge/python-3.10%2B-blue.svg" alt="Python 3.10+"></a>
  <a href="https://opensource.org/licenses/MIT"><img src="https://img.shields.io/badge/License-MIT-yellow.svg" alt="License: MIT"></a>
</p>

<p align="center">
  <strong>Molecular dynamics simulation toolkit for enzyme-polymer systems.</strong>
</p>

<p align="center">
  <a href="https://polyzymd.readthedocs.io">Documentation</a> •
  <a href="#installation">Installation</a> •
  <a href="#quick-start">Quick Start</a>
</p>

---

## Overview

PolyzyMD provides a streamlined workflow for setting up and running MD simulations of enzymes with co-polymer chains. It handles:

- **System Building**: Combine enzyme structures, docked substrates, and random co-polymers
- **Solvation**: Add water, ions, and optional co-solvents with PACKMOL
- **Simulation**: Run equilibration and production with OpenMM
- **HPC Integration**: Daisy-chain job submission for SLURM clusters
- **Configuration**: YAML-based configuration with validation

## Installation

PolyzyMD requires OpenMM and the OpenFF stack, which are only available via conda/mamba.

### User Installation

**Step 1: Create conda environment**

Using mamba (recommended for speed):

```bash
mamba create -n polyzymd-env python=3.11 openmm openff-toolkit openff-interchange \
    openff-nagl openff-nagl-models openff-forcefields openff-units packmol \
    mbuild mdtraj numpy scipy pandas pydantic pyyaml click tqdm -c conda-forge
mamba activate polyzymd-env
```

Or using our environment file:

```bash
mamba env create -f https://raw.githubusercontent.com/joelaforet/polyzymd/main/devtools/conda-envs/polyzymd-env.yml
mamba activate polyzymd-env
```

**Step 2: Install PolyzyMD**

```bash
pip install polyzymd
```

### Developer Installation

For contributing or modifying the source code:

```bash
git clone https://github.com/joelaforet/polyzymd.git
cd polyzymd
mamba env create -f devtools/conda-envs/polyzymd-env.yml
mamba activate polyzymd-env
pip install -e ".[dev]"
```

> **Note for uv users:** OpenMM and the OpenFF stack are not available on PyPI and must be installed via conda/mamba. You can use uv for other packages within an activated conda environment, but the core simulation dependencies require conda.

For detailed installation instructions, including HPC setup and troubleshooting, see the [Installation Guide](https://polyzymd.readthedocs.io/en/latest/tutorials/installation.html).

## Quick Start

### 1. Initialize a Project

```bash
polyzymd init --name my_simulation
cd my_simulation
```

This creates a project directory with a template `config.yaml` and placeholder files.

### 2. Add Your Structure Files

```bash
cp /path/to/enzyme.pdb structures/
cp /path/to/substrate.sdf structures/  # optional
```

### 3. Edit Configuration & Run

```bash
# Edit config.yaml with your settings, then:
polyzymd validate -c config.yaml
polyzymd submit -c config.yaml --replicates 1-5 --preset aa100
```

See the [Quick Start Guide](https://polyzymd.readthedocs.io/en/latest/tutorials/quickstart.html) for a complete walkthrough.

## CLI Commands

| Command | Description |
|---------|-------------|
| `polyzymd init -n my_project` | Initialize a new project directory |
| `polyzymd validate -c config.yaml` | Validate configuration file |
| `polyzymd build -c config.yaml` | Build simulation system |
| `polyzymd run -c config.yaml` | Run simulation locally |
| `polyzymd submit -c config.yaml` | Submit daisy-chain jobs to SLURM |
| `polyzymd continue -w workdir -s 2` | Continue from previous segment |
| `polyzymd info` | Show installation information |

## Documentation

Full documentation is available at **[polyzymd.readthedocs.io](https://polyzymd.readthedocs.io)**.

- [Installation Guide](https://polyzymd.readthedocs.io/en/latest/tutorials/installation.html)
- [Configuration Reference](https://polyzymd.readthedocs.io/en/latest/tutorials/configuration.html)
- [HPC & SLURM Guide](https://polyzymd.readthedocs.io/en/latest/tutorials/hpc_slurm.html)
- [API Reference](https://polyzymd.readthedocs.io/en/latest/api/overview.html)

## License

MIT License - see LICENSE file for details.

## Citation

If you use PolyzyMD in your research, please cite:

```bibtex
@software{polyzymd,
  author = {Laforet Jr., Joseph R.},
  title = {PolyzyMD: MD Simulations for Enzyme-Polymer Systems},
  year = {2026},
  url = {https://github.com/joelaforet/polyzymd}
}
```
