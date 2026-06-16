<p align="center">
  <img src="docs/source/_static/Higher_resolution_PolyzyMD_logo_v1.png" alt="PolyzyMD Logo" width="400">
</p>

<h1 align="center">PolyzyMD</h1>

<p align="center">
  <a href="https://github.com/joelaforet/polyzymd/actions/workflows/ci.yml"><img src="https://github.com/joelaforet/polyzymd/actions/workflows/ci.yml/badge.svg" alt="CI"></a>
  <a href="https://www.python.org/downloads/"><img src="https://img.shields.io/badge/python-3.12%2B-blue.svg" alt="Python 3.12+"></a>
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
- **HPC Integration**: Self-resubmitting job submission for SLURM clusters
- **Configuration**: YAML-based configuration with validation

## Analysis Status

- Stable comparison and plotting workflows: RMSF, contacts, distances, catalytic triad, secondary structure, RMSD, Rg, SASA, and hydrogen bonds
- Experimental analysis features have been removed from the active CLI surface until their definitions and interpretation are finalized
- Analysis supports both OpenMM (DCD) and GROMACS (XTC) trajectories via engine-aware trajectory resolution

## Installation

PolyzyMD uses [pixi](https://pixi.sh) for environment management. Pixi handles
all dependencies (OpenMM, OpenFF stack, CUDA, etc.) from conda-forge with
reproducible lockfiles.

### 1. Install pixi

```bash
curl -fsSL https://pixi.sh/install.sh | sh
```

### 2. Clone and install

```bash
git clone https://github.com/joelaforet/polyzymd.git
cd polyzymd
```

**For local use (building systems, validation, package work):**

```bash
pixi install -e build
pixi shell -e build
```

**For HPC clusters (GPU simulations):**

Pick the simulation runtime environment that matches your cluster's CUDA version:

| Cluster | CUDA | Environment |
|---------|------|-------------|
| CU Boulder Blanca | 12.4 | `sim-cuda-12-4` |
| PSC Bridges2 | 12.6 | `sim-cuda-12-6` |

```bash
# Example for Blanca:
pixi install -e sim-cuda-12-4
pixi shell -e sim-cuda-12-4

# Example for Bridges2:
pixi install -e sim-cuda-12-6
pixi shell -e sim-cuda-12-6
```

After `pixi shell`, the `polyzymd` command is on PATH and works normally.

PolyzyMD v1.3 intentionally splits environments by workflow. Use `build` to
prepare systems, use `sim-cuda-*` only on GPU nodes to execute simulations, and
use `analysis` to compare trajectories and make plots. The same project files,
prepared systems, checkpoints, and trajectories move between these environments.
This split preserves CUDA 12.4 cluster support while keeping package and
analysis work on Python 3.12 with NumPy 2.

macOS is supported for non-CUDA setup, development, and trajectory analysis
where dependencies solve. The CUDA simulation environments are Linux-only and
intended for GPU clusters. AmberTools/AM1-BCC support is Linux-focused in the
pixi workflow; macOS users should use NAGL/OpenFF charging or pre-charged
molecules and open an issue if they need AmberTools support.

### How to find your CUDA version

Run on a GPU node:

```bash
nvidia-smi | head -1
```

The driver version in the top-right maps to a maximum supported CUDA version.
Use the environment whose CUDA version does not exceed your driver.

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
```

Switch to the CUDA simulation environment that matches your cluster before
running OpenMM simulation commands such as `polyzymd submit`:

```bash
# Blanca / CUDA 12.4
pixi shell -e sim-cuda-12-4
polyzymd submit -c config.yaml --replicates 1-5 --preset blanca-shirts

# Bridges2 / CUDA 12.6
pixi shell -e sim-cuda-12-6
polyzymd submit -c config.yaml --replicates 1-5 --preset bridges2
```

The `--preset` flag selects SLURM configuration and automatically picks the
correct pixi environment for generated job scripts (`sim-cuda-12-4` for Blanca,
`sim-cuda-12-6` for Bridges2). You can override with `--pixi-env`:

```bash
polyzymd submit -c config.yaml --replicates 1-5 --preset bridges2 --pixi-env sim-cuda-12-6
```

See the [Quick Start Guide](https://polyzymd.readthedocs.io/en/latest/get_started/quickstart.html) for a complete walkthrough.

## CLI Commands

| Command | Description |
|---------|-------------|
| `polyzymd init -n my_project` | Initialize a new project directory |
| `polyzymd validate -c config.yaml` | Validate configuration file |
| `polyzymd build -c config.yaml` | Build simulation system |
| `polyzymd run -c config.yaml --engine gromacs` | Build and run GROMACS simulation |
| `polyzymd submit -c config.yaml` | Submit self-resubmitting jobs to SLURM |
| `polyzymd run-segment -c config.yaml` | Run a single production segment |
| `polyzymd check-progress -c config.yaml` | Check simulation completion status |
| `polyzymd recover -c config.yaml` | Resume a stalled simulation |
| `polyzymd info` | Show installation information |

## Pixi Environments

PolyzyMD uses [pixi](https://pixi.sh) instead of conda/mamba. Key differences:

- **No `conda activate`** — use `pixi shell -e <env>` instead
- **No environment YAML** — `pixi.toml` + `pixi.lock` are the single source of truth
- **Reproducible** — the lockfile pins every package to exact versions
- **CUDA-aware** — each environment pins the correct CUDA version

| Environment | Use case | Stack policy | Requires GPU? |
|-------------|----------|--------------|---------------|
| `build` | System building, PDB prep, validation, package work | Python 3.12 + NumPy 2 | No |
| `analysis` | Trajectory analysis, comparison, and plotting | Python 3.12 + NumPy 2 | No |
| `test` | CI/local tests and lint/docs tooling | Python 3.12 + NumPy 2 | No |
| `sim-cuda-12-4` | Simulation execution on CUDA 12.4 clusters (Blanca) | Python 3.12 + NumPy 1.x + OpenMM 8.1.x | Yes |
| `sim-cuda-12-6` | Simulation execution on CUDA 12.6 clusters (Bridges2) | Python 3.12 + current OpenMM | Yes |

AmberTools/AM1-BCC support is Linux-focused in PolyzyMD's pixi workflow and is
not part of the default v1.3 NumPy 2 build/analysis/test solve. Use NAGL/OpenFF
charging or provide pre-charged molecules; open an issue if you need a
dedicated AmberTools/AM1-BCC environment for a specific platform.

### Adding support for a new cluster

1. Determine the CUDA version (`nvidia-smi` on a GPU node)
2. Add a new `[feature.sim-cuda-X-Y]` block in `pixi.toml` following the existing pattern
3. Add the corresponding environment in `[environments]`
4. Add the preset mapping in `PRESET_DEFAULT_PIXI_ENV` in `slurm.py`
5. File a PR

## Documentation

Full documentation is available at **[polyzymd.readthedocs.io](https://polyzymd.readthedocs.io)**.

- [Quick Start Guide](https://polyzymd.readthedocs.io/en/latest/get_started/quickstart.html)
- [Installation Guide](https://polyzymd.readthedocs.io/en/latest/get_started/installation.html)
- [Configuration Reference](https://polyzymd.readthedocs.io/en/latest/reference/configuration.html)
- [HPC & SLURM Guide](https://polyzymd.readthedocs.io/en/latest/how_to/hpc_slurm.html)
- [API Reference](https://polyzymd.readthedocs.io/en/latest/api/overview.html)

## License

MIT License - see LICENSE file for details.

## Citation

If you use PolyzyMD in your research, please cite:

```bibtex
@software{polyzymd,
  author = {Laforet Jr., Joseph R.},
  title = {PolyzyMD: Polymer-Enzyme Interactions Studied with Molecular Dynamics},
  year = {2026},
  url = {https://github.com/joelaforet/polyzymd}
}
```
