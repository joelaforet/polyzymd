<p align="center">
  <img src="docs/source/_static/logo.png" alt="PolyzyMD Logo" width="400">
</p>

<h1 align="center">PolyzyMD</h1>

<p align="center">
  <a href="https://github.com/joelaforet/polyzymd/actions/workflows/ci.yml"><img src="https://github.com/joelaforet/polyzymd/actions/workflows/ci.yml/badge.svg" alt="CI"></a>
  <a href="https://www.python.org/downloads/"><img src="https://img.shields.io/badge/python-3.10%2B-blue.svg" alt="Python 3.10+"></a>
  <a href="https://opensource.org/licenses/MIT"><img src="https://img.shields.io/badge/License-MIT-yellow.svg" alt="License: MIT"></a>
</p>

<p align="center">
  <strong>Molecular dynamics simulation toolkit for enzyme-polymer systems.</strong>
</p>

---

## Overview

PolyzyMD provides a streamlined workflow for setting up and running MD simulations of enzymes with co-polymer chains. It handles:

- **System Building**: Combine enzyme structures, docked substrates, and random co-polymers
- **Solvation**: Add water, ions, and optional co-solvents with PACKMOL
- **Simulation**: Run equilibration and production with OpenMM
- **HPC Integration**: Daisy-chain job submission for SLURM clusters
- **Configuration**: YAML-based configuration with validation

## Documentation (web-hosted)

You can find documentation for this package [here!](https://polyzymd.readthedocs.io/en/latest/)

## Installation

### Recommended: Conda Environment

PolyzyMD depends on packages that require conda (OpenMM, OpenFF stack, PACKMOL). We recommend creating a dedicated environment:

```bash
# Clone the repository
git clone https://github.com/joelaforet/polyzymd.git
cd polyzymd

# Create environment from our tested configuration
conda env create -f devtools/conda-envs/test-env.yml
conda activate test-env

# Install polyzymd
pip install -e .
```

### Core Dependencies

These are installed automatically via the conda environment:

| Package | Source | Purpose |
|---------|--------|---------|
| OpenMM >= 8.0 | conda | MD engine |
| OpenFF Toolkit >= 0.16 | conda | Force field parameterization |
| OpenFF Interchange >= 0.4 | conda | System interchange |
| OpenFF NAGL >= 0.3 | conda | ML-based partial charges |
| Polymerist >= 1.0 | pip | Polymer generation |
| PACKMOL | conda | System packing |
| MDTraj | conda | Trajectory analysis |

### Optional Dependencies

For AM1-BCC charge assignment (alternative to NAGL):

```bash
# Option 1: AmberTools (free, open source)
# Note: Only available for Linux; conflicts may require a separate environment
conda install -c conda-forge ambertools

# Option 2: OpenEye Toolkit (commercial, faster)
conda install -c openeye openeye-toolkits
```

## Quick Start

### 1. Initialize a Project

```bash
polyzymd init --name my_simulation
cd my_simulation
```

This creates a project directory with a template `config.yaml` and placeholder files.

### 2. Add Your Structure Files

```bash
# Copy your prepared enzyme PDB
cp /path/to/enzyme.pdb structures/

# Copy your docked substrate SDF (optional)
cp /path/to/substrate.sdf structures/

# Remove placeholder files
rm structures/*.placeholder.txt
```

### 3. Edit the Configuration

Open `config.yaml` and uncomment/customize the sections you need. The template
includes all options with example values - uncomment what you need and update
the values for your system.

Example minimal configuration:

```yaml
name: "LipA_substrate_polymer"

enzyme:
  name: "LipA"
  pdb_path: "structures/enzyme.pdb"

substrate:
  name: "MySubstrate"
  sdf_path: "structures/substrate.sdf"
  charge_method: "nagl"

polymers:
  enabled: true
  type_prefix: "SBMA-EGPMA"
  monomers:
    - label: "A"
      probability: 0.98
      name: "SBMA"
    - label: "B"
      probability: 0.02
      name: "EGPMA"
  length: 5
  count: 2

solvent:
  primary:
    type: "water"
    model: "tip3p"
  ions:
    neutralize: true
    nacl_concentration: 0.15

thermodynamics:
  temperature: 300.0
  pressure: 1.0

simulation_phases:
  equilibration:
    ensemble: "NVT"
    duration: 1.0
  production:
    ensemble: "NPT"
    duration: 100.0
  segments: 10

output:
  projects_directory: "."
  scratch_directory: null
  naming_template: "{enzyme}_{polymer_type}_{temperature}K_run{replicate}"
```

### 4. Validate Configuration

```bash
polyzymd validate -c config.yaml
```

### 5. Submit to HPC Cluster

```bash
# Dry run (generate scripts without submitting)
polyzymd submit -c config.yaml --replicates 1-5 --preset aa100 --dry-run

# Submit for real
polyzymd submit -c config.yaml --replicates 1-5 --preset aa100 --email you@university.edu
```

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

## Directory Structure

PolyzyMD supports separate directories for long-term storage (projects) and high-performance scratch:

```
/projects/user/polyzymd/          # Scripts, logs, configs
├── config.yaml
├── job_scripts/
└── slurm_logs/

/scratch/user/simulations/        # Trajectories, checkpoints
├── LipA_SBMA-EGPMA_300K_run1/
│   ├── equilibration/
│   └── production_seg0/
└── ...
```

## SLURM Presets

| Preset | Partition | Description |
|--------|-----------|-------------|
| `aa100` | aa100 | NVIDIA A100 GPUs (24h limit) |
| `al40` | al40 | NVIDIA L40 GPUs (24h limit) |
| `blanca-shirts` | blanca-shirts | Shirts lab partition (7d limit) |
| `testing` | atesting | Quick tests (1h limit) |

Override time limits with `--time-limit`:

```bash
# Quick 2-minute test
polyzymd submit -c config.yaml --preset testing --time-limit 0:02:00
```

## Documentation (for developers)

Documentation is built with Sphinx and located in the `docs/` directory.

### Building the Docs

```bash
cd docs
make html
```

The built documentation will be in `docs/build/html/`. Open `docs/build/html/index.html` in a browser to view.

### Live Rebuild (Development)

For active documentation work, use auto-rebuild which watches for changes:

```bash
cd docs
make livehtml
```

This starts a local server (typically at `http://127.0.0.1:8000`) that automatically rebuilds when you edit source files.

### Clean Rebuild

If you encounter stale content or build errors:

```bash
cd docs
make clean && make html
```

## License

MIT License - see LICENSE file for details.

## Citation

If you use PolyzyMD in your research, please cite:

```
@software{polyzymd,
  author = {Laforet Jr., Joseph R.},
  title = {PolyzyMD: MD Simulations for Enzyme-Polymer Systems},
  year = {2026},
  url = {https://github.com/joelaforet/polyzymd}
}
```
