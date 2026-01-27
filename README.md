# PolyzyMD

Molecular dynamics simulation toolkit for enzyme-polymer systems.

## Overview

PolyzyMD provides a streamlined workflow for setting up and running MD simulations of enzymes with co-polymer chains. It handles:

- **System Building**: Combine enzyme structures, docked substrates, and random co-polymers
- **Solvation**: Add water, ions, and optional co-solvents with PACKMOL
- **Simulation**: Run equilibration and production with OpenMM
- **HPC Integration**: Daisy-chain job submission for SLURM clusters
- **Configuration**: YAML-based configuration with validation

## Installation

### From Source (Development)

```bash
git clone https://github.com/joelaforet/polyzymd.git
cd polyzymd
pip install -e .
```

### Dependencies

PolyzyMD requires:
- Python >= 3.10
- OpenMM >= 8.0
- OpenFF Toolkit >= 0.14.0
- OpenFF Interchange >= 0.3.0
- Polymerist >= 1.0.0

We recommend using a conda environment with these packages pre-installed.

## Quick Start

### 1. Create a Configuration File

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

### 2. Validate Configuration

```bash
polyzymd validate -c config.yaml
```

### 3. Submit to HPC Cluster

```bash
# Dry run (generate scripts without submitting)
polyzymd submit -c config.yaml --replicates 1-5 --preset aa100 --dry-run

# Submit for real
polyzymd submit -c config.yaml --replicates 1-5 --preset aa100 --email you@university.edu
```

## CLI Commands

| Command | Description |
|---------|-------------|
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

## License

MIT License - see LICENSE file for details.

## Citation

If you use PolyzyMD in your research, please cite:

```
@software{polyzymd,
  author = {Laforet, Joe},
  title = {PolyzyMD: MD Simulations for Enzyme-Polymer Systems},
  year = {2024},
  url = {https://github.com/joelaforet/polyzymd}
}
```
