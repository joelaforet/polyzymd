# Quick Start Guide

This guide walks through running your first PolyzyMD simulation.

## Overview

A typical PolyzyMD workflow:

1. **Prepare input files** (PDB, SDF)
2. **Create configuration** (YAML file)
3. **Validate** the configuration
4. **Submit** jobs to HPC cluster

## Step 1: Prepare Your Input Files

Create a project directory with your input structures:

```bash
mkdir my_simulation
cd my_simulation
mkdir structures
```

### Required Files

| File | Description | Format |
|------|-------------|--------|
| `structures/enzyme.pdb` | Prepared enzyme structure | PDB |
| `structures/substrate.sdf` | Docked substrate (optional) | SDF |

### Enzyme PDB Requirements

Your enzyme PDB must be **simulation-ready**:

- All hydrogens added
- Missing residues/atoms modeled  
- Proper protonation states
- No alternate conformations
- No crystallographic waters/ligands (unless intentional)

**Recommended tools:**
- [PDB2PQR](https://server.poissonboltzmann.org/) - Protonation
- [CHARMM-GUI](https://www.charmm-gui.org/) - Full preparation
- PyMOL/Chimera - Manual inspection

## Step 2: Create Configuration File

Create `config.yaml`:

```yaml
name: "MyEnzyme_simulation"
description: "My first PolyzyMD simulation"

# Enzyme configuration
enzyme:
  name: "MyEnzyme"
  pdb_path: "structures/enzyme.pdb"

# Substrate (optional - set to null for apo simulation)
substrate: null

# Polymers (optional - set to null for enzyme-only)
polymers: null

# Solvent
solvent:
  primary:
    type: "water"
    model: "tip3p"
  ions:
    neutralize: true
    nacl_concentration: 0.15  # 150 mM

# Thermodynamics
thermodynamics:
  temperature: 300.0  # Kelvin
  pressure: 1.0       # atm

# Simulation settings
simulation_phases:
  equilibration:
    ensemble: "NVT"
    duration: 1.0     # ns
  production:
    ensemble: "NPT"
    duration: 10.0    # ns (short for testing)
  segments: 1         # Single segment for testing

# Output
output:
  projects_directory: "."
  scratch_directory: null
  naming_template: "{enzyme}_{temperature}K_run{replicate}"
```

## Step 3: Validate Configuration

Before running, validate your configuration:

```bash
polyzymd validate -c config.yaml
```

Expected output:
```
Validating configuration: config.yaml
Configuration is valid!

Summary:
  Name: MyEnzyme_simulation
  Enzyme: MyEnzyme
  Substrate: None (apo simulation)
  Polymers: Disabled
  Temperature: 300.0 K
  Pressure: 1.0 atm

Simulation phases:
  Equilibration: 1.0 ns (NVT)
  Production: 10.0 ns (NPT)
  Segments: 1
```

## Step 4: Test Build (Dry Run)

Test that the system can be built:

```bash
polyzymd build -c config.yaml --dry-run
```

## Step 5: Run Locally (Optional)

For testing on a local machine with GPU:

```bash
polyzymd run -c config.yaml --replicate 1
```

## Step 6: Submit to HPC

For production runs on an HPC cluster:

```bash
# Generate scripts without submitting (dry run)
polyzymd submit -c config.yaml \
    --replicates 1-3 \
    --preset aa100 \
    --dry-run

# Submit for real
polyzymd submit -c config.yaml \
    --replicates 1-3 \
    --preset aa100 \
    --email your.email@university.edu
```

## Step 7: Monitor Jobs

```bash
# Check job status
squeue -u $USER

# View job output
tail -f slurm_logs/*.out
```

## Output Files

After completion:

```
my_simulation/
├── config.yaml
├── job_scripts/           # SLURM scripts
├── slurm_logs/            # Job output
└── MyEnzyme_300K_run1/    # Simulation output
    ├── system.pdb
    ├── equilibration/
    │   └── trajectory.dcd
    └── production_seg0/
        ├── trajectory.dcd
        ├── checkpoint.chk
        └── state_data.csv
```

## Next Steps

- Add a substrate: See {doc}`configuration`
- Add polymers: See {doc}`polymers`
- Add restraints: See {doc}`restraints`
- Run on HPC: See {doc}`hpc_slurm`
