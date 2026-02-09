# GROMACS Export and Simulation

This guide explains how to run PolyzyMD simulations using GROMACS instead of the default OpenMM backend.

## Overview

PolyzyMD can export your simulation system to GROMACS format, allowing you to:

- Run simulations on HPC clusters where GROMACS is preferred
- Use GROMACS-specific analysis tools
- Leverage GROMACS performance optimizations
- Integrate with existing GROMACS workflows

The GROMACS export maintains **1:1 parity** with OpenMM simulations using OpenFF force field defaults.

---

## When to Use GROMACS vs OpenMM

| Feature | OpenMM | GROMACS |
|---------|--------|---------|
| **GPU acceleration** | Excellent | Excellent |
| **HPC availability** | Less common | Very common |
| **Analysis tools** | MDTraj, MDAnalysis | Built-in + MDAnalysis |
| **Restart/checkpoint** | Native in PolyzyMD | Standard .cpt files |
| **Daisy-chain jobs** | Built-in support | Manual or script-based |
| **Learning curve** | Simpler | More complex |

**Use GROMACS when:**
- Your HPC cluster has optimized GROMACS installations
- You need GROMACS-specific analysis tools
- You want to integrate with existing GROMACS workflows
- You prefer GROMACS file formats (.gro, .xtc, .edr)

**Use OpenMM when:**
- You want the simplest workflow
- You need built-in daisy-chain job management
- You're running locally or on a workstation

---

## Quick Start

### Export and Run with GROMACS

```bash
# Full workflow: build system + run GROMACS simulation
polyzymd run -c config.yaml -r 1 --gromacs

# Export files only (for manual execution or HPC submission)
polyzymd run -c config.yaml -r 1 --gromacs --dry-run

# Build only (export GROMACS files, no simulation)
polyzymd build -c config.yaml -r 1 --gromacs
```

### Using a Custom GROMACS Installation

```bash
# Specify custom GROMACS path
polyzymd run -c config.yaml --gromacs --gmx-path /usr/local/gromacs/bin/gmx

# Or with module system
module load gromacs/2023.3
polyzymd run -c config.yaml --gromacs --gmx-path $(which gmx)
```

---

## Complete Example: Multi-Component System

Here's a complete workflow for running a simulation with enzyme, dynamically generated polymers, organic co-solvent, water, and substrate using GROMACS.

### Step 1: Create Configuration

```yaml
# config_gromacs.yaml
name: "LipA_polymer_gromacs"
description: "Lipase A with SBMA-EGPMA polymers in DMS/water - GROMACS"

# Enzyme
enzyme:
  name: "LipA"
  pdb_path: "structures/enzyme.pdb"

# Substrate (ligand)
substrate:
  name: "ResorufinButyrate"
  sdf_path: "structures/substrate.sdf"
  charge_method: "nagl"
  residue_name: "RBY"

# Dynamic Polymer Generation
polymers:
  enabled: true
  generation_mode: "dynamic"
  type_prefix: "SBMA-EGPMA"
  
  monomers:
    - label: "A"
      probability: 0.7
      name: "SBMA"
      smiles: "[H]C([H])=C(C(=O)OC([H])([H])C([H])([H])[N+](C([H])([H])[H])(C([H])([H])[H])C([H])([H])C([H])([H])C([H])([H])S(=O)(=O)[O-])C([H])([H])[H]"
      residue_name: "SBM"
    - label: "B"
      probability: 0.3
      name: "EGPMA"
      smiles: "[H]C([H])=C(C(=O)OC([H])([H])C([H])([H])Oc1c([H])c([H])c([H])c([H])c1[H])C([H])([H])[H]"
      residue_name: "EGM"
  
  length: 5
  count: 2
  charger: "nagl"

# Solvent: DMS/water mixture
solvent:
  primary:
    type: "water"
    model: "tip3p"
  co_solvents:
    - name: "dmso"
      volume_fraction: 0.30
      residue_name: "DMS"
  ions:
    neutralize: true
    nacl_concentration: 0.0
  box:
    padding: 1.2
    shape: "rhombic_dodecahedron"
    target_density: 1.05
    tolerance: 2.0

# Thermodynamics
thermodynamics:
  temperature: 300.0
  pressure: 1.0

# Simulation phases
simulation_phases:
  equilibration:
    ensemble: "NVT"
    duration: 1.0
    samples: 100
    time_step: 2.0
    thermostat: "LangevinMiddle"
    thermostat_timescale: 1.0
  production:
    ensemble: "NPT"
    duration: 100.0
    samples: 2500
    time_step: 2.0
    thermostat: "LangevinMiddle"
    thermostat_timescale: 1.0
    barostat: "MC"
    barostat_frequency: 25
  segments: 10

# Output
output:
  projects_directory: "."
  naming_template: "{enzyme}_{polymer_type}_gromacs_run{replicate}"
```

### Step 2: Validate Configuration

```bash
polyzymd validate -c config_gromacs.yaml
```

### Step 3: Run with GROMACS

```bash
# Full workflow (build + run)
polyzymd run -c config_gromacs.yaml -r 1 --gromacs

# Or dry-run to just export files
polyzymd run -c config_gromacs.yaml -r 1 --gromacs --dry-run
```

---

## GROMACS Output Files

When using `--gromacs`, files are created in `{projects_dir}/replicate_{N}/gromacs/`:

```
gromacs/
├── {system}.gro              # Initial coordinates
├── {system}.top              # Topology (includes all molecule types)
│
├── em.mdp                    # Energy minimization parameters
├── eq_01_nvt.mdp             # NVT equilibration (stage 1)
├── eq_02_npt.mdp             # NPT equilibration (stage 2, if applicable)
├── prod.mdp                  # Production MD parameters
│
├── posre_Protein.itp         # Position restraints for protein
├── posre_*.itp               # Position restraints for other molecules
│
├── run_{system}_gromacs.sh   # Generated shell script to run everything
│
├── em.tpr                    # Energy minimization run input
├── em.gro                    # Minimized coordinates
├── em.edr                    # Energy file
├── em.log                    # Log file
│
├── eq_01.tpr, eq_01.gro, ... # Equilibration outputs
├── eq_02.tpr, eq_02.gro, ... # (if multiple eq stages)
│
├── prod.tpr                  # Production run input
├── prod.xtc                  # Production trajectory
├── prod.edr                  # Production energies
├── prod.gro                  # Final coordinates
├── prod.cpt                  # Checkpoint for restart
│
├── prod_nojump.xtc           # Trajectory with PBC jumps removed
└── prod_centered.xtc         # Centered trajectory for visualization
```

---

## MDP Parameters

PolyzyMD generates MDP files that match OpenFF/OpenMM defaults for consistency:

### Key Parameters

| Parameter | Value | Notes |
|-----------|-------|-------|
| `cutoff-scheme` | Verlet | Modern GROMACS default |
| `rcoulomb` | 0.9 nm | OpenFF default |
| `rvdw` | 0.9 nm | OpenFF default |
| `coulombtype` | PME | Particle Mesh Ewald |
| `vdwtype` | Cut-off | With potential-shift |
| `constraints` | h-bonds | LINCS for hydrogen bonds |
| `dt` | 0.002 ps | 2 fs timestep |

### Thermostat/Barostat

The MDP parameters are generated based on your config:

```yaml
simulation_phases:
  production:
    thermostat: "LangevinMiddle"    # Maps to sd integrator
    thermostat_timescale: 1.0       # tau_t = 1.0 ps
    barostat: "MC"                  # Maps to Parrinello-Rahman
    barostat_frequency: 25          # nstpcouple
```

---

## Generated Run Script

The `run_{system}_gromacs.sh` script automates the entire workflow:

```bash
#!/bin/bash
# Generated by PolyzyMD
# Run GROMACS simulation for: LipA_SBMA-EGPMA_gromacs_run1

GMX="${GMX:-gmx}"
set -e  # Exit on error

# Energy minimization
$GMX grompp -f em.mdp -c LipA_SBMA-EGPMA_gromacs_run1.gro -p LipA_SBMA-EGPMA_gromacs_run1.top -o em.tpr
$GMX mdrun -deffnm em -v

# Equilibration stage 1 (NVT)
$GMX grompp -f eq_01_nvt.mdp -c em.gro -p LipA_SBMA-EGPMA_gromacs_run1.top -r em.gro -o eq_01.tpr
$GMX mdrun -deffnm eq_01 -v

# Equilibration stage 2 (NPT) - if applicable
$GMX grompp -f eq_02_npt.mdp -c eq_01.gro -p LipA_SBMA-EGPMA_gromacs_run1.top -r em.gro -o eq_02.tpr
$GMX mdrun -deffnm eq_02 -v

# Production
$GMX grompp -f prod.mdp -c eq_02.gro -p LipA_SBMA-EGPMA_gromacs_run1.top -o prod.tpr
$GMX mdrun -deffnm prod -v

# Post-processing
echo "System" | $GMX trjconv -s prod.tpr -f prod.xtc -o prod_nojump.xtc -pbc nojump
echo "System" | $GMX trjconv -s prod.tpr -f prod_nojump.xtc -o prod_centered.xtc -center -pbc mol

echo "Simulation complete!"
```

You can customize this script or use it as a template for HPC submission.

---

## Running on HPC Clusters

### Manual SLURM Submission

After exporting with `--dry-run`, create a SLURM script:

```bash
#!/bin/bash
#SBATCH --job-name=LipA_gmx
#SBATCH --partition=aa100
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --gres=gpu:1
#SBATCH --time=24:00:00
#SBATCH --output=slurm_%j.out

module load gromacs/2023.3

cd /path/to/replicate_1/gromacs

# Set GMX environment variable for the run script
export GMX=$(which gmx)

# Run the generated script
bash run_LipA_SBMA-EGPMA_gromacs_run1.sh
```

### Using GPU Acceleration

GROMACS GPU acceleration is configured via environment variables:

```bash
# Use GPU for non-bonded and PME
export GMX_GPU_DD_COMMS=1
export GMX_GPU_PME_PP_COMMS=1
export GMX_FORCE_UPDATE_DEFAULT_GPU=1

$GMX mdrun -deffnm prod -v -nb gpu -pme gpu -bonded gpu
```

---

## Troubleshooting

### "GROMACS executable not found"

**Cause**: `gmx` command not in PATH.

**Solutions**:
```bash
# Option 1: Load module
module load gromacs

# Option 2: Specify path explicitly
polyzymd run -c config.yaml --gromacs --gmx-path /opt/gromacs/bin/gmx

# Option 3: Add to PATH
export PATH=$PATH:/opt/gromacs/bin
```

### "grompp warnings about charge groups"

**Cause**: OpenFF generates systems without charge groups.

**Solution**: This is expected and safe to ignore. GROMACS works fine with Verlet cutoff scheme.

### "Fatal error: Number of atoms in [file] does not match"

**Cause**: Topology/coordinate mismatch, often from interrupted build.

**Solution**: Rebuild from scratch:
```bash
rm -rf replicate_*/gromacs/
polyzymd run -c config.yaml --gromacs
```

### "Simulation crashes during equilibration"

**Cause**: Initial structure has clashes or bad contacts.

**Solutions**:
1. Use longer/more aggressive energy minimization
2. Increase box padding
3. Check your input structures for issues

### "Trajectory has broken molecules"

**Cause**: PBC artifacts in raw trajectory.

**Solution**: Use the post-processed trajectories:
- `prod_nojump.xtc` - Molecules don't jump across boundaries
- `prod_centered.xtc` - System centered, good for visualization

---

## Post-Processing and Analysis

### Extract Energy Data

```bash
# Extract temperature
echo "Temperature" | gmx energy -f prod.edr -o temperature.xvg

# Extract pressure
echo "Pressure" | gmx energy -f prod.edr -o pressure.xvg

# Extract potential energy
echo "Potential" | gmx energy -f prod.edr -o potential.xvg
```

### Trajectory Analysis

```bash
# RMSD of protein
echo "Backbone Backbone" | gmx rms -s prod.tpr -f prod_centered.xtc -o rmsd.xvg

# RMSF per residue
echo "Backbone" | gmx rmsf -s prod.tpr -f prod_centered.xtc -o rmsf.xvg -res

# Radius of gyration
echo "Protein" | gmx gyrate -s prod.tpr -f prod_centered.xtc -o gyrate.xvg
```

### Visualize with VMD

```bash
vmd replicate_1/gromacs/prod.gro replicate_1/gromacs/prod_centered.xtc
```

---

## Comparison with OpenMM Output

| Metric | OpenMM | GROMACS | Notes |
|--------|--------|---------|-------|
| Coordinates | `.pdb` | `.gro` | Both readable by MDAnalysis |
| Trajectory | `.dcd` | `.xtc` | XTC is compressed |
| Energies | state_data.csv | `.edr` | Use `gmx energy` to extract |
| Checkpoint | `.xml` | `.cpt` | For restarts |
| Topology | `.xml` | `.top` | Different formats |

---

## See Also

- {doc}`cli_reference` - Complete CLI documentation including GROMACS options
- {doc}`dynamic_polymers` - Dynamic polymer generation
- {doc}`configuration` - Configuration file reference
- {doc}`quickstart` - Getting started guide
