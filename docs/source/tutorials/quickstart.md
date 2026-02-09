# Quick Start Guide

This guide walks through running your first PolyzyMD simulation.

## Overview

A typical PolyzyMD workflow:

1. **Initialize project** with `polyzymd init`
2. **Add structure files** (PDB, SDF)
3. **Configure simulation** (edit YAML file)
4. **Validate** the configuration
5. **Submit** jobs to HPC cluster

## Step 1: Initialize Your Project

The easiest way to start is with the `init` command, which creates a complete
project scaffold:

```bash
polyzymd init --name my_simulation
```

This creates:

```
my_simulation/
├── config.yaml              <- Template configuration (edit this)
├── structures/              <- Add your PDB/SDF files here
│   ├── place_protein_here.placeholder.txt
│   └── place_ligand_here.placeholder.txt
├── job_scripts/             <- SLURM scripts will go here
└── slurm_logs/              <- Job output will go here
```

```{tip}
The `config.yaml` template has all sections commented out with example values.
This teaches you the configuration syntax while you customize it for your system.
```

## Step 2: Add Your Structure Files

Copy your prepared structure files into the `structures/` directory:

```bash
cd my_simulation

# Copy your enzyme PDB
cp /path/to/my_enzyme.pdb structures/enzyme.pdb

# Copy your substrate SDF (if using)
cp /path/to/my_substrate.sdf structures/substrate.sdf

# Remove the placeholder files
rm structures/*.placeholder.txt
```

### Enzyme PDB Requirements

Your enzyme PDB must be **simulation-ready**:

- All hydrogens added
- Missing residues/atoms modeled  
- Proper protonation states
- No alternate conformations
- No crystallographic waters/ligands (unless intentional)

**Recommended preparation tools:**
- [PDB2PQR](https://server.poissonboltzmann.org/) - Protonation
- [CHARMM-GUI](https://www.charmm-gui.org/) - Full preparation
- PyMOL/Chimera - Manual inspection

### Substrate SDF Requirements (if using)

- 3D coordinates (from docking or crystal structure)
- Correct protonation state for simulation pH
- Single conformer (or specify `conformer_index` in config)

## Step 3: Configure Your Simulation

Open `config.yaml` in your editor. The template has all sections commented out
with placeholder values. You need to:

1. **Uncomment** the sections you need
2. **Replace** placeholder values with your actual data
3. **Leave commented** sections you don't need

### Minimal Configuration (Enzyme Only)

For a simple enzyme-in-water simulation, uncomment and edit these sections:

```yaml
name: "my_first_simulation"
description: "Testing PolyzyMD with my enzyme"

# Enzyme - REQUIRED
enzyme:
  name: "MyEnzyme"
  pdb_path: "structures/enzyme.pdb"

# Solvent - REQUIRED
solvent:
  primary:
    type: "water"
    model: "tip3p"
  co_solvents: []
  ions:
    neutralize: true
    nacl_concentration: 0.15
  box:
    padding: 1.2
    shape: "rhombic_dodecahedron"
    target_density: 1.0
    tolerance: 2.0

# Thermodynamics - REQUIRED
thermodynamics:
  temperature: 300.0
  pressure: 1.0

# Simulation phases - REQUIRED
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
    duration: 10.0       # Short for testing
    samples: 250
    time_step: 2.0
    thermostat: "LangevinMiddle"
    thermostat_timescale: 1.0
    barostat: "MC"
    barostat_frequency: 25
  segments: 1            # Single segment for testing

# Output - REQUIRED
output:
  projects_directory: "."
  scratch_directory: null
  job_scripts_subdir: "job_scripts"
  slurm_logs_subdir: "slurm_logs"
  naming_template: "{enzyme}_{temperature}K_run{replicate}"
  save_checkpoint: true
  save_state_data: true
  trajectory_format: "dcd"
```

### Adding a Substrate

To include a docked ligand, uncomment the substrate section:

```yaml
substrate:
  name: "MyLigand"
  sdf_path: "structures/substrate.sdf"
  conformer_index: 0
  charge_method: "nagl"
  residue_name: "LIG"
```

### Adding Polymers

To add co-polymer chains around your enzyme:

```yaml
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
  cache_directory: ".polymer_cache"
```

```{warning}
**YAML List Syntax**

When defining monomers (or any list), each `-` starts a **new list item**. All fields for one monomer must be grouped together:

**Incorrect:**
~~~yaml
monomers:
  - label: "A"
  - probability: 0.98
  - name: "SBMA"
~~~

**Correct:**
~~~yaml
monomers:
  - label: "A"
    probability: 0.98
    name: "SBMA"
~~~
```

See {doc}`polymers` for detailed polymer configuration.

### Adding Restraints

To keep a substrate in the active site:

```yaml
restraints:
  - type: "flat_bottom"
    name: "substrate_restraint"
    atom1:
      selection: "resid 77 and name OG"
      description: "Catalytic serine"
    atom2:
      selection: "resname LIG and name C1"
      description: "Substrate carbon"
    distance: 3.5
    force_constant: 10000.0
    enabled: true
```

See {doc}`restraints` for detailed restraint configuration.

## Step 4: Validate Configuration

Before running, validate your configuration:

```bash
polyzymd validate -c config.yaml
```

This checks:
- All required sections are present
- File paths exist
- Values are within valid ranges
- No conflicting settings

Expected output:
```
Validating configuration: config.yaml
Configuration is valid!

Summary:
  Name: my_first_simulation
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

## Step 5: Test Build (Dry Run)

Test that the system can be built without actually building:

```bash
polyzymd build -c config.yaml --dry-run
```

## Step 6: Run Locally (Optional)

For testing on a local machine with GPU:

```bash
polyzymd run -c config.yaml --replicate 1
```

## Step 7: Submit to HPC

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

## Step 8: Monitor Jobs

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
    ├── solvated_system.pdb
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
