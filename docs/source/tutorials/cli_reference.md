# CLI Reference

Complete reference for all PolyzyMD command-line interface commands.

## Global Options

All commands support these global options:

```bash
polyzymd --version        # Show version
polyzymd --help           # Show help
polyzymd -v <command>     # Verbose output (debug logging)
```

---

## polyzymd init

Initialize a new PolyzyMD project directory with template files.

```bash
polyzymd init --name <project_name>
polyzymd init -n <project_name>
```

### Options

| Option | Short | Required | Description |
|--------|-------|----------|-------------|
| `--name` | `-n` | Yes | Name of the project directory to create |

### What It Creates

```
<project_name>/
├── config.yaml              <- Template configuration (edit this)
├── structures/              <- Add your PDB/SDF files here
│   ├── place_protein_here.placeholder.txt
│   └── place_ligand_here.placeholder.txt
├── job_scripts/             <- Generated SLURM scripts go here
└── slurm_logs/              <- SLURM output logs go here
```

### Example

```bash
# Create a new project
polyzymd init --name lipase_dmso_study
cd lipase_dmso_study

# Add your structure files
cp ~/structures/LipA.pdb structures/enzyme.pdb
cp ~/docking/substrate.sdf structures/substrate.sdf

# Remove placeholder files
rm structures/*.placeholder.txt

# Edit the configuration
nano config.yaml

# Validate
polyzymd validate -c config.yaml
```

### Notes

- The command will fail if the directory already exists
- The template `config.yaml` has all sections commented out with example values
- Uncomment and modify only the sections you need

---

## polyzymd validate

Validate a configuration file without building or running.

```bash
polyzymd validate --config <path>
polyzymd validate -c <path>
```

### Options

| Option | Short | Required | Description |
|--------|-------|----------|-------------|
| `--config` | `-c` | Yes | Path to YAML configuration file |

### What It Checks

- YAML syntax validity
- Required fields are present
- Referenced files (PDB, SDF) exist
- Monomer probabilities sum to 1.0
- Valid enum values (water model, ensemble, etc.)
- Co-solvent specification (volume_fraction XOR concentration)

### Example

```bash
polyzymd validate -c config.yaml
```

**Output (success):**
```
Validating configuration: config.yaml
Configuration is valid!

Summary:
  Name: LipA_polymer_simulation
  Enzyme: LipA
  Substrate: ResorufinButyrate
  Polymers: SBMA-EGPMA
    Count: 2
    Length: 5
    Monomer A: 98.0%
    Monomer B: 2.0%
  Temperature: 300.0 K
  Pressure: 1.0 atm

Simulation phases:
  Equilibration: 1.0 ns (NVT)
  Production: 100.0 ns (NPT)
  Segments: 10
```

---

## polyzymd build

Build the simulation system (parameterize, solvate) without running.

```bash
polyzymd build --config <path> [options]
polyzymd build -c <path> -r <replicate>
```

### Options

| Option | Short | Required | Default | Description |
|--------|-------|----------|---------|-------------|
| `--config` | `-c` | Yes | - | Path to YAML configuration file |
| `--replicate` | `-r` | No | 1 | Replicate number (affects polymer random seed) |
| `--scratch-dir` | - | No | from config | Override scratch directory |
| `--projects-dir` | - | No | from config | Override projects directory |
| `--output-dir` | `-o` | No | from config | Alias for --scratch-dir |
| `--dry-run` | - | No | false | Validate only, don't build |

### Example

```bash
# Build replicate 1
polyzymd build -c config.yaml -r 1

# Build with custom output directory
polyzymd build -c config.yaml -r 1 --scratch-dir ./test_output

# Dry run to check configuration
polyzymd build -c config.yaml --dry-run
```

### Output Files

The build command creates:
- `solvated_system.pdb` - Complete system with water and ions
- `system.xml` - OpenMM serialized system
- `topology.json` - Topology information
- `build.log` - Build process log

---

## polyzymd run

Build and run a complete simulation (equilibration + first production segment).

```bash
polyzymd run --config <path> [options]
polyzymd run -c <path> -r <replicate>
```

### Options

| Option | Short | Required | Default | Description |
|--------|-------|----------|---------|-------------|
| `--config` | `-c` | Yes | - | Path to YAML configuration file |
| `--replicate` | `-r` | No | 1 | Replicate number |
| `--scratch-dir` | - | No | from config | Override scratch directory |
| `--projects-dir` | - | No | from config | Override projects directory |
| `--segment-time` | - | No | auto | Override production time per segment (ns) |
| `--segment-frames` | - | No | auto | Override frames per segment |
| `--skip-build` | - | No | false | Skip building (use existing system) |

### Example

```bash
# Run locally (useful for testing)
polyzymd run -c config.yaml -r 1

# Run with shorter segment for testing
polyzymd run -c config.yaml -r 1 --segment-time 0.1 --segment-frames 10
```

### Workflow

1. Load and validate configuration
2. Build system (enzyme + substrate + polymers + solvent)
3. Apply restraints (if configured)
4. Run energy minimization
5. Run equilibration (NVT)
6. Run first production segment (NPT)

---

## polyzymd submit

Submit daisy-chain jobs to SLURM for HPC execution.

```bash
polyzymd submit --config <path> --replicates <range> [options]
polyzymd submit -c <path> -r 1-5 --preset aa100
```

### Options

| Option | Short | Required | Default | Description |
|--------|-------|----------|---------|-------------|
| `--config` | `-c` | Yes | - | Path to YAML configuration file |
| `--replicates` | `-r` | No | "1" | Replicate range (e.g., "1-5", "1,3,5") |
| `--preset` | - | No | aa100 | SLURM partition preset |
| `--email` | - | No | "" | Email for job notifications |
| `--scratch-dir` | - | No | from config | Override scratch directory |
| `--projects-dir` | - | No | from config | Override projects directory |
| `--output-dir` | - | No | auto | Directory for job scripts |
| `--time-limit` | - | No | from preset | Override SLURM time limit (HH:MM:SS) |
| `--dry-run` | - | No | false | Generate scripts without submitting |

### SLURM Presets

| Preset | Partition | Time Limit | Description |
|--------|-----------|------------|-------------|
| `aa100` | aa100 | 24:00:00 | NVIDIA A100 GPUs |
| `al40` | al40 | 24:00:00 | NVIDIA L40 GPUs |
| `blanca-shirts` | blanca-shirts | 7-00:00:00 | Shirts lab partition |
| `testing` | atesting | 01:00:00 | Quick tests |

### Example

```bash
# Dry run first to check scripts
polyzymd submit -c config.yaml -r 1-5 --preset aa100 --dry-run

# Submit for real with email notifications
polyzymd submit -c config.yaml -r 1-5 --preset aa100 --email you@university.edu

# Quick test with short time limit
polyzymd submit -c config.yaml -r 1 --preset testing --time-limit 0:05:00

# Custom directories for HPC
polyzymd submit -c config.yaml -r 1-3 --preset aa100 \
    --scratch-dir /scratch/alpine/$USER/sims \
    --projects-dir /projects/$USER/polyzymd
```

### Daisy-Chain Jobs

The submit command creates a chain of dependent SLURM jobs:

```
Job 1 (initial)  →  Job 2 (continue)  →  Job 3 (continue)  →  ...
   build + eq          segment 1           segment 2
   + segment 0
```

Each job only starts after the previous one completes successfully.

---

## polyzymd continue

Continue a simulation from a previous segment checkpoint.

```bash
polyzymd continue --working-dir <path> --segment <n> --segment-time <ns>
polyzymd continue -w <path> -s <n> -t <ns>
```

### Options

| Option | Short | Required | Default | Description |
|--------|-------|----------|---------|-------------|
| `--working-dir` | `-w` | Yes | - | Working directory with previous segment |
| `--segment` | `-s` | Yes | - | Segment index to run (1-based) |
| `--segment-time` | `-t` | Yes | - | Duration of this segment (ns) |
| `--num-samples` | `-n` | No | 250 | Number of frames to save |

### Example

```bash
# Continue to segment 2 (after segment 1 completed)
polyzymd continue -w /scratch/user/sim/LipA_300K_run1 -s 2 -t 10.0 -n 250
```

### Notes

- This command is typically called by SLURM continuation scripts, not manually
- It loads the checkpoint from the previous segment automatically
- The segment index is 1-based (segment 1 continues from segment 0)

---

## polyzymd info

Display PolyzyMD installation and dependency information.

```bash
polyzymd info
```

### Example Output

```
PolyzyMD - Molecular Dynamics for Enzyme-Polymer Systems
Version: 0.1.0

Dependencies:
  OpenMM: 8.1.1
  OpenFF Toolkit: 0.16.0
  OpenFF Interchange: 0.3.25
  Pydantic: 2.7.1

Example configs: polyzymd/configs/examples/
```

### Use Cases

- Verify installation is complete
- Check dependency versions for troubleshooting
- Confirm GPU-enabled OpenMM is installed

---

## Environment Variables

PolyzyMD expands environment variables in configuration paths:

| Variable | Example | Description |
|----------|---------|-------------|
| `$USER` | jola3134 | Current username |
| `$HOME` | /home/jola3134 | Home directory |
| `~` | /home/jola3134 | Home directory shortcut |
| `${VAR}` | - | Any environment variable |

### Example

```yaml
output:
  projects_directory: "/projects/$USER/polyzymd"
  scratch_directory: "/scratch/alpine/$USER/simulations"
```

---

## Exit Codes

| Code | Meaning |
|------|---------|
| 0 | Success |
| 1 | Error (validation failure, build failure, etc.) |

---

## See Also

- {doc}`quickstart` - Getting started tutorial
- {doc}`configuration` - Configuration file reference
- {doc}`hpc_slurm` - HPC and SLURM guide
