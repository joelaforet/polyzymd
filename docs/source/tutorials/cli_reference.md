# CLI Reference

Complete reference for all PolyzyMD command-line interface commands.

## Global Options

All commands support these global options:

```bash
polyzymd --version        # Show version
polyzymd --help           # Show help
polyzymd -v <command>     # Verbose output (debug logging)
polyzymd --openff-logs <command>  # Enable verbose OpenFF logs
```

### Logging Behavior

By default, PolyzyMD suppresses verbose log messages from OpenFF Interchange and Toolkit libraries. These libraries generate per-atom INFO messages during system building (e.g., "Preset charges applied to atom index 8667" or "Key collision with different parameters"). For large systems with tens of thousands of atoms, this can produce millions of log lines.

**Default behavior:** OpenFF INFO logs are suppressed; only WARNING and ERROR messages are shown.

**To enable OpenFF logs for debugging:**

```bash
polyzymd --openff-logs build -c config.yaml
polyzymd --openff-logs run -c config.yaml
```

This is useful when:
- Debugging force field parameter assignment issues
- Investigating charge assignment problems
- Troubleshooting system building failures

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
polyzymd build -c <path> --gromacs    # Export for GROMACS
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
| `--gromacs` | - | No | false | Export to GROMACS format instead of OpenMM |

### Example

```bash
# Build replicate 1 for OpenMM
polyzymd build -c config.yaml -r 1

# Build with custom output directory
polyzymd build -c config.yaml -r 1 --scratch-dir ./test_output

# Dry run to check configuration
polyzymd build -c config.yaml --dry-run

# Export to GROMACS format
polyzymd build -c config.yaml -r 1 --gromacs
```

### Output Files (OpenMM)

The build command creates:
- `solvated_system.pdb` - Complete system with water and ions
- `system.xml` - OpenMM serialized system
- `topology.json` - Topology information
- `build.log` - Build process log

### Output Files (GROMACS)

With `--gromacs`, the build command creates in `{projects_dir}/replicate_{N}/gromacs/`:
- `{system}.gro` - GROMACS coordinate file
- `{system}.top` - GROMACS topology file
- `em.mdp` - Energy minimization parameters
- `eq_XX_name.mdp` - Equilibration stage parameters
- `prod.mdp` - Production parameters
- `posre_*.itp` - Position restraint files (if configured)
- `run_{system}_gromacs.sh` - Shell script to run the workflow

---

## polyzymd run

Build and run a complete simulation.

By default, runs using OpenMM. Use `--gromacs` to run using GROMACS instead.

```bash
# OpenMM (default)
polyzymd run --config <path> [options]
polyzymd run -c <path> -r <replicate>

# GROMACS
polyzymd run -c <path> --gromacs
polyzymd run -c <path> --gromacs --gmx-path /usr/local/gromacs/bin/gmx
polyzymd run -c <path> --gromacs --dry-run
```

### Options

| Option | Short | Required | Default | Description |
|--------|-------|----------|---------|-------------|
| `--config` | `-c` | Yes | - | Path to YAML configuration file |
| `--replicate` | `-r` | No | 1 | Replicate number |
| `--scratch-dir` | - | No | from config | Override scratch directory |
| `--projects-dir` | - | No | from config | Override projects directory |
| `--segment-time` | - | No | auto | Override production time per segment (ns) [OpenMM only] |
| `--segment-frames` | - | No | auto | Override frames per segment [OpenMM only] |
| `--skip-build` | - | No | false | Skip building (use existing system) [OpenMM only] |
| `--gromacs` | - | No | false | Run using GROMACS instead of OpenMM |
| `--gmx-path` | - | No | "gmx" | Path to GROMACS executable [GROMACS only] |
| `--dry-run` | - | No | false | Export files but don't run simulation [GROMACS only] |

### Example (OpenMM)

```bash
# Run locally (useful for testing)
polyzymd run -c config.yaml -r 1

# Run with shorter segment for testing
polyzymd run -c config.yaml -r 1 --segment-time 0.1 --segment-frames 10

# Skip building (use pre-built system)
polyzymd run -c config.yaml -r 1 --skip-build
```

### Example (GROMACS)

```bash
# Run full GROMACS workflow locally
polyzymd run -c config.yaml -r 1 --gromacs

# Use custom GROMACS installation
polyzymd run -c config.yaml --gromacs --gmx-path /usr/local/gromacs/bin/gmx

# Export files only (for manual execution or HPC)
polyzymd run -c config.yaml --gromacs --dry-run
```

### OpenMM Workflow

1. Load and validate configuration
2. Build system (enzyme + substrate + polymers + solvent)
3. Apply restraints (if configured)
4. Run energy minimization
5. Run equilibration (NVT/NPT stages)
6. Run first production segment (NPT)

### GROMACS Workflow

1. Load and validate configuration
2. Build system (same as OpenMM)
3. Export to GROMACS format (.gro, .top, .mdp files)
4. Run energy minimization (grompp + mdrun)
5. Run equilibration stages (grompp + mdrun for each stage)
6. Run production MD (grompp + mdrun)
7. Post-process trajectory (trjconv for PBC handling)

All GROMACS output is streamed in real-time for familiar user experience.
On any failure, execution stops immediately and all intermediate files are
preserved for debugging.

### GROMACS Notes

- Requires GROMACS to be installed and accessible via PATH
- Use `--gmx-path` to specify a custom GROMACS executable location
- MDP parameters are generated from your config.yaml to match OpenMM settings
- OpenFF force field defaults are used (rcoulomb=0.9, rvdw=0.9, PME) for 1:1 parity with OpenMM
- Position restraints are automatically generated for equilibration stages
- Post-processing creates `prod_nojump.xtc` and `prod_centered.xtc` trajectories

### GROMACS Output Files

Files are created in `{projects_dir}/replicate_{N}/gromacs/`:

```
gromacs/
├── {system}.gro              # Initial coordinates
├── {system}.top              # Topology
├── em.mdp                    # Energy minimization parameters
├── eq_01_nvt.mdp             # Equilibration stage 1 (NVT)
├── eq_02_npt.mdp             # Equilibration stage 2 (NPT)
├── prod.mdp                  # Production parameters
├── posre_*.itp               # Position restraint files
├── run_{system}_gromacs.sh   # Generated run script
├── em.tpr, em.gro, em.edr    # Energy minimization outputs
├── eq_01.*, eq_02.*          # Equilibration outputs
├── prod.tpr, prod.xtc, ...   # Production outputs
├── prod_nojump.xtc           # Trajectory without PBC jumps
└── prod_centered.xtc         # Centered trajectory for visualization
```

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
| `--memory` | - | No | 3G | Override SLURM memory allocation |
| `--openff-logs` | - | No | false | Enable verbose OpenFF logs in job scripts |
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

## polyzymd analyze

Analyze MD trajectories with various metrics.

```bash
polyzymd analyze COMMAND [OPTIONS]

Commands:
  rmsf       Compute RMSF (Root Mean Square Fluctuation) analysis
  distances  Compute inter-atomic distance analysis
```

### polyzymd analyze rmsf

Compute per-residue flexibility from MD trajectories.

```bash
polyzymd analyze rmsf --config <path> --replicates <range> [OPTIONS]
polyzymd analyze rmsf -c <path> -r 1-5 --eq-time 100ns
```

#### Options

| Option | Short | Required | Default | Description |
|--------|-------|----------|---------|-------------|
| `--config` | `-c` | Yes | - | Path to YAML configuration file |
| `--replicates` | `-r` | No | "1" | Replicate specification: "1-5", "1,3,5", or "1" |
| `--eq-time` | - | No | "0ns" | Equilibration time to skip: "100ns", "5000ps" |
| `--selection` | - | No | "protein and name CA" | MDAnalysis selection for RMSF atoms |
| `--reference-mode` | - | No | "centroid" | Reference structure: "centroid", "average", or "frame" |
| `--reference-frame` | - | No | - | Frame index when --reference-mode=frame (1-indexed) |
| `--alignment-selection` | - | No | "protein and name CA" | Selection for trajectory alignment |
| `--centroid-selection` | - | No | "protein" | Selection for centroid finding |
| `--plot` | - | No | false | Generate plot after analysis |
| `--recompute` | - | No | false | Force recompute even if cached |
| `--output-dir` | `-o` | No | auto | Custom output directory |

#### Reference Modes

| Mode | Description | Use Case |
|------|-------------|----------|
| `centroid` | Most populated conformational state (K-Means) | Equilibrium flexibility analysis |
| `average` | Mathematical mean structure | Pure thermal fluctuations |
| `frame` | User-specified frame number | Functional state analysis |

#### Example

```bash
# Basic RMSF analysis
polyzymd analyze rmsf -c config.yaml -r 1 --eq-time 10ns

# Multiple replicates with plot
polyzymd analyze rmsf -c config.yaml -r 1-3 --eq-time 10ns --plot

# Custom reference structure
polyzymd analyze rmsf -c config.yaml -r 1 --reference-mode average

# Specific frame as reference (e.g., catalytically competent)
polyzymd analyze rmsf -c config.yaml -r 1 --reference-mode frame --reference-frame 500
```

#### Output

Results are saved in JSON format:

```
<projects_directory>/
└── analysis/
    └── rmsf/
        ├── run_1/rmsf_eq10ns.json
        ├── run_2/rmsf_eq10ns.json
        └── aggregated/rmsf_reps1-3_eq10ns.json
```

---

## polyzymd compare

Compare analysis results across multiple simulation conditions with statistical testing.

```bash
polyzymd compare COMMAND [OPTIONS]

Commands:
  init      Initialize a new comparison project
  rmsf      Compare RMSF across conditions
  validate  Validate comparison configuration
  show      Display saved comparison results
  plot      Generate comparison plots
```

### polyzymd compare init

Create a new comparison project with template configuration.

```bash
polyzymd compare init NAME [OPTIONS]

Arguments:
  NAME                   Project name (creates directory)

Options:
  --eq-time TEXT         Default equilibration time [default: 10ns]
  -o, --output-dir PATH  Parent directory [default: current]
```

#### Example

```bash
polyzymd compare init polymer_study
cd polymer_study
# Edit comparison.yaml to add your conditions
```

### polyzymd compare rmsf

Run statistical comparison of RMSF across conditions.

```bash
polyzymd compare rmsf [OPTIONS]
```

**Requires** an `rmsf:` section in comparison.yaml. This is the single source
of truth for RMSF settings.

#### Options

| Option | Short | Default | Description |
|--------|-------|---------|-------------|
| `--file` | `-f` | comparison.yaml | Path to comparison config file |
| `--eq-time` | - | from config | Override equilibration time |
| `--override` | - | false | Enable CLI overrides for RMSF settings |
| `--selection` | - | from rmsf config | Override atom selection (requires --override) |
| `--reference-mode` | - | from rmsf config | Override reference mode (requires --override) |
| `--reference-frame` | - | from rmsf config | Override reference frame (requires --override) |
| `--recompute` | - | false | Force recompute RMSF |
| `--format` | - | table | Output format: table, markdown, json |
| `--output` | `-o` | - | Save formatted output to file |
| `--verbose` | `-v` | false | Show detailed logging |

#### Example

```bash
# Run comparison with default settings (uses rmsf: section from YAML)
polyzymd compare rmsf

# Override equilibration time
polyzymd compare rmsf --eq-time 20ns

# Override RMSF settings (requires --override flag)
polyzymd compare rmsf --override --selection "protein and name CA CB"

# Output as markdown
polyzymd compare rmsf --format markdown -o report.md
```

### polyzymd compare validate

Check that a comparison.yaml configuration is valid.

```bash
polyzymd compare validate [OPTIONS]

Options:
  -f, --file PATH    Config file [default: comparison.yaml]
```

### polyzymd compare show

Display a previously saved comparison result.

```bash
polyzymd compare show RESULT_FILE [OPTIONS]

Arguments:
  RESULT_FILE                     Path to saved JSON result

Options:
  --format [table|markdown|json]  Output format [default: table]
```

### polyzymd compare plot

Generate publication-ready plots from comparison results.

```bash
polyzymd compare plot RESULT_FILE [OPTIONS]

Arguments:
  RESULT_FILE                     Path to saved comparison JSON

Options:
  -o, --output-dir PATH           Output directory [default: figures/]
  --format [png|pdf|svg]          Image format [default: png]
  --dpi INTEGER                   Resolution for PNG [default: 150]
  --summary / --no-summary        Generate summary panel [default: yes]
  --show / --no-show              Display interactively [default: no]
```

#### Generated Plots

| File | Description |
|------|-------------|
| `rmsf_comparison.png` | Bar chart of mean RMSF by condition |
| `percent_change.png` | Bar chart of % change vs control |
| `effect_sizes.png` | Forest plot of Cohen's d values |
| `summary_panel.png` | Combined 3-panel summary figure |

#### Example

```bash
# Generate all plots
polyzymd compare plot results/rmsf_comparison_my_study.json

# High resolution for publication
polyzymd compare plot results/rmsf_comparison_my_study.json --dpi 300

# PDF format with interactive preview
polyzymd compare plot results/rmsf_comparison_my_study.json --format pdf --show
```

---

## polyzymd plot

Standalone plotting for analysis results (separate from `compare plot`).

```bash
polyzymd plot COMMAND [OPTIONS]

Commands:
  rmsf       Plot and compare RMSF results
  distances  Plot and compare distance analysis results
```

### polyzymd plot rmsf

Plot RMSF analysis results, optionally comparing multiple conditions.

```bash
polyzymd plot rmsf --inputs <files> [OPTIONS]

Options:
  --inputs PATH...    One or more RMSF result JSON files [required]
  --labels TEXT...    Labels for each input (same order as inputs)
  -o, --output PATH   Output file path
  --format TEXT       Output format: png, pdf, svg [default: png]
```

#### Example

```bash
# Single condition
polyzymd plot rmsf --inputs analysis/rmsf/run_1/rmsf_eq10ns.json

# Compare two conditions
polyzymd plot rmsf \
    --inputs no_polymer/analysis/rmsf/aggregated/rmsf_reps1-3_eq10ns.json \
             with_polymer/analysis/rmsf/aggregated/rmsf_reps1-3_eq10ns.json \
    --labels "No Polymer" "With Polymer" \
    -o comparison.png
```

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
- {doc}`analysis_rmsf_quickstart` - RMSF analysis tutorial
- {doc}`analysis_compare_conditions` - Comparing simulation conditions
