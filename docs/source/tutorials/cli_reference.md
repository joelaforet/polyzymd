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
polyzymd --openff-logs run-gromacs -c config.yaml
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

## polyzymd run-gromacs

Build and run a complete simulation using GROMACS.

Builds the system, exports to GROMACS format (.gro, .top, .mdp), and
executes the full GROMACS workflow locally (energy minimization,
equilibration, production, and trajectory post-processing).

```bash
polyzymd run-gromacs -c <path> [options]
polyzymd run-gromacs -c <path> --gmx-path /usr/local/gromacs/bin/gmx
polyzymd run-gromacs -c <path> --dry-run
```

### Options

| Option | Short | Required | Default | Description |
|--------|-------|----------|---------|-------------|
| `--config` | `-c` | Yes | - | Path to YAML configuration file |
| `--replicate` | `-r` | No | 1 | Replicate number |
| `--scratch-dir` | - | No | from config | Override scratch directory |
| `--projects-dir` | - | No | from config | Override projects directory |
| `--gmx-path` | - | No | "gmx" | Path to GROMACS executable |
| `--dry-run` | - | No | false | Export files but don't run simulation |

### Example

```bash
# Run full GROMACS workflow locally
polyzymd run-gromacs -c config.yaml -r 1

# Use custom GROMACS installation
polyzymd run-gromacs -c config.yaml --gmx-path /usr/local/gromacs/bin/gmx

# Export files only (for manual execution or HPC)
polyzymd run-gromacs -c config.yaml --dry-run
```

### Workflow

1. Load and validate configuration
2. Build system (enzyme + substrate + polymers + solvent)
3. Export to GROMACS format (.gro, .top, .mdp files)
4. Run energy minimization (grompp + mdrun)
5. Run equilibration stages (grompp + mdrun for each stage)
6. Run production MD (grompp + mdrun)
7. Post-process trajectory (trjconv for PBC handling)

All GROMACS output is streamed in real-time for familiar user experience.
On any failure, execution stops immediately and all intermediate files are
preserved for debugging.

### Notes

- Requires GROMACS to be installed and accessible via PATH
- Use `--gmx-path` to specify a custom GROMACS executable location
- MDP parameters are generated from your config.yaml to match OpenMM settings
- OpenFF force field defaults are used (rcoulomb=0.9, rvdw=0.9, PME) for 1:1 parity with OpenMM
- Position restraints are automatically generated for equilibration stages
- Post-processing creates `prod_nojump.xtc` and `prod_centered.xtc` trajectories
- For OpenMM simulations, use `polyzymd run-segment` (for a single segment) or
  `polyzymd submit` (to submit self-resubmitting SLURM jobs)

### Output Files

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

Submit self-resubmitting simulation jobs to SLURM for HPC execution.

Each replicate gets one SLURM script that handles the full simulation
lifecycle: building, equilibration, production segments, interruption
recovery, and resubmission. See {doc}`hpc_slurm` for details.

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

### Self-Resubmitting Jobs

The submit command creates one self-resubmitting SLURM script per replicate:

```
  ┌─────────────────────────┐
  │  Job runs segment       │
  │  Job checks progress    │◄──── resubmits itself
  │  Job resubmits if       │      if work remains
  │  work remains           │
  └─────────────────────────┘
```

Each job is identical and idempotent — it scans the filesystem to determine
what work remains. See {doc}`hpc_slurm` for details.

---

## polyzymd run-segment

Unified entry point for SLURM jobs. Determines what work remains by loading
progress state, then runs the next segment of work.

```bash
polyzymd run-segment -c CONFIG [OPTIONS]
```

### Options

| Option | Short | Required | Default | Description |
|--------|-------|----------|---------|-------------|
| `--config` | `-c` | Yes | - | Path to YAML configuration file |
| `--replicate` | `-r` | No | 1 | Replicate number |
| `--scratch-dir` | - | No | from config | Override scratch directory |
| `--skip-build` | - | No | false | Skip system building for initial segment |

### Behavior

- If no segments exist: builds system, equilibrates, runs production segment 0
- If segments exist but simulation incomplete: continues from last completed segment
- If simulation is already complete: exits 0 immediately

### Exit Codes

| Code | Meaning |
|------|---------|
| 0 | Segment completed successfully |
| 1 | Error |
| 99 | Graceful interruption (wall-time/preemption signal) |

### Notes

- This command is called by the generated SLURM scripts, not typically by users directly
- Progress is tracked in `progress.json` in the working directory

---

## polyzymd check-progress

Check whether a simulation is complete. Used by SLURM resubmission logic
to decide whether to resubmit.

```bash
polyzymd check-progress -c CONFIG [OPTIONS]
```

### Options

| Option | Short | Required | Default | Description |
|--------|-------|----------|---------|-------------|
| `--config` | `-c` | Yes | - | Path to YAML configuration file |
| `--replicate` | `-r` | No | 1 | Replicate number |
| `--scratch-dir` | - | No | from config | Override scratch directory |

### Exit Codes

| Code | Meaning |
|------|---------|
| 0 | Simulation complete — do NOT resubmit |
| 1 | Work remains — resubmit |

### Example

```bash
polyzymd check-progress -c config.yaml -r 1

# Output:
# Progress: 50000000/50000000 steps (100.0%), 10 segment(s)
# Status: COMPLETE
```

### Notes

- This command is called by the generated SLURM scripts, not typically by users directly
- For interactive progress checking, `polyzymd recover` provides a more detailed view

---

(cli-recover)=
## polyzymd recover

Resume a stalled or interrupted simulation. Scans the working directory,
loads progress state, and reports how much work remains. With `--submit`,
generates and submits a self-resubmitting SLURM job that will automatically
continue from the last completed segment.

```bash
polyzymd recover -c CONFIG [OPTIONS]
```

### Options

| Option | Short | Required | Default | Description |
|--------|-------|----------|---------|-------------|
| `--config` | `-c` | Yes | - | Path to YAML configuration file |
| `--replicate` | `-r` | No | 1 | Replicate number |
| `--scratch-dir` | - | No | from config | Override scratch directory |
| `--preset` | - | No | aa100 | SLURM preset for recovery job |
| `--submit / --no-submit` | - | No | --no-submit | Submit a recovery job (default: status only) |
| `--dry-run` | - | No | false | Show what would be submitted without submitting |

### Example

```bash
# Check status only
polyzymd recover -c config.yaml -r 1

# Submit a recovery job
polyzymd recover -c config.yaml -r 1 --submit --preset blanca-shirts

# Dry-run (show what would be submitted)
polyzymd recover -c config.yaml -r 1 --submit --dry-run
```

### Example Output (Status Only)

```
Working directory: /scratch/user/sim/LipA_300K_run1
Progress: 12500000/50000000 steps (25.0%)
Status: in_progress
Segments: 5
  segment 0: completed (100%)
  segment 1: completed (100%)
  segment 2: completed (100%)
  segment 3: completed (100%)
  segment 4: interrupted (50%)

Remaining: 75.000 ns (37500000 steps)

To resume, run:
  polyzymd recover -c config.yaml -r 1 --submit --preset aa100
```

### Notes

- Without `--submit`, this is a read-only status report — useful for inspecting
  simulation health across replicates
- With `--submit`, generates a self-resubmitting SLURM job in
  `{working_dir}/recovery_scripts/` and submits it
- The recovery job is identical to a normal submission job — it uses `run-segment`
  to determine what work remains and continues from there

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
| 99 | Graceful shutdown — simulation was interrupted but interrupted state was saved (see {doc}`hpc_slurm`) |

---

## See Also

- {doc}`quickstart` - Getting started tutorial
- {doc}`configuration` - Configuration file reference
- {doc}`hpc_slurm` - HPC and SLURM guide
