# CLI Reference

Complete reference for all PolyzyMD command-line interface commands.

## Global Options

All commands support these global options (placed **before** the subcommand name):

```bash
polyzymd --version                # Show version and exit
polyzymd --help                   # Show top-level help
polyzymd <command> --help         # Show subcommand help
polyzymd -v <command>             # Enable verbose output
polyzymd --openff-logs <command>  # Show OpenFF toolkit logging
polyzymd --no-color <command>     # Disable colored output
```

> **Note:** Global options must appear *before* the subcommand.
> For example: `polyzymd --no-color check-progress -c config.yaml`
> (not `polyzymd check-progress --no-color -c config.yaml`).
> `--version` prints the installed version and exits immediately.
> `--help` shows top-level help; use `polyzymd <command> --help` for subcommand help.
> `--verbose`/`-v` enables verbose output.
> `--openff-logs` enables OpenFF toolkit logging.
> `--no-color` disables colored output.

### Colored Output

PolyzyMD uses per-module colored logging to help you visually distinguish
which subsystem (building, simulation, workflow, etc.) produced each log
line. Colors are auto-detected based on your terminal capabilities and can
be disabled with `--no-color` or the `NO_COLOR` environment variable.

See the [Colored Logging Guide](../explanation/colored_logging.md) for full details
including the color table, terminal support levels, and HPC notes.

### Logging Behavior

By default, PolyzyMD suppresses verbose log messages from OpenFF Interchange and Toolkit libraries. These libraries generate per-atom INFO messages during system building (e.g., "Preset charges applied to atom index 8667" or "Key collision with different parameters"). For large systems with tens of thousands of atoms, this can produce millions of log lines.

**Default behavior:** OpenFF INFO logs are suppressed; only WARNING and ERROR messages are shown.

**To enable OpenFF logs for debugging:**

```bash
polyzymd --openff-logs build -c config.yaml
polyzymd --openff-logs run -c config.yaml --engine gromacs
```

**OpenFF logs:** OpenFF Interchange and Toolkit libraries are suppressed by default (they generate per-atom INFO messages during system building). Use `--openff-logs` to enable them for debugging force field issues.
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
polyzymd build -c <path> -r <replicates>
polyzymd build -c <path> --format gromacs    # Export for GROMACS
```

### Options

| Option | Short | Required | Default | Description |
|--------|-------|----------|---------|-------------|
| `--config` | `-c` | Yes | - | Path to YAML configuration file |
| `--replicates` | `-r` | No | "1" | Replicate range (for example "1", "1-3", "1,3,5") |
| `--scratch-dir` | - | No | from config | Override scratch directory |
| `--projects-dir` | - | No | from config | Override projects directory |
| `--output-dir` | `-o` | No | from config | Alias for --scratch-dir |
| `--dry-run` | - | No | false | Validate only, don't build |
| `--format` | - | No | OpenMM | Export format (`gromacs`, `lammps` (planned), or `amber` (planned)) |

> **Note:** `--replicate` and `--gromacs` are retained as hidden deprecated aliases.

### Example

```bash
# Build replicate 1 for OpenMM
polyzymd build -c config.yaml -r 1

# Build replicates 1 through 3 for OpenMM
polyzymd build -c config.yaml -r 1-3

# Build with custom output directory
polyzymd build -c config.yaml -r 1-3 --scratch-dir ./test_output

# Dry run to check configuration
polyzymd build -c config.yaml --dry-run

# Export to GROMACS format
polyzymd build -c config.yaml -r 1 --format gromacs
```

### Output Files (OpenMM)

The build command creates:
- `solvated_system.pdb` - Complete system with water and ions
- `system.xml` - OpenMM serialized system with restraints

### Output Files (GROMACS)

With `--format gromacs`, the build command creates in `{projects_dir}/replicate_{N}/gromacs/`:
- `{system}.gro` - GROMACS coordinate file
- `{system}.top` - GROMACS topology file
- `*.itp` - Molecule parameter files (one per component)
- `em.mdp` - Energy minimization parameters
- `eq_XX_name.mdp` - Equilibration stage parameters
- `prod.mdp` - Production parameters
- Position restraints (`#ifdef POSRES_PROTEIN`, etc.) appended into molecule `.itp` files
- `run_{system}_gromacs.sh` - Shell script to run the workflow

---

## polyzymd run

Build and run a complete local simulation with OpenMM or GROMACS.

Builds the system and executes the selected local engine workflow:

- `--engine gromacs` exports GROMACS files and runs the full GROMACS workflow
- `--engine openmm` builds and runs the OpenMM simulation locally

```bash
polyzymd run -c <path> --engine <gromacs|openmm> [options]
polyzymd run -c <path> --engine gromacs --gmx-path /usr/local/gromacs/bin/gmx
polyzymd run -c <path> --engine openmm --dry-run
```

### Options

| Option | Short | Required | Default | Description |
|--------|-------|----------|---------|-------------|
| `--config` | `-c` | Yes | - | Path to YAML configuration file |
| `--replicates` | `-r` | No | "1" | Replicate range (for example "1", "1-3", "1,3,5") |
| `--engine` | - | Yes | - | Local engine to run: `gromacs` or `openmm` |
| `--scratch-dir` | - | No | from config | Override scratch directory |
| `--projects-dir` | - | No | from config | Override projects directory |
| `--gmx-path` | - | No | unset | Path to GROMACS executable (gromacs engine only) |
| `--dry-run` | - | No | false | Validate and preview actions only (writes nothing) |

> **Note:** `--replicate` is retained as a hidden deprecated alias.

### Example

```bash
# Run full GROMACS workflow locally
polyzymd run -c config.yaml -r 1-3 --engine gromacs

# Use custom GROMACS installation
polyzymd run -c config.yaml --engine gromacs --gmx-path /usr/local/gromacs/bin/gmx

# Run a full OpenMM simulation locally
polyzymd run -c config.yaml -r 1 --engine openmm

# Preview without writing files
polyzymd run -c config.yaml -r 1-3 --engine gromacs --dry-run
```

### Workflow

1. Load and validate configuration
2. Build system (enzyme + substrate + polymers + solvent)
3. Run selected engine workflow:
   - GROMACS: export `.gro/.top/.mdp` then run EM/equilibration/production/post-processing
   - OpenMM: run minimization/equilibration/production locally

GROMACS output is streamed in real-time for familiar user experience.
On any failure, execution stops immediately and intermediate files are preserved.

### Notes

- Requires GROMACS only when `--engine gromacs` is selected
- Use `--gmx-path` only with `--engine gromacs`
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
├── *.itp                     # Molecule parameters (one per component)
├── em.mdp                    # Energy minimization parameters
├── eq_01_heating.mdp         # Equilibration stage 1
├── eq_02_free_equilibration.mdp  # Equilibration stage 2
├── prod.mdp                  # Production parameters
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
recovery, and resubmission. See {doc}`../how_to/hpc_slurm` for details.

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
| `--engine` | - | No | from config or openmm | Simulation engine: `gromacs` or `openmm` |
| `--email` | - | No | "" | Email for job notifications |
| `--scratch-dir` | - | No | from config | Override scratch directory |
| `--projects-dir` | - | No | from config | Override projects directory |
| `--output-dir` | - | No | auto | Directory for job scripts |
| `--time-limit` | - | No | from preset | Override SLURM time limit (HH:MM:SS) |
| `--memory` | - | No | 3G | Override SLURM memory allocation |
| `--account` | - | No | - | Override SLURM account / allocation ID |
| `--partition` | - | No | from preset | Override SLURM partition |
| `--qos` | - | No | - | Override SLURM QoS |
| `--gpu-type` | - | No | - | GPU type for GRES (e.g., "a100", "a40", "mi100") |
| `--constraint` | - | No | - | SLURM `--constraint` for node features (e.g., "A40", "A40\|A100") |
| `--nodelist` | - | No | - | SLURM `--nodelist` override (e.g., "gpu-node-001") |
| `--pixi-env` | - | No | from preset | Pixi environment for SLURM jobs (`cuda-12-4` or `cuda-12-6`) |
| `--skip-build` | - | No | false | Skip system building (use pre-built system from `polyzymd build`) |
| `--force` | - | No | false | Skip duplicate-job check |
| `--openff-logs` | - | No | false | Enable verbose OpenFF logs in job scripts |
| `--dry-run` | - | No | false | Preview submission plan only (no files written, no submission) |
| `--generate-only` | - | No | false | Generate SLURM scripts without submitting (the previous `--dry-run` behavior) |

:::{note}
`--dry-run` and `--generate-only` are mutually exclusive. Use `--dry-run` to
preview the submission plan without writing any files. Use `--generate-only` to
generate SLURM scripts for inspection without submitting them to the scheduler.
:::

### SLURM Presets

| Preset | Partition | Time Limit | Description |
|--------|-----------|------------|-------------|
| `aa100` | aa100 | 24:00:00 | NVIDIA A100 GPUs |
| `al40` | al40 | 24:00:00 | NVIDIA L40 GPUs |
| `blanca-shirts` | blanca-shirts | 7-00:00:00 | Blanca condo partition |
| `bridges2` | GPU | 48:00:00 | PSC Bridges2 GPU |
| `testing` | atesting | 01:00:00 | Quick tests |

### Example

```bash
# Preview what would be submitted (no files written)
polyzymd submit -c config.yaml -r 1-5 --preset aa100 --dry-run

# Generate scripts for inspection without submitting
polyzymd submit -c config.yaml -r 1-5 --preset aa100 --generate-only

# Submit for real with email notifications
polyzymd submit -c config.yaml -r 1-5 --preset aa100 --email you@university.edu

# Quick test with short time limit
polyzymd submit -c config.yaml -r 1 --preset testing --time-limit 0:05:00

# Custom directories for HPC
polyzymd submit -c config.yaml -r 1-3 --preset aa100 \
    --scratch-dir /scratch/alpine/$USER/sims \
    --projects-dir /projects/$USER/polyzymd

# GROMACS GPU submission with constraint
polyzymd submit -c config.yaml -r 1-3 \
    --engine gromacs \
    --preset blanca-shirts \
    --constraint "A40" \
    --email you@university.edu

# GROMACS CPU submission
polyzymd submit -c config.yaml -r 1-3 \
    --engine gromacs \
    --preset aa100
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
what work remains. See {doc}`../how_to/hpc_slurm` for details.

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
- For a visual overview of all replicates, use `polyzymd status` instead

---

(cli-status)=
## polyzymd status

Show a compact progress overview for all replicates of a simulation.
Auto-detects replicate directories on disk and displays colored progress
bars with completion percentage, nanoseconds completed, and simulation
status.

### Usage

```bash
polyzymd status -c config.yaml
```

### Options

| Option | Required | Description |
|--------|----------|-------------|
| `-c, --config PATH` | Yes | Path to YAML configuration file |

### Output Format

```
  polyzymd status — fnIII_apo_OEGMA-SBMA_A50_B50_100ns_310K
  ──────────────────────────────────────────────────────

  run1  ██████████████████████████████████████████  100.0%  100.0/100.0 ns  completed
  run2  █████████████████████░░░░░░░░░░░░░░░░░░░░░   50.2%   50.2/100.0 ns  running
  run3  ███████████████░░░░░░░░░░░░░░░░░░░░░░░░░░░   35.0%   35.0/100.0 ns  interrupted
  run4  ░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░    0.0%    0.0/100.0 ns  not_started

  1/4 need attention (recover with: polyzymd recover -c config.yaml -r <N> --submit)
```

### Status Colors

| Status | Color | Meaning |
|--------|-------|---------|
| `completed` | Green | Production run finished |
| `running` | Cyan | Currently executing |
| `interrupted` | Amber | Stalled — needs `polyzymd recover` |
| `failed` | Red | Error occurred |
| `not_started` | Gray | Directory exists but no progress data |
| `not found` | Gray | Expected directory not on disk |

### Notes

- This is a **read-only** command — it only reads `progress.json` files
- Replicate directories are auto-detected via the naming template in the config
- The command is a one-shot snapshot (prints and exits)
- Use `polyzymd recover -c config.yaml -r <N> --submit` to resume interrupted replicates

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
| `--engine` | - | No | from config | Override simulation engine (`gromacs` or `openmm`) |
| `--submit / --no-submit` | - | No | --no-submit | Submit a recovery job (default: status only) |
| `--dry-run` | - | No | false | Show what would be submitted without submitting |
| `--email` | - | No | "" | Email for job notifications |
| `--memory` | - | No | 3G | Override SLURM memory allocation (e.g. '4G', '8G') |
| `--partition` | - | No | from preset | Override SLURM partition |
| `--qos` | - | No | - | Override SLURM QoS |
| `--constraint` | - | No | - | SLURM `--constraint` for node features (e.g., "A40", "A40\|A100") |
| `--nodelist` | - | No | - | SLURM `--nodelist` override (e.g., "gpu-node-001") |
| `--pixi-env` | - | No | from preset | Pixi environment for recovery job |
| `--force` | - | No | false | Skip duplicate-job check |

### Example

```bash
# Check status only
polyzymd recover -c config.yaml -r 1

# Submit a recovery job
polyzymd recover -c config.yaml -r 1 --submit --preset blanca-shirts

# GROMACS recovery with GPU constraint
polyzymd recover -c config.yaml -r 1 \
    --engine gromacs \
    --submit \
    --preset blanca-shirts \
    --constraint "A40"

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

Example configs: polyzymd/templates/examples/
```

### Use Cases

- Verify installation is complete
- Check dependency versions for troubleshooting
- Confirm GPU-enabled OpenMM is installed

---

(polyzymd-new-analysis)=
## polyzymd new-analysis

Scaffold an analysis plugin and matching tests.

```bash
polyzymd new-analysis NAME [OPTIONS]

Options:
  --class-name TEXT                 PascalCase class prefix
  --style [dict]                    Advanced package style
  --advanced                        Request the advanced MDAnalysis-native package scaffold
  --project-root DIRECTORY          Repository root
  --force                           Overwrite existing files
  --dry-run                         Print paths without writing files
```

By default, the command creates a single-file MDAnalysis-native plugin at
`src/polyzymd/analyses/<NAME>.py`. The generated analysis subclasses `Analysis`,
builds an `MDAAnalysisJob.from_function()` job, returns a `ReplicateArtifact`
with explicit `payload["metrics"]`, and relies on the default artifact
aggregation path.

Omit `--style` to use the default single-file scaffold.

`--advanced` and `--style dict` create an advanced package scaffold at
`src/polyzymd/analyses/<NAME>/` with lifecycle wiring in `__init__.py`, a
lazy-imported `AnalysisBase` helper in `_mda.py`, and dict metrics stored in
canonical artifacts.

```bash
polyzymd new-analysis solvent_shell
polyzymd new-analysis solvent_shell --advanced
polyzymd new-analysis solvent_shell --style dict
pixi run -e build pytest tests/analyses/plugins/test_solvent_shell.py -v
```

---

## polyzymd compare

Compare analysis results across multiple simulation conditions with statistical testing.

```bash
polyzymd compare COMMAND [OPTIONS]

Commands:
  init      Initialize a new comparison project
  validate  Validate comparison configuration
  run       Run a comparison by analysis type
  run-all   Run all enabled comparisons
  plot-all  Generate comparison plots from a workspace
  submit    Submit analysis as SLURM job DAG (HPC)
  submit-all Submit all enabled analyses as dependency-ordered SLURM DAGs
  status    Show status of submitted SLURM analysis jobs
  finalize  Run comparison + plotting from aggregated on-disk results
```

### polyzymd compare init

Create a new comparison project with template configuration.

```bash
polyzymd compare init -n NAME [OPTIONS]

Options:
  -n, --name TEXT       Project name (creates directory) [required]
  --eq-time TEXT         Default equilibration time [default: 10ns]
  -o, --output-dir PATH  Parent directory [default: current]
```

#### Example

```bash
polyzymd compare init -n polymer_study
cd polymer_study
# Edit comparison.yaml to add your conditions
```

### polyzymd compare run

Run a single analysis comparison by type. This is a generic command that works
with any discovered analysis plugin.

```bash
polyzymd compare run COMPARISON_TYPE [OPTIONS]

Arguments:
  COMPARISON_TYPE        Analysis plugin name (e.g. rmsf, contacts, distances)

Options:
  -f, --file PATH        Path to comparison.yaml [default: comparison.yaml]
  --eq-time TEXT          Override equilibration time (e.g. '10ns', '5000ps')
  --recompute            Force recompute even if cached results exist
  --format TEXT           Output format: table, markdown, json [default: table]
  -o, --output PATH      Save formatted output to file
  -q, --quiet            Suppress INFO messages
  --debug                Enable DEBUG logging
  --list                 List available comparison types and exit
```

#### Example

```bash
# Run RMSF comparison (uses plugins.rmsf from comparison.yaml)
polyzymd compare run rmsf

# Override equilibration time
polyzymd compare run rmsf --eq-time 20ns

# Run contacts comparison with markdown output
polyzymd compare run contacts --format markdown -o report.md

# List all available analysis types
polyzymd compare run --list
```

### polyzymd compare validate

Validate a comparison.yaml configuration file without running analyses.

```bash
polyzymd compare validate [OPTIONS]

Options:
  -f, --file PATH        Path to comparison.yaml [default: comparison.yaml]
  --format [table|json]  Output format [default: table]
```

#### What It Checks

- YAML syntax and structure
- Required fields present
- At least 1 condition defined
- Condition labels are unique
- Control label matches a condition (if specified)
- Config files exist for each condition

#### Example

```bash
# Basic validation
polyzymd compare validate

# Validate specific file
polyzymd compare validate -f path/to/comparison.yaml

# JSON output for CI integration
polyzymd compare validate --format json
```

**Output (success):**
```
Validating: /path/to/comparison.yaml

✓ Configuration is valid

  Name: polymer_study
  Conditions: 3
    - WT, PEG, SBMA
  Control: WT
  Analysis sections: rmsf, catalytic_triad
```

**Output (errors):**
```
Validating: /path/to/comparison.yaml

✗ Configuration has errors

  • Control 'NoPolymer' not found in conditions: ['WT', 'PEG']
  • Config file not found: /path/to/missing/config.yaml
```

**JSON output:**
```json
{
  "file": "/path/to/comparison.yaml",
  "valid": true,
  "errors": [],
  "summary": {
    "name": "polymer_study",
    "conditions_count": 3,
    "condition_labels": ["WT", "PEG", "SBMA"],
    "control": "WT",
    "sections_configured": ["rmsf", "catalytic_triad"]
  }
}
```

### polyzymd compare plot-all

Generate configured plots from a comparison workspace.

```bash
polyzymd compare plot-all [OPTIONS]

Options:
  -f, --file PATH                 Path to comparison.yaml [default: comparison.yaml]
  -o, --output-dir PATH           Override plot output directory
  -a, --analysis TEXT             Plot one analysis type only
  --list-available                List registered/available plots and exit
  -q, --quiet                     Suppress INFO messages
  --debug                         Enable DEBUG logging
```

#### Example

```bash
# Generate all configured plots
polyzymd compare plot-all

# High-level availability check
polyzymd compare plot-all --list-available

# One analysis only
polyzymd compare plot-all -a rmsf
```

### polyzymd compare submit

Submit a replicate-level SLURM analysis DAG for one plugin. Each replicate runs
as an independent SLURM job, followed by per-condition aggregation jobs and a
final comparison + plotting job.

Before submission, this command runs a dependency preflight check: if the
target plugin declares `dependencies`, required upstream comparison results must
already exist on disk (or use `compare submit-all` instead).

```bash
polyzymd compare submit ANALYSIS [OPTIONS]

Arguments:
  ANALYSIS               Analysis plugin name (e.g. rmsf, contacts)

Options:
  -f, --file PATH        Path to comparison.yaml [default: comparison.yaml]
  --partition TEXT        SLURM partition [default: cluster default]
  --qos TEXT             SLURM QoS
  --account TEXT         SLURM account/allocation
  --pixi-path TEXT       Path to pixi executable [default: pixi]
  --ntasks INT           SLURM ntasks [default: 1]
  --cpus-per-task INT    SLURM cpus-per-task [default: 1]
  --mem TEXT             SLURM memory request [default: 4G]
  --time TEXT            SLURM walltime [default: 01:00:00]
  --max-retries INT      Max retries for failed jobs [default: 3]
  --mail-user TEXT       Email for failure notifications
  --recompute            Force recomputation in workers
  --allow-partial        Allow finalize when some conditions are missing results
  --equilibration TEXT   Override equilibration time
  --dry-run              Generate scripts without submitting jobs
  --job-arrays           Submit one SLURM array job per condition
```

#### Example

```bash
# Submit RMSF analysis to SLURM
polyzymd compare submit rmsf --partition gpu --account my_alloc

# Dry run to inspect generated scripts
polyzymd compare submit contacts --dry-run

# Use job arrays for efficiency
polyzymd compare submit rmsf --job-arrays --partition aa100

# Rely on plugin memory hints and cluster default partition
polyzymd compare submit secondary_structure --qos normal

# Blanca condo nodes at CU Boulder
module load slurm/blanca
polyzymd compare submit sasa \
    -f comparison.yaml \
    --partition blanca-shirts \
    --account blanca-shirts \
    --qos blanca-shirts \
    --mem 8G \
    --time 02:00:00
```

### polyzymd compare submit-all

Submit all enabled analyses from `comparison.yaml` in dependency order with
cross-plugin finalize dependencies.

```bash
polyzymd compare submit-all [OPTIONS]

Options:
  -f, --file PATH         Path to comparison.yaml [default: comparison.yaml]
  --partition TEXT        SLURM partition [default: cluster default]
  --qos TEXT              SLURM QoS
  --account TEXT          SLURM account/allocation
  --pixi-path TEXT        Path to pixi executable [default: pixi]
  --ntasks INT            SLURM ntasks [default: 1]
  --cpus-per-task INT     SLURM cpus-per-task [default: 1]
  --mem TEXT              SLURM memory request [default: 4G]
  --time TEXT             SLURM walltime [default: 01:00:00]
  --max-retries INT       Max retries for failed jobs [default: 3]
  --mail-user TEXT        Email for failure notifications
  --recompute             Force recomputation in workers
  --allow-partial         Allow finalize when some conditions are missing results
  --equilibration TEXT    Override equilibration time
  --dry-run               Generate scripts without submitting jobs
  --exclude TEXT          Exclude one analysis (repeatable)
```

#### Example

```bash
# Submit everything enabled in comparison.yaml
polyzymd compare submit-all -f comparison.yaml --partition aa100 --qos normal

# Skip selected plugins
polyzymd compare submit-all -f comparison.yaml --exclude sasa --exclude hydrogen_bonds

# Dry-run planning only
polyzymd compare submit-all -f comparison.yaml --dry-run
```

### polyzymd compare status

Show the status of a submitted SLURM analysis DAG. Reports counts of pending,
running, succeeded, and failed jobs.

```bash
polyzymd compare status ANALYSIS [OPTIONS]

Arguments:
  ANALYSIS               Analysis plugin name

Options:
  -f, --file PATH        Path to comparison.yaml [default: comparison.yaml]
  --reconcile            Reconcile status files with sacct before reporting
  --json                 Print machine-readable JSON status
```

#### Example

```bash
# Check status of RMSF SLURM jobs
polyzymd compare status rmsf

# Reconcile with SLURM scheduler and get JSON output
polyzymd compare status rmsf --reconcile --json
```

### polyzymd compare finalize

Run comparison and plotting from aggregated on-disk results. Use this after
SLURM jobs complete, or to re-run comparison/plotting without recomputing
per-replicate results.

```bash
polyzymd compare finalize ANALYSIS [OPTIONS]

Arguments:
  ANALYSIS               Analysis plugin name

Options:
  -f, --file PATH        Path to comparison.yaml [default: comparison.yaml]
  --recompute            Retained for CLI compatibility (no effect)
  --allow-partial        Allow finalize when some conditions are missing results
```

#### Example

```bash
# Finalize after SLURM jobs complete
polyzymd compare finalize rmsf

# Allow partial results (some conditions may have failed)
polyzymd compare finalize contacts --allow-partial
```

---

## Plotting comparisons (`polyzymd compare plot-all`)

The standalone `polyzymd plot` command group was removed.
Use `polyzymd compare plot-all` for all comparison plotting workflows.

```bash
polyzymd compare plot-all [OPTIONS]

Options:
  -f, --file PATH                 Path to comparison.yaml [default: comparison.yaml]
  -o, --output-dir PATH           Override plot output directory
  -a, --analysis TEXT             Plot one analysis type only
  --list-available                List registered/available plots and exit
  -q, --quiet                     Suppress INFO messages
  --debug                         Enable DEBUG logging
```

### Example

```bash
# Generate all configured plots
polyzymd compare plot-all

# Show available plots
polyzymd compare plot-all --list-available

# Plot only RMSF
polyzymd compare plot-all -a rmsf
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
| 99 | Graceful shutdown — simulation was interrupted but interrupted state was saved (see {doc}`../how_to/hpc_slurm`) |

---

## See Also

- {doc}`../get_started/quickstart` - Getting started tutorial
- {doc}`configuration` - Configuration file reference
- {doc}`../how_to/hpc_slurm` - HPC and SLURM guide
- {doc}`../how_to/analysis_rmsf_quickstart` - RMSF analysis tutorial
- {doc}`../how_to/analysis_compare_conditions` - Comparing simulation conditions
