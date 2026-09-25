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
- Referenced files (PDB, SDF, cached polymer SDFs, reaction templates) are reported as
  warnings when missing
- Monomer probabilities sum to 1.0
- Valid enum values (water model, ensemble, etc.)
- Co-solvent specification (mole_fraction XOR concentration)

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
| `--dry-run` | - | No | false | Validate only, don't build |
| `--format` | - | No | OpenMM | Export format (`gromacs`, `lammps` (planned), or `amber` (planned)) |

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
- `solvated_system.pdb` - Complete system with water and ions, for viewers
- `system.prmtop` - Analysis topology with every atom, element, charge and
  bond; MDAnalysis reads it directly and it has no atom limit
- `system.xml` - OpenMM serialized system with restraints
- `build_manifest.json` - SHA-256 hashes of the files above, the config
  hash, OpenMM and PolyzyMD versions, and the PACKMOL seeds under `provenance`

PACKMOL is seeded with the replicate number, so replicates start from
independent coordinates. The build aborts with `SolvationClashError` if packed
solvent or polymer atoms overlap the solute (see
{doc}`../how_to/troubleshooting`).

### Output Files (GROMACS)

With `--format gromacs`, the build command creates a build-only handoff in
`{projects_dir}/replicate_{N}/gromacs/`. The core handoff files are:

- `{system}.gro` - GROMACS coordinate file
- `{system}.top` - GROMACS topology file
- `*.itp` - Molecule parameter files (one per component)
- Position restraints (`#ifdef POSRES_PROTEIN`, etc.) appended into molecule `.itp` files

PolyzyMD may also generate convenience defaults:

- `em.mdp` - Energy minimization parameters
- `eq_XX_name.mdp` - Equilibration stage parameters
- `prod.mdp` - Production parameters
- `run_*_gromacs.sh` - Convenience shell script

Once production has run, the directory also holds `prod.tpr`, the compiled
run input. Analyses read it in preference to the PDB or GRO, because it
carries every atom, bond, mass and charge and has no atom limit.

The `.mdp` files and run script are not required to continue outside PolyzyMD;
you may replace them with your own GROMACS workflow. Use
`polyzymd run --engine gromacs` when you want PolyzyMD to perform the full local
build-and-run workflow.

---

## polyzymd analysis-topology

Write `system.prmtop` for runs built before PolyzyMD wrote it, from the
`solvated_system.pdb` and `system.xml` already in each run directory.

```bash
polyzymd analysis-topology RUN_DIR...
```

### Options

| Option | Description |
|--------|-------------|
| `--overwrite` | Rewrite `system.prmtop` where it already exists |

### Example

```bash
polyzymd analysis-topology /scratch/campaign/*/run_*
```

### Notes

- New builds write the file themselves; this command is for existing runs.
- The PDB is read with OpenMM's own reader, which accepts the hex serials it
  writes above 99,999 atoms, so systems of any size convert.
- A PDB and a `system.xml` with different particle counts are refused.
- Exit code 1 if any directory could not be converted; the others are still written.

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
| `--exclude` | - | No | from preset | SLURM `--exclude` override (e.g., "bgpu-g4-u20,bgpu-g4-u24"). Replaces the preset list; pass `""` to exclude nothing |
| `--pixi-env` | - | No | engine-specific | Runtime for generated Slurm jobs; OpenMM `auto` uses a fixed environment for known-site presets, and GROMACS uses `build` |
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

## polyzymd cancel

Stop self-resubmitting simulation chains, and hand them back when you want
to continue.

```bash
polyzymd cancel -c CONFIG [OPTIONS]
```

`scancel` on its own does not stop a chain. SLURM sends `SIGTERM`,
`run-segment` exits 99, and the job wrapper reads that as "interrupted, work
remains" and queues a successor within seconds. `polyzymd cancel` writes a
`STOP` marker into each replicate working directory first — the wrapper
refuses to submit a successor while it exists, and a successor that is
already queued exits before starting a segment — and then cancels the
matching queued and running jobs by job name.

### Options

| Option | Short | Required | Default | Description |
|--------|-------|----------|---------|-------------|
| `--config` | `-c` | Yes | - | Path to YAML configuration file |
| `--replicates` | `-r` | No | `1` | Replicate range (e.g. `1-5`, `1,3,5`) |
| `--scratch-dir` | - | No | from config | Scratch directory override; must match submission |
| `--resume` | - | No | false | Remove the `STOP` marker instead of writing it |
| `--stop-only` | - | No | false | Write the marker but leave running jobs alone |
| `--dry-run` | - | No | false | Report what would happen; write and cancel nothing |

### Example

```bash
# Stop three replicates now
polyzymd cancel -c config.yaml -r 1-3

# Let the current segment finish, then stop
polyzymd cancel -c config.yaml -r 1-3 --stop-only

# Allow the chains to run again, then resubmit
polyzymd cancel -c config.yaml -r 1-3 --resume
polyzymd submit -c config.yaml -r 1-3 --preset blanca-shirts
```

### The STOP marker

`<working_dir>/STOP` is a plain-text file that records who wrote it, on which
host, when, the config path and the replicate, and how to undo it. Deleting
the file by hand is equivalent to `--resume`.

The job wrapper also honours `POLYZYMD_STOP_CHAIN=1` in the job environment
and `POLYZYMD_STOP_FILE=<path>` to relocate the marker.

### Notes

- The marker is written before `scancel` runs, so a successor queued during
  cancellation still sees it.
- `--resume` only removes the marker; it does not resubmit. Use
  `polyzymd submit` afterwards.
- Outside a SLURM environment (no `scancel`) the marker is still written and
  the missing scheduler is reported as a warning.
- Already-running chains keep the job script that was rendered at submission
  time. A chain submitted before this feature existed does not check the
  marker; stop those with `scancel --batch --signal=KILL <job_id>`.

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
- Installs signal handlers before configuration loading or simulation setup
- Skips minimization/equilibration only after an atomic completed phase record
- Restarts an incomplete phase when no synchronized recovery record is valid

### Exit Codes

| Code | Meaning |
|------|---------|
| 0 | Segment completed successfully |
| 1 | Error |
| 99 | Graceful interruption (wall-time/preemption signal) |

### Notes

- This command is called by the generated SLURM scripts, not typically by users directly
- Progress is tracked in `progress.json` in the working directory
- Lifecycle state is tracked in per-phase `phase.json` files; orphan binary
  checkpoints never prove phase completion

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
| `-c, --config PATH` | One of `-c`/`--all` | Path to a YAML configuration file. Repeatable with `--format agent` or `json`. |
| `--all DIR` | One of `-c`/`--all` | Search `DIR` (up to 3 levels deep, hidden directories skipped) for `config.yaml` files. Repeatable. |
| `--format table\|agent\|json` | No | `table` (default) prints progress bars for one config. `agent` prints one compact line per replicate with SLURM state, throughput and ETA. `json` emits the same data as JSON. |
| `--no-slurm` | No | Skip the `squeue` query. Verdicts then rely on `progress.json` alone and cannot separate dead chains from running ones. |
| `--preset NAME` | No | Preset name to print in the resubmit hint for dead chains (`agent` format). |

### Agent format

`--format agent` is designed for scripts and LLM agents: it answers "is
anything still driving this replicate?" and "when will it finish?" in as few
characters as possible. It makes **one** `squeue -u $USER` call regardless of
how many configs and replicates it covers, joins that with each replicate's
`progress.json`, and for chains with no live job reads the newest SLURM log
for the line that explains the death.

```bash
polyzymd status --format agent --all /projects/me/sims --preset blanca-shirts
```

```
# polyzymd status  2026-09-11 16:28 UTC  2 system(s)  8 replicate(s): 5 running, 1 queued, 1 dead, 1 not_started

## CALB_ResorufinButyrate_none_1000ns_343K  (CALB/noPoly_CALB_water_343K/config.yaml)
run1   361.6/1000ns   36%  RUNNING      job 28248421 R 5:27 bgpu-shirts3  176ns/d  eta 3.6d
run4   354.2/1000ns   35%  QUEUED       job 28248423 PD ((Priority))  eta ?
run5   708.2/1000ns   71%  RUNNING      job 28248388 R 39:13 bgpu-shirts2  324ns/d  eta 22h

## RML_ResorufinButyrate_none_1000ns_333K  (RML/noPoly_RML_water_333K/config.yaml)
run3   393.4/1000ns   39%  DEAD         no job  last: FATAL: CUDA routing failed after 3 retries [RML_..._run3.28235208.out]
run4     0.0/1000ns    0%  NOT_STARTED  no job  last: Validation error: 354 polymer atom(s) lie within 1.00 A of the solute [build_r4_28214393.out]

# dead chains — resume from checkpoint with:
polyzymd submit -c RML/noPoly_RML_water_333K/config.yaml -r 3 --preset blanca-shirts
```

Verdict vocabulary (the fourth column) is fixed so callers can branch on it:

| Verdict | Meaning |
|---------|---------|
| `COMPLETED` | `progress.json` reports all production steps done |
| `RUNNING` | A SLURM job with this replicate's job name is in state `R` (or completing/configuring) |
| `QUEUED` | A matching job exists but is pending; the reason is shown in parentheses |
| `DEAD` | Work remains and no matching job is queued or running. Nothing will restart it. The `last:` field is the most informative error line near the end of the newest SLURM log, with the log filename in brackets. |
| `NOT_STARTED` | Directory exists but production never began (typically a failed build; the build log is consulted). |
| `NOT_FOUND` | Expected replicate directory is missing from scratch |

Throughput (`ns/d`) is measured from the newest segment with a known wall
window (a finished or interrupted segment), falling back to the live segment
timed from its start to now. Windows under 10 minutes or 1000 steps are
ignored. The ETA is remaining nanoseconds divided by that rate and is only
printed for `RUNNING` and `QUEUED` replicates.

If `squeue` is unavailable the header carries a warning and every non-complete
replicate falls back to `progress.json` alone; treat `DEAD` as unreliable in
that case.

### Table format

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
| `--pixi-env` | - | No | engine-specific | Runtime for the recovery job; OpenMM `auto` uses a fixed environment for known-site presets, and GROMACS uses `build` |
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
  --builtin                         Write a built-in analysis into the source tree
  --project-root DIRECTORY          Repository root for --builtin
  --force                           Overwrite existing files
  --dry-run                         Print paths without writing files
```

Inside a study, meaning the current directory holds `study.yaml` or sits below
one, the command writes `analyses/<NAME>.py` and `analyses/test_<NAME>.py` in
the study. The plugin holds a pydantic settings model and a class with a
`compute()` that returns observables. NAME may not be a built-in analysis or an
importable Python module, since the test imports the plugin by that name.

```bash
polyzymd new-analysis lid_opening
pytest analyses/test_lid_opening.py -q
```

Outside a study, or with `--builtin`, it writes a built-in analysis to
`src/polyzymd/analyses/<NAME>.py` and its tests to
`tests/analyses/plugins/test_<NAME>.py`.

```bash
polyzymd new-analysis solvent_shell --builtin
PYTHONPATH=$PWD/src pixi run -e test pytest tests/analyses/plugins/test_solvent_shell.py -q
```

---

(polyzymd-study-init)=
## polyzymd study init

Create a study folder: `study.yaml`, a `README.md` describing the layout, a
`.gitignore` for trajectories and checkpoints, and the folders `conditions/`,
`comparisons/`, `analyses/`, `structures/` and `workflows/`. See
{doc}`../how_to/study_layout`.

```bash
polyzymd study init -n NAME [--description TEXT]
```

(polyzymd-study-results)=
## polyzymd study results

Write every comparison's numbers in the study as `conditions.csv`,
`comparisons.csv` and `profiles.csv`. Run it anywhere inside a study.

```bash
polyzymd study results [-o DIRECTORY] [--analysis NAME ...]
```

`-o` defaults to `results/` in the study root. `--analysis` limits the tables to
the named analyses. The Python equivalent is
`polyzymd.analyses.load_results(study_root)`.

(polyzymd-study-export)=
## polyzymd study export

Package the study as a zip for publication, without its trajectories.

```bash
polyzymd study export [-o PATH.zip]
```

Only condition folders listed by a `comparison.yaml` are packaged. Files with
the suffixes `.dcd`, `.xtc`, `.trr`, `.nc`, `.chk` and `.cpt`, and the folders
`slurm_logs`, `.polymer_cache`, `.pixi`, `.git` and `__pycache__`, are left out.
The zip holds `bundle_manifest.json` with the SHA-256 of every file and the
size and fingerprint of every trajectory the results were computed from.

(polyzymd-study-verify)=
## polyzymd study verify

Check an unpacked study against its `bundle_manifest.json`.

```bash
polyzymd study verify [PATH]
```

Exits non-zero when a packaged file is missing or changed, or a downloaded
trajectory does not match the manifest. Trajectories not yet downloaded are
counted, not failed.

---

(cli-analyze)=
## polyzymd analyze

Run one analysis and print a validated result. One `-c` gives a per-condition
summary; two or more give pairwise comparisons with the first config as the
control. The command builds the comparison in memory, so no `comparison.yaml`
has to be written first.

### Usage

```bash
polyzymd analyze NAME -c config.yaml [-c other/config.yaml ...] [OPTIONS]
polyzymd analyze NAME -f comparison.yaml [OPTIONS]
```

List the available analysis names with `polyzymd compare run --list`.

### Options

| Option | Required | Description |
|--------|----------|-------------|
| `NAME` | Yes | Canonical analysis name, for example `rg`. |
| `-c, --config PATH` | One of `-c`/`-f` | Simulation `config.yaml`. Repeatable; the first one is the control. |
| `-f, --file PATH` | One of `-c`/`-f` | Existing `comparison.yaml` to analyze instead of `-c` configs. Cannot be combined with `-c` or `--set`. |
| `--replicates SPEC` | No | Replicates to analyze, for example `1-3`, `1,3,5` or `1-9:2`. Default: the replicate directories found on disk for each condition. |
| `--eq TEXT` | No | Equilibration window discarded from every replicate, for example `10ns`. Default: the comparison default, `10ns`. |
| `--label TEXT` | No | Condition label, one per `-c` in the same order. Default: the name of the directory holding the config. |
| `--run LABEL` | No | Run or pair label to report when the analysis measures one metric on several selections, for example `Protein` or `Polymer Oligomers` for rg. Default: the first one the plugin lists; the rest appear in `all_runs`. |
| `--set KEY=VALUE` | No | Plugin setting. Repeatable. The value is read as YAML, so `--set n_bins=50` gives an integer; a dotted key nests. |
| `--format agent\|json` | No | `agent` (default) prints at most 25 lines; `json` prints the full `ProtocolReport`. For a human-readable table of the same comparison, use `polyzymd compare run --format table`. |
| `-o, --output PATH` | No | Also write the rendered output to this file. |
| `--output-dir PATH` | No | Directory for `analysis/`, `comparison/` and `figures/`. Default: the current directory. |
| `--recompute` | No | Recompute replicates instead of reusing cached results. |

### Agent format

`--format agent` is the default and is designed for scripts and LLM agents. It
prints a header, one line per condition, one line per comparison, any warnings
and one `verdict:` line per comparison, in at most 25 lines with no borders and
no colour.

```bash
polyzymd analyze rg -c noPoly/config.yaml -c SBMA50/config.yaml --eq 10ns
```

```
# polyzymd analyze rg  metric mean_rg  unit A  eq 10ns  conditions 2  replicates 3,3  protocol rg/1
noPoly  n 3  mean 18.42  sem 0.05  ci95 18.2 to 18.64  values 18.4, 18.5, 18.36
SBMA50  n 3  mean 18.73  sem 0.06  ci95 18.47 to 18.99  values 18.71, 18.8, 18.68
noPoly vs SBMA50  delta +0.31  ci95 0.02 to 0.6  p 0.041  p_adj 0.041  test student_t  correction BH  d 1.9  significant
verdict: SBMA50 larger mean_rg than noPoly (delta +0.31 A, 95% CI 0.02 to 0.6, p_adj 0.041, n 3 vs 3)
```

Line shapes:

| Line | Fields |
|---|---|
| header | `# polyzymd analyze <analysis>  metric <key>  unit <unit or none>[  run <label>]  eq <window>  conditions <count>  replicates <n,n,...>  protocol <analysis>/<protocol_version>` |
| condition | `<label>  n <count>  mean <value>  sem <value>  ci95 <low> to <high>  values <per-replicate values>` |
| comparison | `<a> vs <b>  delta <signed>  ci95 <low> to <high>  p <value>  p_adj <value>  test <name>  correction <name>  d <value>  significant\|not_significant\|no_test\|not_testable` |
| warning | `warning: <text>` |
| verdict | `verdict: <sentence>` |

`na` stands for a number that does not exist, such as the standard error of a
single replicate. When a report does not fit in 25 lines, condition and
comparison lines are dropped and the last line says how many.

The verdict vocabulary is fixed so a caller can branch on it:

| Word | Meaning |
|---|---|
| `larger` | The second condition differs from the control after correction and the difference is positive |
| `smaller` | The second condition differs from the control after correction and the difference is negative |
| `no significant difference` | The test ran and the adjusted p value did not clear alpha |
| `no test recorded` | The plugin stored no multiplicity-corrected p value, so the comparison describes a difference without deciding it |
| `changed` | The difference is significant but the two means are equal at the stored precision |
| `not testable` | A condition has fewer than two replicates, so the test is undefined |

`--format json` prints the full report. Every field is documented in
{doc}`analysis_protocol_report`.

### Exit codes

| Code | Meaning |
|------|---------|
| 0 | The analysis ran and the report was printed |
| 2 | A typed analysis error: unknown analysis name, missing config, bad `--replicates` or `--set`, conflicting `-c` and `-f`, or a pipeline failure. The message is printed on one line prefixed `error:` and the fix on the next prefixed `fix:`, both on stderr |

### Example

```bash
# Single condition
polyzymd analyze rg -c enzyme_water/config.yaml --eq 10ns

# Comparison with explicit labels and replicates
polyzymd analyze rmsf -c A/config.yaml -c B/config.yaml \
  --label "no polymer" --label "50% SBMA" --replicates 1-3

# Full JSON record saved to a file
polyzymd analyze sasa -c A/config.yaml -c B/config.yaml \
  --format json -o sasa_report.json

# An existing comparison project
polyzymd analyze rmsf -f comparison.yaml

# Report the polymer selection instead of the protein
polyzymd analyze rg -c A/config.yaml -c B/config.yaml --run "Polymer Oligomers"
```

### Notes

- Run through `pixi run -e analysis`; the default pixi environment has no
  `polyzymd`.
- The replicate is the sampling unit. Every interval and every test uses the
  replicate count as its sample size.
- The run writes and reuses the same cached artifacts as
  `polyzymd compare run`, so the two commands share results.

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
  --format TEXT           Output format: table, markdown, json, agent [default: table]
  -o, --output PATH      Save formatted output to file
  -q, --quiet            Suppress INFO messages
  --debug                Enable DEBUG logging
  --list                 List available comparison types and exit
```

Without `--recompute`, a cached `run_<replicate>/result.json` is reused only
when the trajectory files it records still have the same size and modification
time and the settings fingerprint matches; otherwise that replicate is
recomputed. See the cache reuse section of the comparison reference for the
rules that apply to the other commands.

#### Example

```bash
# Run RMSF comparison (uses plugins.rmsf from comparison.yaml)
polyzymd compare run rmsf

# Override equilibration time
polyzymd compare run rmsf --eq-time 20ns

# Run contacts comparison with markdown output
polyzymd compare run contacts --format markdown -o report.md

# Print the compact agent report instead of the plugin's table
polyzymd compare run rg --format agent

# List all available analysis types
polyzymd compare run --list
```

`--format agent` renders the finished comparison through the same
`ProtocolReport` renderer as {ref}`polyzymd analyze <cli-analyze>`: a header,
one line per condition, one line per comparison, any warnings and one
`verdict:` line per comparison, in at most 25 lines. Use it when a script or an
agent has to read the result; use `--format json` when it needs the full
comparison artifact, which carries more per-plugin detail than the report does.
The line shapes and the verdict vocabulary are documented under
{ref}`polyzymd analyze <cli-analyze>`, and the report fields in
{doc}`analysis_protocol_report`.

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
  --recompute            Regenerate comparison and plot outputs
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
| 2 | Typed analysis error from {ref}`polyzymd analyze <cli-analyze>`; the message and the fix are printed on stderr, one line each |
| 99 | Graceful shutdown — simulation was interrupted but interrupted state was saved (see {doc}`../how_to/hpc_slurm`) |

---

## See Also

- {doc}`../get_started/quickstart` - Getting started tutorial
- {doc}`configuration` - Configuration file reference
- {doc}`../how_to/hpc_slurm` - HPC and SLURM guide
- {doc}`../how_to/analysis_rmsf_quickstart` - RMSF analysis tutorial
- {doc}`../how_to/analysis_compare_conditions` - Comparing simulation conditions
- {doc}`../how_to/analysis_agent_protocol` - Getting a validated number with one command
- {doc}`analysis_protocol_report` - The `ProtocolReport` schema
