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

Global options go *before* the subcommand: `polyzymd --no-color status -c
config.yaml`, not `polyzymd status --no-color -c config.yaml`.

### Colored Output

The color of an INFO or DEBUG log line shows which part of PolyzyMD wrote
it. WARNING lines are always yellow and ERROR lines are always red. The
message text carries the same information, so a log without color loses
nothing.

| Part of PolyzyMD | Messages | Color (16-color fallback) |
|---|---|---|
| `polyzymd.cli`, `polyzymd.workflow` | Commands and job workflow | Lavender (bright blue) |
| `polyzymd.builders` | System build | Sage green (bright green) |
| `polyzymd.simulation.runner`, `.continuation` | OpenMM setup, continuation and recovery | Warm peach (bright magenta) |
| `polyzymd.simulation.progress`, `.signals` | Progress and signal handling | Steel blue (bright cyan) |
| `polyzymd.core` | Shared data structures | Parchment (dark yellow) |
| `polyzymd.data` | Bundled data | Aqua (dark cyan) |
| `polyzymd.exporters` | Export to other formats | Lilac (dark magenta) |
| `polyzymd.utils` | Other utilities | Gray (white) |

PolyzyMD reads the terminal to choose the color depth:

| Depth | When |
|---|---|
| 24-bit color | `COLORTERM` is `truecolor` or `24bit` |
| 256 colors | `TERM` contains `256color` |
| 16 colors | Any other terminal, for example `TERM=xterm-16color` on many cluster login nodes |
| No color | Standard error is not a terminal (a pipe or a file), `TERM=dumb`, `NO_COLOR` is set, or `--no-color` is given |

To turn color off, use one of these:

```bash
polyzymd --no-color build -c config.yaml   # global option, before the command
NO_COLOR=1 polyzymd build -c config.yaml   # any non-empty value (no-color.org)
```

Some messages, such as "Build complete!" and error exits, use Click's green
and red styles instead of the colors above.

### Logging Behavior

By default, PolyzyMD suppresses verbose log messages from OpenFF Interchange and Toolkit libraries. These libraries generate per-atom INFO messages during system building (e.g., "Preset charges applied to atom index 8667" or "Key collision with different parameters"). For large systems with tens of thousands of atoms, this can produce millions of log lines.

**Default behavior:** OpenFF INFO logs are suppressed; only WARNING and ERROR messages are shown.

**To enable OpenFF logs for debugging:**

```bash
polyzymd --openff-logs build -c config.yaml
polyzymd --openff-logs run -c config.yaml --engine gromacs
```

Use `--openff-logs` to investigate charge assignment problems and build
failures.

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

- YAML syntax validity: the file is a mapping and gives each key once
- Required fields are present, numbers are finite and not `true`/`false`
- The enzyme PDB has atoms, the substrate SDF holds the conformer, SMILES parse
  and the charge method can charge each molecule (errors, exit 1)
- Cached polymer SDFs and reaction templates are reported as warnings when missing
- Settings a new build would not run as written: `NVE`, Nose-Hoover or Andersen
  on OpenMM, `samples` above the MD steps, unset `$VAR` in paths (errors)
- Monomer probabilities sum to 1.0
- Valid enum values (water model, ensemble, etc.)
- Co-solvent specification (mole_fraction XOR concentration)

### Example

```bash
cd examples/quickstart
polyzymd validate -c config.yaml
```

**Output (success):**
```
Validating configuration: config.yaml
Configuration is valid!


Summary:
  Name: trpcage_water
  Engine: openmm
  Enzyme: trpcage
  Substrate: None (apo simulation)
  Polymers: Disabled
  Co-solvents: none
  Temperature: 300.0 K
  Pressure: 1.0 atm

Simulation phases:
  Equilibration: 0.002000 ns across 1 stage(s)
    - equil: 0.002 ns (NVT)
  Production: 0.004 ns (NPT)
```

With polymers, the summary also lists `Count`, `Length` and one
`Monomer <label>: <percent>` line per monomer. When a file that the config
names is missing, such as the polymer SDF directory, a `Referenced file
warnings:` block with one `Warning:` line per file prints between
`Configuration is valid!` and `Summary:`.

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
| `--format` | - | No | the config's `engine` | Export format (`gromacs`). Without it, PolyzyMD writes OpenMM files, or GROMACS files when the config's `engine` is `gromacs` |

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
`<replicate folder>/gromacs/`. For `<prefix>` and `<stage>`, see
{ref}`gromacs-output-files`. The core handoff files are:

- `<prefix>.gro` - GROMACS coordinate file
- `<prefix>.top` - GROMACS topology file
- `*.itp` - Molecule parameter files (one per component)
- Position restraints (`#ifdef POSRES_PROTEIN`, etc.) appended into molecule `.itp` files

PolyzyMD may also generate convenience defaults:

- `em.mdp` - Energy minimization parameters
- `eq_01_<stage>.mdp` - Equilibration stage parameters, one file per stage
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

Builds the system and runs it with the config's engine, or with the engine
that `--engine` names:

- `gromacs` exports GROMACS files and runs the full GROMACS workflow
- `openmm` builds and runs the OpenMM simulation locally

```bash
polyzymd run -c <path> [--engine gromacs|openmm] [options]
polyzymd run -c <path> --engine gromacs --gmx-path /usr/local/gromacs/bin/gmx
polyzymd run -c <path> --engine openmm --dry-run
```

### Options

| Option | Short | Required | Default | Description |
|--------|-------|----------|---------|-------------|
| `--config` | `-c` | Yes | - | Path to YAML configuration file |
| `--replicates` | `-r` | No | "1" | Replicate range (for example "1", "1-3", "1,3,5") |
| `--engine` | - | No | the config's `engine` | Local engine to run: `gromacs` or `openmm` |
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
2. Build system (enzyme + substrate + polymers + solvent). If an earlier
   `polyzymd build` wrote a build for this config, reuse it. The log says which
   build is used.
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
- For OpenMM, use `polyzymd run` with `engine: openmm` in the config (or
  `--engine openmm`) to run on this machine, or `polyzymd submit` to submit self-resubmitting SLURM jobs

### Output Files

A GROMACS run writes its files to `<replicate folder>/gromacs/`. For the list
of files, see {ref}`gromacs-output-files`.

---

## polyzymd submit

Submit self-resubmitting simulation jobs to SLURM for HPC execution.

Each replicate gets one SLURM script that handles the simulation lifecycle:
equilibration, production segments, interruption recovery, and resubmission.
See {doc}`../how_to/hpc_slurm` for details.

`submit` does not build. Run `polyzymd build` for each replicate first, in a
compute job. Before it writes any script, `submit` checks each replicate: an
OpenMM build must match its `build_manifest.json` and the config, and a GROMACS
build must have the `.top`, `.gro`, `em.mdp` and `prod.mdp` files. If a
replicate fails the check, `submit` stops with an error that gives the
`polyzymd build` command. `--generate-only` makes the same check.
`--dry-run` does not.

`submit` runs `sbatch --export=NONE`, so a job does not inherit the shell
environment of the submitting host. GROMACS `module_load` runs in the job only.

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
| `--time-limit` | - | No | from preset | Override SLURM time limit: minutes, M:SS, H:MM:SS or D-H:MM:SS; other values are refused |
| `--memory` | - | No | 3G | Override SLURM memory allocation |
| `--account` | - | No | - | Override SLURM account / allocation ID |
| `--partition` | - | No | from preset | Override SLURM partition |
| `--qos` | - | No | - | Override SLURM QoS |
| `--gpu-type` | - | No | - | GPU type for GRES (e.g., "a100", "a40", "mi100") |
| `--constraint` | - | No | - | SLURM `--constraint` for node features (e.g., "A40", "A40\|A100") |
| `--nodelist` | - | No | - | SLURM `--nodelist` override (e.g., "gpu-node-001") |
| `--exclude` | - | No | from preset | SLURM `--exclude` override (e.g., "node01,node02"). Replaces the preset list; pass `""` to exclude nothing |
| `--pixi-env` | - | No | engine-specific | Runtime for generated Slurm jobs; OpenMM `auto` uses a fixed environment for known-site presets, and GROMACS uses `build` |
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

A preset sets the partition, QoS, account and time limit of each job. For the
values of each preset, see the preset table in {doc}`../how_to/hpc_slurm`.

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

# GROMACS GPU submission with constraint (after polyzymd build ... --format gromacs)
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
queued and running jobs that work in that directory. A replicate whose
working directory does not exist and that has no queued job has nothing to
cancel: the command says so and writes nothing.

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

A GROMACS job works in `<working_dir>/gromacs` and checks for `STOP` there
and in `<working_dir>`. For a GROMACS config, `cancel` also writes
`<working_dir>/gromacs/STOP` when that folder exists, so that a job script
written before this check stops too. `--resume` removes both files.

The job wrapper also honours `POLYZYMD_STOP_CHAIN=1` in the job environment
and `POLYZYMD_STOP_FILE=<path>` to relocate the marker.

### Notes

- The marker is written before `scancel` runs, so a successor queued during
  cancellation still sees it.
- `--resume` only removes the marker; it does not resubmit. Use
  `polyzymd submit` afterwards.
- Outside a SLURM environment (no `squeue` or `scancel`) the marker is still
  written. `cancel` warns `squeue not found — skipping duplicate-job check`
  and then prints `no queued or running job` for each replicate. That line
  only means that no job could be listed; check the queue on the cluster.
- A running chain keeps the job script that `submit` wrote. OpenMM and
  GROMACS job scripts check the marker. Job scripts of PolyzyMD 1.2 and
  earlier do not. Stop those chains with
  `scancel --batch --signal=KILL <job_id>`.

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
| `--allow-report-interval-change` | - | No | false | Continue even if the config now gives a different number of steps between trajectory frames than earlier segments used |

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
- When the run has a `progress.json`, it updates it from the files on disk. It
  never creates one, so a run without one (a legacy or downloaded run) is left
  as it is
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
| `--unfinished` | No | `agent` and `json` only. Omit completed replicates; a fully completed system collapses to its header line. |

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

Each system header ends with `done/total completed`, which answers "does every condition have enough finished replicates" without reading the rows. When a `DEAD` replicate's newest log holds no error line, as after a node failure or power loss, its `last:` field shows the SLURM end state from one `sacct` call, for example `slurm: NODE_FAIL exit 0:0`.

Verdict vocabulary (the fourth column) is fixed so callers can branch on it:

| Verdict | Meaning |
|---------|---------|
| `COMPLETED` | `progress.json` reports all production steps done |
| `RUNNING` | A SLURM job that works in this replicate's run directory is in state `R` (or completing/configuring) |
| `QUEUED` | A matching job exists but is pending; the reason is shown in parentheses |
| `DEAD` | Work remains and no matching job is queued or running. Nothing will restart it. The `last:` field is the most informative error line near the end of the newest SLURM log whose `Work dir:` line names this run directory, with the log filename in brackets. |
| `NOT_STARTED` | Directory exists but production never began (typically a failed build; the build log is consulted). |
| `NOT_FOUND` | Expected replicate directory is missing from scratch |
| `CORRUPT` | `progress.json` exists but cannot be read (empty, truncated or of the wrong shape). The `last:` field says why. `status` leaves the file as it is; fix or remove it by hand. |

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

- This is a **read-only** command: it reads `progress.json` files and writes
  nothing. A run without `progress.json` is shown from a scan of its files
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
| `--dry-run` | - | No | false | With `--submit`, print the script path and SLURM settings; write and submit nothing |
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
  `{working_dir}/recovery_scripts/` and submits it. The chain's own script in
  `daisy_chain_scripts/` is left as it was. The job uses `--preset` (default
  `aa100`), so pass the preset the chain was submitted with
- `recover` never writes `progress.json`; the recovery job's runner does
- If `progress.json` cannot be read, `recover` (and `check-progress` in a
  running chain) stops with an error and leaves the file as it is
- The recovery job is identical to a normal submission job — it uses `run-segment`
  to determine what work remains and continues from there

---

## polyzymd clean-pdb

Replace nonstandard residues with their standard residues, and add the
missing hydrogens, with PDBFixer.

```bash
polyzymd clean-pdb -i <input.pdb> [-o <output.pdb>] [--ph 7.4]
```

### Options

| Option | Short | Required | Default | Description |
|--------|-------|----------|---------|-------------|
| `--input` | `-i` | Yes | - | The input PDB file |
| `--output` | `-o` | No | `<input name>_clean.pdb` | The cleaned PDB file |
| `--ph` | - | No | `7.4` | The pH for the protonation states of the added hydrogens |

### Notes

- The output keeps the chain IDs and residue numbers of the input.
- The command does not remove waters or other molecules, does not select one
  copy of a protein, does not set chain IDs, and does not add missing residues
  or heavy atoms. See {doc}`../tutorials/prepare_pdb_for_openff`.
- PDBFixer places the hydrogens with OpenMM on the CPU platform. To use
  another platform, set `OPENMM_DEFAULT_PLATFORM`, for example to `CUDA`.
- Run it in the `build` environment, which holds PDBFixer.

---

## polyzymd info

Display PolyzyMD installation and dependency information.

```bash
polyzymd info
```

### Example Output

After a banner with the PolyzyMD logo:

```
PolyzyMD - Molecular Dynamics of Proteins with Polymers, Ligands and Co-solvents
Version: 1.3.0

Dependencies:
  OpenMM: <version>
  OpenFF Toolkit: <version>
  OpenFF Interchange: <version>
  Pydantic: <version>
  packmol: /path/to/env/bin/packmol
  gmx: NOT ON PATH

Example configs: polyzymd/templates/examples/
```

`packmol` and `gmx` show the executable found on `PATH`, or `NOT ON PATH`.

### Use Cases

- Verify installation is complete
- Check dependency versions for troubleshooting
- Confirm GPU-enabled OpenMM is installed

---

(cli-analyze)=
## polyzymd analyze

Run one analysis and print a validated result. One `-c` gives a per-condition
summary; two or more give pairwise comparisons with the first config as the
control. The command reads the simulation configs given with `-c`.

### Usage

```bash
polyzymd analyze NAME -c config.yaml [-c other/config.yaml ...] [OPTIONS]
polyzymd analyze [RUN] --study study.yaml [OPTIONS]
```

With `--study`, the conditions, equilibration window, stride, replicates and
settings come from the study file, and `RUN` is one of its `analyses:`
entries, a shipped analysis or the study's own function; with no `RUN`,
every entry runs in turn. See {doc}`../how_to/study_yaml`.

`NAME` is one of `rg`, `rmsd`, `rmsf`, `rmsd_per_residue`, `sasa`,
`secondary_structure`, `contacts`, `native_contacts`, `hydrogen_bonds` and
`distances`. `rmsd` is one value per frame, the root mean square over atoms
of their distance from the reference; `rmsd_per_residue` is one value per residue,
the root mean square over frames of each atom's distance from the reference,
as `gmx rmsf -od` reports. An unknown name is refused with the list of all of them.

### Options

| Option | Required | Description |
|--------|----------|-------------|
| `NAME` | Yes, without `--study` | Canonical analysis name, for example `rg`; with `--study`, a run of the study file, and every run when left out. |
| `-c, --config PATH` | Yes, without `--study` | Simulation `config.yaml`. Repeatable; the first one is the control. Refused with `--study`. |
| `--data DIR` | No | Directory holding the run directories of every condition, for this command only, in place of each config's `scratch_directory` and of a study's `data.local.yaml`. The config hash is that of the config as written. |
| `--project PATH` | No | `project.yaml`, or the project folder. Runs `RUN` (or every analysis) in each study of the project that runs it, as `--study` would, one study after another; a study that fails is reported and the next runs, then the command exits 2. Refused with `--study`, `-c` or `--output-dir`. |
| `--study PATH` | No | `study.yaml`, or the folder holding it. Gives the conditions, `--eq`, `--stride`, `--replicates` and the run's settings, and stores the run in `<study>/results/<run>/` with its `report.json`. Options given on the command line override the file; `--label` then picks conditions of the study, and that report is printed but not saved as the run's `report.json`. |
| `--replicates SPEC` | No | Replicates to analyze, for example `1-3`, `1,3,5` or `1-9:2`. Default: the replicate directories found on disk for each condition. |
| `--eq TEXT` | No | Equilibration window discarded from every replicate, for example `10ns`. Default `10ns`. |
| `--label TEXT` | No | Condition label, one per `-c` in the same order. Default: the name of the directory holding the config. |
| `--list` | No | Print every shipped analysis, what it measures, and its settings with their defaults, then exit. |
| `--run LABEL` | No | Run or pair label to report when the analysis measures one metric on several selections, for example one atom pair for distances. Default: the first one the analysis lists; the rest appear in `all_runs`. |
| `--set KEY=VALUE` | No | Top-level analysis setting. Repeatable. The value is read as YAML, so `--set threshold=3.0` gives a number and `--set groups='{protein: chainid A, polymer: chainid C}'` a mapping. A dotted key such as `groups.protein` is refused; give the whole top-level setting as a mapping instead. |
| `--format agent\|json` | No | `agent` (default) prints one line per condition and comparison; `json` prints the full `ProtocolReport`. |
| `-o, --output PATH` | No | Also write the rendered output to this file. |
| `--output-dir PATH` | No | Directory for `polyzymd_results/`, where the measured values are stored, and `figures/`. Default: the current directory. |
| `--until TIME` | No | End of a common analysis window, for example `38ns`: production frames after it are left out for every condition, so conditions simulated for different lengths are compared over the same time. It is recorded with each stored result. A report comparing conditions whose analysed production ends more than 10% apart warns and suggests it. |
| `--stride N` | No | Measure every `N`-th production frame of every replicate, starting with the first after the window. Default `1`. The report header then shows `stride N`. |
| `--recompute` | No | Recompute replicates instead of reusing cached results. |
| `--no-plots` | No | Do not draw figures. By default rg, rmsd, rmsf, rmsd_per_residue, distances, sasa, secondary_structure, contacts, native_contacts and hydrogen_bonds draw theirs to `<output-dir>/figures/<analysis>/`, and the folder is recorded under `output_paths.figures` in the JSON report. |
| `--no-eq-check` | No | Skip the pymbar equilibration diagnostic, which rg, rmsd, native_contacts, distances and the totals of sasa report for each replicate's per-frame series. Values and statistics are the same either way. |
| `--submit` | No | Submit to SLURM instead of running here: one array task per condition and replicate, then a report job; see [SLURM submission](#cli-analyze-submit). |
| `--dry-run` | No | With `--submit`, or alone, write the SLURM scripts and print the two `sbatch` commands without submitting. |
| `--preset NAME` | No | SLURM partition, QoS and account of a cluster for `--submit`: `alpine-cpu`, `blanca-shirts`, `blanca-chbe-rdi` or `bridges2-rm`. |
| `--partition TEXT` | No | SLURM partition for `--submit`, replacing the preset's. |
| `--account TEXT` | No | SLURM account for `--submit`, replacing the preset's. |
| `--qos TEXT` | No | SLURM QoS for `--submit`, replacing the preset's. |
| `--time TEXT` | No | Time limit of each job for `--submit`. Default `12:00:00`. |
| `--mem TEXT` | No | Memory of each job for `--submit`. Default `16G`. |
| `--cpus N` | No | CPUs of each job for `--submit`. Default `2`. |

(cli-analyze-submit)=
### SLURM submission

`--submit` runs the same analysis as a SLURM array instead of in the current
process:

```bash
polyzymd analyze hydrogen_bonds -c A/config.yaml -c B/config.yaml --eq 10ns \
  --submit --preset blanca-shirts
```

It writes `<output-dir>/slurm/<analysis>_<YYYYmmdd-HHMMSS>/` holding
`tasks.tsv` (one config, label and replicate per line), `replicates.sbatch`
(the array: each task measures one replicate with `--no-plots` and stores it
under `polyzymd_results/`), `report.sbatch` (the full command, which reads
every stored result, measures any replicate a task left unmeasured, draws the
figures and writes the report) and `logs/`. The report goes to `-o PATH`, or
by default to `report.txt` (`report.json` with `--format json`) in that
folder. The report job is submitted with `--dependency=afterany:<array>`, so it
starts once every task has ended, whether or not each succeeded.

Both jobs run the Python interpreter and `PYTHONPATH` of the submitting
process, in the submitting directory, so they measure with the same PolyzyMD
and resolve relative paths as the command did. `--dry-run` writes the folder
and prints the two `sbatch` commands without submitting. Every name
`polyzymd analyze` runs can be submitted; see {doc}`../how_to/hpc_execution`.

### Agent format

`--format agent` is the default and is designed for scripts and LLM agents. It
prints a header, one line per condition, one line per comparison, any warnings
and one `verdict:` line per comparison, with no borders and no colour.

```bash
polyzymd analyze rg -c noPoly/config.yaml -c SBMA50/config.yaml --eq 10ns
```

```
# polyzymd analyze rg  metric mean_rg  unit A  eq 10ns  conditions 2  replicates 3,3  protocol rg/2
noPoly  n 3  mean 18.42  sem 0.05  ci95 18.2 to 18.64  values 18.4, 18.5, 18.36  replicates 1,2,3  g 473.7  n_eff 19
SBMA50  n 3  mean 18.73  sem 0.06  ci95 18.47 to 18.99  values 18.71, 18.8, 18.68  replicates 1,2,3  g 402.1  n_eff 22
noPoly vs SBMA50  delta +0.31  ci95 0.02 to 0.6  p 0.041  p_adj 0.041  test welch_t  correction BH  d 1.9  significant
verdict: SBMA50 larger mean_rg than noPoly (delta +0.31 A, 95% CI 0.02 to 0.6, p_adj 0.041, n 3 vs 3)
```

Line shapes:

| Line | Fields |
|---|---|
| header | `# polyzymd analyze <analysis>  metric <key>  unit <unit or none>[  run <label>]  eq <window>  conditions <count>  replicates <n,n,...>  protocol <analysis>/<protocol_version>` |
| condition | `<label>  n <count>  mean <value>  sem <value>  ci95 <low> to <high>  values <per-replicate values>[  replicates <numbers>  g <value>  n_eff <value>  eq_detected <ns>]`; the bracketed fields appear when the per-frame series is stored; `g` and `n_eff` are the statistical inefficiency and effective sample size of each replicate's series, and `eq_detected` is the latest start of an equilibrated region that pymbar detects among the replicates (see {doc}`../explanation/convergence_detection`) |
| comparison | `<a> vs <b>  delta <signed>  ci95 <low> to <high>  p <value>  p_adj <value>  test <name>  correction <name>[  family <m>]  d <value>  significant\|not_significant\|no_test\|not_testable` |
| per-label comparison | For a labelled result such as a per-residue profile, one line per compared condition in place of one line per label: `<a> vs <b>  labels <count>  tested <count>  family <m>  test <name>  correction <name>  lower <count>  higher <count>`, where `family` counts every tested label of every compared condition |
| per-label hits | `<a> vs <b>  lower: <label> delta <signed> p_adj <value>, ...` and the same with `higher:`, listing every label that differs significantly; printed only when there is one |
| note | `note: the <n> per-label condition rows and <m> per-label comparison rows are in the JSON report` |
| warning | `warning: <text>` |
| verdict | `verdict: <sentence>` |

`na` stands for a number that does not exist, such as the standard error of a
single replicate. Every condition and comparison is printed, however many
there are.

The verdict vocabulary is fixed so a caller can branch on it:

| Word | Meaning |
|---|---|
| `larger` | The second condition differs from the control after correction and the difference is positive |
| `smaller` | The second condition differs from the control after correction and the difference is negative |
| `no significant difference` | The test ran and the adjusted p value did not clear alpha |
| `no test recorded` | The comparison stored no multiplicity-corrected p value, so the comparison describes a difference without deciding it |
| `changed` | The difference is significant but the two means are equal at the stored precision |
| `not testable` | A condition has fewer than two replicates, or every replicate of both conditions has the same value, so the test is undefined |

`--format json` prints the full report. Every field is documented in
{doc}`analysis_protocol_report`.

### Exit codes

| Code | Meaning |
|------|---------|
| 0 | The analysis ran and the report was printed |
| 2 | A typed analysis error: unknown analysis name, missing config, bad `--replicates` or `--set`, an invalid study file, or a pipeline failure; with `--study` and no `RUN`, any run that failed. The message is printed on one line prefixed `error:` and the fix on the next prefixed `fix:`, both on stderr |

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

# Fraction of frames the polymer buries each residue, compared residue by residue
polyzymd analyze contacts -c A/config.yaml -c B/config.yaml --stride 10 --run contact_fraction_residues

# Hydrogen bonds between the protein and the polymer, compared residue pair by residue pair
polyzymd analyze hydrogen_bonds -c A/config.yaml -c B/config.yaml --eq 10ns --run protein_polymer_pairs

# Measure the protein backbone instead of all protein atoms
polyzymd analyze rg -c A/config.yaml -c B/config.yaml --set selection='protein and name CA'
```

### Notes

- Run through `pixi run -e analysis`; the default pixi environment has no
  `polyzymd`.
- The replicate is the sampling unit. Every interval and every test uses the
  replicate count as its sample size.
- Every analysis runs through the study API and stores what it measures,
  per condition and replicate, under `polyzymd_results/` in the output
  directory, so a second run with the same inputs and settings reads the
  stored values back; `--recompute` measures again.

---

## polyzymd project

Commands on a project folder: the studies of one paper; see
{doc}`../how_to/project` and {doc}`../explanation/projects`.

### polyzymd project check

```bash
polyzymd project check [PATH] [--production]
```

Reads `project.yaml` (`PATH` is the file or its folder, by default the
current directory) and each study's `study.yaml` without loading any
trajectory. Prints `analysis <run>: studies <labels>` for each project
analysis, `stats <function>: up to date|stale: ...|not run` when it names a
plan, then `== study <label>` and the `polyzymd study check` of each study,
then `== project` with the git line and the project's metadata line, once.
A study's metadata line is printed only when the study has its own metadata.
`--production` adds each condition's production length to each study's
lines, as `study check --production` does. It reads every run's trajectory.
Exits 2 when the project file or a study cannot be read.

### polyzymd project init

```bash
polyzymd project init PATH --study LABEL [--study ...] [--holder NAME] [--no-git]
```

Writes a project folder at `PATH`: `project.yaml` listing the studies,
`analyses/`, `stats/` and `figures/`, licences, a README, and one study folder
per `--study` with a `study.yaml` to fill in. Labels are also folder names
(lower case, digits, `_`). Without `--no-git` the project becomes a git
repository with one commit. To move existing studies in, see
{doc}`../how_to/move_studies_into_project`.

### polyzymd project add-study

```bash
polyzymd project add-study LABEL [--project PATH]
```

Writes the study folder `LABEL` in the project at `PATH` (the
`project.yaml` or its folder, by default the current directory), with a
`study.yaml` to fill in, as `project init --study` does. Adds one line under
`studies:` in `project.yaml` and keeps the rest of the file. Commits nothing.
Exits 2 when the label is not a folder name (lower case, digits, `_`), or the
project already lists the label or holds its folder.

### polyzymd project freeze

```bash
polyzymd project freeze [PATH] [--tag TAG]
```

Freezes every study of the project, then writes the project's
`manifest.json`, `CITATION.cff` and `.zenodo.json`, commits the generated
files and every study's results, tags the project (`project-v<n>` by
default), and lays out `deposit/` with `deposit/UPLOAD.md`. Every gap is a
warning, prefixed with the study's label when it is a study's. Exits 2 only
when a file cannot be read or the tag exists.

---

## polyzymd study

Commands on a study folder; see {doc}`../how_to/study_yaml` and
{doc}`../explanation/study_folders`.

### polyzymd study check

```bash
polyzymd study check [PATH] [--production]
```

`--production` adds each condition's production length, read from every run's
trajectory headers and segments (minutes for long restarted chains); without
it no trajectory is read.

Reads `study.yaml` (`PATH` is the file or its folder, by default the current
directory) without loading any trajectory, and prints:

| Line | Fields |
|---|---|
| header | `study <path>  equilibration <window>  stride <n>[  replicates <list>]` |
| system | `system: <description>`, when the study has one |
| project | `project <project.yaml> as study <label>`, when a project lists the study |
| names | `structure <name>: <path>` and `region <name>: <selection>`, one line each |
| condition | `control\|condition <label>: replicates <numbers> under <directory> (from data.local.yaml\|config)`, with `; production <ns>` (a range when the replicates differ) under `--production`, or `no replicates found under <directory> ...` |
| analysis | `analysis <run>: <settings>; stored results in <folder>[ with its report]`, or `no stored results`; for the study's own function, `analysis <run> (<file>:<function>, <kind>)`, after importing it |
| git | `git: commit <sha>; inputs committed`, `git: commit <sha>; <n> uncommitted inputs: <first five paths> and <m> more`, or `git: not a repository` |
| metadata | `metadata (<file>): complete`, or `metadata (<file>): <n> gaps for publishing: <gap>; <gap> ...` |
| publish | `publish: follow <study>/deposit/UPLOAD.md` after a freeze, or `publish: when the analyses are final, run polyzymd study freeze` |
| reproduce | In a downloaded frozen study (a `manifest.json` but no `deposit/`): how to point it at the trajectories and rerun or redraw |
| citation | `cite: <how to cite PolyzyMD>` |

`analyze` and the `study` commands print only reports and warnings on the
console; the full log, with library messages, goes to `logs/polyzymd-<command>-<time>.log`
(in the study folder, or the output directory), whose path is printed first
as `log: <path>`. `polyzymd -v` keeps INFO on the console.

It exits 2 when the study file, a condition's config or a listed function
cannot be read, and 0 otherwise; missing runs are not errors.

### polyzymd study init

```bash
polyzymd study init DIRECTORY [--condition LABEL=CONFIG ...] [--new-condition LABEL ...] \
  [--equilibration 100ns] [--holder NAME] [--no-git]
```

Creates a study folder: `study.yaml`, `conditions/`, `structures/`,
`analyses/`, `figures/`, `results/`, `environment/`, `README.md`,
`LICENSE-data` (CC-BY-4.0), `LICENSE-code` (MIT), `data.example.yaml` and
`.gitignore`, then makes it a git repository with one commit.

| Option | Description |
|---|---|
| `--condition LABEL=CONFIG` | Copies an existing `config.yaml` to `conditions/<label>/config.yaml`, with every input file it names copied to `conditions/<label>/structures/` and its path made relative. The copy writes new runs into `runs/<label>/` of the study; where the existing runs are goes into `data.local.yaml`. Repeatable; the first is the control |
| `--new-condition LABEL` | Creates `conditions/<label>/` with a template `config.yaml` and `structures/`. Repeatable |
| `--equilibration TEXT` | The study's equilibration window; `0ns` with a note when left out |
| `--holder NAME` | Copyright holder in both LICENSE files. Default: `git config user.name` |
| `--no-git` | Do not make a git repository |

It refuses a folder that already holds a `study.yaml`, and exits 2 on an
unreadable config.

### polyzymd study add-condition

```bash
polyzymd study add-condition LABEL (--config CONFIG | --from OTHER_LABEL | --new) [--study PATH]
```

Adds a condition to an existing study, in `conditions/<label>/`:

- `--new` writes a template `config.yaml` to fill in, and `structures/` with
  placeholder files.
- `--from OTHER_LABEL` copies the config of that condition of the study, with
  the input files it names.
- `--config CONFIG` copies `CONFIG`, with the input files it names, as
  `study init --condition` does. When the config's `scratch_directory` holds
  runs of the config, that folder is written to `data.local.yaml` and printed
  as `data <label>: <path> (from the config's scratch_directory)`.

The new config writes its runs into `runs/<study>/<label>/` of the project,
or into `runs/<label>/` of a study in no project. Git ignores `runs/`, and
freeze never publishes it. The command prints this as a warning, because
trajectories can use a lot of disk space: on a cluster, set
`scratch_directory` in the config to scratch storage. The condition is added
as one line under `conditions:` in `study.yaml`, keeping the rest of the
file. Exits 2 when not exactly one of `--config`, `--from` and `--new` is
given, `OTHER_LABEL` is not a condition of the study, or the label or its
folder is taken.

### polyzymd study locate

```bash
polyzymd study locate DIRECTORY [--study PATH] [--verify]
```

Searches `DIRECTORY` (six levels deep) for each condition's run directories,
named by its config's `naming_template`, and writes the folder holding the
most of them to `data.local.yaml` beside `study.yaml`, keeping entries of
conditions it does not find. With a `manifest.json` from `study freeze`, it
prefers a folder whose files have the recorded sizes (and with `--verify`,
SHA-256), and prints `<label>: <n> files match manifest.json`. When two
conditions name their runs alike, a folder named for the condition
(`no_polymer/`) is chosen before file sizes, and one folder is never written
for two such conditions. Conditions whose runs are named apart can share one
folder. Prints `<label>: replicates [<numbers>] under <folder>` per
condition found. Writes nothing when no condition is found. A condition
with a missing or different file is not written, and its entry in
`data.local.yaml` is removed (with a printed line) when it names that folder. Exits 2 when a
condition is not found, two conditions with alike run names are found in one
folder, or a file is
missing or different. See
{doc}`../how_to/study_folder` and {doc}`../how_to/study_freeze`.

### polyzymd study freeze

```bash
polyzymd study freeze [PATH] [--tag NAME]
```

Prepares the study for publication: checks the metadata, the git state and
whether each analysis's stored results match the study; hashes the
trajectories and writes their engine inputs and final frames to `deposit/`;
writes `manifest.json`, `CITATION.cff`, `.zenodo.json`, `md_checklist.yaml`
and `system_summary.csv`; commits those and `results/`, tags the commit
(`study-v1`, `study-v2`, ... unless `--tag`), and lays out `deposit/`: the
files to add to a Zenodo upload in `deposit/upload/`, the trajectory files in
record-sized batches in `deposit/trajectories.csv`, and the steps in
`deposit/UPLOAD.md`. It uploads and publishes nothing. Every gap is printed as
`warning:` and none stops it. Exits 2 when the study file or its metadata cannot be read, or
the tag exists. See {doc}`../how_to/study_freeze`.

---

## polyzymd hash-trajectories

```bash
polyzymd hash-trajectories (-c CONFIG ... | --study PATH) [--replicates SPEC] [--verify] [--dry-run] [--rehash-changed]
```

Records the SHA-256 and size of every finished trajectory file that has
none yet in `trajectory_hashes.json`, in the run's engine working directory:
for runs whose runner did not record them, such as runs that finished before
PolyzyMD recorded hashes, downsampled copies, and GROMACS runs. Stored
analysis results and frozen studies then identify those trajectories
without reading them again.

It never writes `progress.json`, which only the runner writes, so it cannot
change which segments analyses read. It works for every simulation engine;
the engine of each config says which files are finished trajectories:

| Engine | Trajectory files hashed |
|---|---|
| OpenMM | `production_<n>/production_<n>_trajectory.dcd` of each segment analyses read: with `progress.json`, the completed and interrupted ones; without it, every one |
| GROMACS | `gromacs/prod.xtc`, and `prod_nojump.xtc` and `prod_centered.xtc` when present, once the run has completed or when it has no `progress.json` |

`trajectory_hashes.json` maps each file's path, relative to the engine
working directory, to its SHA-256 and size. A hash the OpenMM runner
recorded for a segment in `progress.json` takes precedence over it.

It is idempotent. A file whose recorded hash has the file's size is left as
it is without reading it, so running the command again changes nothing, and
`trajectory_hashes.json` is written, atomically, only when a hash was added.
A recorded hash is never overwritten: a recorded size that differs from the
file's, or with `--verify` a recomputed hash that differs, is printed as
`conflict:` and exits 2. When the file was legitimately extended, such as a
GROMACS run continued after it was hashed, `--rehash-changed` hashes it
again and replaces the entry in `trajectory_hashes.json`; a hash the runner
recorded is never replaced.

| Option | Description |
|---|---|
| `-c, --config PATH` | A simulation `config.yaml`; every run it finds is hashed. Repeatable |
| `--study PATH` | Every condition of a study, using its `data.local.yaml` |
| `--replicates SPEC` | Only these replicates, for example `1-5` |
| `--verify` | Also rehash files with a recorded hash and compare |
| `--dry-run` | Report what would be hashed, and write nothing |
| `--rehash-changed` | Hash again the entries of `trajectory_hashes.json` whose file size changed |

It prints one line per replicate, such as `SBMA 50% replicate 1 (openmm):
hashed 11, already recorded 1`, and `<label>: no replicates found under <folder>` for
a config whose runs are not where it says. Reading takes about a second per
gigabyte, so on a cluster run it in a batch job.

---

## Retired commands

`polyzymd init ...` is hidden. It accepts any arguments, says to make a project
with `polyzymd project init PATH --study LABEL`, then a condition with
`polyzymd study add-condition LABEL --new`, which writes the template config,
and exits 2.

---

## Environment Variables

PolyzyMD expands environment variables in configuration paths:

| Variable | Example | Description |
|----------|---------|-------------|
| `$USER` | `me` | Current username |
| `$HOME` | `/home/me` | Home directory |
| `~` | `/home/me` | Home directory shortcut |
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
| 2 | Typed analysis error from {ref}`polyzymd analyze <cli-analyze>`, whose message and fix are printed on stderr, one line each |
| 3 | `polyzymd check-progress`: error; the SLURM job does not resubmit |
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
