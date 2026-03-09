# HPC and SLURM Guide

This guide covers running PolyzyMD simulations on HPC clusters using SLURM.

## Overview

Long MD simulations often exceed HPC time limits (typically 24-48 hours). PolyzyMD solves this with **self-resubmitting jobs**: each SLURM job runs a single simulation segment, checks whether work remains, and resubmits itself if not finished. This approach is simpler and more robust than dependency chains — every job is identical, and the simulation resumes correctly after wall-time limits, preemptions, or node failures.

## User Workflow

The typical workflow for running a PolyzyMD simulation on an HPC cluster is:

1. **Create a simulation directory** with your configuration and input files
2. **Write a `config.yaml`** file with your simulation parameters
3. **Generate job scripts** using `polyzymd submit --dry-run`
4. **Review** the generated scripts
5. **Submit the jobs** for real

### Step-by-Step Example

```bash
# 1. Create your simulation directory
mkdir -p my_simulation/structures
cd my_simulation

# 2. Copy your input files
cp /path/to/enzyme.pdb structures/
cp /path/to/substrate.sdf structures/
# If using pre-built polymers:
cp -r /path/to/polymer_sdfs ./ATRP_EGPMA_SBMA_5-mer/

# 3. Create your config.yaml (see Configuration Guide)
# The config should reference paths relative to this directory:
#   enzyme.pdb_path: "structures/enzyme.pdb"
#   substrate.sdf_path: "structures/substrate.sdf"
#   polymers.sdf_directory: "ATRP_EGPMA_SBMA_5-mer"

# 4. Test with a dry run first
polyzymd submit -c config.yaml --preset testing --dry-run

# 5. Review the generated scripts
cat job_scripts/r1_300K_LipA.sh

# 6. Submit for real (quick test first)
polyzymd submit -c config.yaml --preset testing --time-limit 0:05:00 --replicates 1

# 7. Once testing passes, submit production jobs
polyzymd submit -c config.yaml --preset aa100 --replicates 1-5 --email your@email.edu
```

---

## Directory Structure

PolyzyMD supports separating:

- **Projects directory**: Long-term storage for scripts, logs, configs
- **Scratch directory**: High-performance storage for trajectories

```
/projects/$USER/polyzymd/           # Long-term storage
├── my_simulation/                  # Your simulation directory
│   ├── config.yaml                 # Main configuration
│   ├── structures/                 # Input structure files
│   │   ├── enzyme.pdb
│   │   └── substrate.sdf
│   ├── ATRP_EGPMA_SBMA_5-mer/     # Pre-built polymer SDFs (optional)
│   │   ├── EGPMA-SBMA_AAAAA_5-mer_charged.sdf
│   │   └── ...
│   ├── job_scripts/                # Generated SLURM scripts
│   │   ├── r1_300K_LipA.sh        # One script per replicate
│   │   ├── r2_300K_LipA.sh
│   │   └── ...
│   └── slurm_logs/                 # Job output files
│       └── r1_300K_LipA.12345.out

/scratch/alpine/$USER/polyzymd_sims/  # High-performance storage
├── LipA_Substrate_EGPMA-SBMA_10pct_300K_run1/
│   ├── system.pdb
│   ├── progress.json               # Progress tracking file
│   ├── equilibration/
│   │   └── trajectory.dcd
│   ├── production_0/
│   │   ├── trajectory.dcd
│   │   ├── checkpoint.chk
│   │   └── state_data.csv
│   └── production_1/
│       └── ...
└── LipA_Substrate_EGPMA-SBMA_10pct_300K_run2/
```

## Configuring Directories

### In YAML Configuration

Environment variables (`$USER`, `$HOME`, etc.) and `~` are automatically expanded:

```yaml
output:
  projects_directory: "/projects/$USER/polyzymd/my_simulation"
  scratch_directory: "/scratch/alpine/$USER/polyzymd_sims"
  job_scripts_subdir: "job_scripts"
  slurm_logs_subdir: "slurm_logs"
```

You can also use `~` for home directory:

```yaml
output:
  projects_directory: "~/polyzymd/my_simulation"
```

### Via CLI Override

```bash
polyzymd submit -c config.yaml \
    --projects-dir /projects/$USER/polyzymd \
    --scratch-dir /scratch/alpine/$USER/simulations \
    --replicates 1-5
```

---

## SLURM Presets

PolyzyMD includes presets for common HPC configurations:

| Preset | Partition | GPUs | Time Limit | Memory | Description |
|--------|-----------|------|------------|--------|-------------|
| `aa100` | aa100 | 1x A100 | 24h | 3GB | CU Boulder Alpine — NVIDIA A100 |
| `al40` | al40 | 1x L40 | 24h | 3GB | CU Boulder Alpine — NVIDIA L40 |
| `blanca-shirts` | blanca-shirts | 1x | 24h | 3GB | CU Boulder Blanca — Shirts lab partition |
| `testing` | atesting_a100 | 1x | 6min | 3GB | CU Boulder Alpine — quick tests |
| `bridges2` | GPU-shared | 1x V100-32 | 24h | (per-GPU) | PSC Bridges2 — NVIDIA V100 32GB |

### Using Presets

```bash
# Use A100 GPUs
polyzymd submit -c config.yaml --preset aa100

# Use testing partition for quick tests
polyzymd submit -c config.yaml --preset testing
```

### Overriding Time Limit

You can override the preset's time limit using `--time-limit`:

```bash
# Use testing preset with a 2-minute time limit
polyzymd submit -c config.yaml --preset testing --time-limit 0:02:00

# Use A100 with a 12-hour limit instead of 24h
polyzymd submit -c config.yaml --preset aa100 --time-limit 12:00:00
```

**Time format options:**
- `MM:SS` - minutes and seconds (e.g., `2:00` for 2 minutes)
- `HH:MM:SS` - hours, minutes, seconds (e.g., `0:02:00`)
- `D-HH:MM:SS` - days, hours, minutes, seconds (e.g., `1-00:00:00` for 1 day)

This is especially useful for:
- Quick testing with short time limits
- Adjusting for segment duration requirements
- Working within specific QOS constraints

### Custom SLURM Settings

For custom configurations, edit the generated scripts in `job_scripts/` before submitting.

---

## Bridges2 (PSC)

[Bridges2](https://www.psc.edu/resources/bridges-2/) is the Pittsburgh Supercomputing Center (PSC) GPU cluster. It uses slightly different SLURM conventions than CU Boulder Alpine, and polyzymd handles these differences automatically via the `bridges2` preset.

### Key Differences from Alpine

| Feature | Alpine (`aa100`) | Bridges2 (`bridges2`) |
|---------|-----------------|----------------------|
| GPU directive | `--gres=gpu:N` | `--gpus=<type>:N` |
| Nodes/tasks | `--nodes=1` + `--ntasks=1` | `-N 1` (single line) |
| QoS | `--qos=normal` | *(omitted — not used)* |
| Memory | `--mem=3G` | *(omitted — per-GPU allocation)* |
| Account | ucb-group (in preset) | *(omitted — inferred from login)* |
| Module load | `ml miniforge` | `ml anaconda3/2024.10-1` |
| Conda frontend | `mamba activate` | `conda activate` |
| Default time limit | 24h | 24h |

### Account

Bridges2 infers the billing allocation from your login session, so **no `--account` directive is emitted by default**. If you have multiple allocations and need to charge a specific one, pass `--account`:

```bash
polyzymd submit -c config.yaml \
    --preset bridges2 \
    --account chm250017p \
    --replicates 1-3 \
    --email collaborator@pitt.edu
```

```{note}
Unlike Alpine presets, Bridges2 scripts omit the `#SBATCH --account=` line
entirely when no account is specified. The `--account` CLI flag is optional
for Bridges2 (it is required on Alpine where the preset always sets a
group account).
```

### GPU Type Selection

Bridges2 has multiple GPU types available. Use `--gpu-type` to select:

| Flag value | GPU | VRAM |
|------------|-----|------|
| `v100-32` *(default)* | NVIDIA V100 | 32 GB |
| `v100-16` | NVIDIA V100 | 16 GB |
| `l40s-48` | NVIDIA L40S | 48 GB |
| `h100-80` | NVIDIA H100 | 80 GB |

```bash
# Default (V100 32GB) — good balance of availability and memory
polyzymd submit -c config.yaml \
    --preset bridges2 \
    --account abc123_gpu

# High-memory GPU for large systems
polyzymd submit -c config.yaml \
    --preset bridges2 \
    --account abc123_gpu \
    --gpu-type h100-80
```

### Full Bridges2 Workflow

```bash
# 1. Dry run — inspect scripts before submitting
polyzymd submit -c config.yaml \
    --preset bridges2 \
    --account abc123_gpu \
    --replicates 1-3 \
    --dry-run

# 2. Inspect the generated SBATCH directives
head -20 job_scripts/r1_300K_LipA.sh
# You should see:
#   #SBATCH --partition=GPU-shared
#   #SBATCH -N 1                    ← single-line nodes directive
#   #SBATCH --gpus=v100-32:1        ← type-specific GPU directive
#   (no --qos line)
#   (no --mem line)
#   (no --account line — inferred from login)

# 3. Submit for real
polyzymd submit -c config.yaml \
    --preset bridges2 \
    --account abc123_gpu \
    --replicates 1-3 \
    --email collaborator@pitt.edu
```

### Bridges2 Directory Structure

On Bridges2, use Ocean storage for long-term data and local scratch for active simulations:

```
/ocean/projects/abc123_gpu/$USER/polyzymd/   # Long-term storage
├── my_simulation/
│   ├── config.yaml
│   ├── structures/
│   │   ├── enzyme.pdb
│   │   └── substrate.sdf
│   ├── job_scripts/
│   │   ├── r1_300K_LipA.sh
│   │   └── ...
│   └── slurm_logs/

/local/scratch/$USER/polyzymd_sims/          # High-performance local scratch
├── LipA_Substrate_300K_run1/
│   ├── system.pdb
│   ├── progress.json
│   ├── equilibration/
│   └── production_0/
```

Set these paths in your `config.yaml`:

```yaml
output:
  projects_directory: "/ocean/projects/abc123_gpu/$USER/polyzymd/my_simulation"
  scratch_directory: "/local/scratch/$USER/polyzymd_sims"
```

Or override on the CLI:

```bash
polyzymd submit -c config.yaml \
    --preset bridges2 \
    --account abc123_gpu \
    --projects-dir "/ocean/projects/abc123_gpu/$USER/polyzymd/my_simulation" \
    --scratch-dir "/local/scratch/$USER/polyzymd_sims"
```

---

## Self-Resubmitting Job Model

### How It Works

PolyzyMD generates one identical SLURM script per replicate. Each job:

1. Calls `polyzymd run-segment` to run the next segment of work
2. After the segment finishes (or is interrupted), calls `polyzymd check-progress` to see if the simulation is complete
3. If work remains, resubmits itself via `sbatch "$THIS_SCRIPT"`
4. If the simulation is complete, exits cleanly

```
  ┌───────────────────────────┐
  │  Job submits itself       │
  │                           │◄──────────────────┐
  │  1. run-segment           │                   │
  │     (build/eq/prod OR     │                   │
  │      continue from last)  │                   │
  │                           │                   │
  │  2. check-progress        │     resubmit      │
  │     └─ complete? exit 0   │                   │
  │     └─ work remains? ─────┼───────────────────┘
  │     └─ error? exit $RC    │
  └───────────────────────────┘
```

This model is simpler and more fault-tolerant than dependency chains:
- **No `afterany` dependencies** — each job is independent
- **Automatic recovery** — if a job is interrupted by wall-time or preemption, it saves a checkpoint, resubmits, and the next invocation picks up where it left off
- **Idempotent** — each job scans the filesystem to determine what work remains, so the same script can be resubmitted manually at any time

### Progress Tracking

PolyzyMD tracks simulation progress in a `progress.json` file in the working directory. This file records:
- Which segments have completed
- How many steps each segment ran
- Whether any segments were interrupted
- The total steps requested vs. completed

On startup, `run-segment` validates the progress file against the filesystem (checking for `production_N/` directories and their contents) to ensure consistency.

### Configuring Production Duration

Specify the total simulation time in your config. PolyzyMD determines segment boundaries automatically based on wall-time:

```yaml
simulation_phases:
  production:
    duration: 100.0    # 100 ns total production time
```

```{tip}
You don't need to configure segments manually. Each job runs as much
production time as it can before the wall-time limit, checkpoints, and
resubmits. The segment duration is determined at runtime by the
``--segment-time`` and ``--segment-frames`` options passed to ``run``
(which ``submit`` computes automatically from your config).
```

---

(smart-restart)=
## Smart Restart & Fault Tolerance

When you run 60 replicates, each resubmitting multiple times, robustness is
critical. PolyzyMD's **smart restart** system handles interruptions
automatically. The generated scripts already include everything described
below — you do not need to configure anything. This section explains what
happens under the hood so you can debug issues or adapt the approach to
other workflows.

### What Happens Automatically

Every generated SLURM script includes three pieces of fault-tolerance
infrastructure:

1. **Signal handling** — Python-side handlers catch SIGUSR1 (wall-time
   warning) and SIGTERM (preemption), save an interrupted checkpoint, and
   exit with code 99.
2. **Signal forwarding** — Bash trap + background + wait pattern forwards
   signals from the SLURM batch shell to the Python child process.
3. **Progress tracking** — After each segment (whether completed or
   interrupted), the progress file is updated so the next invocation
   knows exactly where to resume.

### The Three Scenarios

| Scenario | Signal | What Happens | Outcome |
|----------|--------|--------------|---------|
| **Wall-time warning** | `SIGUSR1` (5 min before limit) | Interrupted state saved, progress updated, exit 99 | Job resubmits and resumes |
| **Preemption** | `SIGTERM` (120 s grace on Blanca) | Interrupted state saved, progress updated, exit 99 | Job resubmits and resumes |
| **Hard crash** | None (OOM, segfault, node failure) | No state saved | Job resubmits; `run-segment` detects incomplete segment and handles it |

```{note}
The wall-time signal is configured via `#SBATCH --signal=B:USR1@300`, which
tells SLURM to send `SIGUSR1` to the batch shell 300 seconds (5 minutes)
before the time limit expires. This gives the simulation enough time to
save a full OpenMM state (~10-30 seconds on GPU).
```

### Interrupted State Files

When an interrupt is detected, the signal handler writes three files into
the current segment's directory (e.g. `production_3/`):

| File | Purpose |
|------|---------|
| `interrupted_state.xml` | Portable OpenMM state (positions, velocities, forces) |
| `interrupted_system.xml` | Serialized OpenMM System (force field parameters) |
| `INTERRUPTED` | Marker file with step-count metadata for recovery |

The `INTERRUPTED` marker contains the information needed for recovery:

```
segment_index=3
steps_completed=1250000
total_steps=2500000
remaining_steps=1250000
```

### Signal Forwarding: Why trap + background + wait?

SLURM sends signals to the **batch shell process**, not to child processes.
Bash ignores `SIGUSR1` by default, so without explicit forwarding, the
Python simulation never sees the signal. The generated scripts use a
standard pattern to solve this:

```bash
# Background the Python process
polyzymd run-segment -c "$CONFIG_PATH" -r "$REPLICATE" --scratch-dir "$WORKING_DIR" &
CHILD_PID=$!

# Trap signals and forward them to the child
trap 'forward_signal USR1' USR1
trap 'forward_signal TERM' TERM

# Wait in a loop (wait is interrupted by trapped signals)
wait "$CHILD_PID"
RC=$?
while kill -0 "$CHILD_PID" 2>/dev/null; do
    wait "$CHILD_PID"
    RC=$?
done
```

```{warning}
Do not remove the `trap`, backgrounding (`&`), or `wait` loop from the
generated scripts. Without them, signals will not reach the Python process
and graceful shutdown will not work.
```

### Manually Triggering an Interrupt

You can test graceful shutdown or manually stop a running simulation by
sending `SIGUSR1` via `scancel`:

```bash
# Send USR1 to a specific job
scancel --signal=USR1 <job_id>

# The job will save interrupted state, update progress, and exit with code 99
# The resubmission logic then resubmits the job to continue
```

This is useful when you realize a simulation has a problem and want to
stop it cleanly without losing progress.

### Cancelling a Job Permanently

Because `scancel` sends `SIGTERM`, and our scripts treat SIGTERM the same as
SIGUSR1 (save state, exit 99, resubmit), a plain `scancel <job_id>` will
**not** permanently stop a simulation — the job will save its state and
resubmit itself.

To truly cancel a simulation so it does not restart, use `--signal=KILL`:

```bash
# Permanently stop a job (no state saved, no resubmission)
scancel --signal=KILL <job_id>

# Equivalent shorthand
scancel -s KILL <job_id>
```

`SIGKILL` cannot be caught or trapped by bash or Python, so the process
dies immediately and the resubmission logic never runs. The last
*completed* segment's checkpoint is still intact — only in-progress work
since the last checkpoint is lost.

| Command | Saves state? | Resubmits? | Use case |
|---------|-------------|-----------|----------|
| `scancel <job_id>` | Yes | **Yes** | Don't use this to permanently stop a simulation |
| `scancel --signal=USR1 <job_id>` | Yes | **Yes** | Graceful stop (saves progress, continues later) |
| `scancel --signal=KILL <job_id>` | No | **No** | Permanently cancel a simulation |

```{tip}
If you need to cancel **all replicates** of a simulation, cancel them all
at once so that none resubmit before you can cancel the others:

    scancel --signal=KILL <job_id_1> <job_id_2> <job_id_3>

Or cancel all your jobs:

    scancel --signal=KILL -u $USER
```

### Manual Recovery

If a simulation is stalled (e.g., the SLURM job exited without resubmitting),
use the `recover` command to inspect status and optionally resume:

```bash
# Check status only
polyzymd recover -c config.yaml -r 1

# Submit a recovery job
polyzymd recover -c config.yaml -r 1 --submit --preset blanca-shirts

# Dry-run (show what would be submitted)
polyzymd recover -c config.yaml -r 1 --submit --dry-run
```

See the {ref}`CLI Reference <cli-recover>` section below for full option
details.

---

## Submitting Jobs

### Dry Run (Recommended First)

Generate scripts without submitting:

```bash
polyzymd submit -c config.yaml \
    --replicates 1-3 \
    --preset aa100 \
    --dry-run
```

Inspect generated scripts:

```bash
cat job_scripts/r1_300K_LipA.sh
```

### Submit for Real

```bash
polyzymd submit -c config.yaml \
    --replicates 1-3 \
    --preset aa100 \
    --email your.email@university.edu
```

### Replicate Specification

```bash
# Single replicate
--replicates 1

# Range
--replicates 1-5

# Specific replicates
--replicates 1,3,5,7

# Combined
--replicates 1-3,5,7-10
```

---

## Monitoring Jobs

### Check Job Status

```bash
# All your jobs
squeue -u $USER

# Specific job
scontrol show job <job_id>

# Watch jobs update
watch -n 30 'squeue -u $USER'
```

### View Job Output

```bash
# Real-time output
tail -f slurm_logs/r1_300K_LipA*.out

# Check for errors
grep -i error slurm_logs/*.out
grep -i fail slurm_logs/*.out
```

### Check Simulation Progress

Use the `check-progress` command to query the progress file:

```bash
# Check a specific replicate
polyzymd check-progress -c config.yaml -r 1

# Example output:
# Progress: 12500000/50000000 steps (25.0%), 5 segment(s)
# Status: in_progress — 75.000 ns remaining
```

Or inspect files directly:

```bash
# List trajectory files
ls -la /scratch/$USER/polyzymd_sims/*/production_*/trajectory.dcd

# Check trajectory sizes
du -h /scratch/$USER/polyzymd_sims/*/production_*/*.dcd
```

---

## Handling Failures

Most failures are handled automatically by the {ref}`smart restart <smart-restart>`
system. This section covers the cases where manual intervention is needed.

### Wall-Time or Preemption Interrupts

**No action needed.** The smart restart system saves an interrupted checkpoint,
updates the progress file, and the job resubmits itself. Check the SLURM
log to confirm:

```bash
# Look for the graceful shutdown message
grep -i "interrupted\|graceful\|resubmit" slurm_logs/r1_300K_*.out
```

### Hard Crash (OOM, Segfault, Node Failure)

If a job crashes without saving state, the resubmission logic still runs
(since the crash only kills the child Python process, not the bash wrapper).
The resubmitted job's `run-segment` will detect the incomplete segment and
handle it appropriately.

To diagnose issues:

1. **Check the error**:
   ```bash
   cat slurm_logs/r1_300K_*.out
   ```

2. **Fix the issue** (e.g. increase memory with `--memory 8G`)

3. **Resume**:
   ```bash
   # Option 1: Resubmit the existing script
   sbatch job_scripts/r1_300K_LipA.sh

   # Option 2: Use recover to inspect and resume
   polyzymd recover -c config.yaml -r 1 --submit --preset aa100
   ```

### Checking Progress

Use `check-progress` to get a quick status report:

```bash
polyzymd check-progress -c config.yaml -r 1
```

Or for a more detailed view with per-segment status:

```bash
polyzymd recover -c config.yaml -r 1
```

### Start Fresh

To restart a simulation from scratch:

```bash
# Remove old output
rm -rf /scratch/$USER/polyzymd_sims/LipA_*_run1/

# Resubmit
polyzymd submit -c config.yaml --replicates 1 --preset aa100
```

---

## Generated Script Structure

The submit command generates **one script per replicate**. Each script is a
self-resubmitting job that handles the entire simulation lifecycle: building,
equilibration, production segments, interruptions, and resubmission.

```bash
#!/bin/bash
#SBATCH --partition=aa100
#SBATCH --job-name=r1_300K_LipA
#SBATCH --output=slurm_logs/r1_300K_LipA.%A_%a.out
#SBATCH --qos=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=3G
#SBATCH --time=23:59:59
#SBATCH --gres=gpu:1
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=your@email.edu
#SBATCH --account=ucb625_asc1
#SBATCH --signal=B:USR1@300
#SBATCH --no-requeue

# =============================================================================
# PolyzyMD Self-Resubmitting Simulation Job
# Generated by polyzymd — do not edit manually
# =============================================================================

module purge 2>/dev/null || true
ml miniforge 2>/dev/null || true

eval "$(conda shell.bash hook)"
conda activate polyzymd-env

set -e

export INTERCHANGE_EXPERIMENTAL=1

# Resolve this script's path for self-resubmission
# ($SLURM_JOB_SCRIPT is only available in SLURM >= 22.05)
THIS_SCRIPT="${SLURM_JOB_SCRIPT:-$(realpath "$0")}"

CONFIG_PATH="/projects/$USER/polyzymd/my_simulation/config.yaml"
REPLICATE=1
WORKING_DIR="/scratch/alpine/$USER/polyzymd_sims/LipA_300K_run1"

mkdir -p "$WORKING_DIR"

echo "=================================================="
echo "PolyzyMD self-resubmitting job"
echo "Config:    $CONFIG_PATH"
echo "Replicate: $REPLICATE"
echo "Work dir:  $WORKING_DIR"
echo "Job ID:    ${SLURM_JOB_ID:-local}"
echo "Timestamp: $(date)"
echo "=================================================="

# Signal forwarding (see Smart Restart docs)
CHILD_PID=""
forward_signal() {
    if [ -n "$CHILD_PID" ] && kill -0 "$CHILD_PID" 2>/dev/null; then
        echo "Forwarding $1 to Python process (PID $CHILD_PID)"
        kill -"$1" "$CHILD_PID"
    fi
}
trap 'forward_signal USR1' USR1
trap 'forward_signal TERM' TERM

# Run the next segment (backgrounded for signal forwarding)
polyzymd run-segment \
    -c "$CONFIG_PATH" \
    -r "$REPLICATE" \
    --scratch-dir "$WORKING_DIR" &
CHILD_PID=$!

# Wait for the child process
set +e
wait "$CHILD_PID" 2>/dev/null
RC=$?
while kill -0 "$CHILD_PID" 2>/dev/null; do
    wait "$CHILD_PID" 2>/dev/null
    RC=$?
done
set -e

echo "run-segment exited with code $RC at $(date)"

# --- Resubmission logic ---

if [ $RC -ne 0 ] && [ $RC -ne 99 ]; then
    echo "FATAL: run-segment failed (exit code $RC) — NOT resubmitting"
    exit $RC
fi

# Check whether more work remains
set +e
polyzymd check-progress -c "$CONFIG_PATH" -r "$REPLICATE" --scratch-dir "$WORKING_DIR"
PROGRESS_RC=$?
set -e

if [ $PROGRESS_RC -eq 0 ]; then
    echo "Simulation complete — no resubmission needed."
    exit 0
fi

# Work remains — resubmit this same script
echo "Work remains — resubmitting job..."
sbatch "$THIS_SCRIPT"
SUBMIT_RC=$?

if [ $SUBMIT_RC -eq 0 ]; then
    echo "Resubmitted successfully."
else
    echo "WARNING: sbatch resubmission failed (exit code $SUBMIT_RC)"
    echo "You can manually resume with:"
    echo "  sbatch $THIS_SCRIPT"
    exit 1
fi

exit 0
```

### Key Features of the Generated Script

| Feature | How It Works |
|---------|-------------|
| **Signal forwarding** | `trap` + background `&` + `wait` loop ensures SIGUSR1/SIGTERM reach the Python process |
| **Unified entry point** | `polyzymd run-segment` handles both initial (build + eq + seg 0) and continuation segments |
| **Progress checking** | `polyzymd check-progress` returns exit code 0 (complete) or 1 (work remains) |
| **Self-resubmission** | `sbatch "$THIS_SCRIPT"` resubmits the exact same script (path resolved at script start via `$SLURM_JOB_SCRIPT` or `realpath "$0"`) |
| **Error handling** | Non-zero, non-99 exit codes abort without resubmitting |

---

## Best Practices

### 1. Always Test First

```bash
# Generate scripts without submitting (dry run)
polyzymd submit -c config.yaml --preset testing --dry-run

# Quick test with 2-minute time limit
polyzymd submit -c config.yaml \
    --preset testing \
    --time-limit 0:02:00 \
    --replicates 1

# Or a slightly longer test
polyzymd submit -c config.yaml \
    --preset testing \
    --time-limit 0:05:00 \
    --replicates 1
```

### 2. Monitor Early Segments

Watch the first segment complete to catch issues early:

```bash
tail -f slurm_logs/r1_300K_*.out
```

### 3. Back Up Important Data

Scratch is often purged. Copy completed simulations to projects:

```bash
# After simulation completes
cp -r /scratch/$USER/polyzymd_sims/LipA_300K_run1 \
      /projects/$USER/completed_simulations/
```

### 4. Use Email Notifications

```bash
polyzymd submit -c config.yaml --email you@university.edu
```

You'll receive emails when jobs start, end, or fail.

### 5. Segment Duration Guidelines

| Cluster Time Limit | Approximate Production per Segment |
|--------------------|-------------------------------------|
| 1 hour (testing) | 0.5 - 1 ns |
| 24 hours | 8 - 12 ns |
| 48 hours | 20 - 30 ns |
| 7 days | 50 - 100 ns |

---

## CLI Reference

### `polyzymd submit`

Submit self-resubmitting simulation jobs to SLURM.

```bash
polyzymd submit -c CONFIG [OPTIONS]
```

**Required:**
- `-c, --config PATH` - Path to YAML configuration file

**Options:**
- `-r, --replicates RANGE` - Replicate range (e.g., "1-5", "1,3,5"). Default: "1"
- `--preset PRESET` - SLURM preset: aa100, al40, blanca-shirts, testing, bridges2. Default: aa100
- `--account ACCOUNT` - HPC allocation account ID (required for Bridges2)
- `--gpu-type TYPE` - GPU type for Bridges2: v100-16, v100-32, l40s-48, h100-80. Default: v100-32
- `--scratch-dir PATH` - Override scratch directory for simulation output
- `--projects-dir PATH` - Override projects directory for scripts/logs
- `--output-dir PATH` - Directory for job scripts. Default: {projects_dir}/job_scripts
- `--email EMAIL` - Email for job notifications
- `--time-limit TIME` - Override SLURM time limit (HH:MM:SS)
- `--memory SIZE` - Override SLURM memory allocation (e.g., "4G"). Bridges2 omits --mem by default (per-GPU allocation)
- `--openff-logs` - Enable verbose OpenFF logs in generated job scripts (for debugging)
- `--dry-run` - Generate scripts without submitting

### `polyzymd run-segment`

Unified entry point for SLURM jobs. Determines what work remains and runs the next segment.

```bash
polyzymd run-segment -c CONFIG [OPTIONS]
```

**Required:**
- `-c, --config PATH` - Path to YAML configuration file

**Options:**
- `-r, --replicate INT` - Replicate number. Default: 1
- `--scratch-dir PATH` - Override scratch directory for simulation output
- `--skip-build` - Skip system building (use existing) for initial segment

**Behavior:**
- If no segments exist: builds system, equilibrates, runs segment 0
- If segments exist but simulation is incomplete: continues from last completed segment
- If simulation is complete: exits 0 immediately

**Exit codes:**
- `0` - Segment completed successfully
- `1` - Error
- `99` - Graceful interruption (wall-time signal)

### `polyzymd check-progress`

Check whether a simulation is complete. Used by SLURM resubmission logic.

```bash
polyzymd check-progress -c CONFIG [OPTIONS]
```

**Required:**
- `-c, --config PATH` - Path to YAML configuration file

**Options:**
- `-r, --replicate INT` - Replicate number. Default: 1
- `--scratch-dir PATH` - Override scratch directory

**Exit codes:**
- `0` - Simulation complete (do NOT resubmit)
- `1` - Work remains (resubmit)

### `polyzymd run`

Run a complete simulation (build + equilibration + first production segment).

```bash
polyzymd run -c CONFIG [OPTIONS]
```

**Required:**
- `-c, --config PATH` - Path to YAML configuration file

**Options:**
- `-r, --replicate INT` - Replicate number. Default: 1
- `--scratch-dir PATH` - Scratch directory for simulation output
- `--projects-dir PATH` - Projects directory for scripts/logs
- `--segment-time FLOAT` - Override production time per segment (ns)
- `--segment-frames INT` - Override frames per segment
- `--skip-build` - Skip system building (use existing)

### `polyzymd continue`

Continue a simulation from a previous segment checkpoint. Prefer `run-segment` for SLURM workflows.

```bash
polyzymd continue -w WORKING_DIR -s SEGMENT -t TIME [OPTIONS]
```

**Required:**
- `-w, --working-dir PATH` - Working directory with previous segment
- `-s, --segment INT` - Segment index to run (1-based)
- `-t, --segment-time FLOAT` - Duration of this segment (ns)

**Options:**
- `-n, --num-samples INT` - Number of frames to save. Default: 250

(hpc-recover)=
### `polyzymd recover`

Resume a stalled or interrupted simulation. Scans the working directory,
loads progress state, and reports how much work remains. With `--submit`,
generates and submits a self-resubmitting SLURM job that will automatically
continue from the last completed segment.

```bash
polyzymd recover -c CONFIG [OPTIONS]
```

**Required:**
- `-c, --config PATH` - Path to YAML configuration file

**Options:**
- `-r, --replicate INT` - Replicate number. Default: 1
- `--scratch-dir PATH` - Override scratch directory
- `--preset PRESET` - SLURM preset for recovery job. Default: aa100
- `--submit / --no-submit` - Submit a recovery job (default: status only)
- `--dry-run` - Show what would be submitted without submitting

**Examples:**
```bash
# Check status only
polyzymd recover -c config.yaml -r 1

# Submit a recovery job
polyzymd recover -c config.yaml -r 1 --submit --preset blanca-shirts

# Dry-run (show what would be submitted)
polyzymd recover -c config.yaml -r 1 --submit --dry-run
```

**Example output (status only):**
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

---

## Troubleshooting

### "Job pending forever"

```bash
squeue -u $USER
# Check REASON column
```

Common reasons:
- `Resources` - Waiting for GPUs
- `Priority` - Queue is busy

### "Module not found in job"

If you see errors like "Lmod has detected the following error: The following module(s) are unknown", the cluster may have different module configurations on different partitions.

Check available modules:
```bash
module spider miniforge
module spider anaconda
```

You may need to edit the generated scripts to use a different module:
```bash
# Instead of:
module load miniforge
mamba activate polymerist-env

# Try:
module load anaconda
conda activate polymerist-env

# Or source conda directly:
source /curc/sw/anaconda3/latest/etc/profile.d/conda.sh
conda activate polymerist-env
```

### "Out of memory"

There are two types of out-of-memory errors:

**GPU Memory (CUDA OOM):**
```
CUDA out of memory
```

Reduce system size:
- Decrease `box.padding`
- Use fewer polymers
- Use smaller production `samples` (fewer frames saved)

**System Memory (SLURM OOM):**
```
slurmstepd: error: Detected 1 oom_kill event in StepId=...
```

The job exceeded its RAM allocation. This often occurs during energy minimization when loading large systems onto the GPU.

**Solution:** Increase memory with the `--memory` flag:
```bash
# Default is 3G, increase for larger systems
polyzymd submit -c config.yaml --memory 4G

# For very large systems
polyzymd submit -c config.yaml --memory 8G
```

### "GPU not detected"

Check:
```bash
nvidia-smi  # In job script
```

Make sure the GPU directive is present in the generated script:
- Alpine presets: `#SBATCH --gres=gpu:1`
- Bridges2 preset: `#SBATCH --gpus=v100-32:1` (or your selected GPU type)

### "config.yaml not found"

The self-resubmitting script stores an absolute path to the config file.
Make sure the config file has not been moved or deleted since job submission:

```bash
# Check the path in the generated script
grep CONFIG_PATH job_scripts/r1_300K_LipA.sh
```

If the config was moved, either move it back or regenerate and resubmit:

```bash
polyzymd submit -c /new/path/to/config.yaml --preset aa100 --replicates 1
```
