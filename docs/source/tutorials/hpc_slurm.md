# HPC and SLURM Guide

This guide covers running PolyzyMD simulations on HPC clusters using SLURM.

## Overview

Long MD simulations often exceed HPC time limits (typically 24-48 hours). PolyzyMD solves this with **daisy-chaining**: breaking simulations into segments that run as dependent SLURM jobs.

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
cat job_scripts/initial_seg0_rep1.sh

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
│   │   ├── initial_seg0_rep1.sh
│   │   ├── continue_seg1_rep1.sh
│   │   └── ...
│   └── slurm_logs/                 # Job output files
│       └── s0_r1_300K_LipA.out

/scratch/alpine/$USER/polyzymd_sims/  # High-performance storage
├── LipA_Substrate_EGPMA-SBMA_10pct_300K_run1/
│   ├── system.pdb
│   ├── equilibration/
│   │   └── trajectory.dcd
│   ├── production_seg0/
│   │   ├── trajectory.dcd
│   │   ├── checkpoint.chk
│   │   └── state_data.csv
│   └── production_seg1/
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
cat job_scripts/initial_seg0_rep1.sh | head -20
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
│   │   ├── initial_seg0_rep1.sh
│   │   └── ...
│   └── slurm_logs/

/local/scratch/$USER/polyzymd_sims/          # High-performance local scratch
├── LipA_Substrate_300K_run1/
│   ├── system.pdb
│   ├── equilibration/
│   └── production_seg0/
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

## Daisy-Chain Workflow

### How It Works

1. **Initial job**: Builds system, runs equilibration, runs first production segment
2. **Continuation jobs**: Load checkpoint, run next segment
3. **Dependencies**: Each job depends on the previous one via `afterany` (not `afterok`), so the chain continues even if a segment is interrupted
4. **Recovery**: Each continuation job checks whether the previous segment completed normally or was interrupted, and automatically recovers if needed

```
Job 1 (initial)     Job 2 (continue)    Job 3 (continue)
┌─────────────┐     ┌─────────────┐     ┌─────────────┐
│ Build       │     │ Check prev  │     │ Check prev  │
│ Equilibrate │ --> │ Recover?    │ --> │ Recover?    │
│ Run seg 0   │     │ Run seg 1   │     │ Run seg 2   │
└─────────────┘     └─────────────┘     └─────────────┘
     afterany            afterany
```

Using `afterany` instead of `afterok` means that if a segment exits with a non-zero code (e.g. 99 for graceful shutdown), the next segment still starts. The recovery preamble in each continuation script decides whether to proceed, recover, or abort. See {ref}`smart-restart` for details.

### Configuring Segments

```yaml
simulation_phases:
  production:
    duration: 100.0    # 100 ns total
  segments: 10         # 10 segments of 10 ns each
```

```{tip}
**Segment duration** = total duration / segments

Choose segment duration to fit within your cluster's time limit with margin:
- 24h limit → ~8-10 ns segments (2h GPU time + overhead)
- 48h limit → ~20 ns segments
```

---

(smart-restart)=
## Smart Restart & Fault Tolerance

When you run 60 daisy chains with 10 segments each, that is 600 SLURM jobs.
With vanilla `afterok` dependencies, a single failed segment kills the entire
downstream chain — SLURM marks every subsequent job as
`DependencyNeverSatisfied` and you lose the remaining wall-time allocation.

PolyzyMD's **smart restart** system handles this automatically.  The
generated scripts already include everything described below — you do not
need to configure anything.  This section explains what happens under the
hood so you can debug issues or adapt the approach to other workflows.

### What Happens Automatically

Every generated SLURM script (both initial and continuation) includes three
pieces of fault-tolerance infrastructure:

1. **Signal handling** — Python-side handlers catch SIGUSR1 (wall-time
   warning) and SIGTERM (preemption), save an emergency checkpoint, and
   exit with code 99.
2. **Signal forwarding** — Bash trap + background + wait pattern forwards
   signals from the SLURM batch shell to the Python child process.
3. **Recovery preamble** — Each continuation script checks whether the
   previous segment completed, was interrupted (recoverable), or crashed
   (unrecoverable), and takes the appropriate action before starting its
   own segment.

The dependency between segments is `afterany` (not `afterok`), so the
chain continues even when a segment exits with a non-zero code.

### The Three Scenarios

| Scenario | Signal | What Happens | Outcome |
|----------|--------|--------------|---------|
| **Wall-time warning** | `SIGUSR1` (5 min before limit) | Emergency state saved, exit 99 | Next job recovers automatically |
| **Preemption** | `SIGTERM` (120 s grace on Blanca) | Emergency state saved, exit 99 | Next job recovers automatically |
| **Hard crash** | None (OOM, segfault, node failure) | No state saved | Next job detects missing state, exits with error |

```{note}
The wall-time signal is configured via `#SBATCH --signal=B:USR1@300`, which
tells SLURM to send `SIGUSR1` to the batch shell 300 seconds (5 minutes)
before the time limit expires.  This gives the simulation enough time to
save a full OpenMM state (~10-30 seconds on GPU).
```

### Emergency Checkpoint Files

When an interrupt is detected, the signal handler writes three files into
the current segment's directory (e.g. `production_3/`):

| File | Purpose |
|------|---------|
| `emergency_state.xml` | Portable OpenMM state (positions, velocities, forces) |
| `emergency_system.xml` | Serialized OpenMM System (force field parameters) |
| `INTERRUPTED` | Marker file with step-count metadata for recovery |

The `INTERRUPTED` marker contains the information needed for recovery:

```
segment_index=3
steps_completed=1250000
total_steps=2500000
remaining_steps=1250000
```

### Recovery and the Resume Subdirectory

When a continuation job starts and finds an `INTERRUPTED` marker in the
previous segment's directory, it runs `polyzymd recover` before starting
its own segment.  Recovery works as follows:

1. Parse the `INTERRUPTED` marker for step counts
2. Load `emergency_state.xml` and `emergency_system.xml`
3. Run the remaining steps, writing trajectory to a **resume subdirectory**
   (`production_N_resume/`) to avoid corrupting the partial DCD file
4. Write the normal end-of-segment files (`state.xml`, `system.xml`,
   `checkpoint.chk`) to the original `production_N/` directory
5. Remove the `INTERRUPTED` marker

```{tip}
The resume subdirectory strategy means that the partial trajectory from the
interrupted run and the completion trajectory from recovery are in separate
files.  During post-processing, concatenate them:

    production_3/production_3_trajectory.dcd        # first half
    production_3_resume/production_3_resume_trajectory.dcd  # second half
```

After recovery completes, the segment directory looks like any normally-
completed segment to the `ContinuationManager`, so the next segment
proceeds as usual.

### Signal Forwarding: Why trap + background + wait?

SLURM sends signals to the **batch shell process**, not to child processes.
Bash ignores `SIGUSR1` by default, so without explicit forwarding, the
Python simulation never sees the signal.  The generated scripts use a
standard pattern to solve this:

```bash
# Background the Python process
polyzymd continue -w "$SCRATCH_DIR" -s 3 -t 10.0 -n 250 &
CHILD_PID=$!

# Trap signals and forward them to the child
trap 'kill -USR1 $CHILD_PID' USR1
trap 'kill -TERM $CHILD_PID' TERM

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
generated scripts.  Without them, signals will not reach the Python process
and graceful shutdown will not work.
```

### Recovery Preamble Logic

Each continuation script includes a recovery preamble that runs before
the segment's own simulation.  The logic is:

```bash
PREV_DIR="${SCRATCH_DIR}/production_${PREV_SEG}"

if [ -f "${PREV_DIR}/INTERRUPTED" ]; then
    # Previous segment was interrupted — run recovery
    polyzymd recover -w "$SCRATCH_DIR" -s "$PREV_SEG"
elif [ -f "${PREV_DIR}/production_${PREV_SEG}_state.xml" ]; then
    # Previous segment completed normally — proceed
    :
else
    # No state and no marker — unrecoverable crash
    exit 1
fi
```

### Manually Triggering an Interrupt

You can test graceful shutdown or manually stop a running segment by
sending `SIGUSR1` via `scancel`:

```bash
# Send USR1 to a specific job
scancel --signal=USR1 <job_id>

# The job will save emergency state and exit with code 99
# The next segment in the chain will recover automatically
```

This is useful when you realize a simulation has a problem and want to
stop it cleanly without losing progress.

### Manual Recovery Commands

If you need to recover segments outside the automatic daisy-chain workflow
(e.g. the chain has already ended), use the CLI commands directly:

```bash
# Recover a single interrupted segment
polyzymd recover -w /scratch/user/sim/LipA_300K_run1 -s 3

# Preview what would be recovered (no simulation)
polyzymd recover -w /scratch/user/sim/LipA_300K_run1 -s 3 --dry-run

# Scan an entire chain and report status
polyzymd recover-chain -w /scratch/user/sim/LipA_300K_run1 --dry-run

# Scan and recover all interrupted segments
polyzymd recover-chain -w /scratch/user/sim/LipA_300K_run1
```

See the {ref}`CLI Reference <cli-recover>` section below for full option
details.

```{note}
Recovery itself is signal-aware.  If a recovery run is interrupted (e.g.
the recovery job also hits a wall-time limit), the `INTERRUPTED` marker is
updated with the new step count and the next attempt picks up where recovery
left off.
```

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
cat job_scripts/initial_seg0_rep1.sh
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
tail -f slurm_logs/s0_r1_300K_LipA*.out

# Check for errors
grep -i error slurm_logs/*.out
grep -i fail slurm_logs/*.out
```

### Check Simulation Progress

```bash
# List trajectory files
ls -la /scratch/$USER/polyzymd_sims/*/production_*/trajectory.dcd

# Check trajectory sizes
du -h /scratch/$USER/polyzymd_sims/*/production_*/*.dcd
```

---

## Handling Failures

Most failures are handled automatically by the {ref}`smart restart <smart-restart>`
system.  This section covers the cases where manual intervention is needed.

### Wall-Time or Preemption Interrupts

**No action needed.** The smart restart system saves an emergency checkpoint
and the next segment in the chain recovers automatically.  Check the SLURM
log to confirm:

```bash
# Look for the graceful shutdown message
grep -i "interrupted\|graceful\|recovery" slurm_logs/s2_r1_300K_*.out
```

### Hard Crash (OOM, Segfault, Node Failure)

If a job crashes without saving state, the next segment will detect the
missing `state.xml` and exit with an error.  To diagnose and recover:

1. **Check the error**:
   ```bash
   cat slurm_logs/s2_r1_300K_*.out
   ```

2. **Fix the issue** (e.g. increase memory with `--memory 8G`)

3. **Re-run the failed segment manually**, then resubmit the rest of
   the chain:
   ```bash
   # Re-run the failed segment
   polyzymd continue -w /scratch/$USER/sim/LipA_300K_run1 -s 2 -t 10.0 -n 250

   # Resubmit remaining continuation scripts
   sbatch job_scripts/continue_seg3_rep1.sh
   ```

### Checking Chain Health

Use `recover-chain` to get a quick status report across all segments:

```bash
polyzymd recover-chain -w /scratch/$USER/sim/LipA_300K_run1 --dry-run
```

This prints each segment's status (COMPLETED, INTERRUPTED, or MISSING)
without modifying anything.

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

### Initial Script

The initial job script (segment 0) builds the system, runs equilibration, and runs the first production segment.  It includes signal forwarding and exit-code handling for the {ref}`smart restart <smart-restart>` system:

```bash
#!/bin/bash
#SBATCH --partition=aa100
#SBATCH --job-name=i_s0_r1_300K_LipA
#SBATCH --output=slurm_logs/s0_r1_300K_LipA.%A_%a.out
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

module purge 2>/dev/null || true
ml miniforge 2>/dev/null || true

eval "$(conda shell.bash hook)"
conda activate polyzymd-env

set -e

export INTERCHANGE_EXPERIMENTAL=1

PROJECTS_DIR="/projects/$USER/polyzymd/my_simulation"
SCRATCH_DIR="/scratch/alpine/$USER/polyzymd_sims/LipA_300K_run1"
mkdir -p "$SCRATCH_DIR"
cd "$PROJECTS_DIR"

echo "Starting initial simulation segment 0"
echo "Timestamp: $(date)"

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

# Run simulation in background for signal forwarding
polyzymd run -c "config.yaml" \
    --replicate 1 \
    --scratch-dir "$SCRATCH_DIR" \
    --segment-time 10.0 \
    --segment-frames 250 &
CHILD_PID=$!

# Wait loop (handles signal interruption of wait)
set +e
wait "$CHILD_PID" 2>/dev/null
RC=$?
while kill -0 "$CHILD_PID" 2>/dev/null; do
    wait "$CHILD_PID" 2>/dev/null
    RC=$?
done
set -e

if [ $RC -eq 99 ]; then
    echo "Segment 0 interrupted (graceful shutdown) at $(date)"
    exit 99
elif [ $RC -ne 0 ]; then
    echo "Segment 0 FAILED with exit code $RC at $(date)"
    exit $RC
fi

echo "Segment 0 completed successfully at $(date)"
```

### Continuation Script

Continuation scripts include a **recovery preamble** that checks the previous segment's status before starting.  If the previous segment was interrupted, `polyzymd recover` runs automatically:

```bash
#!/bin/bash
#SBATCH --partition=aa100
#SBATCH --job-name=c_s1_r1_300K_LipA
#SBATCH --output=slurm_logs/s1_r1_300K_LipA.%A_%a.out
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

module purge 2>/dev/null || true
ml miniforge 2>/dev/null || true

eval "$(conda shell.bash hook)"
conda activate polyzymd-env

export INTERCHANGE_EXPERIMENTAL=1

PROJECTS_DIR="/projects/$USER/polyzymd/my_simulation"
SCRATCH_DIR="/scratch/alpine/$USER/polyzymd_sims/LipA_300K_run1"
cd "$PROJECTS_DIR"

echo "Starting continuation segment 1"
echo "Timestamp: $(date)"

# Recovery preamble: check previous segment status
PREV_SEG=$(( 1 - 1 ))
PREV_DIR="${SCRATCH_DIR}/production_${PREV_SEG}"

if [ -f "${PREV_DIR}/INTERRUPTED" ]; then
    echo "Previous segment $PREV_SEG was interrupted — running recovery"
    polyzymd recover -w "$SCRATCH_DIR" -s "$PREV_SEG"
    RECOVER_RC=$?
    if [ $RECOVER_RC -ne 0 ]; then
        echo "Recovery of segment $PREV_SEG FAILED (exit code $RECOVER_RC)"
        exit 1
    fi
    echo "Recovery of segment $PREV_SEG completed successfully"
elif [ -f "${PREV_DIR}/production_${PREV_SEG}_state.xml" ]; then
    echo "Previous segment $PREV_SEG completed normally — proceeding"
else
    echo "ERROR: Previous segment $PREV_SEG has no state.xml and no INTERRUPTED marker"
    echo "This segment cannot continue — the previous job likely crashed."
    exit 1
fi

set -e

# Signal forwarding (same pattern as initial script)
CHILD_PID=""
forward_signal() {
    if [ -n "$CHILD_PID" ] && kill -0 "$CHILD_PID" 2>/dev/null; then
        echo "Forwarding $1 to Python process (PID $CHILD_PID)"
        kill -"$1" "$CHILD_PID"
    fi
}
trap 'forward_signal USR1' USR1
trap 'forward_signal TERM' TERM

polyzymd continue \
    -w "$SCRATCH_DIR" \
    -s 1 \
    -t 10.0 \
    -n 250 &
CHILD_PID=$!

set +e
wait "$CHILD_PID" 2>/dev/null
RC=$?
while kill -0 "$CHILD_PID" 2>/dev/null; do
    wait "$CHILD_PID" 2>/dev/null
    RC=$?
done
set -e

if [ $RC -eq 99 ]; then
    echo "Segment 1 interrupted (graceful shutdown) at $(date)"
    exit 99
elif [ $RC -ne 0 ]; then
    echo "Segment 1 FAILED with exit code $RC at $(date)"
    exit $RC
fi

echo "Segment 1 completed successfully at $(date)"
```

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
tail -f slurm_logs/*_s0_*.out
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

| Cluster Time Limit | Recommended Segment Duration |
|--------------------|------------------------------|
| 1 hour (testing) | 0.5 - 1 ns |
| 24 hours | 8 - 12 ns |
| 48 hours | 20 - 30 ns |
| 7 days | 50 - 100 ns |

---

## CLI Reference

### `polyzymd submit`

Submit daisy-chain simulation jobs to SLURM.

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

Continue a simulation from a previous segment checkpoint.

```bash
polyzymd continue -w WORKING_DIR -s SEGMENT -t TIME [OPTIONS]
```

**Required:**
- `-w, --working-dir PATH` - Working directory with previous segment
- `-s, --segment INT` - Segment index to run (1-based)
- `-t, --segment-time FLOAT` - Duration of this segment (ns)

**Options:**
- `-n, --num-samples INT` - Number of frames to save. Default: 250

(cli-recover)=
### `polyzymd recover`

Recover an interrupted simulation segment.  Reads the `INTERRUPTED` marker
and emergency state files, runs the remaining steps in a
`production_N_resume/` subdirectory, then writes end-of-segment files so
the next `ContinuationManager` can proceed.

```bash
polyzymd recover -w WORKING_DIR -s SEGMENT [--dry-run]
```

**Required:**
- `-w, --working-dir PATH` - Working directory containing simulation output
- `-s, --segment INT` - Segment index that was interrupted

**Options:**
- `--dry-run` - Show what would be recovered without running the simulation

**Example:**
```bash
# Preview recovery
polyzymd recover -w /scratch/$USER/sim/LipA_300K_run1 -s 3 --dry-run

# Run recovery
polyzymd recover -w /scratch/$USER/sim/LipA_300K_run1 -s 3
```

### `polyzymd recover-chain`

Scan all segments in a working directory and report their status.  Without
`--dry-run`, also recovers any interrupted segments.

```bash
polyzymd recover-chain -w WORKING_DIR [--dry-run]
```

**Required:**
- `-w, --working-dir PATH` - Working directory containing simulation output

**Options:**
- `--dry-run` - Show chain status without taking action

**Example:**
```bash
# Status report only
polyzymd recover-chain -w /scratch/$USER/sim/LipA_300K_run1 --dry-run

# Recover all interrupted segments
polyzymd recover-chain -w /scratch/$USER/sim/LipA_300K_run1
```

**Example output:**
```
Scanning /scratch/user/sim/LipA_300K_run1
Found 10 segment(s):

  segment   0: COMPLETED
  segment   1: COMPLETED
  segment   2: COMPLETED
  segment   3: INTERRUPTED (50% done, 1250000 remaining)
  segment   4: MISSING (no state.xml, no INTERRUPTED marker — likely crashed)
  ...

Summary: 3 completed, 1 interrupted, 1 missing
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
- `Dependency` - Waiting for previous job

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

Make sure you're running `polyzymd submit` from the directory containing your `config.yaml`, or use an absolute path:

```bash
polyzymd submit -c /full/path/to/config.yaml --preset aa100
```
