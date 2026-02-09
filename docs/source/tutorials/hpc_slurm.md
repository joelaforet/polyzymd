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
| `aa100` | aa100 | 1x A100 | 24h | 3GB | NVIDIA A100 (recommended) |
| `al40` | al40 | 1x L40 | 24h | 3GB | NVIDIA L40 |
| `blanca-shirts` | blanca-shirts | 1x | 24h | 3GB | Shirts lab partition |
| `testing` | atesting_a100 | 1x | 6min | 3GB | Quick tests |

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

## Daisy-Chain Workflow

### How It Works

1. **Initial job**: Builds system, runs equilibration, runs first production segment
2. **Continuation jobs**: Load checkpoint, run next segment
3. **Dependencies**: Each job depends on the previous one completing successfully

```
Job 1 (initial)     Job 2 (continue)    Job 3 (continue)
┌─────────────┐     ┌─────────────┐     ┌─────────────┐
│ Build       │     │ Load chkpt  │     │ Load chkpt  │
│ Equilibrate │ --> │ Run seg 1   │ --> │ Run seg 2   │
│ Run seg 0   │     │ Save chkpt  │     │ Save chkpt  │
└─────────────┘     └─────────────┘     └─────────────┘
```

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

### Job Failed Mid-Segment

If a job fails, the dependent jobs won't start. To restart:

1. **Check the error**:
   ```bash
   cat slurm_logs/s2_r1_300K_*.out
   ```

2. **Fix the issue** (if possible)

3. **Manually continue**:
   ```bash
   # Edit and resubmit the continuation script
   sbatch job_scripts/continue_seg2_rep1.sh
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

### Initial Script

The initial job script (segment 0) builds the system, runs equilibration, and runs the first production segment:

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

# Exit immediately if any command fails
set -e

module purge
module load miniforge
mamba activate polymerist-env

# Projects directory (scripts, configs, logs)
PROJECTS_DIR="/projects/$USER/polyzymd/my_simulation"

# Scratch directory (simulation output)
SCRATCH_DIR="/scratch/alpine/$USER/polyzymd_sims/LipA_300K_run1"

# Ensure scratch directory exists
mkdir -p "$SCRATCH_DIR"

# Change to projects directory where config lives
cd "$PROJECTS_DIR"

echo "Starting initial simulation segment 0"
echo "Projects dir: $PROJECTS_DIR"
echo "Scratch dir: $SCRATCH_DIR"
echo "Config: config.yaml"
echo "Replicate: 1"
echo "Timestamp: $(date)"

# Run the initial simulation using polyzymd CLI
polyzymd run -c "config.yaml" \
    --replicate 1 \
    --scratch-dir "$SCRATCH_DIR" \
    --segment-time 10.0 \
    --segment-frames 250

echo "Segment 0 completed successfully at $(date)"
```

### Continuation Script

Continuation scripts load the checkpoint from the previous segment and continue the simulation:

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

# Exit immediately if any command fails
set -e

module purge
module load miniforge
mamba activate polymerist-env

# Projects directory (scripts, configs, logs)
PROJECTS_DIR="/projects/$USER/polyzymd/my_simulation"

# Scratch directory (simulation output - where previous segment data lives)
SCRATCH_DIR="/scratch/alpine/$USER/polyzymd_sims/LipA_300K_run1"

# Change to projects directory
cd "$PROJECTS_DIR"

echo "Starting continuation segment 1"
echo "Projects dir: $PROJECTS_DIR"
echo "Scratch dir: $SCRATCH_DIR"
echo "Timestamp: $(date)"

# Continue simulation from previous segment using polyzymd CLI
polyzymd continue \
    -w "$SCRATCH_DIR" \
    -s 1 \
    -t 10.0 \
    -n 250

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
- `--preset PRESET` - SLURM preset: aa100, al40, blanca-shirts, testing. Default: aa100
- `--scratch-dir PATH` - Override scratch directory for simulation output
- `--projects-dir PATH` - Override projects directory for scripts/logs
- `--output-dir PATH` - Directory for job scripts. Default: {projects_dir}/job_scripts
- `--email EMAIL` - Email for job notifications
- `--time-limit TIME` - Override SLURM time limit (HH:MM:SS)
- `--memory SIZE` - Override SLURM memory allocation (e.g., "4G", "8G"). Default: 3G
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

Make sure `#SBATCH --gres=gpu:1` is present.

### "config.yaml not found"

Make sure you're running `polyzymd submit` from the directory containing your `config.yaml`, or use an absolute path:

```bash
polyzymd submit -c /full/path/to/config.yaml --preset aa100
```
