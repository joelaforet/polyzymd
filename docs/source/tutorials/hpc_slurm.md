# HPC and SLURM Guide

This guide covers running PolyzyMD simulations on HPC clusters using SLURM.

## Overview

Long MD simulations often exceed HPC time limits (typically 24-48 hours). PolyzyMD solves this with **daisy-chaining**: breaking simulations into segments that run as dependent SLURM jobs.

## Directory Structure

PolyzyMD supports separating:

- **Projects directory**: Long-term storage for scripts, logs, configs
- **Scratch directory**: High-performance storage for trajectories

```
/projects/$USER/polyzymd/           # Long-term storage
├── config.yaml
├── structures/
├── job_scripts/                    # Generated SLURM scripts
│   ├── initial_rep1.sh
│   ├── continue_seg1_rep1.sh
│   └── ...
└── slurm_logs/                     # Job output files
    └── MySimulation_seg0_rep1.out

/scratch/alpine/$USER/simulations/  # High-performance storage
├── MyEnzyme_300K_run1/
│   ├── system.pdb
│   ├── equilibration/
│   │   └── trajectory.dcd
│   ├── production_seg0/
│   │   ├── trajectory.dcd
│   │   ├── checkpoint.chk
│   │   └── state_data.csv
│   └── production_seg1/
│       └── ...
└── MyEnzyme_300K_run2/
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
| `aa100` | aa100 | 1x A100 | 24h | 64GB | NVIDIA A100 (recommended) |
| `al40` | al40 | 1x L40 | 24h | 64GB | NVIDIA L40 |
| `blanca-shirts` | blanca-shirts | 1x | 7d | 64GB | Shirts lab partition |
| `testing` | atesting | 1x | 1h | 16GB | Quick tests |

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
cat job_scripts/initial_rep1.sh
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
tail -f slurm_logs/MySimulation_seg0_rep1.out

# Check for errors
grep -i error slurm_logs/*.out
grep -i fail slurm_logs/*.out
```

### Check Simulation Progress

```bash
# List trajectory files
ls -la /scratch/$USER/simulations/*/production_*/trajectory.dcd

# Check trajectory sizes
du -h /scratch/$USER/simulations/*/production_*/*.dcd
```

---

## Handling Failures

### Job Failed Mid-Segment

If a job fails, the dependent jobs won't start. To restart:

1. **Check the error**:
   ```bash
   cat slurm_logs/MySimulation_seg2_rep1.out
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
rm -rf /scratch/$USER/simulations/MySimulation_run1/

# Resubmit
polyzymd submit -c config.yaml --replicates 1 --preset aa100
```

---

## Generated Script Structure

### Initial Script

```bash
#!/bin/bash
#SBATCH --job-name=MySimulation_rep1
#SBATCH --partition=aa100
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --gres=gpu:1
#SBATCH --time=24:00:00
#SBATCH --output=slurm_logs/MySimulation_seg0_rep1.out

# Environment
module load anaconda
conda activate polymerist-env

# Directories
export PROJECTS_DIR="/projects/$USER/polyzymd"
export SCRATCH_DIR="/scratch/$USER/simulations"

# Create output directory
mkdir -p $SCRATCH_DIR/MySimulation_300K_run1

# Run initial simulation (build + equilibration + seg 0)
polyzymd run -c config.yaml --replicate 1
```

### Continuation Script

```bash
#!/bin/bash
#SBATCH --job-name=MySimulation_seg1_rep1
#SBATCH --dependency=afterok:<previous_job_id>
# ... SLURM headers ...

# Continue from previous segment
polyzymd continue \
    --working-dir $SCRATCH_DIR/MySimulation_300K_run1 \
    --segment 1 \
    --segment-time 10.0 \
    --num-samples 250
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
tail -f slurm_logs/*_seg0_*.out
```

### 3. Back Up Important Data

Scratch is often purged. Copy completed simulations to projects:

```bash
# After simulation completes
cp -r /scratch/$USER/simulations/MySimulation_run1 \
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

Add module loading to your job script or check that your `~/.bashrc` loads modules.

### "Out of memory"

Reduce system size:
- Decrease `box.padding`
- Use fewer polymers
- Use smaller production `samples` (fewer frames saved)

### "GPU not detected"

Check:
```bash
nvidia-smi  # In job script
```

Make sure `#SBATCH --gres=gpu:1` is present.
