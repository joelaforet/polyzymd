# How To: Submit Analysis Jobs to a SLURM Cluster

This guide shows you how to submit PolyzyMD analysis computations as SLURM
jobs. It covers submitting a full analysis DAG (replicate → aggregate →
finalize), monitoring progress, and collecting comparison results — all without
running analysis interactively on a login node.

```{note}
This guide covers **analysis** job submission via `polyzymd compare submit`.
For submitting **simulation** jobs, see {doc}`../tutorials/hpc_slurm`.
```

## Before You Start

You need:

- access to a SLURM cluster with `sbatch` available on PATH
- a working pixi environment on the cluster (`pixi install -e build`)
- completed simulation trajectories for at least two conditions
- a `comparison.yaml` that defines your conditions and analysis settings

If you have not yet created a `comparison.yaml`, start with
{doc}`../tutorials/analysis_compare_conditions`.

## What the DAG Looks Like

When you run `polyzymd compare submit`, the framework generates and submits a
directed acyclic graph (DAG) of SLURM jobs:

```text
                ┌─────────────────────┐
                │  Replicate Jobs     │
                │  (one per replicate │
                │   per condition)    │
                └──────────┬──────────┘
                           │ afterany
                ┌──────────▼──────────┐
                │  Aggregate Jobs     │
                │  (one per condition)│
                └──────────┬──────────┘
                           │ afterany
                ┌──────────▼──────────┐
                │  Finalize Job       │
                │  (compare + plot)   │
                └─────────────────────┘
```

**Replicate jobs** each run `compute_replicate()` for one (condition, replicate)
pair. **Aggregate jobs** wait for all replicates of their condition to finish,
then run `aggregate()`. The **finalize job** waits for all aggregate jobs, then
runs the cross-condition comparison and generates plots.

Each job includes automatic retry logic. If a worker exits with a non-zero
code, it requeues itself up to `--max-retries` times (default: 3) before
marking the task as failed.

## The Example Study

This tutorial uses a CALB enzyme study with three conditions:

```text
calb_study/
├── comparison.yaml
├── noPoly_CALB_pNPB/
│   ├── config.yaml
│   └── scratch/
├── SBMA_100_CALB_pNPB/
│   ├── config.yaml
│   └── scratch/
└── EGMA_100_CALB_pNPB/
    ├── config.yaml
    └── scratch/
```

The `comparison.yaml` defines three conditions with three replicates each,
and enables the SASA analysis plugin:

```yaml
name: "calb_polymer_study"
description: "CALB with SBMA and EGMA polymers"
control: "No Polymer"

conditions:
  - label: "No Polymer"
    config: "../noPoly_CALB_pNPB/config.yaml"
    replicates: [1, 2, 3]

  - label: "SBMA-100"
    config: "../SBMA_100_CALB_pNPB/config.yaml"
    replicates: [1, 2, 3]

  - label: "EGMA-100"
    config: "../EGMA_100_CALB_pNPB/config.yaml"
    replicates: [1, 2, 3]

defaults:
  equilibration_time: "10ns"

plugins:
  sasa:
    runs:
      - label: "protein_isolated"
        target_selection: "protein"
        context_selection: "protein"
      - label: "protein_with_polymer"
        target_selection: "protein"
        context_selection: "protein or chainID C"
    probe_radius_nm: 0.14
    n_sphere_points: 960

plot_settings:
  format: "png"
  dpi: 300
  style: "publication"
```

## Step 1: Dry Run

Before submitting real jobs, generate the scripts without sending them to the
scheduler. This lets you inspect the generated SLURM scripts and verify
that paths, partition names, and resource requests are correct.

```bash
pixi run -e build polyzymd compare submit sasa \
    --comparison-yaml comparison.yaml \
    --partition aa100 \
    --mem 8G \
    --time 02:00:00 \
    --dry-run
```

You will see output like:

```text
Submitted 10 jobs (6 replicate + 3 aggregate + 1 finalize)
Dry run only: no jobs were submitted
```

The generated scripts are written to the HPC artifact directory:

```text
comparison/sasa/_hpc/
├── manifest.json          # Snapshot of analysis inputs
├── scripts/
│   ├── replicate__no_polymer__r1.sh
│   ├── replicate__no_polymer__r2.sh
│   ├── replicate__no_polymer__r3.sh
│   ├── replicate__sbma_100__r1.sh
│   ├── ...
│   ├── aggregate__no_polymer.sh
│   ├── aggregate__sbma_100.sh
│   ├── aggregate__egma_100.sh
│   └── finalize.sh
├── logs/
└── status/
    ├── replicates/
    └── conditions/
```

:::{tip}
Open one of the generated `.sh` scripts and check that:

- the `#SBATCH --partition` matches your cluster
- the `pixi` path resolves correctly (override with `--pixi-path` if needed)
- the `--mem` and `--time` values are appropriate for your system size
:::

## Step 2: Submit the DAG

Once you are satisfied with the dry run, submit for real:

```bash
pixi run -e build polyzymd compare submit sasa \
    --comparison-yaml comparison.yaml \
    --partition aa100 \
    --mem 8G \
    --time 02:00:00
```

Expected output:

```text
Submitted 10 jobs (6 replicate + 3 aggregate + 1 finalize)
```

The framework uses SLURM `--dependency=afterany:...` to wire the DAG. Aggregate
jobs will only start after their replicate jobs finish, and the finalize job
will only start after all aggregate jobs complete.

:::{warning}
The `submit` command requires `sbatch` to be available on your PATH. If you are
on a login node where `sbatch` is not available, the command will fail with a
clear error message telling you to use `polyzymd compare run` for local
execution instead.
:::

## Step 3: Monitor Progress

Check the status of your submitted DAG at any time:

```bash
pixi run -e build polyzymd compare status sasa \
    --comparison-yaml comparison.yaml
```

Sample output:

```text
Analysis: sasa
HPC dir: /path/to/calb_study/comparison/sasa/_hpc
States: pending=0 running=2 retrying=0 succeeded=7 failed=0 unknown=0
```

The status reads JSON files that each worker updates atomically. The states
are:

| State | Meaning |
|-------|---------|
| `pending` | Job has not started yet |
| `running` | Worker is currently executing |
| `retrying` | Worker failed but is requeued for another attempt |
| `succeeded` | Worker completed successfully |
| `failed` | Worker exhausted all retries |
| `unknown` | Status file is corrupted or unreadable |

For machine-readable output (useful in scripts), add `--json`:

```bash
pixi run -e build polyzymd compare status sasa \
    --comparison-yaml comparison.yaml --json
```

You can also use standard SLURM tools alongside PolyzyMD status:

```bash
squeue -u $USER
tail -f comparison/sasa/_hpc/logs/*.out
```

## Step 4: Finalize Results

If the finalize job in the DAG succeeded, your comparison results and plots
are already generated. However, there are cases where you may want to run
finalize manually:

- the finalize SLURM job failed but all aggregates succeeded
- you want to re-plot with different settings
- you submitted with `--allow-partial` and want to finalize with available data

Run finalize manually:

```bash
pixi run -e build polyzymd compare finalize sasa \
    --comparison-yaml comparison.yaml
```

Output:

```text
Saved result: /path/to/calb_study/comparison/sasa/sasa_comparison.json
```

The finalize step runs `compare()` (cross-condition statistics) and `plot()`
(figure generation) using the aggregated results that are already on disk.
It does **not** re-run any replicate or aggregate computations.

:::{tip}
If some conditions failed but you still want partial results, pass
`--allow-partial`:

```bash
pixi run -e build polyzymd compare finalize sasa \
    --comparison-yaml comparison.yaml --allow-partial
```

You will see a warning listing the missing conditions, but the comparison
will proceed with whatever data is available.
:::

## Troubleshooting

### Job fails with `pixi: command not found`

SLURM jobs run in a non-interactive shell that may not have your login-time
PATH. Use `--pixi-path` to provide the absolute path:

```bash
pixi run -e build polyzymd compare submit sasa \
    --comparison-yaml comparison.yaml \
    --pixi-path /home/youruser/.pixi/bin/pixi \
    --partition aa100
```

:::{tip}
Find your pixi path with `which pixi` before submitting.
:::

### Job fails with OOM (out of memory)

SASA computation on large systems can be memory-intensive. Increase the
memory allocation:

```bash
pixi run -e build polyzymd compare submit sasa \
    --comparison-yaml comparison.yaml \
    --mem 16G \
    --partition aa100
```

You can also reduce memory pressure by increasing the `chunk_size` in
your SASA plugin settings (which processes fewer frames at once), although
the default of 100 is already conservative.

### Job times out

Increase the wall-time limit:

```bash
pixi run -e build polyzymd compare submit sasa \
    --comparison-yaml comparison.yaml \
    --time 04:00:00 \
    --partition aa100
```

### A replicate is stuck in `retrying`

Check the SLURM log for that replicate:

```bash
cat comparison/sasa/_hpc/logs/replicate__sbma_100__r2.*.out
```

Common causes: trajectory file not found, selection string matches zero atoms,
or a corrupted DCD file. Fix the underlying issue and resubmit the full DAG.

### Manifest/config drift error

If you change `comparison.yaml` or plugin settings after submitting, the
workers will detect a snapshot hash mismatch and raise a `RuntimeError`. This
is a safety check — it prevents inconsistent results from mixing old and new
configurations.

**Fix:** Resubmit the entire DAG with the updated configuration. The old
HPC artifacts will be overwritten.

### Finalize fails with missing aggregated results

This means one or more conditions did not produce an aggregated result file.
Check `polyzymd compare status` to identify the failed conditions, inspect
their logs, fix the issue, and resubmit. Alternatively, use `--allow-partial`
to finalize with available data.

## Resource Configuration Reference

All resource options are passed as CLI flags to `polyzymd compare submit`:

| Flag | Default | Description |
|------|---------|-------------|
| `--partition` | `aa100` | SLURM partition name |
| `--qos` | *(none)* | Quality of service |
| `--account` | *(none)* | SLURM account/allocation |
| `--mem` | `4G` | Memory per job |
| `--time` | `01:00:00` | Wall-time limit (HH:MM:SS or D-HH:MM:SS) |
| `--max-retries` | `3` | Retry count before marking a task as failed |
| `--ntasks` | `1` | SLURM tasks per job |
| `--cpus-per-task` | `1` | CPUs per task |
| `--pixi-path` | `pixi` | Path to pixi executable |
| `--mail-user` | *(none)* | Email for failure notifications |
| `--recompute` | off | Force recomputation even if cached results exist |
| `--allow-partial` | off | Allow finalize with incomplete condition data |
| `--equilibration` | *(from YAML)* | Override equilibration time |
| `--dry-run` | off | Generate scripts without submitting |

## What You Have Now

After following this guide, you have:

- a submitted SLURM DAG that processes replicates in parallel
- the ability to monitor progress without logging into compute nodes
- comparison results and plots generated automatically by the finalize job
- the knowledge to troubleshoot common failure modes

## See Also

- {doc}`../tutorials/sasa_analysis` — Tutorial for configuring and interpreting SASA analysis
- {doc}`../tutorials/hpc_slurm` — Submitting simulation jobs to SLURM
- {doc}`../tutorials/analysis_compare_conditions` — Setting up comparison.yaml
- {doc}`../tutorials/analysis_complete_workflow` — Full local analysis workflow
