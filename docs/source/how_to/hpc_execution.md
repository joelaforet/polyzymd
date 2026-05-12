# How To: Submit Analysis Jobs to a SLURM Cluster

This guide shows you how to submit PolyzyMD analysis computations as SLURM
jobs. It covers submitting a full analysis DAG (replicate → aggregate →
finalize), monitoring progress, cross-plugin dependency ordering, and
collecting comparison results — all without running analysis interactively on a
login node.

```{note}
This guide covers **analysis** job submission via `polyzymd compare submit`.
For submitting **simulation** jobs, see {doc}`hpc_slurm`.
```

## Before You Start

You need:

- access to a SLURM cluster with `sbatch` available on PATH
- a working pixi environment on the cluster (`pixi install -e build`)
- completed simulation trajectories for at least two conditions
- a `comparison.yaml` that defines your conditions and analysis settings

If you have not yet created a `comparison.yaml`, start with
{doc}`analysis_compare_conditions`.

:::{admonition} Use compute resources, not login nodes
:class: important

Trajectory-analysis jobs can require substantial RAM, CPU/GPU time, and scratch
I/O. On shared HPC systems, run `polyzymd compare submit` from a place where it
can submit jobs to compute nodes, and reserve interactive `polyzymd compare run`
for an allocated compute session or suitably provisioned workstation. Do not run
heavy analysis commands directly on a login node.
:::

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

**Replicate jobs** each run the runner-backed replicate stage for one
(condition, replicate) pair. **Aggregate jobs** wait for all replicates of their
condition to finish, then run `aggregate()`. The **finalize job** waits for all
aggregate jobs, then runs the cross-condition comparison and generates plots.

Each job includes automatic retry logic. If a worker exits with a non-zero
code, it requeues itself up to `--max-retries` times (default: 3) before
marking the task as failed.

## Submitting all enabled analyses with dependency ordering

Use `compare submit-all` to submit every enabled plugin from `comparison.yaml`
in dependency order:

```bash
pixi run -e build polyzymd compare submit-all \
    -f comparison.yaml \
    --partition aa100 \
    --qos normal
```

This command:

- discovers enabled plugins from `plugins:`
- orders plugins by declared `dependencies`
- submits each plugin DAG with cross-plugin finalize dependencies

Exclude one or more analyses with repeatable `--exclude`:

```bash
pixi run -e build polyzymd compare submit-all \
    -f comparison.yaml \
    --exclude sasa \
    --partition aa100
```

Use `--dry-run` to generate all scripts and print the submission summary table
without dispatching jobs.

## Framework finalize-only mode

The submission framework supports plugins that do not implement
compute/aggregate stages and only run comparison/plot logic. For those plugins,
the manifest pipeline mode is `finalize_only`, and submission creates a single
finalize job.

The stable v1.3 analysis plugins use compute and/or aggregate stages before
finalization. Treat `finalize_only` as a framework capability for future or
custom compare-only plugins, not as the normal path for the stable plugins
listed in this guide.

This behavior applies to both `compare submit` and `compare submit-all`.

## The Example Study

This guide uses a CALB enzyme study with three conditions:

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
    config: "noPoly_CALB_pNPB/config.yaml"
    replicates: [1, 2, 3]

  - label: "SBMA-100"
    config: "SBMA_100_CALB_pNPB/config.yaml"
    replicates: [1, 2, 3]

  - label: "EGMA-100"
    config: "EGMA_100_CALB_pNPB/config.yaml"
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

:::{note}
In this guide, `--dry-run` refers to `polyzymd compare submit --dry-run`,
which generates scripts without submitting. This is different from
`polyzymd submit --dry-run` (simulation submission), which in v1.3.0 is
preview-only and writes nothing. See the {doc}`/reference/cli_reference` for details.
:::

Before submitting real jobs, generate the scripts without sending them to the
scheduler. This lets you inspect the generated SLURM scripts and verify
that paths, partition names, and resource requests are correct.

```bash
pixi run -e build polyzymd compare submit sasa \
    -f comparison.yaml \
    --partition aa100 \
    --mem 8G \
    --time 02:00:00 \
    --dry-run
```

You will see output like:

```text
Would submit 13 jobs (9 replicate + 3 aggregate + 1 finalize)
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
    -f comparison.yaml \
    --partition aa100 \
    --mem 8G \
    --time 02:00:00
```

Expected output:

```text
Submitted 13 jobs (9 replicate + 3 aggregate + 1 finalize)
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
    -f comparison.yaml
```

Sample output:

```text
Analysis: sasa
HPC dir: /path/to/calb_study/comparison/sasa/_hpc
States: pending=4 running=2 retrying=0 succeeded=7 failed=0 unknown=0
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
    -f comparison.yaml --json
```

### Reconciling status with SLURM

If workers were preempted or cancelled externally, the local status files may
be stale. Use `--reconcile` to query `sacct` and update status files
atomically:

```bash
pixi run -e build polyzymd compare status sasa \
    -f comparison.yaml --reconcile
```

This checks all `pending`, `running`, and `retrying` tasks against SLURM's
accounting database and updates any that have reached a terminal state.

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
    -f comparison.yaml
```

Output:

```text
Saved result: /path/to/calb_study/comparison/sasa/result.json
```

The finalize step runs `compare()` (cross-condition statistics) and `plot()`
(figure generation) using the aggregated results that are already on disk.
It does **not** re-run any replicate or aggregate computations.

:::{tip}
If some conditions failed but you still want partial results, pass
`--allow-partial`:

```bash
pixi run -e build polyzymd compare finalize sasa \
    -f comparison.yaml --allow-partial
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
    -f comparison.yaml \
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
    -f comparison.yaml \
    --mem 16G \
    --partition aa100
```

You can also reduce memory pressure by increasing the `stride` in
your SASA plugin settings (which analyzes fewer frames), or by
decreasing `chunk_size` (which processes fewer frames per batch
at the cost of more I/O overhead). The default `chunk_size` of 100
is already conservative for most systems.

### Job times out

Increase the wall-time limit:

```bash
pixi run -e build polyzymd compare submit sasa \
    -f comparison.yaml \
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

```{note}
If the configured control condition is intentionally filtered out by a plugin's
`filter_conditions()` (for example, polymer-dependent analyses removing a
no-polymer control), finalize now auto-switches to all-vs-all comparison
without requiring `--allow-partial`.
```

For true runtime data loss (failed/missing condition outputs), strict behavior
still applies and you must pass `--allow-partial`.

### MDAnalysis ChainReader errors (F13)

If a worker fails with MDAnalysis ChainReader exceptions, this is usually an
external reader/data issue rather than a PolyzyMD logic bug.

Recommended checks:

- validate trajectory files are complete and readable
- verify topology/trajectory atom ordering consistency
- re-run with a single replicate to isolate the broken input file

### SLURM association and job-limit failures (for example MaxJobs=0, F15)

If `sbatch` rejects jobs due to association, account, or quota limits, inspect
your scheduler associations:

```bash
sacctmgr show association user=$USER
```

Then resubmit with explicit account/QoS flags as required by your site.

### QoS/account errors

Many clusters require both partition and QoS/account policy settings.

- set `--qos <name>` when jobs remain pending or are rejected
- set `--account <allocation>` when your scheduler enforces account routing

PolyzyMD prints a submit-time tip when `--partition` is set but `--qos` is not.

### Job preempted on a shared partition

If your cluster uses preemption, jobs may be killed and requeued
automatically. PolyzyMD workers handle this via `#SBATCH --requeue` and
built-in retry logic. If a worker exhausts all retries due to repeated
preemption:

1. Increase `--max-retries` (default: 3) when submitting.
2. Use `--reconcile` with `polyzymd compare status` to sync stale status files.
3. Target a non-preemptible partition or QoS if available.

### Why `afterany` instead of `afterok`?

The DAG uses `--dependency=afterany:...` rather than `afterok`. This means
aggregate and finalize jobs **always start** regardless of whether upstream
jobs succeeded or failed. This is intentional:

- Aggregate jobs are designed to handle partial replicate data gracefully.
- The finalize job can produce partial comparison results via `--allow-partial`.
- Using `afterok` would silently block downstream jobs forever when a single
  replicate fails, leaving the pipeline in a zombie state with no output.

## Resource Configuration Reference

All resource options are passed as CLI flags to `polyzymd compare submit`:

| Flag | Default | Description |
|------|---------|-------------|
| `--partition` | *(none)* | SLURM partition name (optional; if omitted, cluster default partition is used) |
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
| `--job-arrays` | off | Submit one SLURM array job per condition instead of individual jobs |

### Resource precedence and plugin hints

For `--mem`, `--time`, and `--cpus-per-task`, submission precedence is:

1. explicit CLI flag
2. plugin `slurm_resource_hint`
3. system default

Current plugin resource hints:

- `sasa`: `8G`, `02:00:00`
- `secondary_structure`: `16G`
- `hydrogen_bonds`: `16G`

This means large plugins get safer defaults, while explicit CLI requests still
override plugin hints.

### Using Job Arrays

For studies with many replicates per condition, job arrays reduce the number
of `sbatch` calls and make the SLURM queue easier to manage. Pass
`--job-arrays` to submit one array job per condition instead of individual
replicate jobs:

```bash
pixi run -e build polyzymd compare submit sasa \
    -f comparison.yaml \
    --partition aa100 \
    --mem 8G \
    --time 02:00:00 \
    --job-arrays
```

With job arrays, the DAG structure changes:

- **Without `--job-arrays`:** One `sbatch` call per replicate (e.g. 9 calls for 3×3).
- **With `--job-arrays`:** One array `sbatch` call per condition (e.g. 3 calls for 3 conditions), each containing replicate tasks as array elements.

Aggregate jobs depend on their condition's array job completing (`afterany`),
and the finalize job depends on all aggregate jobs. Log files for array tasks
include the array task ID in their filename.

:::{note}
Job arrays are opt-in because some clusters have restrictive array size
limits or non-standard array scheduling behavior. Test with `--dry-run
--job-arrays` first to inspect the generated scripts.
:::

## CU Boulder HPC: Alpine and Blanca

CU Boulder operates two SLURM clusters: **Alpine** (shared campus resource) and
**Blanca** (condo model with PI-owned nodes). Both clusters require
`--partition`, `--account`, and `--qos` to be set explicitly for job
submission.

### Switching between clusters

Use environment modules to select which cluster's SLURM scheduler you target:

```bash
# Target Blanca (PI-owned condo nodes)
module load slurm/blanca

# Target Alpine (shared campus resource)
module load slurm/alpine
```

Run the appropriate `module load slurm/<cluster>` command **before** any
`polyzymd` submission command. The module swap updates `sbatch`, `squeue`, and
other SLURM utilities to point at the selected cluster.

### Required SLURM flags

Both clusters require all three scheduling flags. Omitting any of them causes
`sbatch` to reject the job.

| Flag | Alpine (shared) | Blanca (condo) |
|------|-----------------|----------------|
| `--partition` | `amilan` (CPU), `aa100` or `ami100` (GPU) | `blanca-<group>` (e.g. `blanca-shirts`) |
| `--account` | Your allocation (e.g. `ucb625_asc1`) | Same as partition (e.g. `blanca-shirts`) |
| `--qos` | `normal` | Same as partition (e.g. `blanca-shirts`) |

:::{warning}
If you omit `--partition` on Blanca, `sbatch` fails with
*"A partition has not been provided"* — and the error message unhelpfully
references Alpine documentation. This does not mean you need to switch to
Alpine. Add `--partition=blanca-<group>` and resubmit.
:::

### Example: submit all analyses on Alpine

```bash
module load slurm/alpine

pixi run -e build polyzymd compare submit-all \
    -f comparison.yaml \
    --partition amilan \
    --account ucb625_asc1 \
    --qos normal
```

### Example: submit all analyses on Blanca

On Blanca, the partition, account, and QoS are typically the same value —
your PI's condo allocation name:

```bash
module load slurm/blanca

pixi run -e build polyzymd compare submit-all \
    -f comparison.yaml \
    --partition blanca-shirts \
    --account blanca-shirts \
    --qos blanca-shirts
```

### Example: submit a single analysis on Blanca

The same flags work with `compare submit` for individual plugins:

```bash
pixi run -e build polyzymd compare submit sasa \
    -f comparison.yaml \
    --partition blanca-shirts \
    --account blanca-shirts \
    --qos blanca-shirts \
    --mem 8G \
    --time 02:00:00
```

:::{tip}
If you are unsure which accounts and partitions you have access to, run:

```bash
sacctmgr show association user=$USER format=account,partition,qos
```

This lists every account/partition/QoS combination available to your user.
:::

## What You Have Now

After following this guide, you have:

- a submitted SLURM DAG that processes replicates in parallel
- the ability to monitor progress without logging into compute nodes
- comparison results and plots generated automatically by the finalize job
- the knowledge to troubleshoot common failure modes

## See Also

- {doc}`../tutorials/sasa_analysis` — Tutorial for configuring and interpreting SASA analysis
- {doc}`hpc_slurm` — Submitting simulation jobs to SLURM
- {doc}`analysis_compare_conditions` — Setting up comparison.yaml
- {doc}`../tutorials/analysis_complete_workflow` — Full local analysis workflow
