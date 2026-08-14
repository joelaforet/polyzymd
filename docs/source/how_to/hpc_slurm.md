# Run PolyzyMD on SLURM Clusters

Use this guide when you already have a working `config.yaml` and want the
shortest path to a reliable SLURM submission workflow.

PolyzyMD generates self-resubmitting job scripts. Each replicate runs one
segment, checks whether more work remains, and resubmits itself when needed.
That lets long simulations continue across wall-time limits without requiring
manual dependency chains.

## Before you start

- validate your config locally first
- know which `pixi` CUDA environment matches your cluster
- know which SLURM preset you want to use

If you are still setting up the project itself, start with {doc}`../get_started/quickstart`.

:::{admonition} Use compute resources, not login nodes
:class: important

Validation and SLURM script generation are lightweight. System builds and local
simulation commands can require substantial RAM, CPU/GPU time, and scratch I/O.
On shared HPC systems, submit jobs to compute nodes or use an interactive
compute allocation; do not run heavy build or simulation commands directly on a
login node.
:::

## Step 1: validate and dry-run locally

From the repository root or a subdirectory under it:

```bash
pixi run -e build polyzymd validate -c config.yaml
pixi run -e sim-cuda-12-4 polyzymd submit -c config.yaml --preset aa100 --replicates 1 --generate-only
```

The `--generate-only` flag creates a script in `job_scripts/` without submitting it,
so you can inspect it before launching real jobs.

:::{versionchanged} 1.3.0
`--dry-run` is now preview-only (no files written, no submission). Use
`--generate-only` to generate SLURM scripts without submitting — this is
the behavior that `--dry-run` had in earlier versions. The two flags are
mutually exclusive.
:::

## Step 2: pick a preset

PolyzyMD includes presets for common clusters:

| Preset | Cluster style | Typical use |
|--------|---------------|-------------|
| `aa100` | NVIDIA A100 partition | main production runs |
| `al40` | NVIDIA L40 partition | production runs on L40 nodes |
| `blanca-shirts` | Blanca condo partition | preemptable or condo runs |
| `testing` | short queue | smoke tests only |
| `bridges2` | PSC Bridges2 | Bridges2 GPU jobs |

Use `testing` first when you are verifying a new system or workflow.

## Step 3: submit one small test job

Run a short job before launching many replicates:

```bash
pixi run -e sim-cuda-12-4 polyzymd submit \
    -c config.yaml \
    --preset testing \
    --time-limit 0:05:00 \
    --replicates 1
```

This is the fastest way to catch bad paths, scheduler issues, or environment
problems.

## Step 4: submit your real run

Once the short test succeeds, submit production jobs:

```bash
pixi run -e sim-cuda-12-4 polyzymd submit \
    -c config.yaml \
    --preset aa100 \
    --replicates 1-5 \
    --email your.email@university.edu
```

Useful variants:

```bash
# Override storage locations
pixi run -e sim-cuda-12-4 polyzymd submit \
    -c config.yaml \
    --preset aa100 \
    --projects-dir /projects/$USER/polyzymd \
    --scratch-dir /scratch/alpine/$USER/polyzymd_sims

# Give a larger system more RAM
pixi run -e sim-cuda-12-4 polyzymd submit \
    -c config.yaml \
    --preset aa100 \
    --memory 8G
```

## Monitor jobs

Use normal SLURM tools for the scheduler view:

```bash
squeue -u $USER
scontrol show job <job_id>
tail -f slurm_logs/*.out
```

Use PolyzyMD for simulation progress:

```bash
pixi run -e sim-cuda-12-4 polyzymd status -c config.yaml
pixi run -e sim-cuda-12-4 polyzymd check-progress -c config.yaml -r 1
```

## Recover a stalled replicate

If a replicate stops progressing, inspect it first:

```bash
pixi run -e sim-cuda-12-4 polyzymd recover -c config.yaml -r 1
```

If the report shows unfinished work, resubmit a recovery job:

```bash
pixi run -e sim-cuda-12-4 polyzymd recover -c config.yaml -r 1 --submit --preset aa100
```

## Cluster-specific note for Bridges2

Use the `bridges2` preset when running on PSC Bridges2:

```bash
pixi run -e sim-cuda-12-6 polyzymd submit \
    -c config.yaml \
    --preset bridges2 \
    --account abc123_gpu \
    --replicates 1-3
```

Common Bridges2 differences:

- it uses the `sim-cuda-12-6` environment
- you may need `--account` if you want to charge a specific allocation
- GPU selection can be adjusted with `--gpu-type`

## Cluster-specific note for CU Boulder (Alpine and Blanca)

OpenMM jobs should normally use `--pixi-env auto`. The generated script runs
`nvidia-smi` on the allocated node before activating Pixi, chooses the highest
compatible checked-in CUDA 12.0, 12.4, or 12.6 environment, then creates a
minimal explicit CUDA Context. Missing GPUs, unknown drivers, malformed probe
output, incompatible overrides, and CUDA Context failures stop the job. There
is no automatic CPU fallback.

The job atomically records the environment, OpenMM version, CUDA runtime,
platform, precision, driver, device, and compute capability in
`runtime_platform.json`. Resubmission refuses to change environment, OpenMM
version, platform, or precision within a replica.

When temporarily using the old fixed `--pixi-env sim-cuda-12-4` workflow,
exclude `bgpu-g4-u20,bgpu-g4-u24`. Auto routing does not rely on this node
list and removes that workaround once those nodes report a supported runtime.

The repository does not contain claimed Blanca probe results yet. Release
acceptance still requires short allocations on representative CUDA 12.0 and
12.4 nodes, successful explicit CUDA Context creation, a deterministic short
benchmark, and the documented cross-environment trajectory/energy comparison.
Record node, GPU, driver, maximum CUDA, compute capability, OpenMM Context
result, and numerical tolerances from those real jobs; do not infer them from
partition or node names.

CU Boulder runs two SLURM clusters. Switch between them with environment
modules before submitting:

```bash
ml slurm/alpine   # shared campus resource
ml slurm/blanca   # PI-owned condo nodes
```

:::{important}
You must run `module load slurm/blanca` (or `ml slurm/blanca`) before
`sbatch` to see Blanca partitions. Without it, Blanca queues are invisible
to the scheduler.
:::

Both clusters require `--partition`, `--account`, and `--qos` explicitly.
Alpine example:

```bash
pixi run -e sim-cuda-12-4 polyzymd submit \
    -c config.yaml \
    --preset aa100 \
    --account ucb625_asc1 \
    --replicates 1-5
```

Blanca example (partition, account, and QoS are typically the same value):

```bash
pixi run -e sim-cuda-12-4 polyzymd submit \
    -c config.yaml \
    --preset blanca-shirts \
    --replicates 1-5
```

Blanca has 9+ different GPU types (P100, T4, V100, RTX 6000, A40, A100, L40,
H100, RTX Pro 6000). If you are running GPU GROMACS on Blanca, use
`--constraint` to pin your job to a compatible architecture — see
[GROMACS engine](#gromacs-engine) below.

The `blanca-shirts` preset uses `qos=preemptable`, which means jobs can be
preempted by the node owner. GROMACS scripts handle this gracefully via
SIGTERM trapping (see [Preemption resilience](#preemption-resilience)).

:::{tip}
If you are also running analysis jobs via `polyzymd compare submit-all`, see
{doc}`hpc_execution` for detailed CU Boulder cluster configuration including
partition tables and troubleshooting.
:::

## What the generated scripts do

Each generated script follows the same loop:

1. detect GPU capability, select or validate the CUDA environment, activate it,
   and create an explicit CUDA Context
2. run `polyzymd run-segment`
3. call `polyzymd check-progress`
4. resubmit itself if work remains

That is why long runs can continue automatically after wall-time expiry or a
graceful interruption.

For GROMACS jobs the scripts additionally:

- run EM, equilibration stages, and production with checkpoint restart
- pass `-maxh` so GROMACS exits cleanly before the wall-time limit
- trap SIGTERM and forward it to `gmx mdrun`, which flushes a checkpoint
- self-resubmit until the full production duration completes

## GROMACS engine

:::{versionchanged} 1.3.0
`polyzymd submit` now supports `--engine gromacs` for GROMACS SLURM
submission with the same self-resubmitting workflow used by OpenMM.
:::

:::{seealso}
For the complete GROMACS HPC guide including cluster recipes, flag glossary,
config reference, and troubleshooting, see {doc}`gromacs_export`.
:::

### CPU GROMACS

```bash
pixi run -e build polyzymd submit \
    -c config.yaml \
    --engine gromacs \
    --preset aa100 \
    --replicates 1-3
```

:::{note}
GROMACS CPU jobs use `pixi run -e build` (not `sim-cuda-12-4`) because the
GROMACS binary comes from `module load`, not from the pixi environment.
:::

### GPU GROMACS

```bash
pixi run -e build polyzymd submit \
    -c config.yaml \
    --engine gromacs \
    --preset blanca-shirts \
    --constraint "A40|A100" \
    --replicates 1-3
```

### GPU constraints (`--constraint`)

GROMACS uses ahead-of-time compiled CUDA kernels. Unlike OpenMM, which
JIT-compiles kernels at launch, a GROMACS binary compiled for one GPU
architecture may not run on another. If your cluster has mixed GPU types,
`--constraint` ensures your job lands on compatible hardware:

```bash
--constraint "A40"              # single GPU type
--constraint "A40|A100"         # either type (OR)
--constraint "avx2&rh8"         # feature AND (CPU + OS flags)
```

This maps directly to the SLURM `#SBATCH --constraint` directive and works
on any cluster — it is not specific to CU Boulder.

(preemption-resilience)=
### Preemption resilience

GROMACS SLURM scripts trap `SIGTERM` (the signal SLURM sends before
preempting a job). When the trap fires, the script:

1. forwards the signal to `gmx mdrun`, which flushes a `.cpt` checkpoint
2. waits for GROMACS to exit
3. resubmits the job so production resumes from the checkpoint

Combined with `--constraint`, this ensures resumed jobs land on compatible
GPU hardware. This is especially important on clusters with preemptable QoS
(e.g., Blanca `qos=preemptable`).

### Module loading (`gromacs.module_load`)

GROMACS on HPC clusters typically requires loading prerequisite modules
(compiler, MPI) before the GROMACS module itself. Use the `gromacs.module_load`
config field — it is inserted verbatim into the generated SLURM script:

```yaml
gromacs:
  module_load: "module load gcc/11.2.0 openmpi/4.1.1 gromacs/2024.2"
```

List prerequisites before the GROMACS module so dependencies resolve in order.

## Common fixes

### `pixi: command not found`

Make sure `pixi` is available in non-interactive shells, not only your login
shell setup.

### job dies with OOM

Increase `--memory`, reduce system size, or test with fewer polymers.

### config path no longer exists

The generated script stores the config path it was given at submission time. If
you move the config, regenerate the scripts and resubmit.

(hpc-slurm-stop-permanently)=
### need to stop a job permanently

Because standard cancellation can trigger graceful restart behavior, use:

```bash
scancel --signal=KILL <job_id>
```

## Related reference pages

- command details: {doc}`../reference/cli_reference`
- configuration fields: {doc}`../reference/configuration`
- GROMACS HPC guide: {doc}`gromacs_export`
- first-run setup: {doc}`../get_started/quickstart`

<!-- IMAGE OPPORTUNITY: Add a simple lifecycle diagram showing `submit ->
run-segment -> check-progress -> resubmit`, plus a second annotated screenshot
of a generated SLURM script header. -->
