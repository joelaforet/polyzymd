# Run PolyzyMD on SLURM Clusters

Use this guide when you already have a working `config.yaml` and want the
shortest path to a reliable SLURM submission workflow.

PolyzyMD generates self-resubmitting job scripts. Each replicate runs one
segment, checks whether more work remains, and resubmits itself when needed.
That lets long simulations continue across wall-time limits without requiring
manual dependency chains.

## Before you start

- validate your config locally first
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
pixi run -e build polyzymd submit \
    -c config.yaml \
    --preset aa100 \
    --pixi-env auto \
    --replicates 1 \
    --generate-only
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
pixi run -e build polyzymd submit \
    -c config.yaml \
    --preset testing \
    --pixi-env auto \
    --time-limit 0:05:00 \
    --replicates 1
```

This is the fastest way to catch bad paths, scheduler issues, or environment
problems.

## Step 4: submit your real run

Once the short test succeeds, submit production jobs:

```bash
pixi run -e build polyzymd submit \
    -c config.yaml \
    --preset aa100 \
    --pixi-env auto \
    --replicates 1-5 \
    --email your.email@university.edu
```

Useful variants:

```bash
# Override storage locations
pixi run -e build polyzymd submit \
    -c config.yaml \
    --preset aa100 \
    --pixi-env auto \
    --projects-dir /projects/$USER/polyzymd \
    --scratch-dir /scratch/alpine/$USER/polyzymd_sims

# Give a larger system more RAM
pixi run -e build polyzymd submit \
    -c config.yaml \
    --preset aa100 \
    --pixi-env auto \
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
pixi run -e build polyzymd status -c config.yaml
pixi run -e build polyzymd check-progress -c config.yaml -r 1
```

## Recover a stalled replicate

If a replicate stops progressing, inspect it first:

```bash
pixi run -e build polyzymd recover -c config.yaml -r 1
```

If the report shows unfinished work, resubmit a recovery job:

```bash
pixi run -e build polyzymd recover \
    -c config.yaml \
    -r 1 \
    --submit \
    --preset aa100 \
    --pixi-env auto
```

## OpenMM runtime policy

Run `polyzymd submit` from the `build` environment. The `--pixi-env` option
selects the environment that the Slurm job will use. It does not select the
environment that runs the submit command.

For a known site, `--pixi-env auto` resolves to the site policy when PolyzyMD
generates the script. Blanca uses `sim-cuda-12-4`. Bridges-2 uses
`sim-cuda-12-6`. A newer driver does not change this selection.

On the allocated node, PolyzyMD checks the driver and activates the selected
environment. It then creates an explicit CUDA Context, calculates an energy,
and runs one integration step. This test detects an unusable CUDA runtime or
PTX compiler before molecular setup. PolyzyMD does not fall back to CPU.

If the node is not compatible, PolyzyMD submits a replacement job and excludes
that node. It stops after three failed routing attempts.

PolyzyMD records the selected runtime. A replica cannot change its Pixi
environment, OpenMM version, platform, or precision during resubmission. For
supported runtimes and instructions for new hardware, see
{doc}`hardware_platforms`.

## Bridges-2

Use the `bridges2` preset to request PSC Bridges-2 scheduler resources. The
preset resolves `auto` to the checked-in `sim-cuda-12-6` site environment:

```bash
pixi run -e build polyzymd submit \
    -c config.yaml \
    --preset bridges2 \
    --account abc123_gpu \
    --pixi-env auto \
    --replicates 1-3
```

Common Bridges-2 differences:

- you may need `--account` if you want to charge a specific allocation
- GPU selection can be adjusted with `--gpu-type`

The allocated node must be able to create a CUDA Context with this environment.
It does not select a different environment.

## CU Boulder Alpine and Blanca

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
pixi run -e build polyzymd submit \
    -c config.yaml \
    --preset aa100 \
    --account ucb625_asc1 \
    --pixi-env auto \
    --replicates 1-5
```

Blanca example (partition, account, and QoS are typically the same value):

```bash
pixi run -e build polyzymd submit \
    -c config.yaml \
    --preset blanca-shirts \
    --pixi-env auto \
    --replicates 1-5
```

If you run GROMACS on Blanca, use `--constraint` to request hardware that is
compatible with the site GROMACS build. See {doc}`gromacs_export`.

:::{tip}
If you are also running analysis jobs via `polyzymd compare submit-all`, see
{doc}`hpc_execution` for detailed CU Boulder cluster configuration including
partition tables and troubleshooting.
:::

## What the generated scripts do

Each generated NVIDIA OpenMM script follows the same loop:

1. detect GPU capability and validate the environment selected at submission
2. activate the environment, create an explicit CUDA Context, calculate an
   energy, and run one integration step
3. run `polyzymd run-segment`
4. call `polyzymd check-progress`
5. resubmit itself if work remains

This loop lets a simulation continue across wall-time limits.

Version 1.3 OpenMM SLURM scripts require an NVIDIA GPU. The Python simulation
path can use explicit CPU and OpenCL platforms, but these platforms need a
site-managed batch wrapper. See {doc}`hardware_platforms`.

For OpenMM jobs, `SIGTERM` or `SIGUSR1` immediately submits one successor with
an `afterany` dependency on the current job, then forwards the signal to
`run-segment`. A receipt in the replicate directory prevents duplicate traps
and normal exit handling from submitting a second successor. If `sbatch`
fails, the job exits with an error and prints the manual recovery command.

OpenMM records each minimization and equilibration phase in an atomic
`phase.json`. Only `status: completed` permits a successor to skip a phase.
A checkpoint without that record is incomplete: PolyzyMD resumes from a
synchronized portable state when one exists, or restarts the current phase
from the preceding completed phase. Temperature ramps record the exact step
and scheduled temperature. Minimization is restarted when interrupted because
OpenMM minimization itself has no portable mid-minimization resume point.

Every dynamics loop first calibrates with at most 1,000 steps and then targets
five seconds per `Simulation.step()` call. Reporter frequency does not
determine how quickly Python notices a preemption signal.

For GROMACS jobs the scripts additionally:

- run EM, equilibration stages, and production with checkpoint restart
- pass `-maxh` so GROMACS exits cleanly before the wall-time limit
- trap SIGTERM and forward it to `gmx mdrun`, which flushes a checkpoint
- self-resubmit until the full production duration completes

## Submit GROMACS jobs

### CPU GROMACS

```bash
pixi run -e build polyzymd submit \
    -c config.yaml \
    --engine gromacs \
    --preset aa100 \
    --replicates 1-3
```

### GPU GROMACS

```bash
pixi run -e build polyzymd submit \
    -c config.yaml \
    --engine gromacs \
    --preset blanca-shirts \
    --constraint "A40|A100" \
    --replicates 1-3
```

GROMACS uses the site module or container from `config.yaml`. For CPU and GPU
configuration, MPI settings, constraints, and recovery details, see
{doc}`gromacs_export`.

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
- hardware portability and extension: {doc}`hardware_platforms`

<!-- IMAGE OPPORTUNITY: Add a simple lifecycle diagram showing `submit ->
run-segment -> check-progress -> resubmit`, plus a second annotated screenshot
of a generated SLURM script header. -->
# Build integrity and recovery

For prebuilt OpenMM campaigns, wait for `polyzymd build` to finish before
submitting simulation jobs. A successful build contains these three files in
each replicate directory:

- `solvated_system.pdb`
- `system.xml`
- `build_manifest.json`

The manifest is the final commit marker. `--skip-build` verifies its
configuration hash, artifact hashes, and particle count before OpenMM creates a
simulation. A missing or malformed manifest means the build is incomplete.
Older campaigns without a manifest remain recoverable only when the topology
and System particle counts agree; PolyzyMD emits a warning and does not silently
create a manifest.

Build and simulation processes share `.polyzymd.lock` in the replicate
directory. A second process exits instead of writing concurrently. PolyzyMD
also refuses to rebuild a directory after progress, minimization,
equilibration, or production artifacts exist. Use a new output directory for a
new molecular system.

Continuation loads `production_N_topology.pdb` from the predecessor production
segment together with its System and State. The root `solvated_system.pdb` is a
warned legacy fallback only. If a diagnostic names different particle counts,
restore all artifacts from the same build; do not copy individual files until
the counts happen to match.
