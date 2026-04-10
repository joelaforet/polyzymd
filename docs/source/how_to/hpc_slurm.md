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

If you are still setting up the project itself, start with {doc}`../tutorials/quickstart`.

## Step 1: validate and dry-run locally

From the repository root or a subdirectory under it:

```bash
pixi run -e build polyzymd validate -c config.yaml
pixi run -e cuda-12-4 polyzymd submit -c config.yaml --preset aa100 --replicates 1 --generate-only
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
| `aa100` | CU Boulder Alpine A100 | main production runs |
| `al40` | CU Boulder Alpine L40 | production runs on L40 nodes |
| `blanca-shirts` | CU Boulder Blanca | preemptable or lab-specific runs |
| `testing` | short queue | smoke tests only |
| `bridges2` | PSC Bridges2 | Bridges2 GPU jobs |

Use `testing` first when you are verifying a new system or workflow.

## Step 3: submit one small test job

Run a short job before launching many replicates:

```bash
pixi run -e cuda-12-4 polyzymd submit \
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
pixi run -e cuda-12-4 polyzymd submit \
    -c config.yaml \
    --preset aa100 \
    --replicates 1-5 \
    --email your.email@university.edu
```

Useful variants:

```bash
# Override storage locations
pixi run -e cuda-12-4 polyzymd submit \
    -c config.yaml \
    --preset aa100 \
    --projects-dir /projects/$USER/polyzymd \
    --scratch-dir /scratch/alpine/$USER/polyzymd_sims

# Give a larger system more RAM
pixi run -e cuda-12-4 polyzymd submit \
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
pixi run -e cuda-12-4 polyzymd status -c config.yaml
pixi run -e cuda-12-4 polyzymd check-progress -c config.yaml -r 1
```

## Recover a stalled replicate

If a replicate stops progressing, inspect it first:

```bash
pixi run -e cuda-12-4 polyzymd recover -c config.yaml -r 1
```

If the report shows unfinished work, resubmit a recovery job:

```bash
pixi run -e cuda-12-4 polyzymd recover -c config.yaml -r 1 --submit --preset aa100
```

## Cluster-specific note for Bridges2

Use the `bridges2` preset when running on PSC Bridges2:

```bash
pixi run -e cuda-12-6 polyzymd submit \
    -c config.yaml \
    --preset bridges2 \
    --account abc123_gpu \
    --replicates 1-3
```

Common Bridges2 differences:

- it uses the `cuda-12-6` environment
- you may need `--account` if you want to charge a specific allocation
- GPU selection can be adjusted with `--gpu-type`

## What the generated scripts do

Each generated script follows the same loop:

1. activate the selected `pixi` environment
2. run `polyzymd run-segment`
3. call `polyzymd check-progress`
4. resubmit itself if work remains

That is why long runs can continue automatically after wall-time expiry or a
graceful interruption.

## Common fixes

### `pixi: command not found`

Make sure `pixi` is available in non-interactive shells, not only your login
shell setup.

### job dies with OOM

Increase `--memory`, reduce system size, or test with fewer polymers.

### config path no longer exists

The generated script stores the config path it was given at submission time. If
you move the config, regenerate the scripts and resubmit.

### need to stop a job permanently

Because standard cancellation can trigger graceful restart behavior, use:

```bash
scancel --signal=KILL <job_id>
```

## Related reference pages

- command details: {doc}`../reference/cli_reference`
- configuration fields: {doc}`../reference/configuration`
- first-run setup: {doc}`../tutorials/quickstart`

<!-- IMAGE OPPORTUNITY: Add a simple lifecycle diagram showing `submit ->
run-segment -> check-progress -> resubmit`, plus a second annotated screenshot
of a generated SLURM script header. -->
