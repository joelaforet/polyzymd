# Run simulations on SLURM clusters

Submit each replicate as a chain of SLURM jobs. Each job runs one
{term}`segment` of the simulation. It then checks whether work remains, and
submits its successor if it does. The chain continues across wall-time limits
without manual job dependencies.

## Before you start

- Validate your config with `polyzymd validate -c config.yaml`.
- Choose a SLURM preset (see Step 2).

If you do not have a config yet, do {doc}`../get_started/quickstart` first.

:::{admonition} Use compute resources, not login nodes
:class: important

Validation and script generation are light. A system build or a local
simulation needs a lot of memory, CPU or GPU time, and scratch I/O. On a shared
cluster, run these commands in a batch job or an interactive compute
allocation. Do not run them on a login node.
:::

Run `polyzymd submit` from the `build` environment. Prefix each command with
`pixi run -e build`, or activate the environment once with `pixi shell -e build`.

## Step 1: write the job script and read it

```bash
pixi run -e build polyzymd validate -c config.yaml
pixi run -e build polyzymd submit \
    -c config.yaml \
    --preset aa100 \
    --pixi-env auto \
    --replicates 1 \
    --generate-only
```

`--generate-only` writes the script to `job_scripts/` and does not submit it.
Read the `#SBATCH` lines before you submit real jobs.

`--dry-run` prints the submission plan. It writes no file and submits nothing.
You cannot use `--dry-run` and `--generate-only` together.

## Step 2: choose a preset

A preset sets the partition, QoS, account and time limit of each job.

| Preset | Partition | QoS | Time limit | Use |
|--------|-----------|-----|------------|-----|
| `aa100` | `aa100` | `normal` | 23:59:59 | NVIDIA A100 nodes (CU Boulder Alpine) |
| `al40` | `al40` | `normal` | 23:59:59 | NVIDIA L40 nodes (CU Boulder Alpine) |
| `blanca-shirts` | `blanca,blanca-shirts` | `preemptable` | 23:59:59 | CU Boulder Blanca condo nodes |
| `blanca-chbe-rdi` | `blanca,blanca-chbe-rdi` | `preemptable` | 23:59:59 | CU Boulder Blanca condo nodes |
| `bridges2` | `GPU-shared` | none | 24:00:00 | PSC Bridges-2 GPU nodes |
| `testing` | `atesting_a100` | `testing` | 0:05:59 | Short tests on CU Boulder Alpine |

The Alpine presets (`aa100`, `al40` and `testing`) set a CU Boulder
allocation as the account. On another allocation, give your own with
`--account <account>`. To use another cluster, override the preset's fields
with `--partition`, `--account`, `--qos` and `--gpu-type`. See
{doc}`hardware_platforms`.

Use `testing` first when you try a new system or a new workflow.

## Step 3: submit one short test job

Run a short job before you submit many replicates:

```bash
pixi run -e build polyzymd submit \
    -c config.yaml \
    --preset testing \
    --pixi-env auto \
    --time-limit 0:05:00 \
    --replicates 1
```

A short job finds a wrong path, a scheduler problem or a broken environment
in minutes.

## Step 4: submit the production replicates

```bash
pixi run -e build polyzymd submit \
    -c config.yaml \
    --preset aa100 \
    --pixi-env auto \
    --replicates 1-5 \
    --email your.email@university.edu
```

To write the output to other folders, give `--projects-dir` and
`--scratch-dir`:

```bash
pixi run -e build polyzymd submit \
    -c config.yaml \
    --preset aa100 \
    --pixi-env auto \
    --projects-dir /projects/$USER/polyzymd \
    --scratch-dir /scratch/alpine/$USER/polyzymd_sims
```

To give a large system more memory, set `--memory`:

```bash
pixi run -e build polyzymd submit \
    -c config.yaml \
    --preset aa100 \
    --pixi-env auto \
    --memory 8G
```

## Monitor the jobs

Use the SLURM tools to see the scheduler state:

```bash
squeue -u $USER
scontrol show job <job_id>
tail -f slurm_logs/*.out
```

Use PolyzyMD to see the simulation progress:

```bash
pixi run -e build polyzymd status -c config.yaml
pixi run -e build polyzymd check-progress -c config.yaml -r 1
```

`polyzymd status --format agent` prints the SLURM state, the throughput, the
time to completion and the reason a chain stopped, for many configs at once.
See {doc}`monitor_simulations`.

## Recover a stalled replicate

Show what remains of replicate 1:

```bash
pixi run -e build polyzymd recover -c config.yaml -r 1
```

If work remains, submit a recovery job:

```bash
pixi run -e build polyzymd recover \
    -c config.yaml \
    -r 1 \
    --submit \
    --preset aa100 \
    --pixi-env auto
```

(hpc-slurm-stop-a-chain)=
## Stop a chain

`scancel` alone does not stop an OpenMM chain. SLURM sends `SIGTERM` to the
job. `run-segment` then exits with code 99. The job script reads code 99 as
"interrupted, work remains" and submits a successor within seconds. Each new
attempt also makes a new `production_N` folder.

Stop the chain with `polyzymd cancel`:

```bash
pixi run -e build polyzymd cancel -c config.yaml -r 1-3
```

The command does two things:

1. It writes a `STOP` file into each replicate folder.
2. It cancels the queued and running jobs of those replicates.

The job script checks for `STOP` at the start of every job and before it
submits a successor. A successor that is already queued exits before it starts
a segment.

To let the current segment finish and then stop, keep the running job:

```bash
pixi run -e build polyzymd cancel -c config.yaml -r 1-3 --stop-only
```

To start again, remove the marker and submit again:

```bash
pixi run -e build polyzymd cancel -c config.yaml -r 1-3 --resume
pixi run -e build polyzymd submit -c config.yaml -r 1-3 --preset blanca-shirts
```

`STOP` is a text file at `<replicate folder>/STOP`. It names who stopped the
chain, when, and how to undo the stop. To delete it by hand has the same effect
as `--resume`. The job script also stops when the job environment sets
`POLYZYMD_STOP_CHAIN=1`. `POLYZYMD_STOP_FILE=<path>` moves the marker to
another path.

A chain uses the job script that `submit` wrote, and OpenMM and GROMACS job
scripts both check for `STOP`. A GROMACS script written by PolyzyMD 1.2 or
older, and an OpenMM script written before `polyzymd cancel` existed, do
not. Stop such a chain with `SIGKILL`, which the script cannot trap:

```bash
scancel --batch --signal=KILL <job_id>
```

Then cancel any successor that is already queued (`squeue -u $USER`).

## Build integrity and recovery

For a campaign with a separate build step, wait until `polyzymd build` ends
before you submit simulation jobs. A complete build has these files in each
replicate folder:

- `solvated_system.pdb`
- `system.prmtop`, the analysis topology. If ParmEd cannot convert the
  system, the file is missing and the build log says so.
- `system.xml`
- `build_manifest.json`

The build writes `build_manifest.json` last. Before OpenMM makes a simulation,
`--skip-build` checks the manifest's config hash, its file hashes and its
particle count. A missing or damaged manifest means that the build did not
finish. A replicate folder without a manifest can still continue when the
particle counts of the topology and the System agree. PolyzyMD then logs a
warning. It does not write a manifest for that folder.

The build and the simulation share the lock file `.polyzymd.lock` in the
replicate folder. A second process exits and does not write. PolyzyMD refuses
to build again in a folder that has progress, minimization, equilibration or
production files. Use a new output folder for a new molecular system.

A continuation loads `production_N_topology.pdb`, the System and the State of
the previous segment. If the segment has no `production_N_topology.pdb`, the
continuation stops with an error. If an error
names two different particle counts, restore all files from the same build. Do
not copy single files until the counts agree.

## OpenMM runtime policy

`--pixi-env` selects the environment that the SLURM job uses. It does not
select the environment that runs `submit`.

For a known site, `--pixi-env auto` takes the site's environment when
PolyzyMD writes the script. Blanca uses `sim-cuda-12-4`. Bridges-2 uses
`sim-cuda-12-6`. A newer driver does not change this choice.

On the allocated node, the job checks the driver and activates the
environment. It then makes an explicit CUDA Context, calculates an energy and
runs one integration step. This test finds an unusable CUDA runtime or PTX
compiler before the molecular setup. PolyzyMD does not fall back to the CPU.

If the node is not compatible, the job submits a replacement job that
excludes the node. After three failed attempts, the chain stops. The preset
excludes nodes that are known to be incompatible, so these nodes never use an
attempt. See [Excluded Blanca GPU nodes](#excluded-blanca-gpu-nodes).

PolyzyMD records the runtime of each replicate. A replicate cannot change its
pixi environment, OpenMM version, platform or precision when it is submitted
again. For the supported runtimes and for new hardware, see
{doc}`hardware_platforms`.

## Bridges-2

The `bridges2` preset requests PSC Bridges-2 resources. It resolves `auto` to
the `sim-cuda-12-6` environment:

```bash
pixi run -e build polyzymd submit \
    -c config.yaml \
    --preset bridges2 \
    --account <allocation> \
    --pixi-env auto \
    --replicates 1-3
```

- Give `--account` to charge a specific allocation.
- Change the GPU type with `--gpu-type`. The default is `v100-32`.
- The preset sets no `--mem`, because Bridges-2 allocates memory per GPU.

The allocated node must make a CUDA Context with this environment. The job
does not select another environment.

(cu-boulder-site-notes)=
## CU Boulder Alpine and Blanca

CU Boulder runs two SLURM clusters. Load the module of the cluster before you
submit:

```bash
ml slurm/alpine   # shared campus cluster
ml slurm/blanca   # condo nodes owned by research groups
```

:::{important}
Run `ml slurm/blanca` before `sbatch` or `polyzymd submit` to use Blanca. If
you do not, SLURM does not show the Blanca partitions.
:::

Both clusters need a partition, an account and a QoS. The presets set all
three. Alpine example, with your own allocation:

```bash
pixi run -e build polyzymd submit \
    -c config.yaml \
    --preset aa100 \
    --account <account> \
    --pixi-env auto \
    --replicates 1-5
```

On Blanca, the partition, the account and the QoS usually have the same name:

```bash
pixi run -e build polyzymd submit \
    -c config.yaml \
    --preset blanca-shirts \
    --pixi-env auto \
    --replicates 1-5
```

For GROMACS on Blanca, request compatible hardware with `--constraint`. See
{doc}`gromacs_export`.

(excluded-blanca-gpu-nodes)=
### Excluded Blanca GPU nodes

The `blanca-shirts` and `blanca-chbe-rdi` presets give SLURM an `--exclude`
list. Jobs do not start on these nodes:

| Node | Reason |
|------|--------|
| `bgpu-bortz1` | Unreliable node |
| `bgpu-g4-u20` | NVIDIA driver 525.147. The pinned `sim-cuda-12-4` builds need driver 550 or newer. |
| `bgpu-g4-u24` | The same old driver as `bgpu-g4-u20` |

Without the list, a job on `bgpu-g4-u20` or `bgpu-g4-u24` logs
`ROUTING: sim-cuda-12-4 is incompatible with driver 525.147`. It then uses one
of its three routing attempts and submits again. A chain that gets both nodes
in a row uses all its attempts and stops with
`FATAL: CUDA routing failed after 3 retries`.

A running chain keeps the exclude list of its current script until its next
submission.

To use a different list, give `--exclude`. For example, test a node again after
CURC upgrades its driver:

```bash
pixi run -e build polyzymd submit \
    -c config.yaml \
    --preset blanca-shirts \
    --exclude bgpu-bortz1 \
    --replicates 1
```

`--exclude` replaces the preset's list. It does not add to it. Give
`--exclude ""` to exclude no node. Leave out the option to keep the preset's
list.

### GROMACS on Blanca

A GROMACS acceptance test ran on 14 August 2026 with the `blanca-shirts`
account and QoS. Both short jobs used GROMACS 2024.2 from the site module.

| Mode | Node | Hardware | Result |
|------|------|----------|--------|
| CPU | `bgpu-shirts3` | 2 CPU threads | 20-step test completed |
| GPU | `bgpu-shirts1` | NVIDIA A40, driver 550.90.07 | 20-step GPU nonbonded test completed |

The site module reported CUDA GPU support. These results apply only to the
tested module and nodes. Run the short test again after a module or driver
update.

:::{tip}
To run analyses on the cluster, with one job per condition and replicate, see
{doc}`hpc_execution`.
:::

## What the job scripts do

Each OpenMM job script for an NVIDIA GPU does these steps:

1. It exits if a `STOP` file exists (see [Stop a chain](#hpc-slurm-stop-a-chain)).
2. It reads the GPU's capability and checks the environment selected at
   submission.
3. It activates the environment, makes an explicit CUDA Context, calculates an
   energy and runs one integration step.
4. It runs `polyzymd run-segment`.
5. It runs `polyzymd check-progress`.
6. It submits itself again if work remains and no `STOP` file exists.

The OpenMM job scripts need an NVIDIA GPU. A local run with
`polyzymd run` can use the `CPU` or `OpenCL` platform. See
{doc}`hardware_platforms`.

On `SIGTERM` or `SIGUSR1`, an OpenMM job submits one successor at once, with an
`afterany` dependency on the current job. It then passes the signal to
`run-segment`. A receipt file in the replicate folder stops a second trap, or
the normal exit, from submitting a second successor. If `sbatch` fails, the job
exits with an error and prints the command to recover by hand.

The routing state passes from job to job. A rerouted successor inherits
`POLYZYMD_ROUTING_RETRY_COUNT` and `POLYZYMD_ROUTING_FAILED_NODES`. Its
`--exclude` is the failed nodes plus the preset's list. After a segment that
ran, the successor gets a count of zero and an empty node list.

OpenMM records each minimization and equilibration phase in `phase.json`,
which it writes atomically. A successor skips a phase only when the record
says `status: completed`. A checkpoint without that record is incomplete.
PolyzyMD then resumes from a synchronized portable state if one exists.
Otherwise, it starts the phase again from the end of the previous completed
phase. A temperature ramp records the step and the scheduled temperature.
An interrupted minimization starts again, because OpenMM cannot resume a
minimization part way through.

Each dynamics loop first calibrates with at most 1,000 steps. It then runs
about five seconds of steps per `Simulation.step()` call. Python therefore
sees a preemption signal within seconds, whatever the reporter interval.

The GROMACS job scripts do these steps:

- They run minimization, the equilibration stages and production, and restart
  each from its checkpoint.
- They pass `-maxh` to `gmx mdrun`, so GROMACS stops before the wall-time
  limit.
- They pass `SIGTERM` to `gmx mdrun`, which then writes a checkpoint.
- They submit themselves again until production is complete.

## Submit GROMACS jobs

Set `engine: gromacs` in the config. `submit` then writes GROMACS scripts.
The `--engine` option overrides the config.

On CPU nodes:

```bash
pixi run -e build polyzymd submit \
    -c config.yaml \
    --preset aa100 \
    --replicates 1-3
```

On GPU nodes:

```bash
pixi run -e build polyzymd submit \
    -c config.yaml \
    --preset blanca-shirts \
    --constraint "A40|A100" \
    --replicates 1-3
```

GROMACS uses the site module or container that the config's `gromacs:` block
names. For CPU and GPU settings, MPI, constraints and recovery, see
{doc}`gromacs_export`.

## Common fixes

### `pixi: command not found`

Make `pixi` available in non-interactive shells. A setting in your login shell
files alone is not enough.

### The job stops with an out-of-memory error

Increase `--memory`. You can also make the system smaller, or test with fewer
polymer chains.

### The config path no longer exists

The job script keeps the config path that you gave at submission. If you move
the config, write the scripts again and submit again.

(hpc-slurm-stop-permanently)=
### Stop a job permanently

See [Stop a chain](#hpc-slurm-stop-a-chain). `scancel` alone is not enough,
because the chain submits itself again within seconds.

## Related pages

- Command options: {doc}`../reference/cli_reference`
- Configuration keys: {doc}`../reference/configuration`
- GROMACS on a cluster: {doc}`gromacs_export`
- Other hardware: {doc}`hardware_platforms`
- First simulation: {doc}`../get_started/quickstart`
