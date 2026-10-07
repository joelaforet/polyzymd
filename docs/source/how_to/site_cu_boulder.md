(cu-boulder-site-notes)=
# Site notes: CU Boulder

Use these notes when you run PolyzyMD on the CU Boulder Research Computing
clusters, Alpine and Blanca. On another cluster, see {doc}`hpc_slurm` and
{doc}`hardware_platforms`.

:::{admonition} Pixi environments
:class: tip

The commands below use `pixi run -e build` for simulations and
`pixi run -e analysis` for analyses. You can also activate an environment
once with `pixi shell -e build`.
:::

## Choose the cluster

CU Boulder runs two SLURM clusters. Alpine is the shared campus cluster.
Blanca has condo nodes that research groups own. Load the module of the
cluster before you submit:

```bash
ml slurm/alpine   # shared campus cluster
ml slurm/blanca   # condo nodes owned by research groups
```

:::{important}
Run `ml slurm/blanca` before `sbatch`, `polyzymd submit` or
`polyzymd analyze --submit` to use Blanca. If you do not, SLURM does not show
the Blanca partitions. `polyzymd status` and `polyzymd cancel` also need the
module of the cluster that runs the jobs.
:::

## Give a partition, an account and a QoS

Both clusters refuse a job without all three. The presets of
`polyzymd submit` (see {doc}`hpc_slurm`) and of `polyzymd analyze --submit`
(see {doc}`hpc_execution`) set them. In a batch script that you write
yourself, use these values:

| Flag | Alpine | Blanca |
|------|--------|--------|
| `--partition` | `amilan` (CPU), `aa100` or `ami100` (GPU) | `blanca-<group>`, for example `blanca-shirts` |
| `--account` | Your allocation | The same as the partition |
| `--qos` | `normal` | The same as the partition, or `preemptable` |

The simulation presets `blanca-shirts` and `blanca-chbe-rdi` use the
`preemptable` QoS. The owner's jobs can then stop a job. The job script
writes a checkpoint and submits a successor.

For example, on Blanca:

```bash
#SBATCH --partition=blanca-shirts
#SBATCH --account=blanca-shirts
#SBATCH --qos=blanca-shirts
```

If you omit `--partition` on Blanca, `sbatch` fails with
*"A partition has not been provided"*. The message refers to Alpine
documentation. Add `--partition=blanca-<group>` and submit again.

To list the accounts, partitions and QoS values that you can use:

```bash
sacctmgr show association user=$USER format=account,partition,qos
```

## Submit simulations

The Alpine presets (`aa100`, `al40` and `testing`) set a lab allocation as the
account. Give your own allocation with `--account`:

```bash
ml slurm/alpine
pixi run -e build polyzymd submit \
    -c config.yaml \
    --preset aa100 \
    --account <account> \
    --replicates 1-5
```

On Blanca, the preset sets the partition, the account and the QoS:

```bash
ml slurm/blanca
pixi run -e build polyzymd submit \
    -c config.yaml \
    --preset blanca-shirts \
    --replicates 1-5
```

For OpenMM, `--pixi-env auto` selects `sim-cuda-12-4` on Blanca.

(excluded-blanca-gpu-nodes)=
## Excluded Blanca GPU nodes

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

## Run GROMACS

### Tested Blanca nodes

A GROMACS acceptance test ran on 14 August 2026 with the `blanca-shirts`
account and QoS. Both short jobs used GROMACS 2024.2 from the site module.

| Mode | Node | Hardware | Result |
|------|------|----------|--------|
| CPU | `bgpu-shirts3` | 2 CPU threads | 20-step test completed |
| GPU | `bgpu-shirts1` | NVIDIA A40, driver 550.90.07 | 20-step GPU nonbonded test completed |

The site module reported CUDA GPU support. These results apply only to the
tested module and nodes. Run the short test again after a module or driver
update.

### Recipes

Each recipe gives the `gromacs:` block of `config.yaml` and the submit
command. Blanca has mixed GPU types, so give `--constraint` for GPU jobs. See
{doc}`run_gromacs` for the reason.

**Alpine A100, 3 GPUs, real MPI.** If your own script runs
`mpirun -np 3 gmx_mpi mdrun ...` with its own GMXRC or PLUMED set-up, put
those parts in `mpi_launcher_flags` and `setup_commands`.

```yaml
gromacs:
  gmx_binary: "gmx_mpi"
  gpu: true
  gpus: 3
  ntmpi: 3
  ntomp: 4
  memory: "64G"
  module_load: "module load intel/2022.1.2 impi/2021.5.0 gromacs/2023.3"
  mpi_launcher_flags: "-np 3"
  env_exports:
    GMX_GPU_DD_COMMS: "true"
    GMX_GPU_PME_PP_COMMS: "true"
    GMX_FORCE_UPDATE_DEFAULT_GPU: "true"
  mdrun_flags: "-pme gpu -nb gpu -bonded gpu -npme 1 -gpu_id 012 -ntomp 4"
```

```bash
pixi run -e build polyzymd submit \
    -c config.yaml \
    --engine gromacs \
    --preset aa100 \
    --replicates 1-3
```

**Alpine MI100 (AMD), Singularity container.** `slurm_ntasks: 16` sets the
SLURM task count apart from the GROMACS rank count (`ntmpi: 3`).

```yaml
gromacs:
  gmx_binary: "gmx"
  command_prefix: "singularity exec --rocm --bind $PWD /projects/shared/gromacs-rocm.sif"
  gpu: true
  gpus: 3
  ntmpi: 3
  ntomp: 3
  slurm_ntasks: 16
  memory: "64G"
  module_load: "module load singularity"
  mdrun_flags: "-pme gpu -nb gpu -bonded gpu"
```

```bash
pixi run -e build polyzymd submit \
    -c config.yaml \
    --engine gromacs \
    --preset aa100 \
    --gpu-type mi100 \
    --replicates 1-3
```

**Alpine Amilan, CPU only.** `--partition` replaces the preset's partition.

```yaml
gromacs:
  gmx_binary: "gmx_mpi"
  ntmpi: 8
  ntomp: 8
  memory: "16G"
  module_load: "module load gcc/11.2.0 openmpi/4.1.1 gromacs/2024.2"
```

```bash
pixi run -e build polyzymd submit \
    -c config.yaml \
    --engine gromacs \
    --preset aa100 \
    --partition amilan \
    --time-limit 24:00:00 \
    --replicates 1-3
```

**Blanca GPU, Intel MPI, A40 or A100.** To pin the job to one node, add
`--nodelist <node>`.

```yaml
gromacs:
  gmx_binary: "gmx_mpi"
  gpu: true
  gpus: 3
  ntmpi: 3
  ntomp: 4
  memory: "64G"
  module_load: "module load intel/2022.1.2 impi/2021.5.0 gromacs/2023.3"
  mpi_launcher_flags: "-np 3 -genv I_MPI_FABRICS shm:tcp"
  env_exports:
    GMX_GPU_DD_COMMS: "true"
    GMX_GPU_PME_PP_COMMS: "true"
    GMX_FORCE_UPDATE_DEFAULT_GPU: "true"
  mdrun_flags: "-pme gpu -nb gpu -bonded gpu -npme 1 -gpu_id 012 -ntomp 4"
```

```bash
ml slurm/blanca
pixi run -e build polyzymd submit \
    -c config.yaml \
    --engine gromacs \
    --preset blanca-shirts \
    --constraint "A40|A100" \
    --email you@university.edu \
    --replicates 1-3
```

**Blanca CPU only, one CPU architecture.** `--constraint "cascadelake"` keeps
the job on nodes with that instruction set. `--time-limit 7-00:00:00` asks
for seven days, a common limit for long runs on condo partitions.

```yaml
gromacs:
  gmx_binary: "gmx_mpi"
  ntmpi: 8
  ntomp: 8
  memory: "16G"
  module_load: "module load gcc/11.2.0 openmpi/4.1.1 gromacs/2024.2"
```

```bash
ml slurm/blanca
pixi run -e build polyzymd submit \
    -c config.yaml \
    --engine gromacs \
    --preset blanca-shirts \
    --constraint "cascadelake" \
    --time-limit 7-00:00:00 \
    --replicates 1-3
```

## Related pages

- Simulation jobs: {doc}`hpc_slurm`
- Analysis jobs: {doc}`hpc_execution`
- GROMACS: {doc}`run_gromacs`
