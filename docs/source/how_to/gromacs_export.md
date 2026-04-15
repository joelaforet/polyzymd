# Run GROMACS Simulations on HPC Clusters

This guide shows how to configure, submit, and manage GROMACS simulations on
HPC clusters using PolyzyMD's self-resubmitting SLURM workflow.

:::{admonition} Prerequisites
:class: tip

- A working `config.yaml` validated with `polyzymd validate`
- GROMACS installed on your cluster (via `module load` or container)
- Familiarity with your cluster's SLURM partitions and GPU types
- PolyzyMD installed with `pixi` (see {doc}`../tutorials/installation`)

All commands below assume you prefix with `pixi run -e build` or have
activated the environment with `pixi shell -e build`.
:::

---

## Quick Start: CPU

Add a minimal `gromacs:` block to your `config.yaml` and submit:

```yaml
# config.yaml (add this block alongside your existing config)
gromacs:
  module_load: "module load gcc/11.2.0 gromacs/2024.2"
  ntmpi: 1
  ntomp: 8
```

```bash
pixi run -e build polyzymd submit \
    -c config.yaml \
    --engine gromacs \
    --preset aa100 \
    --replicates 1-3
```

This generates self-resubmitting SLURM scripts that run EM, equilibration,
and production with checkpoint-based restart.

---

## Quick Start: GPU

For GPU-accelerated runs, set `gpu: true` and use the thread-MPI `gmx`
binary (not `gmx_mpi`):

```yaml
# config.yaml
gromacs:
  gpu: true
  gpus: 1
  gmx_binary: "gmx"
  ntmpi: 1
  ntomp: 12
  module_load: "module load gcc/11.2.0 gromacs/2024.2"
  mdrun_flags: "-nb gpu -pme gpu -bonded gpu -update gpu -pin on"
```

```bash
pixi run -e build polyzymd submit \
    -c config.yaml \
    --engine gromacs \
    --preset blanca-shirts \
    --constraint "A40" \
    --replicates 1-3
```

:::{important}
Use `--constraint` to pin your job to a compatible GPU architecture. GROMACS
uses ahead-of-time compiled CUDA kernels, so a binary compiled for one GPU
type may not run on another. OpenMM does not need this because it JIT-compiles
kernels at launch.
:::

---

## GROMACS Flag Glossary

These flags are passed to `gmx mdrun` via the `mdrun_flags`,
`mdrun_flags_equilibration`, and `mdrun_flags_production` config fields.

### Parallelism flags

| Flag | Description |
|------|-------------|
| `-ntmpi N` | Thread-MPI ranks. Used by thread-MPI builds (`gmx`). Set via `gromacs.ntmpi` in config. |
| `-ntomp N` | OpenMP threads per rank. Set via `gromacs.ntomp` in config. |
| `-npme N` | Dedicated PME ranks. Useful with multi-GPU runs (e.g., `-npme 1` with 3 GPUs). |

### GPU offload flags

| Flag | Description | Safe for EM? |
|------|-------------|:---:|
| `-nb gpu` | Offload nonbonded forces to GPU | Yes |
| `-pme gpu` | Offload PME electrostatics to GPU | **No** |
| `-bonded gpu` | Offload bonded forces to GPU | **No** |
| `-update gpu` | Offload integration/constraints to GPU | **No** |

### Performance flags

| Flag | Description |
|------|-------------|
| `-pin on` | Pin threads to CPU cores (recommended for performance) |
| `-pinstride N` | Stride between pinned threads |
| `-dlb yes\|auto` | Dynamic load balancing |
| `-gpu_id NNN` | Explicit GPU device IDs (e.g., `012` for 3 GPUs) |

### Energy minimization restrictions

:::{warning}
Only `-nb gpu` is safe during energy minimization. GROMACS will fail
with "Non-dynamical integrator" if `-pme gpu`, `-bonded gpu`, or
`-update gpu` are used during EM.

**PolyzyMD handles this automatically** — it strips unsafe GPU flags from
the `mdrun` command during EM stages. You do not need separate EM-specific
flag configuration.
:::

---

## Understanding MPI vs OpenMP Parallelism in GROMACS

GROMACS supports two parallelism models that are easy to confuse:

| Concept | Thread-MPI | Real MPI |
|---------|-----------|----------|
| **Binary** | `gmx` | `gmx_mpi` |
| **Launch** | Direct (`gmx mdrun -ntmpi N`) | Via launcher (`mpirun -np N gmx_mpi mdrun`) |
| **GPU support** | Yes (CUDA thread-MPI) | Varies by build (often GPU-disabled) |
| **Multi-node** | No (single node only) | Yes |
| **Config field** | `gromacs.ntmpi` | `gromacs.ntmpi` + `gromacs.mpi_launcher_flags` |

### Rule of thumb

```
ntmpi × ntomp = total CPU cores allocated
```

**GPU runs** typically use thread-MPI with 1 MPI rank and many OpenMP threads:

```yaml
gromacs:
  gmx_binary: "gmx"       # thread-MPI build
  ntmpi: 1                 # 1 rank = 1 GPU
  ntomp: 12                # 12 OpenMP threads
```

**CPU runs** can use either model. For multi-node runs, use real MPI:

```yaml
gromacs:
  gmx_binary: "gmx_mpi"   # real MPI build
  ntmpi: 8                 # 8 MPI ranks
  ntomp: 1                 # 1 OpenMP thread per rank
```

:::{tip}
On most clusters, the `gmx` binary (thread-MPI) is the best choice for
single-node GPU runs. Use `gmx_mpi` only when you need multi-node scaling.
:::

---

## Complete `gromacs:` Config Reference

All fields in the `gromacs:` block of your `config.yaml`:

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `gmx_binary` | `str \| null` | `null` | GROMACS binary path or name. Resolved via config > `$GMX_BIN` > PATH if null. |
| `mdrun_flags` | `str` | `""` | Extra flags for `gmx mdrun` (applied to all stages). |
| `mdrun_flags_equilibration` | `str \| null` | `null` | Override `mdrun_flags` for equilibration only. Falls back to `mdrun_flags` if null. |
| `mdrun_flags_production` | `str \| null` | `null` | Override `mdrun_flags` for production only. Falls back to `mdrun_flags` if null. |
| `grompp_flags` | `str` | `"-maxwarn 1"` | Extra flags for `gmx grompp`. |
| `command_prefix` | `str \| null` | `null` | Prefix prepended to all GROMACS commands (e.g., Singularity wrapper). |
| `mpi_launcher_flags` | `str` | `""` | Extra flags for the MPI launcher (`mpirun`). Only used with real-MPI binaries. |
| `module_load` | `str \| null` | `null` | Module load command inserted into SLURM scripts verbatim. |
| `env_exports` | `dict[str, str]` | `{}` | Environment variables exported before GROMACS commands. |
| `setup_commands` | `list[str]` | `[]` | Shell commands run after `module_load` and before GROMACS. |
| `ntmpi` | `int` | `1` | MPI ranks (`-ntmpi`). Also sets SLURM `--ntasks` (unless `slurm_ntasks` overrides). |
| `slurm_ntasks` | `int \| null` | `null` | Override SLURM `--ntasks` independently of `-ntmpi`. For multi-node MPI+GPU workflows. |
| `ntomp` | `int` | `8` | OpenMP threads per rank (`-ntomp`). Sets SLURM `--cpus-per-task`. |
| `gpu` | `bool` | `false` | Request GPU via SLURM. When false, `--gres=gpu` is omitted entirely. |
| `gpus` | `int` | `1` | Number of GPUs to request when `gpu` is true. |
| `memory` | `str` | `"16G"` | SLURM `--mem` allocation for GROMACS jobs. |

### Stage-specific mdrun flags

Use `mdrun_flags_equilibration` and `mdrun_flags_production` to apply
different flags during different simulation phases. When set to `null`
(default), they fall back to `mdrun_flags`.

```yaml
gromacs:
  # Applied to all stages by default
  mdrun_flags: "-pin on"

  # Override for equilibration only (more conservative)
  mdrun_flags_equilibration: "-pin on -dlb yes"

  # Override for production only (full GPU offload)
  mdrun_flags_production: "-pin on -dlb auto -nb gpu -pme gpu -bonded gpu -update gpu"
```

### command_prefix vs mpi_launcher_flags

These two fields serve different purposes. Use **one** approach, not both:

**`command_prefix`** — wraps all GROMACS commands with a prefix. Use for
containers or site-specific launchers:

```yaml
gromacs:
  command_prefix: "singularity exec --rocm --bind $PWD /path/to/gromacs.sif"
```

**`mpi_launcher_flags`** — passes extra flags to the `mpirun` launcher that
PolyzyMD generates for real-MPI binaries:

```yaml
gromacs:
  gmx_binary: "gmx_mpi"
  mpi_launcher_flags: "-genv I_MPI_FABRICS shm:tcp"
  # Result: mpirun -genv I_MPI_FABRICS shm:tcp gmx_mpi mdrun ...
```

:::{note}
When `command_prefix` is set with a real-MPI binary (`gmx_mpi`), PolyzyMD
skips automatic `mpirun` wrapping to avoid double-launching. Any
`mpi_launcher_flags` are ignored (a warning is emitted).
:::

### env_exports

Environment variables exported in the SLURM script before GROMACS commands.
Keys must be valid shell variable names (`[A-Za-z_][A-Za-z0-9_]*`).

```yaml
gromacs:
  env_exports:
    GMX_GPU_DD_COMMS: "true"
    GMX_GPU_PME_PP_COMMS: "true"
    GMX_FORCE_UPDATE_DEFAULT_GPU: "true"
    OMP_PLACES: "cores"
    OMP_PROC_BIND: "close"
```

### setup_commands

Shell commands run in order after `module_load` and `env_exports`, before
any GROMACS commands:

```yaml
gromacs:
  setup_commands:
    - "ulimit -s unlimited"
    - "source /opt/gromacs-2024/bin/GMXRC"
```

---

## SBATCH to PolyzyMD YAML Mapping

If you are translating an existing SLURM batch script to PolyzyMD, use this
table to find where each directive goes:

| SBATCH Directive | PolyzyMD Equivalent | Location |
|-----------------|---------------------|----------|
| `--partition` | `--partition` CLI or `--preset` | CLI override |
| `--qos` | `--qos` CLI | CLI override |
| `--time` | `--time-limit` CLI | CLI override |
| `--gres=gpu:TYPE:N` | `gromacs.gpu` + `gromacs.gpus` + `--gpu-type` CLI | Config + CLI |
| `--constraint` | `--constraint` CLI | CLI override |
| `--nodelist` | `--nodelist` CLI | CLI override |
| `--ntasks` | `gromacs.slurm_ntasks` or `gromacs.ntmpi` | Config |
| `--cpus-per-task` | `gromacs.ntomp` | Config |
| `--mem` | `gromacs.memory` or `--memory` CLI | Config + CLI |
| `--mail-user` | `--email` CLI | CLI override |
| `--account` | `--account` CLI | CLI override |
| `module load ...` | `gromacs.module_load` | Config |
| `export VAR=value` | `gromacs.env_exports` | Config |
| Setup commands | `gromacs.setup_commands` | Config |
| `mpirun` flags | `gromacs.mpi_launcher_flags` | Config |
| `singularity exec` | `gromacs.command_prefix` | Config |

---

## Cluster Recipes

Copy-pasteable configurations for common cluster setups. Each recipe shows
the `gromacs:` YAML block and the CLI submit command.

### Alpine AA100 (NVIDIA A100, 3 GPUs, real MPI)

Multi-GPU run using real MPI (`gmx_mpi`) with Intel compiler stack.

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

:::{note}
If your existing workflow uses `mpirun -np 3 gmx_mpi mdrun ...` with custom
GMXRC sourcing or Plumed integration, translate those into the
`mpi_launcher_flags` and `setup_commands` config fields.
:::

### Alpine AMI100 (AMD MI100, Singularity container)

AMD GPU run using a ROCm container with Singularity.

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

:::{note}
`slurm_ntasks: 16` decouples the SLURM task count from the GROMACS
thread-MPI rank count (`ntmpi: 3`). This is needed when the scheduler
requires more tasks than GROMACS MPI ranks (e.g., for container overhead
or multi-GPU resource allocation).
:::

### Alpine Amilan (CPU-only)

CPU-only run with direct MPI.

```yaml
gromacs:
  gmx_binary: "gmx_mpi"
  ntmpi: 8
  ntomp: 1
  memory: "16G"
  module_load: "module load gcc/11.2.0 openmpi/4.1.1 gromacs/2024.2"
  mdrun_flags: "-ntomp 8"
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

:::{note}
CPU-only runs do not need `gpu: true` or `--constraint`. The `--partition`
CLI flag overrides the preset's default partition.
:::

### Blanca GPU (preemptable, Intel MPI, A40/A100)

Multi-GPU run on Blanca's preemptable QoS with Intel MPI fabric settings.

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
ml slurm/blanca  # required to see Blanca partitions

pixi run -e build polyzymd submit \
    -c config.yaml \
    --engine gromacs \
    --preset blanca-shirts \
    --constraint "A40|A100" \
    --replicates 1-3
```

:::{important}
Blanca uses `qos=preemptable`, meaning jobs can be killed by the node owner.
The generated SLURM scripts handle this gracefully — see
[Preemption resilience](#gromacs-preemption-resilience) below.

You must run `module load slurm/blanca` before `sbatch` to see Blanca
partitions.
:::

### Dedicated GPU Node (nodelist pinning, A40)

Pin a job to a specific GPU node using `--nodelist`. Useful for lab-owned or
reserved nodes where you know the hardware.

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
    --constraint "A40" \
    --nodelist "gpu-node-001" \
    --replicates 1-3
```

:::{note}
`--nodelist` pins the job to a specific node. This is useful when you have
dedicated or reserved nodes with known GPU hardware. Replace `gpu-node-001`
with your actual node hostname.
:::

### CPU-Only with Architecture Constraint (cascadelake)

CPU-only run on Blanca with CPU architecture pinning.

```yaml
gromacs:
  gmx_binary: "gmx_mpi"
  ntmpi: 8
  ntomp: 1
  memory: "16G"
  module_load: "module load gcc/11.2.0 openmpi/4.1.1 gromacs/2024.2"
  mdrun_flags: "-ntomp 8"
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

:::{note}
CPU constraints like `cascadelake` ensure the job runs on nodes with
a compatible instruction set. A 7-day wall time (`7-00:00:00`) is common
for long production runs on condo partitions.
:::

---

## GPU Constraints and Preemption

### Why GPU constraints matter

GROMACS uses ahead-of-time compiled CUDA kernels (unlike OpenMM, which
JIT-compiles at launch). A binary compiled for one GPU architecture may
crash on a different one. If your cluster has mixed GPU types, always use
`--constraint`:

```bash
--constraint "A40"              # single GPU type
--constraint "A40|A100"         # either type (OR)
--constraint "avx2&rh8"         # feature AND (CPU + OS flags)
```

This maps directly to `#SBATCH --constraint` in the generated script.

(gromacs-preemption-resilience)=
### Preemption resilience

GROMACS SLURM scripts trap `SIGTERM` (the signal SLURM sends before
preempting a job). When the trap fires:

1. The script forwards SIGTERM to `gmx mdrun`
2. GROMACS flushes a `.cpt` checkpoint file
3. The script waits for GROMACS to exit
4. The script resubmits itself via `sbatch`

Combined with `--constraint`, this ensures resumed jobs land on compatible
GPU hardware. This is especially important on clusters with `qos=preemptable`
(e.g., Blanca).

The `-maxh` flag is automatically set so GROMACS exits cleanly before the
SLURM wall-time limit.

---

## Recovering Preempted Jobs

If automatic resubmission fails (e.g., `sbatch` was unavailable, or the job
was killed without a grace period), use `polyzymd recover`:

### Check recovery status

```bash
pixi run -e build polyzymd recover \
    -c config.yaml -r 1 --engine gromacs
```

This shows per-stage progress without submitting anything.

### Submit a recovery job

```bash
pixi run -e build polyzymd recover \
    -c config.yaml -r 1 \
    --engine gromacs \
    --submit \
    --preset blanca-shirts \
    --constraint "A40"
```

### How checkpoint resume works

| Stage | Checkpoint | Resume behavior |
|-------|-----------|-----------------|
| Energy minimization | `em.cpt` | Resumes from last EM step |
| Equilibration stage N | `eq_XX.cpt` | Resumes from last equilibration step |
| Production | `state.cpt` | Resumes from last production checkpoint |

Completed stages (those with a `.gro` output file) are automatically skipped
on resubmission. Partially completed stages resume from their checkpoint.

### Dry-run recovery preview

```bash
pixi run -e build polyzymd recover \
    -c config.yaml -r 1 \
    --engine gromacs \
    --submit \
    --dry-run
```

---

## GROMACS Output Files

When a GROMACS job completes, files are located in
`{projects_dir}/replicate_{N}/gromacs/`:

```text
gromacs/
├── {system}.gro              # Initial coordinates
├── {system}.top              # Topology (includes all molecule types)
├── *.itp                     # Molecule parameters (one per component)
│
├── em.mdp                    # Energy minimization parameters
├── eq_01_heating.mdp         # Heating stage
├── eq_02_free_equilibration.mdp
├── prod.mdp                  # Production MD parameters
│
├── run_{system}_gromacs.sh   # Generated shell script
│
├── em.tpr, em.gro, em.edr   # Energy minimization outputs
├── eq_01.*, eq_02.*          # Equilibration outputs
│
├── prod.tpr                  # Production run input
├── prod.xtc                  # Production trajectory
├── prod.edr                  # Production energies
├── prod.gro                  # Final coordinates
├── state.cpt                 # Checkpoint for restart
│
├── prod_nojump.xtc           # Trajectory with PBC jumps removed
└── prod_centered.xtc         # Centered trajectory for visualization
```

Position restraints are appended as `#ifdef POSRES_*` blocks inside the
molecule `.itp` files. MDP files use `-DPOSRES_PROTEIN`, `-DPOSRES_POLYMER`,
etc. to activate them during equilibration stages.

---

## Troubleshooting

### "GROMACS executable not found"

**Cause**: `gmx` command not in PATH after module loading.

**Fix**: Check your `gromacs.module_load` field. List prerequisites
(compiler, MPI) before the GROMACS module:

```yaml
gromacs:
  module_load: "module load gcc/11.2.0 openmpi/4.1.1 gromacs/2024.2"
```

### "Non-dynamical integrator" error during EM

**Cause**: GPU offload flags incompatible with energy minimization.

**Fix**: This should not happen — PolyzyMD automatically strips unsafe GPU
flags during EM. If it does, check that you are using PolyzyMD v1.3.0 or
later.

### "grompp warnings about charge groups"

**Cause**: OpenFF generates systems without charge groups.

**Fix**: This is expected and safe. The `grompp_flags: "-maxwarn 1"` default
suppresses it.

### "Fatal error: Number of atoms does not match"

**Cause**: Topology/coordinate mismatch from an interrupted build.

**Fix**: Rebuild from scratch:

```bash
rm -rf replicate_*/gromacs/
pixi run -e build polyzymd build -c config.yaml --format gromacs
```

### Job dies with OOM

**Fix**: Increase `gromacs.memory` in config or use `--memory` CLI override:

```bash
polyzymd submit -c config.yaml --engine gromacs --memory 64G ...
```

### Trajectory has broken molecules

**Fix**: Use the post-processed trajectories:
- `prod_nojump.xtc` — molecules don't jump across PBC boundaries
- `prod_centered.xtc` — system centered for visualization

---

## See Also

- {doc}`../reference/cli_reference` — CLI options for `submit` and `recover`
- {doc}`../reference/configuration` — Full configuration reference including `gromacs:` block
- {doc}`hpc_slurm` — General SLURM submission workflow (OpenMM and GROMACS)
- {doc}`../tutorials/quickstart` — Getting started guide
