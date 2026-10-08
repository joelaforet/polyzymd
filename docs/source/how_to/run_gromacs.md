# Run GROMACS simulations

Run a PolyzyMD system with GROMACS on your computer, or as a chain of SLURM
jobs on a cluster. PolyzyMD builds the system with OpenFF, writes the GROMACS
files, and runs minimization, equilibration and production.

:::{admonition} Prerequisites
:class: tip

- A working `config.yaml` validated with `polyzymd validate`
- GROMACS on your computer, or on your cluster (from `module load` or a
  container)
- For a cluster: the SLURM partitions and GPU types that you can use
- PolyzyMD installed with `pixi` (see {doc}`../get_started/installation`)

All commands below assume you prefix with `pixi run -e build` or have
activated the environment with `pixi shell -e build`.
:::

---

## Run locally

Set `engine: gromacs` in `config.yaml`. GROMACS does not read the `openmm:`
block of the template, and it writes `.xtc` trajectories whatever
`output.trajectory_format` says. Then run one replicate:

```bash
pixi run -e build polyzymd run -c config.yaml -r 1
```

The run calls `gmx` from `PATH`. To use another GROMACS binary, give
`--gmx-path /path/to/gmx`. To run GROMACS for a config whose `engine` is
`openmm`, give `--engine gromacs`.

The run calls `gmx mdrun -deffnm <stage> -v` for each stage. It does not read
the run settings in the `gromacs:` block of `config.yaml` (`mdrun_flags`,
`ntmpi`, `ntomp` and the other fields except `analysis_topology`); only
[SLURM jobs](#submit-to-a-slurm-cluster) use them. GROMACS then chooses its own
thread counts and GPU use.

PolyzyMD builds the system, writes the GROMACS files to
`<replicate folder>/gromacs/` and runs each stage in order. It prints the
GROMACS output while it runs. If a stage fails, the run stops and keeps the
files. For the list of files, see [Output files](#gromacs-output-files).

To see the plan without writing a file, add `--dry-run`.

## Export the files only

Use `polyzymd build --format gromacs` when you only want PolyzyMD to build and
parameterize the system:

```bash
pixi run -e build polyzymd build -c config.yaml --format gromacs
```

The handoff is the `.gro` coordinate file, the `.top` topology file and the
component `.itp` parameter files. PolyzyMD also writes MDP files and a run
script from `config.yaml`. You do not need them. Use them as a starting point,
or replace them with your own GROMACS workflow.

---

## Submit to a SLURM cluster

`polyzymd submit` writes one self-resubmitting SLURM script for each
replicate. The script runs minimization, equilibration and production, and
restarts each stage from its checkpoint.

`submit` does not build. Build the GROMACS inputs of each replicate first, in
a compute job (see {doc}`hpc_slurm`):

```bash
pixi run -e build polyzymd build -c config.yaml -r 1-3 --format gromacs
```

`--format gromacs` is the default for a config with `engine: gromacs`. Without
the inputs, `submit` stops with an error and writes no script.

### CPU jobs

Add a `gromacs:` block to `config.yaml` and submit:

```yaml
# config.yaml (add this block alongside your existing config)
gromacs:
  module_load: "module load gcc/11.2.0 openmpi/4.1.1 gromacs/2024.2"
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

### GPU jobs

For GPU runs, set `gpu: true` and use the thread-MPI `gmx` binary (not
`gmx_mpi`). `gpu: true` adds `-nb gpu -pme gpu -bonded gpu` to the `mdrun`
flags that you do not set yourself. It does not add `-update gpu`: GROMACS
updates on the GPU only with `integrator = md`, and the Langevin thermostats
write `integrator = sd`. With `-update gpu` and `sd`, GROMACS 2024 stops with
*"Only the md integrator is supported"*.

```yaml
# config.yaml
gromacs:
  gpu: true
  gpus: 1
  gmx_binary: "gmx"
  ntmpi: 1
  ntomp: 12
  module_load: "module load gcc/11.2.0 openmpi/4.1.1 gromacs/2024.2"
  mdrun_flags: "-pin on"
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
Use `--constraint` to keep the job on a compatible GPU type. GROMACS uses
CUDA kernels compiled ahead of time, so a binary compiled for one GPU type
may not run on another. OpenMM compiles its kernels at launch and does not
need this.
:::

For every field of the `gromacs:` block and the `gmx mdrun` flags, see
{ref}`config-gromacs`. For thread-MPI, real MPI and OpenMP threads, see
{doc}`../explanation/gromacs_parallelism`.

## Use GROMACS on different clusters

GROMACS does not use the OpenMM platform router. PolyzyMD does not inspect the
CUDA driver for a GROMACS job. The site GROMACS module or container supplies
the CPU, CUDA, or ROCm runtime.

The shared `--pixi-env auto` option maps to `build` for GROMACS. The generated
script activates `build` for the PolyzyMD commands. It then runs the configured
`module_load`, `env_exports`, and `setup_commands` before it starts GROMACS.
`submit` starts the job with `sbatch --export=NONE`, so the job does not
inherit the environment of the submitting shell. `module_load` runs only in
the job, never on the login node.
If your `module_load` loads the scheduler module, load that module in your
shell before `submit` instead; without `sbatch` on `PATH`,
`submit` stops with an error.

Use these settings to move a job to another SLURM cluster:

| Requirement | PolyzyMD setting |
|-------------|------------------|
| Partition | `--partition` or `--preset` |
| Account | `--account` |
| QoS | `--qos` |
| GPU type | `--gpu-type` |
| Node feature | `--constraint` |
| Site GROMACS installation | `gromacs.module_load` |
| Container runtime | `gromacs.command_prefix` |
| Site environment variables | `gromacs.env_exports` |
| CPU allocation | `gromacs.ntmpi`, `gromacs.ntomp`, and `gromacs.memory` |
| GPU allocation | `gromacs.gpu: true` and `gromacs.gpus` |

An existing preset can supply the SLURM directive style. Override its account,
partition, QoS, and GPU type for the target site. Use `aa100` for sites that
accept `--gres=gpu`. Use `bridges2` for sites that accept `--gpus=<type>:<n>`.
Always inspect a generated script before the first submission on a new site.

```bash
pixi run -e build polyzymd submit \
    -c config.yaml \
    --engine gromacs \
    --preset aa100 \
    --partition <site-partition> \
    --account <site-account> \
    --qos <site-qos> \
    --generate-only
```

The CPU and GPU scripts use the same restart process:

1. GROMACS writes a checkpoint during `mdrun`, every
   `simulation_phases.production.checkpoint_interval` seconds (`mdrun -cpt`).
2. The script passes `-cpi` and `-append` when it restarts production.
3. The script forwards `SIGTERM` to `mdrun` so GROMACS can write a checkpoint.
4. The script submits one successor when work remains.
5. The successor uses the same resource request, module, and GROMACS flags.

For recipes for the CU Boulder clusters, see {doc}`site_cu_boulder`.

### Bridges-2 CPU job

Use a CPU partition that is valid for your Bridges-2 allocation. The
`gromacs.gpu` field must be false:

```yaml
gromacs:
  gpu: false
  gmx_binary: gmx_mpi
  ntmpi: 4
  ntomp: 4
  memory: 16G
  module_load: "module load <site-gromacs-module>"
```

```bash
pixi run -e build polyzymd submit \
    -c config.yaml \
    --engine gromacs \
    --preset bridges2 \
    --partition <bridges2-cpu-partition> \
    --account <allocation> \
    --replicates 1-3
```

PolyzyMD omits the GPU directive when `gromacs.gpu` is false.

### Bridges-2 GPU job

Use the Bridges-2 GPU type that your allocation permits:

```yaml
gromacs:
  gpu: true
  gpus: 1
  gmx_binary: gmx
  ntmpi: 1
  ntomp: 8
  module_load: "module load <site-gromacs-module>"
  mdrun_flags: "-nb gpu -pin on"
```

```bash
pixi run -e build polyzymd submit \
    -c config.yaml \
    --engine gromacs \
    --preset bridges2 \
    --account <allocation> \
    --gpu-type v100-32 \
    --replicates 1-3
```

The site GROMACS build controls GPU compatibility. Test `gmx --version` and a
short `gmx mdrun` in an allocation before a production campaign. Do not use the
OpenMM `sim-cuda-12-*` environments to select a GROMACS runtime.

### Translate an existing SBATCH script

Each directive of a SLURM batch script has a place in PolyzyMD:

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

## GPU constraints and preemption

### Why GPU constraints matter

GROMACS uses CUDA kernels compiled ahead of time. OpenMM compiles its kernels
at launch. A GROMACS binary compiled for one GPU architecture can crash on
another. If your cluster has mixed GPU types, always use `--constraint`:

```bash
--constraint "A40"              # single GPU type
--constraint "A40|A100"         # either type (OR)
--constraint "avx2&rh8"         # feature AND (CPU + OS flags)
```

This maps directly to `#SBATCH --constraint` in the generated script.

(gromacs-preemption-resilience)=
### Preemption

GROMACS SLURM scripts trap `SIGTERM`, the signal that SLURM sends before it
preempts a job. When the trap fires:

1. The script forwards `SIGTERM` to `gmx mdrun`.
2. GROMACS writes a `.cpt` checkpoint file.
3. The script waits for GROMACS to exit.
4. The script submits itself again with `sbatch`.

With `--constraint`, the next job also lands on a compatible GPU. This
matters most with a preemptable QoS.

The script sets `-maxh` to 90 % of the wall time, so GROMACS stops before
the SLURM wall-time limit. Production is complete when the last checkpoint
in `prod.log` is at `nsteps`; until then the script submits a successor.
Progress counts only the steps up to the last checkpoint, because a restart
runs the later steps again.

```{note}
**Stop a GROMACS chain.**
Use `polyzymd cancel`, as for OpenMM. `scancel <job_id>` alone sends
`SIGTERM`, which the job script reads as a preemption: it submits a
successor. `polyzymd cancel` first writes a `STOP` file in the run folder.
The job script checks for it at start and before each successor, so the
chain ends. See {ref}`Stop a chain <hpc-slurm-stop-a-chain>` in the SLURM
guide.
```

---

## Recover a stopped chain

If the automatic resubmission fails (for example, `sbatch` was not available,
or the job was killed without a grace period), use `polyzymd recover`.

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
    --constraint "A40" \
    --email you@university.edu
```

The job script is written to `recovery_scripts/recover_rep<N>.sh`. The
chain's own `daisy_chain_scripts/run_rep<N>.sh` is left as it was.

### How checkpoint resume works

| Stage | Checkpoint | Resume behavior |
|-------|-----------|-----------------|
| Energy minimization | `em.cpt` | Resumes from last EM step |
| Equilibration stage N | `eq_XX.cpt` | Resumes from last equilibration step |
| Production | `state.cpt` | Resumes from last production checkpoint |

These are the checkpoints of the SLURM job script. A local `polyzymd run`
writes the production checkpoint to `prod.cpt`.

Completed stages (those with a `.gro` output file) are skipped on
resubmission. Partially completed stages resume from their checkpoint.

### Dry-run recovery preview

```bash
pixi run -e build polyzymd recover \
    -c config.yaml -r 1 \
    --engine gromacs \
    --submit \
    --dry-run
```

This prints the script path and the SLURM settings, and writes nothing.

---

(gromacs-output-files)=
## Output files

A GROMACS run writes its files to `<replicate folder>/gromacs/`. `<prefix>` is
the enzyme name, followed by `_<polymers.type_prefix>` when the system has
polymers. `<stage>` is the name of an equilibration stage. A local run of the
quickstart writes these files:

```text
gromacs/
├── <prefix>.gro              # Initial coordinates
├── <prefix>.top              # Topology (includes all molecule types)
├── <prefix>_*.itp            # Atom types and molecule parameters
├── <prefix>_pointenergy.mdp  # Single-point energy check
│
├── em.mdp                    # Energy minimization parameters
├── eq_01_<stage>.mdp         # One file per equilibration stage
├── prod.mdp                  # Production parameters
├── mdout.mdp                 # Parameters as grompp read them
│
├── run_<prefix>_gromacs.sh   # Generated shell script
│
├── em.tpr, em.gro, em.edr    # Energy minimization outputs
├── eq_01.tpr, eq_01.gro, eq_01.cpt, eq_01.xtc
│                             # Outputs of equilibration stage 1
│
├── prod.tpr                  # Production run input
├── prod.xtc                  # Production trajectory
├── prod.edr                  # Production energies
├── prod.gro                  # Final coordinates
├── prod.cpt                  # Checkpoint for restart (state.cpt in a SLURM job)
│
├── prod_nojump.xtc           # Trajectory with PBC jumps removed
├── prod_centered.xtc         # Centered trajectory for visualization
└── progress.json             # Stage and segment records: times, seeds
```

Each `mdrun` stage also writes a `.log` file, and a `.trr` file when the stage
writes full-precision coordinates. `solvated_system.pdb` and
`build_manifest.json` (PACKMOL seeds, box and the SHA-256 of each input file)
are in the replicate folder, beside `gromacs/`.

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

**Cause**: GPU offload flags that energy minimization does not accept.

**Fix**: PolyzyMD removes `-pme gpu`, `-bonded gpu` and `-update gpu` from the
minimization `mdrun` command. If you see this error, check that you use
PolyzyMD 1.3.0 or later.

### "Only the md integrator is supported"

**Cause**: `-update gpu` in `mdrun_flags` with a Langevin thermostat, which
GROMACS runs as `integrator = sd`.

**Fix**: Remove `-update gpu` from `mdrun_flags`, or give `-update cpu`.

### "grompp stops with a warning"

**Cause**: PolyzyMD passes no `-maxwarn`, so every `grompp` warning stops the
run. A neutral system built by PolyzyMD gives none. A common one is "You are
using Ewald electrostatics in a system with net charge": the system is not
neutral.

**Fix**: Read the warning. For a net charge, set `solvent.ions.neutralize:
true`. To accept a warning on purpose, set `grompp_flags: "-maxwarn 1"`.
Only the SLURM job scripts from `polyzymd submit` pass `grompp_flags`; a
local `polyzymd run` calls `grompp` without them. See
{doc}`../reference/gromacs_openmm`.

### "Fatal error: Number of atoms does not match"

**Cause**: Topology/coordinate mismatch from an interrupted build.

**Fix**: Delete the GROMACS files and build again:

```bash
rm -rf <replicate folder>/gromacs/
pixi run -e build polyzymd build -c config.yaml --format gromacs
```

### Job dies with OOM

**Fix**: Increase `gromacs.memory` in config or use `--memory` CLI override:

```bash
polyzymd submit -c config.yaml --engine gromacs --memory 64G ...
```

### Trajectory has broken molecules

**Fix**: Use the post-processed trajectories:
- `prod_nojump.xtc`: molecules do not jump across PBC boundaries
- `prod_centered.xtc`: system centered for visualization

---

## See Also

- {doc}`../reference/cli_reference` — CLI options for `run`, `submit` and `recover`
- {doc}`../reference/configuration` — Full configuration reference including `gromacs:` block
- {doc}`../explanation/gromacs_parallelism` — Thread-MPI, real MPI and OpenMP threads
- {doc}`hpc_slurm` — General SLURM submission workflow (OpenMM and GROMACS)
- {doc}`site_cu_boulder` — CU Boulder Alpine and Blanca recipes
- {doc}`../get_started/quickstart` — Getting started guide
