# Run Conjugated-System Workflows on CUDA HPC Environments

Use this guide when you need to prepare enzyme-polymer conjugates on one
machine and run OpenMM production simulations on an HPC cluster whose GPU nodes
provide CUDA 12.4.

The short version is:

1. Use a **full-stack build environment** for `polyzymd build` and conjugation
   preparation.
2. Use `sim-cuda-12-4` only for the prepared simulation runtime on CUDA 12.4
   nodes.
3. Do not try to combine the modern OpenFF build stack with the legacy CUDA
   12.4 OpenMM/NumPy stack in one pixi solve.

:::{important}
`polyzymd build` for conjugated systems is now a full-stack workflow. It creates
OpenFF Interchange objects and runs OpenMM restrained minimization plus vacuum
smoke MD before solvation. A lean simulation-only environment is therefore not
enough for full conjugation builds.
:::

## Choose the right environment

| Goal | Recommended pixi environment | Why |
|------|------------------------------|-----|
| Local CPU build, parameterization, validation, docs, and development | `build` | Full OpenFF/OpenMM/analysis stack without CUDA pinning. If a future or site-specific `build-cpu` environment exists, it should serve this role. |
| Local full-stack CUDA build on a workstation with a CUDA 12.6-compatible driver | `build-cuda-12-6` | Includes modern OpenFF tooling, OpenFF Interchange, NumPy 2, analysis tools, and OpenMM with CUDA 12.6. |
| CUDA 12.4 HPC OpenMM simulation runtime | `sim-cuda-12-4` | Legacy runtime stack: OpenMM 8.1.x, NumPy 1.x, and SciPy `<1.14`, matched to CUDA 12.4 nodes. |
| CUDA 12.6 simulation-only handoff runtime | `sim-cuda-12-6` | Initializes CUDA 12.6 OpenMM, but does not include OpenFF build tooling. Use only after preparation is complete. |

Do not use `sim-cuda-12-4` or `sim-cuda-12-6` for full conjugation builds.
They are intentionally lean simulation runtimes. In particular, `sim-cuda-12-6`
can initialize CUDA, but it does not contain the OpenFF stack needed by
`polyzymd build`.

## Why CUDA 12.4 needs a separate runtime stack

CUDA 12.4 HPC nodes currently require the legacy PolyzyMD simulation stack:

- OpenMM 8.1.x
- NumPy 1.x
- SciPy `<1.14`

The current OpenFF build and parameterization stack resolves against newer
dependencies, including NumPy 2 and OpenMM 8.4-era packages. Those constraints
are not compatible with the CUDA 12.4 simulation runtime in a single modern pixi
environment.

Keep these concerns separate:

- **Build stack:** chemistry perception, OpenFF Interchange creation,
  conjugation, restrained minimization, vacuum smoke MD, solvation, export.
- **Simulation runtime:** load an already prepared simulation artifact, create
  or load the OpenMM objects needed for production, and run on the cluster GPU.

## Recommended two-stage HPC workflow

### Stage 1: build and prepare in a full-stack environment

Run conjugation build on a workstation, CPU node, or local GPU environment that
has the full OpenFF/OpenMM stack.

CPU-oriented preparation:

```bash
pixi run -e build polyzymd validate -c config.yaml
pixi run -e build polyzymd build -c config.yaml
```

Local CUDA 12.6 preparation, when you want OpenMM CUDA available during the
build-time minimization and smoke-MD checks:

```bash
pixi run -e build-cuda-12-6 polyzymd validate -c config.yaml
pixi run -e build-cuda-12-6 polyzymd build -c config.yaml
```

:::{tip}
Use `build-cuda-12-6` rather than the default `build` environment on local
machines where the default solve selects a CUDA version newer than the driver can
run. For example, a driver that supports CUDA 13.1 can still fail if the
environment resolves CUDA 13.2 kernels.
:::

### Stage 2: run prepared inputs on CUDA 12.4 HPC

Use `sim-cuda-12-4` only after Stage 1 has already produced runtime inputs that
the simulation job can load directly. In the current CLI, `--skip-build` tells
generated OpenMM jobs to skip `SystemBuilder` and load pre-built files such as
`solvated_system.pdb` and `system.xml` from the replicate working directory.

If your project has those pre-built runtime inputs in the expected locations,
generate or submit the CUDA 12.4 job with `--skip-build`:

```bash
pixi run -e sim-cuda-12-4 polyzymd submit \
    -c config.yaml \
    --preset aa100 \
    --replicates 1-5 \
    --skip-build
```

:::{warning}
Do not use `sim-cuda-12-4 polyzymd submit` as a full conjugation build command.
With no pre-built runtime inputs, the generated job would need the build stack,
including OpenFF Interchange and build-time OpenMM checks, which this lean
runtime environment intentionally does not provide.
:::

The intended long-term production pattern is to transfer a final prepared
Interchange or equivalent simulation artifact to the cluster and reload it in a
simulation-only runtime. That exact user-facing export/reload handoff may still
need implementation or site-specific scripting. Do not assume a handoff command
exists unless it appears in `polyzymd --help` or {doc}`../reference/cli_reference`.

## Verify each environment

Run these checks in the environment you plan to use. Always run CUDA checks on a
GPU node or an interactive GPU allocation, not on a login node.

### Full-stack CPU build environment

```bash
pixi run -e build python -c "import numpy, openmm, openff.toolkit, openff.interchange; print('full stack ok', numpy.__version__, openmm.version.version)"
pixi run -e build polyzymd validate -c config.yaml
```

### Full-stack CUDA 12.6 build environment

```bash
pixi run -e build-cuda-12-6 python -m openmm.testInstallation
pixi run -e build-cuda-12-6 python -c "import openff.toolkit, openff.interchange, openmm; print('build CUDA stack ok', openmm.version.version)"
```

### CUDA 12.4 HPC simulation runtime

```bash
pixi run -e sim-cuda-12-4 python -m openmm.testInstallation
pixi run -e sim-cuda-12-4 python -c "import numpy, scipy, openmm; print('sim stack ok', numpy.__version__, scipy.__version__, openmm.version.version)"
```

If OpenMM cannot initialize CUDA, fix the environment or node allocation before
submitting production jobs.

## Select the OpenMM platform explicitly

Production configs and CLI workflows should select the intended OpenMM platform
explicitly and fail early when CUDA is not available or not compatible with the
driver. Use `python -m openmm.testInstallation` as the first diagnostic, then
run a short PolyzyMD validation or smoke job before submitting many replicates.

For SLURM workflows, start with a short generated or submitted test:

```bash
pixi run -e sim-cuda-12-4 polyzymd submit \
    -c config.yaml \
    --preset testing \
    --time-limit 0:05:00 \
    --replicates 1 \
    --skip-build \
    --generate-only
```

Inspect the generated script, then submit a small real job once the pixi
environment and scheduler settings are correct.

## Troubleshoot `CUDA_ERROR_UNSUPPORTED_PTX_VERSION`

This error means OpenMM produced or loaded CUDA PTX that the installed NVIDIA
driver cannot JIT-compile. It often appears when the environment resolves a CUDA
toolkit newer than the node or workstation driver supports.

Common fixes:

1. Check the driver on the target GPU node.

   ```bash
   nvidia-smi
   ```

2. Use a pixi environment whose CUDA version does not exceed the driver support.
   For a local workstation where the default `build` environment selected CUDA
   13.2 but the driver supported CUDA 13.1, use the CUDA 12.6 full-stack
   environment instead.

   ```bash
   pixi run -e build-cuda-12-6 python -m openmm.testInstallation
   ```

3. On CUDA 12.4 HPC nodes, use the legacy simulation runtime for simulation
   execution.

   ```bash
   pixi run -e sim-cuda-12-4 python -m openmm.testInstallation
   ```

4. If the test still fails, request a different GPU node or ask the site
   administrators which CUDA driver/runtime combinations are supported.

## See also

- {doc}`hpc_slurm` — submit OpenMM simulation jobs to SLURM
- {doc}`../get_started/installation` — install pixi environments
- {doc}`../reference/cli_reference` — verify which CLI commands and flags exist
- {doc}`troubleshooting` — general PolyzyMD troubleshooting
