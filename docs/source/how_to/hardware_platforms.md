# Run OpenMM on Other Hardware

PolyzyMD does not require a particular cluster. The molecular build files and the
OpenMM simulation code are portable. The automatic SLURM routing in version
1.3 has a smaller scope: it supports NVIDIA GPUs whose drivers can use one of
the checked-in CUDA 12.0, 12.4, or 12.6 environments.

Use this guide when you run on another cluster, use a CPU, or evaluate an AMD
GPU.

## Know the three configuration layers

These settings solve different problems:

| Layer | Purpose | Where to change it |
|-------|---------|--------------------|
| SLURM preset | Requests a partition, account, time limit, and GPU type | `src/polyzymd/workflow/slurm.py` or CLI overrides |
| Pixi environment | Supplies a compatible OpenMM and accelerator runtime | `pixi.toml` and `pixi.lock` |
| OpenMM platform | Selects `CUDA`, `OpenCL`, or `CPU` when PolyzyMD creates a Context | `openmm.platform` in `config.yaml` |

A known-site preset selects one runtime for reproducibility. It does not prove
that each node is compatible. The allocated node must still pass the NVIDIA
driver probe and the explicit CUDA Context preflight.

## Current support limits

| Hardware and launch path | Status in version 1.3 |
|--------------------------|-----------------------|
| NVIDIA GPU with `polyzymd submit` | Supported when the driver is compatible with a checked-in CUDA environment |
| NVIDIA GPU on a SLURM cluster without a preset | Supported with a suitable preset or CLI resource overrides |
| Bridges-2 NVIDIA GPU | The scheduler preset is included; test the selected GPU type and current driver before a campaign |
| NVIDIA GPU on this machine with `polyzymd run` | Supported: build in `build`, then run in a `sim-cuda-*` environment |
| CPU with `polyzymd run` | Supported when `openmm.platform` is `CPU` |
| AMD GPU with `polyzymd run` | Possible through OpenMM `OpenCL` when the site supplies a working OpenCL runtime; not tested by PolyzyMD CI |
| CPU or AMD GPU with generated OpenMM SLURM scripts | Not supported in version 1.3; the generated script requires `nvidia-smi` and performs a CUDA preflight |
| GROMACS CPU or GPU with generated SLURM scripts | Supported through a site module or container; GROMACS does not use OpenMM platform routing |

PolyzyMD never changes from CUDA or OpenCL to CPU after a Context failure. An
unavailable requested platform is a fatal error. This rule prevents a large
GPU job from running slowly on CPU without the user's knowledge.

This OpenMM rule does not apply to the GROMACS backend. GROMACS uses its own
site module or container and its own self-restarting SLURM template. See
{doc}`run_gromacs` for portable CPU and GPU recipes.

## Use another NVIDIA SLURM cluster

First, request one interactive GPU. Run the same probe that the generated job
uses:

```bash
nvidia-smi --query-gpu=driver_version,compute_cap --format=csv,noheader
```

Select one validated environment for the campaign. Then generate one test
script. Use an existing preset and override its Slurm fields when possible:

```bash
pixi run -e build polyzymd submit \
    -c config.yaml \
    --preset aa100 \
    --partition <site-partition> \
    --account <site-account> \
    --qos <site-qos> \
    --pixi-env <validated-sim-cuda-environment> \
    --replicates 1 \
    --generate-only
```

Remove an option when the site does not use it. Inspect the generated `#SBATCH`
lines before submission. Run a short test job before you start a campaign.

For presets without a site policy, `auto` selects the newest compatible
checked-in environment. Use an explicit environment after site validation so
all replicates use the same OpenMM version. The node probe does not use a static
node list.

If a node fails the probe or the Context preflight, the job submits one
successor to the same queue and excludes the failed node. It then exits. The
successor uses an `afterany` dependency on the current job. PolyzyMD stops after
three routing retries so an unsupported partition cannot create an endless
submission loop.

## Use Bridges-2

The `bridges2` preset configures the scheduler request and resolves `auto` to
`sim-cuda-12-6`. Select a GPU type that your allocation can use:

```bash
pixi run -e build polyzymd submit \
    -c config.yaml \
    --preset bridges2 \
    --account <allocation> \
    --gpu-type v100-32 \
    --pixi-env auto \
    --replicates 1
```

The allocated node must pass the driver probe and CUDA execution test. A newer
driver does not cause PolyzyMD to select another environment.

## Run on a local NVIDIA GPU

Run on the GPU of a workstation or a laptop in two steps. Build in the `build`
environment, then run in a `sim-cuda-*` environment.

The `build` environment holds OpenFF, so only it can build. Its OpenMM uses a
recent CUDA version. An older NVIDIA driver cannot run that version, and the
run fails with `CUDA_ERROR_UNSUPPORTED_PTX_VERSION`. A `sim-cuda-*`
environment cannot build, because it has no OpenFF, but its OpenMM uses an
older CUDA version.

1. Find the highest CUDA version of your driver. `nvidia-smi` prints it as
   `CUDA Version` in its first lines.
2. Pick the newest `sim-cuda-*` environment that does not exceed that
   version: `sim-cuda-12-6`, `sim-cuda-12-4` or `sim-cuda-12-0`.
3. Set the CUDA platform in the simulation configuration:

   ```yaml
   openmm:
     platform: CUDA
     precision: mixed
   ```

4. Build in `build`, then run in the CUDA environment:

   ```bash
   pixi run -e build polyzymd build -c config.yaml -r 1
   pixi run -e sim-cuda-12-6 polyzymd run -c config.yaml -r 1
   ```

`polyzymd run` finds the finished build and prints
`Reusing the build in <replicate folder>`. It then minimizes, equilibrates
and runs production on the GPU. Give the same `-r` replicates to both
commands. Without a build, `polyzymd run` in a `sim-cuda-*` environment tries
to build and fails, because OpenFF is missing.

## Run on a CPU

Set the platform in the simulation configuration:

```yaml
openmm:
  platform: CPU
  precision: mixed
```

Build and run one replicate on this machine with `polyzymd run`. The `build`
environment holds PolyzyMD and OpenMM:

```bash
pixi run -e build polyzymd run -c config.yaml -r 1
```

`polyzymd run` builds the system, minimizes it, runs the equilibration stages
and runs production, in one process. The shipped example
`examples/quickstart/config.yaml` uses the CPU platform. See
{doc}`../get_started/quickstart`.

The CPU platform uses `SLURM_CPUS_PER_TASK` as the OpenMM thread count when
the variable exists. `polyzymd submit` does not write CPU batch scripts for
OpenMM. To run on CPU nodes of a cluster, write a batch script that requests
CPU resources, activates the environment and runs `polyzymd run`.

## Evaluate an AMD GPU

OpenMM can expose an `OpenCL` platform when the operating system, device
driver, OpenCL loader, and OpenMM package are compatible. Set:

```yaml
openmm:
  platform: OpenCL
  precision: mixed
```

PolyzyMD passes the explicit platform to initial production and continuation.
It fails if OpenMM cannot create that platform. The repository does not include
an AMD Pixi environment, an AMD device probe, or an AMD SLURM template. A site
maintainer must supply and test these parts before production use.

Test the site environment before you use PolyzyMD:

```python
import openmm
from openmm import unit

system = openmm.System()
system.addParticle(1.0)
integrator = openmm.VerletIntegrator(1.0 * unit.femtosecond)
platform = openmm.Platform.getPlatformByName("OpenCL")
context = openmm.Context(system, integrator, platform)
print(context.getPlatform().getName())
```

After this preflight succeeds, run `polyzymd run` as in the CPU section above. Run a short
simulation and compare its initial potential energy with a known result, within
a stated tolerance, before a full campaign.

## Add a supported hardware environment

To add a new NVIDIA driver or GPU cohort to PolyzyMD, follow the checklist in
{ref}`add-hardware-environment` of the contributor guide.

## Respond to a driver update

Do not replace an environment while active replicates use it. Progress metadata
records the Pixi environment, OpenMM version, CUDA runtime, platform, precision,
driver, and device. PolyzyMD rejects a changed simulation environment during a
resubmission.

Use this sequence after a driver update:

1. Run the driver and compute-capability probe on an allocated node.
2. Test each existing Pixi environment with an explicit CUDA Context.
3. Keep an old environment when it still works and active campaigns use it.
4. Add a new rich platform and environment when the new driver needs one.
5. Update the routing threshold only after the Context and benchmark tests pass.
6. Start new campaigns with the new environment. Finish each existing replicate
   with its recorded environment and OpenMM version.

This process makes hardware support explicit and testable. It also keeps the
simulation results independent of cluster names.
