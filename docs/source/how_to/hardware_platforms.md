# Run OpenMM on Other Hardware

PolyzyMD does not require a Blanca cluster. The molecular build files and the
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
| NVIDIA GPU on a non-Blanca SLURM cluster | Supported with a suitable preset or CLI resource overrides |
| Bridges-2 NVIDIA GPU | The scheduler preset is included; test the selected GPU type and current driver before a campaign |
| CPU with the Python API | Supported when `openmm.platform` is `CPU` |
| AMD GPU with the Python API | Possible through OpenMM `OpenCL` when the site supplies a working OpenCL runtime; not tested by PolyzyMD CI |
| CPU or AMD GPU with generated OpenMM SLURM scripts | Not supported in version 1.3; the generated script requires `nvidia-smi` and performs a CUDA preflight |
| GROMACS CPU or GPU with generated SLURM scripts | Supported through a site module or container; GROMACS does not use OpenMM platform routing |

PolyzyMD never changes from CUDA or OpenCL to CPU after a Context failure. An
unavailable requested platform is a fatal error. This rule prevents a large
GPU job from running slowly on CPU without the user's knowledge.

This OpenMM rule does not apply to the GROMACS backend. GROMACS uses its own
site module or container and its own self-restarting SLURM template. See
{doc}`gromacs_export` for portable CPU and GPU recipes.

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
all replicas use the same OpenMM version. The node probe does not use a static
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

## Run with CPU through the API

Set the platform in the simulation configuration:

```yaml
openmm:
  platform: CPU
  precision: mixed
```

Run the OpenMM engine from a site-managed environment that contains PolyzyMD
and OpenMM:

```python
from pathlib import Path

from polyzymd.config.schema import SimulationConfig
from polyzymd.engines import create_engine

config = SimulationConfig.from_yaml("config.yaml")
engine = create_engine(config)
replicate = 1
working_dir = Path(config.get_working_directory(replicate))
engine.run_local(replicate, working_dir, skip_build=True)
```

The CPU platform uses `SLURM_CPUS_PER_TASK` as the OpenMM thread count when the
variable exists. The Python API does not generate or submit a CPU batch script.
Write a site wrapper that requests CPU resources, activates the site
environment, and runs this Python entry point.

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

After this preflight succeeds, use the Python API example above. Run a short,
deterministic simulation and compare its energies with a known result before a
full campaign.

## Add a supported hardware environment

Use this checklist when a site has a new NVIDIA driver or GPU cohort:

1. Record the node, GPU model, driver, maximum CUDA version, and compute
   capability from an allocated node.
2. Add a named rich platform to the `workspace.platforms` list in `pixi.toml`.
3. Add or update a simulation feature, environment, and solve group in
   `pixi.toml`. Keep the OpenMM version fixed for an active campaign.
4. Regenerate `pixi.lock` with Pixi 0.72.2 or newer. Install the environment on
   an allocated node.
5. Add the validated driver threshold and environment name in
   `src/polyzymd/workflow/cuda_routing.py` and in the OpenMM Slurm template.
6. Add routing tests in `tests/workflow/test_cuda_routing.py` and script tests
   in `tests/workflow/test_slurm.py`.
7. Create an explicit CUDA Context on each hardware cohort. Then run a short,
   deterministic benchmark without CPU fallback.
8. Compare particle identity and energy behavior with an existing supported
   environment. Record the test tolerance and result.
9. Add the verified hardware to this guide. Do not infer support from a node
   name or a CUDA toolkit module.

To add a reusable scheduler preset, edit `SlurmConfig.from_preset()` and
`PRESET_DEFAULT_PIXI_ENV` in `src/polyzymd/workflow/slurm.py`. Assign the
validated environment as the preset default. Add its CLI choice and tests. Use
CLI overrides when the difference is only an account, partition, QoS,
constraint, or node selection.

## Respond to a driver update

Do not replace an environment while active replicas use it. Progress metadata
records the Pixi environment, OpenMM version, CUDA runtime, platform, precision,
driver, and device. PolyzyMD rejects a changed simulation environment during a
resubmission.

Use this sequence after a driver update:

1. Run the driver and compute-capability probe on an allocated node.
2. Test each existing Pixi environment with an explicit CUDA Context.
3. Keep an old environment when it still works and active campaigns use it.
4. Add a new rich platform and environment when the new driver needs one.
5. Update the routing threshold only after the Context and benchmark tests pass.
6. Start new campaigns with the new environment. Finish each existing replica
   with its recorded environment and OpenMM version.

This process makes hardware support explicit and testable. It also keeps the
simulation results independent of cluster names.
