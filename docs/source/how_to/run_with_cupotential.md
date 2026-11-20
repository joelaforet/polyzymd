# Run OpenMM dynamics on NVIDIA cuPotential

NVIDIA's cuPotential plugin (`nv_openmm_plugin`) replaces OpenMM's nonbonded
force, LangevinMiddle integrator and Monte Carlo barostat with faster CUDA
versions. PolyzyMD can run equilibration and production dynamics on it.

## Enable it

```yaml
engine: openmm
openmm:
  platform: CUDA
  precision: mixed   # or single; double is refused
  backend: openmm_cupotential
```

The backend is a runtime choice: it does not change the build or the config
hash, so an existing prebuilt bundle runs unchanged. The default, `openmm`,
runs stock OpenMM for everything.

## Requirements

- The `CUDA` platform with `mixed` or `single` precision.
- One `NonbondedForce` with PME or LJPME. PolyzyMD sets its Ewald tolerance
  to 5e-4 and lets the plugin choose the PME parameters, as the plugin requires.
- The `LangevinMiddle` thermostat and, under NPT, the isotropic Monte Carlo
  barostat. Other thermostats and barostats are refused.
- glibc 2.38 or newer and a CUDA 13.4 driver API. On CURC, jobs run inside
  the `openmm-cupotential` Apptainer image (see [GPU support](#gpu-support)).

## Get the plugin and the image

The plugin runs in its own pixi workspace, `envs/openmm-cupotential/pixi.toml`,
which installs PolyzyMD from your checkout. It is not part of the main
`pixi.toml`, so `pixi lock` keeps working on machines without the plugin.

1. Download `nv-openmm-plugin-conda-channel-0.3.5.tar.gz` from the
   [cuPotential releases](https://github.com/NVIDIA/cuPotential/releases) and
   unpack it. The unpacked folder holds `linux-64/repodata.json`; it is a
   local conda channel.
2. Point PolyzyMD at it and at the Apptainer image:

   ```bash
   export POLYZYMD_CUPOTENTIAL_CHANNEL=/path/to/nv-openmm-plugin-conda-channel-0.3.5
   export POLYZYMD_CUPOTENTIAL_IMAGE=/path/to/openmm-cupotential_0.3.5.sif
   ```

   With `--preset blanca-shirts`, both default to the Shirts lab copies under
   `/pl/active/shirts_archive/LaforetJoe/containers/openmm-cupotential/`.
   Without an image or a channel, `polyzymd submit` stops before it submits
   anything and names the variable to set.

`polyzymd submit` links the channel to `envs/openmm-cupotential/nv-channel`
(a pixi manifest cannot read environment variables). To lock the environment
by hand on a login node:

```bash
ln -sfn "$POLYZYMD_CUPOTENTIAL_CHANNEL" envs/openmm-cupotential/nv-channel
pixi lock --manifest-path envs/openmm-cupotential/pixi.toml
```

The lock file is not committed: pixi records the channel by the absolute path
of your checkout, so a lock from one machine does not work on another. The
environment installs only inside the image, which has glibc 2.39.

Without the plugin, a `openmm_cupotential` run stops with an error that
names these steps.

## What stays the same

- `system.xml` and every segment and interrupt `*_system.xml` hold the stock
  OpenMM System. PolyzyMD converts a private copy only when it creates a
  Context. Energy minimization always runs on stock forces.
- The plugin integrator removes centre-of-mass motion at the frequency of
  the stock System's `CMMotionRemover` (every step for Interchange systems),
  never at the plugin default of 100 steps.
- A run never mixes backends. Each phase records `backend` in its provenance,
  and a later job with another backend stops. A record without `backend` ran
  `openmm`.

## GPU support

`polyzymd submit` runs an `openmm_cupotential` job in its own environment,
`openmm-cupotential`, and records that name in `runtime_platform.json`. The
CUDA 12 environments (`sim-cuda-12-*`) and their driver-based routing are not
used, and a run never moves between them. `--pixi-env` with another
environment is refused.

Each job:

1. Asks Blanca for capable GPUs only: `#SBATCH --constraint=A40|A100|L40|h100|rtx6000|T4`
   on the Blanca presets. A `--constraint` you pass replaces this list. V100
   and P100 nodes are never requested: CUDA 13 does not support them. One L40S
   node has no feature tag and is missed.
2. Runs `cupo-run` in the image. It uses NVIDIA's CUDA 13.4 forward-compatible
   libraries when the host driver accepts them, else the host driver. When
   neither gives a CUDA 13.4 driver API, it exits 86, and the job excludes the
   node and submits itself again, up to three times, like any other routing
   failure. A cuPotential job never falls back to stock OpenMM: it waits for a
   capable node instead.
3. Installs `envs/openmm-cupotential` inside the image on first use. Jobs that
   start together wait on a lock (`.install.lock`). The pixi cache goes to
   `envs/openmm-cupotential/.pixi-cache`, or `$POLYZYMD_CUPOTENTIAL_PIXI_CACHE`,
   never `$HOME`: the home quota on CURC is 2 GB. Keep the checkout in
   `/projects`, not `$HOME`, for the same reason.
4. Runs the CUDA check, `run-segment` and `check-progress` with
   `apptainer exec --nv`. `/projects/$USER`, `/scratch/alpine/$USER`, the run
   folder, the config folder, the checkout and the channel are bound; set
   `APPTAINER_BIND` for any other folder the job reads. Apptainer passes
   `SIGUSR1` and `SIGTERM` on to PolyzyMD, so wall-time stops and preemption
   work as in a stock job.

Which NVIDIA driver branches work:

| Host driver | Result |
|-------------|--------|
| R615 or newer | Works with the host driver (tested: 620 on Blanca's stateless image) |
| R570, R580, R590 | Works through the forward-compatible libraries (tested on Blanca and Alpine) |
| R595, R610 | Should work through the forward-compatible libraries (listed by NVIDIA; untested here) |
| R575 | Exit 86: not a forward-compatibility branch; OpenMM cannot load its kernels |
| R550, R525 | Exit 86: the node is excluded |
| Any driver on V100 or P100 | Never requested |

Measured on Blanca (2026-11-20; 55k-atom enzyme in water, NPT, 2 fs, mixed precision),
stock OpenMM 8.4 vs `openmm_cupotential` in ns/day: L40/L40S (R590) 632–646 → 1,798–1,836
(2.8x); H100 NVL (R580) 498 → 1,336 (2.7x); H100 PCIe (R570) 455 → 1,199 (2.6x);
A100 40GB (R590) 404 → 1,018 (2.5x); Quadro RTX 6000 (R620) 299 → 563 (1.9x).

Other branches are untested; `cupo-run` decides on each node when the job
starts.

Apptainer runs only on compute nodes, not on login nodes.
