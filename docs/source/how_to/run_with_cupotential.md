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
- glibc 2.38 or newer and a CUDA 13.4 driver. On CURC, run inside the
  `openmm-cupotential` Apptainer image.

## Get the plugin

Download `nv-openmm-plugin-conda-channel-0.3.5.tar.gz` from the
[cuPotential releases](https://github.com/NVIDIA/cuPotential/releases), unpack
it, and install `nv-openmm-plugin 0.3.5` from the unpacked folder as a local
conda channel. Without the plugin, a `openmm_cupotential` run stops with an error
that names these steps.

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

## GPU support and benchmarks

To be filled.
