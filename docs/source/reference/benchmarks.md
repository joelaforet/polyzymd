# Performance benchmarks

These tables give the measured simulation speed of PolyzyMD systems on cluster
GPUs. Use them to estimate how long a simulation takes.

## How to read the tables

Each GPU has one table. Each row gives a system size in atoms and the speed in
ns/day. Find the row nearest to your system size.

```{tip}
**Estimating your system size:**
- Enzyme only: ~5,000-10,000 atoms (protein) + ~50,000-100,000 atoms (water/ions)
- With polymers: Add ~500-2,000 atoms per polymer chain
- After `polyzymd build`, count the atoms of the built system:
  `grep -c '^\(ATOM\|HETATM\)' solvated_system.pdb`
```

## CU Boulder Alpine Cluster

### NVIDIA A100 (aa100 partition)

| System Size (atoms) | ns/day | Notes |
|--------------------:|-------:|:------|
| ~18,000 | 442 | Measured on Alpine. Blanca also has A100 GPUs. |

### NVIDIA L40 (al40 partition)

No measurement.

## Blanca Cluster

### NVIDIA A40 (blanca-shirts partition)

| System Size (atoms) | ns/day | Notes |
|--------------------:|-------:|:------|
| ~70,000 | 159 | Enzyme + substrate + 2 polymers, NPT ensemble, 2fs timestep |

## Benchmark Conditions

Unless otherwise noted in the "Notes" column, all benchmarks use:

- **Ensemble**: NPT (production phase)
- **Timestep**: 2 fs
- **Thermostat**: Langevin Middle Integrator
- **Barostat**: Monte Carlo (25 step frequency)
- **Precision**: mixed (the PolyzyMD default, `openmm.precision`)
- **Nonbonded method**: PME with 1.0 nm cutoff

## Contributing Benchmarks

To add a benchmark entry:

1. Run a simulation segment (at least 1 ns recommended for accurate timing)
2. Check the SLURM output or simulation logs for performance (ns/day)
3. Count the atoms of `solvated_system.pdb` in the replicate folder
4. Add a row to the appropriate GPU table above

Example log output showing performance:
```
Performance: 159.234 ns/day
```

## Estimating Simulation Time

To estimate how long your simulation will take:

```
Wall time (days) = Simulation length (ns) / Performance (ns/day)
```

**Example**: 100 ns simulation at 159 ns/day:
```
100 ns / 159 ns/day = 0.63 days = ~15 hours
```

```{note}
Actual wall time will be slightly longer due to:
- Energy minimization (~5-10 minutes)
- Equilibration phase
- Checkpoint saving overhead
- Job startup/teardown
```
