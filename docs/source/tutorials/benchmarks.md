# Performance Benchmarks

This page contains performance benchmarks for PolyzyMD simulations on various GPU hardware available at CU Boulder. Use these benchmarks to estimate simulation throughput and plan your computational resource needs.

## How to Read These Tables

Each GPU has a table mapping **system size (atoms)** to **simulation speed (ns/day)**. Find your approximate system size and read across to estimate performance.

```{tip}
**Estimating your system size:**
- Enzyme only: ~5,000-10,000 atoms (protein) + ~50,000-100,000 atoms (water/ions)
- With polymers: Add ~500-2,000 atoms per polymer chain
- Use `polyzymd build --dry-run` to see exact atom counts before running
```

## CU Boulder Alpine Cluster

### NVIDIA A100 (aa100 partition)

| System Size (atoms) | ns/day | Notes |
|--------------------:|-------:|:------|
| | | |

### NVIDIA L40 (al40 partition)

| System Size (atoms) | ns/day | Notes |
|--------------------:|-------:|:------|
| | | |

## Blanca Cluster (Shirts Lab)

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
- **Precision**: Mixed precision (OpenMM default)
- **Nonbonded method**: PME with 1.0 nm cutoff

## Contributing Benchmarks

To add a benchmark entry:

1. Run a simulation segment (at least 1 ns recommended for accurate timing)
2. Check the SLURM output or simulation logs for performance (ns/day)
3. Note the system size from the build output or `system.pdb`
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
