# Performance benchmarks

These tables give measured simulation speeds of PolyzyMD systems. Use them to
estimate how long a simulation takes. Each row says what was measured, on
which hardware, with which software, and when, as far as it was recorded.

## How to read the tables

Find the row nearest to your system size and hardware. The speed falls
roughly in proportion to the atom count.

```{tip}
**Estimating your system size:**
- Enzyme only: ~5,000-10,000 atoms (protein) + ~50,000-100,000 atoms (water/ions)
- With polymers: Add ~500-2,000 atoms per polymer chain
- After `polyzymd build`, count the atoms of the built system:
  `grep -c '^\(ATOM\|HETATM\)' solvated_system.pdb`
```

## Cluster GPUs (OpenMM)

| Hardware | System | Atoms | ns/day | Software | Date |
|---|---|---:|---:|---|---|
| NVIDIA A100, CU Boulder Alpine (`aa100`) | not recorded | ~18,000 | 442 | not recorded | not recorded |
| NVIDIA A40, CU Boulder Blanca (`blanca-shirts`) | enzyme, substrate and 2 polymer chains, NPT, 2 fs | ~70,000 | 159 | not recorded | not recorded |

No measurement exists yet for the L40 nodes (`al40`) or for Bridges-2.

## Laptop (OpenMM and GROMACS)

Measured on 2026-10-25 on one laptop: Intel Core Ultra 9 275HX (24 cores),
NVIDIA GeForce RTX 5070 Ti Laptop GPU (driver 591.91), Linux under WSL2,
PolyzyMD 1.3.0rc5. Other runs shared the machine at the same time, so treat
the numbers as lower bounds.

**GROMACS 2026.0** (conda-forge, mixed precision, CPU only, 1 thread-MPI
rank and 2 OpenMP threads, the `gromacs:` settings of the shipped
`config_gromacs.yaml`). The speed is the `Performance` line of `prod.log`.

| System | Atoms | Production steps | ns/day |
|---|---:|---:|---:|
| Quickstart: Trp-cage in water with 0.15 M NaCl | 6,138 | 2,000 | 125 |
| Trp-cage with 2 SBMA chains, water, NaCl (3 replicates) | 11,274 | 25,000 | 68 to 70 |

**OpenMM 8.4** (`polyzymd run -r 1` of the quickstart, 6,138 atoms). The
quickstart production is only 4 ps, too short for an ns/day value, so the
table gives wall times from the run logs.

| Platform | Environment | Build and physics audit | Run |
|---|---|---:|---:|
| CPU | `build` | about 15 s | 42 to 44 s for the whole command, build included (4 runs) |
| CUDA, mixed precision | `sim-cuda-12-6` | 11 s (`polyzymd build`, in `build`) | 6 s (`polyzymd run` after the build) |

The build and audit time runs from the start of `polyzymd run` to the start
of minimization. For a large system the audit takes longer: about one minute
for 90,000 atoms with 78 polymer chains on a cluster node (see
{doc}`../explanation/physics_audit`).

## Benchmark conditions

Unless a row says otherwise, the runs use:

- **Ensemble**: NPT (production)
- **Time step**: 2 fs, bonds to hydrogen constrained
- **Thermostat**: Langevin (`LangevinMiddle`; `sd` on GROMACS)
- **Barostat**: Monte Carlo every 25 steps (stochastic cell rescaling on GROMACS)
- **Precision**: mixed (the PolyzyMD default, `openmm.precision`)
- **Nonbonded**: PME, 0.9 nm cutoff, Lennard-Jones switched from 0.8 nm,
  with dispersion correction. These settings come from the OpenFF force
  field, and every build records them in `build_provenance.md`
  (`nonbonded.settings`).

## Add a benchmark

1. Run a production segment of at least 1 ns.
2. Read the speed in ns/day: `polyzymd status` prints the throughput of a
   running chain, and GROMACS prints a `Performance` line at the end of
   `prod.log`.
3. Count the atoms of `solvated_system.pdb` in the replicate folder.
4. Add a row with the hardware, the system, the atom count, the speed, the
   PolyzyMD and engine versions (from `build_provenance.md`) and the date.

## Estimate the simulation time

```
Wall time (days) = Simulation length (ns) / Performance (ns/day)
```

**Example**: 100 ns at 159 ns/day takes 100 / 159 = 0.63 days, about 15 hours.

```{note}
The wall time is a little longer than this estimate, because of:
- energy minimization (5 to 10 minutes for a large system)
- the equilibration stages
- checkpoint writing
- job startup and the wait between chained SLURM jobs
```
