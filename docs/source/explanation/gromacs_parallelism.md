# MPI and OpenMP threads in GROMACS

GROMACS has two kinds of MPI build. The `gromacs:` block of `config.yaml`
sets each one in a different way, so it helps to know which build your site
gives you.

| Concept | Thread-MPI | Real MPI |
|---------|-----------|----------|
| **Binary** | `gmx` | `gmx_mpi` |
| **Launch** | Direct (`gmx mdrun -ntmpi N`) | Through a launcher (`mpirun -np N gmx_mpi mdrun`) |
| **GPU support** | Yes (CUDA thread-MPI) | Depends on the build (often no GPU) |
| **Multi-node** | No (one node only) | Yes |
| **Config field** | `gromacs.ntmpi` | `gromacs.ntmpi` + `gromacs.mpi_launcher_flags` |

PolyzyMD reads the binary name. A name that contains `_mpi` is a real-MPI
build. PolyzyMD starts it with `mpirun`, unless `command_prefix` is set. Any
other name is a thread-MPI build. PolyzyMD passes `-ntomp` to both builds,
and `-ntmpi` only to a thread-MPI build. It does not add a flag that
`mdrun_flags` already has.

## How many ranks and threads

Each MPI rank runs `ntomp` OpenMP threads, so:

```
ntmpi × ntomp = total CPU cores allocated
```

PolyzyMD asks SLURM for `ntmpi` tasks (or `slurm_ntasks`) and `ntomp` CPUs
per task.

**GPU runs** usually use thread-MPI with one rank for each GPU and many
OpenMP threads:

```yaml
gromacs:
  gmx_binary: "gmx"       # thread-MPI build
  ntmpi: 1                 # 1 rank = 1 GPU
  ntomp: 12                # 12 OpenMP threads
```

**CPU runs** can use either build. A run on more than one node needs real
MPI:

```yaml
gromacs:
  gmx_binary: "gmx_mpi"   # real MPI build
  ntmpi: 8                 # 8 MPI ranks
  ntomp: 1                 # 1 OpenMP thread per rank
```

On most clusters, the thread-MPI `gmx` binary is the best choice for a GPU
run on one node. Use `gmx_mpi` only when you need more than one node.

## See also

- {doc}`../how_to/run_gromacs`: run and submit GROMACS jobs
- {ref}`config-gromacs`: every field of the `gromacs:` block
