# How To: Run Analysis Jobs on a SLURM Cluster

This guide shows you how to run `polyzymd analyze` as a SLURM batch job, so
that trajectory analysis runs on a compute node instead of a login node, and
how to rerun it cheaply once replicate results are stored.

```{note}
This guide covers **analysis** jobs. For submitting **simulation** jobs, see
{doc}`hpc_slurm`.
```

## Before You Start

You need:

- access to a SLURM cluster with `sbatch` available on PATH
- a working analysis pixi environment on the cluster (`pixi install -e analysis`)
- completed simulation trajectories, with one simulation `config.yaml` per
  condition

:::{admonition} Use compute resources, not login nodes
:class: important

Trajectory-analysis jobs can require substantial RAM, CPU time, and scratch
I/O. Run `polyzymd analyze` inside a batch job or an allocated interactive
session, not directly on a login node.
:::

## Write a batch script

`polyzymd analyze` measures every replicate of every condition given with
`-c`, control first, in one process. Put the command in a batch script:

```bash
#!/bin/bash
#SBATCH --job-name=hbonds
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=24:00:00
#SBATCH --output=slurm_logs/%x.%j.out

pixi run -e analysis polyzymd analyze hydrogen_bonds \
    -c noPoly_CALB_pNPB/config.yaml \
    -c SBMA_100_CALB_pNPB/config.yaml \
    --label "No Polymer" --label "SBMA-100" \
    --eq 10ns \
    --format json -o hbonds_report.json \
    --output-dir analysis_results
```

Submit it from the directory that holds the condition folders:

```bash
mkdir -p slurm_logs
sbatch hbonds.sbatch
```

The job writes each replicate's values and their record under
`analysis_results/polyzymd_results/`, the figures under
`analysis_results/figures/hydrogen_bonds/`, and the full `ProtocolReport` to
`hbonds_report.json`. The printed report also goes to the job's log; leave
out `--format json` to print the short agent-format lines instead.

## Rerun from the stored results

A replicate result is read back instead of measured when its record matches
the new call: the same function, selections and settings, the same config,
the same input files, the same equilibration window and the same frames. For
example, `--run protein_polymer_any_fraction` reads back the values that a
run of the default `protein_polymer_mean_hbonds` stored, because both come
from one `functions.hydrogen_bonds` call per replicate. Add `--recompute` to
measure every replicate again.

A replicate whose input files changed since its result was stored, for
example because its trajectory grew, is measured again, so rerunning the same
script while a campaign is still producing trajectories measures only the
replicates that changed.

## One job per analysis

Each analysis is one command, so submit one script per analysis, or run
several commands one after another in the same script:

```bash
pixi run -e analysis polyzymd analyze rmsf \
    -c noPoly_CALB_pNPB/config.yaml -c SBMA_100_CALB_pNPB/config.yaml \
    --eq 10ns --output-dir analysis_results
pixi run -e analysis polyzymd analyze sasa \
    -c noPoly_CALB_pNPB/config.yaml -c SBMA_100_CALB_pNPB/config.yaml \
    --eq 10ns --output-dir analysis_results \
    --set "contexts={isolated: protein, with_polymer: protein or chainid C}" \
    --run with_polymer
```

Give every job the same `--output-dir`, so that each one reuses what the
others stored. `--stride N` measures every N-th production frame, which
shortens a first look at a long campaign; a stored result is reused only for
the same stride.

## Cluster-specific notes (CU Boulder CURC)

CU Boulder's Research Computing provides two SLURM clusters: **Alpine**
(shared campus resource) and **Blanca** (condo model with PI-owned nodes).
Both require `--partition`, `--account`, and `--qos` to be set for job
submission.

### Switching between clusters

Use environment modules to select which cluster's SLURM scheduler you target:

```bash
# Target Blanca (PI-owned condo nodes)
module load slurm/blanca

# Target Alpine (shared campus resource)
module load slurm/alpine
```

Run the appropriate `module load slurm/<cluster>` command **before** `sbatch`.
The module swap points `sbatch`, `squeue`, and the other SLURM utilities at the
selected cluster.

### Required SLURM flags

Both clusters require all three scheduling flags. Omitting any of them causes
`sbatch` to reject the job.

| Flag | Alpine (shared) | Blanca (condo) |
|------|-----------------|----------------|
| `--partition` | `amilan` (CPU), `aa100` or `ami100` (GPU) | `blanca-<group>` (e.g. `blanca-shirts`) |
| `--account` | Your allocation (e.g. `ucb625_asc1`) | Same as partition (e.g. `blanca-shirts`) |
| `--qos` | `normal` | Same as partition (e.g. `blanca-shirts`) |

Add them to the batch script, for example on Blanca:

```bash
#SBATCH --partition=blanca-shirts
#SBATCH --account=blanca-shirts
#SBATCH --qos=blanca-shirts
```

:::{warning}
If you omit `--partition` on Blanca, `sbatch` fails with
*"A partition has not been provided"*, and the error message references
Alpine documentation. This does not mean you need to switch to Alpine. Add
`--partition=blanca-<group>` and resubmit.
:::

:::{tip}
If you are unsure which accounts and partitions you have access to, run:

```bash
sacctmgr show association user=$USER format=account,partition,qos
```

This lists every account/partition/QoS combination available to your user.
:::

## Analysis plugins you register yourself

`polyzymd compare submit` and `submit-all` submit one SLURM job per
replicate, condition and comparison for an analysis plugin registered with the
plugin framework, and `polyzymd compare status` and `finalize` follow and
finish those jobs. No shipped analysis is such a plugin,
so these commands run only plugins you register yourself, and they are being
removed. See `polyzymd compare submit --help` for their options.

## See Also

- {doc}`analysis_agent_protocol` — Every `polyzymd analyze` option, and how to read the report
- {doc}`hpc_slurm` — Submitting simulation jobs to SLURM
- {doc}`analysis_compare_conditions` — Comparing conditions with `polyzymd analyze`
- {doc}`../tutorials/analysis_complete_workflow` — Full local analysis workflow
