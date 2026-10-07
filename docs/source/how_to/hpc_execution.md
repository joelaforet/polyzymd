# How To: Run Analysis Jobs on a SLURM Cluster

This guide shows you how to run `polyzymd analyze` as a SLURM batch job, so
that trajectory analysis runs on a compute node instead of a login node, and
how to rerun it cheaply once replicate results are stored, and how to measure
the replicates in parallel with `--submit`.

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

Run the same script again. PolyzyMD reads back each stored replicate result
whose record matches the new call, and measures only the replicates whose
inputs changed, for example a trajectory that grew. For the rule, see "Why a
cached result is checked against its inputs" in
{doc}`../explanation/analysis_concepts`. Add `--recompute` to measure every
replicate again.

For example, `--run protein_polymer_any_fraction` reads back the values that a
run of the default `protein_polymer_mean_hbonds` stored. Both come from one
`functions.hydrogen_bonds` call per replicate.

## Measure replicates in parallel

A long campaign takes too long to measure in one job. Add `--submit` and a
cluster preset to the command, and PolyzyMD submits one SLURM array task per
condition and replicate and a report job that runs after them:

```bash
pixi run -e analysis polyzymd analyze hydrogen_bonds \
    -c noPoly_CALB_pNPB/config.yaml \
    -c SBMA_100_CALB_pNPB/config.yaml \
    --label "No Polymer" --label "SBMA-100" \
    --eq 10ns --output-dir analysis_results \
    --submit --preset <preset>
```

```text
  condition × replicate jobs (SLURM array)          report job
 ┌──────────────────────────────────────┐
 │ No Polymer, replicate 1  ──┐         │
 │ No Polymer, replicate 2  ──┤         │     ┌──────────────────────────┐
 │ ...                        ├─ store ─┼────▶│ polyzymd analyze, all -c │
 │ SBMA-100,  replicate 1   ──┤ results │     │ reads every stored result│
 │ SBMA-100,  replicate 2   ──┘         │     │ statistics + figures     │
 └──────────────────────────────────────┘     └──────────────────────────┘
                                   --dependency=afterany
```

On a site with several SLURM clusters, load the module of the cluster first.
For CU Boulder, see {doc}`site_cu_boulder`.

The command writes `analysis_results/slurm/hydrogen_bonds_<YYYYmmdd-HHMMSS>/`:

| File | What it holds |
|---|---|
| `tasks.tsv` | One line per array task: config, label and replicate |
| `replicates.sbatch` | The array. Each task runs `polyzymd analyze hydrogen_bonds -c <config> --label <label> --replicates <n>` with the same `--eq`, `--stride`, `--set`, `--run` and `--output-dir` and `--no-plots`, and stores its replicate's result under `analysis_results/polyzymd_results/` |
| `report.sbatch` | The full command with every `-c`, which reads every stored result, measures any replicate a task left unmeasured, draws the figures and writes the report |
| `logs/` | `replicate.<array>_<task>.out` and `report.<job>.out` |
| `report.txt` | The report, written by the report job; `report.json` with `--format json`, or the path given with `-o` |

It submits the array, then the report job with
`--dependency=afterany:<array>`, so the report starts once every task has
ended, whether or not each one succeeded, and prints both job IDs. A task that
fails leaves its replicate unmeasured, and the report job measures it itself;
check `logs/` for the failure.

When a condition still cannot be measured, the report job does not fail with
it: it reports the conditions that work, compared with the control when the
control is among them, and marks the report `status partial`, with one
`problem:` line naming each condition left out and its error. `report.json`
records the same as `status` and `problems`, and `polyzymd study check` shows
a partial report with its problems. An array task, which measures one
condition and replicate, stores its result but never writes the run's
`report.json`.

- **`--dry-run`** writes the folder and prints the two `sbatch` commands
  without submitting, so you can read the scripts first.
- **Environment.** Both jobs run the Python interpreter and `PYTHONPATH` of
  the submitting process, in the directory you submitted from, so they measure
  with the same PolyzyMD and packages and resolve relative `-c` paths as the
  command did. Submit through `pixi run -e analysis` and the jobs use the
  analysis environment.
- **Presets.** `--preset` sets the partition, QoS and account:

  | Preset | Partition | QoS | Account |
  |---|---|---|---|
  | `alpine-cpu` | `acpu` | `cpu-normal` | none |
  | `blanca-shirts` | `blanca-shirts` | `blanca-shirts` | `blanca-shirts` |
  | `blanca-chbe-rdi` | `blanca-chbe-rdi` | `blanca-chbe-rdi` | `blanca-chbe-rdi` |
  | `bridges2-rm` | `RM-shared` | none | none |

  `--partition`, `--account` and `--qos` replace the preset's values, or give
  all three without a preset. Each job gets `--time 12:00:00`, `--mem 16G`
  and `--cpus 2` unless you pass others.

Run on Blanca with two LipA conditions of five replicates each, the submitted
hydrogen-bond analysis gave the same report as one job measuring every
replicate in turn.

### What `--submit` writes, and other schedulers

The submitted scripts do what the hand-written ones below do, and these work
on a cluster no preset covers or adapt to another scheduler. List the
conditions in a file, one config and its label per line, the labels exactly
as in the report script, and write the array script, whose task number picks
a condition and a replicate:

```bash
# conditions.txt
noPoly_CALB_pNPB/config.yaml No Polymer
SBMA_100_CALB_pNPB/config.yaml SBMA-100
```

```bash
#!/bin/bash
#SBATCH --job-name=hbonds-replicates
#SBATCH --array=0-9            # 2 conditions x 5 replicates
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=12:00:00
#SBATCH --output=slurm_logs/%x.%A_%a.out

REPLICATES=5
line=$(( SLURM_ARRAY_TASK_ID / REPLICATES + 1 ))
replicate=$(( SLURM_ARRAY_TASK_ID % REPLICATES + 1 ))
read -r config label < <(sed -n "${line}p" conditions.txt)

pixi run -e analysis polyzymd analyze hydrogen_bonds \
    -c "$config" --label "$label" --replicates "$replicate" \
    --eq 10ns --output-dir analysis_results --no-plots
```

Then submit the report job, the script of [Write a batch script](#write-a-batch-script),
to start when every array task has ended:

```bash
mkdir -p slurm_logs
array=$(sbatch --parsable hbonds_replicates.sbatch)
sbatch --dependency=afterany:$array hbonds.sbatch
```

A stored result is read back only when its record matches, so give every
job the same analysis, `--set` settings, `--label` per condition,
`--eq` and `--stride`. Each task loads its replicate as a single run would,
so an analysis's own requirements hold for every task, such as the
force-field files that {doc}`hydrogen_bonds` needs beside each trajectory
(`<segment>_system.xml` for OpenMM, `prod.tpr` and `<prefix>.top` for
GROMACS).

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

## Cluster settings

Every cluster needs its own partition, account and QoS flags in the
`#SBATCH` lines. `polyzymd analyze --submit --preset` sets them for the
clusters in the preset table under
[Measure replicates in parallel](#measure-replicates-in-parallel). For
simulation jobs, `polyzymd submit --preset` sets them, see {doc}`hpc_slurm`.
For the values on the CU Boulder clusters, see {doc}`site_cu_boulder`.

## See Also

- {doc}`analysis_agent_protocol` — Every `polyzymd analyze` option, and how to read the report
- {doc}`hpc_slurm` — Submitting simulation jobs to SLURM
- {doc}`analysis_compare_conditions` — Comparing conditions with `polyzymd analyze`
- {doc}`../tutorials/analysis_complete_workflow` — Full local analysis workflow
