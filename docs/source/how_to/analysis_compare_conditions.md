# How to Compare Simulation Conditions

Use this guide when you already have completed PolyzyMD simulations and want to
know whether conditions differ, for example an enzyme with and without a
polymer.

You will:

- pick the simulation `config.yaml` of every condition, control first
- run `polyzymd analyze NAME -c ... -c ...` for each analysis
- read the per-condition values and the comparisons against the control
- find the stored results and the figures

```{note}
If you have not yet run a full analysis workflow, start with
[Tutorial: Analyze a Study from Finished Simulations](../tutorials/analysis_complete_workflow.md).
Each analysis has its own quick start: {doc}`analysis_rmsd_quickstart`,
{doc}`analysis_rg_quickstart`, {doc}`analysis_rmsf_quickstart`,
{doc}`analysis_distances_quickstart`, {doc}`analysis_secondary_structure_quickstart`,
{doc}`analysis_sasa_quickstart`, {doc}`analysis_contacts_quickstart`,
{doc}`analysis_native_contacts_quickstart` and {doc}`hydrogen_bonds`. The
catalytic triad is a routine on the analysis API: {doc}`analysis_triad_quickstart`.
```

:::{admonition} Environment Setup
:class: tip

All analysis commands below assume you have activated the PolyzyMD analysis
pixi environment:

```bash
pixi shell -e analysis
```

Alternatively, prefix each command with `pixi run -e analysis`.
:::

:::{admonition} Resource requirements
:class: important

`polyzymd analyze` loads trajectories and can require substantial RAM, CPU
time, and scratch I/O. On shared HPC systems, run it inside an allocated job or
interactive compute session, not on a login node; {doc}`hpc_execution` shows a
batch script.
:::

## Before You Start

Make sure each condition already has:

- a simulation `config.yaml`
- finished trajectories for the replicates you want to compare

`polyzymd analyze` stores each replicate's result under `polyzymd_results/`
and reads it back on a later run when the function, settings, config, input
files, equilibration window and frames are unchanged. A replicate whose
trajectory grew since its result was stored is measured again. Add
`--recompute` to measure every replicate again.

## Analyze a campaign that is still running

Analyses leave out any OpenMM production segment that `progress.json` records
as running or failed, because its trajectory file ends wherever the last flush
landed. A warning names the segments that were dropped.

- If the warning says the window ends before the excluded segments, the results
  describe the completed part of the run. Wait for the run to finish and rerun
  with `--recompute` when you want the full window.
- If it says the exclusion left a gap, an excluded segment sits between two
  that were kept and the lineage check will refuse the replicate. Wait for the
  run to finish, or load it deliberately with `require_complete=False`:

  ```python
  from polyzymd.analyses.shared.loader import TrajectoryLoader
  u = TrajectoryLoader(config).load_universe(replicate=1, require_complete=False)
  ```

  The incomplete segments are then read as they stand and listed in
  `incomplete_segments` on the layout.

On GROMACS there is no per-segment status to consult, so a production XTC that
is still being written is read as-is. Check that the job has finished before
analyzing a live GROMACS run.

## Step 1: Run One Comparison

Give one `-c` per condition, the control first, and name the conditions with
`--label` in the same order. Protein-polymer hydrogen bonds exist only where
there is a polymer, so this example compares three polymer conditions:

```bash
polyzymd analyze hydrogen_bonds \
  -c ../SBMA_100_enzyme_DMSO/config.yaml \
  -c ../EGMA_100_enzyme_DMSO/config.yaml \
  -c ../SBMA_50_enzyme_DMSO/config.yaml \
  --label "100% SBMA" --label "100% EGMA" --label "50% SBMA" \
  --replicates 1-3 --eq 10ns --set d_a_cutoff=3.0
```

This command:

- builds every condition from its `config.yaml` and the replicates given with
  `--replicates`, or every replicate found on disk without it
- discards the first 10 ns of every replicate
- measures, on every production frame, the hydrogen bonds between the
  protein (`chainid A`) and the polymer (`chainid C`), the default groups
- prints each condition's replicate mean with its 95% interval, and a Welch's
  t test of every condition against the control with the Benjamini-Hochberg
  correction

Without `--label`, each condition is named after the directory holding its
config. {doc}`analysis_agent_protocol` explains every line of the output, and
{doc}`hydrogen_bonds` lists the other settings, such as named groups and
summaries.

## Step 2: Pick Another Result

An analysis that reports several results reports one at a time and lists the
others in `all_runs` of the JSON report. Pick one with `--run`:

```bash
polyzymd analyze hydrogen_bonds -c A/config.yaml -c B/config.yaml --eq 10ns \
  --run protein_polymer_residues
```

The replicate results stored by the first command are read back, so only
results that were not measured yet load the trajectories.

## Step 3: Save the Full Report

`--format json` prints the whole `ProtocolReport`, and `-o` also writes it to a
file:

```bash
polyzymd analyze hydrogen_bonds -c A/config.yaml -c B/config.yaml --eq 10ns \
  --format json -o reports/hydrogen_bonds.json
```

## Step 4: Run the Other Analyses

Every analysis takes the same `-c`, `--label`, `--replicates` and `--eq`
options, so run each one on the same conditions:

```bash
polyzymd analyze rmsf -c A/config.yaml -c B/config.yaml --eq 10ns
polyzymd analyze contacts -c A/config.yaml -c B/config.yaml --eq 10ns \
  --set polymer_selection='chainid C' --set protein_selection='chainid A'
polyzymd analyze distances -c A/config.yaml -c B/config.yaml --eq 10ns \
  --set pairs=pairs.yaml
```

For the pair distances of a catalytic triad, write the pairs to `pairs.yaml`;
for its hydrogen bonds, follow {doc}`analysis_triad_quickstart`.

## Step 5: Check the Outputs

After a successful run in `polymer_stability_study/`, expect files like these:

```text
polymer_stability_study/
├── polyzymd_results/
│   └── hydrogen_bonds_protein_polymer/
│       ├── 100_SBMA/
│       │   ├── replicate_1/
│       │   │   ├── record.json
│       │   │   └── values.npz
│       │   └── ...
│       └── 100_EGMA/
│           └── ...
└── figures/
    └── hydrogen_bonds/
        └── hbonds_protein_polymer_mean_hbonds_comparison.png
```

Each folder name is the result or condition label with every run of
characters other than letters, digits, `.`, `+` and `-` replaced by `_`.
`record.json` holds what was measured and on which inputs, and `values.npz`
the replicate's values. `--output-dir` puts `polyzymd_results/` and
`figures/` in another directory, and `--no-plots` draws no figures.

## Programmatic Use

In Python, `analyze` returns the same report:

```python
from polyzymd.analyses import analyze

report = analyze(
    "hydrogen_bonds",
    ["A/config.yaml", "B/config.yaml"],
    labels=["100% SBMA", "100% EGMA"],
    equilibration="10ns",
)
print(report.to_agent_text())
```

To measure a function of your own on the same conditions, build a
`pz.Study` and call `study.timeseries` or `study.per_replicate`; see
{doc}`../reference/analysis_functions` and {doc}`custom_artifact_plotting`.

## If you have a `comparison.yaml`

`polyzymd analyze` does not read `comparison.yaml`. `polyzymd analyze NAME -f
comparison.yaml` exits with an error that prints the equivalent
`polyzymd analyze NAME -c <config> --label <label> ... --replicates ... --eq ...`
command built from the file's conditions, labels, replicates and
equilibration window. `polyzymd compare`, with any arguments, exits 2 and
prints the `polyzymd analyze` and `polyzymd analyze ... --submit` commands
that replace it.

## Troubleshooting

### `config` path not found

Relative `-c` paths are resolved from your current shell directory.

### `no run directory under the scratch directory`

The replicates given with `--replicates` have no run directory for that
config. Leave out `--replicates` to use every replicate found on disk.

### `the control ... has no replicate where every selection matches atoms`

A condition without a polymer has no `chainid C`, so `hydrogen_bonds` and
`contacts` leave its replicates out, and when it is the control the other
conditions are only summarised. Put a condition with a polymer first to compare
against it, or compare a group that every condition has, such as
`--set "summaries={protein: {within: protein}}"`.

## See Also

- [Tutorial: Analyze a Study from Finished Simulations](../tutorials/analysis_complete_workflow.md)
- [Get a validated number with one command](analysis_agent_protocol.md)
- [Statistical Best Practices for Analysis](../explanation/analysis_statistics_best_practices.md)
