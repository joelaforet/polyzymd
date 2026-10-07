# Analyze the replicates of a study

In this tutorial you run two more replicates of the quickstart simulation and
measure how much each residue of Trp-cage fluctuates (RMSF). You run the
analysis on the {term}`study`, read the report, and find the stored results
and figures. At the end you see when a quick look with `-c` is enough.

You learn these steps:

1. Run more replicates of a condition with `polyzymd run -r`.
2. List an analysis in `project.yaml`.
3. Run it on the study with `polyzymd analyze NAME --study`.
4. Report another result of the same analysis with `--run`.
5. Find the stored results and the figures.

## Before you start

Do {doc}`../get_started/quickstart` first. This tutorial continues in its
project folder, `~/pz_quickstart`.

:::{admonition} Environment Setup
:class: tip

Run every command of this tutorial in the `build` environment. It includes
the analysis tools. From the repository root, activate it once:

```bash
pixi shell -e build
```
:::

## Step 1: Run two more replicates

A mean needs more than one replicate to have an interval. Run replicates 2
and 3 of the `Water` condition:

```bash
cd ~/pz_quickstart
polyzymd run -c trpcage/conditions/water/config.yaml -r 2-3
```

The replicate number seeds the starting structure, so each replicate is a
different simulation. The two runs take about three minutes. The last line
is:

```
All 2 replicate(s) completed successfully.
```

`runs/trpcage/water/` now holds `trpcage_300K_run1/`, `trpcage_300K_run2/`
and `trpcage_300K_run3/`.

## Step 2: List the analysis

Open `project.yaml` and add `rmsf: {}` under `analyses:`:

```yaml
analyses:
  rg: {}
  rmsf: {}
```

`rmsf: {}` runs RMSF with its default settings in every study of the
project. Commit the change:

```bash
git add -A
git commit -m "Add the rmsf analysis"
```

## Step 3: Run RMSF on the study

```bash
polyzymd analyze rmsf --study trpcage
```

`--study` names the study folder. The study gives the conditions, the
control, the replicates and the equilibration window. The project gives the
settings of `rmsf`. The analysis superposes every production frame on a
reference frame of its replicate. Then it measures how far each C-alpha atom
moves about its mean position. The output is:

```
log: /home/me/pz_quickstart/trpcage/logs/polyzymd-analyze-20261006-205942-pid12810.log
# polyzymd analyze rmsf  metric core_rmsf  unit A  run core_rmsf  eq 0ns  conditions 1  replicates 3  protocol rmsf/2
Water  n 3  mean 0.3463  sem 0.01661  ci95 0.2748 to 0.4178  values 0.3669, 0.3585, 0.3134
verdict: Water core_rmsf 0.3463 A (95% CI 0.2748 to 0.4178, n 3)
```

Your values differ a little, because the CPU threads add up the forces in a
different order on each run.

Read the report in this order:

1. The header line starts with `#`. It names the analysis, the result
   (`core_rmsf`), its unit, the equilibration window (`eq 0ns`) and the
   number of replicates.
2. The `Water` line gives the number of replicates (`n 3`), their mean, the
   standard error (`sem`), the 95 % confidence interval (`ci95`) and the
   value of each replicate.
3. `verdict:` gives the result in one line.

`core_rmsf` is one number for each replicate: the root of the mean square
fluctuation of the residues. The simulations are only picoseconds long, so
the values are small. They show the steps, not the physics of Trp-cage.

## Step 4: Report another result

RMSF has several results. To see the value of each residue, report the
`rmsf` result:

```bash
polyzymd analyze rmsf --study trpcage --run rmsf
```

The output has one line for each residue, starting with the residue number.
These are the first lines and the end:

```
# polyzymd analyze rmsf  metric rmsf  unit A  run rmsf  eq 0ns  conditions 1  replicates 3  protocol rmsf/2
1  Water  n 3  mean 0.4587  sem 0.0461  ci95 0.2604 to 0.6571  values 0.4763, 0.3716, 0.5283
2  Water  n 3  mean 0.3254  sem 0.05403  ci95 0.0929 to 0.5578  values 0.4322, 0.2581, 0.2858
...
20  Water  n 3  mean 0.5097  sem 0.05915  ci95 0.2552 to 0.7643  values 0.5066, 0.6137, 0.4089
warning: the 95 percent interval of condition Water at 8 extends past the bounds 0 to inf of rmsf, where a t interval is not reliable
warning: the 95 percent interval of condition Water at 16 extends past the bounds 0 to inf of rmsf, where a t interval is not reliable
verdict: Water rmsf over 20 labels, label means from 0.2236 to 0.5097 A
```

The `warning:` lines are part of the result. With three replicates, the
interval of residues 8 and 16 reaches below 0, which an RMSF cannot be.

This command takes a few seconds. It does not read the trajectories again. It
reuses the values that step 3 stored, because the inputs and the settings are
the same.

## Step 5: Find the results

The results of the study are in `trpcage/results/rmsf/`:

```text
trpcage/results/rmsf/
├── report.json                    # the full report of the last command
├── figures/
│   └── rmsf/
│       ├── rmsf_profile.png       # the RMSF of each residue
│       ├── offset_profile.png
│       ├── rmsd_per_residue_profile.png
│       ├── rms_decomposition.png
│       └── rmsf_comparison.png
└── polyzymd_results/
    └── rms_decomposition/
        └── Water/
            ├── replicate_1/
            │   ├── record.json    # what was measured, on which inputs, with which settings
            │   ├── parts.json
            │   └── values.npz     # the values of this replicate
            ├── replicate_2/
            └── replicate_3/
```

`record.json` holds the function and the hash of its code, the settings, the
hash of the config, the size and SHA-256 of the topology and of each
trajectory file, the frames used and the software versions. For every field
of the report, see {doc}`../reference/analysis_protocol_report`.

## When a quick look with `-c` is enough

You can also give a config instead of a study:

```bash
mkdir ~/quick_look
cd ~/quick_look
polyzymd analyze rmsf -c ~/pz_quickstart/trpcage/conditions/water/config.yaml --label Water --eq 0ns
```

It prints the same report. The results go into the current folder, in
`polyzymd_results/` and `figures/`. No study records them.

Use `-c` for a quick look: to check one simulation while it runs, or to try a
setting before you add it to a study. Use `--study` or `--project` for the
results that you keep. The study then holds the conditions, the control, the
window, the settings and the results together, and `polyzymd project freeze`
can publish them.

## What you did

You ran three replicates of one condition, measured their RMSF on the study,
read the report and found the stored results. Next, compare two conditions in
{doc}`analysis_complete_workflow`.
