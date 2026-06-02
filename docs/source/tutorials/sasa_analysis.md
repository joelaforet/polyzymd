# Tutorial: Measure Polymer Shielding with SASA

This tutorial walks through one guided SASA (Solvent Accessible Surface Area)
workflow: compare an enzyme without polymer to polymer-conjugated conditions and
interpret whether the polymer shields the protein surface.

By the end, you will have:

- configured paired SASA runs for an isolated and polymer-aware protein surface,
- run the `sasa` analysis plugin on a comparison study,
- checked the expected output files, and
- interpreted the main shielding signal in the comparison plots.

For task recipes, use {doc}`../how_to/analysis_sasa_quickstart`. For field and
artifact lookup, use {doc}`../reference/analysis_sasa_reference`.

## Prerequisites

Before starting, make sure you have:

1. a working PolyzyMD pixi environment,
2. completed production trajectories for at least two conditions,
3. a `comparison.yaml` with those conditions, and
4. the PolyzyMD chain convention in mind: A = protein, B = substrate,
   C = polymer, D+ = solvent/ions/other.

If you have not run a comparison analysis before, complete
{doc}`first_analysis` first.

## What SASA will tell us

SASA reports how much molecular surface is accessible to solvent. The SASA
plugin lets each run define two selections:

| Selection | Role in this tutorial |
|-----------|-----------------------|
| `target_selection` | the atoms whose SASA is reported |
| `context_selection` | the atoms allowed to block the target surface |

The shielding idea is simple:

1. Compute protein SASA with only protein atoms in the context.
2. Compute protein SASA with protein plus polymer atoms in the context.
3. Compare the two values.

If the polymer-aware run has lower protein SASA, the polymer is shielding part
of the protein surface.

## Step 1: Add paired SASA runs

Add a `sasa` section to your `comparison.yaml` under `plugins:`.

```yaml
plugins:
  sasa:
    runs:
      - label: "protein_isolated"
        target_selection: "protein"
        context_selection: "protein"

      - label: "protein_with_polymer"
        target_selection: "protein"
        context_selection: "protein or chainid C"

    probe_radius_nm: 0.14
    n_sphere_points: 960
```

The first run is the baseline. The second run lets the polymer on chain C block
protein surface points.

:::{tip}
For a no-polymer control condition, `protein_with_polymer` should be nearly the
same as `protein_isolated` because there are no chain C atoms. Treat this as a
useful consistency check.
:::

## Step 2: Use the runs in a comparison file

Here is a minimal two-condition example. Replace the config paths and labels
with your own completed simulations.

```yaml
name: "calb_sasa_shielding"
description: "SASA shielding tutorial"
control: "No Polymer"

conditions:
  - label: "No Polymer"
    config: "../no_polymer/config.yaml"
    replicates: [1, 2, 3]

  - label: "With Polymer"
    config: "../with_polymer/config.yaml"
    replicates: [1, 2, 3]

defaults:
  equilibration_time: "10ns"

plugins:
  sasa:
    runs:
      - label: "protein_isolated"
        target_selection: "protein"
        context_selection: "protein"
      - label: "protein_with_polymer"
        target_selection: "protein"
        context_selection: "protein or chainid C"
    probe_radius_nm: 0.14
    n_sphere_points: 960

plot_settings:
  format: "png"
  dpi: 300
  style: "compact"
```

## Step 3: Run the SASA analysis

Run the plugin from the directory containing `comparison.yaml`:

```bash
pixi run -e build polyzymd compare run sasa -f comparison.yaml
```

PolyzyMD computes each replicate, aggregates results by condition, compares the
conditions, and writes plots.

:::{note}
SASA can be slower than RMSD or Rg because each frame must evaluate surface
points around atoms in the context selection. For large studies, use the SLURM
workflow in {doc}`../how_to/hpc_execution` rather than copying HPC commands
from this tutorial.
:::

## Step 4: Check that the run succeeded

You should see these outputs:

```text
analysis/<condition>/sasa/run_<replicate>/result.json
analysis/<condition>/sasa/aggregated/result.json
comparison/sasa/result.json
figures/sasa/
```

The comparison file is the first place to check after a successful run:

```bash
python - <<'PY'
import json
from pathlib import Path

result = json.loads(Path("comparison/sasa/result.json").read_text())
print(result["run_labels"])
for condition in result["conditions"]:
    print(condition["label"])
    for run in condition["run_summaries"]:
        print(f"  {run['label']}: {run['mean_sasa']:.1f} ± {run['sem_sasa']:.1f} A^2")
PY
```

Success looks like one summary for each condition and one entry for each SASA
run label.

## Step 5: Interpret the shielding signal

Start with the `protein_isolated` run. It asks whether the protein itself has a
similar accessible surface across conditions before polymer blocking is counted.
Large differences here can mean the protein compactness or conformation differs
between conditions.

Then inspect `protein_with_polymer`. In polymer conditions, a lower mean SASA
relative to the no-polymer control is the main shielding signal.

The strongest evidence for polymer shielding is:

1. `protein_isolated` remains similar across conditions, and
2. `protein_with_polymer` decreases in polymer-conjugated conditions.

You can also compare the two runs within a polymer condition:

```text
shielded area = protein_isolated mean - protein_with_polymer mean
```

A larger positive difference indicates that more protein surface is blocked by
the polymer context.

## Step 6: Use the plots

Open the SASA plots in the configured plot output directory, usually
`figures/sasa/`. The most useful first checks are:

- `sasa_comparison_protein_with_polymer.png` — mean SASA by condition, with
  replicate scatter points.
- `sasa_normalized_comparison_protein_with_polymer.png` — percent change for
  each non-control condition relative to the configured control.

Negative normalized values indicate lower SASA than the control, which is
consistent with shielding for this tutorial's `protein_with_polymer` run.

## What you have now

You have completed a guided SASA shielding analysis and can now answer:

- Did the polymer reduce protein solvent-accessible surface area?
- Was the reduction specific to the polymer-aware context?
- Which plots should you inspect first for replicate consistency and effect
  direction?

## Next steps

- Use {doc}`../how_to/analysis_sasa_quickstart` for active-site,
  monomer-specific, stride, and quick-command recipes.
- Use {doc}`../reference/analysis_sasa_reference` for settings, canonical
  artifact paths, plot filenames, and programmatic artifact loading.
- Use {doc}`../how_to/hpc_execution` when the SASA workload needs SLURM.
