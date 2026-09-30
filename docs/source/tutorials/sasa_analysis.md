# Tutorial: Measure Polymer Shielding with SASA

This tutorial walks through one guided SASA (solvent-accessible surface area)
workflow: compare an enzyme without polymer to polymer-conjugated conditions and
interpret whether the polymer shields the protein surface.

By the end, you will have:

- measured the protein's SASA with and without polymer in the calculation,
- compared each condition with the no-polymer control,
- computed the area the polymer covers, frame by frame, and compared it, and
- found the figures to check first.

For task recipes, use {doc}`../how_to/analysis_sasa_quickstart`. For what each
function measures, use {doc}`../reference/analysis_functions`.

## Prerequisites

Before starting, make sure you have:

1. a working PolyzyMD pixi environment,
2. completed production trajectories for at least two conditions,
3. the `config.yaml` of each condition, and
4. the residue names of the polymer's monomers in your topology.

If you have not run an analysis before, complete {doc}`first_analysis` first.

:::{admonition} Environment Setup
:class: tip

All analysis commands below assume you have activated the PolyzyMD analysis
pixi environment:

```bash
pixi shell -e analysis
```

Alternatively, prefix each command with `pixi run -e analysis`.
:::

## What SASA will tell us

SASA is the area of a molecule's surface that a solvent-sized probe can reach.
Each SASA calculation has two selections:

| Selection | Role in this tutorial |
|-----------|-----------------------|
| `target` | the atoms whose SASA is reported, here the protein |
| context | the atoms present in the calculation, which can cover the target's surface |

The shielding idea is simple:

1. Compute protein SASA with only protein atoms in the context.
2. Compute protein SASA with protein and polymer atoms in the context.
3. Compare the two values.

If the protein's SASA is lower with polymer in the context, the polymer is
covering part of the protein surface.

## Step 1: Measure the protein on its own

```bash
polyzymd analyze sasa \
  -c ../noPoly_enzyme/config.yaml -c ../SBMA_100_enzyme/config.yaml \
  --label "No Polymer" --label "100% SBMA" --eq 200ns \
  --set "contexts={isolated: protein, with_polymer: protein or resname SBM EGM}" \
  --run isolated
```

Replace `SBM EGM` with the residue names of your polymer. The report gives each
condition's mean protein SASA with its 95 percent interval and every replicate
value, and compares `100% SBMA` with `No Polymer`.

`isolated` asks whether the protein itself has a similar accessible surface
across conditions before polymer covering is counted. A large difference here
can mean the protein's compactness or conformation differs between conditions.

## Step 2: Measure the protein with polymer present

Run the same command with `--run with_polymer`:

```bash
polyzymd analyze sasa \
  -c ../noPoly_enzyme/config.yaml -c ../SBMA_100_enzyme/config.yaml \
  --label "No Polymer" --label "100% SBMA" --eq 200ns \
  --set "contexts={isolated: protein, with_polymer: protein or resname SBM EGM}" \
  --run with_polymer
```

In the no-polymer condition the context adds no atoms, so `with_polymer` equals
`isolated` there. In polymer conditions, a lower mean than the control is the
shielding signal.

The strongest evidence for polymer shielding is:

1. `isolated` stays similar across conditions, and
2. `with_polymer` decreases in the polymer conditions.

## Step 3: Compare the covered area

The area the polymer covers is `isolated` minus `with_polymer`, frame by frame.
In Python, `Timeseries.transform` computes it from the two stored series without
reading the trajectories again:

```python
import numpy as np
import polyzymd as pz
from polyzymd.analyses.functions import sasa

study = pz.Study.from_configs(
    {"No Polymer": "../noPoly_enzyme/config.yaml", "100% SBMA": "../SBMA_100_enzyme/config.yaml"},
    equilibration="200ns",
)
protein = pz.select("protein")
isolated = study.timeseries(sasa, protein, protein, unit="A^2", name="sasa_isolated")
with_polymer = study.timeseries(
    sasa, protein, pz.select("protein or resname SBM EGM"), unit="A^2", name="sasa_with_polymer"
)
covered = isolated.transform(np.subtract, with_polymer, unit="A^2", name="sasa_covered")
print(covered.reduce("mean").compare(control="No Polymer").to_agent_text())
```

`sasa_isolated` and `sasa_with_polymer` are the names `polyzymd analyze sasa`
stores its results under, so run in the same folder, the two series are read
back rather than measured again. A larger positive covered area means more
protein surface is covered by polymer.

## Step 4: Use the figures

`polyzymd analyze sasa` writes its figures to `figures/sasa/`. The most useful
first checks are:

- `sasa_comparison_with_polymer` — each condition's mean SASA with its interval
  and every replicate value.
- `sasa_timeseries_with_polymer` — every replicate's SASA against time, to see
  whether it settles after the equilibration window.

To see which residues are covered, run `--run with_polymer_residues`, which
compares each residue's SASA with the control and draws
`sasa_profile_with_polymer` and `sasa_difference_with_polymer`.

## What you have now

You have completed a guided SASA shielding analysis and can now answer:

- Did the polymer reduce the protein's solvent-accessible surface area?
- Was the reduction specific to the polymer-aware context?
- How much area does the polymer cover, and does it differ between conditions?

## Next steps

- Use {doc}`../how_to/analysis_sasa_quickstart` for active-site targets,
  per-residue results and the settings.
- Use {doc}`../explanation/analysis_sasa_verification` for how the values were
  checked.
- Use {doc}`../explanation/analysis_api` to run your own functions on every
  replicate.
