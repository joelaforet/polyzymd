# Measure polymer shielding with SASA

In this tutorial you measure how much of the surface of Trp-cage the SBMA
chains cover. You compute the solvent-accessible surface area (SASA) of the
protein twice: once with only the protein in the calculation, and once with
the polymer too. Then you compute the covered area in Python and compare it
between the conditions.

You learn these steps:

1. Give an analysis settings in `project.yaml`.
2. Report two results of SASA with `--run`.
3. Compute a new quantity from two time series with `Timeseries.transform`.
4. Find the figures to check first.

## Before you start

Do {doc}`analysis_complete_workflow` first. This tutorial continues in its
project folder, `~/pz_quickstart`, with the conditions `Water` and `SBMA`.

:::{admonition} Environment Setup
:class: tip

Run every command of this tutorial in the `build` environment. From the
repository root, activate it once:

```bash
pixi shell -e build
```
:::

## What SASA tells you

SASA is the area of the surface of a molecule that a probe the size of a
water molecule can reach. Each SASA calculation has two selections:

| Selection | Role |
|---|---|
| target | The atoms whose SASA is reported, here the protein |
| context | The atoms in the calculation. They can cover the surface of the target |

If the protein has a smaller SASA with the polymer in the context, the
polymer covers part of the protein surface.

## Step 1: List SASA with two contexts

Open `project.yaml` and add `sasa` under `analyses:`:

```yaml
analyses:
  rg: {}
  rmsf: {}
  contacts: {}
  hydrogen_bonds: {}
  sasa:
    contexts:
      isolated: protein
      with_polymer: protein or resname SBM
```

`contexts` names two calculations. `isolated` has only the protein.
`with_polymer` adds the monomers, whose residue name is `SBM`. In `Water`
there is no `SBM`, so both contexts are the same there. Commit the change:

```bash
git add -A
git commit -m "Add the sasa analysis"
```

## Step 2: Measure the protein on its own

```bash
polyzymd analyze sasa --study trpcage --run isolated
```

The output is:

```
log: /home/me/pz_quickstart/trpcage/logs/polyzymd-analyze-20261006-210823-pid23165.log
# polyzymd analyze sasa  metric mean_sasa  unit A^2  run isolated  eq 0ns  conditions 2  replicates 3,3  protocol sasa/2
Water  n 3  mean 1927  sem 10.85  ci95 1881 to 1974  values 1937, 1906, 1939  replicates 1, 2, 3  g 1, 1, 1  n_eff 4, 4, 4  eq_detected 0.001 ns
SBMA  n 3  mean 1924  sem 22.86  ci95 1826 to 2022  values 1912, 1968, 1892  replicates 1, 2, 3  g 1, 1, 1  n_eff 4, 4, 4  eq_detected 0.001 ns
Water vs SBMA  delta -3.236  ci95 -86.08 to 79.6  p 0.9067  p_adj 0.9067  test welch_t  correction BH  family 1  d -0.1044  not_significant
warning: condition Water, condition SBMA: replicates 1, 2, 3 have fewer than 20 effective samples, so the start of an equilibrated region cannot be detected reliably; values and statistics are unaffected
verdict: no significant difference in mean_sasa between Water and SBMA (delta -3.236 A^2, 95% CI -86.08 to 79.6, p_adj 0.9067, p 0.9067, n 3 vs 3)
```

Your values differ a little, because each run adds up the forces in a
different order.

These runs are three replicates of a few picoseconds. The chains start within
0.5 nm of the protein, so they cover it from the start. The numbers show the
method, not a property of SBMA.

The protein alone has about 1930 Å² of SASA in both conditions. At this
sample size no difference shows between them. That does not prove that the
polymer leaves the surface of the protein unchanged.

## Step 3: Measure the protein with the polymer

```bash
polyzymd analyze sasa --study trpcage --run with_polymer
```

The output is:

```
log: /home/me/pz_quickstart/trpcage/logs/polyzymd-analyze-20261006-210828-pid23332.log
# polyzymd analyze sasa  metric mean_sasa  unit A^2  run with_polymer  eq 0ns  conditions 2  replicates 3,3  protocol sasa/2
Water  n 3  mean 1927  sem 10.85  ci95 1881 to 1974  values 1937, 1906, 1939  replicates 1, 2, 3  g 1, 1, 1  n_eff 4, 4, 4  eq_detected 0.001 ns
SBMA  n 3  mean 1610  sem 43.08  ci95 1425 to 1795  values 1526, 1637, 1668  replicates 1, 2, 3  g 1, 1, 1  n_eff 4, 4, 4  eq_detected 0.001 ns
Water vs SBMA  delta -317.1  ci95 -489.1 to -145.1  p 0.01371  p_adj 0.01371  test welch_t  correction BH  family 1  d -5.829  significant
warning: condition Water, condition SBMA: replicates 1, 2, 3 have fewer than 20 effective samples, so the start of an equilibrated region cannot be detected reliably; values and statistics are unaffected
verdict: SBMA smaller mean_sasa than Water (delta -317.1 A^2, 95% CI -489.1 to -145.1, p_adj 0.01371, p 0.01371, n 3 vs 3)
```

With the polymer in the context, the protein in `SBMA` has about 320 Å² less
accessible surface. `Water` gives the same value as in step 2.

The two steps together show the pattern that shielding gives:

1. `isolated` shows no difference between the conditions.
2. `with_polymer` is smaller in the polymer condition.

## Step 4: Compute the covered area

The covered area is `isolated` minus `with_polymer`, frame by frame.
`Timeseries.transform` computes it from two time series. Save this script as
`covered.py` in the project folder:

```python
import numpy as np
import polyzymd as pz
from polyzymd.analyses.functions import sasa

study = pz.Study("trpcage")
folder = study.root / "results" / "sasa_covered"
protein = pz.select("protein")
with_polymer = pz.select("protein or resname SBM")
isolated = study.timeseries(
    sasa, protein, protein, unit="A^2", name="isolated", output_dir=folder
)
sasa_with_polymer = study.timeseries(
    sasa, protein, with_polymer, unit="A^2", name="with_polymer", output_dir=folder
)
covered = isolated.transform(np.subtract, sasa_with_polymer, unit="A^2", name="covered")
print(covered.reduce("mean").compare(control="Water").to_agent_text())
```

`pz.Study("trpcage")` reads the conditions and the window from
`trpcage/study.yaml`. `study.timeseries` runs the shipped `sasa` function on
every frame of every replicate. `reduce("mean")` gives one value per
replicate, and `compare` tests each condition against `Water`. Run it:

```bash
python covered.py
```

It prints:

```
# polyzymd analyze covered  metric mean_covered  unit A^2  eq 0ns  conditions 2  replicates 3,3  protocol covered/2
Water  n 3  mean 0  sem 0  ci95 na  values 0, 0, 0  replicates 1, 2, 3  g 1, 1, 1  n_eff 4, 4, 4  eq_detected 0.001 ns
SBMA  n 3  mean 313.9  sem 47.48  ci95 109.6 to 518.2  values 386, 331.3, 224.3  replicates 1, 2, 3  g 1, 1, 1  n_eff 4, 4, 4  eq_detected 0.001 ns
Water vs SBMA  delta +313.9  ci95 109.6 to 518.2  p 0.02213  p_adj 0.02213  test welch_t  correction BH  family 1  d 5.397  significant
warning: condition Water, condition SBMA: replicates 1, 2, 3 have fewer than 20 effective samples, so the start of an equilibrated region cannot be detected reliably; values and statistics are unaffected
warning: condition Water has the same mean_covered in every replicate, so its interval is not estimable
verdict: SBMA larger mean_covered than Water (delta +313.9 A^2, 95% CI 109.6 to 518.2, p_adj 0.02213, p 0.02213, n 3 vs 3; little power: Water has the same value in every replicate)
```

Python also prints warnings of MDAnalysis and one that says the bytecode of
`np.subtract` is hashed, because it has no source file. They do not change
the values.

`SBMA` covers about 310 Å² of the protein surface. Step 3 compared the SASA
of `SBMA` with that of `Water`. The covered area subtracts the SASA of the
same `SBMA` replicate, so the two numbers differ by a few Å². The values are
stored in `trpcage/results/sasa_covered/polyzymd_results/`. A second run
reads the two SASA series back instead of measuring them again.

## Step 5: Check the figures

`polyzymd analyze sasa` wrote its figures to `trpcage/results/sasa/figures/sasa/`:

| Figure | What it shows |
|---|---|
| `sasa_comparison_with_polymer.png` | The mean SASA of each condition, with its interval and every replicate value |
| `sasa_timeseries_with_polymer.png` | The SASA of each replicate against time |
| `sasa_distribution_with_polymer.png` | The distribution of the SASA values of each condition |

Each figure has an `isolated` version too. In a real study, check in the time
series that the SASA settles after the equilibration window.

## What you did

You measured the SASA of a protein with and without a polymer in the
calculation, and computed the covered area from the two series in Python. To
see which residues the polymer covers, report `--run with_polymer_residues`.
For the settings, such as an active-site target, see
{doc}`../how_to/analysis_sasa_quickstart`. To run your own function on every
replicate, see {doc}`../how_to/study_api`.
