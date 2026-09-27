# Catalytic triad analysis: quick start

Measure the distances that define a catalytic triad on every production frame
of every replicate, and report how often each contact, and the whole triad, is
formed.

```{note}
This guide focuses on getting results quickly. For interpretation, see
{doc}`../explanation/analysis_triad_best_practices`. For what each shipped
function measures, see {doc}`../reference/analysis_functions`.
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

## Define the triad

Write the pairs to `triad.yaml`:

```yaml
- label: "Ser77-His156"
  selection_a: "protein and resid 76 and name OG"
  selection_b: "protein and resid 155 and name NE2"
- label: "His156-Asp133"
  selection_a: "protein and resid 155 and name ND1"
  selection_b: "protein and resid 132 and name OD2"
```

## Run it

```bash
polyzymd analyze catalytic_triad -c noPoly/config.yaml -c SBMA50/config.yaml \
  --label "No polymer" --label "SBMA 50%" --eq 200ns --set pairs=triad.yaml
```

For each replicate this reports:

| Result | Meaning |
|--------|---------|
| `simultaneous` | Fraction of frames in which every pair is below its threshold at once |
| `<label>` | Mean distance of that pair, in Å |
| `<label> below <threshold> A` | Fraction of frames in which that pair is below its threshold |

The threshold is 3.5 Å unless you set `--set threshold=...` or give a pair its
own `threshold`. The report shows `simultaneous`; pick another result with
`--run "<label>"`. The simultaneous fraction is computed from the stored pair
distances, so the trajectory is read once per pair. Conditions are compared by
Welch's t test with the Benjamini-Hochberg correction.

```{important}
The simultaneous contact fraction is the primary metric.
It estimates the fraction of analyzed frames where the full triad geometry is
simultaneously compatible with your contact threshold definition. It is a
fraction from 0 to 1.
```

## Selection Rules That Prevent Most Errors

Use chain-aware selections for catalytic residues.

```yaml
# Avoid: can match multiple chains
selection_a: "resid 77 and name OG"

# Prefer: constrained to the protein
selection_a: "protein and resid 77 and name OG"
```

Supported selection forms:

| Syntax | Description | Example |
|--------|-------------|---------|
| Standard | MDAnalysis selection | `protein and resid 77 and name OG` |
| `midpoint()` | Geometric midpoint | `midpoint(protein and resid 133 and name OD1 OD2)` |
| `com()` | Center of mass | `com(protein and resid 133 and name OD1 OD2)` |

```{tip}
`midpoint()` is often a good choice for Asp/Glu carboxylate acceptors.
```

## From Python

```python
import polyzymd as pz
from polyzymd.analyses.functions import all_below, pair_distance

study = pz.Study.from_configs(
    {"No polymer": "noPoly/config.yaml", "SBMA 50%": "SBMA50/config.yaml"},
    equilibration="200ns",
)
ser_his = study.timeseries(pair_distance, pz.select("protein and resid 76 and name OG"),
                           pz.select("protein and resid 155 and name NE2"), unit="Å")
his_asp = study.timeseries(pair_distance, pz.select("protein and resid 155 and name ND1"),
                           pz.select("protein and resid 132 and name OD2"), unit="Å")
both = ser_his.transform(all_below, his_asp, thresholds=[3.5, 3.5], unit=None)
print(both.reduce("fraction").compare(control="No polymer").to_agent_text())
```

## Next Steps

- **Interpret triad results**: {doc}`../explanation/analysis_triad_best_practices`
- **Measure other atom pairs**: {doc}`analysis_distances_quickstart`
- **Run RMSD as a complementary global metric**: {doc}`analysis_rmsd_quickstart`
