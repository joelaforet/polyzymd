# Distance analysis: quick start

Measure the distance between named atom pairs on every production frame of
every replicate, report each pair's mean distance and the fraction of frames
below a threshold, and compare conditions with the replicate as the sampling
unit.

```{note}
**Want to understand the statistics?** This guide focuses on getting results
quickly. For uncertainty and replicate-level comparison, see
{doc}`../explanation/analysis_statistics_best_practices`. For what each shipped
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

```{tip}
**When to use distances vs. contacts vs. triad:**
- **Distances**: Specific atom pairs with continuous distance values
- **Contacts**: All residue-residue contacts at an interface (binary count)
- **Triad**: a routine on the analysis API, {doc}`analysis_triad_quickstart`, that
  counts the triad's hydrogen bonds on every frame and combines them with these distances
```

## Define the pairs

Write the pairs to a YAML (or JSON) file, for example `pairs.yaml`:

```yaml
- label: "Ser77(OG)-Substrate(carbonyl C)"
  selection_a: "protein and resid 76 and name OG"
  selection_b: "resname RBY and name C13x"
  threshold: 10.0
  below_label: "Within 10 Angstrom"
- label: "Met78(N)-Substrate(carbonyl C)"
  selection_a: "protein and resid 77 and name N"
  selection_b: "resname RBY and name C13x"
```

Each pair needs `label`, `selection_a` and `selection_b`. `threshold` sets the
cutoff for that pair's fraction of frames below it; pairs without one use the
analysis threshold, 3.5 Å unless you set `--set threshold=...`. `below_label`
names that fraction in the report.

## Run it

```bash
polyzymd analyze distances -c noPoly/config.yaml -c SBMA50/config.yaml \
  --label "No polymer" --label "SBMA 50%" --eq 200ns --set pairs=pairs.yaml
```

The first `-c` is the control. Every pair is measured on every production frame
of every replicate, and two results are reported per pair: `<label>`, the mean
distance, and `<label> <below_label>` (or `<label> below <threshold> A`), the
fraction of frames below the threshold. The report shows the first pair's mean
distance; pick another result with `--run`, for example
`--run "Ser77(OG)-Substrate(carbonyl C) Within 10 Angstrom"`. The measured
distances are stored, so a second `--run` does not read the trajectory again.
Conditions are compared by Welch's t test with the Benjamini-Hochberg
correction across the conditions compared.

## Write Robust Selections

PolyzyMD supports standard MDAnalysis selections plus helper syntax like
`midpoint(...)`, `com(...)`, and `pdbindex N`.

```{warning}
**Chain-aware selections are required**

Residue numbers restart by chain in PolyzyMD systems. A selection like
`resid 141-148` can match multiple chains.

For protein residues, include `protein and ...`:

```yaml
# Incorrect
selection_a: "com(resid 141-148)"

# Correct
selection_a: "com(protein and resid 141-148)"
```
```

Common patterns:

```yaml
# Midpoint of Asp carboxylate oxygens
selection_a: "midpoint(protein and resid 133 and name OD1 OD2)"

# Center of mass of ligand
selection_b: "com(resname LIG)"

# Single atom
selection_a: "protein and resid 77 and name OG"

# Atom by PDB serial number
selection_a: "pdbindex 2740"
```

A plain selection must match exactly one atom.

## Keep PBC on

Distances use the minimum image convention by default (`use_pbc`), with the
box stored in each frame. Leave it on unless you know your trajectory is
already unwrapped, because it is what keeps a pair from being measured the
long way around the box. A frame with no valid box is measured without the
minimum image, and a warning says so. Distances are never aligned: a distance
does not change when the whole system is rotated or translated.

## From Python

```python
import polyzymd as pz
from polyzymd.analyses.functions import pair_distance

study = pz.Study.from_configs(
    {"No polymer": "noPoly/config.yaml", "SBMA 50%": "SBMA50/config.yaml"},
    equilibration="200ns",
)
d = study.timeseries(
    pair_distance,
    pz.select("protein and resid 76 and name OG"),
    pz.select("resname RBY and name C13x"),
    unit="Å",
)
print(d.reduce("mean").compare(control="No polymer").to_agent_text())
below = d.transform(lambda x, cutoff: x < cutoff, cutoff=10.0, unit=None)
print(below.reduce("fraction").compare(control="No polymer").to_agent_text())
```

`transform` builds the per-frame fraction from the stored distances without
reading the trajectory again. Pass values such as the cutoff as keyword
arguments, as above, so they are recorded with the result.

## Next Steps

- **Catalytic triad analysis**: {doc}`analysis_triad_quickstart`
- **Understand statistics**: {doc}`../explanation/analysis_statistics_best_practices`
- **Contact analysis**: {doc}`analysis_contacts_quickstart`
