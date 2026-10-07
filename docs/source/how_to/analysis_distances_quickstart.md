# Run distance analysis

Measure the distance between named atom pairs on each production frame of each
replicate. For each pair, PolyzyMD reports the mean distance and the fraction
of frames below a threshold. Then it compares the conditions, with one value
per replicate.

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
Use the analysis that fits the question:

- **Distances**: the distance between specific atoms, as a continuous value.
- **Contacts**: the fraction of frames in which the polymer buries or touches
  each protein residue. See {doc}`analysis_contacts_quickstart`.
- **Catalytic triad**: the hydrogen bonds of the triad on each frame, combined
  with these distances. See {doc}`analysis_triad_quickstart`.
```

## Define the pairs

Write the pairs to a YAML (or JSON) file, for example `pairs.yaml`:

```yaml
- label: "Ser76(OG)-Substrate(carbonyl C)"
  selection_a: "protein and resid 76 and name OG"
  selection_b: "resname RBY and name C13x"
  threshold: 10.0
  below_label: "Within 10 Angstrom"
- label: "Met77(N)-Substrate(carbonyl C)"
  selection_a: "protein and resid 77 and name N"
  selection_b: "resname RBY and name C13x"
```

Each pair needs `label`, `selection_a` and `selection_b`. Each pair label must
be unique. The optional keys are:

- `threshold`: the cutoff for the fraction of frames below it. A pair without
  one uses the analysis threshold. The default is 3.5 Å. Change it with
  `--set threshold=...`.
- `below_label`: the name of that fraction in the report.

## Run it

In a {term}`study`, give the study folder. The study names the conditions,
the control and the equilibration window. The results go to
`<study>/results/distances/`:

```bash
polyzymd analyze distances --study my_study --set pairs=pairs.yaml
```

For a quick look without a study, give the `config.yaml` of each condition
instead. The results then go to the current folder:

```bash
polyzymd analyze distances -c noPoly/config.yaml -c SBMA50/config.yaml \
  --label "No polymer" --label "SBMA 50%" --eq 200ns --set pairs=pairs.yaml
```

The first `-c` is the control. PolyzyMD measures each pair on each production
frame of each replicate. It reports two results for each pair:

| Result | Meaning |
|---|---|
| `<label>` | The mean distance |
| `<label> <below_label>`, or `<label> below <threshold> A` | The fraction of frames below the threshold |

The report shows the mean distance of the first pair. To report a different
result, use `--run`, for example
`--run "Ser76(OG)-Substrate(carbonyl C) Within 10 Angstrom"`. PolyzyMD stores
the measured distances, so a second `--run` does not read the trajectory again.
PolyzyMD compares the conditions by Welch's t test. It corrects the p values
over the compared conditions with the {term}`Benjamini-Hochberg` method.

## Write selections

A selection is an MDAnalysis selection string. PolyzyMD adds three forms:
`midpoint(...)`, `com(...)` and `pdbindex N`.

````{warning}
**Include the chain in each selection**

Residue numbers restart in each chain of a PolyzyMD system. A selection such as
`resid 141-148` can match atoms in more than one chain.

For protein residues, start the selection with `protein and`:

```yaml
# Incorrect
selection_a: "com(resid 141-148)"

# Correct
selection_a: "com(protein and resid 141-148)"
```
````

Common patterns:

```yaml
# Midpoint of Asp carboxylate oxygens
selection_a: "midpoint(protein and resid 132 and name OD1 OD2)"

# Center of mass of ligand
selection_b: "com(resname LIG)"

# Single atom
selection_a: "protein and resid 76 and name OG"

# The 2740th atom of the system, counted from 1
selection_a: "pdbindex 2740"
```

A plain selection must match exactly one atom. `pdbindex N` selects the N-th
atom of the system, as in restraints. It is the same as the MDAnalysis
selection `bynum N`.

## Keep PBC on

By default (`use_pbc: true`), PolyzyMD measures distances with the minimum
image convention and the box of each frame. This stops PolyzyMD from measuring
a pair the long way around the box. Keep it on, unless the trajectory is
already unwrapped. If a frame has no valid box, PolyzyMD measures it without
the minimum image and prints a warning. PolyzyMD does not align frames for
distances, because a distance does not change when the system moves or
rotates.

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
    unit="A",
)
print(d.reduce("mean").compare(control="No polymer").to_agent_text())
below = d.transform(lambda x, cutoff: x < cutoff, cutoff=10.0, unit=None)
print(below.reduce("fraction").compare(control="No polymer").to_agent_text())
```

`transform` computes the per-frame fraction from the stored distances. It does
not read the trajectory again. Give values such as the cutoff as keyword
arguments, as above. PolyzyMD then records them with the result. To measure
your own quantity, see {doc}`study_api`.

## Next steps

- **Catalytic triad analysis**: {doc}`analysis_triad_quickstart`
- **Understand statistics**: {doc}`../explanation/analysis_statistics_best_practices`
- **Contact analysis**: {doc}`analysis_contacts_quickstart`
