# Catalytic triad: a routine on the analysis API

Measure whether the hydrogen bonds that make a catalytic triad work are
formed, pair by pair and all at once, on every production frame of every
replicate, and compare conditions with the replicate as the sampling unit.
PolyzyMD has no separate triad analysis: the triad is a few lines on the
analysis API, which shows how to combine the shipped functions for a question
of your own.

```{note}
For interpretation, see {doc}`../explanation/analysis_triad_best_practices`.
For what each shipped function measures, see
{doc}`../reference/analysis_functions`.
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

## The hydrogen bonds of the triad

In a serine hydrolase the serine hydroxyl donates to the histidine NE2, and
the histidine ND1-H donates to the aspartate carboxylate. The numbering below
is that of *B. subtilis* lipase A in a PolyzyMD topology, Ser76, His155 and
Asp132; use your enzyme's residues.

`hbond_count(group_a, group_b)` counts, at one frame, the hydrogen bonds
between two groups with MDAnalysis `HydrogenBondAnalysis`, the donor within
3.5 Å of the acceptor and a donor-hydrogen-acceptor angle of at least 150° by
default. Give each group the atoms of one side of the bond, the donor with its
hydrogen:

```python
import polyzymd as pz
from polyzymd.analyses.functions import hbond_count

study = pz.Study.from_configs(
    {"No polymer": "noPoly/config.yaml", "SBMA 50%": "SBMA50/config.yaml"},
    equilibration="200ns",
)
ser_his = study.timeseries(
    hbond_count,
    pz.select("protein and resid 76 and name OG HG"),
    pz.select("protein and resid 155 and name NE2"),
    unit=None,
    name="ser_his_hbonds",
)
his_asp = study.timeseries(
    hbond_count,
    pz.select("protein and resid 155 and name ND1 HD1"),
    pz.select("protein and resid 132 and name OD1 OD2"),
    unit=None,
    name="his_asp_hbonds",
)
both = ser_his.transform(
    lambda a, b: (a > 0) & (b > 0), his_asp, unit=None, name="triad_hbonds", bounds=(0.0, 1.0)
)
print(both.reduce("mean").compare(control="No polymer").to_agent_text())
```

`both` is 1 on the frames where both hydrogen bonds are formed, so its mean
per replicate is the fraction of frames with the whole triad hydrogen-bonded.
`ser_his.transform(lambda a: a > 0, ...)` gives one bond's fraction the same
way. The donors and acceptors are chosen by the element and valency rule of
{doc}`hydrogen_bonds`, and every result is stored under `polyzymd_results/`
with its record, so a second run reads it back. Pass `d_a_cutoff=3.0` or
`d_h_a_angle_cutoff=...` to `study.timeseries` to change the geometry. The
histidine tautomer matters: with the proton on NE2 instead of ND1, the
serine-histidine bond cannot form as written, and the counts are 0 rather than
an error.

On a 363 K lipase A replicate with a 50:50 SBMA-EGMA polymer, 716 production
frames took 10 s: the serine-histidine bond was formed on 80 percent of the
frames, the histidine-aspartate bond on 56 percent, and both on 44 percent.

## The triad distances from the command line

Heavy-atom distances need no hydrogens and run through `polyzymd analyze
distances`. Write the pairs to `triad.yaml`:

```yaml
- label: "Ser76-His155"
  selection_a: "protein and resid 76 and name OG"
  selection_b: "protein and resid 155 and name NE2"
- label: "His155-Asp132"
  selection_a: "protein and resid 155 and name ND1"
  selection_b: "midpoint(protein and resid 132 and name OD1 OD2)"
```

```bash
polyzymd analyze distances -c noPoly/config.yaml -c SBMA50/config.yaml \
  --label "No polymer" --label "SBMA 50%" --eq 200ns --set pairs=triad.yaml
```

For each pair this reports its mean distance and the fraction of frames below
the pair's threshold, 3.5 Å by default; see {doc}`analysis_distances_quickstart`.
To have every pair below its threshold at once, combine the distance series
with `all_below`:

```python
from polyzymd.analyses.functions import all_below, pair_distance

ser_his_d = study.timeseries(pair_distance, pz.select("protein and resid 76 and name OG"),
                             pz.select("protein and resid 155 and name NE2"), unit="A")
his_asp_d = study.timeseries(pair_distance, pz.select("protein and resid 155 and name ND1"),
                             pz.select("protein and resid 132 and name OD2"), unit="A")
close = ser_his_d.transform(all_below, his_asp_d, thresholds=[3.5, 3.5], unit=None)
print(close.reduce("mean").compare(control="No polymer").to_agent_text())
```

## Selection rules that prevent most errors

Keep catalytic residue selections inside the protein, since polymer and
substrate residues can share residue numbers:

```yaml
# Avoid: can match other chains
selection_a: "resid 76 and name OG"

# Prefer: constrained to the protein
selection_a: "protein and resid 76 and name OG"
```

For `pair_distance`, `midpoint(...)` and `com(...)` place one point at the
geometric center or center of mass of several atoms, which suits a
carboxylate acceptor.

## Next steps

- **Interpret triad results**: {doc}`../explanation/analysis_triad_best_practices`
- **Count hydrogen bonds between groups**: {doc}`hydrogen_bonds`
- **Measure other atom pairs**: {doc}`analysis_distances_quickstart`
- **Write a measurement of your own**: {doc}`../explanation/analysis_api`
