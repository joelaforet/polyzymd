# Run polymer Rg analysis

Measure the radius of gyration (Rg) of each polymer chain on each production
frame of each replicate, and compare how compact the chains are between
conditions, with one value per replicate.

:::{admonition} Environment Setup
:class: tip

All analysis commands below assume you have activated the PolyzyMD analysis
pixi environment:

```bash
pixi shell -e analysis
```

Alternatively, prefix each command with `pixi run -e analysis`.
:::

## From the command line

```bash
polyzymd analyze polymer_rg --study my_study
polyzymd analyze polymer_rg -c SBMA25/config.yaml -c SBMA50/config.yaml --eq 200ns
```

The command does these steps:

1. It splits the polymer, `chainid C`, into chains by its bonds: a PolyzyMD
   build puts every chain in chain C, so the bonds are what tell them apart.
2. On each frame after the equilibration window, it makes each chain whole
   across the periodic box and measures its mass-weighted Rg.
3. `--run mean_rg`, the default, averages the chains on each frame and takes
   the mean over frames of each replicate. `--run chain_rg` reports each
   chain's mean over frames instead, labelled 1, 2, ... by its position in the
   topology.
4. It compares each condition with the control, with the replicate as the
   sampling unit.

A condition without polymer, such as a no-polymer control, is left out with a
warning, so put a condition with polymer first to compare against it. To
measure other chains, such as a surfactant, add `--set selection='resname SDS'`.

## From Python

```python
import polyzymd as pz
from polyzymd.analyses.functions import polymer_rg

study = pz.Study.from_configs(
    {"SBMA 25%": "SBMA25/config.yaml", "SBMA 50%": "SBMA50/config.yaml"},
    equilibration="200ns",
)
rg = study.timeseries(polymer_rg, pz.select("chainid C"), unit="A")
print(rg.reduce("mean").compare().to_agent_text())
```

## Before you interpret the numbers

- The topology must carry bonds for every polymer atom; PolyzyMD reads them
  from the run's OpenMM system or GROMACS run input. Without them the command
  stops with an error rather than reporting a number.
- `chain_rg` labels chains by position. Chain 2 of one condition is the same
  molecule as chain 2 of another only when the builds place the same chain
  there, so compare per-chain values across conditions with care; `mean_rg`
  does not depend on the order.
- The mean over chains counts every chain once, whatever its length.

See {doc}`../reference/analysis_functions` for the functions.
