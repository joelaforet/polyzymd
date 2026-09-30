# Which Analysis Should I Run?

You already have trajectories. This page helps you choose the analysis plugins that match the question you want to answer.

:::{admonition} Environment Setup
:class: tip

All analysis commands below assume you have activated the PolyzyMD analysis
pixi environment:

```bash
pixi shell -e analysis
```

Alternatively, prefix each command with `pixi run -e analysis`.
:::

## Quick Recommendation

If you are new to PolyzyMD or doing routine characterization, start with:

- `rmsd` — overall structural stability over time
- `rmsf` — per-residue flexibility
- `contacts` — polymer-protein interactions (if polymer is present)
- `secondary_structure` — fold integrity

This set gives a useful first pass before you move to more specialized analyses.

## Question-to-plugin chooser

| Question | Recommended plugins |
|----------|---------------------|
| Is my protein stable? | `rmsd`, `rmsf`, `secondary_structure` |
| Does the polymer interact with the protein? | `contacts`, `hydrogen_bonds` |
| Which residues interact with the polymer? | `contacts`, `hydrogen_bonds` |
| Is the active site accessible? | `sasa`, `catalytic_triad` |
| Does the polymer shield the protein surface? | `sasa` |
| Are catalytic residues properly positioned? | `catalytic_triad`, `distances` |
| Is a specific atom-pair distance maintained? | `distances` |
| How compact is the protein? | `rg` |

## Analysis quick reference

| Analysis | Runs through | Cost | Input needed |
|----------|--------------|------|--------------|
| `rmsd` | `polyzymd analyze rmsd` | Low | Atom selection |
| `rg` | `polyzymd analyze rg` | Low | Atom selection |
| `rmsf` | `polyzymd analyze rmsf` | Low | Atom selection |
| `secondary_structure` | `polyzymd analyze secondary_structure` | Low | (uses protein by default) |
| `contacts` | `polyzymd analyze contacts` | High (SASA with and without the polymer on every frame); lower with `--set method=distance` | Polymer + protein selections |
| `distances` | `polyzymd analyze distances` | Low | Atom pairs |
| `catalytic_triad` | `polyzymd analyze catalytic_triad` | Low | Residue pairs + threshold |
| `sasa` | `polyzymd analyze sasa` | High | Target + context selections |
| `hydrogen_bonds` | `polyzymd compare run hydrogen_bonds` (or `polyzymd analyze hydrogen_bonds -f comparison.yaml`) | High | Groups + summaries |

## Run several analyses

`polyzymd analyze` runs one analysis over every condition given with `-c`, so
run it once per analysis:

```bash
polyzymd analyze rmsf -c A/config.yaml -c B/config.yaml --eq 10ns
polyzymd analyze contacts -c A/config.yaml -c B/config.yaml --eq 10ns
```

`polyzymd compare run-all` runs every plugin enabled under `plugins:` in
`comparison.yaml`, such as `hydrogen_bonds`.

## See Also

- {doc}`analysis_compare_conditions` — Setting up comparison.yaml
- {doc}`../explanation/analysis_concepts` — How the analysis pipeline works
- {doc}`../reference/analysis_comparison_reference` — Full plugin listing
