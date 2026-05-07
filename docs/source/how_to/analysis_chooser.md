# Which Analysis Should I Run?

You already have trajectories. This page helps you choose the analysis plugins that match the question you want to answer.

:::{admonition} Environment Setup
:class: tip

All commands below assume you have activated the PolyzyMD pixi environment:

```bash
pixi shell -e build
```

Alternatively, prefix each command with `pixi run -e build`.
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
| Does the polymer bridge different protein regions? | `polymer_bridging` (experimental) |

## Plugin quick reference

| Plugin | Stable? | Cost | Input needed |
|--------|---------|------|--------------|
| `rmsd` | ✓ | Low | Atom selection |
| `rg` | ✓ | Low | Atom selection |
| `rmsf` | ✓ | Low | Atom selection |
| `secondary_structure` | ✓ | Low | (uses protein by default) |
| `contacts` | ✓ | Medium | Polymer + protein selections |
| `distances` | ✓ | Low | Atom pairs |
| `catalytic_triad` | ✓ | Low | Residue pairs + threshold |
| `sasa` | ✓ | High | Target + context selections |
| `hydrogen_bonds` | ✓ | High | Groups + summaries |
| `polymer_bridging` | Experimental | Medium | Polymer + protein selections |

## Run multiple plugins in one pass

You can enable multiple plugins in the same `comparison.yaml` file:

```yaml
plugins:
  rmsd:
    runs:
      - label: "Protein Backbone"
        selection: "protein and name CA"
  rmsf:
    selection: "protein and name CA"
  contacts: {}
  secondary_structure: {}
```

Then run all configured plugins at once:

```bash
polyzymd compare run-all
```

## See Also

- {doc}`analysis_compare_conditions` — Setting up comparison.yaml
- {doc}`../explanation/analysis_concepts` — How the analysis pipeline works
- {doc}`../reference/analysis_comparison_reference` — Full plugin listing
