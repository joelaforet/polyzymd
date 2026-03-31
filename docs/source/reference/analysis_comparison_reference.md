# Comparison and Plotting Reference

Use this page when you need quick lookup information for `polyzymd compare`,
`comparison.yaml`, output paths, or plotting behavior.

## Comparison Project Layout

`polyzymd compare init -n my_study` creates a workspace like this:

```text
my_study/
├── comparison.yaml
├── comparison/
├── figures/
└── structures/
```

## Core `comparison.yaml` Fields

```yaml
name: "polymer_stability_study"
description: "Optional human-readable summary"
control: "No Polymer"  # optional

conditions:
  - label: "No Polymer"
    config: "../no_polymer/config.yaml"
    replicates: [1, 2, 3]

defaults:
  equilibration_time: "10ns"

plugins:
  rmsf:
    selection: "protein and name CA"
```

## Stable Plugin Keys

Stable release workflows for `v1.2.0`:

- `rmsf`
- `contacts`
- `distances`
- `catalytic_triad`
- `secondary_structure`

Experimental but still available:

- binding preference through `contacts`
- `exposure`
- `binding_free_energy`
- `polymer_affinity`
- `polymer_bridging` (alias: `bridging`)

## Path Rules

- relative paths in `config:` are resolved relative to `comparison.yaml`
- absolute paths are used as-is
- `replicates` must be an explicit list such as `[1, 2, 3]`

## Commands

| Command | Purpose |
|---------|---------|
| `polyzymd compare init -n NAME` | Create a comparison workspace |
| `polyzymd compare validate` | Check `comparison.yaml` before running |
| `polyzymd compare run TYPE` | Run one analysis plugin |
| `polyzymd compare run --list` | List available comparison types and aliases |
| `polyzymd compare run-all` | Run every enabled plugin in one pass |
| `polyzymd compare plot-all` | Generate configured figures |
| `polyzymd compare plot-all --list-available` | List available plots and experimental labels |

## Common Stable Commands

```bash
polyzymd compare run rmsf
polyzymd compare run contacts
polyzymd compare run distances
polyzymd compare run triad
polyzymd compare run-all
polyzymd compare plot-all
```

## Experimental Commands

```bash
polyzymd compare run exposure
polyzymd compare run binding_free_energy
polyzymd compare run polymer_affinity
polyzymd compare run polymer_bridging   # alias: polyzymd compare run bridging
```

These remain callable, but PolyzyMD labels them as experimental in CLI output,
docs, and generated figures.

## Output Locations

- comparison JSON files are written to `comparison/<analysis>/result.json`
- figures are written under the configured `plot_settings.output_dir`
- default project scaffolds create a `figures/` directory next to
  `comparison.yaml`

Typical comparison cache paths:

```text
comparison/rmsf/result.json
comparison/contacts/result.json
comparison/distances/result.json
comparison/catalytic_triad/result.json
```

## Plotting Smoke Test

For a final smoke test after comparisons finish:

```bash
polyzymd compare plot-all --list-available
polyzymd compare plot-all
```

## Statistical Terms

- `p-value`: significance of the observed difference under the null hypothesis
- `Cohen's d`: effect size magnitude
- `ANOVA`: omnibus test across multiple conditions
- `SEM`: standard error of the mean across replicates

For interpretation guidance rather than lookup, see:

- [Statistical Best Practices for Analysis](../tutorials/analysis_statistics_best_practices.md)
- [How to Compare Simulation Conditions](../tutorials/analysis_compare_conditions.md)
