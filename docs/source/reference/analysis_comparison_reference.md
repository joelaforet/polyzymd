# Comparison and Plotting Reference

Use this page when you need lookup information for `polyzymd compare`,
`comparison.yaml`, output files, or the plotting workflow.

## Comparison Project Layout

`polyzymd compare init -n my_study` creates a workspace like this:

```text
my_study/
├── comparison.yaml
├── results/
├── figures/
└── structures/
```

<!-- IMAGE OPPORTUNITY: Add an annotated project-tree graphic showing where
comparison configuration, results JSON, and figures are written. -->

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

analysis_settings:
  rmsf:
    selection: "protein and name CA"

comparison_settings:
  rmsf: {}
```

## Stable Analysis Keys

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

## Path Rules

- Relative paths in `config:` are resolved relative to `comparison.yaml`
- Absolute paths are used as-is
- `replicates` must be an explicit list such as `[1, 2, 3]`

## Commands

| Command | Purpose |
|---------|---------|
| `polyzymd compare init -n NAME` | Create a comparison workspace |
| `polyzymd compare validate` | Check `comparison.yaml` before running |
| `polyzymd compare run TYPE` | Run a single comparison through the registry |
| `polyzymd compare run-all` | Run every enabled comparison in one pass |
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
polyzymd compare exposure
polyzymd compare binding-free-energy
polyzymd compare polymer-affinity
```

These remain callable, but PolyzyMD labels them as experimental in CLI output,
docs, and generated figures.

## Output Locations

- comparison JSON files are written to `results/`
- figures are written under the configured `plot_settings.output_dir`
- default project scaffolds create a `figures/` directory next to
  `comparison.yaml`

Typical result filenames:

```text
results/rmsf_comparison_<study>.json
results/contacts_comparison_<study>.json
results/distances_comparison_<study>.json
results/triad_comparison_<study>.json
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
