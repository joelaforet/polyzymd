# Catalytic Triad Analysis: Quick Start

Analyze catalytic triad geometry and integrity in under 5 minutes.

```{note}
This guide focuses on running the analysis quickly.
For interpretation and statistical guidance, see
{doc}`../explanation/analysis_triad_best_practices`.
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

## TL;DR

```bash
# Run catalytic triad comparison (all conditions in comparison.yaml)
polyzymd compare run catalytic_triad -f comparison.yaml --eq-time 100ns

# Equivalent full plugin name
polyzymd compare run catalytic_triad -f comparison.yaml --eq-time 100ns

# Run all enabled analyses in one workflow
polyzymd compare run-all -f comparison.yaml --eq-time 100ns

# Force recompute and machine-readable output
polyzymd compare run catalytic_triad -f comparison.yaml --eq-time 100ns --recompute --format json
```

## Prerequisites

Before running catalytic triad analysis, you need:

1. Completed production simulation output
2. A `comparison.yaml` file with a `plugins.catalytic_triad` section
3. One or more replicates per condition; use at least two replicates for SEM
   and robust inferential comparisons

## What Catalytic Triad Analysis Provides

The catalytic triad plugin computes:

| Feature | Description |
|---------|-------------|
| Per-pair distance | Mean distance for each configured pair |
| Per-pair contact fraction | Fraction below the configured threshold |
| Simultaneous contact fraction | Fraction where all pairs are in contact together |
| SEM | Autocorrelation-aware uncertainty estimates |
| Replicate aggregation | Condition-level means and SEM across replicates |
| Condition comparison | Ranking, pairwise tests, and effect sizes |

```{important}
The simultaneous contact fraction is the primary metric.
It estimates the fraction of analyzed frames where the full triad geometry is
simultaneously compatible with your contact threshold definition.
```

## Basic Usage

`````{tab-set}

````{tab-item} YAML (recommended)
Define triad settings in `comparison.yaml`:

```yaml
name: "polymer_triad_study"
control: "No Polymer"

conditions:
  - label: "No Polymer"
    config: "../projects/no_polymer/config.yaml"
    replicates: [1, 2, 3]

  - label: "With Polymer"
    config: "../projects/with_polymer/config.yaml"
    replicates: [1, 2, 3]

defaults:
  equilibration_time: "100ns"

plugins:
  catalytic_triad:
    name: "LipA_catalytic_triad"
    description: "Ser-His-Asp catalytic triad"
    threshold: 3.5
    pairs:
      - label: "Asp133-His156"
        selection_a: "midpoint(protein and resid 133 and name OD1 OD2)"
        selection_b: "protein and resid 156 and name ND1"
      - label: "His156-Ser77"
        selection_a: "protein and resid 156 and name NE2"
        selection_b: "protein and resid 77 and name OG"
```

Then run:

```bash
polyzymd compare run catalytic_triad -f comparison.yaml --eq-time 100ns
polyzymd compare run catalytic_triad -f comparison.yaml --eq-time 100ns --recompute
```
````

````{tab-item} CLI
Run the configured comparison:

```bash
polyzymd compare run catalytic_triad -f comparison.yaml --eq-time 100ns
```

Expected output pattern:

```text
Comparison: polymer_triad_study
Type: catalytic_triad
Conditions: 2
Equilibration: 100ns

Comparison result: comparison/catalytic_triad/result.json

Condition         Simultaneous Contact   SEM
------------------------------------------------
No Polymer        49.9%                  27.3%
With Polymer      87.3%                   2.2%
```

Replicates and conditions always come from `comparison.yaml`.
````

````{tab-item} Python
Use the API for scripted workflows:

```python
from polyzymd.analyses.discovery import get_analysis
from polyzymd.analyses.orchestrator import run_comparison
from polyzymd.config.comparison import ComparisonConfig

config = ComparisonConfig.from_yaml("comparison.yaml")
analysis = get_analysis("catalytic_triad")()

pipeline_result = run_comparison(
    analysis,
    config,
    recompute=False,
    equilibration="100ns",
)

result = pipeline_result["comparison"]
print(f"Ranking (best triad first): {result.ranking}")

for condition in result.conditions:
    contact_pct = condition.mean_simultaneous_contact * 100
    sem_pct = condition.sem_simultaneous_contact * 100
    print(f"{condition.label}: {contact_pct:.1f} ± {sem_pct:.1f}%")
```
````

`````

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

## Comparing Conditions

Use one `comparison.yaml` with all conditions, then run one command:

```bash
polyzymd compare run catalytic_triad -f comparison.yaml --eq-time 100ns
```

The comparison includes:
- Ranking by simultaneous contact fraction (higher is better)
- Pairwise tests with p-values and effect sizes
- Percent change relative to control

For broader comparison workflow guidance, see
{doc}`analysis_compare_conditions`.

## Reference and Troubleshooting

For full field tables, output file structures, JSON schemas, plot descriptions,
CLI option lookup, and troubleshooting fixes, see
{doc}`../reference/analysis_triad_reference`.

For statistical interpretation and best practices, see
{doc}`../explanation/analysis_triad_best_practices`.

## Next Steps

- Compare across conditions: {doc}`analysis_compare_conditions`
- Understand triad statistics: {doc}`../explanation/analysis_triad_best_practices`
- Run RMSD as a complementary global metric: {doc}`analysis_rmsd_quickstart`
- Run RMSF for per-residue flexibility: {doc}`analysis_rmsf_quickstart`
