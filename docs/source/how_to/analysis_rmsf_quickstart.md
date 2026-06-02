# RMSF Analysis: Quick Start

Run RMSF (Root Mean Square Fluctuation) in minutes and compare flexibility
across conditions.

```{note}
Need field-by-field settings, output schemas, plot options, or troubleshooting?
Use {doc}`../reference/analysis_rmsf_reference`.
```

:::{admonition} Environment Setup
:class: tip

All commands below assume you have activated the PolyzyMD pixi environment:

```bash
pixi shell -e build
```

Alternatively, prefix each command with `pixi run -e build`.
:::

## TL;DR

```bash
# Run RMSF for conditions in comparison.yaml
polyzymd compare run rmsf -f comparison.yaml --eq-time 10ns

# Run all enabled analyses
polyzymd compare run-all -f comparison.yaml --eq-time 10ns

# Force recompute and write machine-readable output
polyzymd compare run rmsf -f comparison.yaml --eq-time 10ns --recompute --format json
```

## Prerequisites

Before you run RMSF, confirm you have:

1. Completed production trajectories for at least one replicate
2. A `comparison.yaml` with one or more conditions
3. The original simulation `config.yaml` paths still valid

Quick check:

```bash
ls $(polyzymd info -c config.yaml --scratch-dir)/production_*/
```

## What RMSF gives you

RMSF reports per-residue flexibility plus condition-level summary statistics.

| Output | Use |
|--------|-----|
| Mean RMSF ± SEM | Compare overall flexibility between conditions |
| Per-residue RMSF profile | Find flexible loops and rigid core regions |
| Pairwise statistics | Quantify condition differences |
| Ranking | Order conditions by stability (lower RMSF first) |

## Set up `comparison.yaml`

Use `plugins.rmsf` for RMSF-specific controls.

```yaml
# comparison.yaml
name: "rmsf_quickstart"
control: "No Polymer"

conditions:
  - label: "No Polymer"
    config: "configs/no_polymer.yaml"
    replicates: [1, 2, 3]
  - label: "With Polymer"
    config: "configs/with_polymer.yaml"
    replicates: [1, 2, 3]

defaults:
  equilibration_time: "10ns"

plugins:
  rmsf:
    selection: "protein and name CA"
    reference_mode: "centroid"
```

## Run RMSF

```bash
# RMSF only
polyzymd compare run rmsf -f comparison.yaml --eq-time 10ns

# RMSF + all other enabled plugins
polyzymd compare run-all -f comparison.yaml --eq-time 10ns

# Generate plots in the same run
polyzymd compare run-all -f comparison.yaml --eq-time 10ns --plot
```

## Compare conditions

By default, `compare run rmsf` produces:

- Condition means with SEM
- Pairwise tests and effect sizes
- Ranking (lowest RMSF = most stable)

Example:

```bash
polyzymd compare run rmsf -f comparison.yaml --format table
```

For broader multi-plugin workflows, see {doc}`analysis_compare_conditions`.

## Comparing Two Conditions

Use this minimal two-condition template:

```yaml
name: "polymer_stabilization"
control: "No Polymer"

conditions:
  - label: "No Polymer"
    config: "../no_polymer/config.yaml"
    replicates: [1, 2, 3]
  - label: "With Polymer"
    config: "../with_polymer/config.yaml"
    replicates: [1, 2, 3]

plugins:
  rmsf:
    selection: "protein and name CA"
    reference_mode: "centroid"
```

Then run:

```bash
polyzymd compare run rmsf -f comparison.yaml --eq-time 10ns
```

## Common task patterns

### Change equilibration cutoff

```bash
polyzymd compare run rmsf -f comparison.yaml --eq-time 100ns
```

### Restrict analysis to a region

```yaml
plugins:
  rmsf:
    selection: "protein and name CA and resid 50-150"
```

### Trim flexible termini

```yaml
plugins:
  rmsf:
    selection: "protein and name CA and resid 5:175"
```

### Use an external reference structure

```yaml
plugins:
  rmsf:
    selection: "protein and name CA"
    reference_mode: "external"
    reference_file: "/path/to/crystal_structure.pdb"
```

### Re-run from scratch

```bash
polyzymd compare run rmsf -f comparison.yaml --eq-time 10ns --recompute
```

### Save results for scripting

```bash
polyzymd compare run rmsf -f comparison.yaml --format json -o rmsf_result.json
```

## Python API (programmatic run)

```python
from polyzymd.analyses.discovery import get_analysis
from polyzymd.analyses.orchestrator import run_comparison
from polyzymd.config.comparison import ComparisonConfig

config = ComparisonConfig.from_yaml("comparison.yaml")
analysis = get_analysis("rmsf")()

pipeline_result = run_comparison(
    analysis,
    config,
    equilibration="10ns",
)

result = pipeline_result["comparison"]
print(result.ranking)
print(pipeline_result["comparison_path"])
```

## Where to find details

- Full settings tables: {doc}`../reference/analysis_rmsf_reference`
- Output file tree and JSON schemas: {doc}`../reference/analysis_rmsf_reference`
- Plot types and plot settings: {doc}`../reference/analysis_rmsf_reference`
- CLI options and troubleshooting: {doc}`../reference/analysis_rmsf_reference`
- Statistical interpretation: {doc}`../explanation/analysis_rmsf_best_practices`
- Reference mode guidance: {doc}`../explanation/analysis_reference_selection`

## Next steps

- Compare multiple plugins together: {doc}`analysis_compare_conditions`
- Run RMSD alongside RMSF: {doc}`analysis_rmsd_quickstart`
- Add distance checks: {doc}`analysis_distances_quickstart`
