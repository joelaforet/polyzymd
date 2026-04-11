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
    config: "no_polymer/config.yaml"
    replicates: [1, 2, 3]

defaults:
  equilibration_time: "10ns"

plugins:
  rmsf:
    selection: "protein and name CA"
```

## Per-Plugin Statistical Settings

Some plugins support per-plugin statistical settings configured under the
`plugins:` block in `comparison.yaml`. These control false discovery rate
correction, effect-size filtering, and output truncation for cross-condition
comparisons.

### Canonical YAML Example

```yaml
plugins:
  contacts:
    cutoff: 4.5
    fdr_alpha: 0.05
    min_effect_size: 0.5
    top_residues: 10

  binding_free_energy:
    units: "kcal/mol"
    fdr_alpha: 0.05

  polymer_affinity:
    surface_exposure_threshold: 0.2
    fdr_alpha: 0.05
```

### Settings Support Matrix

| Setting | contacts | binding_free_energy | polymer_affinity | Default |
|---------|----------|---------------------|------------------|---------|
| `fdr_alpha` | ✓ | ✓ | ✓ | 0.05 |
| `min_effect_size` | ✓ | — | — | 0.5 |
| `top_residues` | ✓ | — | — | 10 |

### Setting Descriptions

- **`fdr_alpha`** — Significance threshold for the Benjamini-Hochberg (BH)
  false discovery rate correction applied to pairwise p-values. Lower values
  are more conservative.
- **`min_effect_size`** — Minimum Cohen's d required for practical
  significance. Pairs that meet or exceed this threshold are highlighted with
  "†" in formatted output; all pairs are shown regardless.
- **`top_residues`** — Maximum number of contacted residues shown per
  condition, ranked by aggregated `contact_fraction_mean`. Affects both saved
  JSON and CLI output.

## Stable Plugin Keys

Stable analysis plugins:

- `rmsd`
- `rg`
- `rmsf`
- `contacts`
- `distances`
- `catalytic_triad`
- `secondary_structure`
- `sasa`
- `hydrogen_bonds` (aliases: `hbonds`, `hbond`)

Experimental but still available:

- binding preference through `contacts`
- `exposure`
- `binding_free_energy`
- `polymer_affinity`
- `polymer_bridging` (alias: `bridging`)

## Plugin Summary Table

| Plugin | Default compare? | Primary metric | Key feature | Statistical method |
|--------|-----------------|----------------|-------------|-------------------|
| `rmsd` | Yes | `mean_rmsd` | Backbone stability over time | FDR-corrected pairwise t-tests + ANOVA |
| `rg` | Yes | `mean_rg` | Protein compactness | FDR-corrected pairwise t-tests + ANOVA |
| `rmsf` | Yes | `mean_rmsf` | Per-residue flexibility | FDR-corrected pairwise t-tests + ANOVA |
| `contacts` | No (custom) | Coverage + contact fraction | Per-residue contact mapping | FDR-corrected pairwise t-tests per residue |
| `distances` | No (custom) | Multiple distance metrics | Named distance pairs | Per-distance t-tests + ANOVA |
| `catalytic_triad` | Yes | `mean_triad_proximity` | Active-site geometry | FDR-corrected pairwise t-tests + ANOVA |
| `secondary_structure` | Yes | `helix_fraction` | Secondary structure content | FDR-corrected pairwise t-tests + ANOVA |
| `sasa` | No (custom) | Per-run mean SASA | Multi-run target/context model | Per-run pairwise t-tests + ANOVA |
| `hydrogen_bonds` | Yes | `mean_hbonds_per_frame` per summary | Flexible named groups + summaries + composition analysis | FDR-corrected pairwise t-tests + ANOVA |
| `exposure` | No (custom) | Exposure dynamics metrics | Time-resolved surface exposure | Custom statistical pipeline |
| `binding_free_energy` | No (custom) | Per-contact ΔG_sel | Free energy decomposition | Custom statistical pipeline |
| `polymer_affinity` | No (custom) | Total interaction score | Combined contact + energetic scoring | Custom statistical pipeline |
| `polymer_bridging` | No (custom) | Bridging event counts | Polymer-mediated inter-chain contacts | Custom statistical pipeline |

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
| `polyzymd compare submit ANALYSIS` | Submit a SLURM DAG for one analysis plugin |
| `polyzymd compare status ANALYSIS` | Show status of a submitted SLURM DAG |
| `polyzymd compare finalize ANALYSIS` | Run comparison + plotting from on-disk aggregated results |

## Common Stable Commands

All commands below assume you are inside the pixi environment
(`pixi shell -e build`) or are prefixed with `pixi run -e build`.

```bash
polyzymd compare run rmsd
polyzymd compare run rg
polyzymd compare run rmsf
polyzymd compare run contacts
polyzymd compare run distances
polyzymd compare run triad
polyzymd compare run sasa
polyzymd compare run hydrogen_bonds  # aliases: hbonds, hbond
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
comparison/rmsd/result.json
comparison/rg/result.json
comparison/rmsf/result.json
comparison/contacts/result.json
comparison/distances/result.json
comparison/catalytic_triad/result.json
comparison/sasa/result.json
comparison/hydrogen_bonds/result.json
```

## Plotting Smoke Test

For a final smoke test after comparisons finish:

```bash
polyzymd compare plot-all --list-available
polyzymd compare plot-all
```

## Plugin-Specific Metadata Fields

Some plugins include additional metadata in their comparison output beyond the
standard ranking and statistical fields. These fields are **additive
diagnostics** — they do not affect rankings, p-values, or effect sizes.

### RMSD Convergence Output

The RMSD plugin includes per-run convergence diagnostics generated by the
sliding-window convergence heuristic in `analyses/shared/convergence.py`.
These fields appear in the per-condition summaries within
`comparison/rmsd/result.json`:

| Field | Type | Description |
|-------|------|-------------|
| `convergence_fraction` | `float` | Fraction of replicates that converged (0.0–1.0) |
| `n_converged_replicates` | `int` | Count of replicates where sustained convergence was detected |
| `mean_convergence_time_ns` | `float \| null` | Mean convergence time across converged replicates (ns) |
| `median_convergence_time_ns` | `float \| null` | Median convergence time across converged replicates (ns) |

:::{note}
Convergence metadata is purely informational. It does not influence the RMSD
ranking, pairwise t-tests, ANOVA, or effect-size calculations. Use it to
identify conditions where one or more replicates failed to reach a stable
plateau, which may warrant longer production runs or additional replicates.
:::

## Statistical Terms

- `p-value`: significance of the observed difference under the null hypothesis
- `Cohen's d`: effect size magnitude
- `ANOVA`: omnibus test across multiple conditions
- `SEM`: standard error of the mean across replicates
- `Benjamini-Hochberg (BH)`: step-up procedure for controlling the false
  discovery rate across multiple hypothesis tests
- `Adjusted p-value (p_adj)`: p-value corrected for multiple comparisons via
  the BH procedure
- `False Discovery Rate (FDR)`: expected proportion of false positives among
  rejected hypotheses
- `Effect size threshold`: minimum Cohen's d required for a pairwise difference
  to be considered practically significant

For interpretation guidance rather than lookup, see:

- [Statistical Best Practices for Analysis](../explanation/analysis_statistics_best_practices.md)
- [How to Compare Simulation Conditions](../how_to/analysis_compare_conditions.md)
- [Post-Hoc Testing Reference](posthoc_testing.md) — full post-hoc method details, output fields, and edge cases

:::{admonition} Version 1.3.0 migration
:class: note

In versions before 1.3.0, `fdr_alpha`, `min_effect_size`, and `top_residues`
were accepted by the YAML schema but silently ignored at runtime. As of 1.3.0,
these settings actively affect comparison results: pairwise p-values are
BH-corrected, effect-size thresholds annotate pairwise output, and
`top_residues` limits both CLI output and saved JSON.
:::
