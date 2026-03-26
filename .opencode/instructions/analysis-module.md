# Analysis Module Rules

## Overview

PolyzyMD analysis is split across two packages:

| Package | Role |
|---------|------|
| `analysis/` | Per-condition calculators, results, aggregation — the **compute** layer |
| `analyses/` | **Plugin system** — unified lifecycle (compute → aggregate → compare → plot → format) |

New analysis types are added as **plugins in `analyses/`**. The `analysis/`
package provides the underlying computation that plugins delegate to.

## Plugin System (`analyses/`)

### Module Structure

```
src/polyzymd/analyses/
├── __init__.py             # Public API: get_analysis, list_analyses, run_comparison
├── base.py                 # Analysis ABC, context objects, result models
├── discovery.py            # pkgutil-based auto-discovery
├── orchestrator.py         # Framework engine: compute → aggregate → compare → plot
├── stats.py                # default_scalar_comparison(), format_scalar_comparison()
├── rmsf.py                 # RMSF plugin (simplest default-compare example)
├── catalytic_triad.py      # Catalytic triad plugin (default compare)
├── secondary_structure.py  # Secondary structure plugin (default compare)
├── distances.py            # Distances plugin (custom compare)
├── contacts.py             # Contacts plugin (custom compare)
├── exposure.py             # Exposure dynamics plugin (custom compare)
├── binding_free_energy.py  # Binding free energy plugin (custom compare)
└── polymer_affinity.py     # Polymer affinity plugin (custom compare)
```

### How to Add a New Analysis

1. Create `src/polyzymd/analyses/<name>.py`
2. Define a `Settings` inner class (Pydantic v2 `BaseModel`)
3. Subclass `Analysis` and implement `compute_replicate()` and `aggregate()`
4. Done — framework discovers it via `pkgutil` (no registries, no imports)

### Required Class Variables

```python
class MyAnalysis(Analysis):
    name: ClassVar[str] = "my_analysis"       # Used in CLI and config
    Settings: ClassVar[type] = MySettings      # Or as inner class
```

### Required Methods

| Method | When Called | Signature |
|--------|-----------|-----------|
| `compute_replicate()` | Once per replicate per condition | `(ctx: ReplicateContext, replicate: int) -> Any` |
| `aggregate()` | Once per condition (after all replicates) | `(ctx: AggregateContext, results: Sequence[Any]) -> Any` |

### Optional Overrides

| Method | Default Behavior | Override When |
|--------|-----------------|---------------|
| `extract_metrics()` | Returns `{}` | Using default `compare()` — return `dict[str, MetricValue]` for automatic t-tests/ANOVA |
| `compare()` | Calls `extract_metrics()` + `default_scalar_comparison()` | Multi-metric or entry-table comparisons (e.g. contacts, distances) |
| `filter_conditions()` | Keeps all conditions | Excluding conditions (e.g. no-polymer for contacts) |
| `plot()` | Returns `[]` (no plots) | Custom matplotlib visualizations |
| `format()` | JSON dump or legacy formatter | Custom CLI output formatting |

### Two Comparison Paths

**Simple path** (RMSF, catalytic_triad, secondary_structure):
- Implement `extract_metrics()` to return `dict[str, MetricValue]`
- Framework handles t-tests, ANOVA, ranking, and formatting automatically
- Also implement `_deserialize_result()` so the default `compare()` can load
  aggregated JSON results

**Custom path** (contacts, distances, exposure, BFE, polymer_affinity):
- Override `compare()` entirely — return your own Pydantic model with `.save()`
- The returned model must be saveable/loadable

### Context Objects

Plugins receive framework-provided context objects — never load configs yourself:

| Context | Passed To | Key Attributes |
|---------|-----------|----------------|
| `ReplicateContext` | `compute_replicate()` | `.sim_config`, `.settings`, `.output_dir`, `.replicate`, `.recompute` |
| `AggregateContext` | `aggregate()` | `.condition`, `.replicates`, `.output_dir`, `.settings` |
| `ComparisonContext` | `compare()` | `.conditions`, `.analysis_dirs`, `.results_dir`, `.effective_control`, `.settings` |
| `PlotContext` | `plot()` | `.conditions`, `.analysis_dirs`, `.output_dir`, `.settings`, `.plot_settings` |

### Important: `_deserialize_result()` Footgun

If you use the **default `compare()` path**, the framework calls
`_load_aggregated_result()` → `_deserialize_result()` to load your aggregated
JSON. The base implementation raises `NotImplementedError`. You MUST either:

1. Implement `_deserialize_result(path)` to load your aggregated result model, OR
2. Override `compare()` entirely (then `_deserialize_result()` is never called)

## MetricType System (Autocorrelation Handling)

Based on LiveCoMS best practices (Grossfield et al., 2018), metrics are classified
by how they should handle autocorrelated MD data:

| MetricType | Frame Strategy | Uncertainty Strategy | Examples |
|------------|---------------|---------------------|----------|
| **MEAN_BASED** | Use ALL frames | Correct SEM with N_eff = N/g | Contact fraction, triad proximity |
| **VARIANCE_BASED** | Subsample by 2τ | Standard formula on independent samples | RMSF, fluctuation metrics |

This classification lives in `analysis/core/metric_type.py` and is used by the
underlying calculators in `analysis/`, not directly by the plugin system.

## Compute Layer (`analysis/`)

The `analysis/` package provides calculators that plugins delegate to:

### Key Classes

- **`RMSFCalculator`** — Per-residue RMSF calculation
- **`DistanceCalculator`** — Inter-group distance tracking
- **`CatalyticTriadAnalyzer`** — Catalytic triad integrity
- **`ParallelContactAnalyzer`** — Polymer-protein contacts
- **`SecondaryStructureCalculator`** — DSSP secondary structure
- **`MolecularSelector`** — Strategy for atom group selection
- **`ContactCriteria`** — Strategy for contact definitions
- **`MetricType`** — Autocorrelation handling strategy

### Caching

Analysis results are cached to avoid recomputation. The cache uses file hashes
to detect changes. **Known issue:** the hash mismatch warning currently prints
66+ times instead of once (see known-issues.md).

## Existing Plugins

| Plugin | Default compare? | Primary metric |
|--------|-----------------|----------------|
| `rmsf` | Yes | `mean_rmsf` (lower = more stable) |
| `catalytic_triad` | Yes | `mean_triad_proximity` (lower = closer) |
| `secondary_structure` | Yes | `helix_fraction` (higher = more structured) |
| `distances` | No (custom) | Multiple distance metrics |
| `contacts` | No (custom) | Coverage + contact fraction |
| `exposure` | No (custom) | Exposure dynamics metrics |
| `binding_free_energy` | No (custom) | Per-contact ΔG_sel |
| `polymer_affinity` | No (custom) | Total interaction score |

## Issue #20 — Remaining TODOs

GitHub Issue #20 tracks the analysis module roadmap. The plugin system
(Phases 1-5) is complete. Remaining items:

1. Fix cache hash mismatch warning (prints 66x)
2. Resolve contact criteria cutoff disagreement (4.0A vs 4.5A)
3. Migrate legacy formatters into plugin `format()` methods
4. Add test templates for new contributor plugins
