# Analysis Module Rules

## Overview

The analysis system lives in the `analyses/` package:

```
src/polyzymd/analyses/
├── base.py                 # Analysis ABC, context objects, result models
├── discovery.py            # pkgutil-based auto-discovery
├── orchestrator.py         # Framework engine: compute → aggregate → compare → plot
├── stats.py                # default_scalar_comparison(), format_scalar_comparison()
├── shared/                 # Reusable utilities (TrajectoryLoader, alignment, statistics, etc.)
│   ├── binding_preference.py   # Cross-plugin compute (contacts, BFE, polymer_affinity)
│   ├── binding_preference_helpers.py  # Orchestration helpers for binding preference
│   └── surface_exposure.py     # SASA-based surface exposure utility
├── _results_base.py        # Base result model (shared, stays at top level)
├── rmsf/                   # RMSF plugin sub-package
│   ├── __init__.py         #   RMSFAnalysis plugin class
│   ├── _calculator.py      #   RMSFCalculator
│   ├── _results.py         #   Result models
│   └── _plotting.py        #   Legacy plot helpers
├── distances/              # Distances plugin sub-package
├── catalytic_triad/        # Catalytic triad plugin sub-package
├── secondary_structure/    # Secondary structure plugin sub-package
├── contacts/               # Contacts plugin sub-package
├── exposure/               # Exposure dynamics plugin sub-package
├── binding_free_energy/    # Binding free energy plugin (custom compare)
└── polymer_affinity/       # Polymer affinity plugin (custom compare)
```

The private `_calculator.py` modules (inside each sub-package) provide underlying computation that
some plugins delegate to. New plugins can compute directly in
`compute_replicate()` — there is no separate `analysis/` package.

### How to Add a New Analysis

1. Run `polyzymd new-analysis <name>` to scaffold the plugin package and tests, OR
   create `src/polyzymd/analyses/<name>/` sub-package manually
2. Define a `Settings` class (Pydantic v2 `BaseModel`)
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

**Simple path** (rg, RMSF, catalytic_triad, secondary_structure):
- Implement `extract_metrics()` to return `dict[str, MetricValue]`
- Implement `_deserialize_result()` so the default `compare()` can load
  aggregated JSON results
- Framework handles t-tests, ANOVA, ranking, and formatting automatically

**Custom path** (contacts, distances, exposure, BFE, polymer_affinity):
- Override `compare()` entirely — return your own Pydantic model with `.save()`
- The returned model must be saveable/loadable

### Context Objects

Plugins receive framework-provided context objects — never load configs yourself:

| Context | Passed To | Key Attributes |
|---------|-----------|----------------|
| `ReplicateContext` | `compute_replicate()` | `.sim_config`, `.settings`, `.output_dir`, `.replicate`, `.recompute`, `.result_path` |
| `AggregateContext` | `aggregate()` | `.condition`, `.replicates`, `.output_dir`, `.settings`, `.result_path` |
| `ComparisonContext` | `compare()` | `.conditions`, `.analysis_dirs`, `.results_dir`, `.effective_control`, `.settings` |
| `PlotContext` | `plot()` | `.conditions`, `.analysis_dirs`, `.output_dir`, `.settings`, `.plot_settings` |

### Result Saving

Existing plugins save results to disk explicitly — using custom filenames
for per-replicate caching (e.g. `rmsf_eq10ns.json`) and `ctx.result_path`
for aggregated results. The orchestrator has a **fallback** auto-save that
writes to `ctx.result_path` only if the file doesn't already exist.

Simple plugins can skip manual saves and rely on the fallback. Plugins that
want equilibration-aware caching should save explicitly (see `rmsf/` for the
pattern).

### Return Types: Dicts vs Pydantic Models

Both plain dicts and Pydantic `BaseModel` instances are supported as return
types from `compute_replicate()` and `aggregate()`. Dicts are recommended
for new plugins; Pydantic models are useful for complex results that need
validation or NPZ sidecar storage.

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

This classification lives in `analyses/shared/metric_type.py` and is used by
the underlying calculators (`_calculator.py` inside each sub-package), not directly by the plugin
system.

## Compute Layer (Private `_calculator.py` Modules)

The private modules inside `analyses/` provide calculators that some plugins
delegate to:

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
3. ~~Migrate legacy formatters into plugin `format()` methods~~ (DONE — Phase 7)
4. Add test templates for new contributor plugins
