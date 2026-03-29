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
│   ├── _plotters.py        #   Plotting functions (extracted from __init__.py)
│   ├── _results.py         #   Result models
│   └── _comparison_results.py  # Comparison result model
├── distances/              # Distances plugin sub-package
│   ├── __init__.py         #   DistancesAnalysis + DistanceCalculator
│   └── _plotters.py        #   Plotting functions
├── catalytic_triad/        # Catalytic triad plugin sub-package (default-compare lifecycle)
│   ├── __init__.py         #   CatalyticTriadAnalysis plugin class
│   ├── _plotters.py        #   Plotting functions (KDE panels, threshold bars)
│   └── _results.py         #   Result models
├── secondary_structure/    # Secondary structure plugin sub-package
│   ├── __init__.py         #   SecondaryStructureAnalysis plugin class
│   └── _plotters.py        #   Plotting functions
├── contacts/               # Contacts plugin sub-package
│   ├── __init__.py         #   ContactsAnalysis plugin class
│   ├── _plotters.py        #   Plotting functions (19 functions)
│   ├── _results.py         #   Result models (widely imported — do not modify)
│   ├── _comparison_results.py  # Comparison result model
│   ├── _aggregator.py      #   Aggregation logic
│   └── _formatters.py      #   CLI formatting
├── exposure/               # Exposure dynamics plugin sub-package
├── binding_free_energy/    # Binding free energy plugin (custom compare)
└── polymer_affinity/       # Polymer affinity plugin (custom compare)
```

Each plugin is a self-contained package. All established plugins extract
plotting functions into a `_plotters.py` module to keep `__init__.py` focused
on the Analysis lifecycle. New plugins can start with all logic in
`__init__.py` and extract plotting later as complexity grows.

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

### Extracting Plotting to `_plotters.py`

Established plugins separate plotting functions into a `_plotters.py` module
within the plugin package. This keeps `__init__.py` focused on the lifecycle
methods (`compute_replicate`, `aggregate`, `compare`) while plotting functions
live in a dedicated file.

**Current state** (all 5 established plugins have `_plotters.py`):
- `rmsf/_plotters.py` — 8 plotting functions
- `secondary_structure/_plotters.py` — 7 plotting functions + 3 constants
- `distances/_plotters.py` — 10 plotting functions
- `contacts/_plotters.py` — 18 functions (8 data loaders + 10 plotters)
- `catalytic_triad/_plotters.py` — 7 functions (KDE panels, threshold bars, 2D KDE)

**Pattern:** The plugin's `plot()` method calls functions from `_plotters.py`:

```python
# In __init__.py
from ._plotters import plot_comparison_bars, plot_per_residue

def plot(self, ctx: PlotContext) -> list[Path]:
    data, labels = self._build_plot_data(ctx)
    paths = []
    paths.append(plot_comparison_bars(data, labels, ctx.output_dir, ctx.plot_settings))
    paths.append(plot_per_residue(data, labels, ctx.output_dir, ctx.plot_settings))
    return [p for p in paths if p is not None]
```

**New plugins** should start with plotting inline in `plot()` and extract to
`_plotters.py` when the module exceeds ~500 lines or has 3+ plot functions.

### Two Comparison Paths

**Simple path** (rmsf, catalytic_triad, secondary_structure):
- Implement `extract_metrics()` to return `dict[str, MetricValue]`
- Set `AggregatedResultClass` if using Pydantic models (dict plugins work
  automatically via `json.loads()`)
- Framework handles loading, t-tests, ANOVA, ranking, and formatting automatically

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

### Return Types: Pydantic Models vs Dicts

**Recommended:** Return typed Pydantic `BaseModel` instances from
`compute_replicate()` and `aggregate()`. This gives you validation, type
safety, and IDE autocomplete. Plain dicts are also supported for rapid
prototyping.

The orchestrator issues a warning if return types are not `BaseModel` or `dict`.

### Result Deserialization

The default `compare()` path calls `_load_aggregated_result()` →
`_deserialize_result()` to load aggregated results from disk. The base
implementation handles both paths automatically:

- If `AggregatedResultClass` is set: uses `.load(path)` or `.model_validate_json()`
- Otherwise: uses `json.loads()` for plain-dict results

You only need to override `_deserialize_result()` for non-standard
deserialization (e.g. NPZ sidecars, custom migrations).

### PlotContext.plot_settings Guarantee

The orchestrator guarantees that `PlotContext.plot_settings` is always a valid
`PlotSettings` instance — never `None`. Plugins should NOT guard against `None`:

```python
# CORRECT — trust the guarantee
def plot(self, ctx: PlotContext) -> list[Path]:
    theme = get_theme(ctx.plot_settings)

# WRONG — unnecessary guard (removed from all plugins)
def plot(self, ctx: PlotContext) -> list[Path]:
    if ctx.plot_settings is None:
        ...
```

### Loading Results in `plot()`

Use `_build_plot_data()` to collect per-condition paths, then
`_load_aggregated_result()` to load each condition's result:

```python
data, labels = self._build_plot_data(ctx)
for label in labels:
    if label not in data or label == "__meta__":
        continue
    agg_dir = data[label]["aggregated_dir"]
    summary = self._load_aggregated_result(agg_dir)
    if summary is not None:
        # ... plot data ...
```

## MetricType System (Autocorrelation Handling)

Based on LiveCoMS best practices (Grossfield et al., 2018), metrics are classified
by how they should handle autocorrelated MD data:

| MetricType | Frame Strategy | Uncertainty Strategy | Examples |
|------------|---------------|---------------------|----------|
| **MEAN_BASED** | Use ALL frames | Correct SEM with N_eff = N/g | Contact fraction, triad proximity |
| **VARIANCE_BASED** | Subsample by 2τ | Standard formula on independent samples | RMSF, fluctuation metrics |

This classification lives in `analyses/shared/metric_type.py` and is used
internally by the compute logic within each plugin.

## Key Classes

- **`Analysis`** — Plugin ABC (in `analyses/base.py`)
- **`ParallelContactAnalyzer`** — Polymer-protein contacts (in `contacts/`)
- **`DistanceCalculator`** — Inter-group distance tracking (in `distances/__init__.py`)
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
