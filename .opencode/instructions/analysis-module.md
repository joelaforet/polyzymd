# Analysis Module Rules

## Overview

The analysis system lives in the `analyses/` package:

```
src/polyzymd/analyses/
├── base.py                 # Public facade: Analysis, contexts, result models
├── _analysis_runner.py     # Internal runner-backed replicate dispatch
├── _analysis_compare.py    # Internal default comparison implementation
├── _analysis_io.py         # Internal result paths, cache, serialization helpers
├── _analysis_contract.py   # Internal plugin contract validation
├── _contexts.py            # Internal definitions re-exported by base.py
├── _comparison_models.py   # Internal models re-exported by base.py
├── discovery.py            # pkgutil-based auto-discovery
├── orchestrator.py         # Framework engine: compute → aggregate → compare → plot
├── stats.py                # default_scalar_comparison(), format_scalar_comparison()
├── shared/                 # Reusable utilities (TrajectoryLoader, alignment, statistics, etc.)
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
│   ├── __init__.py         #   Public ContactsAnalysis facade
│   ├── _cache.py           #   Internal cache helpers
│   ├── _lifecycle.py       #   Internal lifecycle helpers
│   ├── _filters.py         #   Internal condition filtering
│   ├── _comparison.py      #   Internal custom comparison implementation
│   ├── _plotting.py        #   Internal plot lifecycle orchestration
│   ├── _plotters.py        #   Internal plotting functions
│   ├── _results.py         #   Internal result models
│   ├── _comparison_results.py  # Internal comparison result model
│   ├── _aggregator.py      #   Internal aggregation logic
│   ├── _formatters.py      #   Internal CLI formatting
│   ├── _identity.py        #   Internal settings fingerprints
│   ├── _runner.py          #   Internal trajectory runner and analyzer
│   └── _paths.py           #   Internal result path helpers
```

`polyzymd.analyses.base` is the supported public import surface. It re-exports
`Analysis`, lifecycle contexts, metric descriptors, and comparison models from
private framework modules. Contributors should import from
`polyzymd.analyses.base` and should not import `_analysis_*`, `_contexts.py`, or
`_comparison_models.py` directly.

Each plugin is a self-contained package. All established plugins extract
plotting functions into a `_plotters.py` module to keep `__init__.py` focused
on Analysis lifecycle wiring. Contacts is larger, so `ContactsAnalysis` stays
public in `contacts/__init__.py` while cache, filtering, lifecycle, comparison,
plotting, result, and runner helpers live in private `contacts/_*.py` modules.
These modules are implementation details, not
contributor API. Public contributor plugins use the runner-backed lifecycle;
trajectory-native plugins should isolate MDAnalysis trajectory logic in a
dedicated module such as `_runner.py`.

### How to Add a New Analysis

1. Run `polyzymd new-analysis <name>` to scaffold the plugin package and tests, OR
   create `src/polyzymd/analyses/<name>/` sub-package manually
2. Define a `Settings` class (Pydantic v2 `BaseModel`)
3. Subclass `Analysis` and choose the lifecycle mode for your plugin
4. When `has_compute_stage=True`, implement `build_runner()` +
   `summarize_replicate()` for the MDAnalysis-first runner path; keep lifecycle
   wiring in `__init__.py` and put runner logic in a dedicated module such as
   `_runner.py`
5. If the plugin is compare-only, set `has_compute_stage=False`
6. Implement `aggregate()` only when `has_aggregate_stage=True`
7. Done — framework discovers it via `pkgutil` (no registries, no imports)

### Required Class Variables

```python
class MyAnalysis(Analysis):
    name: ClassVar[str] = "my_analysis"       # Used in CLI and config
    Settings: ClassVar[type] = MySettings      # Or as inner class
```

### Lifecycle Hooks

Required hooks depend on the plugin mode:

| Hook | When Used | Signature / Notes |
|------|-----------|-------------------|
| `build_runner()` + `summarize_replicate()` | Public compute-stage plugins with `has_compute_stage=True` | MDAnalysis owns per-trajectory iteration while PolyzyMD owns caching, ensemble aggregation, and comparison workflow |
| `run_replicate()` | Advanced/internal cache or sidecar wrappers only | Framework entry point used by the orchestrator; new contributor plugins should not override it |
| `aggregate()` | Only when `has_aggregate_stage=True` | `(ctx: AggregateContext, results: Sequence[Any]) -> Any` |

Compare-only or no-compute plugins set `has_compute_stage=False` and skip the
runner-backed compute path.

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
methods for the chosen plugin mode while plotting functions live in a dedicated
file.

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

**Custom path** (contacts, distances):
- Override `compare()` entirely — return your own Pydantic model with `.save()`
- The returned model must be saveable/loadable

### Context Objects

Plugins receive framework-provided context objects — never load configs yourself:

| Context | Passed To | Key Attributes |
|---------|-----------|----------------|
| `ReplicateContext` | `build_runner()`, `summarize_replicate()`, advanced `run_replicate()` wrappers | `.sim_config`, `.settings`, `.output_dir`, `.replicate`, `.recompute`, `.result_path` |
| `AggregateContext` | `aggregate()` | `.condition`, `.replicates`, `.output_dir`, `.settings`, `.result_path` |
| `ComparisonContext` | `compare()` | `.conditions`, `.analysis_dirs`, `.results_dir`, `.effective_control`, `.settings` |
| `PlotContext` | `plot()` | `.conditions`, `.analysis_dirs`, `.output_dir`, `.settings`, `.plot_settings` |

### Result Saving

Existing plugins save results to disk explicitly — using custom filenames
for per-replicate caching (e.g. `rmsf_eq10ns.json`) and `ctx.result_path`
for aggregated results. The orchestrator also provides fallback result saving
for plugins that return serializable outputs without writing them manually.

Simple plugins can skip manual saves and rely on the fallback. Plugins that
want equilibration-aware caching should save explicitly (see `rmsf/` for the
pattern).

### Return Types: Pydantic Models vs Dicts

Lifecycle hooks must return a Pydantic `BaseModel` instance or a plain `dict`
(`summarize_replicate()`, plus `aggregate()` when used). `compare()` may also
return `None` when no comparison result is produced. Invalid return types are
enforced as plugin contract failures and raise `PluginContractError`.

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

This section is a planned future enhancement and is not implemented in the
current codebase.

Based on LiveCoMS best practices (Grossfield et al., 2018), metrics may be
classified by how they should handle autocorrelated MD data.

Planned (not yet implemented):

| MetricType | Frame Strategy | Uncertainty Strategy | Examples |
|------------|---------------|---------------------|----------|
| **MEAN_BASED** | Use ALL frames | Correct SEM with N_eff = N/g | Contact fraction, triad proximity |
| **VARIANCE_BASED** | Subsample by 2τ | Standard formula on independent samples | RMSF, fluctuation metrics |

## Key Classes

- **`Analysis`** — Plugin ABC (in `analyses/base.py`)
- **`ParallelContactAnalyzer`** — Polymer-protein contacts (in `contacts/`)
- **`DistanceCalculator`** — Inter-group distance tracking (in `distances/__init__.py`)
- **`MolecularSelector`** — Strategy for atom group selection
- **`ContactCriteria`** — Strategy for contact definitions
- **`MetricType`** — Planned autocorrelation handling strategy (future)

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
| `rmsd` | No (custom) | Per-run mean RMSD (lower = more stable) |
| `rg` | No (custom) | Per-run mean Rg (lower = more compact) |
| `sasa` | No (custom) | Solvent-accessible surface area metrics |
| `distances` | No (custom) | Multiple distance metrics |
| `contacts` | No (custom) | Coverage + contact fraction |
| `hydrogen_bonds` | No (custom) | H-bond occupancy + lifetime |

## Issue #20 — Remaining TODOs

GitHub Issue #20 tracks the analysis module roadmap. The plugin system
(Phases 1-5) is complete. Remaining items:

1. Fix cache hash mismatch warning (prints 66x)
2. Resolve contact criteria cutoff disagreement (4.0A vs 4.5A)
3. ~~Migrate legacy formatters into plugin `format()` methods~~ (DONE — Phase 7)
4. Add test templates for new contributor plugins
