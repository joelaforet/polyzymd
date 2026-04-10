# Extending the Analysis Framework

This guide shows you how to add a new analysis type to PolyzyMD.
By the end, you will have a working analysis plugin that can compute metrics
per replicate, aggregate across replicates, compare conditions, generate
figures, and display results on the CLI.

```{note}
This guide covers the **plugin system** in ``analyses/``.
Each plugin is a self-contained package — you implement the science
(compute, aggregate) and the framework handles discovery, iteration,
result saving, comparison statistics, and CLI wiring.
```

## Prerequisites

- A working pixi environment: `pixi install -e build`
- Basic familiarity with Pydantic v2 `BaseModel`
- Understanding of the PolyzyMD chain convention (A=protein, B=substrate,
  C=polymer)

## How the Plugin System Works

Every analysis in PolyzyMD is a single class that subclasses `Analysis` and
lives in a package under `src/polyzymd/analyses/<name>/`. The framework discovers
plugins automatically — no registry edits, no imports, no boilerplate. You drop
a package, and the CLI can run it.

The lifecycle for each analysis is:

```text
For each condition:
  For each replicate:
    result = compute_replicate(ctx, replicate)       # your code
  aggregated = aggregate(ctx, [replicate_results])   # your code

filter_conditions(conditions)  →  filtered list
compare(ctx)                   →  ComparisonResult
plot(ctx)                      →  [figure paths]
format(result, "text")         →  CLI output string
```

You implement the science (compute, aggregate). The framework provides
discovery, replicate iteration, result saving, comparison statistics, and CLI
wiring.

## Anatomy of a Plugin

Before seeing the full code, here is what each piece does and why it exists.

### `name` (required)

A unique, lowercase string identifier. Used in CLI commands
(`polyzymd compare run solvent_contacts`), config files, and output directories.

### `Settings` (required)

A Pydantic v2 `BaseModel` inner class. Users configure your analysis through
the `plugins:` section of `comparison.yaml`, and the framework deserializes
those settings into your `Settings` model. Provide sensible defaults for every
field so the analysis works out of the box.

### `compute_replicate(ctx, replicate)` (required)

Runs once per replicate per condition. The framework passes a
`ReplicateContext` containing everything you need: the loaded
`SimulationConfig`, your `Settings`, output paths, and equilibration time. Use
`TrajectoryLoader` from `analyses/shared/` to load trajectories — it handles
topology discovery, trajectory segment daisy-chaining, timestep detection, and
equilibration skipping.

Return a dict or Pydantic model. The framework collects these for `aggregate()`.

### `aggregate(ctx, results)` (required)

Runs once per condition after all replicates complete. Receives the list of
objects your `compute_replicate()` returned. Compute summary statistics
(mean, SEM) and return the aggregated result. The framework auto-saves this
to `ctx.result_path`.

### `extract_metrics(summary)` (optional, but recommended)

Exposes scalar metrics for the **default comparison path**. Return a dict of
`MetricValue` objects, and the framework automatically runs t-tests, computes
Cohen's d, runs ANOVA (≥3 conditions), and ranks conditions. This is the
simplest way to get a fully statistical comparison — no need to implement
`compare()`.

### `compare(ctx)` (optional override)

Override only for complex multi-metric or entry-table comparisons. If you
implement `extract_metrics()`, you typically do **not** need to override
`compare()`.

### `plot(ctx)` (optional)

Generate matplotlib figures. Use `self._build_plot_data(ctx)` to collect
per-condition paths, then `self._load_aggregated_result(agg_dir)` to load each
condition's data. Use shared plotting helpers (`get_theme`, `apply_axis_style`,
`save_figure`, `get_colors`) for consistent theming. Return a list of figure
paths. The scaffold generates a working `plot()` method automatically.

### `format(result, output_format)` (optional)

Custom CLI display. Use `format_scalar_comparison()` from `analyses.stats`
for formatted tables with rankings, effect sizes, and significance stars.

## The Complete Example

Here is a single, complete plugin file that uses real framework utilities.
This is the pattern you should follow for new analyses.

```{important}
The `SolventContactsAnalysis` implementation below is a walkthrough example for
a hypothetical `solvent_contacts` plugin name that does not collide with
built-in plugins.
```

Create `src/polyzymd/analyses/solvent_contacts/` as a package in your own branch or downstream
project, or use the scaffold command to generate the boilerplate automatically
(the scaffold includes compute, aggregate, comparison, plotting, and tests):

```bash
# Default: dict-based results (simplest starting point)
polyzymd new-analysis solvent_contacts

# Typed Pydantic result models (better for complex analyses)
polyzymd new-analysis solvent_contacts --style pydantic
```

The `--style` flag controls how results are structured:

- **`dict`** (default) — Uses plain dicts for replicate/aggregated results.
  No result model classes needed. Best for getting started quickly.
- **`pydantic`** — Generates typed `ReplicateResult` and `AggregatedResult`
  Pydantic models with `AggregatedResultClass`. Better when you need
  validation, IDE autocomplete, or NPZ sidecar storage.

If creating manually:

```python
"""Solvent contacts analysis plugin."""
from __future__ import annotations

import logging
from pathlib import Path
from typing import Any, ClassVar, Sequence

import numpy as np
from pydantic import BaseModel, Field

from polyzymd.analyses.base import (
    AggregateContext,
    Analysis,
    MetricValue,
    PlotContext,
    ReplicateContext,
)
from polyzymd.analyses.shared import TrajectoryLoader, convert_time, parse_time_string
from polyzymd.analyses.stats import format_scalar_comparison

logger = logging.getLogger(__name__)


class SolventContactsAnalysis(Analysis):
    """Solvent contacts around a selected protein region."""

    name: ClassVar[str] = "solvent_contacts"

    class Settings(BaseModel):
        """Settings for solvent contact-fraction analysis.

        Attributes
        ----------
        selection : str
            MDAnalysis selection string for target atoms
        solvent_selection : str
            MDAnalysis selection for solvent atoms used in contacts
        cutoff : float
            Contact cutoff distance in Angstrom
        """

        selection: str = Field(
            default="protein and name CA",
            description="MDAnalysis atom selection for target atoms",
        )
        solvent_selection: str = Field(
            default="resname HOH and name O",
            description="MDAnalysis atom selection for solvent contact atoms",
        )
        cutoff: float = Field(
            default=4.5,
            gt=0,
            description="Distance cutoff in Angstrom for a contact",
        )

    def compute_replicate(self, ctx: ReplicateContext, replicate: int) -> dict[str, Any]:
        """Compute solvent contact fraction for one replicate.

        Uses TrajectoryLoader for topology discovery, trajectory segment
        daisy-chaining, and equilibration frame skipping
        """
        # Lazy-import heavy third-party dependency
        import MDAnalysis as mda  # noqa: F401
        from MDAnalysis.analysis.distances import distance_array

        loader = TrajectoryLoader(ctx.sim_config)
        u = loader.load_universe(replicate)

        eq_value, eq_unit = parse_time_string(ctx.equilibration)
        timestep_ps = loader.get_timestep(replicate, unit="ps")
        eq_time_ps = convert_time(eq_value, eq_unit, "ps")
        start_frame = int(eq_time_ps / timestep_ps)

        target_atoms = u.select_atoms(ctx.settings.selection)
        solvent_atoms = u.select_atoms(ctx.settings.solvent_selection)
        contact_fractions = []
        for _ts in u.trajectory[start_frame:]:
            if len(target_atoms) == 0 or len(solvent_atoms) == 0:
                contact_fractions.append(0.0)
                continue
            dists = distance_array(target_atoms.positions, solvent_atoms.positions)
            contact_pairs = np.sum(dists <= ctx.settings.cutoff)
            normalizer = len(target_atoms) * len(solvent_atoms)
            contact_fractions.append(float(contact_pairs / normalizer))

        return {
            "mean_contact_fraction": float(np.mean(contact_fractions)),
            "std_contact_fraction": float(np.std(contact_fractions)),
            "n_frames": len(contact_fractions),
            "replicate": replicate,
        }

    def aggregate(self, ctx: AggregateContext, results: Sequence[Any]) -> dict[str, Any]:
        """Average contact fraction across replicates with SEM."""
        values = [r["mean_contact_fraction"] for r in results]
        return {
            "mean_contact_fraction": float(np.mean(values)),
            "sem_contact_fraction": float(np.std(values, ddof=1) / np.sqrt(len(values))),
            "replicate_values": values,
        }

    def extract_metrics(self, summary: Any) -> dict[str, MetricValue]:
        """Expose contact fraction as scalar metric for default comparison."""
        return {
            "mean_contact_fraction": MetricValue(
                name="mean_contact_fraction",
                mean=summary["mean_contact_fraction"],
                sem=summary["sem_contact_fraction"],
                replicate_values=summary["replicate_values"],
                higher_is_better=False,
                direction_labels=("more exposed", "unchanged", "more shielded"),
            )
        }

    def plot(self, ctx: PlotContext) -> list[Path]:
        """Generate a contact-fraction bar chart comparing conditions."""
        import matplotlib.pyplot as plt

        from polyzymd.analyses.shared import apply_axis_style, get_colors, get_theme, save_figure

        # Use the base class helper to gather per-condition paths
        data, labels = self._build_plot_data(ctx)
        if not labels:
            return []

        # Load aggregated results for each condition
        plot_labels = []
        means = []
        sems = []
        for label in labels:
            if label not in data or label == "__meta__":
                continue
            agg_dir = data[label]["aggregated_dir"]
            summary = self._load_aggregated_result(agg_dir)
            if summary is None:
                continue
            plot_labels.append(label)
            means.append(summary["mean_contact_fraction"])
            sems.append(summary["sem_contact_fraction"])

        if not plot_labels:
            return []

        # Use shared plotting utilities for consistent theming
        plot_settings = ctx.plot_settings
        theme = get_theme(plot_settings)
        colors = get_colors(len(plot_labels), plot_settings)

        fig, ax = plt.subplots(figsize=(8, 5))
        ax.bar(plot_labels, means, yerr=sems, capsize=5, color=colors[: len(plot_labels)])
        apply_axis_style(
            ax,
            plot_settings,
            title="Solvent Contact Fraction Across Conditions",
            ylabel="Contact fraction",
        )

        output_path = ctx.output_dir / "solvent_contacts_comparison.png"
        output_path.parent.mkdir(parents=True, exist_ok=True)
        save_figure(fig, output_path, plot_settings)
        plt.close(fig)
        return [output_path]

    def format(self, result: Any, output_format: str = "text") -> str:
        """Format comparison result for CLI display."""
        from polyzymd.analyses.base import ComparisonResult

        if isinstance(result, ComparisonResult):
            return format_scalar_comparison(
                result,
                title="Solvent Contacts Comparison",
                metric_label="Mean contact fraction",
                metric_unit="",
                metric_key="mean_contact_fraction",
                output_format=output_format,
                higher_is_better=False,
            )
        return super().format(result, output_format)
```

That file is now discoverable, runnable through compare workflows, and
compatible with default framework formatting/statistics.

## Base Class Helpers

The `Analysis` base class provides three utility methods that reduce
boilerplate across the lifecycle. The `SolventContactsAnalysis` example above uses
`_build_plot_data()` in its `plot()` method, but all three are available
to every plugin.

### `_build_plot_data(ctx, *, include_replicates=False)`

Collects per-condition paths from a `PlotContext` into a dict ready for
plotter functions. Returns `(data, labels)` where `data[label]` contains
`"analysis_dir"` and `"aggregated_dir"` keys:

```python
def plot(self, ctx: PlotContext) -> list[Path]:
    # One line replaces the 8-12 line loop through ctx.conditions
    data, labels = self._build_plot_data(ctx, include_replicates=True)
    if not labels:
        return []

    for label in labels:
        agg_dir = data[label]["aggregated_dir"]
        summary = self._load_aggregated_result(agg_dir)
        # ... plot summary data ...
```

```{admonition} Manual alternative
:class: dropdown

Without `_build_plot_data()`, you write the loop yourself — which is
fine for understanding the data flow:

~~~python
def plot(self, ctx: PlotContext) -> list[Path]:
    for cond in ctx.conditions:
        agg_dir = ctx.analysis_dirs[cond.label] / "aggregated"
        summary = self._load_aggregated_result(agg_dir)
        if summary is None:
            continue
        # ... plot summary data ...
~~~
```

### `_check_cache(result_cls, cache_path, *, recompute, sim_config=None)`

Consolidates the cache-checking pattern used in `compute_replicate()`.
Loads a Pydantic result from disk if the cache file exists and
`recompute` is False, optionally validating the config hash:

```python
def compute_replicate(self, ctx: ReplicateContext, replicate: int) -> dict:
    result_file = ctx.output_dir / "solvent_contacts_result.json"

    # One line replaces the 5-line cache check pattern
    cached = self._check_cache(
        SolventContactsReplicateResult, result_file,
        recompute=ctx.recompute, sim_config=ctx.sim_config,
    )
    if cached is not None:
        return cached

    # ... expensive computation ...
```

```{note}
`_check_cache()` requires a Pydantic model class with a `.load(path)`
method (i.e., a `BaseAnalysisResult` subclass). Dict-based plugins that
want caching should use explicit `json.loads()` / `Path.exists()` checks
instead.
```

### `_format_replicate_range(replicates)`

Formats replicate tuples into compact strings for log messages and
filenames. Contiguous ranges are collapsed:

```python
cls._format_replicate_range((1, 2, 3))    # → "reps1-3"
cls._format_replicate_range((1, 3, 5))    # → "reps1_3_5"
```

Useful in `aggregate()` for logging or naming output files:

```python
def aggregate(self, ctx: AggregateContext, results):
    rep_str = self._format_replicate_range(ctx.replicates)
    logger.info("Aggregating %s across %s", self.name, rep_str)
    # ...
```

### Summary of helpers

| Helper | Used in | Purpose |
|--------|---------|---------|
| `_build_plot_data()` | `plot()` | Collect per-condition paths from `PlotContext` |
| `_check_cache()` | `compute_replicate()` | Load cached Pydantic result if valid |
| `_format_replicate_range()` | `aggregate()`, logging | Compact replicate string for messages/filenames |

## Import Rules

The import structure in the example above is intentional. Getting imports
wrong causes subtle test failures, so follow these rules:

**Module-level imports** (top of file):

- Standard library (`logging`, `pathlib`)
- NumPy — imported at module level because it is available in project
  environments
- Framework utilities from `analyses/shared/` (`TrajectoryLoader`,
  `parse_time_string`, `convert_time`, `AlignmentConfig`, etc.)
- Framework base classes from `analyses/base` and `analyses/stats`

**Lazy imports inside methods** (inside `compute_replicate`, `plot`, etc.):

- `MDAnalysis` / `mdtraj` — heavy simulation libraries
- `matplotlib` / `seaborn` — plotting libraries
- Any package that may not be installed in all environments

```python
# At module level — always available, needed for patch targets in tests
from polyzymd.analyses.shared import TrajectoryLoader, convert_time, parse_time_string


# Inside methods — heavy deps that may not be installed
def compute_replicate(self, ctx, replicate):
    import MDAnalysis as mda
    ...


def plot(self, ctx):
    import matplotlib.pyplot as plt
    ...
```

**Why this matters for testing:** when you mock `TrajectoryLoader` in tests,
you write `@patch("polyzymd.analyses.solvent_contacts.TrajectoryLoader")`. This replaces
the module-level name. If `TrajectoryLoader` were imported lazily inside the
method, the patch target would be different and easier to get wrong.

## Context Objects Reference

The framework provides context objects so plugins never need to load configs
or discover paths themselves:

| Context | Passed To | Key Attributes |
|---------|-----------|----------------|
| `ReplicateContext` | `compute_replicate()` | `.sim_config`, `.settings`, `.output_dir`, `.replicate`, `.recompute`, `.equilibration`, `.result_path` |
| `AggregateContext` | `aggregate()` | `.condition`, `.replicates`, `.output_dir`, `.settings`, `.equilibration`, `.result_path` |
| `ComparisonContext` | `compare()` | `.name`, `.conditions`, `.analysis_dirs`, `.results_dir`, `.effective_control`, `.control_label`, `.settings`, `.equilibration`, `.recompute`, `.fdr_alpha`, `.result_path`, `.aggregated_results` |
| `PlotContext` | `plot()` | `.conditions`, `.analysis_dirs`, `.output_dir`, `.results_dir`, `.settings`, `.plot_settings`, `.comparison_path`, `.control_label` |

(plotcontext-plot-settings)=
## `PlotContext.plot_settings`

`PlotContext.plot_settings` carries global plotting configuration from
`polyzymd.compare.config.PlotSettings`. Use this object instead of hard-coding
figure style, DPI, or output behavior.

### Type and global fields

`ctx.plot_settings` is always an instance of `PlotSettings` (the orchestrator
guarantees this — it is never `None`).

Global fields available on `PlotSettings` include:

- `output_dir`
- `format`
- `dpi`
- `style`
- `color_palette`
- `theme`

### Theme system

`PlotSettings.theme` uses a `PlotTheme` model that includes 28 visual style
fields (for example: `title_fontsize`, `label_fontsize`, `tick_fontsize`,
`dot_size`, line widths, transparency, and grid controls). The theme system has
three presets:

- `PlotTheme.publication()`
- `PlotTheme.presentation()`
- `PlotTheme.minimal()`

### Per-analysis settings

Per-analysis settings are accessed as attributes on `PlotSettings`, for example:

- `plot_settings.contacts`
- `plot_settings.rmsf`

These fields use a `__getattr__` fallback pattern, so access is stable even
when a specific analysis block is not configured; defaults are returned.

### Standard plugin pattern

The orchestrator guarantees that `ctx.plot_settings` is always a valid
`PlotSettings` instance — never `None`. Use it directly:

```python
plot_settings = ctx.plot_settings
```

### Example using shared plotting helpers

```python
from pathlib import Path

from polyzymd.analyses.base import PlotContext
from polyzymd.analyses.shared import (
    apply_axis_style,
    get_colors,
    get_theme,
    save_figure,
)


def plot(self, ctx: PlotContext) -> list[Path]:
    import matplotlib.pyplot as plt

    plot_settings = ctx.plot_settings
    theme = get_theme(plot_settings)
    colors = get_colors(len(ctx.conditions), plot_settings)

    fig, ax = plt.subplots(figsize=(8, 5))
    _ = theme
    _ = colors
    # ... plotting code ...
    apply_axis_style(ax, plot_settings, title="My Plot", ylabel="Value (Å)")
    out = ctx.output_dir / "my_plot.png"
    save_figure(fig, out, plot_settings)
    plt.close(fig)
    return [out]
```

```{tip}
Prefer shared plotting helpers over custom style code in each plugin. This
keeps visual output consistent across analyses and makes global theming work.
```

## Choosing Your Comparison Path

### Path A: Default Scalar Comparison (recommended)

Implement `extract_metrics()` and you get t-tests, Cohen's d, ANOVA, and
ranking for free.

The framework loads aggregated results from disk automatically using:

- `json.loads()` for dict-based plugins
- your `AggregatedResultClass.load(path)` when `AggregatedResultClass` is set

You generally do **not** need to override `_deserialize_result()`.

### Path B: Custom Comparison

Override `compare()` for multi-metric or entry-table analyses:

```python
from typing import Any

from polyzymd.analyses.base import ComparisonContext


def compare(self, ctx: ComparisonContext) -> Any:
    # Load results, compute custom statistics, return your result model
    # The returned object must support save(path)
    ...
```

See `src/polyzymd/analyses/contacts/` or `src/polyzymd/analyses/distances/`
for full examples.

## Loading Results in `plot()`

`PlotContext` provides paths and conditions but does **not** carry pre-loaded
aggregated results. The recommended approach is to use `_build_plot_data()` to
collect per-condition paths, then `_load_aggregated_result()` to load each
condition's data:

```python
from pathlib import Path

from polyzymd.analyses.base import PlotContext


def plot(self, ctx: PlotContext) -> list[Path]:
    import matplotlib.pyplot as plt

    data, labels = self._build_plot_data(ctx)
    for label in labels:
        if label not in data or label == "__meta__":
            continue
        agg_dir = data[label]["aggregated_dir"]
        summary = self._load_aggregated_result(agg_dir)
        if summary is not None:
            # summary is a dict or AggregatedResultClass instance
            # ... plot data from summary ...
            pass

    return []
```

You can also loop through `ctx.conditions` directly for fine-grained control:

```python
def plot(self, ctx: PlotContext) -> list[Path]:
    import matplotlib.pyplot as plt

    for cond in ctx.conditions:
        agg_dir = ctx.analysis_dirs[cond.label] / "aggregated"
        summary = self._load_aggregated_result(agg_dir)
        if summary is not None:
            # ... plot data from summary ...
            pass

    return []
```

`_load_aggregated_result()` looks for `result.json` in the aggregated
directory, falling back to the most recent `*.json` file. It deserializes
using `AggregatedResultClass` if set, otherwise `json.loads()`. Returns `None`
if no result file is found.

```{warning}
Do not bypass `_load_aggregated_result()` with custom file parsing unless your
plugin requires a non-standard storage format. The helper already handles both
dict and model-based deserialization behavior.
```

## Formatting with `format_scalar_comparison()`

For analyses using the default comparison path, `format_scalar_comparison()`
from `analyses.stats` produces formatted tables with rankings, effect sizes,
and significance stars:

```python
from typing import Any

from polyzymd.analyses.base import ComparisonResult
from polyzymd.analyses.stats import format_scalar_comparison


def format(self, result: Any, output_format: str = "text") -> str:
    if isinstance(result, ComparisonResult):
        return format_scalar_comparison(
            result,
            title="My Analysis Comparison",
            metric_label="My Metric",
            metric_unit="Å",
            metric_key="my_metric",
            output_format=output_format,
            higher_is_better=True,
        )
    return super().format(result, output_format)
```

**Parameters:**

| Parameter | Type | Description |
|-----------|------|-------------|
| `result` | `ComparisonResult` | The result to format |
| `title` | `str` | Display title (for example, `"RMSF Comparison"`) |
| `metric_label` | `str` | Human-readable metric name for table headers |
| `metric_unit` | `str` | Unit string appended to values (for example, `"Å"`) |
| `metric_key` | `str \| None` | Key prefix in `ConditionSummary` extra fields; auto-detected if `None` |
| `output_format` | `str` | `"text"`, `"markdown"`, or `"json"` |
| `higher_is_better` | `bool` | Affects interpretation wording |

(return-types)=
## Return Types: Dicts vs Pydantic Models

Plugins can return either plain dicts or typed Pydantic models from
`compute_replicate()` and `aggregate()`. Both paths are fully supported.

**For new plugins, start with dicts.** They require no extra imports, no extra
classes, and they work with the default comparison pipeline out of the box.
You can always migrate to Pydantic models later if your plugin grows in
complexity.

### Dict path (recommended starting point)

Dict plugins are the simplest way to get a working analysis:

- `compute_replicate()` returns `dict`
- `aggregate()` returns `dict`
- No `AggregatedResultClass` class variable needed
- Framework saves/loads JSON using standard `json.dumps()` / `json.loads()`

```python
class MyAnalysis(Analysis):
    name = "my_analysis"

    class Settings(BaseModel):
        selection: str = "protein and name CA"

    def compute_replicate(self, ctx, replicate):
        return {"value": 1.23, "replicate": replicate}

    def aggregate(self, ctx, results):
        vals = [r["value"] for r in results]
        return {"mean_value": sum(vals) / len(vals), "replicate_values": vals}

    def extract_metrics(self, summary):
        return {
            "value": MetricValue(
                name="value",
                mean=summary["mean_value"],
                sem=0.1,
                replicate_values=summary["replicate_values"],
            )
        }
```

This is the path used in the complete `SolventContactsAnalysis` example earlier in this
tutorial. The framework handles serialization, deserialization, and the full
comparison lifecycle without any Pydantic result models.

### What is Pydantic, and when would you use it?

[Pydantic](https://docs.pydantic.dev/) is a Python data validation library.
A Pydantic `BaseModel` subclass is like a dataclass with automatic type
checking: if you declare a field as `float`, Pydantic raises a validation
error if someone passes a string. PolyzyMD already uses Pydantic for its
config system (`Settings` classes), so it is always available in the
environment.

For **result models**, Pydantic adds:

- **Validation on construction** — catches type errors early (e.g., passing a
  string where a float is expected)
- **Typed attribute access** — `result.mean_rg` instead of
  `result["mean_rg"]`, with IDE autocomplete
- **Structured extension path** — start with scaffold-generated `BaseModel`
  classes, then upgrade to `BaseAnalysisResult` if you need helper I/O methods
  or NPZ sidecar storage
- **Nested structure** — complex result hierarchies with validated sub-models

### Pydantic model path

Use models when you need validation, nested structure, strict typing, or more
explicit result contracts.

Key pieces:

1. The scaffold `--style pydantic` path generates result models inheriting
   from Pydantic `BaseModel`
2. Plugin class sets `AggregatedResultClass` to the aggregated model class
3. Framework deserializes aggregated files through
   `AggregatedResultClass.model_validate_json(...)` (or `.load(path)` when the
   class provides it)
4. `BaseModel` classes are the simple default and do not include `.save()` /
   `.load()` helpers
5. If you need NPZ sidecars or model-level I/O helpers, manually upgrade your
   result models to
   `polyzymd.analyses._results_base.BaseAnalysisResult`

```python
from pathlib import Path
from typing import ClassVar

from pydantic import BaseModel

from polyzymd.analyses.base import Analysis


class ReplicateResult(BaseModel):
    mean_rg: float
    replicate: int


class AggregatedResult(BaseModel):
    mean_rg: float
    sem_rg: float
    replicate_values: list[float]


class SolventContactsAnalysis(Analysis):
    name: ClassVar[str] = "solvent_contacts"
    AggregatedResultClass = AggregatedResult

    class Settings(BaseModel):
        selection: str = "protein and name CA"

    def compute_replicate(self, ctx, replicate) -> ReplicateResult:
        return ReplicateResult(mean_rg=15.0, replicate=replicate)

    def aggregate(self, ctx, results) -> AggregatedResult:
        values = [r.mean_rg for r in results]
        return AggregatedResult(
            mean_rg=sum(values) / len(values),
            sem_rg=0.3,
            replicate_values=values,
        )
```

If you later need NPZ sidecar storage or `.save()` / `.load()` helper methods,
replace `BaseModel` with `BaseAnalysisResult` in those result classes.

```{important}
**Gotcha: `AggregatedResultClass` requires model returns**

When you set `AggregatedResultClass`, your `compute_replicate()` and
`aggregate()` should return instances of the corresponding Pydantic result
models. Do not return plain dicts in this mode.

Why: `_load_aggregated_result()` / `_deserialize_result()` uses the declared
class for deserialization. If on-disk JSON represents dict-shaped data that does
not match the model contract expected by `AggregatedResultClass.load(path)`,
deserialization fails.
```

### Comparison: dicts vs Pydantic models

| Aspect | Dicts | Pydantic models |
|--------|-------|-----------------|
| **Setup effort** | None — just return `{}` | Define result classes (scaffold uses `BaseModel` by default) |
| **Validation** | None — typos in keys are silent | Automatic type checking on construction |
| **Attribute access** | `result["mean_rg"]` | `result.mean_rg` with IDE autocomplete |
| **Large array support** | Manual — save NumPy arrays yourself | BaseModel: manual; BaseAnalysisResult: NPZ sidecar helpers |
| **Nested structure** | Manual — dicts of dicts | Validated sub-models with typed fields |
| **Serialization** | `json.dumps()` / `json.loads()` | BaseModel JSON by framework; optional `BaseAnalysisResult.save()` / `.load()` |
| **Framework wiring** | No `AggregatedResultClass` needed | Set `AggregatedResultClass` on plugin class |
| **Best for** | Scalar analyses, quick prototypes, new contributors | Complex results, per-residue arrays, production plugins |

**Rule of thumb:** Start with dicts. If you find yourself writing defensive
`result.get("key", default)` checks, duplicating key strings, or needing to
store large NumPy arrays alongside JSON, migrate to Pydantic models.

## A Note on the `compare/` Package

You may notice a `compare/` package with `compare/results/` and
`compare/settings.py`. Existing plugins reference these because they were
designed in an earlier architecture stage.

**You do NOT need to create files in `compare/` for a new plugin.** Keep your
plotting logic in your plugin's `plot()` method and your formatting in
`format()`. Start with all logic in `__init__.py`; as your plugin grows, extract
plotting functions into a `_plotters.py` module within the package (see
{ref}`plotters-extraction` below).

(shared-utilities)=
## Shared Utilities (`analyses/shared/`)

The `analyses/shared/` package provides reusable infrastructure for plugin
authors. The example above uses `TrajectoryLoader`, `parse_time_string`, and
`convert_time`. The package now also re-exports plotting and config-hash helper
symbols for direct plugin use.

### TrajectoryLoader

The most important shared utility. Handles topology discovery, daisy-chain
trajectory segments, timestep detection, and equilibration frame skipping:

```python
from polyzymd.analyses.shared import TrajectoryLoader, convert_time, parse_time_string

loader = TrajectoryLoader(ctx.sim_config)
u = loader.load_universe(replicate)

# Convert equilibration time to frame offset
eq_value, eq_unit = parse_time_string(ctx.equilibration)
timestep_ps = loader.get_timestep(replicate, unit="ps")
eq_time_ps = convert_time(eq_value, eq_unit, "ps")
start_frame = int(eq_time_ps / timestep_ps)
```

### Available re-exports

The following symbols are re-exported from `polyzymd.analyses.shared` for
convenient one-line imports:

| Source module | Re-exported symbols |
|---------------|---------------------|
| `shared.loader` | `TrajectoryLoader`, `TrajectoryInfo`, `parse_time_string`, `convert_time`, `time_to_frame` |
| `shared.alignment` | `AlignmentConfig`, `ReferenceMode`, `align_trajectory`, `get_alignment_description` |
| `shared.statistics` | `compute_sem`, `aggregate_per_residue_stats`, `aggregate_region_stats`, `weighted_mean_with_sem`, `StatResult`, `PerResidueStats` |
| `shared.autocorrelation` | `compute_acf`, `estimate_correlation_time`, `statistical_inefficiency`, `statistical_inefficiency_multiple`, `n_effective`, `get_independent_indices`, `check_statistical_reliability`, `ACFResult`, `CorrelationTimeResult` |
| `shared.pbc` | `minimum_image_distance`, `pairwise_distances_pbc` |
| `shared.constants` | `DEFAULT_CONTACT_CUTOFF`, `DEFAULT_DISTANCE_THRESHOLD`, `DEFAULT_SURFACE_EXPOSURE_THRESHOLD` |
| `shared.defaults` | `AnalysisDefaults` |
| `shared.config_hash` | `compute_config_hash`, `validate_config_hash` |
| `shared.plotting` | `get_theme`, `apply_axis_style`, `apply_legend`, `get_colors`, `get_output_path`, `save_figure`, `grouped_bars` |

#### Convergence detection

`analyses/shared/convergence.py` provides `find_convergence_time()`, a
sliding-window slope heuristic for detecting sustained convergence in any
timeseries. The RMSD plugin uses it to annotate per-run convergence
diagnostics, but it is generic enough for any signal (e.g., Rg, energy).

```python
from polyzymd.analyses.shared.convergence import find_convergence_time

result = find_convergence_time(time_ns, rmsd_values)
if result.converged:
    print(f"Converged at {result.convergence_time_ns:.1f} ns")
```

Other submodules (`selections`, `aggregation`, `diagnostics`, `aa_classification`,
`logging_utils`) are available via direct import but are not re-exported from
the package root:

The `metric_type` classification is a planned future enhancement and is not
yet implemented in the current codebase.

```python
# Direct submodule import (not re-exported)
from polyzymd.analyses.shared.selections import get_positions
from polyzymd.analyses.shared.aggregation import collect_replicate_results
```

Plugins can now import plotting and config hash utilities directly from the
package root:

```python
from polyzymd.analyses.shared import (
    compute_config_hash,
    get_output_path,
    save_figure,
)
```

```{note}
Re-export imports are preferred for plugin code because they keep import paths
stable even if internal module locations change.
```

## Contributor Checklist

When creating a new analysis plugin:

- [ ] File: `src/polyzymd/analyses/<name>/` (use `polyzymd new-analysis <name>` to scaffold; add `--style pydantic` for typed result models)
- [ ] `name: ClassVar[str]` — unique, lowercase, used in CLI
- [ ] `Settings` class — Pydantic v2 BaseModel with defaults
- [ ] `compute_replicate()` — uses `TrajectoryLoader`, returns dict or Pydantic model
- [ ] `aggregate()` — returns dict or Pydantic model; framework auto-saves to `ctx.result_path`
- [ ] Choose comparison path:
  - [ ] Default: implement `extract_metrics()`
  - [ ] Custom: override `compare()` returning a model with `.save()`
- [ ] `plot()` — use `_build_plot_data()` + shared plotting helpers (`get_theme`, `apply_axis_style`, `save_figure`, `get_colors`)
- [ ] (Optional) `format()` — CLI display via `format_scalar_comparison()`
- [ ] (Optional) Use `_check_cache()` in `compute_replicate()` for result caching
- [ ] (Optional) Use `_format_replicate_range()` for compact log messages
- [ ] If using `AggregatedResultClass`, return matching result model instances
- [ ] Tests: `tests/test_<name>_plugin.py` (scaffold generates plot tests too)
- [ ] Verify: `pixi run -e build pytest tests/test_<name>_plugin.py -v`

## What a Good Community PR Looks Like

A strong contribution usually includes:

- one new plugin package in `src/polyzymd/analyses/<name>/`
- one focused test file
- one short docs update showing the `plugins:` config block
- figures only if they add scientific value

Aim for a plugin that another MDAnalysis user can understand in one read.

(plotters-extraction)=
## Extracting Plotting to `_plotters.py`

As a plugin grows, its `plot()` method may accumulate many helper functions.
Established plugins extract these into a `_plotters.py` module within the
plugin package:

```text
analyses/rmsf/
├── __init__.py      # Analysis lifecycle (compute, aggregate, compare, plot)
├── _plotters.py     # Plotting functions called by plot()
├── _results.py      # Result models
└── ...
```

The `plot()` method in `__init__.py` delegates to functions in `_plotters.py`:

```python
from ._plotters import plot_comparison_bars, plot_per_residue

def plot(self, ctx: PlotContext) -> list[Path]:
    data, labels = self._build_plot_data(ctx)
    paths = []
    paths.append(plot_comparison_bars(data, labels, ctx.output_dir, ctx.plot_settings))
    paths.append(plot_per_residue(data, labels, ctx.output_dir, ctx.plot_settings))
    return [p for p in paths if p is not None]
```

**When to extract:** Consider extracting when your plugin has 3+ plot
functions or `__init__.py` exceeds ~500 lines. For simple plugins (like
`catalytic_triad/`), keeping plotting inline is fine.

**All 6 established plugins have `_plotters.py`:** rmsf, secondary_structure,
distances, contacts, catalytic_triad, hydrogen_bonds.

## Anti-Patterns to Avoid

### Do not bypass the context

```python
# Wrong
def compute_replicate(self, ctx, replicate):
    config = SimulationConfig.from_yaml("/my/path/config.yaml")
```

```python
# Right
def compute_replicate(self, ctx, replicate):
    sim_config = ctx.sim_config  # already loaded by framework
```

### Understand the result-saving convention

The framework auto-saves return values from `compute_replicate()` and
`aggregate()` to `ctx.result_path`. For dict plugins, returning your result is
usually sufficient.

For caching with explicit target naming, save explicitly:

```python
import json


def aggregate(self, ctx, results):
    agg = {"mean": 15.0}
    target_path = ctx.result_path
    if target_path is None:
        target_path = ctx.output_dir / "result.json"
    target_path.parent.mkdir(parents=True, exist_ok=True)
    target_path.write_text(json.dumps(agg, indent=2), encoding="utf-8")
    return agg
```

### Do not import heavy deps at module level

```python
# Wrong — breaks environments without MDAnalysis
import MDAnalysis as mda


# Right — lazy import inside method
def compute_replicate(self, ctx, replicate):
    import MDAnalysis as mda
```

Framework utilities from `analyses/shared/` are imported at module level —
they are always available and patchable in tests. Only heavy third-party
dependencies (MDAnalysis, mdtraj, matplotlib) should be imported inside
methods.

### Do not create files in `compare/`

New plugins should keep all logic in the plugin package. The `compare/results/`
directory is primarily used by existing plugins designed under older
architecture constraints.

## Existing Plugins to Study

| Plugin | Comparison path | Good example of |
|--------|------------------|-----------------|
| `catalytic_triad/` | Default (`extract_metrics`) | Simplest lifecycle — minimal compute logic with `_plotters.py` for KDE/bar plots |
| `rmsf/` | Default (`extract_metrics`) | Default scalar comparison with `_plotters.py` extraction |
| `secondary_structure/` | Default (`extract_metrics`) | Wrapping existing calculators with `_plotters.py` and typed result classes |
| `contacts/` | Custom (`compare` override) | Multi-metric comparison, condition filtering, large `_plotters.py` |
| `distances/` | Custom (`compare` override) | Entry-table style comparison with `_plotters.py` extraction |

```{important}
Do not interpret plugin differences as "simple" versus "complex" quality.
Several existing plugins reflect earlier architectural phases. Use them as
pattern references for specific tasks (result modeling, comparison style,
plotting), not as maturity rankings.
```

## Testing Your Plugin

### Test file naming

Name your test file `tests/test_<name>_plugin.py` (for example,
`tests/test_solvent_contacts_plugin.py`).

### Testing strategy

A reliable test suite for plugins usually combines:

- discovery and class-contract tests (`name`, `Settings`, discovery)
- pure logic tests (`extract_metrics`, `aggregate`)
- isolated compute tests (`compute_replicate` with mocked trajectory loader)
- filesystem-backed integration tests for `compare()` and `plot()`

This layering keeps tests fast while still validating framework integration.

### Basic tests (discovery, settings, metrics)

These tests use no trajectory data and verify the plugin integrates with the
framework:

```python
from polyzymd.analyses import get_analysis, list_analyses
from polyzymd.analyses.base import MetricValue


class TestSolventContactsDiscovery:
    """Test that the plugin is discoverable."""

    def test_discovered(self):
        assert "solvent_contacts" in list_analyses()

    def test_class_variables(self):
        cls = get_analysis("solvent_contacts")
        assert cls.name == "solvent_contacts"
        assert hasattr(cls, "Settings")

    def test_settings_defaults(self):
        cls = get_analysis("solvent_contacts")
        settings = cls.Settings()
        assert settings.selection == "protein and name CA"


class TestSolventContactsMetrics:
    """Test metric extraction for default comparison."""

    def test_extract_metrics(self):
        cls = get_analysis("solvent_contacts")
        analysis = cls()
        summary = {
            "mean_contact_fraction": 0.142,
            "sem_contact_fraction": 0.010,
            "replicate_values": [0.136, 0.148],
        }
        metrics = analysis.extract_metrics(summary)
        assert "mean_contact_fraction" in metrics
        assert isinstance(metrics["mean_contact_fraction"], MetricValue)
        assert metrics["mean_contact_fraction"].higher_is_better is False
```

### Testing `compute_replicate()` with mocked trajectories

`compute_replicate()` requires MDAnalysis and trajectory files. Mock
`TrajectoryLoader` at the module level to test compute logic without real
trajectory data:

```python
from pathlib import Path
from unittest.mock import MagicMock, patch

from polyzymd.analyses import get_analysis
from polyzymd.analyses.base import Condition, ReplicateContext


class TestSolventContactsCompute:
    """Test compute_replicate with mocked trajectories."""

    @patch("polyzymd.analyses.solvent_contacts.TrajectoryLoader")
    def test_computes_metric(self, mock_loader_cls, tmp_path):
        cls = get_analysis("solvent_contacts")
        analysis = cls()
        settings = cls.Settings()

        mock_loader = MagicMock()
        mock_loader_cls.return_value = mock_loader

        mock_universe = MagicMock()
        mock_target = MagicMock()
        mock_target.__len__.return_value = 10
        mock_target.positions = [[0.0, 0.0, 0.0]] * 10
        mock_solvent = MagicMock()
        mock_solvent.__len__.return_value = 20
        mock_solvent.positions = [[1.0, 1.0, 1.0]] * 20
        mock_universe.select_atoms.side_effect = [mock_target, mock_solvent]

        mock_trajectory = MagicMock()
        mock_trajectory.__getitem__.return_value = range(50)
        mock_universe.trajectory = mock_trajectory

        mock_loader.load_universe.return_value = mock_universe
        mock_loader.get_timestep.return_value = 10.0

        condition = Condition(
            label="Test",
            config_path=Path("/fake/config.yaml"),
            replicates=(1,),
            sim_config=MagicMock(),
        )
        ctx = ReplicateContext(
            condition=condition,
            replicate=1,
            sim_config=condition.sim_config,
            output_dir=tmp_path / "run_1",
            equilibration="0ns",
            recompute=True,
            settings=settings,
        )

        result = analysis.compute_replicate(ctx, replicate=1)
        assert "mean_contact_fraction" in result
        assert result["mean_contact_fraction"] >= 0
        assert result["replicate"] == 1
```

Note the patch target:
`@patch("polyzymd.analyses.solvent_contacts.TrajectoryLoader")`. This works
because `TrajectoryLoader` is imported at module level in
`solvent_contacts/__init__.py`.

### Testing `aggregate()`

Aggregation takes plain data and needs no trajectory mocks:

```python
from pathlib import Path
from unittest.mock import MagicMock

import numpy as np

from polyzymd.analyses import get_analysis
from polyzymd.analyses.base import AggregateContext, Condition


class TestSolventContactsAggregate:
    """Test aggregation across replicates."""

    def test_aggregate(self, tmp_path):
        cls = get_analysis("solvent_contacts")
        analysis = cls()
        settings = cls.Settings()

        condition = Condition(
            label="Test",
            config_path=Path("/fake/config.yaml"),
            replicates=(1, 2, 3),
            sim_config=MagicMock(),
        )
        ctx = AggregateContext(
            condition=condition,
            replicates=(1, 2, 3),
            output_dir=tmp_path / "aggregated",
            equilibration="10ns",
            settings=settings,
        )

        results = [
            {"mean_contact_fraction": 0.12, "replicate": 1},
            {"mean_contact_fraction": 0.10, "replicate": 2},
            {"mean_contact_fraction": 0.11, "replicate": 3},
        ]

        agg = analysis.aggregate(ctx, results)
        expected_mean = np.mean([0.12, 0.10, 0.11])
        assert abs(agg["mean_contact_fraction"] - expected_mean) < 1e-10
        assert agg["replicate_values"] == [0.12, 0.10, 0.11]
```

### Testing `compare()` with filesystem-backed aggregated results

`compare()` reads condition-specific result files from analysis directories, so
tests should construct realistic filesystem layout with `result.json` files and
a `ComparisonContext`:

```python
import json
from pathlib import Path
from unittest.mock import MagicMock

from polyzymd.analyses import get_analysis
from polyzymd.analyses.base import ComparisonContext, Condition


class TestCompare:
    def test_compare_two_conditions(self, tmp_path):
        cls = get_analysis("solvent_contacts")
        analysis = cls()

        # Set up filesystem with aggregated results
        for label, value in [("WT", 0.12), ("Mutant", 0.18)]:
            agg_dir = tmp_path / label / "solvent_contacts" / "aggregated"
            agg_dir.mkdir(parents=True)
            (agg_dir / "result.json").write_text(
                json.dumps(
                    {
                        "mean_contact_fraction": value,
                        "sem_contact_fraction": 0.01,
                        "replicate_values": [value - 0.01, value + 0.01],
                    }
                ),
                encoding="utf-8",
            )

        conditions = [
            Condition(
                label="WT",
                config_path=Path("/fake"),
                replicates=(1, 2),
                sim_config=MagicMock(),
            ),
            Condition(
                label="Mutant",
                config_path=Path("/fake"),
                replicates=(1, 2),
                sim_config=MagicMock(),
            ),
        ]

        ctx = ComparisonContext(
            name="test",
            conditions=conditions,
            excluded_conditions=[],
            control_label="WT",
            analysis_dirs={
                "WT": tmp_path / "WT" / "solvent_contacts",
                "Mutant": tmp_path / "Mutant" / "solvent_contacts",
            },
            results_dir=tmp_path / "results",
            equilibration="10ns",
            settings=cls.Settings(),
            recompute=True,
        )

        result = analysis.compare(ctx)
        assert result is not None
```

```{tip}
This style of test is robust because it validates the same on-disk contract
used by real compare runs (`<condition>/<analysis>/aggregated/result.json`).
```

### Matplotlib backend in CI/headless tests

When testing plot code in CI or headless environments, force a non-interactive
backend before importing `matplotlib.pyplot`:

```python
import matplotlib

matplotlib.use("Agg")  # headless backend for CI / testing
import matplotlib.pyplot as plt
```

This avoids backend/display errors and keeps plotting tests deterministic.

Run tests with:

```bash
pixi run -e build pytest tests/test_solvent_contacts_plugin.py -v
```

## See Also

- {doc}`../explanation/architecture` — Overall system architecture
- {doc}`../explanation/analysis_statistics_best_practices` — Autocorrelation and uncertainty
- {doc}`../how_to/analysis_compare_conditions` — User guide for running comparisons
- API reference: `polyzymd.analyses.base.Analysis`

(comparison-yaml)=
## Appendix: The `comparison.yaml` File

The `comparison.yaml` file is the single configuration file that drives the
entire comparison workflow. It defines which simulation conditions to compare,
which analyses to run, and how to generate plots. Understanding this file is
important for plugin authors because your `Settings` class maps directly to
the `plugins:` section.

### Generating a template

The fastest way to create a `comparison.yaml` is with the scaffold command:

```bash
polyzymd compare init -n my_study
```

This creates a project directory with a fully commented template. You can
also specify a custom equilibration time:

```bash
polyzymd compare init -n my_study --eq-time 20ns
```

### File structure

A `comparison.yaml` has four main sections:

```yaml
# Project metadata
name: "my_study"
description: "Comparison of polymer conditions"
control: "Condition A"  # or null if no control

# Conditions — each points to a simulation's config.yaml
conditions:
  - label: "Condition A"
    config: "../path/to/condition_a/config.yaml"
    replicates: [1, 2, 3]

  - label: "Condition B"
    config: "../path/to/condition_b/config.yaml"
    replicates: [1, 2, 3]

# Defaults
defaults:
  equilibration_time: "10ns"

# Plugins — which analyses to run and their settings
plugins:
  rmsf:
    selection: "protein and name CA"
    reference_mode: "centroid"

  # Your plugin's Settings fields go here:
  my_analysis:
    example_parameter: "value"

# Plot settings — visual customization
plot_settings:
  output_dir: "figures/"
  format: "png"
  dpi: 300
  style: "publication"
```

### How `plugins:` connects to your `Settings` class

When the framework runs your plugin, it reads the `plugins:` section of
`comparison.yaml` and deserializes the matching key into your `Settings`
model. For example, if your plugin has:

```python
class MySettings(BaseModel):
    cutoff: float = 4.5
    selection: str = "protein and name CA"
```

Then users configure it in `comparison.yaml` as:

```yaml
plugins:
  my_analysis:
    cutoff: 5.0
    selection: "chainID A"
```

Any fields not specified in the YAML use the defaults from your `Settings`
class. This is why providing sensible defaults for every field is important —
it means the analysis works without any explicit configuration.

### Running comparisons

After editing `comparison.yaml`:

```bash
# Validate the configuration
polyzymd compare validate

# Run a specific analysis
polyzymd compare run rmsf

# Run all enabled analyses
polyzymd compare run-all

# Generate plots
polyzymd compare plot-all
```

See {doc}`../how_to/analysis_compare_conditions` for a complete user guide on running
comparisons.
