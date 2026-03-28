# Extending the Analysis Framework

This guide shows you how to add a new analysis type to PolyzyMD in one file.
By the end, you will have a working analysis plugin that can compute metrics
per replicate, aggregate across replicates, compare conditions, generate
figures, and display results on the CLI.

```{note}
This guide covers the **plugin system** in ``analyses/``.
The underlying per-condition calculators (private ``_calculator.py`` modules
inside each analysis sub-package in ``analyses/``) are implementation details — you do not need to touch them.
```

## Prerequisites

- A working pixi environment: `pixi install -e build`
- Basic familiarity with Pydantic v2 `BaseModel`
- Understanding of the PolyzyMD chain convention (A=protein, B=substrate,
  C=polymer)

## How the Plugin System Works

Every analysis in PolyzyMD is a single class that subclasses `Analysis` and
lives in `src/polyzymd/analyses/`. The framework discovers plugins
automatically — no registry edits, no imports, no boilerplate. You drop a file,
and the CLI can run it.

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
(`polyzymd compare run rg`), config files, and output directories.

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

Generate matplotlib figures. Use `self._load_aggregated_result(agg_dir)` to
load each condition's data. Return a list of figure paths.

### `format(result, output_format)` (optional)

Custom CLI display. Use `format_scalar_comparison()` from `analyses.stats`
for formatted tables with rankings, effect sizes, and significance stars.

## The Complete Example

Here is a single, complete plugin file that uses real framework utilities.
This is the pattern you should follow for new analyses.

```{important}
The `RgAnalysis` implementation below is a **tutorial illustration** in docs.
Treat it as a reference template for your own plugin file. It is not presented
here as a required built-in deployed analysis.
```

Create `src/polyzymd/analyses/rg.py` in your own branch or downstream project:

```python
"""Radius of gyration analysis plugin."""
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


class RgAnalysis(Analysis):
    """Radius of gyration: compactness of the protein structure."""

    name: ClassVar[str] = "rg"

    class Settings(BaseModel):
        """Settings for radius of gyration analysis.

        Attributes
        ----------
        selection : str
            MDAnalysis selection string for atoms to include in Rg
            calculation. Defaults to C-alpha atoms
        """

        selection: str = Field(
            default="protein and name CA",
            description="MDAnalysis atom selection for Rg calculation",
        )

    def compute_replicate(self, ctx: ReplicateContext, replicate: int) -> dict[str, Any]:
        """Compute mean Rg for one replicate.

        Uses TrajectoryLoader for topology discovery, trajectory segment
        daisy-chaining, and equilibration frame skipping
        """
        # Lazy-import heavy third-party dependency
        import MDAnalysis as mda  # noqa: F401

        loader = TrajectoryLoader(ctx.sim_config)
        u = loader.load_universe(replicate)

        eq_value, eq_unit = parse_time_string(ctx.equilibration)
        timestep_ps = loader.get_timestep(replicate, unit="ps")
        eq_time_ps = convert_time(eq_value, eq_unit, "ps")
        start_frame = int(eq_time_ps / timestep_ps)

        atoms = u.select_atoms(ctx.settings.selection)
        rg_values = []
        for _ts in u.trajectory[start_frame:]:
            rg_values.append(atoms.radius_of_gyration())

        return {
            "mean_rg": float(np.mean(rg_values)),
            "std_rg": float(np.std(rg_values)),
            "n_frames": len(rg_values),
            "replicate": replicate,
        }

    def aggregate(self, ctx: AggregateContext, results: Sequence[Any]) -> dict[str, Any]:
        """Average Rg across replicates with SEM."""
        values = [r["mean_rg"] for r in results]
        return {
            "mean_rg": float(np.mean(values)),
            "sem_rg": float(np.std(values, ddof=1) / np.sqrt(len(values))),
            "replicate_values": values,
        }

    def extract_metrics(self, summary: Any) -> dict[str, MetricValue]:
        """Expose Rg as a scalar metric for automatic t-tests and ANOVA."""
        return {
            "mean_rg": MetricValue(
                name="mean_rg",
                mean=summary["mean_rg"],
                sem=summary["sem_rg"],
                replicate_values=summary["replicate_values"],
                higher_is_better=False,
                direction_labels=("compacting", "unchanged", "expanding"),
            )
        }

    def plot(self, ctx: PlotContext) -> list[Path]:
        """Generate an Rg bar chart comparing conditions."""
        import matplotlib.pyplot as plt

        labels = []
        means = []
        sems = []
        for cond in ctx.conditions:
            agg_dir = ctx.analysis_dirs[cond.label] / "aggregated"
            summary = self._load_aggregated_result(agg_dir)
            if summary is None:
                continue
            labels.append(cond.label)
            means.append(summary["mean_rg"])
            sems.append(summary["sem_rg"])

        if not labels:
            return []

        fig, ax = plt.subplots(figsize=(8, 5))
        ax.bar(labels, means, yerr=sems, capsize=5, color="steelblue")
        ax.set_ylabel("Radius of Gyration (Å)")
        ax.set_title("Rg Comparison Across Conditions")

        output_path = ctx.output_dir / "rg_comparison.png"
        output_path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(output_path, dpi=150, bbox_inches="tight")
        plt.close(fig)
        return [output_path]

    def format(self, result: Any, output_format: str = "text") -> str:
        """Format comparison result for CLI display."""
        from polyzymd.analyses.base import ComparisonResult

        if isinstance(result, ComparisonResult):
            return format_scalar_comparison(
                result,
                title="Radius of Gyration Comparison",
                metric_label="Mean Rg",
                metric_unit="Å",
                metric_key="mean_rg",
                output_format=output_format,
                higher_is_better=False,
            )
        return super().format(result, output_format)
```

That file is now discoverable, runnable through compare workflows, and
compatible with default framework formatting/statistics.

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
you write `@patch("polyzymd.analyses.rg.TrajectoryLoader")`. This replaces
the module-level name. If `TrajectoryLoader` were imported lazily inside the
method, the patch target would be different and easier to get wrong.

## Context Objects Reference

The framework provides context objects so plugins never need to load configs
or discover paths themselves:

| Context | Passed To | Key Attributes |
|---------|-----------|----------------|
| `ReplicateContext` | `compute_replicate()` | `.sim_config`, `.settings`, `.output_dir`, `.replicate`, `.recompute`, `.equilibration`, `.result_path` |
| `AggregateContext` | `aggregate()` | `.condition`, `.replicates`, `.output_dir`, `.settings`, `.equilibration`, `.result_path` |
| `ComparisonContext` | `compare()` | `.conditions`, `.analysis_dirs`, `.results_dir`, `.effective_control`, `.settings`, `.recompute` |
| `PlotContext` | `plot()` | `.conditions`, `.analysis_dirs`, `.output_dir`, `.settings`, `.plot_settings` |

(plotcontext-plot-settings)=
## `PlotContext.plot_settings`

`PlotContext.plot_settings` carries global plotting configuration from
`polyzymd.compare.config.PlotSettings`. Use this object instead of hard-coding
figure style, DPI, or output behavior.

### Type and global fields

`ctx.plot_settings` is either:

- an instance of `PlotSettings`, or
- `None` (for backward compatibility)

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

All plugins should follow this guard pattern:

```python
plot_settings = ctx.plot_settings
if plot_settings is None:
    from polyzymd.compare.config import PlotSettings

    plot_settings = PlotSettings()
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
    if plot_settings is None:
        from polyzymd.compare.config import PlotSettings

        plot_settings = PlotSettings()

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
aggregated results. Use `self._load_aggregated_result(agg_dir)` — inherited
from `Analysis` — to load each condition's data:

```python
from pathlib import Path

from polyzymd.analyses.base import PlotContext


def plot(self, ctx: PlotContext) -> list[Path]:
    import matplotlib.pyplot as plt

    for cond in ctx.conditions:
        agg_dir = ctx.analysis_dirs[cond.label] / "aggregated"
        summary = self._load_aggregated_result(agg_dir)
        if summary is not None:
            # summary is a dict or AggregatedResultClass instance
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
`compute_replicate()` and `aggregate()`. Both are supported, but they follow
different deserialization paths.

### Dict path

Dict plugins are straightforward:

- `compute_replicate()` returns `dict`
- `aggregate()` returns `dict`
- no `AggregatedResultClass` class variable is required
- framework saves/loads JSON using standard `json.dumps()` / `json.loads()`

This path is often best for scalar analyses and quick plugin development.

```python
class MyAnalysis(Analysis):
    name = "my_analysis"

    def compute_replicate(self, ctx, replicate):
        return {"value": 1.23, "replicate": replicate}

    def aggregate(self, ctx, results):
        vals = [r["value"] for r in results]
        return {"mean_value": sum(vals) / len(vals), "replicate_values": vals}
```

### Pydantic model path

Use models when you want validation, nested structure, strict typing, and
sidecar storage for large arrays.

Key pieces:

1. Result models inherit from `BaseAnalysisResult`
   (`polyzymd.analyses._results_base.BaseAnalysisResult`)
2. `BaseAnalysisResult` provides `save()` and `load()`
3. `BaseAnalysisResult` supports NPZ sidecar storage when models include large
   NumPy data
4. Plugin class sets `AggregatedResultClass` to the aggregated model class
5. Framework deserializes aggregated files through
   `AggregatedResultClass.load(path)` (or model-JSON validation fallback)

```python
from pathlib import Path
from typing import ClassVar

from pydantic import BaseModel

from polyzymd.analyses._results_base import BaseAnalysisResult
from polyzymd.analyses.base import Analysis


class ReplicateResult(BaseAnalysisResult):
    mean_rg: float
    replicate: int


class AggregatedResult(BaseAnalysisResult):
    mean_rg: float
    sem_rg: float
    replicate_values: list[float]


class RgAnalysis(Analysis):
    name: ClassVar[str] = "rg"
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

### Quick decision guide

- Choose **dicts** when your output is small, flat, and naturally JSON
- Choose **Pydantic models** when you need validation, richer structure,
  predictable typing, or NPZ sidecar support

## A Note on the `compare/` Package

You may notice a `compare/` package with `compare/results/` and
`compare/formatters.py`. Existing plugins reference these because they were
designed in an earlier architecture stage.

**You do NOT need to create files in `compare/` for a new plugin.** Keep your
plotting logic in your plugin's `plot()` method and your formatting in
`format()`. Everything can live in a single file.

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

Other submodules (`selections`, `aggregation`, `diagnostics`, `aa_classification`,
`metric_type`, `logging_utils`) are available via direct import but are not
re-exported from the package root:

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

- [ ] File: `src/polyzymd/analyses/<name>.py`
- [ ] `name: ClassVar[str]` — unique, lowercase, used in CLI
- [ ] `Settings` inner class — Pydantic v2 BaseModel with defaults
- [ ] `compute_replicate()` — uses `TrajectoryLoader`, returns dict or Pydantic model
- [ ] `aggregate()` — returns dict or Pydantic model; framework auto-saves to `ctx.result_path`
- [ ] Choose comparison path:
  - [ ] Default: implement `extract_metrics()`
  - [ ] Custom: override `compare()` returning a model with `.save()`
- [ ] (Optional) `plot()` — matplotlib figures, load data via `self._load_aggregated_result()`
- [ ] (Optional) `format()` — CLI display via `format_scalar_comparison()`
- [ ] If using `AggregatedResultClass`, return matching result model instances
- [ ] Tests: `tests/test_<name>_plugin.py`
- [ ] Verify: `pixi run -e build pytest tests/test_<name>_plugin.py -v`

## What a Good Community PR Looks Like

A strong contribution usually includes:

- one new plugin file in `src/polyzymd/analyses/`
- one focused test file
- one short docs update showing the `plugins:` config block
- figures only if they add scientific value

Aim for a plugin that another MDAnalysis user can understand in one read.

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

New plugins should keep all logic in a single file. The `compare/results/`
directory is primarily used by existing plugins designed under older
architecture constraints.

## Existing Plugins to Study

| Plugin | Comparison path | Good example of |
|--------|------------------|-----------------|
| `secondary_structure/` | Default (`extract_metrics`) | Wrapping existing calculators and using typed result classes |
| `rmsf/` | Default (`extract_metrics`) | Default scalar comparison plus plotting/formatting workflow |
| `catalytic_triad/` | Default (`extract_metrics`) | Custom geometric computation inside plugin lifecycle |
| `contacts/` | Custom (`compare` override) | Multi-metric comparison and condition filtering |
| `distances/` | Custom (`compare` override) | Entry-table style comparison output |

```{important}
Do not interpret plugin differences as "simple" versus "complex" quality.
Several existing plugins reflect earlier architectural phases. Use them as
pattern references for specific tasks (result modeling, comparison style,
plotting), not as maturity rankings.
```

## Testing Your Plugin

### Test file naming

Name your test file `tests/test_<name>_plugin.py` (for example,
`tests/test_rg_plugin.py`).

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


class TestRgDiscovery:
    """Test that the plugin is discoverable."""

    def test_discovered(self):
        assert "rg" in list_analyses()

    def test_class_variables(self):
        cls = get_analysis("rg")
        assert cls.name == "rg"
        assert hasattr(cls, "Settings")

    def test_settings_defaults(self):
        cls = get_analysis("rg")
        settings = cls.Settings()
        assert settings.selection == "protein and name CA"


class TestRgMetrics:
    """Test metric extraction for default comparison."""

    def test_extract_metrics(self):
        cls = get_analysis("rg")
        analysis = cls()
        summary = {
            "mean_rg": 15.2,
            "sem_rg": 0.3,
            "replicate_values": [14.9, 15.5],
        }
        metrics = analysis.extract_metrics(summary)
        assert "mean_rg" in metrics
        assert isinstance(metrics["mean_rg"], MetricValue)
        assert metrics["mean_rg"].higher_is_better is False
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


class TestRgCompute:
    """Test compute_replicate with mocked trajectories."""

    @patch("polyzymd.analyses.rg.TrajectoryLoader")
    def test_computes_metric(self, mock_loader_cls, tmp_path):
        cls = get_analysis("rg")
        analysis = cls()
        settings = cls.Settings()

        mock_loader = MagicMock()
        mock_loader_cls.return_value = mock_loader

        mock_universe = MagicMock()
        mock_atoms = MagicMock()
        mock_atoms.radius_of_gyration = MagicMock(side_effect=[15.0 + i * 0.01 for i in range(50)])
        mock_universe.select_atoms.return_value = mock_atoms

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
        assert "mean_rg" in result
        assert result["mean_rg"] > 0
        assert result["replicate"] == 1
```

Note the patch target:
`@patch("polyzymd.analyses.rg.TrajectoryLoader")`. This works because
`TrajectoryLoader` is imported at module level in `rg.py`.

### Testing `aggregate()`

Aggregation takes plain data and needs no trajectory mocks:

```python
from pathlib import Path
from unittest.mock import MagicMock

import numpy as np

from polyzymd.analyses import get_analysis
from polyzymd.analyses.base import AggregateContext, Condition


class TestRgAggregate:
    """Test aggregation across replicates."""

    def test_aggregate(self, tmp_path):
        cls = get_analysis("rg")
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
            {"mean_rg": 15.0, "replicate": 1},
            {"mean_rg": 15.5, "replicate": 2},
            {"mean_rg": 14.8, "replicate": 3},
        ]

        agg = analysis.aggregate(ctx, results)
        expected_mean = np.mean([15.0, 15.5, 14.8])
        assert abs(agg["mean_rg"] - expected_mean) < 1e-10
        assert agg["replicate_values"] == [15.0, 15.5, 14.8]
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
        cls = get_analysis("rg")
        analysis = cls()

        # Set up filesystem with aggregated results
        for label, value in [("WT", 15.0), ("Mutant", 16.5)]:
            agg_dir = tmp_path / label / "rg" / "aggregated"
            agg_dir.mkdir(parents=True)
            (agg_dir / "result.json").write_text(
                json.dumps(
                    {
                        "mean_rg": value,
                        "sem_rg": 0.3,
                        "replicate_values": [value - 0.2, value + 0.2],
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
                "WT": tmp_path / "WT" / "rg",
                "Mutant": tmp_path / "Mutant" / "rg",
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
pixi run -e build pytest tests/test_rg_plugin.py -v
```

## See Also

- {doc}`architecture` — Overall system architecture
- {doc}`analysis_statistics_best_practices` — Autocorrelation and uncertainty
- {doc}`analysis_compare_conditions` — User guide for running comparisons
- API reference: `polyzymd.analyses.base.Analysis`
