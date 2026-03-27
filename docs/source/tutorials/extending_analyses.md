# Extending the Analysis Framework

This guide shows you how to add a new analysis type to PolyzyMD in one file.
By the end, you will have a working analysis plugin that can compute metrics
per replicate, aggregate across replicates, compare conditions, generate
figures, and display results on the CLI.

```{note}
This guide covers the **plugin system** in ``analyses/``.
The underlying per-condition calculators (private ``_calculator_*.py`` modules
inside ``analyses/``) are implementation details — you do not need to touch them.
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
    result = compute_replicate(ctx, replicate)    # your code
  aggregated = aggregate(ctx, [replicate_results])  # your code

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
implement `extract_metrics()`, you do **not** need to touch `compare()`.

### `plot(ctx)` (optional)

Generate matplotlib figures. Use `self._load_aggregated_result(agg_dir)` to
load each condition's data. Return a list of figure paths.

### `format(result, output_format)` (optional)

Custom CLI display. Use `format_scalar_comparison()` from `analyses.stats`
for formatted tables with rankings, effect sizes, and significance stars.

## The Complete Example

Here is a single, complete plugin file that uses real framework utilities.
This is the pattern you should follow for new analyses.

Create `src/polyzymd/analyses/rg.py`:

```python
"""Radius of gyration analysis plugin."""
from __future__ import annotations

import json
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
from polyzymd.analyses.shared import TrajectoryLoader, parse_time_string, convert_time
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
            calculation. Defaults to C-alpha atoms.
        """

        selection: str = Field(
            default="protein and name CA",
            description="MDAnalysis atom selection for Rg calculation",
        )

    # --- Required methods ---

    def compute_replicate(
        self, ctx: ReplicateContext, replicate: int
    ) -> dict[str, Any]:
        """Compute mean Rg for one replicate.

        Uses TrajectoryLoader for topology discovery, trajectory segment
        daisy-chaining, and equilibration frame skipping.
        """
        # Lazy-import heavy third-party dependency
        import MDAnalysis as mda

        # Use TrajectoryLoader — handles topology, trajectory segments,
        # timestep detection, and equilibration offset
        loader = TrajectoryLoader(ctx.sim_config)
        u = loader.load_universe(replicate)

        # Convert equilibration time to frame offset
        eq_value, eq_unit = parse_time_string(ctx.equilibration)
        timestep_ps = loader.get_timestep(replicate, unit="ps")
        eq_time_ps = convert_time(eq_value, eq_unit, "ps")
        start_frame = int(eq_time_ps / timestep_ps)

        # Compute Rg over production frames
        atoms = u.select_atoms(ctx.settings.selection)
        rg_values = []
        for ts in u.trajectory[start_frame:]:
            rg_values.append(atoms.radius_of_gyration())

        return {
            "mean_rg": float(np.mean(rg_values)),
            "std_rg": float(np.std(rg_values)),
            "n_frames": len(rg_values),
            "replicate": replicate,
        }

    def aggregate(
        self, ctx: AggregateContext, results: Sequence[Any]
    ) -> dict[str, Any]:
        """Average Rg across replicates with SEM.

        The framework auto-saves the return value to ctx.result_path.
        """
        values = [r["mean_rg"] for r in results]
        return {
            "mean_rg": float(np.mean(values)),
            "sem_rg": float(np.std(values, ddof=1) / np.sqrt(len(values))),
            "replicate_values": values,
        }

    # --- Default comparison path ---

    def extract_metrics(self, summary: Any) -> dict[str, MetricValue]:
        """Expose Rg as a scalar metric for automatic t-tests and ANOVA.

        Returning MetricValue with higher_is_better=False tells the
        framework that lower Rg ranks better (more compact protein).
        """
        return {
            "mean_rg": MetricValue(
                name="mean_rg",
                mean=summary["mean_rg"],
                sem=summary["sem_rg"],
                replicate_values=summary["replicate_values"],
                higher_is_better=False,
                direction_labels=("compacting", "unchanged", "expanding"),
            ),
        }

    # --- Optional: plots ---

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

    # --- Optional: CLI formatting ---

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

That file is the complete plugin. It is now:

- Discovered by `polyzymd.analyses.list_analyses()`
- Runnable via `polyzymd compare run rg -f comparison.yaml`
- Automatically compared with t-tests, ANOVA, and ranking
- Plotted and formatted on the CLI

## Import Rules

The import structure in the example above is intentional. Getting imports
wrong causes subtle test failures, so follow these rules:

**Module-level imports** (top of file):

- Standard library (`json`, `logging`, `pathlib`)
- NumPy — imported at module level because it's available in all environments
- Framework utilities from `analyses/shared/` (`TrajectoryLoader`,
  `parse_time_string`, `convert_time`, `AlignmentConfig`, etc.)
- Framework base classes from `analyses/base` and `analyses/stats`

**Lazy imports inside methods** (inside `compute_replicate`, `plot`, etc.):

- `MDAnalysis` / `mdtraj` — heavy simulation libraries
- `matplotlib` / `seaborn` — plotting libraries
- Any package that may not be installed in all environments

```python
# At module level — always available, needed for @patch targets in tests
from polyzymd.analyses.shared import TrajectoryLoader, parse_time_string, convert_time

# Inside methods — heavy deps that may not be installed
def compute_replicate(self, ctx, replicate):
    import MDAnalysis as mda
    ...

def plot(self, ctx):
    import matplotlib.pyplot as plt
    ...
```

**Why this matters for testing:** When you mock `TrajectoryLoader` in tests,
you write `@patch("polyzymd.analyses.rg.TrajectoryLoader")`. This replaces
the module-level name. If `TrajectoryLoader` were imported lazily inside the
method, the patch target would be different and harder to mock correctly.

## Context Objects Reference

The framework provides context objects so plugins never need to load configs
or discover paths themselves:

| Context | Passed To | Key Attributes |
|---------|-----------|----------------|
| `ReplicateContext` | `compute_replicate()` | `.sim_config`, `.settings`, `.output_dir`, `.replicate`, `.recompute`, `.equilibration`, `.result_path` |
| `AggregateContext` | `aggregate()` | `.condition`, `.replicates`, `.output_dir`, `.settings`, `.equilibration`, `.result_path` |
| `ComparisonContext` | `compare()` | `.conditions`, `.analysis_dirs`, `.results_dir`, `.effective_control`, `.settings`, `.recompute` |
| `PlotContext` | `plot()` | `.conditions`, `.analysis_dirs`, `.output_dir`, `.settings`, `.plot_settings` |

## Choosing Your Comparison Path

### Path A: Default Scalar Comparison (recommended)

Implement `extract_metrics()` and you get t-tests, Cohen's d, ANOVA, and
ranking for free. This is what the example above does.

The framework loads aggregated results from disk automatically using
`json.loads()` (for dict results) or your `AggregatedResultClass` (for
Pydantic model results). You do **not** need to implement
`_deserialize_result()`.

### Path B: Custom Comparison

Override `compare()` entirely for multi-metric or entry-table analyses:

```python
def compare(self, ctx: ComparisonContext) -> MyComparisonResult:
    # Load results, compute custom statistics, return your Pydantic model
    # The returned object must have a .save(path) method
    ...
```

See `src/polyzymd/analyses/contacts.py` or `src/polyzymd/analyses/distances.py`
for full examples.

## Loading Results in `plot()`

`PlotContext` provides paths and conditions but does **not** carry pre-loaded
aggregated results. Use `self._load_aggregated_result(agg_dir)` — inherited
from `Analysis` — to load each condition's data:

```python
def plot(self, ctx: PlotContext) -> list[Path]:
    import matplotlib.pyplot as plt

    for cond in ctx.conditions:
        agg_dir = ctx.analysis_dirs[cond.label] / "aggregated"
        summary = self._load_aggregated_result(agg_dir)
        if summary is not None:
            # summary is a dict (or your AggregatedResultClass instance)
            # ... plot data from summary ...
```

`_load_aggregated_result()` looks for `result.json` in the aggregated
directory, falling back to the most recent `*.json` file. It deserializes
using `AggregatedResultClass` if set, otherwise `json.loads()`. Returns
`None` if no file is found.

## Formatting with `format_scalar_comparison()`

For analyses using the default comparison path, `format_scalar_comparison()`
from `analyses.stats` produces formatted tables with rankings, effect sizes,
and significance stars:

```python
from polyzymd.analyses.stats import format_scalar_comparison

def format(self, result: Any, output_format: str = "text") -> str:
    from polyzymd.analyses.base import ComparisonResult

    if isinstance(result, ComparisonResult):
        return format_scalar_comparison(
            result,
            title="My Analysis Comparison",
            metric_label="My Metric",      # column header in output tables
            metric_unit="Å",               # appended to values
            metric_key="my_metric",        # key prefix in ConditionSummary extra fields
            output_format=output_format,   # "text", "markdown", or "json"
            higher_is_better=True,         # affects interpretation wording
        )
    return super().format(result, output_format)
```

**Parameters:**

| Parameter | Type | Description |
|-----------|------|-------------|
| `result` | `ComparisonResult` | The result to format |
| `title` | `str` | Display title (e.g. `"RMSF Comparison"`) |
| `metric_label` | `str` | Human-readable metric name for table headers |
| `metric_unit` | `str` | Unit string appended to values (e.g. `"Å"`) |
| `metric_key` | `str \| None` | Key prefix in `ConditionSummary` extra fields; auto-detected if `None` |
| `output_format` | `str` | `"text"`, `"markdown"`, or `"json"` |
| `higher_is_better` | `bool` | Affects interpretation wording |

## Return Types: Dicts vs Pydantic Models

The example above returns plain Python dicts from `compute_replicate()`
and `aggregate()`. This is the **recommended approach for new plugins** — it
is simpler, requires less code, and the framework handles serialization
automatically via `json.dumps()`.

**When to use dicts (recommended for most plugins):**
- Single-metric scalar analyses
- Results that serialize naturally to JSON (numbers, strings, lists)
- When you want to minimize boilerplate

**When to graduate to Pydantic models:**
- You need validation on result fields
- You have complex nested results (per-residue arrays, distance matrices)
- You want IDE autocompletion on result fields
- You need NPZ sidecar storage for large numpy arrays

### Graduating to Pydantic Models

Existing plugins like `rmsf.py` and `secondary_structure.py` use typed Pydantic
result models defined in private `_results_*.py` files. These inherit from
`BaseAnalysisResult` in `analyses/_results_base.py`, which provides
`save()`/`load()` methods and optional NPZ sidecar support for large arrays.

If you follow this pattern:

1. Define your result models in `analyses/_results_<name>.py`
2. Inherit from `BaseAnalysisResult`
3. Return model instances from `compute_replicate()` and `aggregate()`
4. Set `AggregatedResultClass` on your plugin class:

```python
from polyzymd.analyses._results_rg import RgAggregatedResult

class RgAnalysis(Analysis):
    name: ClassVar[str] = "rg"
    AggregatedResultClass = RgAggregatedResult
    ...
```

The framework uses `AggregatedResultClass` to deserialize aggregated results
from disk — trying `.load(path)` first (if available), then
`.model_validate_json()`. For plain-dict plugins, no class variable is needed;
the framework falls back to `json.loads()`.

## A Note on the `compare/` Package

You may notice a `compare/` package with `compare/plotters/`,
`compare/results/`, and `compare/formatters.py`. Existing plugins reference
these because they were written during an earlier version of the architecture.

**You do NOT need to create files in `compare/` for a new plugin.** Keep your
plotting logic in your plugin's `plot()` method and your formatting in
`format()`. Everything can live in a single file.

(shared-utilities)=
## Shared Utilities (`analyses/shared/`)

The `analyses/shared/` package provides reusable infrastructure that all
plugins can use. The example above uses `TrajectoryLoader`, `parse_time_string`,
and `convert_time`. Here is the full set of available utilities:

### TrajectoryLoader

The most important shared utility. Handles topology discovery, daisy-chain
trajectory segments, timestep detection, and equilibration frame skipping:

```python
from polyzymd.analyses.shared import TrajectoryLoader, parse_time_string, convert_time

loader = TrajectoryLoader(ctx.sim_config)
u = loader.load_universe(replicate)

# Convert equilibration time to frame offset
eq_value, eq_unit = parse_time_string(ctx.equilibration)
timestep_ps = loader.get_timestep(replicate, unit="ps")
eq_time_ps = convert_time(eq_value, eq_unit, "ps")
start_frame = int(eq_time_ps / timestep_ps)
```

### Other Available Utilities

| Module | What It Provides |
|--------|-----------------|
| `shared.loader` | `TrajectoryLoader`, `parse_time_string`, `convert_time`, `time_to_frame` |
| `shared.alignment` | `AlignmentConfig`, `align_trajectory` — trajectory alignment with multiple reference modes |
| `shared.statistics` | `compute_sem`, `aggregate_per_residue_stats`, `weighted_mean_with_sem` |
| `shared.autocorrelation` | `compute_acf`, `estimate_correlation_time`, `statistical_inefficiency` |
| `shared.selections` | Extended selection syntax: `midpoint(...)`, `com(...)` |
| `shared.aggregation` | `collect_replicate_results`, `aggregate_distance_pair_stats` |
| `shared.pbc` | `minimum_image_distance`, `pairwise_distances_pbc` |
| `shared.constants` | `DEFAULT_CONTACT_CUTOFF`, `DEFAULT_DISTANCE_THRESHOLD` |

Import directly or via the convenience re-exports:

```python
from polyzymd.analyses.shared import TrajectoryLoader, compute_sem, AlignmentConfig
```

## Contributor Checklist

When creating a new analysis plugin:

- [ ] File: `src/polyzymd/analyses/<name>.py`
- [ ] `name: ClassVar[str]` — unique, lowercase, used in CLI
- [ ] `Settings` inner class — Pydantic v2 BaseModel with defaults
- [ ] `compute_replicate()` — uses `TrajectoryLoader`, returns dict or Pydantic model
- [ ] `aggregate()` — returns dict or Pydantic model; framework auto-saves to `ctx.result_path`
- [ ] Choose comparison path:
  - [ ] Default: implement `extract_metrics()` — framework handles deserialization
  - [ ] Custom: override `compare()` returning a model with `.save()`
- [ ] (Optional) `plot()` — matplotlib figures, load data via `self._load_aggregated_result()`
- [ ] (Optional) `format()` — CLI display via `format_scalar_comparison()`
- [ ] Tests: `tests/test_<name>_plugin.py`
- [ ] Verify: `pixi run -e build pytest tests/test_<name>_plugin.py -v`

## What a Good Community PR Looks Like

A strong contribution usually includes:

- one new plugin file in `src/polyzymd/analyses/`
- one focused test file
- one short docs update showing the `plugins:` config block
- figures only if they add real scientific value

Aim for a plugin that another MDAnalysis user can understand in one read.

## Anti-Patterns to Avoid

### Don't bypass the context

```python
# WRONG
def compute_replicate(self, ctx, replicate):
    config = SimulationConfig.from_yaml("/my/path/config.yaml")
```

```python
# RIGHT
def compute_replicate(self, ctx, replicate):
    sim_config = ctx.sim_config  # Already loaded by framework
```

### Understand the result-saving convention

The framework auto-saves return values from `compute_replicate()` and
`aggregate()` to `ctx.result_path`. For simple plugins, just return your dict
and the framework handles the rest.

For caching with equilibration-aware filenames (the established pattern in
existing plugins), save explicitly:

```python
def aggregate(self, ctx, results):
    agg = {"mean": 15.0}
    target_path = ctx.result_path
    if target_path is None:
        target_path = ctx.output_dir / "result.json"
    target_path.parent.mkdir(parents=True, exist_ok=True)
    target_path.write_text(json.dumps(agg, indent=2))
    return agg
```

### Don't import heavy deps at module level

```python
# WRONG — breaks environments without MDAnalysis
import MDAnalysis as mda

# RIGHT — lazy import inside method
def compute_replicate(self, ctx, replicate):
    import MDAnalysis as mda
```

Framework utilities from `analyses/shared/` are imported at **module level** —
they are always available and must be patchable in tests. Only heavy
*third-party* dependencies (MDAnalysis, mdtraj, matplotlib) go inside methods.

### Don't create files in `compare/`

New plugins should keep all logic in a single file. The `compare/plotters/`
and `compare/results/` directories are used by existing plugins for historical
reasons. New plugins should implement `plot()` and `format()` inline.

## Existing Plugins to Study

| Plugin | Complexity | Comparison | Good Example Of |
|--------|-----------|------------|-----------------|
| `secondary_structure.py` | Simple | Default (extract_metrics) | Wraps existing calculator, uses AggregatedResultClass |
| `rmsf.py` | Simple | Default (extract_metrics) | Default compare with formatting |
| `catalytic_triad.py` | Medium | Default (extract_metrics) | Custom compute with triad geometry |
| `contacts.py` | Complex | Custom (override compare) | Multi-metric comparison, condition filtering |
| `distances.py` | Complex | Custom (override compare) | Entry-table comparison results |

## Testing Your Plugin

### Test file naming

Name your test file `tests/test_<name>_plugin.py` (e.g.,
`tests/test_rg_plugin.py`).

### Basic tests (discovery, settings, metrics)

These tests use no mocks and verify the plugin integrates with the framework:

```python
import json
from pathlib import Path
from unittest.mock import MagicMock

import numpy as np
import pytest

from polyzymd.analyses import get_analysis, list_analyses
from polyzymd.analyses.base import (
    AggregateContext,
    Condition,
    MetricValue,
)


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
`TrajectoryLoader` at the **module level** to test compute logic without
real data:

```python
from unittest.mock import MagicMock, patch

from polyzymd.analyses.base import ReplicateContext


class TestRgCompute:
    """Test compute_replicate with mocked trajectories."""

    @patch("polyzymd.analyses.rg.TrajectoryLoader")
    def test_computes_metric(self, MockLoader, tmp_path):
        cls = get_analysis("rg")
        analysis = cls()
        settings = cls.Settings()

        # Mock the TrajectoryLoader
        mock_loader = MagicMock()
        MockLoader.return_value = mock_loader

        # Mock a universe with 50 frames
        mock_universe = MagicMock()
        mock_atoms = MagicMock()
        mock_atoms.__len__ = MagicMock(return_value=100)
        # Return deterministic Rg values
        mock_atoms.radius_of_gyration = MagicMock(
            side_effect=[15.0 + i * 0.01 for i in range(50)]
        )
        mock_universe.select_atoms.return_value = mock_atoms

        # Mock trajectory slicing
        mock_trajectory = MagicMock()
        mock_trajectory.__len__ = MagicMock(return_value=50)
        mock_trajectory.__getitem__ = MagicMock(return_value=range(50))
        mock_universe.trajectory = mock_trajectory

        mock_loader.load_universe.return_value = mock_universe
        mock_loader.get_timestep.return_value = 10.0  # 10 ps

        # Build context
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

Note the mock target: `@patch("polyzymd.analyses.rg.TrajectoryLoader")`. This
works because `TrajectoryLoader` is imported at module level in `rg.py`. If it
were imported lazily inside the method, you would need to patch
`polyzymd.analyses.shared.loader.TrajectoryLoader` instead — which is
fragile and non-obvious.

### Testing `aggregate()`

Aggregation takes plain data and needs no mocks:

```python
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

Run with:

```bash
pixi run -e build pytest tests/test_rg_plugin.py -v
```

## See Also

- {doc}`architecture` — Overall system architecture
- {doc}`analysis_statistics_best_practices` — Autocorrelation and uncertainty
- {doc}`analysis_compare_conditions` — User guide for running comparisons
- API reference: `polyzymd.analyses.base.Analysis`
