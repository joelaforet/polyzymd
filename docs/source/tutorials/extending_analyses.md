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

## Overview

Every community extension in PolyzyMD should start as a single class in
`src/polyzymd/analyses/` that subclasses `Analysis`.

The goal is simple contributor ergonomics:

- one file to review
- one settings model to configure
- one class to discover from the CLI
- one PR that adds compute, comparison, and plots together

The framework handles:

- **Discovery**: `pkgutil`-based auto-discovery — no registries or imports
- **Replicate iteration**: Calls your `compute_replicate()` for each run
- **Result saving**: Fallback auto-save if your plugin doesn't save itself
- **Aggregation**: Calls your `aggregate()` for each condition
- **Comparison**: Default t-tests/ANOVA or your custom `compare()` override
- **Plotting**: Your optional `plot()` override
- **CLI wiring**: `polyzymd compare run <name>` works automatically

You implement the science. The framework does the plumbing.

If you already have working MDAnalysis code for a metric, the usual path is:

1. move that logic into `compute_replicate()`
2. summarize replicates in `aggregate()`
3. either expose scalar metrics with `extract_metrics()` or write a custom `compare()`
4. add `plot()` if the analysis needs figures

That is the whole extension model.

## Quick Start: One-File Plugin

Create `src/polyzymd/analyses/rg.py`:

```python
"""Radius of gyration analysis plugin."""
from __future__ import annotations

import json
from pathlib import Path
from typing import Any, ClassVar, Sequence

from pydantic import BaseModel, Field

from polyzymd.analyses.base import (
    AggregateContext,
    Analysis,
    MetricValue,
    ReplicateContext,
)


class RgAnalysis(Analysis):
    """Radius of gyration: compactness of the protein structure."""

    name: ClassVar[str] = "rg"

    class Settings(BaseModel):
        """Settings for radius of gyration analysis."""
        selection: str = Field(
            default="protein and name CA",
            description="MDAnalysis selection string for atoms to include",
        )

    # --- Required methods ---

    def compute_replicate(
        self, ctx: ReplicateContext, replicate: int
    ) -> dict[str, Any]:
        """Compute mean Rg for one replicate."""
        import MDAnalysis as mda
        import numpy as np

        sim_config = ctx.sim_config
        run_dir = sim_config.get_working_directory(replicate)
        topology = run_dir / "solvated_system.pdb"
        trajs = sorted(run_dir.glob("production_*/*_trajectory.dcd"))

        u = mda.Universe(str(topology), [str(t) for t in trajs])
        atoms = u.select_atoms(ctx.settings.selection)

        rg_values = [atoms.radius_of_gyration() for _ in u.trajectory]
        return {
            "mean_rg": float(np.mean(rg_values)),
            "replicate": replicate,
        }

    def aggregate(
        self, ctx: AggregateContext, results: Sequence[Any]
    ) -> dict[str, Any]:
        """Average Rg across replicates with SEM."""
        import numpy as np

        values = [r["mean_rg"] for r in results]
        return {
            "mean_rg": float(np.mean(values)),
            "sem_rg": float(np.std(values, ddof=1) / np.sqrt(len(values))),
            "replicate_values": values,
        }

    # --- Required for default comparison path ---

    def extract_metrics(self, summary: Any) -> dict[str, MetricValue]:
        """Expose Rg as a scalar metric for automatic t-tests and ANOVA."""
        return {
            "mean_rg": MetricValue(
                name="mean_rg",
                mean=summary["mean_rg"],
                sem=summary["sem_rg"],
                replicate_values=summary["replicate_values"],
                higher_is_better=False,  # lower Rg = more compact
                direction_labels=("compacting", "unchanged", "expanding"),
            ),
        }

    def _deserialize_result(self, path: Path) -> dict[str, Any]:
        """Load an aggregated result from JSON."""
        return json.loads(path.read_text())
```

That's it. The plugin is now:

- Discovered by `polyzymd.analyses.list_analyses()`
- Runnable via `polyzymd compare run rg -f comparison.yaml`
- Automatically compared with t-tests, ANOVA, and ranking

```{important}
The quick-start above is **complete and working**. Every method shown is
required for the default comparison path. If you remove `_deserialize_result()`,
the default `compare()` will raise `NotImplementedError` when it tries to load
your aggregated results from disk. If you remove `extract_metrics()`, the
default comparison will produce no results.
```

## Mental Model

Think of one plugin file as a complete analysis package:

- `Settings` answers "what can users configure?"
- `compute_replicate()` answers "how do I analyze one trajectory?"
- `aggregate()` answers "how do I summarize replicates for one condition?"
- `compare()` answers "how do I compare conditions?"
- `plot()` answers "what figures should this analysis produce?"
- `format()` answers "what should the CLI print?"

For many analyses, you only need:

- `Settings`
- `compute_replicate()`
- `aggregate()`
- `extract_metrics()` + `_deserialize_result()`

The framework handles comparison, result saving, and CLI wiring automatically.

## Step-by-Step Guide

### Step 1: Choose Your Analysis Name

The `name` class variable is used in CLI commands, config files, and output
directories. Pick a short, descriptive, lowercase string:

```python
name: ClassVar[str] = "rg"  # Used as: polyzymd compare run rg
```

You can also set aliases:

```python
aliases: ClassVar[tuple[str, ...]] = ("radius_of_gyration",)
```

### Step 2: Define Your Settings

Settings are a Pydantic v2 `BaseModel` defined as an inner class. They are
parsed from the `plugins:` section of `comparison.yaml`:

```python
class Settings(BaseModel):
    """Settings for Rg analysis."""
    selection: str = Field(
        default="protein and name CA",
        description="MDAnalysis atom selection",
    )
    skip_first_ns: float = Field(
        default=0.0,
        description="Skip this many ns from start of trajectory",
    )
```

In `comparison.yaml`:

```yaml
plugins:
  rg:
    selection: "protein and name CA"
    skip_first_ns: 10.0
```

### Step 3: Implement `compute_replicate()`

This runs once per replicate per condition. Use the context object — never
load configs yourself:

```python
def compute_replicate(
    self, ctx: ReplicateContext, replicate: int
) -> dict[str, Any]:
    # Use framework-provided context
    sim_config = ctx.sim_config        # SimulationConfig (already loaded)
    settings = ctx.settings            # Your Settings model
    output_dir = ctx.output_dir        # Where to write per-replicate results
    equilibration = ctx.equilibration  # e.g. "10ns"

    # Lazy-import heavy dependencies
    import MDAnalysis as mda

    # Load trajectory
    run_dir = sim_config.get_working_directory(replicate)
    topology = run_dir / "solvated_system.pdb"
    trajs = sorted(run_dir.glob("production_*/*_trajectory.dcd"))

    u = mda.Universe(str(topology), [str(t) for t in trajs])
    # ... compute your metric ...

    return {"my_metric": value, "replicate": replicate}
```

**Key rules:**
- Always use `ctx.sim_config`, never `SimulationConfig.from_yaml(...)` yourself
- Always use `ctx.settings` for your analysis parameters
- Lazy-import MDAnalysis, OpenMM, numpy inside the method
- Return a dict or Pydantic model — the framework collects these for `aggregate()`

```{tip}
For production plugins, use `TrajectoryLoader` from `analyses/shared/` instead
of loading trajectories manually. It handles topology discovery, daisy-chain
trajectory segments, timestep parsing, and equilibration skipping. See
{ref}`shared-utilities` below.
```

### Step 4: Implement `aggregate()`

This runs once per condition after all replicates complete:

```python
def aggregate(
    self, ctx: AggregateContext, results: Sequence[Any]
) -> dict[str, Any]:
    import numpy as np

    values = [r["my_metric"] for r in results]
    return {
        "mean_value": float(np.mean(values)),
        "sem_value": float(np.std(values, ddof=1) / np.sqrt(len(values))),
        "replicate_values": values,
    }
```

The `results` list contains exactly the objects your `compute_replicate()`
returned for each successful replicate.

```{note}
**Result saving convention.** Existing plugins save results to disk
themselves — typically to a custom-named file for per-replicate caching
(e.g. `rmsf_eq10ns.json`) and to `ctx.result_path` for aggregated results.
The orchestrator has a **fallback** auto-save that writes to
`ctx.result_path` only if the file doesn't already exist, so if your
plugin saves manually the fallback is a no-op.

For simple plugins that return plain dicts, you *can* skip manual saves
and let the fallback handle it.  But if you want per-replicate caching
with custom filenames (the established pattern), save explicitly.
```

### Step 5: Choose Your Comparison Path

You have two options:

#### Option A: Default Comparison (Recommended for Single-Metric Analyses)

Implement `extract_metrics()` to expose scalar metrics. The framework
automatically runs t-tests, ANOVA, Cohen's d, and ranking:

```python
def extract_metrics(self, summary: Any) -> dict[str, MetricValue]:
    return {
        "my_metric": MetricValue(
            name="my_metric",
            mean=summary["mean_value"],
            sem=summary["sem_value"],
            replicate_values=summary["replicate_values"],
            higher_is_better=True,  # or False if lower is better
            direction_labels=("decreased", "unchanged", "increased"),
        ),
    }
```

If using this path, you **must** also implement `_deserialize_result()` so the
framework can load your aggregated results from disk during comparison:

```python
def _deserialize_result(self, path: Path) -> dict[str, Any]:
    import json
    return json.loads(path.read_text())
```

```{warning}
This is the most common contributor mistake. If you implement
`extract_metrics()` but forget `_deserialize_result()`, the default
`compare()` will raise `NotImplementedError` at runtime with a helpful error
message telling you exactly what to add.
```

#### Option B: Custom Comparison

Override `compare()` entirely for multi-metric or complex analyses:

```python
def compare(self, ctx: ComparisonContext) -> MyComparisonResult:
    # Load results, compute custom statistics, return your Pydantic model
    # The returned object must have a .save(path) method
    ...
```

See `src/polyzymd/analyses/contacts.py` or `src/polyzymd/analyses/distances.py`
for full examples.

### Step 6: Add Plots (Optional)

Override `plot()` to generate matplotlib figures from your plugin:

```python
def plot(self, ctx: PlotContext) -> list[Path]:
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(10, 6))

    # Load aggregated results for each condition
    for cond in ctx.conditions:
        agg_dir = ctx.analysis_dirs[cond.label] / "aggregated"
        summary = self._load_aggregated_result(agg_dir)
        if summary is not None:
            # ... plot data from summary ...
            pass

    output_path = ctx.output_dir / f"{self.name}_comparison.png"
    fig.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return [output_path]
```

```{note}
`PlotContext` provides paths and conditions but does not carry pre-loaded
aggregated results. Use `self._load_aggregated_result(agg_dir)` (inherited
from `Analysis`) to load each condition's data from the `aggregated/`
directory. This keeps memory usage low for analyses with large results.
```

### Step 7: Add CLI Formatting (Optional)

Override `format()` for custom CLI output:

```python
def format(self, result: Any, output_format: str = "text") -> str:
    from polyzymd.analyses.base import ComparisonResult
    from polyzymd.analyses.stats import format_scalar_comparison

    if isinstance(result, ComparisonResult):
        return format_scalar_comparison(
            result,
            title="My Analysis Comparison",
            metric_label="My Metric",
            metric_unit="Å",
            metric_key="my_metric",
            output_format=output_format,
        )
    return super().format(result, output_format)
```

For the default comparison path, `format_scalar_comparison()` already produces
well-formatted tables with rankings, effect sizes, and significance stars.

## Context Objects Reference

The framework provides context objects so plugins never need to load configs
or discover paths themselves:

| Context | Passed To | Key Attributes |
|---------|-----------|----------------|
| `ReplicateContext` | `compute_replicate()` | `.sim_config`, `.settings`, `.output_dir`, `.replicate`, `.recompute`, `.equilibration`, `.result_path` |
| `AggregateContext` | `aggregate()` | `.condition`, `.replicates`, `.output_dir`, `.settings`, `.equilibration`, `.result_path` |
| `ComparisonContext` | `compare()` | `.conditions`, `.analysis_dirs`, `.results_dir`, `.effective_control`, `.settings`, `.recompute` |
| `PlotContext` | `plot()` | `.conditions`, `.analysis_dirs`, `.output_dir`, `.settings`, `.plot_settings` |

## Method Lifecycle

```text
For each condition:
  For each replicate:
    result = compute_replicate(ctx, replicate)
    if plugin didn't save: framework saves to ctx.result_path  ← fallback
  aggregated = aggregate(ctx, [replicate_results])
  if plugin didn't save: framework saves to ctx.result_path    ← fallback

filter_conditions(conditions)         →  filtered conditions list
compare(ctx)                          →  ComparisonResult
plot(ctx)                             →  [figure paths]
format(result, "text")                →  CLI output string
```

(shared-utilities)=
## Shared Utilities (`analyses/shared/`)

The `analyses/shared/` package provides reusable infrastructure that all
plugins can use. You are not required to use these — the quick-start example
works without them — but they handle common tasks correctly and save you
from re-implementing boilerplate.

### TrajectoryLoader

The most important shared utility. Handles topology discovery, daisy-chain
trajectory segments, timestep detection, and equilibration frame skipping:

```python
from polyzymd.analyses.shared import TrajectoryLoader, parse_time_string, convert_time

def compute_replicate(self, ctx, replicate):
    loader = TrajectoryLoader(ctx.sim_config)
    u = loader.load_universe(replicate)

    # Parse equilibration time and convert to frames
    eq_value, eq_unit = parse_time_string(ctx.equilibration)
    timestep_ps = loader.get_timestep(replicate, unit="ps")
    eq_time_ps = convert_time(eq_value, eq_unit, "ps")
    start_frame = int(eq_time_ps / timestep_ps)

    atoms = u.select_atoms(ctx.settings.selection)
    for ts in u.trajectory[start_frame:]:
        # ... per-frame computation ...
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

## Return Types: Dicts vs Pydantic Models

The quick-start example returns plain Python dicts from `compute_replicate()`
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

If you want to follow this pattern:

1. Define your result models in `analyses/_results_<name>.py`
2. Inherit from `BaseAnalysisResult` (for per-replicate) and
   `AggregatedResultMixin` (for aggregated results)
3. Return model instances from `compute_replicate()` and `aggregate()`
4. Your `_deserialize_result()` becomes:
   ```python
   def _deserialize_result(self, path: Path) -> MyAggregatedResult:
       return MyAggregatedResult.model_validate_json(path.read_text())
   ```

This is an **optional upgrade** — dicts work perfectly for the default
comparison pipeline.

## A Note on the `compare/` Package

You may notice a `compare/` package with `compare/plotters/`,
`compare/results/`, and `compare/formatters.py`. Existing plugins reference
these because they were written during an earlier version of the architecture.

**You do NOT need to create files in `compare/` for a new plugin.** Keep your
plotting logic in your plugin's `plot()` method and your formatting in
`format()`. Everything can live in a single file.

## Contributor Checklist

When creating a new analysis plugin:

- [ ] File: `src/polyzymd/analyses/<name>.py`
- [ ] `name: ClassVar[str]` — unique, lowercase, used in CLI
- [ ] `Settings` inner class — Pydantic v2 BaseModel with defaults
- [ ] `compute_replicate()` — returns dict or Pydantic model per replicate
- [ ] `aggregate()` — returns dict or Pydantic model per condition; save to `ctx.result_path` or let fallback handle it
- [ ] Choose comparison path:
  - [ ] Default: implement `extract_metrics()` **and** `_deserialize_result()`
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

### Don't forget `_deserialize_result()` when using default compare

If you implement `extract_metrics()` but not `_deserialize_result()`, the
default `compare()` will raise `NotImplementedError` when it tries to load
your aggregated results. The error message tells you exactly what to add.

### Understand the result-saving convention

Existing plugins save results to disk explicitly — using custom filenames for
per-replicate caching and `ctx.result_path` for aggregated results. The
orchestrator has a **fallback** that writes to `ctx.result_path` only if
the file doesn't already exist.

For a minimal plugin, you can skip manual saves and rely on the fallback:

```python
# Simple — let the fallback auto-save handle it
def aggregate(self, ctx, results):
    return {"mean": 15.0}
```

For caching with equilibration-aware filenames (the established pattern),
save explicitly like existing plugins do:

```python
# Explicit save — matches rmsf.py, contacts.py, etc.
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

### Don't create files in `compare/`

New plugins should keep all logic in a single file. The `compare/plotters/`
and `compare/results/` directories are used by existing plugins for historical
reasons. New plugins should implement `plot()` and `format()` inline.

## Existing Plugins to Study

| Plugin | Complexity | Comparison | Good Example Of |
|--------|-----------|------------|-----------------|
| `rg.py` | Simple | Default (extract_metrics) | Minimal complete plugin with plots |
| `secondary_structure.py` | Simple | Default (extract_metrics) | Wraps existing calculator |
| `rmsf.py` | Simple | Default (extract_metrics) | Default compare with formatting |
| `catalytic_triad.py` | Simple | Default (extract_metrics) | Custom compute with triad geometry |
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


class TestMyPluginDiscovery:
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


class TestMyPluginMetrics:
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

`compute_replicate()` requires MDAnalysis and trajectory files. Mock the
`TrajectoryLoader` (or MDAnalysis Universe) to test compute logic without
real data:

```python
from unittest.mock import MagicMock, patch

from polyzymd.analyses.base import ReplicateContext


class TestMyPluginCompute:
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

### Testing `aggregate()`

Aggregation takes plain data and needs no mocks:

```python
class TestMyPluginAggregate:
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
