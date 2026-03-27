# Extending the Analysis Framework

This guide shows you how to add a new analysis type to PolyzyMD in one file.
By the end, you will have a working analysis plugin that can compute metrics
per replicate, aggregate across replicates, compare conditions, generate
figures, and display results on the CLI.

```{note}
This guide covers the **plugin system** in ``analyses/``.
For the underlying per-condition calculators, see the ``analysis/`` package
API docs.
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

    # --- Optional: enable default comparison ---

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
```

That's it. The plugin is now:

- Discovered by `polyzymd.analyses.list_analyses()`
- Runnable via `polyzymd compare run rg -f comparison.yaml`
- Automatically compared with t-tests, ANOVA, and ranking

## Step-by-Step Guide

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
- `extract_metrics()`

The framework can do the rest.

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
) -> Any:
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

### Step 4: Implement `aggregate()`

This runs once per condition after all replicates complete:

```python
def aggregate(
    self, ctx: AggregateContext, results: Sequence[Any]
) -> Any:
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

If using this path, you must also implement `_deserialize_result()` so the
framework can load your aggregated results from disk:

```python
def _deserialize_result(self, path: Path) -> Any:
    import json
    return json.loads(path.read_text())
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
    # ... plotting logic using ctx.analysis_dirs, ctx.conditions ...

    output_path = ctx.output_dir / f"{self.name}_comparison.png"
    fig.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return [output_path]
```

### Step 7: Add CLI Formatting (Optional)

Override `format()` for custom CLI output:

```python
def format(self, result: Any, output_format: str = "text") -> str:
    from polyzymd.analyses.stats import format_scalar_comparison
    return format_scalar_comparison(result, output_format)
```

For the default comparison path, `format_scalar_comparison()` already produces
well-formatted tables with rankings, effect sizes, and significance stars.

## Context Objects Reference

The framework provides context objects so plugins never need to load configs
or discover paths themselves:

| Context | Passed To | Key Attributes |
|---------|-----------|----------------|
| `ReplicateContext` | `compute_replicate()` | `.sim_config`, `.settings`, `.output_dir`, `.replicate`, `.recompute`, `.equilibration` |
| `AggregateContext` | `aggregate()` | `.condition`, `.replicates`, `.output_dir`, `.settings`, `.equilibration` |
| `ComparisonContext` | `compare()` | `.conditions`, `.analysis_dirs`, `.results_dir`, `.effective_control`, `.settings`, `.recompute` |
| `PlotContext` | `plot()` | `.conditions`, `.analysis_dirs`, `.output_dir`, `.settings`, `.plot_settings` |

## Method Lifecycle

```text
For each condition:
  For each replicate:
    compute_replicate(ctx, replicate)  →  per-replicate result
  aggregate(ctx, [replicate_results]) →  condition summary

filter_conditions(conditions)         →  filtered conditions list
compare(ctx)                          →  ComparisonResult
plot(ctx)                             →  [figure paths]
format(result, "text")                →  CLI output string
```

## Contributor Checklist

When creating a new analysis plugin:

- [ ] File: `src/polyzymd/analyses/<name>.py`
- [ ] `name: ClassVar[str]` — unique, lowercase, used in CLI
- [ ] `Settings` inner class — Pydantic v2 BaseModel with defaults
- [ ] `compute_replicate()` — returns dict/BaseModel per replicate
- [ ] `aggregate()` — returns dict/BaseModel per condition
- [ ] Choose comparison path:
  - [ ] Default: implement `extract_metrics()` + `_deserialize_result()`
  - [ ] Custom: override `compare()` returning a model with `.save()`
- [ ] (Optional) `plot()` — matplotlib figures
- [ ] (Optional) `format()` — CLI display
- [ ] Tests: `tests/test_analyses_<name>.py`
- [ ] Verify: `pixi run -e build pytest tests/ -v -k <name>`

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
your aggregated results.

### Don't import heavy deps at module level

```python
# WRONG — breaks environments without MDAnalysis
import MDAnalysis as mda

# RIGHT — lazy import inside method
def compute_replicate(self, ctx, replicate):
    import MDAnalysis as mda
```

## Existing Plugins to Study

| Plugin | Complexity | Comparison | Good Example Of |
|--------|-----------|------------|-----------------|
| `secondary_structure.py` | Simple | Default (extract_metrics) | Minimal plugin, wraps existing calculator |
| `rmsf.py` | Simple | Default (extract_metrics) | Default compare with formatting |
| `catalytic_triad.py` | Simple | Default (extract_metrics) | Custom compute with triad geometry |
| `contacts.py` | Complex | Custom (override compare) | Multi-metric comparison, condition filtering |
| `distances.py` | Complex | Custom (override compare) | Entry-table comparison results |

## Testing Your Plugin

Create `tests/test_analyses_<name>.py`:

```python
import pytest
from polyzymd.analyses import get_analysis, list_analyses
from polyzymd.analyses.base import MetricValue


class TestRgAnalysis:
    """Tests for rg analysis plugin."""

    def test_discovered(self):
        """Plugin is found by discovery."""
        assert "rg" in list_analyses()

    def test_class_variables(self):
        cls = get_analysis("rg")
        assert cls.name == "rg"
        assert hasattr(cls, "Settings")

    def test_settings_defaults(self):
        cls = get_analysis("rg")
        settings = cls.Settings()
        assert settings.selection == "protein and name CA"

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

Run with:

```bash
pixi run -e build pytest tests/test_analyses_rg.py -v
```

## See Also

- {doc}`architecture` — Overall system architecture
- {doc}`analysis_statistics_best_practices` — Autocorrelation and uncertainty
- {doc}`analysis_compare_conditions` — User guide for running comparisons
- API reference: `polyzymd.analyses.base.Analysis`
