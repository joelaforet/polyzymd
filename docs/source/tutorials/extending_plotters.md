# Extending the Plotter Framework

This guide shows developers how to create custom plotters using PolyzyMD's registry-based framework. The framework follows the **Open-Closed Principle** (open for extension, closed for modification) and provides automatic discovery of plotter implementations.

```{note}
This guide covers the **plotter** subsystem for generating comparison visualizations.
For statistical comparisons, see {doc}`extending_comparators`.
```

## Overview

The plotter framework provides:

- **BasePlotter**: Abstract base class with shared utilities (`_save_figure`, `_get_output_path`)
- **PlotterRegistry**: Registry for auto-discovery of plotter types via `@register` decorator
- **PlotSettings**: Configuration model for plot appearance (DPI, colors, formats)
- **Automatic discovery**: `plot_all()` finds all registered plotters for an analysis type

### Architecture

```
compare/
├── plotter.py              # BasePlotter, PlotterRegistry, ComparisonPlotter
├── plotters/
│   ├── triad.py            # TriadKDEPanelPlotter, TriadThresholdBarsPlotter
│   ├── contacts.py         # BindingPreferenceHeatmapPlotter, BindingPreferenceBarPlotter
│   ├── rmsf.py             # RMSFBarPlotter, RMSFLinePlotter
│   └── __init__.py
└── config.py               # PlotSettings, PlotSettingsTriad, PlotSettingsContacts
```

## Quick Start: Minimal Plotter

Here's the minimal code to create a new plotter:

```python
from pathlib import Path
from typing import Any, Sequence

from polyzymd.compare.plotter import BasePlotter, PlotterRegistry


@PlotterRegistry.register("my_custom_plot")
class MyCustomPlotter(BasePlotter):
    """Generate custom visualization for my analysis type."""

    @classmethod
    def plot_type(cls) -> str:
        return "my_custom_plot"

    def can_plot(self, comparison_config, analysis_type: str) -> bool:
        """Return True if this plotter handles the analysis type."""
        return analysis_type == "my_analysis"

    def plot(
        self,
        data: dict[str, Any],
        labels: Sequence[str],
        output_dir: Path,
        **kwargs,
    ) -> list[Path]:
        """Generate and save plot(s)."""
        import matplotlib.pyplot as plt

        # Load your data from each condition's analysis_dir
        results = self._load_results(data, labels)

        if not results:
            return []

        # Create your visualization
        fig, ax = plt.subplots(figsize=(10, 6), dpi=self.settings.dpi)
        # ... plotting logic ...

        # Save using helper method
        output_path = self._get_output_path(output_dir, "my_custom_plot")
        return [self._save_figure(fig, output_path)]

    def _load_results(self, data, labels):
        """Load analysis results from filesystem."""
        from my_module.results import MyResult

        results = {}
        for label in labels:
            cond_data = data.get(label)
            if cond_data is None:
                continue

            analysis_dir = Path(cond_data["analysis_dir"])
            result_file = analysis_dir / "my_result.json"

            if result_file.exists():
                results[label] = MyResult.load(result_file)

        return results
```

## Critical: The Data Contract

```{warning}
The most common mistake when implementing plotters is expecting data to be
passed via `kwargs`. **This does not work!** The orchestrator only provides
filesystem paths—plotters must load their own data.
```

### What `plot()` Receives

The `ComparisonPlotter.plot_analysis()` method calls your plotter with:

```python
plotter.plot(
    data=data,           # Dict of condition metadata (see structure below)
    labels=labels,       # List of condition labels in display order
    output_dir=output_dir,  # Where to save plots
    # **kwargs is reserved for future use—DO NOT rely on it
)
```

### The `data` Dictionary Structure

The `data` dict contains metadata for each condition, **not analysis results**:

```python
data = {
    "SBMA_75_25": {
        "condition": ConditionConfig(...),      # Condition metadata
        "sim_config": SimulationConfig(...),    # Full simulation config
        "analysis_dir": Path("path/to/analysis/contacts/"),   # CRITICAL!
        "aggregated_dir": Path("path/to/analysis/contacts/aggregated/"),
        "replicates": [1, 2, 3],                # Replicate numbers
    },
    "EGMA_75_25": {
        # Same structure...
    },
}
```

### Correct Pattern: Load From Filesystem

`````{tab-set}
````{tab-item} Correct ✓
:sync: correct

```python
def plot(self, data, labels, output_dir, **kwargs):
    # Load data from filesystem paths
    for label in labels:
        analysis_dir = Path(data[label]["analysis_dir"])
        result_file = analysis_dir / "my_result_aggregated.json"
        result = MyAggregatedResult.load(result_file)
        # ... use result for plotting ...
```
````

````{tab-item} Incorrect ✗
:sync: incorrect

```python
def plot(self, data, labels, output_dir, **kwargs):
    # WRONG: Expecting pre-loaded result in kwargs
    comparison_result = kwargs.get("comparison_result")
    if comparison_result is None:
        return []  # Always returns empty!
```
````
`````

## Step-by-Step Guide

### Step 1: Choose Your Analysis Type

Each plotter handles one analysis type (e.g., `"contacts"`, `"catalytic_triad"`, `"rmsf"`).
The `can_plot()` method determines which analysis types your plotter supports:

```python
def can_plot(self, comparison_config: "ComparisonConfig", analysis_type: str) -> bool:
    """Check if this plotter can handle the analysis type."""
    if analysis_type != "contacts":
        return False

    # Optionally check settings to enable/disable
    return self.settings.contacts.generate_my_plot
```

### Step 2: Register Your Plotter

Use the `@PlotterRegistry.register()` decorator with a unique key:

```python
@PlotterRegistry.register("binding_preference_heatmap")
class BindingPreferenceHeatmapPlotter(BasePlotter):
    ...
```

The registry key should be:
- Lowercase with underscores
- Descriptive of the plot type
- Unique across all plotters

### Step 3: Implement Required Methods

#### `plot_type()` (classmethod)

Return the registry key:

```python
@classmethod
def plot_type(cls) -> str:
    return "binding_preference_heatmap"
```

#### `can_plot()`

Check if this plotter should be called:

```python
def can_plot(self, comparison_config: "ComparisonConfig", analysis_type: str) -> bool:
    # Only handle contacts analysis
    if analysis_type != "contacts":
        return False

    # Check if this specific plot type is enabled in settings
    return self.settings.contacts.generate_enrichment_heatmap
```

#### `plot()`

Generate and save the visualization. Key responsibilities:

1. **Load data from filesystem** (not kwargs!)
2. **Create matplotlib figure(s)**
3. **Save using `_save_figure()`**
4. **Return list of output paths**

```python
def plot(
    self,
    data: dict[str, Any],
    labels: Sequence[str],
    output_dir: Path,
    **kwargs,
) -> list[Path]:
    import matplotlib.pyplot as plt

    # 1. Load data from each condition
    results = {}
    for label in labels:
        analysis_dir = Path(data[label]["analysis_dir"])
        result_file = analysis_dir / "binding_preference_aggregated_reps1-3.json"
        if result_file.exists():
            results[label] = AggregatedBindingPreferenceResult.load(result_file)

    if not results:
        return []

    # 2. Create figure
    fig, ax = plt.subplots(figsize=self.settings.contacts.figsize_enrichment_heatmap)

    # ... plotting logic ...

    # 3. Save and return
    output_path = self._get_output_path(output_dir, "my_plot_name")
    return [self._save_figure(fig, output_path)]
```

### Step 4: Access Plot Settings

The `self.settings` attribute provides access to `PlotSettings`:

```python
# Global settings
self.settings.dpi              # int, default 150
self.settings.format           # str, "png" or "svg"
self.settings.color_palette    # str, default "Set2"

# Analysis-specific settings
self.settings.contacts.figsize_enrichment_heatmap  # tuple
self.settings.contacts.enrichment_colormap         # str
self.settings.contacts.show_enrichment_error       # bool

self.settings.triad.figsize_kde_panel   # tuple
self.settings.triad.kde_fill_alpha      # float
```

To add new settings, extend the settings models in `compare/config.py`.

### Step 5: Handle Multiple Plot Files

If your plotter generates multiple files (e.g., one per polymer type), return all paths:

```python
def plot(self, data, labels, output_dir, **kwargs) -> list[Path]:
    output_paths = []

    for polymer_type in polymer_types:
        fig, ax = plt.subplots(...)
        # ... plot for this polymer type ...

        output_path = self._get_output_path(
            output_dir, f"binding_bars_{polymer_type.lower()}"
        )
        output_paths.append(self._save_figure(fig, output_path))

    return output_paths
```

## Complete Example: Bar Chart Plotter

Here's a complete example of a plotter that generates grouped bar charts:

```python
"""Example plotter for binding preference bar charts."""

from __future__ import annotations

import logging
from pathlib import Path
from typing import TYPE_CHECKING, Any, Sequence

import numpy as np

from polyzymd.compare.plotter import BasePlotter, PlotterRegistry

if TYPE_CHECKING:
    from polyzymd.analysis.contacts.binding_preference import (
        AggregatedBindingPreferenceResult,
    )
    from polyzymd.compare.config import ComparisonConfig

logger = logging.getLogger(__name__)


@PlotterRegistry.register("binding_preference_bars")
class BindingPreferenceBarPlotter(BasePlotter):
    """Generate grouped bar chart of binding preference enrichment.

    Creates a figure showing enrichment ratios as grouped bars with:
    - Groups: Protein groups (e.g., aromatic, polar, charged)
    - Bars within group: One per condition
    - Error bars: SEM across replicates
    - Reference line at 1.0 (neutral enrichment)

    One plot is generated per polymer type.
    """

    @classmethod
    def plot_type(cls) -> str:
        return "binding_preference_bars"

    def can_plot(self, comparison_config: "ComparisonConfig", analysis_type: str) -> bool:
        """Check if this plotter can handle the analysis type."""
        if analysis_type != "contacts":
            return False
        return self.settings.contacts.generate_enrichment_bars

    def plot(
        self,
        data: dict[str, Any],
        labels: Sequence[str],
        output_dir: Path,
        **kwargs,
    ) -> list[Path]:
        """Generate enrichment bar chart from filesystem data."""
        import matplotlib.pyplot as plt

        # Load binding preference results from each condition
        binding_results = self._load_binding_preference_results(data, labels)

        if not binding_results:
            logger.info("No binding preference data found - skipping bar plots")
            return []

        # Collect all polymer types and protein groups
        all_polymer_types: set[str] = set()
        all_protein_groups: set[str] = set()
        for result in binding_results.values():
            all_polymer_types.update(result.polymer_types())
            all_protein_groups.update(result.protein_groups())

        polymer_types = sorted(all_polymer_types)
        protein_groups = sorted(all_protein_groups)

        if not polymer_types or not protein_groups:
            return []

        # Generate one plot per polymer type
        output_paths: list[Path] = []
        valid_labels = [label for label in labels if label in binding_results]

        for poly_type in polymer_types:
            fig, ax = plt.subplots(
                figsize=self.settings.contacts.figsize_enrichment_bars,
                dpi=self.settings.dpi,
            )

            n_groups = len(protein_groups)
            n_conditions = len(valid_labels)
            bar_width = 0.8 / n_conditions
            x = np.arange(n_groups)

            for i, cond_label in enumerate(valid_labels):
                result = binding_results[cond_label]
                means = []
                sems = []

                for prot_group in protein_groups:
                    entry = result.get_entry(poly_type, prot_group)
                    if entry and entry.mean_enrichment is not None:
                        means.append(entry.mean_enrichment)
                        sems.append(entry.sem_enrichment or 0.0)
                    else:
                        means.append(0.0)
                        sems.append(0.0)

                offset = (i - n_conditions / 2 + 0.5) * bar_width
                ax.bar(x + offset, means, bar_width, yerr=sems, label=cond_label)

            ax.axhline(y=1.0, color="black", linestyle="--", label="Neutral")
            ax.set_xlabel("Protein Group")
            ax.set_ylabel("Enrichment Ratio")
            ax.set_title(f"Binding Preference: {poly_type}")
            ax.set_xticks(x)
            ax.set_xticklabels(protein_groups, rotation=45, ha="right")
            ax.legend()
            plt.tight_layout()

            output_path = self._get_output_path(
                output_dir, f"binding_preference_bars_{poly_type.lower()}"
            )
            output_paths.append(self._save_figure(fig, output_path))

        return output_paths

    def _load_binding_preference_results(
        self,
        data: dict[str, Any],
        labels: Sequence[str],
    ) -> dict[str, "AggregatedBindingPreferenceResult"]:
        """Load aggregated binding preference results for each condition."""
        from polyzymd.analysis.contacts.binding_preference import (
            AggregatedBindingPreferenceResult,
        )

        results: dict[str, AggregatedBindingPreferenceResult] = {}

        for label in labels:
            cond_data = data.get(label)
            if cond_data is None:
                continue

            analysis_dir = Path(cond_data.get("analysis_dir", ""))
            if not analysis_dir:
                continue

            # Find aggregated file
            agg_files = list(analysis_dir.glob("binding_preference_aggregated*.json"))
            if not agg_files:
                continue

            result_file = sorted(agg_files)[-1]  # Most recent
            try:
                results[label] = AggregatedBindingPreferenceResult.load(result_file)
            except Exception as e:
                logger.warning(f"Failed to load {result_file}: {e}")

        return results
```

## Testing Your Plotter

### Unit Test

```python
import pytest
from pathlib import Path
from unittest.mock import MagicMock, patch


class TestMyCustomPlotter:
    """Tests for MyCustomPlotter."""

    def test_plot_type_returns_registry_key(self):
        plotter = MyCustomPlotter(settings=MagicMock())
        assert plotter.plot_type() == "my_custom_plot"

    def test_can_plot_returns_true_for_correct_type(self):
        settings = MagicMock()
        settings.my_analysis.generate_custom_plot = True
        plotter = MyCustomPlotter(settings=settings)

        assert plotter.can_plot(MagicMock(), "my_analysis") is True
        assert plotter.can_plot(MagicMock(), "other_type") is False

    def test_plot_returns_empty_when_no_data(self):
        plotter = MyCustomPlotter(settings=MagicMock())
        result = plotter.plot({}, [], Path("/tmp"))
        assert result == []
```

### Integration Test

```bash
# Run plot-all to test discovery and execution
mamba run -n polyzymd-env polyzymd compare plot-all -f comparison.yaml
```

## Common Patterns

### Pattern: Load Aggregated Results

Most plotters load aggregated (replicate-averaged) results:

```python
def _load_aggregated(self, data, labels):
    from my_module import MyAggregatedResult

    results = {}
    for label in labels:
        agg_dir = Path(data[label].get("aggregated_dir", ""))
        result_file = agg_dir / "my_aggregated.json"
        if result_file.exists():
            results[label] = MyAggregatedResult.load(result_file)
    return results
```

### Pattern: Pool Per-Replicate Data

For distribution plots (KDEs, histograms), pool raw data across replicates:

```python
def _pool_replicate_data(self, data, labels):
    pooled = {}
    for label in labels:
        analysis_dir = Path(data[label]["analysis_dir"])
        replicates = data[label]["replicates"]

        all_values = []
        for rep in replicates:
            rep_file = analysis_dir / f"run_{rep}" / "result.json"
            if rep_file.exists():
                result = MyResult.load(rep_file)
                all_values.extend(result.values)

        if all_values:
            pooled[label] = np.array(all_values)

    return pooled
```

### Pattern: Conditional Plot Generation

Only generate plots when data is available:

```python
def plot(self, data, labels, output_dir, **kwargs):
    results = self._load_results(data, labels)

    if not results:
        logger.info("No data found - skipping plot")
        return []  # Return empty list, not raise exception

    # Continue with plotting...
```

## Adding Plot Settings

To add configurable settings for your plotter:

### 1. Extend Settings Model

In `compare/config.py`:

```python
class PlotSettingsMyAnalysis(BaseModel):
    """Plot settings for my_analysis type."""

    generate_custom_plot: bool = True
    figsize_custom: tuple[int, int] = (10, 6)
    custom_colormap: str = "viridis"


class PlotSettings(BaseModel):
    """Root plot settings."""

    # Existing fields...
    my_analysis: PlotSettingsMyAnalysis = Field(default_factory=PlotSettingsMyAnalysis)
```

### 2. Access in Plotter

```python
def plot(self, data, labels, output_dir, **kwargs):
    fig, ax = plt.subplots(figsize=self.settings.my_analysis.figsize_custom)
    ax.imshow(matrix, cmap=self.settings.my_analysis.custom_colormap)
```

## Troubleshooting

### "No data found" but data exists

Check that your file pattern matches:
```python
# Verify file exists
agg_files = list(analysis_dir.glob("binding_preference_aggregated*.json"))
print(f"Found files: {agg_files}")  # Debug
```

### Plot not appearing in `plot-all`

1. Verify registration:
   ```python
   from polyzymd.compare.plotter import PlotterRegistry
   print(PlotterRegistry.list())  # Should include your key
   ```

2. Check `can_plot()` returns `True`:
   ```python
   assert plotter.can_plot(config, "my_analysis") is True
   ```

3. Verify settings enable the plot:
   ```yaml
   plot_settings:
     my_analysis:
       generate_custom_plot: true
   ```

### Type errors in matplotlib

Use tuples for `add_axes` and `tight_layout`:
```python
# Correct
cbar_ax = fig.add_axes((0.92, 0.15, 0.02, 0.7))  # tuple
plt.tight_layout(rect=(0, 0, 0.9, 0.95))         # tuple

# Incorrect (causes type errors)
cbar_ax = fig.add_axes([0.92, 0.15, 0.02, 0.7])  # list
```

## See Also

- {doc}`extending_comparators` — Creating custom statistical comparators
- {doc}`analysis_binding_preference` — Binding preference analysis details
- {doc}`architecture` — Overall system architecture
