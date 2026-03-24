# Extending the Comparison Framework

This guide shows developers how to create custom comparators and analyzers using PolyzyMD's object-oriented framework. The framework follows the **Open-Closed Principle** (open for extension, closed for modification) and uses the **Template Method** design pattern for DRY comparison logic.

```{note}
This guide covers the **compare** module for cross-condition statistical comparisons.
For single-condition analysis, see the analysis module documentation.
```

## Overview

The comparison framework provides:

- **BaseComparator**: Abstract base class implementing the Template Method pattern
- **ComparatorRegistry**: Registry for auto-discovery of comparator types
- **MetricType**: Classification system for correct autocorrelation handling
- **Shared statistical methods**: T-tests, ANOVA, Cohen's d (implemented once)

### Architecture

```
compare/
├── core/
│   ├── base.py          # BaseComparator, BaseConditionSummary, BaseComparisonResult
│   ├── registry.py      # ComparatorRegistry with @register decorator
│   └── __init__.py
├── comparators/
│   ├── rmsf.py          # RMSFComparator (VARIANCE_BASED)
│   ├── triad.py         # TriadComparator (MEAN_BASED)
│   ├── contacts.py      # ContactsComparator (MEAN_BASED)
│   └── __init__.py
├── results/
│   ├── rmsf.py          # RMSFConditionSummary, RMSFComparisonResult
│   ├── triad.py         # TriadConditionSummary, TriadComparisonResult
│   ├── contacts.py      # ContactsConditionSummary, ContactsComparisonResult
│   └── __init__.py
├── settings.py          # Pydantic settings for each analysis type
├── config.py            # ComparisonConfig, ConditionConfig
├── statistics.py        # Shared statistical functions
└── cli.py               # CLI commands
```

## Quick Start: Minimal Comparator

Here's the minimal code to create a new comparator:

```python
from polyzymd.analysis.core.metric_type import MetricType
from polyzymd.compare.core.base import BaseComparator, ANOVASummary
from polyzymd.compare.core.registry import ComparatorRegistry

# Step 1: Define your result models (see "Result Models" section)
from my_module.results import MyConditionSummary, MyComparisonResult


@ComparatorRegistry.register("my_metric")
class MyComparator(BaseComparator[MySettings, dict, MyConditionSummary, MyComparisonResult]):
    """Compare my_metric across simulation conditions."""

    @classmethod
    def comparison_type_name(cls) -> str:
        return "my_metric"

    @property
    def metric_type(self) -> MetricType:
        # Choose based on what you measure (see MetricType section)
        return MetricType.MEAN_BASED

    def _load_or_compute(self, cond, recompute) -> dict:
        # Load cached results or compute analysis
        ...

    def _build_condition_summary(self, cond, data) -> MyConditionSummary:
        # Convert raw data to structured summary
        ...

    def _build_result(self, summaries, comparisons, anova, ranking, effective_control, excluded) -> MyComparisonResult:
        # Construct final result object
        ...

    def _get_replicate_values(self, summary) -> list[float]:
        # Extract per-replicate values for statistical tests
        return summary.replicate_values

    def _get_mean_value(self, summary) -> float:
        # Get the mean value for percent change calculations
        return summary.mean_value

    def _interpret_direction(self, percent_change) -> str:
        # Return "improving", "worsening", "stabilizing", etc.
        if percent_change > 0:
            return "increasing"
        elif percent_change < 0:
            return "decreasing"
        return "unchanged"

    def _rank_summaries(self, summaries) -> list:
        # Sort summaries (best first)
        return sorted(summaries, key=lambda s: s.mean_value, reverse=True)
```

## Step-by-Step Guide

### Step 1: Choose Your MetricType

The `metric_type` property determines how autocorrelation is handled in the underlying analysis. This is **critical** for correct uncertainty quantification.

#### MEAN_BASED Metrics

Use for metrics that compute **averages** over frames:

- Contact fraction (average fraction of frames in contact)
- Average distance
- Average angle
- Catalytic triad proximity

**Frame strategy**: Use ALL frames for computation (maximizes precision)

**Uncertainty strategy**: Correct SEM using N_eff = N/g where g is the statistical inefficiency

```python
@property
def metric_type(self) -> MetricType:
    return MetricType.MEAN_BASED
```

#### VARIANCE_BASED Metrics

Use for metrics that compute **fluctuations** or variances:

- RMSF (root-mean-square fluctuation)
- Standard deviation of distances
- B-factors

**Frame strategy**: Subsample to independent frames separated by 2τ (correlation time)

**Uncertainty strategy**: Use standard formula on subsampled (independent) data

```python
@property
def metric_type(self) -> MetricType:
    return MetricType.VARIANCE_BASED
```

```{warning}
Declaring the wrong MetricType leads to incorrect uncertainty estimates!
- MEAN_BASED with variance metric → underestimated uncertainty
- VARIANCE_BASED with mean metric → unnecessarily large uncertainty (discards data)
```

**Reference**: Grossfield et al. (2018) LiveCoMS 1:5067 — Best Practices for Quantification of Uncertainty

### Step 2: Define Result Models

Create Pydantic models for your condition summary and comparison result.

#### Condition Summary

Inherit from `BaseConditionSummary` and add metric-specific fields:

```python
from pydantic import Field
from polyzymd.compare.core.base import BaseConditionSummary


class MyConditionSummary(BaseConditionSummary):
    """Summary statistics for one condition in a MyMetric comparison.

    Attributes
    ----------
    label : str
        Display name for this condition (inherited).
    config_path : str
        Path to the simulation config file (inherited).
    n_replicates : int
        Number of replicates included (inherited).
    replicate_values : list[float]
        Per-replicate values of the primary metric (inherited).
    mean_value : float
        Mean of my metric across replicates.
    sem_value : float
        Standard error of the mean.
    """

    mean_value: float = Field(..., description="Mean of my metric")
    sem_value: float = Field(..., description="Standard error of mean")

    @property
    def primary_metric_value(self) -> float:
        """Return the primary metric value for ranking."""
        return self.mean_value

    @property
    def primary_metric_sem(self) -> float:
        """Return the SEM of the primary metric."""
        return self.sem_value
```

#### Comparison Result

Inherit from `BaseComparisonResult` and add analysis-specific fields:

```python
from datetime import datetime
from typing import ClassVar
from polyzymd.compare.core.base import (
    BaseComparisonResult,
    PairwiseComparison,
    ANOVASummary,
)


class MyComparisonResult(BaseComparisonResult[MyConditionSummary, PairwiseComparison]):
    """Complete comparison result for MyMetric.

    Attributes
    ----------
    metric : str
        Always "my_metric".
    name : str
        Name of the comparison project.
    conditions : list[MyConditionSummary]
        Summary for each condition.
    pairwise_comparisons : list[PairwiseComparison]
        Statistical comparisons.
    anova : ANOVASummary, optional
        ANOVA result if 3+ conditions.
    ranking : list[str]
        Labels sorted by metric (best first).
    my_custom_field : str
        Any additional fields specific to your analysis.
    """

    comparison_type: ClassVar[str] = "my_metric"

    # Override with specific types
    conditions: list[MyConditionSummary]
    pairwise_comparisons: list[PairwiseComparison]
    anova: ANOVASummary | None = None

    # Add your custom fields
    my_custom_field: str = Field(..., description="Custom field for my analysis")
```

### Step 3: Implement Abstract Methods

The `BaseComparator` class requires these abstract methods:

| Method | Purpose | Return Type |
|--------|---------|-------------|
| `comparison_type_name()` | Registry key (e.g., "rmsf") | `str` |
| `metric_type` | Autocorrelation strategy | `MetricType` |
| `_load_or_compute()` | Load cached or compute analysis | Your data dict |
| `_build_condition_summary()` | Convert raw data to summary | Your summary type |
| `_build_result()` | Construct final result | Your result type |
| `_get_replicate_values()` | Extract values for stats | `list[float]` |
| `_get_mean_value()` | Get mean for percent change | `float` |
| `_interpret_direction()` | Interpret change direction | `str` |
| `_rank_summaries()` | Sort summaries (best first) | `list[Summary]` |

#### Example: `_load_or_compute()`

This method should load cached results if available, or compute the analysis:

```python
def _load_or_compute(
    self,
    cond: "ConditionConfig",
    recompute: bool,
) -> dict:
    """Load existing results or compute analysis.

    Parameters
    ----------
    cond : ConditionConfig
        Condition to analyze (has label, config path, replicates).
    recompute : bool
        If True, ignore cache and recompute.

    Returns
    -------
    dict
        Raw analysis data with keys needed by _build_condition_summary.
    """
    from polyzymd.config.schema import SimulationConfig

    logger.info(f"Processing condition: {cond.label}")

    # Load simulation config
    sim_config = SimulationConfig.from_yaml(cond.config)

    # Try to find cached result
    result_path = self._find_cached_result(sim_config, cond.replicates)

    if result_path and result_path.exists() and not recompute:
        logger.info(f"  Loading cached result: {result_path}")
        cached = MyResult.load(result_path)
        return {
            "mean_value": cached.mean,
            "sem_value": cached.sem,
            "n_replicates": cached.n_replicates,
            "replicate_values": cached.per_replicate_values,
        }

    # Compute analysis
    logger.info(f"  Computing for replicates {cond.replicates}...")
    result = self._compute_analysis(sim_config, cond.replicates)

    return {
        "mean_value": result.mean,
        "sem_value": result.sem,
        "n_replicates": result.n_replicates,
        "replicate_values": result.per_replicate_values,
    }
```

### Step 4: Register with the Registry

The `@ComparatorRegistry.register()` decorator automatically registers your comparator:

```python
@ComparatorRegistry.register("my_metric")
class MyComparator(BaseComparator[...]):
    ...
```

After registration, your comparator is available via:

```python
# List all registered comparators
ComparatorRegistry.list_available()
# ['contacts', 'my_metric', 'rmsf', 'triad']

# Create instance via registry
comparator = ComparatorRegistry.create(
    "my_metric",
    config=comparison_config,
    analysis_settings=my_settings,
)

# Run comparison
result = comparator.compare()
```

### Step 5: Add to CLI (Optional)

To make your comparator available via the CLI, add it to `compare/cli.py`:

```python
@compare.command("my-metric")
@click.option("-c", "--config", required=True, type=click.Path(exists=True))
@click.option("--recompute", is_flag=True, help="Force recompute")
def compare_my_metric(config: str, recompute: bool):
    """Compare my_metric across conditions."""
    from polyzymd.compare.config import ComparisonConfig
    from polyzymd.compare.comparators.my_metric import MyComparator

    cfg = ComparisonConfig.from_yaml(config)
    settings = cfg.analysis_settings.get("my_metric")

    comparator = MyComparator(cfg, settings)
    result = comparator.compare(recompute=recompute)

    # Print results...
```

Or use the generic command (no code changes needed):

```bash
polyzymd compare run my_metric -c comparison.yaml
```

## Advanced Patterns

### Custom Pairwise Comparison

If your metric requires custom pairwise comparison logic (e.g., multiple sub-metrics), override `_compare_pair()`:

```python
def _compare_pair(
    self,
    cond_a: MyConditionSummary,
    cond_b: MyConditionSummary,
) -> MyPairwiseComparison:
    """Compare two conditions with custom logic."""
    # Get replicate values
    values_a = self._get_replicate_values(cond_a)
    values_b = self._get_replicate_values(cond_b)

    # Custom statistical tests
    from polyzymd.compare.statistics import independent_ttest, cohens_d

    ttest = independent_ttest(values_a, values_b)
    effect = cohens_d(values_a, values_b)

    # Return your custom pairwise comparison type
    return MyPairwiseComparison(
        condition_a=cond_a.label,
        condition_b=cond_b.label,
        custom_metric=compute_custom_metric(values_a, values_b),
        t_statistic=ttest.t_statistic,
        p_value=ttest.p_value,
        cohens_d=effect.cohens_d,
        ...
    )
```

### Filtering Conditions

Override `_filter_conditions()` to exclude certain conditions (e.g., "No Polymer" for contacts):

```python
def _filter_conditions(
    self,
) -> tuple[list["ConditionConfig"], list["ConditionConfig"]]:
    """Filter conditions before analysis.

    Returns
    -------
    tuple[list[ConditionConfig], list[ConditionConfig]]
        (valid_conditions, excluded_conditions)
    """
    valid = []
    excluded = []

    for cond in self.config.conditions:
        if self._should_include(cond):
            valid.append(cond)
        else:
            excluded.append(cond)
            logger.info(f"Excluding '{cond.label}': reason...")

    return valid, excluded
```

### Multiple ANOVA Metrics

If your analysis has multiple metrics (like contacts with coverage and contact_fraction), override the entire `compare()` method:

```python
def compare(self, recompute: bool = False) -> MyComparisonResult:
    """Run comparison with custom multi-metric logic."""
    # ... custom implementation
    # See ContactsComparator for a full example
```

### Custom Cohen's d Interpretation

For RMSF, negative Cohen's d indicates stabilization (desirable). Override the hook:

```python
def _use_rmsf_mode_for_cohens_d(self) -> bool:
    """Use RMSF-specific Cohen's d interpretation.

    When True, negative d is interpreted as "stabilizing" rather than
    the default "decreasing".
    """
    return True  # Only RMSFComparator returns True
```

## Complete Example: GyrationComparator

Here's a complete example implementing a radius of gyration comparator:

```python
"""Radius of gyration comparator."""

from __future__ import annotations

import logging
from datetime import datetime
from pathlib import Path
from typing import TYPE_CHECKING, Any, ClassVar

from pydantic import Field

from polyzymd import __version__
from polyzymd.analysis.core.metric_type import MetricType
from polyzymd.compare.core.base import (
    ANOVASummary,
    BaseComparator,
    BaseComparisonResult,
    BaseConditionSummary,
    PairwiseComparison,
)
from polyzymd.compare.core.registry import ComparatorRegistry

if TYPE_CHECKING:
    from polyzymd.compare.config import ComparisonConfig, ConditionConfig

logger = logging.getLogger("polyzymd.compare")


# =============================================================================
# Result Models
# =============================================================================


class GyrationConditionSummary(BaseConditionSummary):
    """Summary for one condition in a gyration comparison."""

    mean_rg: float = Field(..., description="Mean radius of gyration (Angstroms)")
    sem_rg: float = Field(..., description="Standard error of mean Rg")

    @property
    def primary_metric_value(self) -> float:
        return self.mean_rg

    @property
    def primary_metric_sem(self) -> float:
        return self.sem_rg


class GyrationComparisonResult(
    BaseComparisonResult[GyrationConditionSummary, PairwiseComparison]
):
    """Complete gyration comparison result."""

    comparison_type: ClassVar[str] = "gyration"

    conditions: list[GyrationConditionSummary]
    pairwise_comparisons: list[PairwiseComparison]
    anova: ANOVASummary | None = None
    selection: str = Field(..., description="Atom selection used")


# =============================================================================
# Settings
# =============================================================================


class GyrationAnalysisSettings:
    """Settings for gyration analysis."""

    def __init__(self, selection: str = "protein and name CA"):
        self.selection = selection


# =============================================================================
# Comparator
# =============================================================================


@ComparatorRegistry.register("gyration")
class GyrationComparator(
    BaseComparator[
        GyrationAnalysisSettings,
        dict[str, Any],
        GyrationConditionSummary,
        GyrationComparisonResult,
    ]
):
    """Compare radius of gyration across simulation conditions.

    Lower Rg indicates a more compact structure.

    Parameters
    ----------
    config : ComparisonConfig
        Comparison configuration.
    analysis_settings : GyrationAnalysisSettings
        Gyration analysis settings.
    equilibration : str, optional
        Equilibration time override.
    """

    comparison_type: ClassVar[str] = "gyration"

    @classmethod
    def comparison_type_name(cls) -> str:
        return "gyration"

    @property
    def metric_type(self) -> MetricType:
        """Rg is a mean-based metric (average over frames)."""
        return MetricType.MEAN_BASED

    def _load_or_compute(
        self,
        cond: "ConditionConfig",
        recompute: bool,
    ) -> dict[str, Any]:
        """Load or compute Rg analysis."""
        import numpy as np

        from polyzymd.config.schema import SimulationConfig

        logger.info(f"Processing condition: {cond.label}")
        sim_config = SimulationConfig.from_yaml(cond.config)

        # Compute Rg for each replicate
        rg_values = []
        for rep in cond.replicates:
            rg = self._compute_rg_for_replicate(sim_config, rep)
            if rg is not None:
                rg_values.append(rg)

        return {
            "mean_rg": float(np.mean(rg_values)),
            "sem_rg": float(np.std(rg_values, ddof=1) / np.sqrt(len(rg_values))),
            "n_replicates": len(rg_values),
            "replicate_values": rg_values,
        }

    def _compute_rg_for_replicate(self, sim_config, replicate: int) -> float | None:
        """Compute mean Rg for a single replicate."""
        import MDAnalysis as mda
        import numpy as np

        run_dir = sim_config.get_working_directory(replicate)
        topology = run_dir / "solvated_system.pdb"

        if not topology.exists():
            logger.warning(f"  Replicate {replicate} not found")
            return None

        # Find trajectory files
        import glob

        trajs = sorted(glob.glob(str(run_dir / "production_*" / "*_trajectory.dcd")))
        if not trajs:
            return None

        u = mda.Universe(str(topology), trajs)
        atoms = u.select_atoms(self.analysis_settings.selection)

        rg_per_frame = []
        for ts in u.trajectory:
            rg_per_frame.append(atoms.radius_of_gyration())

        return float(np.mean(rg_per_frame))

    def _build_condition_summary(
        self,
        cond: "ConditionConfig",
        data: dict[str, Any],
    ) -> GyrationConditionSummary:
        return GyrationConditionSummary(
            label=cond.label,
            config_path=str(cond.config),
            n_replicates=data["n_replicates"],
            mean_rg=data["mean_rg"],
            sem_rg=data["sem_rg"],
            replicate_values=data["replicate_values"],
        )

    def _build_result(
        self,
        summaries: list[GyrationConditionSummary],
        comparisons: list[Any],
        anova: ANOVASummary | None,
        ranking: list[str],
        effective_control: str | None,
        excluded_conditions: list["ConditionConfig"],
    ) -> GyrationComparisonResult:
        return GyrationComparisonResult(
            metric="gyration",
            name=self.config.name,
            control_label=effective_control,
            conditions=summaries,
            pairwise_comparisons=comparisons,
            anova=anova,
            ranking=ranking,
            equilibration_time=self.equilibration,
            selection=self.analysis_settings.selection,
            created_at=datetime.now(),
            polyzymd_version=__version__,
        )

    def _get_replicate_values(self, summary: GyrationConditionSummary) -> list[float]:
        return summary.replicate_values

    def _get_mean_value(self, summary: GyrationConditionSummary) -> float:
        return summary.mean_rg

    def _interpret_direction(self, percent_change: float) -> str:
        # Lower Rg = more compact = potentially more stable
        if percent_change < 0:
            return "compacting"
        elif percent_change > 0:
            return "expanding"
        return "unchanged"

    def _rank_summaries(
        self, summaries: list[GyrationConditionSummary]
    ) -> list[GyrationConditionSummary]:
        # Lower Rg is better (more compact)
        return sorted(summaries, key=lambda s: s.mean_rg)
```

## Testing Your Comparator

```python
import pytest
from unittest.mock import MagicMock, patch


class TestGyrationComparator:
    """Tests for GyrationComparator."""

    def test_comparison_type_name(self):
        """comparison_type_name should return 'gyration'."""
        assert GyrationComparator.comparison_type_name() == "gyration"

    def test_metric_type_is_mean_based(self):
        """Gyration should be MEAN_BASED."""
        comparator = GyrationComparator(
            config=MagicMock(),
            analysis_settings=GyrationAnalysisSettings(),
        )
        assert comparator.metric_type == MetricType.MEAN_BASED

    def test_interpret_direction_compacting(self):
        """Negative percent change should be 'compacting'."""
        comparator = GyrationComparator(
            config=MagicMock(),
            analysis_settings=GyrationAnalysisSettings(),
        )
        assert comparator._interpret_direction(-5.0) == "compacting"

    def test_registered_in_registry(self):
        """Comparator should be auto-registered."""
        from polyzymd.compare.core.registry import ComparatorRegistry

        assert ComparatorRegistry.is_registered("gyration")
```

## Summary Checklist

When creating a new comparator:

- [ ] Choose the correct `MetricType` (MEAN_BASED or VARIANCE_BASED)
- [ ] Create `ConditionSummary` model inheriting from `BaseConditionSummary`
- [ ] Create `ComparisonResult` model inheriting from `BaseComparisonResult`
- [ ] Implement all abstract methods
- [ ] Add `@ComparatorRegistry.register()` decorator
- [ ] Write unit tests
- [ ] (Optional) Add CLI command or use generic `compare run <type>`

## See Also

- {doc}`analysis_statistics_best_practices` — Autocorrelation and uncertainty
- {doc}`analysis_compare_conditions` — User guide for running comparisons
- {doc}`architecture` — Overall project architecture
