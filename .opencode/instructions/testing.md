# Testing Rules

## Current State

The test suite covers the full analysis plugin system and core infrastructure:

- **Test directory:** `tests/` with comprehensive coverage of `analyses/`, `compare/`, and `config/`
- **Test count:** 900+ tests passing
- **Fixtures:** `tests/conftest.py` with shared fixtures for comparison configs, mock data, etc.
- **Markers:** `@pytest.mark.slow` for tests requiring simulation data

## Running Tests

```bash
# Full test suite
pixi run -e build pytest tests/ -v

# Specific test file
pixi run -e build pytest tests/test_rmsf_plugin.py -v

# Run tests matching a pattern
pixi run -e build pytest tests/ -v -k "rmsf"

# Run tests for a specific plugin
pixi run -e build pytest tests/ -v -k "secondary_structure"
```

## Writing New Tests

When adding tests:

1. Place unit tests in `tests/` directory
2. Name files `test_<name>_plugin.py` (e.g., `tests/test_rg_plugin.py`)
3. Use pytest conventions (`test_` prefix for functions/methods)
4. Mock heavy dependencies (OpenMM, MDAnalysis) for unit tests
5. Use `@pytest.mark.slow` for tests requiring simulation data
6. Use existing fixtures from `tests/conftest.py`

### Testing an Analysis Plugin

When adding a new analysis plugin in `analyses/`, write tests that cover:

1. **Discovery**: Plugin is found by `list_analyses()` and `get_analysis()`
2. **Class variables**: `name` and `Settings` are set correctly
3. **Settings validation**: Pydantic model validates/rejects correctly
4. **compute_replicate()**: Returns expected structure (mock TrajectoryLoader / MDAnalysis)
5. **aggregate()**: Combines replicate results correctly (no mocks needed)
6. **extract_metrics()**: Returns correct `MetricValue` instances (if using default compare)
7. **_deserialize_result()**: Loads JSON back correctly (if using default compare)
8. **compare()**: Produces a valid `ComparisonResult`
9. **format()**: Generates readable CLI output

Example test structure for a plugin:

```python
import json
from pathlib import Path
from unittest.mock import MagicMock, patch

import numpy as np
import pytest

from polyzymd.analyses import get_analysis, list_analyses
from polyzymd.analyses.base import (
    AggregateContext,
    Condition,
    MetricValue,
    ReplicateContext,
)


class TestMyPluginDiscovery:
    """Tests for discovery and class-level attributes."""

    def test_discovered(self):
        """Plugin should be discovered automatically."""
        analyses = list_analyses()
        assert "my_analysis" in analyses

    def test_class_variables(self):
        cls = get_analysis("my_analysis")
        assert cls.name == "my_analysis"
        assert hasattr(cls, "Settings")

    def test_settings_defaults(self):
        cls = get_analysis("my_analysis")
        settings = cls.Settings()
        assert settings.selection == "protein and name CA"


class TestMyPluginCompute:
    """Test compute_replicate with mocked trajectories."""

    @patch("polyzymd.analyses.my_analysis.TrajectoryLoader")
    def test_computes_metric(self, MockLoader, tmp_path):
        cls = get_analysis("my_analysis")
        analysis = cls()
        settings = cls.Settings()

        # Mock TrajectoryLoader and Universe
        mock_loader = MagicMock()
        MockLoader.return_value = mock_loader
        mock_universe = MagicMock()
        mock_atoms = MagicMock()
        mock_atoms.__len__ = MagicMock(return_value=100)
        mock_universe.select_atoms.return_value = mock_atoms
        mock_trajectory = MagicMock()
        mock_trajectory.__len__ = MagicMock(return_value=50)
        mock_trajectory.__getitem__ = MagicMock(return_value=range(50))
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
        assert isinstance(result, dict)


class TestMyPluginAggregate:
    """Test aggregation — no mocks needed."""

    def test_aggregate(self, tmp_path):
        cls = get_analysis("my_analysis")
        analysis = cls()
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
            settings=cls.Settings(),
        )

        results = [
            {"my_metric": 15.0, "replicate": 1},
            {"my_metric": 15.5, "replicate": 2},
            {"my_metric": 14.8, "replicate": 3},
        ]

        agg = analysis.aggregate(ctx, results)
        assert "replicate_values" in agg


class TestMyPluginMetrics:
    """Test metric extraction for default comparison."""

    def test_extract_metrics(self):
        cls = get_analysis("my_analysis")
        analysis = cls()
        summary = {"mean_value": 1.5, "sem_value": 0.1, "replicate_values": [1.4, 1.6]}
        metrics = analysis.extract_metrics(summary)
        assert isinstance(metrics, dict)
        for v in metrics.values():
            assert isinstance(v, MetricValue)
```

## Test Data

- Real simulation data lives at `../testing_analysis/` (outside repo)
- Do NOT commit trajectory files (.xtc, .dcd) or large PDB files to the repo
- For unit tests, use small synthetic data or mock objects
- For integration tests, use the `@pytest.mark.slow` marker and document
  the required data paths
