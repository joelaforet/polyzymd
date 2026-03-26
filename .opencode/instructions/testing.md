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
pixi run -e build pytest tests/test_analyses_rmsf.py -v

# Run tests matching a pattern
pixi run -e build pytest tests/ -v -k "rmsf"

# Run tests for a specific plugin
pixi run -e build pytest tests/ -v -k "secondary_structure"
```

## Writing New Tests

When adding tests:

1. Place unit tests in `tests/` directory
2. Name files `test_<module>.py` (e.g., `tests/test_analyses_my_plugin.py`)
3. Use pytest conventions (`test_` prefix for functions/methods)
4. Mock heavy dependencies (OpenMM, MDAnalysis) for unit tests
5. Use `@pytest.mark.slow` for tests requiring simulation data
6. Use existing fixtures from `tests/conftest.py`

### Testing an Analysis Plugin

When adding a new analysis plugin in `analyses/`, write tests that cover:

1. **Discovery**: Plugin is found by `list_analyses()` and `get_analysis()`
2. **Class variables**: `name` and `Settings` are set correctly
3. **Settings validation**: Pydantic model validates/rejects correctly
4. **compute_replicate()**: Returns expected structure (mock MDAnalysis)
5. **aggregate()**: Combines replicate results correctly
6. **extract_metrics()**: Returns correct `MetricValue` instances (if using default compare)
7. **compare()**: Produces a valid comparison result
8. **format()**: Generates readable CLI output

Example test structure for a plugin:

```python
import pytest
from unittest.mock import MagicMock, patch

from polyzymd.analyses import get_analysis, list_analyses
from polyzymd.analyses.base import MetricValue


class TestMyPlugin:
    """Tests for my_analysis plugin."""

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

    def test_extract_metrics_returns_metric_values(self):
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
