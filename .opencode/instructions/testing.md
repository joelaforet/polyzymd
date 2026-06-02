# Testing Rules

## Current State

The test suite covers the full analysis plugin system and core infrastructure:

- **Test directory:** `tests/` with subdirectories mirroring the source tree
- **Test count:** 1,470 collected (1,464 passed, 6 skipped)
- **Fixtures:** `tests/conftest.py` with shared fixtures for comparison configs, mock data, etc.
- **Markers:** `@pytest.mark.slow` for tests requiring simulation data

### Directory Structure

```
tests/
├── conftest.py                  # Shared fixtures
├── _support/                    # Test utilities (analysis_testkit.py)
├── analyses/                    # Tests for analyses/ source tree
│   ├── test_base.py             # analyses/base.py
│   ├── test_discovery.py        # analyses/discovery.py
│   ├── test_orchestrator.py     # analyses/orchestrator.py
│   ├── test_orchestrator_workers.py
│   ├── test_orchestrator_cost_hints.py
│   ├── test_stats.py            # analyses/stats.py
│   ├── shared/                  # analyses/shared/ utilities
│   │   ├── test_convergence.py
│   │   ├── test_defaults.py
│   │   ├── test_inferential_statistics.py
│   │   ├── test_multi_run_comparison.py
│   │   ├── test_multi_run_formatting.py
│   │   ├── test_paths.py
│   │   ├── test_result_io.py
│   │   └── test_sasa.py
│   ├── plugins/                 # One file per analysis plugin
│   │   ├── test_catalytic_triad.py
│   │   ├── test_contacts.py
│   │   ├── test_distances.py
│   │   ├── test_hydrogen_bonds.py
│   │   ├── test_rg.py
│   │   ├── test_rmsd.py
│   │   ├── test_rmsf.py
│   │   ├── test_sasa.py
│   │   └── test_secondary_structure.py
│   └── integration/             # Cross-plugin integration tests
│       ├── test_fdr_plugin_wiring.py
│       └── test_zero_control_regression.py
├── cli/                         # Tests for cli/ source tree
│   ├── test_main.py
│   ├── test_main_recover.py
│   ├── test_main_status.py
│   ├── test_compare.py
│   ├── test_scaffold.py
│   └── test_colors.py
├── config/                      # Tests for config/ source tree
│   ├── test_schema.py
│   └── test_loader.py
├── simulation/                  # Tests for simulation/ source tree
│   ├── test_continuation.py
│   ├── test_progress.py
│   ├── test_runner.py
│   └── test_signals.py
├── workflow/                    # Tests for workflow/ source tree
│   ├── test_analysis_slurm.py
│   ├── test_slurm.py
│   └── test_daisy_chain.py
├── exporters/                   # Tests for exporters/ source tree
│   ├── test_gromacs.py
│   └── test_interchange.py
└── utils/                       # Tests for utils/ source tree
    ├── test_packmol.py
    └── test_replicates.py
```

## Running Tests

```bash
# Full test suite
pixi run -e build pytest tests/ -v

# Specific test file
pixi run -e build pytest tests/analyses/plugins/test_rmsf.py -v

# Run tests matching a pattern
pixi run -e build pytest tests/ -v -k "rmsf"

# Run tests for a specific plugin
pixi run -e build pytest tests/ -v -k "secondary_structure"
```

## Writing New Tests

When adding tests:

1. Place tests in the subdirectory matching the source module (e.g., `tests/analyses/plugins/`)
2. Name files `test_<source_module>.py` with 1:1 correspondence to source files
3. Use pytest conventions (`test_` prefix for functions/methods)
4. Mock heavy dependencies (OpenMM, MDAnalysis) for unit tests
5. Use `@pytest.mark.slow` for tests requiring simulation data
6. Use existing fixtures from `tests/conftest.py`

### Testing an Analysis Plugin

When adding a new analysis plugin in `analyses/`, write tests that cover:

1. **Discovery**: Plugin is found by `list_analyses()` and `get_analysis()`
2. **Class variables**: `name` and `Settings` are set correctly
3. **Settings validation**: Pydantic model validates/rejects correctly
4. **MDAnalysis job stage**: `build_mda_jobs()` constructs `MDAAnalysisJob`
   objects using `FrameSelection` and fake `AnalysisBase`-compatible work.
5. **Collector artifacts**: `build_mda_collector()` maps completed jobs to a
   valid `ReplicateArtifact` without serializing raw MDAnalysis `Results`.
6. **Artifact aggregation**: `aggregate()` or the default artifact aggregator
   combines replicate artifacts into a `ConditionArtifact` without loading
   trajectories.
7. **Comparison metrics**: `extract_metrics()` or custom `compare()` consumes
   condition artifacts and uses replicate-level statistics.
8. **Artifact-only plots**: `plot()` reads cached artifacts/sidecars only.
9. **format()**: Generates readable CLI output from comparison artifacts/results.

Example test structure for a plugin:

```python
from pathlib import Path
from unittest.mock import MagicMock

from polyzymd.analyses import get_analysis, list_analyses
from polyzymd.analyses.base import (
    AggregateContext,
    Condition,
    MDAReplicateJobContext,
    MetricValue,
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


class TestMyPluginMDAJobs:
    """Test the MDAnalysis job stage with small fakes."""

    class FakeTrajectory:
        def __len__(self) -> int:
            return 50

    class FakeUniverse:
        trajectory = FakeTrajectory()

    def test_builds_mda_jobs(self, tmp_path):
        cls = get_analysis("my_analysis")
        analysis = cls()
        condition = Condition(
            label="Test",
            config_path=Path("/fake/config.yaml"),
            replicates=(1,),
            sim_config=object(),
        )
        ctx = MDAReplicateJobContext(
            condition=condition,
            replicate=1,
            sim_config=condition.sim_config,
            output_dir=tmp_path / "run_1",
            equilibration="0ns",
            recompute=True,
            settings=cls.Settings(),
        )

        jobs = analysis.build_mda_jobs(ctx)
        assert jobs
        assert all(job.name for job in jobs)


class TestMyPluginArtifacts:
    """Test artifact aggregation — no trajectory mocks needed."""

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
