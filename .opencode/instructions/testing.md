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

1. **What compute() measures**: run it on the `synthetic_universe` fixture and
   assert the known answer, the frame count over the production window, and the
   unit of every observable.
2. **What the condition reports**: run `run_contract_analysis` over three
   replicates and assert on the `ObservableAggregate` fields, not on the raw
   observable.

That is the whole plugin test surface. Discovery, settings validation,
persistence, aggregation, hypothesis testing and plotting belong to the
framework and are tested once in `tests/analyses/`, not per plugin.

Example test structure for a plugin:

```python
import pytest

from polyzymd.analyses.contract import ObservableAggregate
from polyzymd.analyses.mda.frame_selection import FrameSelection
from polyzymd.analyses.my_analysis import MyAnalysis, My, MySettings


def test_compute_reports_one_value_per_selected_frame(synthetic_universe):
    """The observable covers the production window and states its unit."""
    frames = FrameSelection(start=0, stop=5, step=1, n_frames_total=5)
    settings = MySettings()

    observables = My().compute(synthetic_universe, frames, settings)

    assert [observable.unit for observable in observables] == ["A"]
    assert len(observables[0].values) == 5


def test_aggregates_over_three_replicates(synthetic_universe, run_contract_analysis):
    """Three identical replicates give the known mean and a zero SEM."""
    settings = MySettings()

    artifact = run_contract_analysis(MyAnalysis, settings, synthetic_universe)

    aggregate = ObservableAggregate.model_validate(artifact.payload["observables"][0])
    assert aggregate.n_replicates == 3
    assert aggregate.mean == pytest.approx(1.0)
    assert aggregate.sem == pytest.approx(0.0)
```

## Test Data

- Real simulation data lives at `../testing_analysis/` (outside repo)
- Do NOT commit trajectory files (.xtc, .dcd) or large PDB files to the repo
- For unit tests, use small synthetic data or mock objects
