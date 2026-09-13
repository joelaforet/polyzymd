# Testing Rules

## Current State

The test suite covers the full analysis plugin system and core infrastructure:

- **Test directory:** `tests/` with subdirectories mirroring the source tree
- **Test count:** 2,284 collected. Run `pytest tests -q --collect-only | tail -1`
  rather than trusting this number; update it when it drifts
- **Fixtures:** `tests/conftest.py` with shared fixtures for comparison configs, mock data, etc.
- **Markers:** `@pytest.mark.slow` for tests requiring simulation data

### Directory Structure

```
tests/
├── conftest.py                  # Shared fixtures
├── _support/                    # Test utilities
├── analyses/                    # Tests for analyses/ source tree
│   ├── conftest.py              # synthetic_universe, run_contract_analysis,
│   │                            # and the plot uncertainty audit
│   ├── test_contract.py         # the observable contract and the cache identity
│   ├── test_contract_plots.py   # one figure per kind
│   ├── test_base.py             # the Analysis lifecycle
│   ├── test_analysis_lifecycle.py
│   ├── test_discovery.py
│   ├── test_orchestrator.py, test_orchestrator_workers.py,
│   │   test_orchestrator_cost_hints.py
│   ├── test_protocols.py, test_protocols_contract_artifacts.py
│   ├── test_aggregate_validation.py, test_cache_identity.py
│   ├── test_stats.py            # interpret_direction and format_pct
│   ├── mda/                     # artifacts, store, frame selection, universe,
│   │                            # the replicate lifecycle, the import surface
│   ├── shared/                  # alignment, autocorrelation, loader, paths,
│   │                            # plotting, selections, statistics, window
│   ├── plugins/                 # one file per analysis plugin
│   ├── scientific/              # known-answer and best-practice checks
│   ├── parity/                  # frozen numbers on the real campaign data
│   └── integration/             # cross-plugin checks
├── cli/                         # Tests for cli/ source tree
├── config/                      # Tests for config/ source tree
├── engines/, builders/, core/, simulation/, workflow/, utils/
└── data/                        # small fixtures, no trajectories
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
