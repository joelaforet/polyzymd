# Testing Rules

## Current State

The test suite covers the study API, the shipped analysis functions, `polyzymd analyze` and core infrastructure:

- **Test directory:** `tests/` with subdirectories mirroring the source tree
- **Test count:** run `pytest tests --collect-only -q | tail -1`; do not copy a number here
- **Fixtures:** `tests/conftest.py` with shared fixtures for simulation configs, mock data, etc.
- **Markers:** `@pytest.mark.slow` for tests requiring simulation data

### Directory Structure

```
tests/
├── conftest.py                  # Shared fixtures
├── _support/                    # Test utilities (analysis_testkit.py)
├── analyses/                    # Tests for analyses/ source tree
│   ├── test_study_timeseries.py # study.py, timeseries.py: Study, timeseries, per_replicate, stride
│   ├── test_reference.py        # reference.py
│   ├── test_figures.py          # figures.py
│   ├── test_protocols.py        # protocols.py: polyzymd analyze, ProtocolReport
│   ├── test_rmsf.py             # rmsf, rmsd_per_residue, rms_decomposition
│   ├── test_sasa.py             # sasa, residue_sasa
│   ├── test_secondary_structure.py  # dssp_occupancy
│   ├── test_residue_contacts.py     # residue_contacts, polyzymd analyze contacts
│   ├── test_residue_occlusion.py    # residue_occlusion
│   ├── test_contact_lifetimes.py    # contact_events, contact_lifetimes, --run mean_lifetime
│   ├── test_native_contacts.py      # native_contacts, polyzymd analyze native_contacts
│   ├── test_hydrogen_bonds_functions.py  # hbond_atoms, hydrogen_bonds, lifetimes, occupancies, hbond_count
│   ├── test_hydrogen_bonds_analyze.py    # polyzymd analyze hydrogen_bonds on OpenMM run directories
│   ├── test_segment_join.py     # loader repairs of restart-chain boundaries
│   ├── test_empty_segments.py   # loader skips of empty segments
│   ├── test_universe.py         # universe.py: UniverseProvider and input file records
│   ├── test_identity.py         # identity.py: the config hash must not change
│   ├── shared/                  # analyses/shared/ utilities (loader, window, statistics, ...)
│   ├── scientific/              # statistical and uncertainty contract tests
│   └── integration/             # Cross-analysis integration tests
├── cli/                         # Tests for cli/ source tree
│   ├── test_main.py
│   ├── test_main_recover.py
│   ├── test_main_status.py
│   ├── test_analyze.py          # polyzymd analyze options, -f refusal
│   ├── test_retired_commands.py # hidden compare and new-analysis stubs
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
│   ├── test_analysis_submit.py  # polyzymd analyze --submit scripts
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
pixi run -e build pytest tests/analyses/test_sasa.py -v

# Run tests matching a pattern
pixi run -e build pytest tests/ -v -k "rmsf"

# Run the hydrogen-bond tests
pixi run -e build pytest tests/ -v -k "hydrogen_bonds"
```

## Writing New Tests

When adding tests:

1. Place tests in the subdirectory matching the source module (e.g., `tests/analyses/` for study API functions)
2. Name files `test_<source_module>.py` with 1:1 correspondence to source files
3. Use pytest conventions (`test_` prefix for functions/methods)
4. Mock heavy dependencies (OpenMM, MDAnalysis) for unit tests
5. Use `@pytest.mark.slow` for tests requiring simulation data
6. Use existing fixtures from `tests/conftest.py`

### Testing an analysis function

A new analysis is a function in `analyses/functions.py`; test it in `tests/analyses/test_<name>.py`:

1. **Known answer on placed geometry**: build a small `MDAnalysis.Universe`
   with atoms at chosen positions and check the function's value by hand
   arithmetic, as `test_hydrogen_bonds_functions.py` does.
2. **Through the study**: write OpenMM run directories with
   `tests/_support/analysis_testkit.py` (`write_simulation_config`,
   `write_openmm_replicate`, `write_openmm_frames`), run
   `Study.timeseries` or `Study.per_replicate`, and check the replicate values
   and the record under `polyzymd_results/`.
3. **polyzymd analyze**: for a shipped analysis, run `analyze(name, configs)` or
   the CLI and check `metric`, `unit`, `all_runs`, the replicate values and the
   figures written.


## Test Data

- Real simulation data lives at `../testing_analysis/` (outside repo)
- Do NOT commit trajectory files (.xtc, .dcd) or large PDB files to the repo
- For unit tests, use small synthetic data or mock objects
- For integration tests, use the `@pytest.mark.slow` marker and document
  the required data paths
