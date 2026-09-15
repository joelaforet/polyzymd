# Contributing to PolyzyMD

Thank you for your interest in contributing to PolyzyMD! This guide covers
everything you need to get started.

## Setting Up Your Development Environment

PolyzyMD uses [pixi](https://pixi.sh) for environment management. Pixi resolves
the full scientific and simulation stack from conda-forge with reproducible
lockfiles.

```bash
# 1. Install pixi (if you don't have it)
curl -fsSL https://pixi.sh/install.sh | sh

# 2. Clone and set up
git clone https://github.com/joelaforet/polyzymd.git
cd polyzymd
pixi install -e build
pixi shell -e build
```

After `pixi shell`, the `polyzymd` CLI and all development tools are on your
PATH. Pixi/conda-forge is the recommended and supported route for full
scientific, simulation, and CUDA workflows. A best-effort pip path is available
for analysis-only workflows with `pip install "polyzymd[analysis]"`; that extra
may install MDAnalysis and MDTraj. Do not use system-pip installs for OpenMM,
OpenFF, AmberTools, RDKit, PACKMOL, CUDA, or full simulation-stack workflows.

## Code Quality Checks

Run these before every commit:

```bash
# Lint
ruff check src/polyzymd/

# Format check
black src/ --check

# Auto-format
black src/

# Type check (config module)
pixi run -e build mypy src/polyzymd/config/

# Tests
pixi run -e build pytest tests/ -v
```

CI runs lint, format check, and type check on every push. All checks must pass
before a PR can merge.

## Git Workflow

- **`main`** — stable releases
- **`dev`** — integration branch
- **Feature branches** — `feature/<short-description>`, branched from `dev`

Commit messages use imperative mood with a 50-character subject line:

```
Add radius of gyration analysis plugin

Report the radius of gyration of the protein selection as one
mean_of_timeseries observable in angstrom. The framework owns the aggregation,
the testing and the figure.

Closes #42
```

Never force-push to `main` or `dev`.

## How to Contribute a New Analysis Plugin

This is the most common type of contribution. PolyzyMD's plugin system is
designed so that adding a new analysis requires one module and no changes to
core code. A plugin is about forty lines.

### Step-by-Step

1. **Read the guide**: `docs/source/contributor_guide/analysis_plugins/index.md`
   walks through the whole contract and has a complete working example.

2. **Scaffold it**: `polyzymd new-analysis <name>` writes
   `src/polyzymd/analyses/<name>.py` and
   `tests/analyses/plugins/test_<name>.py`. Edit those two files and nothing
   else.

3. **Study an existing plugin**: `src/polyzymd/analyses/rg.py` is a short one,
   and `src/polyzymd/analyses/contacts/` is one that also returns an extra
   sidecar array.

4. **Declare the plugin** as a plain class with:
   - `name`, a unique lowercase string identifier
   - `Settings`, a pydantic v2 `BaseModel` with sensible defaults
   - `references`, the citations for the method
   - `compute(universe, frames, settings)`, returning a sequence of
     `Observable`

5. **Pick a kind for each observable**: `mean_of_timeseries`, `fluctuation`,
   `fraction` or `profile`. The kind decides how the replicate reduces, how
   conditions are compared, and which figure is drawn. Never aggregate across
   replicates, run a test, write a file or import matplotlib inside
   `compute()`. End the module with
   `NameAnalysis = contract_analysis(Name)`, which is what discovery finds.

6. **Write two tests**: one asserting what `compute()` measures on the
   `synthetic_universe` fixture, and one running `run_contract_analysis` over
   three replicates and asserting on the `ObservableAggregate` fields.

7. **Run the test suite**: `pixi run -e build pytest tests/ -v`

### Key Rules

- **Use `TrajectoryLoader`** from `analyses/shared/` for trajectory loading —
  it handles topology and trajectory discovery, segment daisy-chaining, and
  timestep access. Equilibration-aware frame slicing belongs in the shared
  window helpers.

- **Import rules**: Import framework utilities (TrajectoryLoader, etc.) at
  module level. Import heavy third-party packages (MDAnalysis, matplotlib,
  mdtraj) lazily inside methods. This matters for testability — `@patch`
  targets must be importable at the module level.

- **Take what you are handed**: `compute()` receives the loaded universe, the
  resolved production window and your settings model. Never load a config or
  resolve a file path yourself.

- **Result serialization**: there is nothing to do. The framework writes the
  replicate artifact, the per-frame series in an NPZ sidecar, the condition
  aggregate and the comparison. A plugin that returns a bulky extra array
  returns it as a second element of a `(observables, extras)` pair and the
  framework writes that sidecar too.

- **Chain convention**: A=protein, B=substrate, C=polymer, D+=solvent.

### Checklist Before Opening a PR

- [ ] Plugin module at `src/polyzymd/analyses/<name>.py`
- [ ] `name` set (lowercase, unique), `Settings` a pydantic model with defaults,
  `references` naming the method paper
- [ ] `compute()` returns observables and does nothing else
- [ ] Every observable states a real `unit`, not the scaffold's `"TODO"`
- [ ] The module ends with `NameAnalysis = contract_analysis(Name)`
- [ ] Test file in `tests/analyses/plugins/test_<name>.py`
- [ ] `ruff check src/polyzymd/` passes
- [ ] `black src/ --check` passes
- [ ] `pixi run -e build pytest tests/ -v` passes

## Other Types of Contributions

### Bug Fixes

1. Open an issue describing the bug with reproduction steps
2. Branch from `dev`: `git checkout -b fix/<short-description> dev`
3. Write a failing test, then fix the bug
4. Run the full test suite
5. Open a PR referencing the issue

### Documentation

Docs are built with Sphinx and MyST-Parser. Source lives in `docs/source/`.

```bash
# Build docs locally
pixi run -e build make -C docs clean html

# View in browser
open docs/build/html/index.html
```

Tutorials follow the [Diataxis](https://diataxis.fr/) framework — check whether
your content is a tutorial (learning-oriented), how-to guide (task-oriented),
reference (information-oriented), or explanation (understanding-oriented), and
place it accordingly.

### Configuration or Build System

PolyzyMD uses:
- **hatchling** for the Python build backend
- **pixi** for environment management (`pixi.toml`)
- **ruff** for linting
- **black** for formatting
- **pytest** for testing

## Project Layout

```
src/polyzymd/
├── cli/          # Click CLI entry point
├── config/       # Pydantic v2 config models, YAML loading
├── builders/     # System construction (PDB to parameterized topology)
├── simulation/   # OpenMM simulation runners
├── workflow/     # Orchestration (build, simulate, analyze)
├── core/         # Base classes, shared types
├── analyses/     # Plugin system + comparison framework — primary extension point
│   ├── shared/   #   Reusable utilities (TrajectoryLoader, alignment, statistics)
│   └── <name>/   #   One sub-package per analysis type (all plugins are packages)
├── exporters/    # GROMACS/other format exporters
├── data/         # Bundled data files (force fields, templates)
└── utils/        # Shared utilities
```

The `analyses/` directory is the primary extension point. Each sub-package
(`<name>.py`) is one analysis plugin. There is no plugin package holding
private helper modules: persistence, statistics and figures belong to the
framework.

## Getting Help

- **Guide**: `docs/source/contributor_guide/analysis_plugins/index.md`
- **Plugin contract**: `src/polyzymd/analyses/contract.py` (module docstring)
- **Issues**: https://github.com/joelaforet/polyzymd/issues

## License

By contributing, you agree that your contributions will be licensed under the
MIT License.
