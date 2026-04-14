# Contributing to PolyzyMD

Thank you for your interest in contributing to PolyzyMD! This guide covers
everything you need to get started.

## Setting Up Your Development Environment

PolyzyMD uses [pixi](https://pixi.sh) for environment management. Pixi resolves
OpenMM, OpenFF, MDAnalysis, and all other heavy dependencies from conda-forge
with reproducible lockfiles.

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
PATH. **Do not `pip install` OpenMM, MDAnalysis, or other simulation-stack
packages outside the pixi environment.**

## Code Quality Checks

Run these before every commit:

```bash
# Lint
ruff check src/polyzymd/

# Format check (add --diff to see what would change)
ruff format --check src/polyzymd/

# Auto-format
ruff format src/polyzymd/

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

Implement compute_replicate using TrajectoryLoader for trajectory
handling. Aggregate with SEM across replicates. Wire into default
scalar comparison path via extract_metrics.

Closes #42
```

Never force-push to `main` or `dev`.

## How to Contribute a New Analysis Plugin

This is the most common type of contribution. PolyzyMD's plugin system is
designed so that adding a new analysis requires **one package** and **no changes
to core code**.

### Step-by-Step

1. **Read the tutorial**: `docs/source/tutorials/extending_analyses.md` — it
   walks through every component of a plugin and has a complete working example.

2. **Study an existing plugin**: Start with `src/polyzymd/analyses/secondary_structure/`
   (simplest real plugin) or `src/polyzymd/analyses/rmsf/` (simple with default
   comparison path).

3. **Create your plugin package**: Use `polyzymd new-analysis <name>` to
   scaffold automatically, or create `src/polyzymd/analyses/<name>/` with an
   `__init__.py` manually

4. **Subclass `Analysis`** and implement the required pieces:
   - `name` — unique lowercase string identifier
   - `Settings` — Pydantic v2 `BaseModel` with sensible defaults
   - `compute_replicate(ctx, replicate)` — per-replicate computation
   - `aggregate(ctx, results)` — summarize across replicates

5. **Choose your comparison path**:
   - **Default (recommended)**: Implement `extract_metrics()` returning
     `MetricValue` objects. The framework handles t-tests, Cohen's d, ANOVA,
     and ranking automatically.
   - **Custom**: Override `compare()` entirely for multi-metric or entry-table
     analyses.

6. **Write tests**: Create `tests/test_<name>_plugin.py` following the pattern
   of existing plugin tests. The standard test structure covers:
   - Discovery and class attributes
   - Settings validation
   - `compute_replicate` with mocked trajectories
   - `aggregate` with sample data
   - `extract_metrics` (if applicable)
   - Plot generation
   - Full lifecycle integration

7. **Run the test suite**: `pixi run -e build pytest tests/ -v`

### Key Rules

- **Use `TrajectoryLoader`** from `analyses/shared/` for trajectory loading —
  it handles topology discovery, segment daisy-chaining, timestep detection,
  and equilibration skipping.

- **Import rules**: Import framework utilities (TrajectoryLoader, etc.) at
  module level. Import heavy third-party packages (MDAnalysis, matplotlib,
  mdtraj) lazily inside methods. This matters for testability — `@patch`
  targets must be importable at the module level.

- **Use the context objects**: `ReplicateContext` and `AggregateContext` provide
  everything you need (sim config, settings, output paths). Never construct
  configs manually.

- **Result serialization**: If your aggregated result is a Pydantic model that
  inherits from `BaseAnalysisResult`, set `AggregatedResultClass = YourModel`
  on the plugin class. The framework handles serialization automatically. For
  dict results, the framework falls back to `json.loads()`.

- **Chain convention**: A=protein, B=substrate, C=polymer, D+=solvent.

### Checklist Before Opening a PR

- [ ] Plugin package in `src/polyzymd/analyses/<name>/`
- [ ] `name` class variable set (lowercase, unique)
- [ ] `Settings` inner class with default values for all fields
- [ ] `compute_replicate` and `aggregate` implemented
- [ ] `extract_metrics` implemented (or `compare` overridden for custom path)
- [ ] `AggregatedResultClass` set if using a Pydantic result model
- [ ] Test file in `tests/test_<name>_plugin.py`
- [ ] `ruff check src/polyzymd/` passes
- [ ] `ruff format --check src/polyzymd/` passes
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
- **ruff** for linting and formatting
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
├── analyses/     # Plugin system — primary extension point
│   ├── shared/   #   Reusable utilities (TrajectoryLoader, alignment, statistics)
│   └── <name>/   #   One sub-package per analysis type (all plugins are packages)
├── compare/      # Cross-condition statistics, formatters, config, IO
├── exporters/    # GROMACS/other format exporters
├── data/         # Bundled data files (force fields, templates)
└── utils/        # Shared utilities
```

The `analyses/` directory is the primary extension point. Each sub-package
(`<name>/`) is one analysis plugin. Private `_*.py` modules inside each package
are internal implementation details (calculators, result models, formatters).

## Getting Help

- **Tutorial**: `docs/source/tutorials/extending_analyses.md`
- **Plugin contract**: `src/polyzymd/analyses/base.py` (class docstring)
- **Issues**: https://github.com/joelaforet/polyzymd/issues

## License

By contributing, you agree that your contributions will be licensed under the
MIT License.
