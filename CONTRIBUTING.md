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
Add a radius of gyration analysis function

Measure the radius of gyration of an atom group at one frame with
MDAnalysis, run it per replicate through Study.timeseries, and add it to
polyzymd analyze with its figures.

Closes #42
```

Never force-push to `main` or `dev`.

## How to Contribute a New Analysis

This is the most common type of contribution. A new analysis is a function of
MDAnalysis atom groups; the study API supplies the replicate universes, the
production frames, the stored records and the cross-replicate statistics, so
adding one needs **one function** and **no changes to the study code**.

### Step-by-Step

1. **Read** `docs/source/explanation/analysis_api.md` and
   `docs/source/reference/analysis_functions.md`, which lists every shipped
   function and what it measures.

2. **Write the function** in `src/polyzymd/analyses/functions.py`:
   - a per-frame function takes atom groups at the current frame and returns
     one number, and runs through `Study.timeseries`;
   - a per-replicate function also takes `frames` and returns a number or one
     value per label (such as per residue), and runs through
     `Study.per_replicate`.
   Call MDAnalysis, MDTraj or another package for the measurement; do not
   reimplement it. The docstring says concretely what is measured and in which
   unit.

3. **Test it** in `tests/analyses/test_<name>.py`: a known answer on a small
   universe with placed atoms, then a run through the study on the synthetic
   OpenMM run directories of `tests/_support/analysis_testkit.py`.

4. **To ship it in `polyzymd analyze`**, add its settings to
   `FUNCTION_ANALYSES` and an `_analyze_<name>` in
   `src/polyzymd/analyses/protocols.py`, with its figures, and document it in
   `docs/source/reference/analysis_functions.md` and a how-to page.

5. **Run the test suite**: `pixi run -e build pytest tests/ -v`

The plugin framework (`Analysis`, `polyzymd new-analysis`, `polyzymd compare`)
still exists, but no shipped analysis uses it and it is being removed. Do not
add plugins.

### Key Rules

- **Use the study API** for trajectories: `Study.from_configs` resolves the
  topology and trajectory files, joins restart segments and applies the
  equilibration window. Do not build a `Universe` by hand in an analysis.

- **Import rules**: Import heavy third-party packages (MDAnalysis, matplotlib,
  mdtraj) lazily inside functions.

- **The replicate is the sampling unit** for every uncertainty and every
  comparison; frames never are.

- **Chain convention**: A=protein, B=substrate, C=polymer, D+=solvent.

### Checklist Before Opening a PR

- [ ] Function in `src/polyzymd/analyses/functions.py` with a docstring that
  says what it measures and in which unit
- [ ] Known-answer test and a study test in `tests/analyses/test_<name>.py`
- [ ] For a shipped analysis: `FUNCTION_ANALYSES` entry, `_analyze_<name>`,
  figures and docs
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
├── analyses/     # Study API and analysis functions — primary extension point
│   └── shared/   #   Reusable utilities (TrajectoryLoader, alignment, statistics)
├── exporters/    # GROMACS/other format exporters
├── data/         # Bundled data files (force fields, templates)
└── utils/        # Shared utilities
```

The `analyses/` directory is the primary extension point: `functions.py` holds
the shipped measurements, `study.py` and `timeseries.py` the study API, and
`protocols.py` the `polyzymd analyze` analyses. `base.py`, `orchestrator.py`,
`discovery.py`, `mda/` and `_framework/` are the plugin framework, which no
shipped analysis uses and which is being removed.

## Getting Help

- **Study API**: `docs/source/explanation/analysis_api.md`
- **Shipped functions**: `docs/source/reference/analysis_functions.md`
- **Issues**: https://github.com/joelaforet/polyzymd/issues

## License

By contributing, you agree that your contributions will be licensed under the
MIT License.
