# Package Structure and Release Guide

Use this page when you need to understand how PolyzyMD is packaged, why the
project is pixi-first, and what release automation expects from contributors.

## Packaging model

PolyzyMD uses standard Python packaging for metadata and distributions, but the
full molecular-simulation stack is managed by pixi. A wheel can be built with
`pyproject.toml`, but pip alone is not the recommended full install path because
OpenMM, OpenFF, MDAnalysis, AmberTools, RDKit, PACKMOL, and related scientific
packages are resolved most reliably from conda-forge.

For contributor and user workflows, prefer:

```bash
pixi install -e build
pixi run -e build polyzymd --help
```

Use pip-only installs only for package metadata checks, lightweight import
checks, or release validation jobs that intentionally avoid the simulation
stack.

## Repository layout

PolyzyMD uses the `src/` layout so tests import the installed package instead of
accidentally importing files from the repository root.

```text
polyzymd/
├── .github/workflows/        # CI, full-test, and release workflows
├── docs/source/              # Sphinx + MyST documentation
│   ├── get_started/          # Installation and quickstart pages
│   ├── tutorials/            # Learning-oriented walkthroughs
│   ├── how_to/               # Task-oriented guides
│   ├── reference/            # CLI, schema, and data-layout reference
│   ├── explanation/          # Architecture and design rationale
│   ├── contributor_guide/    # Contributor documentation
│   └── api/                  # Autodoc API pages
├── src/polyzymd/
│   ├── analyses/             # Analysis plugin system and built-in plugins
│   ├── builders/             # Molecular system construction
│   ├── cli/                  # Click command-line interface
│   ├── config/               # Pydantic configuration schema and loading
│   ├── core/                 # Shared core types
│   ├── data/                 # Bundled data and solvent/co-solvent libraries
│   ├── engines/              # Engine-specific execution/export helpers
│   ├── exporters/            # Format exporters
│   ├── simulation/           # Simulation runners and continuation support
│   ├── templates/            # Project and config templates
│   ├── utils/                # Shared utilities
│   └── workflow/             # SLURM and orchestration helpers
├── tests/
│   ├── analyses/             # Plugin and analysis-framework tests
│   ├── cli/                  # CLI tests
│   ├── config/               # Configuration tests
│   ├── engines/              # Engine integration/export tests
│   ├── exporters/            # Exporter tests
│   ├── regression/           # Regression coverage
│   ├── simulation/           # Runner and continuation tests
│   ├── utils/                # Utility tests
│   └── workflow/             # SLURM/workflow tests
├── pixi.toml                 # Conda-forge environments and tasks
├── pyproject.toml            # Python package metadata and tool config
└── README.md
```

## `pyproject.toml`

`pyproject.toml` defines the Python package metadata, build backend, CLI entry
point, and development-tool configuration.

```toml
[build-system]
requires = ["hatchling"]
build-backend = "hatchling.build"

[project]
name = "polyzymd"
version = "1.3.0"
requires-python = ">=3.12"

[project.scripts]
polyzymd = "polyzymd.cli.main:main"
```

The project currently keeps the release version in three places:

- `pyproject.toml`
- `pixi.toml`
- `src/polyzymd/__init__.py`

When preparing a release, update all three to the same version, for example
`1.3.1`.

## Pixi environments

Current pixi environments are:

| Environment | Use |
|-------------|-----|
| `build` | Local contributor environment with docs, tests, packaging tools, and flexible non-CUDA OpenMM. Use this for editing, docs, configuration, system-building checks, and most tests. |
| `cuda-12-4` | Cluster/GPU workflow for systems whose driver supports CUDA 12.4. |
| `cuda-12-6` | Cluster/GPU workflow for systems whose driver supports CUDA 12.6. |

Install the local contributor environment with:

```bash
pixi install -e build
```

For GPU cluster work, choose the CUDA environment that does not exceed the
cluster driver version:

```bash
pixi install -e cuda-12-6
pixi shell -e cuda-12-6
```

Heavy scientific dependencies must remain pixi-managed. Do not add OpenMM,
OpenFF, MDAnalysis, AmberTools, PACKMOL, or RDKit installation instructions that
ask contributors to install them into a system Python.

## Lazy imports for heavy dependencies

PolyzyMD keeps package imports lightweight by importing heavy dependencies only
inside functions or methods that need them. This lets CI, docs, and metadata
checks import modules that do not actually run simulations.

```python
def load_trajectory(topology_path, trajectory_path):
    import MDAnalysis as mda

    return mda.Universe(str(topology_path), str(trajectory_path))
```

Avoid module-level imports of OpenMM, OpenFF, MDAnalysis, ParmEd, PDBFixer, or
other simulation-stack packages in new code.

## CI workflows

### `ci.yml`

The main CI workflow runs on pushes and pull requests to `main` and `dev`.

- **Lint:** installs `black` and `ruff`, then runs `ruff check`, `black --check`,
  and `ruff format --check` against the source tree.
- **Type check:** installs lightweight Python dependencies and runs mypy against
  `src/polyzymd/config/`. This job is informational and currently uses
  `continue-on-error` because many simulation-stack packages do not provide full
  type stubs.
- **Build package:** builds the sdist/wheel with `python -m build` and validates
  distributions with `twine check`.
- **Docs build:** uses the pixi `build` environment and treats Sphinx warnings as
  errors with `SPHINXOPTS="-W --keep-going"`.
- **Test installation:** installs the built wheel, checks basic imports, and
  verifies `polyzymd --help` and `polyzymd --version`.

### `full-test.yml`

The full test workflow runs on pushes and pull requests to `main` and `dev`, on a
weekly schedule, and by manual dispatch. It uses the pixi `build` environment on
Ubuntu and macOS with Python 3.12, then runs:

```bash
pixi run -e build pytest -v --cov=polyzymd --cov-report=xml --color=yes tests/
```

Coverage upload is attempted with Codecov but does not fail the workflow if the
upload fails.

### `release.yml`

The release workflow runs when a tag matching `v*.*.*` is pushed. It:

1. builds and `twine check`s the distributions,
2. publishes non-prerelease tags to PyPI through trusted publishing, and
3. creates a GitHub release with generated notes and attached distributions.

## Release checklist

For a normal version bump:

1. Update the version in `pyproject.toml`, `pixi.toml`, and
   `src/polyzymd/__init__.py`.
2. Run focused checks for the changed files, then at minimum run the docs build
   before handing off for release review:

   ```bash
   pixi run -e build make -C docs clean html
   ```

3. Commit the version bump with a focused message, such as
   `chore(release): bump version to 1.3.1`.
4. After review, tag the release, for example `v1.3.1`. The tag triggers
   `release.yml`.

## Best practices

- Keep `__init__.py` lightweight and preserve lazy imports.
- Keep the full simulation stack in pixi environments.
- Use `src/` layout conventions for new package modules.
- Update docs, tests, CLI help, and configuration reference pages when behavior
  changes.
- Do not commit generated caches, trajectories, build artifacts, or local
  environment files.

## See also

- {doc}`setup` — contributor environment setup
- {doc}`contributing` — contribution workflow and review expectations
- {doc}`../get_started/installation` — user installation guide
- {doc}`../explanation/architecture` — architecture overview
