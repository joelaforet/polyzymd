# Package Structure and Development Guide

This guide explains how PolyzyMD is structured as a Python package, the design decisions behind it, and how to contribute or create similar packages.

## Overview

PolyzyMD is designed to be:

- **pip-installable**: `pip install polyzymd`
- **pixi-first for heavy scientific dependencies**: use `pixi install -e <env>`
  for the reproducible OpenMM/OpenFF stack
- **Testable in CI**: Works with GitHub Actions for automated testing
- **Developer-friendly**: Clear structure for contributors

This document explains the key concepts that make this work.

---

## Package Layout

### The `src/` Layout

PolyzyMD uses a **src layout**, which is considered best practice for pip-installable packages:

```
polyzymd/                    # Repository root
├── src/
│   └── polyzymd/            # The actual Python package
│       ├── __init__.py
│       ├── config/
│       ├── builders/
│       ├── simulation/
│       └── ...
├── tests/                   # Test suite (outside the package)
├── docs/                    # Documentation
├── devtools/                # Development tools and helper scripts
├── pyproject.toml           # Package metadata and build config
├── README.md
└── LICENSE
```

### Why `src/` Layout?

The `src/` layout has several advantages over a "flat" layout:

| Aspect | src/ Layout | Flat Layout |
|--------|-------------|-------------|
| **Import safety** | Forces testing against installed package | Can accidentally import local directory |
| **Editable installs** | Works correctly with `pip install -e .` | May have import conflicts |
| **Clear separation** | Package code is isolated | Package mixed with repo files |

**Example of the problem with flat layout:**

```python
# With flat layout, this might import local files instead of installed package:
import polyzymd  # Which polyzymd? Local or installed?
```

With `src/` layout, you must install the package (even in editable mode) before imports work, ensuring you always test the real package.

### Alternative: Flat Layout

Some packages (like Polymerist) use a flat layout:

```
polymerist/                  # Repository root
├── polymerist/              # Package has same name as repo
│   ├── __init__.py
│   └── ...
├── pyproject.toml
└── ...
```

This is simpler but requires more care to avoid import issues.

---

## The `pyproject.toml` File

The `pyproject.toml` file is the modern standard for Python package configuration. It replaces the older `setup.py` approach.

### Key Sections

#### Build System

```toml
[build-system]
requires = ["hatchling"]
build-backend = "hatchling.build"
```

This tells pip which tool to use for building the package. Common options:
- **hatchling**: Modern, fast, good defaults (what PolyzyMD uses)
- **setuptools**: Traditional, widely supported
- **flit**: Simple, minimal configuration

#### Project Metadata

```toml
[project]
name = "polyzymd"
version = "1.0.0"
description = "Molecular dynamics simulation toolkit for enzyme-polymer systems"
readme = "README.md"
license = "MIT"
authors = [
    { name = "Joseph R. Laforet Jr.", email = "jola3134@colorado.edu" }
]
requires-python = ">=3.10"
```

#### Dependencies

```toml
dependencies = [
    # Core dependencies that pip can install
    "pydantic>=2.0.0",
    "pyyaml>=6.0",
    "click>=8.0.0",
    "numpy>=1.21.0,<2.0.0",
]

[project.optional-dependencies]
# Groups of optional dependencies
dev = ["pytest>=7.0.0", "ruff>=0.1.0"]
docs = ["sphinx>=6.0.0", "myst-parser>=1.0.0"]
```

#### Entry Points (CLI Commands)

```toml
[project.scripts]
polyzymd = "polyzymd.cli.main:main"
```

This creates the `polyzymd` command that users can run from terminal.

---

## Handling Heavy Dependencies

### The Challenge

PolyzyMD depends on packages that are difficult or impossible to install via pip:

- **OpenMM**: GPU-accelerated MD engine (resolved from conda-forge via pixi)
- **OpenFF Toolkit**: Force field tools (resolved from conda-forge via pixi)
- **RDKit**: Chemistry toolkit (managed in the pixi environment)
- **PACKMOL**: Molecular packing (resolved from conda-forge via pixi)
- **AmberTools**: Optional, for AM1-BCC charging backend (resolved from conda-forge via pixi)

If we list these in `dependencies`, `pip install polyzymd` will fail.

### The Solution: Lazy Imports

We use **lazy imports** so the package can be imported without these heavy dependencies:

```python
# polyzymd/__init__.py

__version__ = "1.0.0"

def __getattr__(name):
    """Lazy import heavy modules only when accessed."""
    if name == "SystemBuilder":
        from polyzymd.builders.system_builder import SystemBuilder
        return SystemBuilder
    raise AttributeError(f"module 'polyzymd' has no attribute {name!r}")
```

This means:

```python
import polyzymd                    # Works without OpenFF installed
print(polyzymd.__version__)        # Works - just returns "1.0.0"

from polyzymd import SystemBuilder  # NOW imports OpenFF (fails if not installed)
```

### Why This Matters

1. **CI can run basic tests** without a full simulation environment
2. **Users can install** the package and see helpful error messages
3. **Documentation builds** don't require simulation dependencies

### Recommended Installation

We recommend users install via pixi so the heavy scientific dependencies are
resolved reproducibly from `pixi.toml`:

```bash
git clone https://github.com/joelaforet/polyzymd.git
cd polyzymd
pixi install -e build
pixi shell -e build
```

### Optional Charging Backends

PolyzyMD defaults to **NAGL** for partial charge assignment, which is fast and
doesn't require additional dependencies beyond the OpenFF stack.

For AM1-BCC charges, you can optionally install additional backends:

#### AmberTools (sqm)

AmberTools is already included in the PolyzyMD pixi environments, so no extra
installation step is usually required.

Then use in your code:

```python
from polyzymd.utils.charging import get_charger

# Use AmberTools for AM1-BCC
charger = get_charger("am1bcc", toolkit="ambertools")
charged_mol = charger.charge_molecule(molecule)
```

#### OpenEye Toolkit (Commercial)

If you have an OpenEye license:

```python
charger = get_charger("am1bcc", toolkit="openeye")
```

Note: NAGL is recommended for most use cases as it's faster and produces
comparable results to AM1-BCC.

---

## Continuous Integration (CI)

### Two-Tier CI Strategy

PolyzyMD uses a two-tier CI approach:

#### Tier 1: Basic CI (Runs on Every PR)

Fast checks that don't require conda:

- **Linting**: Code style with `ruff`
- **Type checking**: Static analysis with `mypy`
- **Build verification**: Ensure package builds correctly
- **Import tests**: Test that config module imports

```yaml
# .github/workflows/ci.yml
jobs:
  lint:
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v4
      - uses: actions/setup-python@v5
      - run: pip install ruff
      - run: ruff check src/polyzymd/
      
  build:
    runs-on: ubuntu-latest
    steps:
      - run: pip install build
      - run: python -m build
```

#### Tier 2: Full CI (Runs Weekly or on Release)

Comprehensive tests with the full pixi-managed simulation environment:

```yaml
# .github/workflows/full-test.yml
jobs:
  test:
    runs-on: ubuntu-latest
    steps:
      - uses: prefix-dev/setup-pixi@v0
      - run: pixi install -e build
      - run: pixi run -e build pip install . --no-deps
      - run: pixi run -e build pytest tests/
```

### The Pixi Manifest

Environment definitions live in `pixi.toml`:

```toml
[workspace]
name = "polyzymd"
channels = ["conda-forge"]

[environments]
build = { features = ["base", "build", "test", "docs"] }
cuda-12-4 = { features = ["base", "cuda-12-4", "test", "docs"] }
cuda-12-6 = { features = ["base", "cuda-12-6", "test", "docs"] }
```

---

## Version Management

### Simple Approach (Current)

PolyzyMD uses a simple hardcoded version:

```python
# polyzymd/__init__.py
__version__ = "1.0.0"
```

```toml
# pyproject.toml
[project]
version = "1.0.0"
```

When releasing a new version, update both files.

### Dynamic Versioning (Alternative)

Some packages use tools like `versioningit` or `setuptools-scm` to derive version from git tags:

```toml
[project]
dynamic = ["version"]

[tool.versioningit]
default-version = "1+unknown"

[tool.versioningit.write]
file = "polymerist/_version.py"
```

```python
from importlib.metadata import version
__version__ = version(__name__)
```

**Pros**: Version always matches git tags, no manual updates
**Cons**: More complex, requires understanding of the tooling

---

## Creating a Release

### Steps for a New Release

1. **Update version numbers**
   ```bash
   # Edit pyproject.toml and src/polyzymd/__init__.py
   # Change version = "1.0.0" to version = "1.1.0"
   ```

2. **Commit the version bump**
   ```bash
   git add pyproject.toml src/polyzymd/__init__.py
   git commit -m "Bump version to 1.1.0"
   git push origin main
   ```

3. **Create and push a git tag**
   ```bash
   git tag -a v1.1.0 -m "Release 1.1.0"
   git push origin v1.1.0
   ```

4. **GitHub Actions handles the rest**
   - Builds the package
   - Publishes to PyPI
   - Creates GitHub release

### PyPI Publishing

The release workflow uses "trusted publishing" - no API tokens needed:

```yaml
# .github/workflows/release.yml
- name: Publish to PyPI
  uses: pypa/gh-action-pypi-publish@release/v1
  # Uses OIDC authentication - configure at pypi.org
```

To set up trusted publishing:
1. Go to https://pypi.org/manage/project/polyzymd/settings/publishing/
2. Add a new "pending publisher" with:
   - Owner: `joelaforet`
   - Repository: `polyzymd`
   - Workflow: `release.yml`

---

## Directory Structure Reference

```
polyzymd/
├── .github/
│   └── workflows/
│       ├── ci.yml           # Basic CI (lint, build, import tests)
│       ├── release.yml      # Publish to PyPI on tag
│       └── full-test.yml    # Full pixi test suite
│
├── devtools/                # Optional helper scripts / dev tooling
│
├── docs/
│   ├── source/
│   │   ├── conf.py          # Sphinx configuration
│   │   ├── index.rst        # Documentation index
│   │   └── tutorials/       # User guides
│   └── requirements.txt     # Docs build dependencies
│
├── src/
│   └── polyzymd/
│       ├── __init__.py      # Package init (lazy imports)
│       ├── config/          # Configuration (pydantic schemas)
│       ├── builders/        # System building
│       ├── simulation/      # MD runners
│       ├── workflow/        # SLURM/HPC integration
│       ├── exporters/       # GROMACS export
│       ├── cli/             # Command-line interface
│       └── data/            # Bundled data files
│
├── tests/
│   ├── __init__.py
│   ├── test_imports.py      # Basic import tests
│   └── test_config.py       # Configuration tests
│
├── .gitignore
├── .readthedocs.yaml        # ReadTheDocs configuration
├── LICENSE
├── README.md
└── pyproject.toml           # Package metadata and build config
```

---

## Best Practices Summary

### Do

- Use `src/` layout for clear package boundaries
- Use lazy imports for heavy dependencies
- Keep `__init__.py` lightweight
- Separate fast CI (lint) from slow CI (full tests)
- Use `pyproject.toml` for all configuration
- Document installation requirements clearly

### Don't

- Don't assume pip alone can provision the full OpenMM/OpenFF stack
- Don't eagerly import heavy modules at package init
- Don't mix test artifacts with source code
- Don't hardcode paths - use `importlib.resources` for data files
- Don't commit generated files (`.pyc`, `__pycache__`, etc.)

---

## Further Reading

- [Python Packaging User Guide](https://packaging.python.org/)
- [PyPA Sample Project](https://github.com/pypa/sampleproject)
- [Hatchling Documentation](https://hatch.pypa.io/)
- [Scientific Python Library Development Guide](https://learn.scientific-python.org/development/)
- [Polymerist Repository](https://github.com/timbernat/polymerist) - Reference implementation

---

## See Also

- {doc}`contributing` - How to contribute to PolyzyMD
- {doc}`architecture` - Internal architecture overview
- {doc}`installation` - User installation guide
