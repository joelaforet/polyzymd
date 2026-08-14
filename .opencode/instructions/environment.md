# Environment Rules

## Pixi Environment (Critical)

The project uses **pixi** for environment management. All simulation and
analysis dependencies that are not available via pip are resolved from
conda-forge through pixi environments.

**Always prefix commands with:**
```bash
pixi run -e build <command>
```

**Or activate a shell:**
```bash
pixi shell -e build
```

**Never `pip install` these packages directly:**
- OpenMM (`openmm`)
- OpenFF Toolkit (`openff-toolkit`)
- OpenFF Interchange (`openff-interchange`)
- MDAnalysis (`MDAnalysis`)
- AmberTools (`ambertools`)
- parmed

These must come from conda-forge via the pixi environment.

## Prebuilt simulation artifacts

`polyzymd build` publishes `solvated_system.pdb`, `system.xml`, and
`build_manifest.json` as one integrity bundle. The manifest is written last and
records configuration identity, particle count, hashes, OpenMM version, and
canonical paths. Do not stage only part of this bundle between filesystems.
Legacy directories without a manifest are recovery inputs only and are accepted
only when all available particle counts agree.

## Environment Setup

The environment specification is managed by `pixi.toml` at the repo root.

To install:
```bash
pixi install -e build
```

To install polyzymd in editable mode within the env:
```bash
pixi run -e build pip install -e ".[dev,docs]"
```

## Dependency Management

- **Runtime deps** go in `pyproject.toml` under `[project.dependencies]`
- **Dev deps** go in `pyproject.toml` under `[project.optional-dependencies.dev]`
- **Docs deps** go in `pyproject.toml` under `[project.optional-dependencies.docs]`
- **Conda-only deps** are resolved by pixi from conda-forge
- Never add OpenMM/OpenFF/MDAnalysis to `pyproject.toml` — they are conda-only

## CI Workflows

GitHub Actions workflows in `.github/workflows/`:

| Workflow | Trigger | Description |
|----------|---------|-------------|
| `ci.yml` | Push/PR to main, dev | Lint + basic tests |
| `full-test.yml` | Manual / scheduled | Full test suite with pixi env |
| `release.yml` | Tag push | Build and publish to PyPI |

## Lazy Imports

Because not all environments have the full simulation stack, heavy dependencies
must be lazily imported inside functions:

```python
def analyze_trajectory(topology_path, trajectory_path):
    import MDAnalysis as mda  # Lazy import
    u = mda.Universe(str(topology_path), str(trajectory_path))
    ...
```

This allows the config/CLI modules to be imported without the full stack
installed (useful for documentation builds, CI lint jobs, etc.).
