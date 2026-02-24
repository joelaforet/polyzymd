# PolyzyMD — Agent Instructions

> Computational toolkit for enzyme-polymer conjugate MD simulations.
> Python >=3.10 | MIT License | hatchling build | src layout

## Environment

**All simulation-stack commands MUST use the conda environment:**

```bash
mamba run -n polyzymd-env <command>
```

The `polyzymd-env` environment contains OpenMM, OpenFF, MDAnalysis and other
heavy dependencies that are conda-only. Never `pip install` these.

**Quick commands:**

| Task | Command |
|------|---------|
| Install (editable) | `mamba run -n polyzymd-env pip install -e ".[dev,docs]"` |
| Run tests | `mamba run -n polyzymd-env pytest tests/ -v` |
| Lint | `ruff check src/` |
| Format | `black src/ --check` (or `black src/` to fix) |
| Build docs | `mamba run -n polyzymd-env make -C docs clean html` |
| Type check | `mamba run -n polyzymd-env mypy src/polyzymd` |

## Git Workflow

- **Branches:** `main` (stable), `dev` (integration), `feature/*` (work)
- Currently on `feature/analysis-module` (25 commits ahead of main)
- Commit messages: imperative mood, 50-char subject, reference issues (`#20`)
- Run `ruff check` and `black --check` before committing
- Never force-push to `main` or `dev`

## Architecture Quick Reference

```
src/polyzymd/
├── cli/          # Click CLI (main.py = entry point)
├── config/       # Pydantic v2 models (schema.py), YAML loading
├── builders/     # System construction (PDB → parameterized topology)
├── simulation/   # OpenMM simulation runners
├── workflow/     # Orchestration (build → simulate → analyze)
├── core/         # Base classes, shared types
├── analysis/     # Post-simulation analysis (contacts, RMSD, etc.)
├── compare/      # Multi-condition comparison engine
├── exporters/    # GROMACS/other format exporters
├── data/         # Bundled data files (force fields, templates)
├── utils/        # Shared utilities
└── configs/      # Default YAML configs
```

## Key Patterns

- **Chain convention:** A=protein, B=substrate, C=polymer, D+=solvent
- **Factory pattern:** `ClassName.from_config(config)` or `ClassName.from_yaml(path)`
- **Lazy imports:** Heavy deps (OpenMM, MDAnalysis) imported inside functions/methods
- **ABC + Strategy:** `ContactCriteria`, `MolecularSelector`, `MoleculeCharger`
- **Registry pattern:** `ComparatorRegistry`, `PlotterRegistry` for extensibility
- **Config:** Pydantic v2 `BaseModel` subclasses with `model_validator`

### Contributor Entry Points for Analysis

When adding a new analysis type or result class, start with these two base
class docstrings — they contain step-by-step instructions:

| Base Class | Location | What It Documents |
|------------|----------|-------------------|
| `BaseAnalyzer` | `analysis/core/registry.py` | How to add a new analyzer (compute, aggregate, register, CLI) |
| `BaseAnalysisResult` | `analysis/results/base.py` | Serialization contract (save/load, field conventions, migration) |

Key rules from those docstrings:

- **Results**: Inherit `BaseAnalysisResult`, set `analysis_type` as `ClassVar[str]`,
  implement `summary()`. Do NOT reimplement `save()`/`load()`.
- **Analyzers**: Implement `analysis_type()`, `from_config()`, `compute()`,
  `compute_aggregated()`, and a `label` property. Register settings with
  `@AnalysisSettingsRegistry.register()`.
- **Nested data objects** (e.g., per-residue stats) inherit `BaseModel`, not
  `BaseAnalysisResult`.
- **Large binary data** (e.g., per-frame SASA) uses NPZ + JSON sidecar instead
  of `BaseAnalysisResult`.

## Design Principles (Critical for Contributors)

This project prioritizes **extensibility** so users can contribute new analyses,
comparators, and plotters without modifying core code. Follow these principles:

### Open-Closed Principle (OCP)

Classes should be **open for extension, closed for modification**. Use:
- Abstract base classes (`ABC`) with well-defined contracts
- Registry patterns for runtime discovery of new implementations
- Strategy pattern for swappable algorithms

**Example:** To add a new plotter, inherit from `BasePlotter` and register
with `@PlotterRegistry.register()`. No changes to `ComparisonPlotter` needed.

### Follow Established Contracts

When extending a registry-based system, **study existing implementations first**:

1. **Read the base class docstrings** — they define the contract
2. **Study 2-3 existing implementations** — understand the expected data flow
3. **Match the pattern exactly** — don't invent new data passing mechanisms

**Anti-pattern to avoid:**
```python
# WRONG: Expecting custom kwargs that the orchestrator doesn't provide
def plot(self, data, labels, output_dir, **kwargs):
    result = kwargs.get("comparison_result")  # Orchestrator never passes this!
    ...
```

**Correct pattern:**
```python
# RIGHT: Load data from filesystem paths provided in data dict
def plot(self, data, labels, output_dir, **kwargs):
    for label in labels:
        analysis_dir = data[label]["analysis_dir"]
        result = MyResult.load(analysis_dir / "my_result.json")
    ...
```

### Registry Pattern Contracts

| Registry | Base Class | Data Source | Key Contract |
|----------|------------|-------------|--------------|
| `ComparatorRegistry` | `BaseComparator` | Load from filesystem via `_load_or_compute()` | Returns structured result object |
| `PlotterRegistry` | `BasePlotter` | Receives `data` dict with `analysis_dir` paths | Load your own data from `analysis_dir` |

### When Adding New Features

1. **Identify the extension point** — which registry/ABC to use?
2. **Read the base class** — understand required methods and their signatures
3. **Study existing implementations** — at least 2-3 similar classes
4. **Follow the data flow** — how does data get to your code?
5. **Test with the orchestrator** — don't just test in isolation

## Code Style

- **Formatter:** Black, line-length=100
- **Linter:** Ruff (see `pyproject.toml` for rule selection)
- **Docstrings:** NumPy style preferred (Google style exists in older modules)
- **Type hints:** `X | None` (3.10+ union syntax) in new code
- **Imports:** stdlib → third-party → local, lazy-import heavy deps

## Known Issues

1. **Config hash mismatch warning** prints 66+ times — should print once
2. **Contacts criteria mismatch** — cached 4.0A vs 4.5A cutoff disagreement
3. **Docs sidebar** — after adding toctree entries, run `make clean html` (not just `make html`)
4. **GitHub Issue #20** — tracks remaining analysis module TODOs

## Modular Instructions

See `.opencode/instructions/` for detailed rules on specific topics:
- `code-style.md` — formatting, linting, import conventions
- `architecture.md` — module structure, design patterns, extension points
- `environment.md` — conda setup, dependency management, CI
- `testing.md` — test infrastructure, running tests, writing new tests
- `analysis-module.md` — analysis-specific patterns and Issue #20 scope
- `documentation.md` — Sphinx/MyST conventions, API docs
- `known-issues.md` — detailed bug descriptions and workarounds
