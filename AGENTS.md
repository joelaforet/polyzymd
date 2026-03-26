# PolyzyMD — Agent Instructions

> Computational toolkit for enzyme-polymer conjugate MD simulations.
> Python >=3.10 | MIT License | hatchling build | src layout

## Environment

**All simulation-stack commands MUST use a PolyzyMD pixi environment:**

```bash
pixi run -e <env> <command>
```

The PolyzyMD pixi environments contain OpenMM, OpenFF, MDAnalysis and other
heavy dependencies resolved from conda-forge. Never `pip install` these
outside the managed pixi environment.

**Quick commands:**

| Task | Command |
|------|---------|
| Install env | `pixi install -e build` |
| Activate shell | `pixi shell -e build` |
| Run tests | `pixi run -e build pytest tests/ -v` |
| Lint | `ruff check src/` |
| Format | `black src/ --check` (or `black src/` to fix) |
| Build docs | `pixi run -e build make -C docs clean html` |
| Type check | `pixi run -e build mypy src/polyzymd` |

## Git Workflow

- **Branches:** `main` (stable), `dev` (integration), `feature/*` (work)
- Currently on `refactor/analysis-ocp-compliance` (checked out from `v1.2.1`)
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
├── analysis/     # Per-condition analysis calculators (RMSF, contacts, etc.)
├── analyses/     # ★ Plugin system — unified analysis lifecycle (primary extension point)
├── compare/      # Statistics, formatters, plotters, config, IO
├── exporters/    # GROMACS/other format exporters
├── data/         # Bundled data files (force fields, templates)
├── utils/        # Shared utilities
└── configs/      # Default YAML configs
```

### `analysis/` vs `analyses/` (important distinction)

| Package | Role |
|---------|------|
| `analysis/` | Per-condition calculators, results, aggregation — the **compute** layer |
| `analyses/` | Plugin system — wraps `analysis/` calculators into a unified lifecycle (compute → aggregate → compare → plot → format) |

New analysis types are added as **plugins in `analyses/`**. The `analysis/`
package provides the underlying computation that plugins delegate to.

## Key Patterns

- **Chain convention:** A=protein, B=substrate, C=polymer, D+=solvent
- **Factory pattern:** `ClassName.from_config(config)` or `ClassName.from_yaml(path)`
- **Lazy imports:** Heavy deps (OpenMM, MDAnalysis) imported inside functions/methods
- **ABC + Strategy:** `ContactCriteria`, `MolecularSelector`, `MoleculeCharger`
- **Plugin discovery:** `pkgutil`-based auto-discovery in `analyses/` — no registries
- **Config:** Pydantic v2 `BaseModel` subclasses with `model_validator`

### Contributor Entry Points for Analysis

To add a new analysis type, create ONE file in `src/polyzymd/analyses/` and
subclass `Analysis`:

| Resource | Location | What It Documents |
|----------|----------|-------------------|
| `Analysis` base class | `analyses/base.py` | Full contract: required methods, optional overrides, context objects |
| Plugin discovery | `analyses/discovery.py` | How auto-discovery works, naming rules |
| Orchestrator | `analyses/orchestrator.py` | How the framework runs your plugin |
| Simplest example | `analyses/secondary_structure.py` | Uses default `compare()` — minimal override |
| Stats utilities | `analyses/stats.py` | `default_scalar_comparison()`, `format_scalar_comparison()` |

Key rules:

- **Required class variables**: `name` (str) and `Settings` (Pydantic BaseModel)
- **Required methods**: `compute_replicate(ctx, replicate)` and `aggregate(ctx, results)`
- **Optional overrides**: `compare()`, `plot()`, `format()`, `extract_metrics()`, `filter_conditions()`
- **Default compare path**: Implement `extract_metrics()` to return `dict[str, MetricValue]` — the framework does t-tests, ANOVA, ranking automatically
- **Custom compare path**: Override `compare()` entirely for multi-metric or entry-table analyses
- **Auto-discovery**: Drop a `.py` file in `analyses/` — no imports, no registries, no bootstrap

### Quick Example — Minimal Plugin

```python
"""Radius of gyration analysis plugin."""
from __future__ import annotations

from typing import Any, ClassVar, Sequence

from pydantic import BaseModel

from polyzymd.analyses.base import (
    AggregateContext, Analysis, MetricValue, ReplicateContext,
)


class RgAnalysis(Analysis):
    name: ClassVar[str] = "rg"

    class Settings(BaseModel):
        selection: str = "protein and name CA"

    def compute_replicate(self, ctx: ReplicateContext, replicate: int) -> Any:
        import MDAnalysis as mda
        import numpy as np

        sim_config = ctx.sim_config
        run_dir = sim_config.get_working_directory(replicate)
        topology = run_dir / "solvated_system.pdb"
        trajs = sorted(run_dir.glob("production_*/*_trajectory.dcd"))

        u = mda.Universe(str(topology), [str(t) for t in trajs])
        atoms = u.select_atoms(ctx.settings.selection)

        rg_values = [atoms.radius_of_gyration() for _ in u.trajectory]
        return {"mean_rg": float(np.mean(rg_values)), "replicate": replicate}

    def aggregate(self, ctx: AggregateContext, results: Sequence[Any]) -> Any:
        import numpy as np

        values = [r["mean_rg"] for r in results]
        return {"mean_rg": float(np.mean(values)), "sem_rg": float(np.std(values, ddof=1) / np.sqrt(len(values))), "replicate_values": values}

    def extract_metrics(self, summary: Any) -> dict[str, MetricValue]:
        return {
            "mean_rg": MetricValue(
                name="mean_rg",
                mean=summary["mean_rg"],
                sem=summary["sem_rg"],
                replicate_values=summary["replicate_values"],
                higher_is_better=False,  # lower Rg = more compact
                direction_labels=("compacting", "unchanged", "expanding"),
            )
        }
```

## Design Principles (Critical for Contributors)

This project prioritizes **extensibility** so users can contribute new analyses
without modifying core code. Follow these principles:

### Open-Closed Principle (OCP)

Classes should be **open for extension, closed for modification**. The plugin
system achieves this:
- Subclass `Analysis` and drop a file in `analyses/` — no core changes needed
- Framework discovers plugins automatically via `pkgutil`
- Default implementations (compare, format, plot) are overridable

### Follow Established Contracts

When writing a new analysis plugin, **study existing implementations first**:

1. **Read `analyses/base.py`** — it defines the full contract
2. **Study `analyses/secondary_structure.py`** or `analyses/rmsf.py` — simplest plugins
3. **Match the context pattern** — use `ctx.settings`, `ctx.sim_config`, etc.

**Anti-pattern to avoid:**
```python
# WRONG: Inventing custom data passing, bypassing the context
def compute_replicate(self, ctx, replicate):
    config = SimulationConfig.from_yaml(self.custom_config_path)  # Don't do this!
```

**Correct pattern:**
```python
# RIGHT: Use the framework-provided context
def compute_replicate(self, ctx, replicate):
    sim_config = ctx.sim_config  # Already loaded by framework
    settings = ctx.settings       # Your Settings model, resolved from YAML
```

### Plugin System Contracts

| Method | When Called | Input | Output |
|--------|-----------|-------|--------|
| `compute_replicate()` | Once per replicate per condition | `ReplicateContext` + replicate int | Pydantic model or dict |
| `aggregate()` | Once per condition (after all replicates) | `AggregateContext` + list of replicate results | Aggregated model or dict |
| `extract_metrics()` | During default `compare()` | Aggregated result | `dict[str, MetricValue]` |
| `compare()` | Once per analysis (cross-condition) | `ComparisonContext` | `ComparisonResult` or custom Pydantic model |
| `plot()` | Once per analysis | `PlotContext` | `list[Path]` of figures |
| `format()` | CLI display | Comparison result + format string | Formatted string |

### When Adding New Features

1. **Start with `analyses/base.py`** — read the class docstring
2. **Pick your complexity level**: simple (use default compare) or custom (override compare)
3. **Study a matching example**: `secondary_structure.py` for simple, `contacts.py` for custom
4. **Write your plugin** in `analyses/<name>.py`
5. **Test**: `pixi run -e build pytest tests/ -v -k <name>`

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
- `environment.md` — pixi environment setup, dependency management, CI
- `testing.md` — test infrastructure, running tests, writing new tests
- `analysis-module.md` — analysis plugin system patterns and contracts
- `documentation.md` — Sphinx/MyST conventions, API docs
- `known-issues.md` — detailed bug descriptions and workarounds
