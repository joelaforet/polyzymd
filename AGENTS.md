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
├── cli/          # Click CLI (main.py = entry point, scaffold.py = new-analysis generator)
├── config/       # Pydantic v2 models (schema.py), YAML loading
├── builders/     # System construction (PDB → parameterized topology)
├── simulation/   # OpenMM simulation runners
├── workflow/     # Orchestration (build → simulate → analyze)
├── core/         # Base classes, shared types
├── analyses/     # ★ Plugin system — unified analysis lifecycle (primary extension point)
│   ├── shared/   #   Reusable utilities (TrajectoryLoader, alignment, statistics, etc.)
│   └── <name>/   #   One package per analysis type (all plugins are packages)
├── exporters/    # GROMACS/other format exporters
├── data/         # Bundled data files (force fields, templates)
├── utils/        # Shared utilities
└── templates/    # Example YAML configs and project templates
```

### Inside `analyses/`

| Layer | Files | Role |
|-------|-------|------|
| **Plugins** (public) | `rmsf/`, `contacts/`, `distances/`, etc. | One class per analysis type — the **extension point** for contributors |
| **Private modules** | `<name>/_plotters.py`, `<name>/_results.py`, etc. | Plotting functions, result models, formatters, aggregators used internally by plugins |
| **Shared utilities** | `shared/loader.py`, `shared/alignment.py`, etc. | `TrajectoryLoader`, alignment, statistics, autocorrelation — reusable across plugins |
| **Shared compute** | `shared/binding_preference.py`, `shared/surface_exposure.py` | Cross-plugin compute (used by contacts, BFE, polymer_affinity) |
| **Framework** | `base.py`, `discovery.py`, `orchestrator.py`, `stats.py` | Plugin ABC, auto-discovery, lifecycle runner, default comparison utilities |

New analysis types are added as **packages in `analyses/<name>/`**. All 13
existing plugins are packages (no single-file plugins exist). Established
plugins extract plotting into `_plotters.py` modules; new plugins can keep
plotting inline in `plot()` or extract it as complexity grows.

## Key Patterns

- **Chain convention:** A=protein, B=substrate, C=polymer, D+=solvent
- **Factory pattern:** `ClassName.from_config(config)` or `ClassName.from_yaml(path)`
- **Lazy imports:** Heavy deps (OpenMM, MDAnalysis) imported inside functions/methods
- **ABC + Strategy:** `ContactCriteria`, `MolecularSelector`, `MoleculeCharger`
- **Plugin discovery:** `pkgutil`-based auto-discovery in `analyses/` — no registries
- **Config:** Pydantic v2 `BaseModel` subclasses with `model_validator`

### Contributor Entry Points for Analysis

To add a new analysis type, use the scaffold command or create a package
in `src/polyzymd/analyses/<name>/` and subclass `Analysis`:

| Resource | Location | What It Documents |
|----------|----------|-------------------|
| **Scaffold CLI** | `polyzymd new-analysis <name>` | Generates plugin package + tests automatically |
| `Analysis` base class | `analyses/base.py` | Full contract: required methods, optional overrides, context objects |
| Plugin discovery | `analyses/discovery.py` | How auto-discovery works, naming rules |
| Orchestrator | `analyses/orchestrator.py` | How the framework runs your plugin |
| Shared utilities | `analyses/shared/` | `TrajectoryLoader`, alignment, statistics, autocorrelation |
| Scaffold output | `polyzymd new-analysis <name>` | Simplest working plugin — start here |
| Richer example | `analyses/catalytic_triad/` | Default-compare lifecycle with DistanceCalculator + complex plotting |
| Stats utilities | `analyses/stats.py` | `default_scalar_comparison()`, `format_scalar_comparison()` |
| Contributor tutorial | `docs/source/contributor_guide/extending_analyses.md` | Step-by-step guide with test examples |

Key rules:

- **Required class variables**: `name` (str) and `Settings` (Pydantic BaseModel)
- **Lifecycle contract**: when `has_compute_stage=True`, use runner-backed `build_runner()` + `summarize_replicate()`; this keeps MDAnalysis-first per-trajectory iteration while PolyzyMD owns caching, ensemble aggregation, and comparison workflow. For compare-only plugins, set `has_compute_stage=False`. `aggregate(ctx, results)` is required only when `has_aggregate_stage=True`
- **Optional overrides**: `compare()`, `plot()`, `format()`, `extract_metrics()`, `filter_conditions()`
- **Default compare path**: Implement `extract_metrics()` — the framework loads results automatically (via `AggregatedResultClass` or `json.loads()`) and does t-tests, ANOVA, ranking
- **Custom compare path**: Override `compare()` entirely for multi-metric or entry-table analyses
- **Auto-discovery**: Drop a package in `analyses/<name>/` — no imports, no registries, no bootstrap
- **Result saving**: Existing plugins save results explicitly; the orchestrator has a fallback auto-save if the plugin doesn't
- **No `compare/` files needed**: Keep comparison and formatting logic in plugin packages, with shared helpers in `analyses/stats.py` and `analyses/shared/inferential_statistics.py`

### Planned MetricType classification (future)

The `MetricType` classification (`MEAN_BASED` and `VARIANCE_BASED`) is a
planned enhancement and is not implemented in the current codebase.

For now, reviewers should treat any `metric_type` checks as aspirational
guidance rather than a required plugin contract.

## Design Principles (Critical for Contributors)

This project prioritizes **extensibility** so users can contribute new analyses
without modifying core code. Follow these principles:

### Open-Closed Principle (OCP)

Classes should be **open for extension, closed for modification**. The plugin
system achieves this:
- Subclass `Analysis` and drop a package in `analyses/<name>/` — no core changes needed
- Framework discovers plugins automatically via `pkgutil`
- Default implementations (compare, format, plot) are overridable

### Follow Established Contracts

When writing a new analysis plugin, **study existing implementations first**:

1. **Read `analyses/base.py`** — it defines the full contract
2. **Start with the scaffold output** — `polyzymd new-analysis <name>` generates a complete working plugin with compute, aggregate, comparison, plotting, and tests
3. **Study `analyses/rmsf/`** for default compare with plots, or **`analyses/catalytic_triad/`** for default-compare lifecycle

**Anti-pattern to avoid:**
```python
# WRONG: Inventing custom data passing, bypassing the context
def build_runner(self, ctx, replicate, universe, window):
    config = SimulationConfig.from_yaml(self.custom_config_path)  # Don't do this!
```

**Correct pattern:**
```python
# RIGHT: Use the framework-provided context
def build_runner(self, ctx, replicate, universe, window):
    sim_config = ctx.sim_config  # Already loaded by framework
    settings = ctx.settings       # Your Settings model, resolved from YAML
```

### Plugin System Contracts

| Method | When Called | Input | Output |
|--------|-----------|-------|--------|
| `build_runner()` + `summarize_replicate()` | Once per replicate per condition on the runner-backed compute path; MDAnalysis owns per-trajectory iteration there while PolyzyMD owns ensemble workflow | `ReplicateContext` + replicate int | Runner + summarized replicate result |
| `aggregate()` | Once per condition after all replicates, only when `has_aggregate_stage=True` | `AggregateContext` + list of replicate results | Aggregated model or dict |
| `extract_metrics()` | During default `compare()` | Aggregated result | `dict[str, MetricValue]` |
| `compare()` | Once per analysis (cross-condition) | `ComparisonContext` | `ComparisonResult` or custom Pydantic model |
| `plot()` | Once per analysis | `PlotContext` | `list[Path]` of figures |
| `format()` | CLI display | Comparison result + format string | Formatted string |

### When Adding New Features

1. **Read the tutorial**: `docs/source/contributor_guide/extending_analyses.md`
2. **Read `analyses/base.py`** — the class docstring defines the full contract
3. **Pick your complexity level**: simple (use default compare) or custom (override compare)
4. **Study a matching example**: start with scaffold output (`polyzymd new-analysis <name>`), then use `rmsf/` for default compare with plots or `contacts/` for custom compare
5. **Write your plugin** in `analyses/<name>/` — keep plugin and lifecycle wiring in `__init__.py`; for trajectory-native runner-backed plugins, isolate MDAnalysis runner logic in a dedicated module such as `_runner.py`, and extract plotting to `_plotters.py` as complexity grows
6. **Test**: `pixi run -e build pytest tests/analyses/plugins/test_<name>.py -v`

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
5. **Pre-existing LSP type errors** — Pyright/Pylance reports false positives in `config/schema.py`, `builders/system_builder.py`, etc. due to missing type stubs for OpenMM/OpenFF. Does NOT affect runtime.

## Modular Instructions

See `.opencode/instructions/` for detailed rules on specific topics:
- `code-style.md` — formatting, linting, import conventions
- `architecture.md` — module structure, design patterns, extension points
- `environment.md` — pixi environment setup, dependency management, CI
- `testing.md` — test infrastructure, running tests, writing new tests
- `analysis-module.md` — analysis plugin system patterns and contracts
- `documentation.md` — Sphinx/MyST conventions, API docs, zero-warning build gate, `:no-index:` rules
- `known-issues.md` — detailed bug descriptions and workarounds
