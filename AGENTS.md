# PolyzyMD — Agent Instructions

> Computational toolkit for enzyme-polymer conjugate MD simulations.
> Python >=3.12 | MIT License | hatchling build | src layout

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

- **Branches:** `main` (released), `release/*` (release integration), short-lived
  `fix/*` and `feature/*` branches (atomic work)
- Resolve the current branch dynamically with `git branch --show-current`; do
  not rely on a static branch name in this file.
- Until the branch migration is complete, treat `feature/v1.3.0-rc5` as
  `release/1.3` and `conjugation-engine-refactor` as `release/1.4`.
- Mark beta snapshots with immutable SemVer prerelease tags such as
  `v1.3.0-rc.1`; do not create a new mutable branch for each snapshot.
- Forward-integrate stabilized v1.3 fixes into the v1.4 integration branch so
  conjugation work does not drift. Prefer merges for published/shared stacks
  and rebases for unpublished local work.
- Use conventional commits with an imperative subject near 50 characters.
  Keep each commit atomic and explain important reasoning in the body.
- Run `ruff check` and `black --check` before committing
- Never force-push to `main`, `dev`, or `release/*`.
- Never merge a pull request. Joe reviews every PR and merges it manually.

See `.opencode/instructions/development-workflow.md` for collaboration,
scope, validation, authorship, push, and PR rules.

## Harness Capabilities and Living Guidance

The GitHub connector lacked PR-write permission, so I used the authenticated gh fallback. I will not merge it.

Treat `AGENTS.md`, `.opencode/instructions/`, and the personal PolyzyMD Skills
as living guidance. When repository behavior, scientific contracts, branch
names, tools, permissions, or recurring workflows change, update the affected
guidance in the same atomic task so future agents do not follow stale rules.
Validate edited Skills and run the relevant repository documentation or static
checks before committing.

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
│   └── <name>/   #   Analysis plugins (single-file simple modules or packages)
├── exporters/    # GROMACS/other format exporters
├── data/         # Bundled data files (force fields, templates)
├── utils/        # Shared utilities
└── templates/    # Example YAML configs and project templates
```

### Inside `analyses/`

| Layer | Files | Role |
|-------|-------|------|
| **Plugins** (public) | `rmsf/`, `contacts/`, `distances/`, etc. | One class per analysis type — the **extension point** for contributors |
| **Private modules** | `_framework/`, `<name>/_*.py`, etc. | Internal framework and plugin implementation details; not contributor import targets |
| **Shared utilities** | `shared/loader.py`, `shared/alignment.py`, etc. | `TrajectoryLoader`, alignment, statistics, autocorrelation — reusable across plugins |
| **Framework** | `base.py`, `discovery.py`, `orchestrator.py`, `stats.py`, `mda/` | Stable public facade, auto-discovery, artifact lifecycle, default comparison utilities |

New analysis types may be simple single-file modules or packages under
`analyses/`. The scaffold defaults to a compact single-file plugin; advanced
and built-in analyses commonly use packages with `_mda.py` and `_plotters.py`
helpers as complexity grows.

`polyzymd.analyses.base` is the stable public facade for contributor imports.
It re-exports context objects and comparison models while delegating
implementation to private `_framework/` modules such as `compare.py`,
`io.py`, `contract.py`, `contexts.py`, and `comparison_models.py`. Contributors should import
`Analysis`, `ReplicateContext`, `AggregateContext`, `ComparisonContext`,
`PlotContext`, `MetricValue`, and `ComparisonResult` from
`polyzymd.analyses.base`, not from private modules.

## Key Patterns

- **Chain convention:** A=protein, B=substrate, C=polymer, D+=solvent
- **OpenFF PDB ingestion:** New diagnosed protein/PDB ingestion failures must be
  documented in `docs/source/how_to/troubleshoot_openff_pdb_ingestion.md` and
  `docs/source/reference/openff_pdb_ingestion.md` before the task is closed,
  unless the user explicitly defers the durable documentation update.
- **Factory pattern:** `ClassName.from_config(config)` or `ClassName.from_yaml(path)`
- **Lazy imports:** Heavy deps (OpenMM, MDAnalysis) imported inside functions/methods
- **ABC + Strategy:** `ContactCriteria`, `MolecularSelector`, `MoleculeCharger`
- **Plugin discovery:** `pkgutil`-based auto-discovery in `analyses/` — no registries
- **Config:** Pydantic v2 `BaseModel` subclasses with `model_validator`

### Contributor Entry Points for Analysis

To add a new analysis type, use the scaffold command or create a module/package
under `src/polyzymd/analyses/` and subclass `Analysis`:

| Resource | Location | What It Documents |
|----------|----------|-------------------|
| **Scaffold CLI** | `polyzymd new-analysis <name>` | Generates a simple plugin, or an advanced package with tests |
| `Analysis` base class | `analyses/base.py` | Stable public facade for the full contract, required methods, optional overrides, and context objects |
| Plugin discovery | `analyses/discovery.py` | How auto-discovery works, naming rules |
| Orchestrator | `analyses/orchestrator.py` | How the framework runs your plugin |
| Shared utilities | `analyses/shared/` | `TrajectoryLoader`, alignment, statistics, autocorrelation |
| Scaffold output | `polyzymd new-analysis <name>` | Simplest working plugin — start here |
| Richer example | `analyses/catalytic_triad/` | Default-compare lifecycle with DistanceCalculator + complex plotting |
| Stats utilities | `analyses/stats.py` | `default_scalar_comparison()`, `format_scalar_comparison()` |
| Contributor tutorial | `docs/source/contributor_guide/extending_analyses.md` | Step-by-step guide with test examples |

Key rules:

- **Required class variables**: `name` (str) and `Settings` (Pydantic BaseModel)
- **Lifecycle contract**: trajectory-native analyses build `MDAAnalysisJob` objects wrapping `AnalysisBase`-compatible work and map completed jobs through collectors into `ReplicateArtifact` objects. PolyzyMD owns `ArtifactStore`, `ConditionArtifact`, `ComparisonArtifact`, ensemble aggregation, statistics, and plotting. For compare-only plugins, set `has_compute_stage=False`. `aggregate(ctx, results)` is required only when `has_aggregate_stage=True`
- **Optional overrides**: `compare()`, `plot()`, `format()`, `extract_metrics()`, `filter_conditions()`
- **Default compare path**: Implement `extract_metrics()` — the framework loads results automatically (via `AggregatedResultClass` or `json.loads()`) and does t-tests, ANOVA, ranking
- **Custom compare path**: Override `compare()` entirely for multi-metric or entry-table analyses
- **Auto-discovery**: Drop a module or package in `analyses/` — no imports, no registries, no bootstrap
- **Result saving**: Prefer canonical artifacts through `ArtifactStore`; do not introduce plugin-specific cache filename schemes
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
- Subclass `Analysis` and drop a module or package in `analyses/` — no core changes needed
- Framework discovers plugins automatically via `pkgutil`
- Default implementations (compare, format, plot) are overridable

### Follow Established Contracts

When writing a new analysis plugin, **study existing implementations first**:

1. **Read `analyses/base.py`** — it defines the full contract
2. **Start with the scaffold output** — `polyzymd new-analysis <name>` generates a complete working plugin with MDAnalysis jobs, artifacts, aggregation, comparison, plotting, and tests
3. **Study `analyses/rmsf/`** for default compare with plots, or **`analyses/catalytic_triad/`** for default-compare lifecycle

**Anti-pattern to avoid:**
```python
# WRONG: Inventing custom data passing, bypassing the context
def build_mda_jobs(self, ctx):
    config = SimulationConfig.from_yaml(self.custom_config_path)  # Don't do this!
```

**Correct pattern:**
```python
# RIGHT: Use the framework-provided context
def build_mda_jobs(self, ctx):
    sim_config = ctx.sim_config  # Already loaded by framework
    settings = ctx.settings       # Your Settings model, resolved from YAML
```

### Plugin System Contracts

| Method | When Called | Input | Output |
|--------|-----------|-------|--------|
| `build_mda_jobs()` + `build_mda_collector()` | Once per replicate per condition on the MDAnalysis job compute path; MDAnalysis owns per-trajectory iteration there while PolyzyMD owns artifacts and ensemble workflow | `MDAReplicateJobContext` / collector context | `ReplicateArtifact` |
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
5. **Write your plugin** as a simple module or package in `analyses/`; for advanced trajectory-native packages, isolate MDAnalysis job helpers in `_mda.py`, and extract plotting to `_plotters.py` as complexity grows
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
- `development-workflow.md` — release flow, atomic scope, agent
  collaboration, commits, pushes, and PR handoff
- `analysis-module.md` — analysis plugin system patterns and contracts
- `documentation.md` — Sphinx/MyST conventions, API docs, zero-warning build gate, `:no-index:` rules
- `openff-pdb-ingestion.md` — OpenFF protein/PDB ingestion troubleshooting and living error-log rules
- `known-issues.md` — detailed bug descriptions and workarounds
