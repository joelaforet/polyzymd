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
| **Plugins** (public) | `rg.py`, `rmsf/`, `contacts/`, etc. | One module per analysis, holding a settings model and a `compute()`. The extension point for contributors |
| **Private modules** | `_framework/`, `<name>/_*.py`, etc. | Internal framework and plugin implementation details; not contributor import targets |
| **Shared utilities** | `shared/loader.py`, `shared/alignment.py`, etc. | `TrajectoryLoader`, alignment, statistics, autocorrelation — reusable across plugins |
| **Framework** | `base.py`, `discovery.py`, `orchestrator.py`, `stats.py`, `mda/` | Stable public facade, auto-discovery, artifact lifecycle, default comparison utilities |

New analysis types may be simple single-file modules or packages under
`analyses/`. Every analysis is one module: a settings model, a `compute()` that
returns observables, and one call to `contract_analysis()`.

A plugin imports `Observable`, `iter_frames` and `contract_analysis` from
`polyzymd.analyses.contract`, and anything else it needs from
`polyzymd.analyses.shared`. `polyzymd.analyses.base` holds `Analysis` and the
four framework contexts (`ReplicateContext`, `AggregateContext`,
`ComparisonContext`, `PlotContext`), which a plugin receives rather than builds.
`polyzymd.analyses._framework` is private.

## Key Patterns

- **Chain convention:** A=protein, B=substrate, C=polymer, D+=solvent
- **OpenFF PDB ingestion:** New diagnosed protein/PDB ingestion failures must be
  documented in `docs/source/how_to/troubleshoot_openff_pdb_ingestion.md` and
  `docs/source/reference/openff_pdb_ingestion.md` before the task is closed,
  unless the user explicitly defers the durable documentation update.
- **OpenMM build identity:** A completed prebuild is committed by `build_manifest.json`. Never copy `solvated_system.pdb`,
  `system.xml`, or segment State/topology files independently. Continuation prefers the predecessor topology; root PDB fallback is legacy-only and must pass count validation.
- **OpenMM site runtime:** Run submission commands from `build`. Known-site
  presets pin one simulation environment per campaign: Blanca uses
  `sim-cuda-12-4`, and Bridges-2 uses `sim-cuda-12-6`. A node probe can reject
  and exclude a node, but it must not upgrade the environment on a newer driver.
- **Factory pattern:** `ClassName.from_config(config)` or `ClassName.from_yaml(path)`
- **Lazy imports:** Heavy deps (OpenMM, MDAnalysis) imported inside functions/methods
- **Preemption lifecycle:** Install handlers at `run-segment` entry; only atomic
  `phase.json` records with `status: completed` may skip OpenMM phases. Keep
  reporter intervals independent from bounded signal-check chunks.
- **ABC + Strategy:** `ContactCriteria`, `MolecularSelector`, `MoleculeCharger`
- **Plugin discovery:** `pkgutil`-based auto-discovery in `analyses/` — no registries
- **Config:** Pydantic v2 `BaseModel` subclasses with `model_validator`

### Contributor Entry Points for Analysis

To add a new analysis type, use the scaffold command or create a module/package
under `src/polyzymd/analyses/` and subclass `Analysis`:

| Resource | Location | What It Documents |
|----------|----------|-------------------|
| **Scaffold CLI** | `polyzymd new-analysis <name>` | Writes one plugin module and its two tests |
| `Analysis` base class | `analyses/base.py` | Stable public facade for the full contract, required methods, optional overrides, and context objects |
| Plugin discovery | `analyses/discovery.py` | How auto-discovery works, naming rules |
| Orchestrator | `analyses/orchestrator.py` | How the framework runs your plugin |
| Shared utilities | `analyses/shared/` | `TrajectoryLoader`, alignment, statistics, autocorrelation |
| Scaffold output | `polyzymd new-analysis <name>` | A working plugin to start from |
| Richer example | `analyses/catalytic_triad/` | Default-compare lifecycle with DistanceCalculator + complex plotting |
| Stats utilities | `analyses/stats.py` | `interpret_direction()`, `format_pct()` |
| Contributor guide | `docs/source/contributor_guide/analysis_plugins/index.md` | Write an analysis plugin, with the checklist |

Key rules:

A plugin declares `name` (str), `Settings` (a pydantic model), `references` (a
tuple of citations) and `compute(universe, frames, settings)`, and nothing else.

`contract_analysis()` builds the `Analysis` subclass that runs it. The framework
owns the universe, the production window, the replicate cache, `ArtifactStore`,
`ConditionArtifact`, `ComparisonArtifact`, aggregation, hypothesis testing,
plotting and formatting. There are no lifecycle hooks to override.

The `kind` of each observable does the work. `mean_of_timeseries`,
`fluctuation`, `fraction` and `profile` decide how a replicate reduces, how
conditions compare, and which figure is drawn.

Discovery walks `analyses/`, so a new module needs no import and no registry
entry. Artifacts persist through `ArtifactStore`; never introduce a
plugin-specific cache filename scheme.

## Design Principles (Critical for Contributors)

This project prioritizes **extensibility** so users can contribute new analyses
without modifying core code. Follow these principles:

### Open-Closed Principle (OCP)

Classes should be **open for extension, closed for modification**. The plugin
system achieves this:
- Write a plugin class and drop a module in `analyses/`, no core changes needed
- Framework discovers plugins automatically via `pkgutil`
- Behaviour comes from the observable `kind`, not from overriding hooks

### Follow Established Contracts

When writing a new analysis plugin, **study existing implementations first**:

1. Read `analyses/contract.py`, which defines the whole contract.
2. Start with the scaffold output. `polyzymd new-analysis <name>` writes a
   working plugin and two tests.
3. Study `analyses/rg.py` for a short plugin, or `analyses/contacts/` for one
   that also returns an extra sidecar.

**Anti-pattern to avoid:**
```python
# WRONG: loading config or averaging replicates inside compute()
def compute(self, universe, frames, settings):
    config = SimulationConfig.from_yaml(self.custom_config_path)  # Don't do this
```

**Correct pattern:**
```python
# RIGHT: measure the frames you were handed and return raw values
def compute(self, universe, frames, settings):
    group = universe.select_atoms(settings.selection)
    values = [measure(group) for _ in iter_frames(universe, frames)]
    return [Observable(name="my", kind="mean_of_timeseries", unit="A", values=values)]
```

### Plugin System Contracts

| What a plugin provides | When called | Input | Output |
|--------|-----------|-------|--------|
| `compute()` | Once per replicate per condition | universe, `FrameSelection`, settings | sequence of `Observable`, optionally with extra sidecar arrays |
| `identity_files()` (optional) | When the replicate cache is checked | settings | paths whose contents belong to the replicate identity |
| `PlotSettings` (optional) | When figures are drawn | n/a | a `ContractPlotSettings` subclass |

Everything else (persistence, aggregation, testing, plotting, formatting) is the
framework's, on `Analysis` in `analyses/base.py`.

### When Adding New Features

1. **Read the guide**: `docs/source/contributor_guide/analysis_plugins/index.md`
2. Read `analyses/contract.py`, whose module docstring defines the contract
3. **Pick your complexity level**: simple (use default compare) or custom (override compare)
4. Study a matching example: start with the scaffold output, then read `analyses/rg.py` for a short plugin or `analyses/contacts/` for one that returns an extra sidecar
5. Write your plugin as one module in `analyses/`. There is no plugin package and no place to put plotting or persistence code, because the framework owns both
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
