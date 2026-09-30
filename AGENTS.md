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
├── analyses/     # ★ Study API: functions of a replicate Universe, Study, polyzymd analyze (primary extension point)
│   ├── shared/   #   Reusable utilities (TrajectoryLoader, alignment, statistics, etc.)
│   └── _framework/, mda/  # Plugin framework: no shipped analysis uses it; being removed
├── exporters/    # GROMACS/other format exporters
├── data/         # Bundled data files (force fields, templates)
├── utils/        # Shared utilities
└── templates/    # Example YAML configs and project templates
```

### Inside `analyses/`

| Layer | Files | Role |
|-------|-------|------|
| **Function analyses** (the extension point) | `functions.py`, `study.py`, `timeseries.py`, `protocols.py`, `figures.py` | rg, rmsd, rmsf, rmsd_per_residue, distances, sasa, secondary_structure, contacts, native_contacts and hydrogen_bonds: functions of a replicate `Universe`, measured by `Study.timeseries` or `Study.per_replicate` and run by `polyzymd analyze NAME -c config.yaml` (`protocols.FUNCTION_ANALYSES`). Before choosing among `rmsd`, `rmsf`, `offset` and `rmsd_per_residue`, read the "Fluctuation, offset and deviation" section of `docs/source/explanation/analysis_rmsf_best_practices.md`: `rmsd` is per frame, the others per residue, and `rmsd_per_residue² = rmsf² + offset²` |
| **Shared utilities** | `shared/loader.py`, `shared/window.py`, etc. | `TrajectoryLoader`, frame windows, statistics, autocorrelation |
| **Plugin framework** | `base.py`, `discovery.py`, `orchestrator.py`, `stats.py`, `mda/`, `_framework/` | The `Analysis` plugin lifecycle behind `polyzymd compare`. No shipped analysis is a plugin any more, and the framework is being removed; do not add plugins |

The catalytic triad is not an analysis: it is a routine on the study API
(`docs/source/how_to/analysis_triad_quickstart.md`), and `polyzymd analyze
catalytic_triad` exits with a pointer to it. `polyzymd analyze` does not read
`comparison.yaml`: `-f comparison.yaml` exits with the equivalent `-c` command.

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
- **ABC + Strategy:** `MolecularSelector`, `MoleculeCharger`
- **Analyses are functions:** a measurement is a function of an MDAnalysis `Universe`; `Study` supplies the replicate universes, records and cross-replicate statistics
- **Config:** Pydantic v2 `BaseModel` subclasses with `model_validator`

### Contributor Entry Points for Analysis

To add a measurement, write a function of an MDAnalysis `Universe` (per frame,
or per replicate with a `frames` argument) and run it with `Study.timeseries`
or `Study.per_replicate`. Add it to `polyzymd analyze` only when it is a
shipped analysis, as a `FUNCTION_ANALYSES` entry and an `_analyze_<name>` in
`protocols.py`.

| Resource | Location | What It Documents |
|----------|----------|-------------------|
| Shipped functions | `analyses/functions.py`, `docs/source/reference/analysis_functions.md` | What each function measures, its arguments and units |
| Study API | `analyses/study.py`, `analyses/timeseries.py` | `Study.from_configs`, `timeseries`, `per_replicate`, `transform`, `reduce`, `compare`, records under `polyzymd_results/` |
| API explanation | `docs/source/explanation/analysis_api.md` | How the study API supplies universes, records and statistics |
| `polyzymd analyze` | `analyses/protocols.py`, `docs/source/how_to/analysis_agent_protocol.md` | `FUNCTION_ANALYSES`, the `_analyze_<name>` functions, `ProtocolReport` |
| Figures | `analyses/figures.py` | `ReplicateValues.plot`, profiles, differences, uncertainty footnotes |
| Worked routine | `docs/source/how_to/analysis_triad_quickstart.md` | Combining shipped functions for a question of your own |

The plugin framework (`Analysis`, `polyzymd new-analysis`, `polyzymd compare`)
still exists but no shipped analysis uses it, and it is being removed. Do not
write new plugins; `docs/source/contributor_guide/analysis_plugins/` documents
it only until the removal.

## Design Principles (Critical for Contributors)

This project prioritizes **extensibility** so users can contribute new analyses
without modifying core code. Follow these principles:

### Open-Closed Principle (OCP)

A new measurement is a new function passed to `Study.timeseries` or
`Study.per_replicate`; the study, its caching, its records and its statistics
need no change. Do not reimplement what MDAnalysis, MDTraj or another package
already computes: wrap it in a function.

### Follow Established Contracts

When writing a new analysis, **study existing implementations first**:

1. **Read `analyses/functions.py`**: `radius_of_gyration` is a per-frame function, `hydrogen_bonds` a per-replicate one that takes `frames`
2. **Read the matching `_analyze_<name>` in `analyses/protocols.py`** to see how a function becomes a `polyzymd analyze` result, with its `--run` names and figures
3. **Read `docs/source/how_to/analysis_triad_quickstart.md`** for a routine that combines shipped functions with `Timeseries.transform`

**Anti-pattern to avoid:**
```python
# WRONG: your own loop over replicates and frames
for replicate in range(1, 4):
    universe = mda.Universe(...)
    values = [my_measure(universe) for ts in universe.trajectory]
```

**Correct pattern:**
```python
# RIGHT: the study supplies the universes, the production frames and the records
series = study.timeseries(my_measure, pz.select("protein"), unit="A", name="my_measure")
print(series.reduce("mean").compare().to_agent_text())
```

### When Adding New Features

1. **Read** `docs/source/explanation/analysis_api.md` and `docs/source/reference/analysis_functions.md`
2. **Write the function** in `analyses/functions.py`, with a docstring that says concretely what it measures and in which unit
3. **Test it on placed geometry** with a known answer, then on the synthetic OpenMM run directories of `tests/_support/analysis_testkit.py` (`write_simulation_config`, `write_openmm_replicate`)
4. **For a shipped analysis**, add the `FUNCTION_ANALYSES` entry and `_analyze_<name>` in `protocols.py`, and document it in `analysis_functions.md` and a how-to page
5. **Test**: `pixi run -e test pytest tests/analyses -v -k <name>`

## Frame counting in analyses

The study API (`Study`, `Replicate.frames`, `pz.reference`) counts frames in
two ways, and mixing them up picks the wrong structure:

- `Replicate.frames` holds trajectory frame indices, counted from 0 from the
  first loaded frame, for the production frames left after the equilibration
  window.
- A frame the user names, such as `pz.reference("frame", ..., frame=N)` or
  `--set reference_frame=N` for `polyzymd analyze rmsd`, is a production frame
  counted from 1 after the window. The default, 1, is the first production
  frame, so the same N points at a different structure when the window
  changes. Each record stores the trajectory frame actually used under
  `chosen`.
- The removed plugins counted differently: rmsd took trajectory frames from 0
  (default 0, inside the equilibration window), and rmsf took trajectory
  frames from 1. With the first production frame's 0-based trajectory index,
  convert an old rmsd value with new = old - (first production frame) + 1 and
  an old rmsf value with new = old - (first production frame).
- A frame inside the equilibration window, such as the starting structure,
  cannot be named as a `frame` reference. Use `external` mode with that
  structure's file, which is also hashed into the record.

## Code Style

- **Formatter:** Black, line-length=100
- **Linter:** Ruff (see `pyproject.toml` for rule selection)
- **Docstrings:** NumPy style preferred (Google style exists in older modules)
- **Type hints:** `X | None` (3.10+ union syntax) in new code
- **Imports:** stdlib → third-party → local, lazy-import heavy deps

## Known Issues

1. **Config hash mismatch warning** prints 66+ times — should print once
2. **Docs sidebar** — after adding toctree entries, run `make clean html` (not just `make html`)
3. **GitHub Issue #20** — tracks remaining analysis module TODOs
4. **Pre-existing LSP type errors** — Pyright/Pylance reports false positives in `config/schema.py`, `builders/system_builder.py`, etc. due to missing type stubs for OpenMM/OpenFF. Does NOT affect runtime.

## Modular Instructions

See `.opencode/instructions/` for detailed rules on specific topics:
- `code-style.md` — formatting, linting, import conventions
- `architecture.md` — module structure, design patterns, extension points
- `environment.md` — pixi environment setup, dependency management, CI
- `testing.md` — test infrastructure, running tests, writing new tests
- `development-workflow.md` — release flow, atomic scope, agent
  collaboration, commits, pushes, and PR handoff
- `analysis-module.md` — study API patterns, `polyzymd analyze` protocols and the plugin framework being removed
- `documentation.md` — Sphinx/MyST conventions, API docs, zero-warning build gate, `:no-index:` rules
- `openff-pdb-ingestion.md` — OpenFF protein/PDB ingestion troubleshooting and living error-log rules
- `known-issues.md` — detailed bug descriptions and workarounds
