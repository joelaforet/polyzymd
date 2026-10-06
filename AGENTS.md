# PolyzyMD — Agent Instructions

> Computational toolkit for enzyme-polymer conjugate MD simulations.
> Python >=3.12 | MIT License | hatchling build | src layout

## Environment

Run every command that needs the simulation or analysis stack in a PolyzyMD
pixi environment: `pixi run -e <env> <command>`. Pixi 0.72.2 or newer reads
`pixi.toml`. The default environment has no numpy, so analysis code and tests
fail there.

| Task | Environment | Command |
|------|-------------|---------|
| Run tests | `test` | `pixi run -e test pytest tests/analyses -q` |
| Run analysis or library code | `analysis` | `pixi run -e analysis python -c "import polyzymd"` |
| Lint | `build` | `pixi run -e build ruff check src tests` |
| Format check | `build` | `pixi run -e build black --check <changed files>` |
| Build docs | `build` | `pixi run -e build make -C docs clean html` |
| Submit OpenMM jobs | `build` | see "OpenMM site runtime" below |

- Do not run `pixi install` or edit `pixi.toml` during a task. Re-solving the
  environments is slow and is the maintainer's job.
- Never `pip install` OpenMM, OpenFF, AmberTools, parmed or PDBFixer. They come
  from conda-forge through pixi. Do not add OpenMM or OpenFF to
  `pyproject.toml`; pip-installable analysis extras go in
  `[project.optional-dependencies]`.

### Worktrees

A git worktree shares the main checkout's pixi environments, and their
editable install points at the main checkout's `src`. Set `PYTHONPATH` to the
worktree's `src` and call the environment's interpreter directly:

```bash
MAIN=$(dirname "$(git rev-parse --path-format=absolute --git-common-dir)")
PYTHONPATH=$PWD/src $MAIN/.pixi/envs/test/bin/python -m pytest tests/analyses -q
$MAIN/.pixi/envs/build/bin/ruff check src tests
```

Before you trust a test run from a worktree, print `polyzymd.__file__` and
check that it points inside the worktree.

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
- Commit as Joe (`git config user.name`/`user.email`). Do not add
  generated-by text or co-author trailers.
- Run `ruff check` and `black --check` before committing.
- Never force-push to `main`, `dev`, or `release/*`. Force-push any other
  shared branch only with Joe's approval, and then with `--force-with-lease`.
- Do not commit or push failing, incomplete or scientifically unresolved work.
- Open pull requests as drafts. Joe reviews and merges them; never merge
  yourself. State the user goal, scientific assumptions, validation and
  limitations in the description.
- Never merge a pull request. Joe reviews every PR and merges it manually.

### Priorities and scope

Priorities, in order: scientific correctness and reproducibility, low-friction
user experience, OpenMM and GROMACS parity, small reviewable changes, docs that
match the code. Ask Joe before you choose a chemistry, topology,
parameterization, force-field or interpretation assumption.

- One end-user goal per PR. Above about 300 changed lines (tests included),
  look for scope creep and split independent behavior. Above 500, stop and
  explain why the change cannot be split. An "and" in a subject is a signal
  to check for two changes.
- Use one writing agent per worktree. Parallel writing tasks use separate
  worktrees.
- Run cheap checks first: static checks, focused tests, subsystem tests, then
  integration and scientific checks. GPU simulations are for feature
  acceptance and release gates.
- Engine-facing features need OpenMM and GROMACS coverage with the same
  particles, masses, bonds, constraints, charges, parameters, exclusions,
  coordinates and box. Document the tolerance for any engine difference, such
  as small PME energy differences.

## Living Guidance

`AGENTS.md`, `CLAUDE.md`, `.claude/skills/` and the personal PolyzyMD Skills
are living guidance. When code, scientific contracts, commands, branch names,
tools or recurring workflows change, update the affected guidance in the same
task. Remove obsolete rules instead of adding contradicting ones. If you cannot
update it, give Joe the stale file and the replacement text.

## Architecture Quick Reference

```
src/polyzymd/
├── cli/          # Click CLI (main.py = entry point, analyze.py = polyzymd analyze, retired.py = hidden compare/new-analysis stubs)
├── config/       # Pydantic v2 models (schema.py), YAML loading
├── builders/     # System construction (PDB → parameterized topology)
├── simulation/   # OpenMM simulation runners
├── workflow/     # Orchestration (build → simulate), SLURM scripts, analysis_submit.py = analyze --submit
├── core/         # Base classes, shared types
├── analyses/     # ★ Study API: functions of a replicate Universe, Study, polyzymd analyze (primary extension point)
│   ├── universe.py, identity.py  #   UniverseProvider (replicate universes, input file records), compute_config_hash
│   └── shared/   #   Reusable utilities (TrajectoryLoader, alignment, statistics, etc.)
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
| **Loading and identity** | `universe.py`, `identity.py` | `UniverseProvider` loads each replicate's `Universe` and records its input files; `compute_config_hash` is recorded by every stored result and must not change |

The catalytic triad is not an analysis: it is a routine on the study API
(`docs/source/how_to/analysis_triad_quickstart.md`), and `polyzymd analyze
catalytic_triad` exits with a pointer to it. `polyzymd analyze` does not read
`comparison.yaml`: `-f comparison.yaml` exits with the equivalent `-c` command.
An unknown name exits 2 listing `FUNCTION_ANALYSES`. `polyzymd analyze NAME -c
... --submit --preset <cluster>` runs one SLURM array task per condition and
replicate and a report job (`workflow/analysis_submit.py`). The hidden
`polyzymd compare` and `polyzymd new-analysis` commands (`cli/retired.py`) only
exit 2 and name their replacements.

A study folder's `study.yaml` (`analyses/study_file.py`) holds the analysis
protocol: conditions, the equilibration window (an entry may set its own
`equilibration:`/`until:`, `StudyFile.window`) and each analysis run's
settings, including the study's own functions (`function: file.py:name`,
`analyses/user_functions.py`, keyed on the whole file's hash).
A study is one protein (or system); a project (`project.yaml`,
`analyses/project_file.py`, `analyses/project.py`, `cli/project.py`) lists a
paper's studies and the analyses each runs, with each study's `regions:` and
`structures:` resolved into `region <name>` / `structure <name>`
(`resolve_names`); `analyze --project`, `project check` and
`pz.Project(...).results(run)` (a `study` column, factor columns).
Statistics (`analyses/study_statistics.py`): `replicate_table(run)` (one row
per replicate, the test unit) and a slope test per numeric condition factor
in every `--study` report (`TrendReport`). A paper's own statistics are
scripts in the project's `stats/` folder; freeze hashes and publishes them.
`polyzymd project init` (`analyses/project_scaffold.py`) writes a project
with empty studies (existing studies are moved in by hand, following
`docs/source/how_to/move_studies_into_project.md`), and `polyzymd project freeze`
(`analyses/project_freeze.py`) freezes every study with
`freeze(..., publish=False)` and publishes the project once. File arguments
are recorded by name and SHA-256, never location, so moved studies reuse
results.
`polyzymd study check` reads it without trajectories (`cli/study.py`);
`polyzymd analyze [RUN] --study study.yaml` runs one run or all of them into
`<study>/results/<run>/`; `pz.Study("study.yaml").results(run)` reads them back
without trajectories (`analyses/results.py`). `polyzymd study init` writes a
study folder and commits it (`analyses/study_scaffold.py`); a gitignored
`data.local.yaml`, written by `polyzymd study locate DIR`, or `--data DIR` says
where this machine keeps each condition's runs without changing its config
hash; reports record the study's git state (`analyses/study_git.py`).
`polyzymd study add-condition` adds a condition to an existing study. Reports
warn when conditions were analysed up to different times (`until`/`--until`
gives a common window) and name segments `progress.json` records but the disk
lacks; `analyze` and the `study` commands print only reports and warnings and
write the full log to `logs/` (`polyzymd -v` for more). A study's own function
with `allow_empty: true` gets an empty AtomGroup where a selection matches
nothing (a polymer selection in a no-polymer control), so it must handle one.
Config hashes identify input structures by content and leave out the
projects and scratch directories (`analyses/identity.py`); stored records
identify trajectory and topology files by SHA-256 and size
(`analyses/shared/file_hashes.py`, from `progress.json`'s `trajectory_sha256`
when the runner recorded it), so moved or downloaded data reuses results.
`polyzymd hash-trajectories` (`cli/hashes.py`) records missing trajectory
hashes for older runs of any engine in `trajectory_hashes.json`,
idempotently, without ever overwriting a recorded hash, and never writing
`progress.json` (only the runner writes it). The logic is
`SimulationEngine.record_trajectory_hashes` in `engines/base.py`; each engine
overrides `trajectory_files` with the finished files analyses read (and may
extend `recorded_trajectory_hashes`, as OpenMM does so the runner's segment
hash wins), and records and freeze
read hashes through `UniverseProvider.recorded_trajectory_hashes`. Anything
new about a run's outputs belongs on the engine in the same way.
`polyzymd study freeze` (`analyses/study_freeze.py`, `analyses/study_metadata.py`)
checks staleness and metadata, writes `manifest.json`, `CITATION.cff` (citing
PolyzyMD through `polyzymd/citation.py`), `.zenodo.json`, `md_checklist.yaml` and
`system_summary.csv`, commits and tags them with `results/`, and lays out the
gitignored `deposit/`, with `deposit/upload/`, `deposit/trajectories.csv` and
the step-by-step `deposit/UPLOAD.md` (`analyses/study_upload_guide.py`).
PolyzyMD never uploads or publishes: that is the author's step. The design of study folders and
the slices still to come are in `docs/source/explanation/study_folders.md`.

## Key Patterns

- **Chain convention:** A=protein, B=substrate, C=polymer, D+=solvent
- **OpenFF PDB ingestion:** New diagnosed protein/PDB ingestion failures must be
  documented in `docs/source/how_to/troubleshoot_openff_pdb_ingestion.md` and
  `docs/source/reference/openff_pdb_ingestion.md` before the task is closed,
  unless the user explicitly defers the durable documentation update.
  Reproduce a failure with `openff.toolkit.Topology.from_pdb(path)` first.
  Treat a charge mismatch as a chemistry or connectivity error: fix the
  structure upstream. Never monkeypatch OpenFF or weaken its validation.
- **OpenMM build identity:** A completed prebuild is committed by `build_manifest.json`. Never copy `solvated_system.pdb`,
  `system.xml`, or segment State/topology files independently. Continuation prefers the predecessor topology; root PDB fallback is legacy-only and must pass count validation.
- **OpenMM site runtime:** Run submission commands from `build`. Known-site
  presets pin one simulation environment per campaign: Blanca uses
  `sim-cuda-12-4`, and Bridges-2 uses `sim-cuda-12-6`. A node probe can reject
  and exclude a node, but it must not upgrade the environment on a newer driver.
- **Factory pattern:** `ClassName.from_config(config)` or `ClassName.from_yaml(path)`
- **Trajectory loading:** `analyses/shared/loader.py` (`TrajectoryLoader`,
  `open_universe`) finds and opens each replicate's files, and
  `analyses/universe.py` (`UniverseProvider`) wraps it. Do not build a
  replicate `Universe` or its file paths by hand elsewhere.
- **Preemption lifecycle:** Install handlers at `run-segment` entry; only atomic
  `phase.json` records with `status: completed` may skip OpenMM phases. Keep
  reporter intervals independent from bounded signal-check chunks.
- **ABC + Strategy:** `MoleculeCharger`
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
| Study API | `analyses/study.py`, `analyses/timeseries.py`, `docs/source/how_to/study_api.md`, `docs/source/reference/study_api.md` | `Study.from_configs`, `timeseries`, `per_replicate`, `transform`, `reduce`, `compare`, records under `polyzymd_results/` |
| API explanation | `docs/source/explanation/analysis_api.md` | How the study API supplies universes, records and statistics |
| `polyzymd analyze` | `analyses/protocols.py`, `docs/source/how_to/analysis_agent_protocol.md` | `FUNCTION_ANALYSES`, the `_analyze_<name>` functions, `ProtocolReport` |
| Study folders | `analyses/study_file.py`, `analyses/results.py`, `analyses/study_scaffold.py`, `analyses/study_git.py`, `docs/source/how_to/study_yaml.md`, `docs/source/how_to/study_folder.md` | `study.yaml`, `--study`, `polyzymd study check/init/locate/freeze`, `data.local.yaml`, `--data`, `Study.results`, `docs/source/how_to/study_freeze.md` |
| Figures | `analyses/figures.py` | `ReplicateValues.plot`, profiles, differences, uncertainty footnotes |
| Worked routine | `docs/source/how_to/analysis_triad_quickstart.md` | Combining shipped functions for a question of your own |

`docs/source/contributor_guide/adding_an_analysis.md` gives the contributor
steps. There is no plugin class, registry or scaffold.

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

1. **Read** `docs/source/reference/study_api.md`, `docs/source/explanation/analysis_api.md` and `docs/source/reference/analysis_functions.md`
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
- A frame inside the equilibration window, such as the starting structure,
  cannot be named as a `frame` reference. Use `external` mode with that
  structure's file, which is also hashed into the record.

## Code Style

- **Formatter:** Black, line-length 100. **Linter:** Ruff (rules in
  `pyproject.toml`). CI also runs `ruff format --check src/polyzymd/`.
- Run black only on the files your change touches. Never run `black src/` or
  format whole files that the change does not otherwise edit: it buries the
  change in unrelated diffs.
- **Docstrings:** NumPy style in new code. Older modules (`builders/`,
  `simulation/`, `workflow/`) use Google style; convert one only when you
  rewrite most of that function.
- **Type hints:** `X | None` in new code. Do not change `Optional[X]` in code
  you do not otherwise touch.
- **Imports:** stdlib, third-party, local. Import `openmm`, `openff.*`,
  `MDAnalysis`, `parmed` and `pdbfixer` inside functions, never at module
  level, so the config, CLI and docs build work without the full stack.

## Testing

- `tests/` mirrors `src/polyzymd/`: put the tests for `analyses/foo.py` in
  `tests/analyses/test_foo.py`. Shared fixtures are in `tests/conftest.py`;
  synthetic OpenMM run directories come from `tests/_support/analysis_testkit.py`.
- A test must not leave files in the repository. Write to `tmp_path`, and use
  `monkeypatch.chdir(tmp_path)` for code that writes to the working directory.
- Do not commit trajectories (`.xtc`, `.dcd`) or large PDB files. Use small
  synthetic data or placed geometry with a known answer.
- Prefer the smallest test set that proves the behavior. Do not add layers of
  fixtures or mocks without a concrete failure they catch.
- Get the test count from `pytest tests --collect-only -q | tail -1`; do not
  copy a number into docs.

## Documentation

Sphinx with MyST Markdown, sources in `docs/source/`. Update the docs in the
same change when code changes a command, config field, output or scientific
meaning.

- Keep each page to one Diátaxis role: `get_started/` and `tutorials/` teach,
  `how_to/` gives task steps, `reference/` lists facts, `explanation/` gives
  concepts and reasons. Link between pages instead of duplicating content.
- Keep the admonitions that tell readers which pixi environment to use.
- Use `docs/source/get_started/quickstart.md` as the voice reference: short
  sentences, exact commands, expected outputs, next steps.
- Name repository paths and public symbols, or use MyST roles
  (`{func}`, `{class}`, `{doc}`). Do not link to source line numbers.
- CI builds the docs with `-W`: zero warnings. After you add or remove a
  `toctree` entry, run `make clean html`, not `make html`, or the sidebar is
  stale. Every new page goes in a `toctree`.
- Add `:no-index:` to an `automodule` directive for a module with Pydantic
  models or dataclasses, or Sphinx reports duplicate object descriptions.
- Heading levels must be consecutive. In docstrings, put a blank line before a
  list, and do not list Pydantic fields under `Attributes:` (autodoc-pydantic
  documents them).
- A new module-level import that the docs build cannot install goes in
  `autodoc_mock_imports` in `docs/source/conf.py`.

## Known Issues

- **Pre-existing LSP type errors:** Pyright/Pylance reports false positives in
  `config/schema.py`, `builders/system_builder.py` and others, from missing
  OpenMM/OpenFF type stubs and Pydantic `default_factory` inference. They do
  not affect runtime.
