# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.3.0] - 2026-04-09 — Analysis Plugin System & OCP Compliance

Large-scale refactoring of the analysis and comparison subsystems to achieve
Open-Closed Principle compliance.  New analysis types can now be added by
dropping a package in `analyses/<name>/` with no modifications to core code.

### Breaking Changes

- **`compare/` public API removed.**  The lazy-export facade in
  `compare/__init__.py` (44 symbols) and the re-export layer in
  `compare/results/__init__.py` (8 classes) have been deleted.  All consumers
  must import from concrete submodules (`compare.statistics`, `compare.config`,
  `compare.io`, `compare.core.base`).  No external downstream consumers exist.
- **Three dead registries removed.**  `AnalyzerRegistry`,
  `AnalysisSettingsRegistry`, and `ComparisonSettingsRegistry` (plus all 16
  `@register` decorators) were deleted from `compare/registries.py` during
  the initial OCP refactor.
- **`compare/plotters/` directory removed.**  All per-analysis plotting code has
  been inlined into each plugin's package.  `compare/plotter.py` (plot-all
  orchestrator) has been rewritten to delegate to plugins.
- **`compare/formatters.py` removed.**  RMSF-specific formatters moved to
  `analyses/rmsf/`.
- **`compare/comparators/` directory removed.**  All 11 old comparator files
  deleted; comparison logic lives in each plugin's `compare()` method.
- **Contacts plugin slimmed.**  `_binding_preference.py` (3,938 lines),
  `_surface_exposure.py` (369 lines), and `_helpers.py` (369 lines) extracted
  to `analyses/shared/` as reusable cross-plugin compute utilities.

### Added

- **`polyzymd new-analysis <name>` scaffold CLI command.**  Generates a complete
  plugin package (`__init__.py`, `_comparison_results.py`) and test file
  (`test_<name>_plugin.py`) with working discovery, compute, aggregate, and
  extract_metrics implementations.  Supports `--class-name`, `--force`,
  `--dry-run`, and `--project-root` options.  Validates plugin names
  (snake_case, no collisions) and class names (PascalCase, valid identifiers).
  (`cli/scaffold.py`, `cli/main.py`)
- **Analysis plugin framework.**  `Analysis` ABC (`analyses/base.py`),
  `pkgutil`-based auto-discovery (`analyses/discovery.py`), lifecycle
  orchestrator (`analyses/orchestrator.py`), and default scalar comparison
  pipeline (`analyses/stats.py`).  Eight analysis types migrated as plugins:
  RMSF, contacts, distances, catalytic triad, secondary structure, exposure,
  binding free energy, polymer affinity.
- **Shared analysis utilities.**  `analyses/shared/` package with reusable
  `TrajectoryLoader`, alignment helpers, autocorrelation, plot utilities,
  `surface_exposure.py`, `binding_preference.py`, and
  `binding_preference_helpers.py`.
- **Framework hardening.**  `compute_replicate()` return values are validated
  (None rejected).  `ComparisonContext` gains `failed_conditions` and
  `aggregated_results` fields.  Orchestrator warns on undeclared plugin
  dependencies.
- **`cohens_d()` default corrected.**  `rmsf_mode` default changed from `True`
  to `False` in `compare/statistics.py` so the general-purpose function is not
  biased toward one analysis type.
- **Unified `run` command.**  `polyzymd run --engine gromacs|openmm` replaces the
  old `run-gromacs` command.  OpenMM engine runs simulations locally using the
  existing runner infrastructure.  (`cli/main.py`)
- **`submit --generate-only` flag.**  Generates SLURM scripts without submitting,
  replacing the previous `--dry-run` behavior for script generation.
  (`cli/main.py`, `workflow/daisy_chain.py`)
- **RMSD convergence detection.**  Always-on sliding-window slope diagnostic
  integrated into the RMSD analysis plugin.  Detects when RMSD timeseries have
  plateaued.  Configurable via `convergence_window_size_ns`,
  `convergence_step_size_ns`, `convergence_slope_threshold`, and
  `convergence_sustained_for_ns` settings.  (`analyses/rmsd/`)
- **Shared convergence utility.**  `analyses/shared/convergence.py` provides
  `find_convergence_time()` for any timeseries — reusable across plugins.
  (`analyses/shared/convergence.py`)
- **Convergence diagnostic plots.**  Dual-axis panels showing RMSD timeseries
  with sliding slope and convergence markers.  Controlled via
  `show_convergence_plots` plot setting.  (`analyses/rmsd/_plotters.py`)
- **"Establishing Convergence" documentation.**  New Explanation page covering
  convergence concepts, the sliding-window algorithm, parameter tuning, and
  limitations.  (`docs/source/explanation/convergence_detection.md`)

### Changed

- **All analysis plugins are now packages.**  Each plugin lives in
  `analyses/<name>/` with `__init__.py` (plugin class), `_comparison_results.py`
  (result models), and optional `_plotting.py`, `_formatters.py` modules.
  Single-file plugins no longer exist.
- **Plugin-specific code co-located.**  Result models, formatters, plotters, and
  settings that were previously scattered across `compare/results/`,
  `compare/formatters.py`, `compare/plotters/`, and `compare/comparators/` now
  live inside each plugin's package.
- **`compare/` is now shared infrastructure only.**  Contains statistics,
  ComparisonConfig/PlotSettings/PlotTheme, IO helpers, CLI commands, and
  base result classes — no analysis-specific code.
- **`submit --dry-run` is now preview-only.**  Validates configuration and
  previews the submission plan without writing any files.  Use `--generate-only`
  for script generation.  (`cli/main.py`, `workflow/daisy_chain.py`)
- **RMSD cache identity.**  Per-replicate and aggregated cache filenames now
  include a settings fingerprint, preventing stale results when analysis
  parameters change.  (`analyses/rmsd/`)
- **RMSD `overall_median` computation.**  Fixed to use `np.median()` instead of
  `np.mean()` for correctness.  (`analyses/rmsd/`)
- **RMSD comparison statistics.**  t-test and ANOVA are now guarded for
  insufficient replicates (n < 2), reporting "not testable" instead of NaN
  statistics.  (`analyses/rmsd/`)
- **RMSD comparison resilience.**  Comparison now gracefully handles missing
  run-condition combinations instead of raising `KeyError`.  (`analyses/rmsd/`)
- **Convergence input validation.**  `find_convergence_time()` rejects NaN/inf
  inputs with clear `ValueError` messages.  (`analyses/shared/convergence.py`)
- **Job naming consistency.**  `create_job_name()` now uses rounding (matching
  directory naming) instead of truncation for polymer percentages.
  (`workflow/slurm.py`)
- **`SubmissionResult` state model.**  Added `is_generated_only` field to
  distinguish generate-only from dry-run submissions.
  (`workflow/daisy_chain.py`)
- **CLI exception handling.**  `run` and `submit` commands now catch specific
  operational exceptions instead of broad `except Exception`.  (`cli/main.py`)

### Removed

- **`run-gromacs` command.**  Replaced by `polyzymd run --engine gromacs`.  All
  references updated across CLI, docs, and configuration.  (`cli/main.py`)

### Fixed

- **`"default"` reaction template paths no longer rejected during config
  loading.**  When a user specified `initiation: "default"` (or
  `polymerization` / `termination`) in the YAML config, `_expand_paths()`
  unconditionally treated the string as a relative filesystem path and
  prepended the config directory, producing e.g.
  `/ocean/projects/.../default`.  The Pydantic field validator in
  `ReactionConfig` then rejected it because the path lacked a `.rxn`
  extension.  The loader now recognises `"default"` as a sentinel value
  and passes it through untouched so the validator can resolve it to the
  bundled ATRP reaction templates.  (`config/loader.py`)
- **Dead code removed.**  `catalytic_triad/_formatters.py` (427 lines, zero
  imports), `rmsf/_comparison_results_legacy.py`, ghost `compare/plotters/`
  directory, stale `__pycache__/` files.
- **`validate_name()` exception handling narrowed.**  `except Exception:` changed
  to `except ImportError:` in the scaffold's collision-check path.
- **NPZ sidecar lookup.**  RMSD plotters now load NPZ paths from result metadata
  instead of reconstructing filenames, preventing stale sidecar binding.
  (`analyses/rmsd/_plotters.py`)
- **File descriptor leaks.**  NPZ file loading in RMSD plotters now uses context
  managers.  (`analyses/rmsd/_plotters.py`)

### Documentation

- Updated `docs/source/tutorials/extending_analyses.md` with scaffold command
  as the primary entry point for contributors.
- Updated `AGENTS.md` and `.opencode/instructions/` for refactored layout.
- Updated `docs/source/api/compare.md` to remove stale module references.

### Tests

- 940 tests passing (6 skipped), up from 908 at branch start.
- Added 44 scaffold tests (name validation, class-name validation, file
  generation, code quality, CliRunner integration).
- Removed 12 obsolete registry and smoke tests.

### Known Limitations

- **Contacts cache hash does not include analysis parameters such as cutoff.**
  Cached contact results may be reused across runs with different cutoff values
  (e.g., 4.0A vs 4.5A). Until cache-key hashing is expanded, clear the contacts
  analysis cache when changing contact criteria.

## [1.2.1] - 2026-04-01

### Fixed

- **`resid` restraint selections now match PDB residue numbers correctly.**
  The restraint atom selection parser (`_parse_selection`) used OpenMM's
  sequential 0-based `Residue.index` and assumed a simple `resid - 1`
  conversion, which only worked when PDB numbering started at 1 with no gaps.
  For structures with non-standard starting residue numbers (e.g., PDB 4TGL
  starting at residue 5), the offset was wrong and selections like
  `resid 144 and name OG` would silently target the wrong residue or fail
  with "No atoms match selection". Now uses `Residue.id` (the PDB `resSeq`
  number) for direct matching without any offset arithmetic.
  (`core/restraints.py`)

## [1.2.0] - 2026-03-24

### Added

- **Analysis workflow and CLI.** Added the `polyzymd analyze` command family with
  YAML-driven setup and execution for RMSF, contacts, distances, catalytic triad,
  and secondary-structure analyses, plus shared loading, alignment, PBC handling,
  result models, and aggregation helpers. (`src/polyzymd/analyses/`)
- **Comparison engine for multi-condition studies.** Added registry-based
  comparators, typed comparison result models, shared statistical utilities, and
  a generic `polyzymd compare run` workflow for RMSF, contacts, distances,
  catalytic triad, exposure dynamics, binding free energy, polymer affinity,
  and secondary structure. (`src/polyzymd/compare/`)
- **Config-driven plotting stack.** Added `polyzymd compare plot-all`,
  registry-based plot discovery, shared plot themes, and publication-oriented
  plotters for RMSF, contacts, distances, catalytic triad, secondary structure,
  binding free energy, exposure, and polymer affinity. (`src/polyzymd/compare/plotter.py`,
  `src/polyzymd/compare/plotters/`)
- **Secondary-structure comparison support.** Added DSSP-backed secondary
  structure analysis, comparison results, and plotting so secondary structure is
  part of the stable release analysis stack. (`src/polyzymd/analyses/secondary_structure.py`,
  `src/polyzymd/compare/comparators/secondary_structure.py`,
  `src/polyzymd/compare/plotters/secondary_structure.py`)
- **Comprehensive analysis documentation.** Added end-to-end tutorials,
  cookbook-style guides, API pages, and extension guides covering analysis,
  comparison, and plotting workflows. (`docs/source/tutorials/`,
  `docs/source/api/compare.md`)

### Changed

- **Release presentation labeling for debated metrics.** Binding preference,
  exposure dynamics, binding free energy, and polymer affinity remain available
  from the CLI and plotting pipeline, but PolyzyMD now marks them explicitly as
  experimental in command output, plot listings, generated text reports, figure
  annotations, config templates, and user-facing docs. (`src/polyzymd/core/experimental.py`,
  `src/polyzymd/compare/cli.py`, `src/polyzymd/compare/plotter.py`, `README.md`)
- **Stable release scope for analysis demos.** The presentation-ready stable
  comparison stack is now RMSF, contacts, distances, catalytic triad, and
  secondary structure, while the debated science-facing metrics remain visible
  but clearly labeled as experimental. (`README.md`,
  `docs/source/tutorials/analysis_compare_conditions.md`)

### Fixed

- **Comparison result and plotting reliability.** Fixed multiple comparison and
  plotting issues uncovered while building the release branch, including cached
  result discovery, condition-specific result paths, partition-aware BFE plots,
  shared-path bugs in the plot orchestrator, contacts/distances comparison edge
  cases, and corrupted-trajectory handling in contacts/exposure workflows.
   (`src/polyzymd/compare/`, `src/polyzymd/analyses/contacts.py`,
    `src/polyzymd/analyses/distances.py`)

### Known Limitations

- **Analysis supports OpenMM trajectories only.** The `polyzymd analyze`
  commands expect DCD trajectories in PolyzyMD's standard directory layout.
  GROMACS XTC trajectory support is planned for v1.2.1 (#47). Users running
  GROMACS simulations should use native GROMACS analysis tools or MDAnalysis
  directly until then.

## [1.1.1] - 2026-03-23

### Fixed

- **Applied `ruff format` to 5 source files** that failed the CI formatting check
  (`cli/main.py`, `config/loader.py`, `config/schema.py`, `simulation/runner.py`,
  `workflow/slurm.py`).

## [1.1.0] - 2026-03-23

### Breaking Changes

- **Pixi replaces conda/mamba for environment management.** SLURM job scripts
  now use `pixi shell-hook` instead of `module load` + `conda activate`.
  Existing conda environments still work for local use, but HPC job submission
  requires a pixi installation. See the updated
  [Installation Guide](https://polyzymd.readthedocs.io/en/latest/tutorials/installation.html).
- **Removed `polyzymd-submit` and `polyzymd-continue` entry points.** These
  console-script aliases were broken and unused. Use `polyzymd submit` and
  `polyzymd run-segment` instead.
- **Removed deprecated GROMACS exporter API.** `PositionRestraintGenerator.generate()`,
  `generate_all_from_config()`, and the `TopologyModifier` class have been
  removed. These were dead code never used externally. Use
  `PositionRestraintGenerator.add_posres_to_itp_files()` instead.

### Added

- **`polyzymd status` CLI command.** Displays a compact progress overview
  for all replicates of a simulation with colored Unicode progress bars,
  completion percentages, nanosecond progress, and per-replicate status.
  Auto-detects replicate directories via the naming template. Read-only
  (uses `load_progress()` only). (`cli/main.py`, `cli/colors.py`,
  `config/schema.py`)
- **Wall-time restart checkpoints for SLURM preemption resilience.** The
  simulation loop now saves portable `restart_state.xml` +
  `restart_system.xml` at a configurable wall-time interval
  (`checkpoint_interval` in config, default 60s). On SLURM preemption,
  the loop detects SIGTERM within ~15s (via adaptive sub-chunking) and
  saves an interrupted state before the grace period expires. Previously,
  a single `simulation.step(200000)` call could block for ~2 minutes,
  leaving no time for graceful shutdown within a 120s grace period.
  (`simulation/runner.py`, `simulation/continuation.py`,
  `simulation/signals.py`, `config/schema.py`)
- **Adaptive sub-chunk sizing.** After the first checkpoint interval, the
  loop measures actual steps/second and adjusts the sub-chunk size to
  target ~15s between interrupt checks. This ensures responsive signal
  handling regardless of system size or hardware speed.
  (`simulation/runner.py`, `simulation/continuation.py`)
- **Portable recovery path priority.** Continuation recovery now prefers
  portable state XML files (`interrupted_state.xml`, `restart_state.xml`)
  over binary `.chk` checkpoints, which are not portable across
  heterogeneous GPU clusters. Binary `.chk` is only used as a last-resort
  fallback for legacy interrupted segments or hard-killed segments.
  (`simulation/continuation.py`)
- **Per-module colored logging.** Each module group (builders, simulation,
  workflow, exporters, etc.) gets a distinct near-white tinted color for
  INFO/DEBUG log messages, making it easy to visually distinguish which
  subsystem produced each log line. WARNING stays amber yellow, ERROR stays
  red. Colors auto-detect terminal capability (truecolor > 256-color > basic
  > none) and respect the `NO_COLOR` environment variable.
  (`cli/colors.py` — new module)
- **`--no-color` CLI flag.** Disables all ANSI color output for logging and
  `colored_echo` messages. Added to the top-level `polyzymd` command group.
  (`cli/main.py`)
- **Hard-kill recovery.** If a SLURM job is killed without a graceful signal
  (e.g., node failure, `scancel`, OOM), the next job detects the incomplete
  segment via checkpoint + CSV analysis and resumes from the last checkpoint
  rather than re-running the entire segment. (`simulation/runner.py`,
  `simulation/progress.py`)
- **Interruptible equilibration.** Equilibration stages now run in chunked
  steps and respond to SIGUSR1/SIGTERM, allowing graceful shutdown during
  long equilibration phases. Previously, interruption during equilibration
  meant the entire stage had to be re-run. (`simulation/runner.py`)
- **Checkpoint-based continuation.** `run-segment` automatically determines
  whether to build, continue from checkpoint, or skip based on filesystem
  state. No manual segment tracking required.
- **FAILED segment cleanup.** `run-segment` detects and removes incomplete
  `FAILED`-state segments before retrying, preventing permanent stuck states.
  (`cli/main.py`)
- **`CheckpointReporter` for production segments.** Segment 0 now saves
  `system.xml` early and uses a dedicated checkpoint reporter, enabling
  recovery even if the simulation crashes before the first trajectory frame.
  (`simulation/runner.py`)
- **`--pixi-env` option** for `polyzymd submit` and `polyzymd recover`.
  Overrides the default pixi environment name in generated SLURM scripts
  (default is auto-selected based on the SLURM preset). (`cli/main.py`)
- **`--memory` option for `polyzymd recover`.**  Overrides the SLURM memory
  allocation in recovery job scripts, matching the existing `--memory` flag
  on `polyzymd submit`. Useful when a job OOM-killed and needs to be resumed
  with more RAM. (`cli/main.py`)
- **squeue-based duplicate detection at submission time.** Both `polyzymd
  submit` and `polyzymd recover --submit` now query `squeue` for
  RUNNING/PENDING jobs with the same job name before submitting.  If a
  duplicate is found, submission is blocked with a clear error message.
  The check is best-effort: if `squeue` is unavailable (non-SLURM
  environment, CI), a warning is logged and submission proceeds normally.
  A new `--force` flag on both commands allows explicit override.
  (`workflow/daisy_chain.py`, `cli/main.py`)
- **`_estimate_steps_from_csv` helper.** Estimates completed steps from
  `state_data.csv` when progress.json is missing or stale, enabling accurate
  progress reporting after hard kills. (`simulation/progress.py`)

### Changed

- All `click.echo()` calls in the CLI migrated to `colored_echo()` with
  phase-aware coloring (e.g., build commands use sage green, workflow
  commands use lavender). Success messages (`click.style(fg="green")`) and
  error messages (`click.style(fg="red")`) are preserved as-is.
- All `print()` calls in production code (`workflow/daisy_chain.py`,
  `exporters/gromacs.py`, `data/solvents/_generator.py`) migrated to
  `LOGGER.info()` so they flow through the colored logging formatter.
- SLURM job scripts now activate the environment via
  `pixi shell-hook -e <env> --manifest-path <path>` instead of
  `module load` + `conda activate`. The manifest path is auto-detected
  from the `polyzymd` binary location at submission time.
- `pixi.toml` trimmed to actual runtime dependencies with three environments:
  `build` (no CUDA), `cuda-12-4` (CU Boulder Blanca), `cuda-12-6`
  (PSC Bridges2).
- Added `openbabel` to `pixi.toml` conda dependencies. Required at import
  time by `polymerist.polymers.building.mbconvert`, which `polyzymd build`
  triggers unconditionally.

### Fixed

- **Concurrency guard prevents duplicate segment execution.** When SLURM
  requeues a preempted job while a recovery script also resubmits, two jobs
  can race to start the next segment. `run-segment` now checks for any
  segment with a recently-modified checkpoint file (< 600s) classified as
  RUNNING and exits with code 2 (`EXIT_CODE_CONCURRENT`) instead of
  launching a concurrent segment. The SLURM bash wrapper intercepts exit
  code 2 and terminates cleanly without resubmitting, breaking infinite
  submit-cancel-resubmit loops that occurred when a job was accidentally
  double-submitted. (`cli/main.py`, `simulation/signals.py`,
  `workflow/slurm.py`)
- **Overall status now reflects the most recent segment, not `any()`.**
  When a simulation had mixed segment statuses (e.g., segment 0 INTERRUPTED,
  segment 1 FAILED — as seen with the CALB replicate 2 infinite-loop bug),
  the `any(INTERRUPTED)` check in the status cascade fired before
  `any(FAILED)`, setting `progress.status` to `"interrupted"`. This misled
  the user into thinking auto-resume would handle recovery when the
  simulation actually needed manual resubmission. The status derivation
  logic now uses the highest-index segment's status to determine the
  overall state: if the latest segment is FAILED, overall = FAILED. A new
  `_derive_overall_status()` helper centralises this logic (previously
  duplicated in `scan_filesystem` and `validate_progress`). Additionally,
  cleanup blocks in `run-segment` (FAILED segment removal, hard-kill
  cleanup) now recompute `progress.status` before saving, preventing stale
  status values from persisting in `progress.json`.
  (`simulation/progress.py`, `cli/main.py`)
- **Hard-killed segments now retry in-place instead of advancing.**
  When SLURM kills a job without a grace period (SIGKILL, node failure,
  OOM), no `INTERRUPTED` marker file is written. Previously, the next job
  would classify the segment as INTERRUPTED via the stale-checkpoint
  heuristic and advance to a new segment index, loading from
  `restart_state.xml` — potentially losing all work done after that
  periodic checkpoint. Now, `run-segment` detects this case (highest-index
  segment classified INTERRUPTED but missing the `INTERRUPTED` marker file),
  cleans up the incomplete directory, and removes it from progress. This
  causes `get_next_segment_info()` to reassign the same index, retrying
  from the previous completed segment's state with no data loss.
  (`cli/main.py`)
- **`_estimate_steps_from_csv` now returns per-segment step counts.**
  Previously, the function returned the raw cumulative step number from the
  last CSV row. OpenMM's `StateDataReporter` writes cumulative integrator
  steps (from time=0 including equilibration and all prior segments), so for
  continuation segments this massively overcounted progress — causing
  `validate_progress()` to mark in-progress simulations as completed.  The
  function now computes `last_step - first_step` for the correct per-segment
  delta.  For single-row CSVs it returns 0 (safe undercount).
  (`simulation/progress.py`)
- **Equilibration `finished_at` timestamps are now populated.**
  `EquilibrationStageRecord.finished_at` existed in the Pydantic model but
  was always null. `scan_equilibration_stages()` now sets it from the
  checkpoint file's mtime, and `_run_initial_segment()` sets it to the
  current time during live runs. (`simulation/progress.py`, `cli/main.py`)
- **`polyzymd status` ns calculation now includes interrupted segments.**
  Previously, the status command used `time_completed_ns()` which only
  counts COMPLETED segments, showing 0.000 ns even when millions of steps
  had been simulated across interrupted replicates. Now uses
  `total_steps_completed * timestep_fs / 1e6` for accurate progress
  display. (`cli/main.py`)
- **`polyzymd status` summary no longer falsely reports "All completed".**
  The summary line only counted `interrupted`, `failed`, `not_started`,
  and `not_found` statuses as needing attention, so replicates with
  `running` status (including stale-running jobs that were killed without
  graceful shutdown) fell through to the "All N replicates completed!"
  message. The summary now tracks `completed_count`, `running_count`, and
  `need_attention` separately, and only shows the green completion message
  when every replicate has `status == completed`. (`cli/main.py`)
- **`polyzymd status` now detects stale "running" replicates.** Switched
  from `load_progress()` (raw JSON read) to `load_or_scan_progress()`
  which validates progress against the filesystem. If a checkpoint file
  is older than 10 minutes, the segment is reclassified from `running` to
  `interrupted`, matching what `polyzymd recover` already sees. The
  corrected status is saved back to `progress.json`. (`cli/main.py`)
- **Position restraints now applied to all polymer ITP files.** Previously,
  only the first polymer ITP (`_MOL1.itp`) received `#ifdef POSRES_POLYMER`
  blocks. With random copolymers, OpenFF Interchange generates a separate
  molecule type (and ITP file) per unique polymer sequence, leaving most
  polymer chains unrestrained. The rewritten `PositionRestraintGenerator`
  discovers all polymer ITPs, parses each one's `[ atoms ]` section to
  identify heavy atoms by atom name (HMR-safe), and appends position
  restraint blocks to every polymer ITP. (`exporters/gromacs.py`)
- **Residue numbering in .gro files is now globally sequential.** OpenFF
  Interchange's GRO writer computes `(residue_index + copy_index) % 100_000`,
  which creates a sliding +1 offset for multi-residue molecules (polymers).
  A new post-processing step (`_fix_gro_residue_numbering`) assigns globally
  sequential residue numbers across all multi-residue molecule copies, enabling
  unique residue-based selection in MDAnalysis (e.g., `resid 11:15` for the
  third polymer chain). Single-residue molecules (water, ions) are left
  unchanged. (`exporters/gromacs.py`)
- Segment 0 progress loss: previously, a hard kill during segment 0 could
  lose all progress because no checkpoint existed. Now `system.xml` is saved
  at the start of production, enabling checkpoint-based recovery.
- Stale `.pyc` files from feature branches no longer cause import errors
  after branch switching (resolved by merging both feature branches).
- **Removed `espaloma-charge` from `pixi.toml`** to prevent a broken import
  chain. `polymerist` eagerly imports `espaloma_charge` at module level
  (in `_toolkits.py`), which pulls in `dgl`, which fails to load
  `libgraphbolt` when dgl and PyTorch versions are mismatched. Since
  polyzymd uses NAGL (not espaloma) for charge assignment, and NAGL >=0.2
  has a pure-PyTorch fallback that works without dgl, removing
  `espaloma-charge` eliminates the crash with no loss of functionality.
- **Fixed indentation bug in generated GROMACS run script.** The
  post-processing section had a misindented `echo` line that would cause
  the script to fail under `set -e`. (`exporters/gromacs.py`)
- **`recover --submit` no longer rebuilds system when equilibration is
  complete.**  When a replicate had completed equilibration but no production
  segments, `polyzymd recover --submit` generated a SLURM script that re-ran
  the full build routine.  Since polymer packing and solvation are
  non-deterministic, this produced a different atom count, causing
  `loadCheckpoint` to crash with `"wrong number of particles"`.  The
  `recover` command now detects pre-built system files
  (`solvated_system.pdb`, `system.xml`) and passes `--skip-build` to the
  generated script.  Additionally, `_run_initial_segment` now skips
  minimization and equilibration when `--skip-build` is active and
  equilibration is already recorded as complete in `progress.json`, jumping
  directly to production segment 0.
  (`cli/main.py`)
- **Co-solvent volume fraction validator no longer crashes with concentration-based
  co-solvents.** `validate_volume_fractions` called `sum()` over `volume_fraction`
  fields without filtering `None` values, raising `TypeError` when any co-solvent
  used `concentration` instead of `volume_fraction`. (`config/schema.py`)
- **Equilibration stages now honour `thermostat_timescale` for integrator friction.**
  The `thermostat_timescale` field was read from the stage config but never used;
   the integrator always received the default friction of 1.0/ps. Friction is now
  computed as `1.0 / thermostat_timescale`. (`simulation/runner.py`)
- **Barostat temperature now tracks the integrator during NPT temperature ramps.**
  When an equilibration stage used NPT ensemble with temperature ramping, the
  `MonteCarloBarostat` was initialized at the starting temperature but never
  updated as the ramp progressed. This caused the barostat to evaluate
  volume-move acceptance at the wrong temperature throughout the entire ramp,
  leading to incorrect pressure coupling. The ramp loop and final-temperature
  section now call `context.setParameter(MonteCarloBarostat.Temperature(), ...)`
  to keep the barostat in sync with the integrator. (`simulation/runner.py`)
- **EQ_INTERRUPTED marker now records the correct temperature during ramps.**
  The temperature ramp loop incremented `current_temp` *before* the interrupt
  check, so the `EQ_INTERRUPTED` marker saved a temperature one increment
  higher than what was actually simulated. The increment is now moved to
  *after* the interrupt check, and log messages report the correct
  temperature. (`simulation/runner.py`)
- **Temperature ramp resume no longer double-counts fast-forwarded chunks.**
  On resume, `current_temp` was initialized from `resume_temperature` (the
  value saved in the marker) and then the fast-forward skip loop *also*
  incremented `current_temp` for each skipped chunk, causing the simulation
  to jump ahead in temperature. The ramp loop now always starts from
  `stage.temperature_start` and lets the fast-forward loop reconstruct the
  correct temperature by skipping completed chunks. (`simulation/runner.py`)
- **Unrecoverable hard-kill state (checkpoint without system.xml) now raises
  immediately.** When a segment was hard-killed and only the periodic
  checkpoint existed (no `system.xml`), the continuation manager logged an
  error but fell through silently, returning paths to non-existent files.
  This caused confusing downstream `FileNotFoundError` messages.  Case 5b
  now raises `FileNotFoundError` immediately with a clear message.
  (`simulation/continuation.py`)
- **`check-progress` errors no longer trigger infinite SLURM resubmission.**
  Errors in `check-progress` (config load failure, missing progress file)
  exited with code 1 — the same code used for "work remains."  The SLURM
  bash wrapper interpreted any non-zero exit as "resubmit," causing an
  infinite loop on persistent errors.  Error conditions now exit with code 3
  (`EXIT_CODE_CHECK_ERROR`), and the SLURM template only resubmits on exit
  code 1. (`cli/main.py`, `simulation/signals.py`, `workflow/slurm.py`)
- **Progress file writes are now crash-safe with `fsync` before rename.**
  `save_progress()` already used atomic write-to-temp-then-rename, but did not
  call `fsync` on the temporary file before `os.replace()`.  On power loss or
  kernel panic the rename could be durable while the file contents were not,
  leaving a zero-length or corrupt `progress.json`.  The function now calls
  `f.flush()` and `os.fsync(f.fileno())` before the rename.
  (`simulation/progress.py`)
- **`SlurmConfig.from_preset()` now raises `ValueError` for unknown preset names.**
  Previously, an unrecognised preset name silently fell back to the `aa100` preset,
  masking typos in config files or CLI arguments.  The error message lists all valid
  presets. (`workflow/slurm.py`)
- **`save_config` no longer mutates the global `yaml.Dumper` representer registry.**
  The custom multiline-string representer was registered via `yaml.add_representer()`,
  which permanently alters `yaml.Dumper` for the entire process.  Now uses a local
  `Dumper` subclass so other YAML consumers are unaffected.  (`config/loader.py`)
- **`build --dry-run --gromacs` now shows the actual output path.**  The GROMACS
  dry-run summary printed the literal string `{projects_dir}/{replicate}/gromacs/`
  instead of interpolating the real directory.  (`cli/main.py`)
- **Reaction template paths (`initiation`, `polymerization`, `termination`) are now
  included in path resolution.**  `_expand_paths` and `_convert_paths_to_relative`
  only knew about `pdb_path`, `sdf_path`, `sdf_directory`, `cache_directory`, and
  `base_directory`.  Relative `.rxn` paths in the `reactions:` config block were
  passed through as-is, causing `FileNotFoundError` when the config file lived in
  a different   directory from the CWD.  (`config/loader.py`)
- **`to_signac_statepoint()` no longer crashes with concentration-based co-solvents.**
  The statepoint export unconditionally accessed `cosolvent.volume_fraction`, which
  is `None` for concentration-based co-solvents.  Now exports `_fraction` or
  `_molarity` depending on which is set.  (`config/schema.py`)
- **`load_checkpoint` now restores velocities, not just positions.**
  `getState()` was called with `getPositions=True` only, so
  `_current_velocities` remained stale (or `None`) after loading a
  checkpoint.  If equilibration stages subsequently checked
  `_current_velocities`, they could incorrectly re-randomize velocities
  instead of continuing from the checkpoint's kinetic state.
  (`simulation/runner.py`)
- **Sorted import block in `run-segment` handler.**  ruff I001 (import sort)
  violation in the equilibration progress save block.  (`cli/main.py`)
- **Signal handler no longer calls `LOGGER` (async-signal-unsafe).**  The
  `_handler()` function used `LOGGER.warning()`, which acquires Python's
  logging lock internally.  If the signal arrives while application code
  already holds that lock, the handler deadlocks.  Replaced with
  `os.write(2, ...)` which is async-signal-safe.  (`simulation/signals.py`)
- **Cross-check INTERRUPTED markers against CSV data to detect stale markers.**
  If a segment was gracefully interrupted, then restarted in-place and ran much
  further before being hard-killed, the old INTERRUPTED marker would persist
  with the original (too-low) step count while the CSV reflected all the work
  actually done.  `_scan_segment_dir` now compares the marker's
  `steps_completed` against the CSV delta; if the CSV shows more than 2× the
  marker value and exceeds 1 million steps, the stale marker is overridden with
  the CSV estimate and a warning is logged.  This prevents undercounting
  completed steps, which would inflate the "remaining" calculation and cause the
  simulation to overshoot its target duration.
  (`simulation/progress.py`)

### Documentation

- Rewrote installation guide for pixi (replaces conda/mamba instructions).
- Removed phantom `polyzymd run` and `polyzymd continue` CLI references
  (these commands never existed in the code).
- Added `polyzymd run` section to CLI reference.
- Updated HPC guide with pixi shell-hook activation examples.
- Fixed stale `polymerist-env` environment name references.
- Updated troubleshooting guide for pixi workflow.

## [1.0.4] - 2025-01-15

- Initial public release on PyPI.
