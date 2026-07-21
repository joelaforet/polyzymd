# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

- Added strict mixed GLYCAM/OpenFF overlay diagnostics, ownership artifacts, and
  an opt-in `final_e2e` acceptance test gate for the compact periodic 1UBQ
  glycan/polymer system.

### Changed

- Canonical `moiety.force_field: glycam06` now selects GLYCAM; omitted moiety
  force fields inherit `force_field.small_molecule`, supported `.offxml` sources
  remain generic, and unknown labels fail without fallback.
- Mixed overlay now preserves baseline nonbonded global settings and force
  metadata while rejecting unsupported baseline custom forces that touch
  GLYCAM-owned atoms.

## [1.3.0] - 2026-04-09 — Analysis Plugin System & OCP Compliance

Large-scale refactoring of the analysis and comparison subsystems to achieve
Open-Closed Principle compliance.  New analysis types can now be added by
dropping a module or package in `analyses/` with no modifications to core code.

### Breaking Changes

- **Python 3.12 is now required.**  Python 3.10 and 3.11 support has been
  dropped for v1.3 and later.  Project metadata, CI, Read the Docs, legacy conda
  environment files, and the build, analysis, test, and CUDA 12.6 pixi
  environments now target Python 3.12 with NumPy 2. The `sim-cuda-12-4`
  environment intentionally uses NumPy 1.26 and OpenMM 8.1.x for CUDA 12.4
  compatibility.
- **Scalar Measurement API removed.**  The obsolete alternate scalar
  measurement abstraction (`MetricSpec`, `CacheIdentity`, `Measurement`,
  `ScalarMeasurement`, and `ScalarMeasurementAnalysis`) was deleted now that
  the MDAnalysis job/artifact lifecycle is the true contributor extension path.
  Catalytic triad now inherits directly from `Analysis` and exposes its primary
  comparison metric through `extract_metrics()`.
- **`compare/` package removed.**  The entire `compare/` package (lazy-export
  facade, re-export layer, registries, plotters, formatters, comparators) has
  been deleted.  Analysis-specific code now lives in each plugin's package under
  `analyses/<name>/`.  Shared statistics live in `analyses/stats.py` and
  `analyses/shared/inferential_statistics.py`.  No external downstream consumers
  existed.
- **Three dead registries removed.**  `AnalyzerRegistry`,
  `AnalysisSettingsRegistry`, and `ComparisonSettingsRegistry` (plus all 16
  `@register` decorators) were deleted from `compare/registries.py` during
  the initial OCP refactor.
- **`compare/plotters/` directory removed.**  All per-analysis plotting code
  now lives in each plugin's `_plotters.py` module.  The orchestrator
  (`analyses/orchestrator.py`) delegates plotting to each plugin's `plot()`
  method.
- **`compare/formatters.py` removed.**  RMSF-specific formatters moved to
  `analyses/rmsf/`.
- **`compare/comparators/` directory removed.**  All 11 old comparator files
  deleted; comparison logic lives in each plugin's `compare()` method.
- **Contacts plugin slimmed.**  Obsolete binding-preference and surface-exposure
  experiments were archived out of the active runtime instead of being exposed
  as reusable shared compute utilities.  Contacts-specific implementation
  details now remain private to the contacts plugin, while cross-plugin helpers
  live only in the documented shared utility modules.

### Added

- **`polyzymd new-analysis <name>` scaffold CLI command.**  Generates a complete
  single-file MDAnalysis-native plugin at `src/polyzymd/analyses/<name>.py` and
  test file (`tests/analyses/plugins/test_<name>.py`) with working discovery,
  compute, aggregate, extract_metrics, format, and plot implementations.
  Advanced/package scaffolds are available with `--advanced`, `--style dict`,
  or `--style pydantic` and create package files such as
  `src/polyzymd/analyses/<name>/__init__.py`.  Supports `--class-name`,
  `--force`, `--dry-run`, `--style`, and `--project-root` options.  Validates
  plugin names (snake_case, no collisions) and class names (PascalCase, valid
  identifiers).  (`cli/scaffold.py`, `cli/main.py`)
- **Analysis plugin framework.**  `Analysis` ABC (`analyses/base.py`),
  `pkgutil`-based auto-discovery (`analyses/discovery.py`), lifecycle
  orchestrator (`analyses/orchestrator.py`), and default scalar comparison
  pipeline (`analyses/stats.py`).  Active v1.3 plugins are contacts,
  distances, hydrogen bonds, RMSD, RMSF, radius of gyration, SASA, secondary
  structure, and catalytic triad.  Archived experiments such as exposure,
  binding free energy, and polymer affinity are not active runtime plugins.
- **Shared analysis utilities.**  `analyses/shared/` package with reusable
  `TrajectoryLoader`, alignment helpers, autocorrelation/statistics helpers,
  selection utilities, path helpers, SASA compatibility utilities, multi-run
  comparison helpers, convergence helpers, and plot utilities used by the
  active plugin set.
- **Framework hardening.**  `compute_replicate()` return values are validated
  (None rejected).  `ComparisonContext` gains `failed_conditions` and
  `aggregated_results` fields.  Orchestrator warns on undeclared plugin
  dependencies.
- **`cohens_d()` default corrected.**  `rmsf_mode` default changed from `True`
  to `False` in `analyses/shared/inferential_statistics.py` so the
  general-purpose function is not biased toward one analysis type.
- **Multi-engine simulation architecture.**  New `SimulationEngine` ABC
  (`engines/base.py`) and `match/case` dispatch (`engines/__init__.py`) allow
  multiple MD engines to coexist behind a shared interface.  Each engine
  provides trajectory layout resolution, SLURM script generation, submission,
  progress tracking, and recovery.  OpenMM is wrapped as the default engine
  (`engines/openmm/`).
- **GROMACS engine with GPU support.**  Full GROMACS simulation lifecycle —
  energy minimization, multi-stage equilibration with checkpoint resume,
  production MD, and trajectory post-processing — driven entirely from the
  existing `config.yaml`.  GPU acceleration is enabled via `gpu: true` in
  `GromacsEngineConfig`, which auto-composes `-nb gpu -pme gpu -update gpu`
  offload flags and requests GPU GRES in SLURM scripts.  Supports both
  thread-MPI (`gmx`) and real-MPI (`gmx_mpi`) binaries with automatic flag
  adjustment.  (`engines/gromacs/engine.py`, `engines/gromacs/slurm.py`,
  `engines/gromacs/binary.py`, `engines/gromacs/progress.py`)
- **GROMACS configuration model.**  `GromacsEngineConfig` Pydantic v2 model
  with 16 validated fields covering binary selection, GPU/MPI hardware,
  thread counts, mdrun flags (global, equilibration-only, production-only),
  environment exports, setup commands, and `command_prefix` for custom
  launchers.  Includes cross-field validators that warn on contradictory
  settings (e.g., `gpu: true` with `gmx_mpi`, explicit `ntmpi` conflicting
  with `mdrun_flags`).  (`config/schema.py`)
- **GROMACS SLURM script generator.**  Self-resubmitting Bash scripts with
  `SIGTERM` trap for preemption resilience.  On receiving SIGTERM the script
  forwards the signal to `gmx mdrun`, waits for checkpoint flush, then
  resubmits itself via `sbatch`.  Generates EM, equilibration, and production
  stages with checkpoint-aware resume.  Supports `env_exports`,
  `setup_commands`, `module_load`, and `-maxh` wall-time safety.
  (`engines/gromacs/slurm.py`)
- **`--constraint` CLI option.**  Maps to `#SBATCH --constraint` in generated
  scripts, supporting single values (`"A40"`), OR expressions (`"A40|A100"`),
  and AND expressions (`"avx2&rh8"`).  Available on both `submit` and
  `recover`.  (`cli/main.py`, `workflow/slurm.py`)
- **`--nodelist` CLI option.**  Pins jobs to specific SLURM nodes.  Validates
  bracket-style hostlists (e.g., `node[01-04]`).  (`cli/main.py`,
  `workflow/slurm.py`)
- **`--partition`, `--qos`, `--email` CLI overrides.**  Runtime overrides for
  SLURM partition, QoS, and email notifications on `submit` and `recover`
  commands.  (`cli/main.py`)
- **Engine-aware CLI commands.**  `submit`, `recover`, `status`, and
  `check-progress` all accept `--engine gromacs|openmm` and resolve
  engine-specific working directories, progress files, and recovery logic
  automatically.  (`cli/main.py`)
- **Unified `run` command.**  `polyzymd run --engine gromacs|openmm` replaces the
  old `run-gromacs` command.  OpenMM engine runs simulations locally using the
  existing runner infrastructure.  (`cli/main.py`)
- **`TrajectoryLayout` model.**  Engine-neutral description of topology and
  trajectory file locations, consumed by `TrajectoryLoader` for analysis.
  Each engine provides its own resolver.  (`engines/base.py`,
  `analyses/shared/loader.py`)
- **Shared `run_sbatch()` helper.**  Centralizes `sbatch` invocation with
  optional `module load` preloading (e.g., `slurm/blanca`), used by both
  `submit` and `recover`.  (`workflow/slurm_submit.py`)
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
- **Benjamini-Hochberg FDR correction.**  Pairwise t-tests are now corrected
  for multiple comparisons across the full family of tests.  Configurable via
  `fdr_alpha` in `comparison.yaml` defaults.  (`analyses/stats.py`,
  `analyses/shared/inferential_statistics.py`)
- **Configurable t-test method.**  `ttest_method` in `comparison.yaml` defaults
  selects Student's (`"student"`, default) or Welch's (`"welch"`) t-test for
  pairwise comparisons.  (`analyses/shared/defaults.py`, `analyses/stats.py`)
- **Tukey's HSD post-hoc testing.**  `posthoc_method` in `comparison.yaml`
  defaults selects either `"ttest_bh"` (default, pairwise t-tests with
  Benjamini-Hochberg correction) or `"tukey_hsd"` (Tukey's Honestly
  Significant Difference).  Tukey HSD uses a single studentized range
  distribution for FWER control without requiring separate FDR correction.
  (`analyses/shared/defaults.py`, `analyses/shared/inferential_statistics.py`,
  `analyses/stats.py`, `analyses/base.py`)
- **Post-hoc testing reference page.**  New documentation page at
  `docs/source/reference/posthoc_testing.md` covering both post-hoc methods,
  configuration, output fields, CLI significance markers, and edge cases.
  Cross-linked from comparison YAML reference, statistics best practices,
  comparison tutorial, and analysis comparison reference.
- **Structured exception hierarchy.**  `AnalysisError`, `ReplicateError`,
  `AggregationError`, `ComparisonError`, `PlotError`, and
  `PluginContractError` in `analyses/exceptions.py` replace generic
  `Exception` raises throughout the orchestrator.

### Changed

- **Scientific dependency stack updated for NumPy 2.**  Runtime metadata now
  requires `numpy>=2,<3` with Python 3.12-compatible lower bounds for SciPy,
  pandas, scikit-learn, matplotlib, seaborn, MDAnalysis, and MDTraj.  Pixi CUDA
  environments retain separate CUDA 12.4 and CUDA 12.6 targets. `sim-cuda-12-6`
  uses the current NumPy 2/OpenMM stack, while `sim-cuda-12-4` intentionally
  stays on NumPy 1.26 and OpenMM 8.1.x for CUDA 12.4 compatibility.
- **Bundled v1.3 analyses are package-organized.**  Active bundled plugins live
  in `analyses/<name>/` with `__init__.py` (plugin class) and optional private
  modules (`_plotters.py`, `_results.py`, `_comparison_results.py`,
  `_formatters.py`, `_aggregator.py`) as needed.  Contributor plugins and
  discovery still support both single-file modules and packages.
- **Plugin contract enforcement hardened.**  `compute_replicate()` returning
  `None` is now a hard `PluginContractError` (not a soft skip).
  `PluginContractError` propagates without being wrapped.  `compare()` and
  `plot()` return values are validated (dict/BaseModel/None for compare;
  `list[Path]` for plot).  `run_comparison()` fails fast on contract
  violations.  (`analyses/orchestrator.py`)
- **Discovery skip logic checks all path components.**  `_should_skip_module()`
  now checks ALL path components (not just the last one), so shared utility
  modules and private helper packages are skipped instead of being probed as
  plugins.  (`analyses/discovery.py`)
- **Settings fingerprint canonicalized.**  `settings_fingerprint()` now uses
  `json.dumps(settings.model_dump(mode="json"), sort_keys=True)` for
  deterministic hashing across Pydantic versions and Python dict ordering.
  (`analyses/shared/config_hash.py`)
- **`PlotContext.plot_settings` always materialised.**  `__post_init__`
  now creates a real `PlotSettings()` when `None` is passed, eliminating
  null-guard boilerplate in every plugin's `plot()` method.
  (`analyses/base.py`)
- **Plugin-specific code co-located.**  Result models, formatters, plotters, and
  settings that were previously scattered across `compare/results/`,
  `compare/formatters.py`, `compare/plotters/`, and `compare/comparators/` now
  live inside each plugin's package.
- **`compare/` package dissolved.**  Statistics moved to `analyses/stats.py` and
  `analyses/shared/inferential_statistics.py`.  ComparisonConfig/PlotSettings/
  PlotTheme moved to `config/comparison.py`.  CLI commands remain in `cli/`.
  Per-analysis result models, formatters, and plotters now live in each plugin's
  package.
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
- **`fdr_alpha` unified across both post-hoc methods.**  `fdr_alpha` is now
  used as the significance threshold for both `"ttest_bh"` (BH FDR) and
  `"tukey_hsd"` (FWER), and is also threaded to ANOVA significance.
  Significance markers and formatter footnotes dynamically use
  `result.fdr_alpha` instead of hardcoded `0.05`.  (`analyses/stats.py`)
- **`ttest_method` scoped to `ttest_bh` only.**  `ttest_method` is now
  explicitly documented and enforced as only relevant when `posthoc_method` is
  `"ttest_bh"`.  Tukey HSD ignores `ttest_method` entirely.
- **Tukey `p_value_adjusted` mirrors `p_value`.**  For Tukey HSD results,
  `p_value_adjusted` is set to the same value as `p_value` (Tukey p-values are
  already family-wise corrected) instead of being left `null`.
  (`analyses/stats.py`)
- **`PlotContext.plot_settings` typing strengthened.**  Field type changed from
  `PlotSettings | None = None` to `PlotSettings = field(default_factory=...)`.
  A `__post_init__` backstop materializes a default `PlotSettings()` if
  explicit `None` is passed, guaranteeing non-null for all plugin code.
  (`analyses/base.py`)
- **`PluginContractError` propagation hardened.**  `run_all_comparisons()` now
  has an explicit `except PluginContractError: raise` before the broader
  `except AnalysisError` catch, ensuring contract violations fail fast.
  (`analyses/orchestrator.py`)
- **`save_result()` errors reclassified.**  All 5 `save_result()` call sites in
  the orchestrator now wrap `OSError` as the appropriate lifecycle error
  (`ReplicateError`, `AggregationError`, or `ComparisonError`).
  (`analyses/orchestrator.py`)
- **Foreign condition filtering in orchestrator.**  `_prepare_conditions_with_filter()`
  now discards (not just warns about) foreign conditions returned by a plugin's
  `filter_conditions()` method.  (`analyses/orchestrator.py`)
- **Discovery import classification.**  `_OPTIONAL_HEAVY_DEPS` allowlist
  (openmm, MDAnalysis, parmed, etc.) governs skip-vs-reraise for
  `ImportError` during plugin discovery.  Known heavy deps are skip+warn;
  unknown `ImportError`s re-raise immediately.  (`analyses/discovery.py`)
- **Cache identity with settings fingerprint.**  Contact artifact path
  resolution now includes settings fingerprints where applicable.
  Fingerprinted files are searched first, legacy files as fallback.
  `find_enzyme_pdb()` sorts glob results and warns on ambiguity.
  (`analyses/contacts/_paths.py`)

### Removed

- **Welch ANOVA removed entirely.**  Only classical one-way ANOVA
  (`scipy.stats.f_oneway`) remains.  The `anova_method` field was removed
  from all models and configuration.  Welch's ANOVA was deemed unnecessary
  given that ANOVA is used only as an omnibus gate before post-hoc tests.

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
- **Scaffold test path corrected.**  `new-analysis` now generates test files at
  `tests/analyses/plugins/test_<name>.py` matching the project's test layout.
  (`cli/scaffold.py`)
- **Cache ambiguity detection.**  `contacts/_paths.py` now raises `ValueError`
  on ambiguous glob matches (>1 file) instead of silently picking the last
  alphabetical match.  (`analyses/contacts/_paths.py`)
- **GROMACS trajectory resolution order corrected.**  `resolve_trajectory_layout()`
  now prefers post-processed trajectories: `prod_centered.xtc` (whole molecules,
  centered protein) before `prod_nojump.xtc` before raw `prod.xtc`.
  Previously the raw trajectory could be selected first, breaking analyses that
  expect unwrapped coordinates.  (`engines/gromacs/engine.py`)
- **GROMACS topology resolution skips build artifacts.**  Generic `*.pdb` glob
  removed from topology search; resolver now uses explicit candidates
  (`solvated_system.pdb`, `<prefix>.pdb`, `<prefix>.gro`) to avoid selecting
  `_PACKING_MOLECULE*.pdb` build artifacts left by Packmol.
  (`engines/gromacs/engine.py`)
- **TrajectoryLoader resolves engine subdirectory.**  `_resolve_layout()` now
  calls `engine.resolve_engine_working_directory()` before layout resolution,
  so GROMACS files in the `gromacs/` subdirectory are found correctly.
  Previously the loader passed the run-level directory directly, missing the
  engine-specific subdirectory.  (`analyses/shared/loader.py`)
- **Unused imports removed.**  Cleaned up stale imports across `cli/main.py`,
  `workflow/daisy_chain.py`, and several analysis plugins.
- **`Optional[X]` → `X | None` migration.**  Updated type annotations in
  `cli/main.py` and `workflow/daisy_chain.py` to use Python 3.10+ union syntax.
- **Tukey HSD graceful degradation on n<2.**  `tukey_hsd()` now returns an
  empty list instead of raising `ValueError` when fewer than 2 groups or
  fewer than 2 observations per group are provided, matching the NaN-return
  pattern used by other statistical functions.
  (`analyses/shared/inferential_statistics.py`)
- **Cohen's d direction for d==0.0.**  Now returns `direction="unchanged"`
  instead of `"lower"` when the effect size is exactly zero.
  (`analyses/shared/inferential_statistics.py`)
- **ANOVA <2 groups guard.**  Explicit guard returns NaN result when fewer
  than 2 groups are provided, instead of letting SciPy raise.
  (`analyses/shared/inferential_statistics.py`)
- **`ConditionSummary.n_replicates` validated across metrics.**  Now uses
  the minimum replicate count across all metrics and warns on mismatch.
  (`analyses/stats.py`)
- **Control label log messages corrected.**  Orchestrator now uses
  `original_control` in log messages before reassignment to prevent
  misleading "not found" messages that name the wrong label.
  (`analyses/orchestrator.py`)
- **`fdr_alpha` validation for all post-hoc methods.**  `pairwise_comparisons()`
  validates `fdr_alpha` is in `(0, 1]` and not `NaN` for both `"ttest_bh"`
  and `"tukey_hsd"`, raising `ValueError` on invalid input.
  (`analyses/stats.py`)

### Documentation

- Corrected `posthoc_testing.md` — Tukey HSD significance now documented as
  using `fdr_alpha` (not hardcoded 0.05); `p_value_adjusted` accurately
  described as mirroring `p_value` for Tukey results; ANOVA significance
  documented as using `fdr_alpha`; Tukey n<2 edge case documented as graceful
  empty-result return.
- Corrected `comparison_yaml.md` — `fdr_alpha` description updated to reflect
  unified behavior across both post-hoc methods.
- Corrected `analysis_comparison_reference.md` — `rmsd` and `rg` correctly
  listed as custom compare (not default); `fdr_alpha` description updated for
  both post-hoc methods.
- Corrected `analysis_statistics_best_practices.md` — added Tukey HSD
  subsection alongside BH; `fdr_alpha` described as configurable at both
  `defaults:` and per-plugin level; multiple comparisons table updated to
  reference both correction methods.

- Added `docs/source/reference/posthoc_testing.md` — post-hoc testing methods
  reference page covering t-test+BH vs Tukey HSD, configuration, output fields,
  CLI significance markers, and edge cases.
- Added `posthoc_method` and `ttest_method` fields to `comparison_yaml.md`
  defaults reference table.
- Cross-linked post-hoc testing page from `analysis_comparison_reference.md`,
  `analysis_statistics_best_practices.md`, and `analysis_compare_conditions.md`.
- Updated `docs/source/contributor_guide/extending_analyses.md` with scaffold
  command as the primary entry point for contributors.
- Updated `AGENTS.md` and `.opencode/instructions/` for refactored layout.
- Removed stale `api/compare.md` API page (module no longer exists).
- Added `docs/source/how_to/gromacs_export.md` — comprehensive GROMACS HPC
  how-to covering GPU/CPU configuration, SLURM submission, preemption
  resilience, constraint-based GPU targeting, recovery, and monitoring.
  (~700 lines with copy-paste recipes for common HPC workflows.)
- Added GROMACS engine config and CLI options to `cli_reference.md` and
  `configuration.md` reference pages.
- Updated `quickstart.md` with GROMACS-tabbed submission and monitoring steps.
- Added `scancel --signal=KILL` warning to both `gromacs_export.md` and
  `hpc_slurm.md` for stopping preemption-resilient jobs permanently.

### Tests

- Expanded regression coverage substantially from the branch baseline.
- Added GROMACS engine tests: binary resolution (`test_gromacs_binary.py`),
  engine adapter (`test_gromacs_engine.py`), SLURM script generation
  (`test_gromacs_slurm.py`), progress tracking (`test_gromacs_progress.py`),
  trajectory layout (`test_gromacs_layout.py`), and engine dispatch
  (`test_dispatch.py`, `test_base.py`), with expanded engine/workflow/cli/config
  coverage.  (`tests/engines/`)
- Added GROMACS smoke tests for the active analysis plugin set, verifying that
  plugins resolve trajectory layouts from the GROMACS engine correctly.
  (`tests/analyses/plugins/test_*_gromacs_smoke.py`)
- Added `run_sbatch()` tests (`test_slurm_submit.py`) and engine-aware CLI
  tests for submit, recover, status, and check-progress.
- Added scaffold tests (name validation, class-name validation, file
  generation, code quality, CliRunner integration).
- Added Tukey HSD tests (basic operation, single-condition edge case,
  insufficient-replicate edge case, format integration).
- Added plugin contract enforcement tests (`TestContractEnforcement` in
  `test_orchestrator.py`).
- Added discovery robustness tests (shared descendant skip, getattr logging).
- Added contacts cache identity/path regression tests and archived setting
  rejection coverage.
- Added settings fingerprint canonicalization tests (`test_config_hash.py`).
- Added Phase 14 test hardening: exact Tukey pair assertions, BH/Tukey boundary
  tests, ANOVA alpha threading tests, `fdr_alpha` validation tests, cache
  precedence tests, fail-fast verification, and import classification tests.
  (`tests/analyses/shared/`, `tests/analyses/test_stats.py`,
  `tests/analyses/test_orchestrator.py`, `tests/analyses/test_discovery.py`,
  `tests/analyses/plugins/test_contacts.py`)
- Full test suite reorganized to mirror source layout (`tests/analyses/plugins/`,
  `tests/analyses/shared/`, `tests/cli/`, etc.).
- Removed obsolete registry and smoke tests.

### Known Limitations

- **Config hash mismatch warning prints 66+ times.**  The warning should print
  once per analysis run but currently fires per-frame.  Does not affect results.

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
  catalytic triad, secondary structure, and historical experimental metrics
  that were later archived from the v1.3 active runtime. (`src/polyzymd/compare/`)
- **Config-driven plotting stack.** Added `polyzymd compare plot-all`,
  registry-based plot discovery, shared plot themes, and publication-oriented
  plotters for RMSF, contacts, distances, catalytic triad, secondary structure,
  and historical experimental metrics that were later archived from the v1.3
  active runtime. (`src/polyzymd/compare/plotter.py`,
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
  exposure dynamics, binding free energy, and polymer affinity were historical
  v1.2 experimental metrics; they were later archived from the v1.3 active
  runtime. At the time, PolyzyMD marked them explicitly as experimental in
  command output, plot listings, generated text reports, figure annotations,
  config templates, and user-facing docs. (`src/polyzymd/core/experimental.py`,
  `src/polyzymd/compare/cli.py`, `src/polyzymd/compare/plotter.py`, `README.md`)
- **Stable release scope for analysis demos.** The presentation-ready stable
  comparison stack is now RMSF, contacts, distances, catalytic triad, and
  secondary structure, while the debated science-facing metrics remain visible
  but clearly labeled as experimental. (`README.md`,
  `docs/source/tutorials/analysis_compare_conditions.md`)

### Fixed

- **Comparison result and plotting reliability.** Fixed multiple comparison and
  plotting issues uncovered while building the release branch, including cached
  result discovery, condition-specific result paths, partition-aware historical
  BFE plots, shared-path bugs in the plot orchestrator, contacts/distances
  comparison edge cases, and corrupted-trajectory handling in historical
  contacts/exposure workflows.
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
- `pixi.toml` trimmed to actual runtime dependencies with separate environments:
  `build` (no CUDA), `sim-cuda-12-4` (CU Boulder Blanca), and `sim-cuda-12-6`
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
