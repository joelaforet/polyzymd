# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.1.0] - Unreleased

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

- **`polyzymd status` ns calculation now includes interrupted segments.**
  Previously, the status command used `time_completed_ns()` which only
  counts COMPLETED segments, showing 0.000 ns even when millions of steps
  had been simulated across interrupted replicates. Now uses
  `total_steps_completed * timestep_fs / 1e6` for accurate progress
  display. (`cli/main.py`)
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

### Documentation

- Rewrote installation guide for pixi (replaces conda/mamba instructions).
- Removed phantom `polyzymd run` and `polyzymd continue` CLI references
  (these commands never existed in the code).
- Added `polyzymd run-gromacs` section to CLI reference.
- Updated HPC guide with pixi shell-hook activation examples.
- Fixed stale `polymerist-env` environment name references.
- Updated troubleshooting guide for pixi workflow.

## [1.0.4] - 2025-01-15

- Initial public release on PyPI.
