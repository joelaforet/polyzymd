# Scaffold your first analysis plugin

This tutorial gives you one safe first success with the PolyzyMD analysis
scaffold. You will generate one default single-file plugin named
`solvent_shell`, inspect the files the scaffold creates, and run the generated
tests.

The scaffolded analysis is intentionally a placeholder. The goal here is not to
design solvent-shell science yet. The goal is to learn where a new analysis
lives, which public APIs it imports, and how to verify that the generated plugin
starts from a working state.

## Before you start

You need:

- A working contributor environment from {doc}`../setup`.
- A feature branch or scratch checkout where it is safe to create throwaway
  files.
- The existing architecture orientation from {doc}`architecture`.

Use the throwaway plugin name `solvent_shell` exactly as written. If you keep
working on a real contribution later, choose a name that describes your actual
analysis.

```{important}
Do not commit the generated `solvent_shell` files unless you are intentionally
turning this tutorial scaffold into a real plugin. For this tutorial, treat them
as scratch files and remove them before committing unrelated documentation or
code changes.
```

## 1. Check that the scaffold command is available

From the repository root, ask the CLI for the command help:

```bash
pixi run -e build polyzymd new-analysis --help
```

You should see help for `polyzymd new-analysis`, including the default files it
creates and the generated-test command. If this command is not available, return
to {doc}`../setup` and make sure you are using the `build` pixi environment.

## 2. Generate the default scaffold

Before generating files, check whether the tutorial paths already exist:

```bash
git status --short src/polyzymd/analyses/solvent_shell.py tests/analyses/plugins/test_solvent_shell.py
```

If either path is already present, stop and inspect it before continuing. The
scaffold may fail, or you may need to understand the scaffold's `--force` option
before replacing anything. Do not overwrite work you intend to keep.

Run the scaffold with the default options:

```bash
pixi run -e build polyzymd new-analysis solvent_shell
```

The command creates a default MDAnalysis-native, single-file plugin and matching
tests.

## 3. Confirm the generated files

Confirm that the scaffold created exactly these two tutorial files:

```text
src/polyzymd/analyses/solvent_shell.py
tests/analyses/plugins/test_solvent_shell.py
```

That is the expected outcome for this tutorial: one analysis module and one test
module.

You can also check the same paths with Git:

```bash
git status --short src/polyzymd/analyses/solvent_shell.py tests/analyses/plugins/test_solvent_shell.py
```

```{warning}
These generated files are scratch content. Do not commit them unless you are
intentionally continuing from this scaffold as a real `solvent_shell` plugin.
```

## 4. Inspect the generated plugin

Open `src/polyzymd/analyses/solvent_shell.py` and skim it from top to bottom.
You do not need to edit anything yet.

### Settings

Find the generated settings model. It should include placeholder options for the
first working scaffold:

- `selection` with default value `"protein and name CA"`
- `scale` with default value `1.0`

These settings are intentionally generic. In a real plugin, replace or extend
them with options that describe your analysis.

### Placeholder measurement

Find the small placeholder measurement function. It receives a loaded universe
and frame-selection arguments, then computes a scaffold metric rather than a
scientific solvent-shell measurement.

The generated metric is named `solvent_shell_value`.

```{note}
Treat the placeholder measurement as a wiring check. Replace it with
domain-specific MDAnalysis logic before developing a scientific plugin.
```

### MDAnalysis job and collector

Look for the MDAnalysis job wiring:

- `MDAAnalysisJob` wraps the per-trajectory work.
- Collector logic converts completed job output into a `ReplicateArtifact`.
- The collector keeps raw trajectory handling out of later aggregation and
  comparison steps.

This is the scaffold's main lesson: trajectory-native plugins should hand
per-replicate results to PolyzyMD's artifact lifecycle instead of inventing a
separate cache format.

### Analysis class

Find the `SolventShellAnalysis` class. It should define:

- `name = "solvent_shell"`
- a `Settings` class
- `build_mda_jobs()`
- `build_mda_collector()`
- `extract_metrics()` when metric metadata needs customization

The generated `extract_metrics()` reads canonical condition-artifact payloads for
the default artifact comparison contract.

### Metrics and aggregation

The scaffold includes small private helpers that keep the generated example
complete enough to test:

- count selected frames for the placeholder calculation
- validate JSON-compatible metrics
- read default aggregate summaries

The generated plugin exposes the scalar metric through `MetricValue`, and the
default aggregation path can combine replicate artifacts into condition-level
results.

### Public import boundaries

The most important breadcrumb is the import boundary. The generated plugin uses
public contributor-facing imports such as:

- `polyzymd.analyses.base` for the `Analysis` subclass and scalar
  `MetricValue` descriptors.
- `polyzymd.analyses.mda` for `MDAAnalysisJob`, MDAnalysis job contexts,
  collectors, and artifact objects.

It should not import from `polyzymd.analyses._framework`. That package is an
internal implementation detail behind the public facades described in
{doc}`architecture`.

## 5. Inspect the generated tests

Open `tests/analyses/plugins/test_solvent_shell.py` next. The tests are part of
the scaffold because a new plugin contribution should start with executable
contract checks, not only a Python module.

### Discovery and settings tests

Skim the first tests for discovery and settings behavior. They prove that:

- the plugin can be discovered by the name `solvent_shell`
- the discovered class subclasses `Analysis`
- generated settings defaults and validation rules work

### MDAnalysis job tests

Find the tests for `build_mda_jobs()`. They prove that the generated plugin can
return an `MDAAnalysisJob` and that the placeholder measurement can run against
small fake universe objects.

### Collector and aggregation tests

Find the tests for collector and aggregation behavior. They prove that:

- the generated collector returns a `ReplicateArtifact`
- default aggregation can read replicate artifacts
- aggregated output exposes scalar metrics for comparison

### Why the tests use fakes

The tests use small fakes instead of real trajectories. That keeps the first
success fast and reviewable.

The fakes are enough for this scaffold because the tutorial checks plugin
wiring, artifact flow, and default aggregation behavior. Real analysis tests can
add trajectory-backed cases later when the scientific implementation exists.

## 6. Run the generated tests

Run only the generated test file:

```bash
pixi run -e build pytest tests/analyses/plugins/test_solvent_shell.py -v
```

A successful run confirms that the scaffolded plugin is discoverable, can build
its MDAnalysis job, can collect a replicate artifact, and can participate in the
default aggregation path.

Your exact pytest timing may differ, but the important success state is that the
generated `test_solvent_shell.py` tests pass.

## 7. Clean up the scratch scaffold

For this tutorial, `solvent_shell` is a throwaway name. Before you commit any
real work, either continue developing it into an intentional plugin or remove the
generated files.

If you are done with the tutorial, remove only these two files:

```bash
rm src/polyzymd/analyses/solvent_shell.py
rm tests/analyses/plugins/test_solvent_shell.py
```

Then confirm your checkout no longer contains the scratch scaffold:

```bash
git status --short src/polyzymd/analyses/solvent_shell.py tests/analyses/plugins/test_solvent_shell.py
```

These are the only files this tutorial asked you to create:

```text
src/polyzymd/analyses/solvent_shell.py
tests/analyses/plugins/test_solvent_shell.py
```

The scaffold is valuable because it gives you a known-good starting point. A
real contribution still needs domain-specific settings, trajectory logic,
artifact payloads, metrics, plots if needed, and tests that match the scientific
question.

## Common mistakes

- **Running the command outside the pixi environment.** Use
  `pixi run -e build` for every PolyzyMD scaffold and test command.
- **Starting from a real analysis name too early.** Use `solvent_shell` as the
  throwaway learning name, then choose an intentional name for real work.
- **Importing from private framework modules.** Contributor plugins should use
  `polyzymd.analyses.base` and `polyzymd.analyses.mda`, not
  `polyzymd.analyses._framework`.
- **Treating the placeholder metric as science.** The generated atom-frame count
  is only a scaffold. Replace it before opening a scientific plugin PR.
- **Committing scratch scaffold files.** If this was only a tutorial run, remove
  `solvent_shell.py` and `test_solvent_shell.py` before committing other work.

## What to read next

- {doc}`architecture` for why the scaffold uses MDAnalysis jobs, collectors, and
  artifacts.
- {doc}`../extending_analyses` for the full current implementation guide when
  you are ready to turn a scaffold into a real analysis.
