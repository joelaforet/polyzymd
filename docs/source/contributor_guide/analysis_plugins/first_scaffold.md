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

Run the scaffold with the default options:

```bash
pixi run -e build polyzymd new-analysis solvent_shell
```

The command creates a default MDAnalysis-native, single-file plugin and matching
tests. The generated paths are:

```text
src/polyzymd/analyses/solvent_shell.py
tests/analyses/plugins/test_solvent_shell.py
```

That is the expected outcome for this tutorial: one analysis module and one test
module.

## 3. Tour the generated plugin file

Open `src/polyzymd/analyses/solvent_shell.py` and skim it from top to bottom.
You do not need to edit anything yet.

Notice these pieces:

| Generated piece | What to look for |
| --- | --- |
| Settings model | A Pydantic settings class with placeholder options such as an MDAnalysis atom selection and a scale factor. |
| Measurement function | A small placeholder function that receives a loaded universe and frame-selection arguments. This is where domain-specific MDAnalysis logic will eventually replace the placeholder atom-frame count. |
| Collector | A collector class that converts one completed MDAnalysis job into a `ReplicateArtifact`. |
| Analysis subclass | A `SolventShellAnalysis` class with `name = "solvent_shell"`, a `Settings` class, `build_mda_jobs()`, `build_mda_collector()`, and `extract_metrics()`. The generated `extract_metrics()` is compatibility/fallback support for non-artifact summaries; default artifact comparison reads condition-artifact metrics first. |
| Helper functions | Small private helpers used by the scaffold to count selected frames, validate JSON-compatible metrics, and read default aggregate summaries. |

The most important breadcrumb is the import boundary. The generated plugin uses
public contributor-facing imports such as:

- `polyzymd.analyses.base` for the `Analysis` subclass and scalar
  `MetricValue` descriptors.
- `polyzymd.analyses.mda` for `MDAAnalysisJob`, MDAnalysis job contexts,
  collectors, and artifact objects.

It should not import from `polyzymd.analyses._framework`. That package is an
internal implementation detail behind the public facades described in
{doc}`architecture`.

## 4. Tour the generated tests

Open `tests/analyses/plugins/test_solvent_shell.py` next. The tests are part of
the scaffold because a new plugin contribution should start with executable
contract checks, not only a Python module.

Skim for these test groups:

| Test group | What it proves |
| --- | --- |
| Discovery tests | The plugin can be discovered by name and subclasses `Analysis`. |
| Settings tests | The generated settings defaults and validation rules work. |
| MDAnalysis job tests | `build_mda_jobs()` returns an `MDAAnalysisJob` and the placeholder function can run against fake universe objects. |
| Collector tests | The generated collector returns a `ReplicateArtifact`. |
| Aggregation tests | Default aggregation can read replicate artifacts and expose scalar metrics for comparison. |

The tests use small fakes instead of real trajectories. That keeps the first
success fast and reviewable.

## 5. Run the generated tests

Run only the generated test file:

```bash
pixi run -e build pytest tests/analyses/plugins/test_solvent_shell.py -v
```

A successful run confirms that the scaffolded plugin is discoverable, can build
its MDAnalysis job, can collect a replicate artifact, and can participate in the
default aggregation path.

Your exact pytest timing may differ, but the important success state is that the
generated `test_solvent_shell.py` tests pass.

## 6. Leave your checkout clean

For this tutorial, `solvent_shell` is a throwaway name. Before you commit any
real work, either continue developing it into an intentional plugin or remove the
generated files:

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
  `polyzymd.analyses.base`, `polyzymd.analyses.mda`, and documented shared
  utilities, not `polyzymd.analyses._framework`.
- **Treating the placeholder metric as science.** The generated atom-frame count
  is only a scaffold. Replace it before opening a scientific plugin PR.
- **Committing scratch scaffold files.** If this was only a tutorial run, remove
  `solvent_shell.py` and `test_solvent_shell.py` before committing other work.

## What to read next

- {doc}`architecture` for why the scaffold uses MDAnalysis jobs, collectors, and
  artifacts.
- {doc}`../extending_analyses` for the full current implementation guide when
  you are ready to turn a scaffold into a real analysis.
