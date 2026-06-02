# Analysis plugin contribution checklist

Use this page to copy a review checklist into a GitHub issue or pull request
body. The rendered documentation shows the instructions, while the checklist
itself lives in a copy-pasteable template so contributors can check items in
GitHub.

## How to use this checklist

Replace `PLUGIN_NAME` with your plugin's stable lowercase name, paste the
template into a GitHub issue or PR body, and check boxes in GitHub as you work.
Delete irrelevant items or mark them `N/A` when a lifecycle stage does not apply
to your plugin.

## What this checklist covers

The template covers the expected contributor contract for analysis plugins:
source placement, public imports, lifecycle hooks, MDAnalysis jobs, artifacts,
aggregation, comparison, plotting, formatting, tests, user-facing documentation,
and pre-PR commands.

## Copy-paste checklist template

````markdown
## Analysis plugin checklist

Plugin: `PLUGIN_NAME`

Use this checklist before opening or reviewing an analysis plugin pull request.
Delete items that do not apply, or mark them `N/A` with a short reason.

### Placement and shape

- [ ] The plugin lives under `src/polyzymd/analyses/` as either one readable
  module or one package with a public `__init__.py`.
- [ ] The public plugin class subclasses `Analysis` and defines `name` and
  `Settings`.
- [ ] Advanced helper files such as `_mda.py`, `_plotters.py`, `_models.py`, and
  `_formatters.py` exist only when they have clear responsibilities.
- [ ] The plugin name is stable, lowercase, and matches CLI/config usage.

### Public imports

- [ ] Contributor-facing code imports lifecycle classes from
  `polyzymd.analyses.base`.
- [ ] MDAnalysis job, artifact, sidecar, and artifact-store types come from
  `polyzymd.analyses.mda`.
- [ ] Shared helpers come only from documented `polyzymd.analyses.shared`
  utilities when they fit the task.
- [ ] The plugin does not import from `polyzymd.analyses._framework`; that
  package is private behind the public facades.
- [ ] Heavy dependencies such as MDAnalysis and matplotlib are imported lazily
  inside functions or methods that need them.

### Plugin contract

- [ ] `Settings` is a Pydantic model with defaults, field descriptions, and
  validation for invalid values.
- [ ] Compute-stage plugins implement `build_mda_jobs()` and, when needed,
  `build_mda_collector()`.
- [ ] Compare-only plugins explicitly set the appropriate lifecycle flags and do
  not pretend to run trajectory work.
- [ ] Lifecycle methods use framework-provided context objects; they do not reload
  YAML configuration on their own.
- [ ] MDAnalysis selections use project conventions such as `chainid A` for
  protein, `chainid B` for substrate, and `chainid C` for polymer.

### MDAnalysis jobs

- [ ] Jobs are `MDAAnalysisJob` objects built from the public MDAnalysis layer.
- [ ] Function-adapter jobs respect `frames` and `start`/`stop`/`step` frame
  selection kwargs.
- [ ] `AnalysisBase`-compatible jobs can be tested with small fake workers.
- [ ] Job metadata records enough provenance to explain selections, frame windows,
  settings fingerprints, and relevant algorithm choices.
- [ ] Job results are converted to JSON-compatible values or sidecars before they
  become artifacts.

### Artifacts and sidecars

- [ ] Collectors return `ReplicateArtifact` objects.
- [ ] Scalar values intended for default aggregation are finite numbers under
  `payload["metrics"]`.
- [ ] Payloads remain compact and JSON-compatible.
- [ ] Large arrays, event tables, and frame-by-frame data are stored as registered
  sidecars, not inline JSON.
- [ ] Sidecar paths are store-relative and are loaded through `ArtifactStore`
  validation helpers.
- [ ] Raw MDAnalysis `Results` objects are never serialized.

### Aggregation and comparison

- [ ] Simple scalar artifact plugins rely on default aggregation to build a
  `ConditionArtifact` when that contract is sufficient.
- [ ] Custom aggregation is used only when condition-level data needs a custom
  schema or sidecar handling.
- [ ] `extract_metrics()` implementations read canonical `ConditionArtifact`
  payloads and return `MetricValue` objects for default scalar comparison.
- [ ] Custom `compare()` methods return a saveable comparison result or artifact
  and have tests for missing or incomplete conditions.
- [ ] Metric names, units, and higher/lower-is-better interpretation are explicit
  where scientific comparison depends on them.

### Plotting and formatting

- [ ] `plot()` reads cached artifacts and registered sidecars only.
- [ ] Plotting does not load trajectories, run MDAnalysis jobs, or recompute
  analysis results.
- [ ] Plot helpers write figures under the provided `PlotContext.output_dir`.
- [ ] `format()` renders existing comparison output; it does not compute new
  scientific quantities.
- [ ] Plot and format behavior is covered by focused tests or intentionally left
  absent for plugins that do not provide them.

### Tests

- [ ] Discovery test covers `list_analyses()` and `get_analysis()`.
- [ ] Settings tests cover defaults and invalid values.
- [ ] `build_mda_jobs()` tests use fakes or mocks and verify frame selection and
  settings wiring.
- [ ] Collector tests assert that a `ReplicateArtifact` is returned with expected
  payload, provenance, warnings, and sidecars.
- [ ] Aggregation or metric-extraction tests cover condition-level outputs.
- [ ] Comparison tests cover the default metric path or custom comparison path.
- [ ] Plot tests are artifact-only and do not require trajectory data.
- [ ] Sidecar tests validate registered sidecars through `ArtifactStore` when the
  plugin writes arrays or tables.
- [ ] Tests that require real trajectory data are marked `@pytest.mark.slow`.

### Documentation and user-facing surfaces

- [ ] New plugins, settings, public APIs, CLI options, or user workflows update
  the relevant docs, reference pages, CLI help, and configuration examples.
- [ ] New analysis settings are documented where users look up
  `comparison.yaml` plugin options.
- [ ] Stable versus experimental behavior is labeled clearly when the feature is
  not yet a settled default workflow.

### Commands before PR

Run the focused checks for your plugin and the project-level checks relevant to
analysis contributions:

```bash
PLUGIN_NAME=my_plugin

pixi run -e build pytest "tests/analyses/plugins/test_${PLUGIN_NAME}.py" -v
pixi run -e build pytest tests/analyses/ -v
pixi run -e build ruff check src/ "tests/analyses/plugins/test_${PLUGIN_NAME}.py"
pixi run -e build black src/ "tests/analyses/plugins/test_${PLUGIN_NAME}.py" --check
pixi run -e build make -C docs clean html
```

If you changed documentation, the clean Sphinx build must complete with zero
warnings before review.
````
