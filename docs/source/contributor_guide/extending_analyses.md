# Extending Analyses

This guide shows the supported contributor workflow for adding a new analysis
plugin. The default path is now a **single-file scalar measurement plugin**:
your code measures one scalar value per replicate, and PolyzyMD handles the
runner bridge, aggregation, cache identity, default scalar comparison, CLI
wiring, and discovery.

Use the advanced runner-backed scaffold only when one scalar measurement is not
enough for the science you need to express.

## Start with the measurement scaffold

Create a working plugin and tests with:

```bash
pixi run -e build polyzymd new-analysis solvent_shell
pixi run -e build pytest tests/analyses/plugins/test_solvent_shell.py -v
```

The default scaffold creates:

```text
src/polyzymd/analyses/solvent_shell.py
tests/analyses/plugins/test_solvent_shell.py
```

The plugin is discovered automatically. Single-file plugins named
`src/polyzymd/analyses/<name>.py` and package plugins under
`src/polyzymd/analyses/<name>/` both participate in discovery; no registry edits,
decorators, or bootstrap imports are needed.

## Public imports

Contributor plugins should import the measurement API from the public facade:

```python
from polyzymd.analyses.base import MetricSpec, ScalarMeasurement, ScalarMeasurementAnalysis
```

`polyzymd.analyses.base` is the stable contributor surface. It re-exports the
measurement API, the `Analysis` base class, lifecycle contexts, metric models,
and comparison result models while implementation details remain private. Do not
import from modules such as `_measurement.py`, `_analysis_runner.py`,
`_contexts.py`, or `_comparison_models.py` in contributor plugins.

## Minimal scalar measurement plugin

A scalar measurement plugin has three pieces:

- a Pydantic settings model
- a `ScalarMeasurement` that implements `measure()`
- a `ScalarMeasurementAnalysis` that binds the analysis name, settings, and
  measurement

```python
from __future__ import annotations

from typing import Any, ClassVar

from pydantic import BaseModel, Field

from polyzymd.analyses.base import MetricSpec, ScalarMeasurement, ScalarMeasurementAnalysis


class SolventShellSettings(BaseModel):
    """Settings for solvent shell analysis."""

    selection: str = Field(default="protein and name CA")
    cutoff: float = Field(default=5.0, gt=0.0)


class SolventShellMeasurement(ScalarMeasurement):
    """Measure one scalar value from one replicate trajectory."""

    name: ClassVar[str] = "solvent_shell_measurement"
    version: ClassVar[str] = "1"
    metric: ClassVar[MetricSpec] = MetricSpec(
        name="mean_shell_count",
        higher_is_better=False,
        label="Mean shell count",
        unit="atoms",
        direction_labels=("decreased", "unchanged", "increased"),
    )

    def measure(
        self,
        universe: Any,
        settings: BaseModel,
        *,
        start: int | None = None,
        stop: int | None = None,
        step: int | None = None,
    ) -> float:
        """Return one scalar value for the requested trajectory window."""
        if not isinstance(settings, SolventShellSettings):
            raise TypeError("settings must be a SolventShellSettings instance.")

        atoms = universe.select_atoms(settings.selection)
        trajectory = universe.trajectory[start:stop:step]
        values = []
        for _frame in trajectory:
            values.append(float(len(atoms)))
        return sum(values) / len(values) if values else 0.0


class SolventShellAnalysis(ScalarMeasurementAnalysis):
    """Solvent shell analysis backed by one scalar measurement."""

    name: ClassVar[str] = "solvent_shell"
    Settings: ClassVar[type[BaseModel]] = SolventShellSettings
    measurement: ClassVar[type[ScalarMeasurement]] = SolventShellMeasurement
```

Replace the placeholder measurement logic with your scientific calculation.
`measure()` receives the MDAnalysis `Universe`, resolved settings, and the frame
window chosen by the comparison workflow. It must return a single `float` for
that replicate.

## What the framework handles

`ScalarMeasurementAnalysis` adapts the simple `measure()` method to the full
analysis lifecycle:

| Lifecycle concern | Handled by PolyzyMD |
|-------------------|---------------------|
| Trajectory loading | The framework creates the `TrajectoryLoader` from the run context |
| Runner bridge | The scalar adapter wraps `measure()` in a runner-compatible object |
| Replicate result | The adapter stores analysis, measurement, metric, replicate, value, and cache identity |
| Aggregation | Replicate scalar values become mean, SEM, replicate values, and counts |
| Cache identity | Measurement name, version, metric, settings, and analysis name are fingerprinted |
| Comparison | The default scalar comparison path builds rankings, pairwise tests, and ANOVA where applicable |

You usually do not implement `build_runner()`, `summarize_replicate()`,
`aggregate()`, or `extract_metrics()` for scalar measurement plugins.

## Measurement metadata and cache identity

Use `MetricSpec` to describe the scalar produced by the measurement:

- `name` is the stable serialized metric key
- `higher_is_better` controls ranking direction
- `label` and `unit` support human-readable output
- `direction_labels` describe decreased, unchanged, and increased effects

Set `measurement.version` when changing the meaning of the measurement. The
default cache identity includes the measurement name, version, metric name,
analysis name, and settings payload, so stale scalar replicate results are
rejected when those inputs change.

## Trajectory access and imports

Heavy dependencies such as MDAnalysis, OpenMM, OpenFF, ParmEd, and PDBFixer must
be imported lazily inside functions or methods. Do not import them at module
level in plugin files.

Selection strings are passed to MDAnalysis `Universe.select_atoms()` unless your
plugin documents a wrapper such as `com(...)` for center-of-mass selections.
Common examples are:

- `protein`
- `protein and name CA`
- `chainid A`
- `chainid C`
- `resname SBM`
- `protein and (resid 77 or resid 156)`

Validate selections on representative topologies because available attributes
depend on the trajectory topology.

## Advanced runner-backed scaffolds

Use the advanced scaffold when your plugin needs multi-metric outputs,
per-frame/per-residue tables, custom result models, sidecar files, or a
trajectory runner that cannot be expressed as one scalar `measure()` call.

```bash
pixi run -e build polyzymd new-analysis solvent_shell --advanced
pixi run -e build polyzymd new-analysis solvent_shell --style dict
pixi run -e build polyzymd new-analysis solvent_shell --style pydantic
```

`--advanced` creates the runner-backed package scaffold. If no style is given,
it uses the plain-dict advanced style. Passing `--style dict` or
`--style pydantic` also selects advanced package scaffolding for backward
compatibility.

Advanced packages create files such as:

```text
src/polyzymd/analyses/<name>/__init__.py
src/polyzymd/analyses/<name>/_runner.py
src/polyzymd/analyses/<name>/_results.py  # pydantic style only
```

Advanced plugins implement the runner-backed lifecycle:

| Hook | Purpose |
|------|---------|
| `build_runner(ctx, replicate, universe, window)` | Construct a trajectory runner for one replicate |
| `summarize_replicate(ctx, replicate, runner, window)` | Convert executed runner output into a dict or Pydantic result |
| `aggregate(ctx, results)` | Combine replicate results for one condition |
| `extract_metrics(summary)` | Expose scalar metrics for the default comparison path |

Keep the inherited per-replicate dispatch unless you need plugin-specific cache
or sidecar behavior. Compare-only plugins can set `has_compute_stage = False`.

## Pydantic result models

Use `--style pydantic` for advanced plugins that need typed replicate or
aggregated result models. The scaffold writes models to `_results.py` and sets
`ReplicateResultClass` and `AggregatedResultClass` so serialization and default
deserialization use those models.

Avoid manually documenting Pydantic fields in class docstrings; Sphinx renders
field documentation automatically.

## Custom comparison guidance

Use the default scalar comparison path whenever your aggregate result can expose
one or more scalar `MetricValue` objects. For scalar measurement plugins, this is
already implemented by `ScalarMeasurementAnalysis`.

Override `compare()` only when the default scalar path cannot represent the
result, such as:

- entry tables or per-residue hypothesis families
- multiple named runs per condition
- custom statistical models or output schemas
- comparisons that need to load sidecar files in addition to `result.json`

Keep custom comparison logic in the plugin package, return a saveable Pydantic
model or framework `ComparisonResult`, and implement `format()` if CLI output
needs more than the default JSON-style rendering.

(plotters-extraction)=
## Plotting guidance

Small plugins can keep plotting directly in `plot(ctx)`. Use the context-provided
output directory, load cached aggregate or comparison results through the plugin
helpers, and return the list of figure paths written.

Extract helpers into `_plotters.py` only when plotting grows enough to justify a
separate module, for example when a plugin has several figure types or the
analysis file becomes hard to review. Established package plugins may also use
`_formatters.py` for custom CLI text. These private helper modules are optional
organization tools, not the default contributor pattern.

## Testing checklist

Generated scalar measurement tests focus on the contributor contract:

- discovery and class variables
- settings defaults and validation
- direct `measure()` behavior with a fake Universe
- base scalar runner dispatch through the framework bridge
- aggregation and default metric extraction inherited from
  `ScalarMeasurementAnalysis`

Run plugin tests through the pixi environment:

```bash
pixi run -e build pytest tests/analyses/plugins/test_<name>.py -v
```

## Style checklist

- Use NumPy-style docstrings for new classes and methods
- Keep imports ordered stdlib, third-party, local
- Keep heavy scientific dependencies lazy
- Use `X | None` annotations rather than `Optional[X]`
- Run Ruff and Black checks on modified Python files

```bash
pixi run -e build ruff check src/polyzymd/analyses/<name>.py tests/analyses/plugins/test_<name>.py
pixi run -e build black src/polyzymd/analyses/<name>.py tests/analyses/plugins/test_<name>.py --check
```
