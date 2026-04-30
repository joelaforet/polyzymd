# Extending Analyses

This guide shows the supported contributor workflow for adding a new analysis
plugin. New public plugins use the **runner-backed** lifecycle: PolyzyMD owns
replicate discovery, caching, aggregation, comparison, plotting, and CLI wiring;
your plugin supplies a small trajectory runner plus summary and aggregation
logic.

## Start with the scaffold

Use the scaffold command for a working package and tests:

```bash
pixi run -e build polyzymd new-analysis solvent_shell
pixi run -e build pytest tests/analyses/plugins/test_solvent_shell.py -v
```

The generated plugin is a small package under `src/polyzymd/analyses/<name>/`:

- `__init__.py` contains the `Analysis` facade and lifecycle wiring
- `_runner.py` contains the trajectory runner
- `_results.py` is added for `--style pydantic` result models

The package is discovered automatically by `pkgutil`; no registry edits are
needed.

Choose the result-container style with `--style`:

```bash
pixi run -e build polyzymd new-analysis solvent_shell --style dict
pixi run -e build polyzymd new-analysis solvent_shell --style pydantic
```

Both styles use the same lifecycle. The style only changes whether replicate
and aggregated results are plain dicts or Pydantic models stored in
`_results.py`.

## The public lifecycle

Implement these methods on your `Analysis` subclass:

| Hook | Purpose |
|------|---------|
| `build_runner(ctx, replicate, universe, window)` | Construct a trajectory runner for one replicate |
| `summarize_replicate(ctx, replicate, runner, window)` | Convert executed runner output into a dict or Pydantic result |
| `aggregate(ctx, results)` | Combine replicate results for one condition |
| `extract_metrics(summary)` | Expose scalar metrics for the default comparison path |
| `plot(ctx)` | Generate optional figures |

Keep the base per-replicate dispatch inherited. The base class loads the
Universe, resolves the trajectory window, calls your runner's `run(...)`, and
passes the executed runner to `summarize_replicate()`.

## Minimal plugin shape

```python
from __future__ import annotations

from typing import Any, ClassVar, Sequence

from pydantic import BaseModel, Field

from polyzymd.analyses.base import AggregateContext, Analysis, MetricValue, ReplicateContext


class MyMetricRunner:
    """Trajectory runner for one replicate."""

    def __init__(self, universe: Any, settings: "MyMetricAnalysis.Settings") -> None:
        self.universe = universe
        self.settings = settings
        self.results: dict[str, Any] = {}

    def run(self, start: int, stop: int | None, step: int = 1) -> "MyMetricRunner":
        stop_frame = len(self.universe.trajectory) if stop is None else stop
        n_frames = len(range(start, stop_frame, step))
        self.results = {"value": float(n_frames), "n_frames": n_frames}
        return self


class MyMetricAnalysis(Analysis):
    """Example runner-backed analysis."""

    name: ClassVar[str] = "my_metric"

    class Settings(BaseModel):
        """Settings for the analysis."""

        example_parameter: str = Field(default="dummy")

    def build_runner(
        self,
        ctx: ReplicateContext,
        replicate: int,
        universe: Any,
        window: Any,
    ) -> MyMetricRunner:
        """Build a runner for one replicate."""
        del replicate, window
        return MyMetricRunner(universe=universe, settings=ctx.settings)

    def summarize_replicate(
        self,
        ctx: ReplicateContext,
        replicate: int,
        runner: Any,
        window: Any,
    ) -> dict[str, Any]:
        """Summarize an executed runner."""
        del ctx, window
        return {
            "replicate": replicate,
            "value": float(runner.results["value"]),
            "n_frames": int(runner.results["n_frames"]),
        }

    def aggregate(self, ctx: AggregateContext, results: Sequence[Any]) -> dict[str, Any]:
        """Aggregate replicate values for one condition."""
        del ctx
        values = [float(result["value"]) for result in results]
        mean = sum(values) / len(values)
        sem = 0.0
        if len(values) > 1:
            variance = sum((value - mean) ** 2 for value in values) / (len(values) - 1)
            sem = (variance**0.5) / (len(values) ** 0.5)
        return {
            "mean_value": mean,
            "sem_value": sem,
            "replicate_values": values,
            "n_replicates": len(values),
        }

    def extract_metrics(self, summary: Any) -> dict[str, MetricValue]:
        """Expose metrics to the default comparison workflow."""
        return {
            "value": MetricValue(
                name="value",
                mean=float(summary["mean_value"]),
                sem=float(summary["sem_value"]),
                replicate_values=[float(value) for value in summary["replicate_values"]],
                higher_is_better=True,
            )
        }
```

The scaffold includes plotting and formatting examples as well as tests for the
base dispatch path.

## Pydantic result models

For more structured output, use `--style pydantic`. The scaffold writes result
models to `_results.py` and re-exports them from the package facade. The
generated class sets `ReplicateResultClass` and `AggregatedResultClass` so
framework serialization and deserialization use your Pydantic models.

Avoid manually documenting Pydantic fields in class docstrings; Sphinx renders
field documentation automatically.

## Trajectory access and imports

The base dispatch creates a `TrajectoryLoader` from `ctx.sim_config` and calls
`loader.load_universe(replicate)`. Heavy dependencies such as MDAnalysis,
OpenMM, OpenFF, ParmEd, and PDBFixer must be imported lazily inside functions
or methods, including functions and methods defined in `_runner.py`. Do not
import these packages at module level in any plugin module.

Selection strings are passed to MDAnalysis `Universe.select_atoms()` unless a
plugin explicitly documents a plugin-specific wrapper such as `com(...)` for
distance pairs or `midpoint(...)` for catalytic-triad pairs. Common structural
examples are:

- `protein`
- `protein and name CA`
- `chainID C`
- `resname SBM`
- `protein and (resid 77 or resid 156)`

Validate selections on representative topology files, because available
attributes depend on the trajectory topology.

## Condition filtering

Use `filter_conditions()` when an entire analysis is not applicable to a
condition, for example a polymer-only analysis on a no-polymer control.
Otherwise, handle empty selections in the runner or summary stage and return a
well-documented empty or NaN result.

## Comparison choices

Use `extract_metrics()` when your aggregated result exposes scalar metrics. The
default comparison workflow loads each condition, performs the configured
pairwise tests and ANOVA where applicable, builds rankings, and saves a
`ComparisonResult`.

Override `compare()` only when the default scalar path cannot represent your
result, such as entry tables, multiple nested run labels, or custom output
models.

(plotters-extraction)=
## Plotting organization

Small plugins can keep plotting in `plot(ctx)`. As plotting grows, move helper
functions into `src/polyzymd/analyses/<name>/_plotters.py` and keep `plot(ctx)`
as the lifecycle wrapper that loads results and calls those helpers. This is the
established pattern for plugins with several figure types.

## Testing checklist

Generated tests cover the standard contributor contract:

- discovery and class variables
- settings defaults and validation
- `build_runner()` smoke behavior
- `summarize_replicate()` with a fake executed runner
- base per-replicate dispatch with a fake loader, Universe, and trajectory window
- `aggregate()` and `extract_metrics()`
- `plot()` with cached aggregated results

Run plugin tests through the pixi environment:

```bash
pixi run -e build pytest tests/analyses/plugins/test_<name>.py -v
```

## Style checklist

- Use NumPy-style docstrings for new classes and methods
- Keep imports ordered stdlib, third-party, local
- Keep heavy scientific dependencies lazy
- Use `X | None` annotations rather than `Optional[X]`
- Run Ruff and Black checks on modified files

```bash
pixi run -e build ruff check src/polyzymd/analyses/<name> tests/analyses/plugins/test_<name>.py
pixi run -e build black src/polyzymd/analyses/<name> tests/analyses/plugins/test_<name>.py --check
```
