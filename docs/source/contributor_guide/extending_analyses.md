# Extending Analyses

This guide shows the supported contributor workflow for adding a new analysis
plugin. The default path is now a **single-file MDAnalysis-native plugin**: your
code builds an `MDAAnalysisJob`, converts completed job results into a canonical
`ReplicateArtifact`, and lets PolyzyMD handle cache, aggregation, comparison,
CLI wiring, and discovery.

Advanced scaffolds generate packages with lifecycle wiring in `__init__.py` and
lazy MDAnalysis `AnalysisBase` helpers in `_mda.py`.

## Start with the MDAnalysis-native scaffold

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

Contributor plugins should import lifecycle classes from the public facade and
MDA extension-layer helpers from `polyzymd.analyses.mda`:

```python
from polyzymd.analyses.base import Analysis, MetricValue
from polyzymd.analyses.mda import MDAAnalysisJob, ReplicateArtifact
```

`polyzymd.analyses.base` is the stable contributor surface. It re-exports the
`Analysis` base class, lifecycle contexts, metric models, and comparison result
models while implementation details remain private. Do not import from modules
such as `_contexts.py` or `_comparison_models.py` in contributor plugins.

## Minimal MDAnalysis-native plugin

The default scaffold has four pieces:

- a Pydantic settings model
- a small function wrapped by `MDAAnalysisJob.from_function()`
- a collector that returns a `ReplicateArtifact` with `payload["metrics"]`
- an `Analysis` subclass implementing `build_mda_jobs()`,
  `build_mda_collector()`, and `extract_metrics()`

```python
from __future__ import annotations

from collections.abc import Sequence
from typing import Any, ClassVar

from pydantic import BaseModel, Field

from polyzymd.analyses.base import Analysis, MetricValue
from polyzymd.analyses.mda import MDAAnalysisJob, MDACollectorContext, ReplicateArtifact

METRIC_NAME = "mean_shell_count"


class SolventShellSettings(BaseModel):
    """Settings for solvent shell analysis."""

    selection: str = Field(default="protein and name CA")
    cutoff: float = Field(default=5.0, gt=0.0)


def measure_solvent_shell(universe: Any, *, settings: SolventShellSettings, **_kwargs: Any):
    """Return strict JSON-compatible job results."""
    atoms = universe.select_atoms(settings.selection)
    return {"metrics": {METRIC_NAME: float(len(atoms))}}


class SolventShellAnalysis(Analysis):
    """Solvent shell analysis backed by an MDAnalysis job."""

    name: ClassVar[str] = "solvent_shell"
    Settings: ClassVar[type[BaseModel]] = SolventShellSettings

    def build_mda_jobs(self, ctx):
        return [
            MDAAnalysisJob.from_function(
                name="solvent_shell",
                function=measure_solvent_shell,
                universe=ctx.universe,
                frame_selection=ctx.frame_selection,
                universe_policy=ctx.universe_policy,
                function_kwargs={"settings": ctx.settings},
            )
        ]

    def build_mda_collector(self, ctx: MDACollectorContext):
        def collect(ctx: MDACollectorContext, completed_jobs: Sequence[Any]) -> ReplicateArtifact:
            metrics = dict(completed_jobs[0].results["metrics"])
            return ReplicateArtifact(
                analysis_name=ctx.analysis_name,
                condition_label=ctx.condition_label,
                replicate=ctx.replicate,
                payload={"metrics": metrics},
                provenance={"source": "solvent_shell_scaffold"},
            )

        return collect

    def extract_metrics(self, summary: Any) -> dict[str, MetricValue]:
        metric = summary.payload["metrics"][METRIC_NAME]
        return {
            METRIC_NAME: MetricValue(
                name=METRIC_NAME,
                mean=metric["mean"],
                sem=metric["sem"],
                replicate_values=metric["values"],
                higher_is_better=False,
            )
        }
```

Replace the placeholder function and collector logic with your scientific
calculation. The collector should map raw job results to JSON-compatible payloads
or sidecars; do not persist files directly from the collector.

## What the framework handles

The generated `Analysis` subclass uses the MDAnalysis job lifecycle:

| Lifecycle concern | Handled by PolyzyMD |
|-------------------|---------------------|
| Trajectory loading | The framework creates the `TrajectoryLoader` from the run context |
| Replicate execution | The framework executes each `MDAAnalysisJob` with the resolved frame window |
| Replicate result | The collector returns a `ReplicateArtifact` persisted through `ArtifactStore` |
| Aggregation | Explicit `payload["metrics"]` values become mean, SEM, replicate values, and counts |
| Cache identity | Settings fingerprints are stored in artifact metadata |
| Comparison | The default scalar comparison path builds rankings, pairwise tests, and ANOVA where applicable |

The base `aggregate()` implementation handles `ReplicateArtifact` inputs when
each artifact declares explicit scalar metrics in `payload["metrics"]`.

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

## Advanced MDAnalysis-native plugins

Use an advanced package when your plugin needs multi-metric outputs,
per-frame/per-residue tables, custom result models, sidecar files, or
MDAnalysis-compatible jobs that cannot be expressed as one scalar `measure()`
call.

```bash
pixi run -e build polyzymd new-analysis solvent_shell --advanced
pixi run -e build polyzymd new-analysis solvent_shell --style dict
pixi run -e build polyzymd new-analysis solvent_shell --style pydantic
```

`--advanced`, `--style dict`, and `--style pydantic` generate MDAnalysis-native
packages. They do not generate or depend on the removed runner package.

Advanced packages use files such as:

```text
src/polyzymd/analyses/<name>/__init__.py
src/polyzymd/analyses/<name>/_mda.py
src/polyzymd/analyses/<name>/_results.py  # pydantic style only
```

Advanced plugins implement the MDAnalysis job lifecycle:

| Hook | Purpose |
|------|---------|
| `build_mda_jobs(ctx)` | Construct MDAnalysis-compatible jobs for one replicate |
| `build_mda_collector(ctx)` | Convert completed job outputs into a canonical artifact |
| `aggregate(ctx, results)` | Combine replicate results for one condition |
| `extract_metrics(summary)` | Expose scalar metrics for the default comparison path |

Keep the inherited per-replicate dispatch unless you need an explicit
`run_replicate()` override for advanced/internal behavior. Compare-only plugins
can set `has_compute_stage = False`.

## Pydantic result models

Use `--style pydantic` when the generated advanced package should validate its
replicate-level metric payload before storing it in a canonical artifact. The
scaffold writes a small helper model to `_results.py`, but it still persists
replicate and condition outputs through `ReplicateArtifact`, `ConditionArtifact`,
and `ArtifactStore`. It does not set `ReplicateResultClass` or
`AggregatedResultClass`.

Avoid manually documenting Pydantic fields in class docstrings; Sphinx renders
field documentation automatically.

## Custom comparison guidance

Use the default scalar comparison path whenever your aggregate result can expose
one or more scalar `MetricValue` objects. The generated scaffold implements
`extract_metrics()` for the artifact aggregate produced by default MDA
aggregation.

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

Generated MDAnalysis-native tests currently focus on the scaffold contract:

- discovery and class variables
- settings defaults and validation
- `build_mda_jobs()` behavior with fake universes or fake AnalysisBase objects
- collector output as `ReplicateArtifact`
- default artifact aggregation from explicit `payload["metrics"]`

For production plugins, consider adding targeted tests for optional behavior you
customize or rely on heavily:

- frame-selection behavior, Universe loading, and replicate cache identity when
  relevant
- aggregation and default metric extraction for any custom metric policies
- custom plotting, formatting, filtering, or cache identity logic

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
