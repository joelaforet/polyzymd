# Convert a plugin into an advanced package

This how-to is for contributors whose analysis plugin has outgrown a single
module. Use it when the plugin already works and the next task is to make its
structure easier to maintain, test, and review.

An advanced package is still one plugin. The public entry point remains one
`Analysis` subclass discovered by PolyzyMD. The package split only moves helper
code into focused private modules inside that plugin package.

## Decide whether to split

Do **not** split a plugin just because it might grow later. Keep a plugin as a
single module when:

- it has one small `Settings` model and one readable `Analysis` subclass;
- `build_mda_jobs()` and the collector fit on the same screen;
- plotting is absent or limited to one short helper;
- results are plain artifacts with compact payloads and no custom output model;
- formatting uses the default formatter or a short `format()` method; and
- tests can find the behavior without navigating several files.

Split into a package when the current module is becoming hard to review. Good
signals include:

- MDAnalysis job setup or collectors need several helper classes or functions;
- the plugin uses array or table sidecars and needs loading helpers;
- plotting has several artifact-only figures or repeated data preparation;
- aggregation, comparison, or CLI output uses structured Pydantic models;
- formatting has enough logic to distract from the lifecycle methods; or
- maintainers are repeatedly asking for clearer separation of responsibilities.

The goal is not abstraction for its own sake. The goal is a smaller public
facade in `__init__.py` plus private modules that each have one reason to
change.

## Use built-in packages as examples, not import targets

Built-in plugins such as `rmsf`, `distances`, `contacts`, `sasa`, `rg`, `rmsd`,
`hydrogen_bonds`, `secondary_structure`, and `catalytic_triad` show real package
shapes. Their private helper modules are useful examples of organization:

- `_mda.py` for trajectory-native job and collector helpers;
- `_plotters.py` for artifact-only plotting helpers;
- `_results.py` for structured result models;
- `_formatters.py` for substantial CLI formatting; and
- plugin-specific modules such as `_comparison.py` or `_filters.py` only when a
  plugin has a genuine need.

Do not import private modules from another plugin. For contributor plugins, use
the public facades only:

```python
from polyzymd.analyses.base import Analysis, PlotContext
from polyzymd.analyses.mda import MDAAnalysisJob, MDAReplicateJobContext
```

Documented utilities from `polyzymd.analyses.shared` are also acceptable when an
existing shared helper fits the task. Do not import from
`polyzymd.analyses._framework`; it is a private implementation detail behind the
public facades.

## Target package layout

For a plugin named `solvent_shell`, convert the single module into a directory
like this:

```text
src/polyzymd/analyses/solvent_shell/
├── __init__.py
├── _mda.py
├── _plotters.py
├── _results.py
└── _formatters.py
```

Create only the files that have a current job. If your plugin has no custom
formatting yet, skip `_formatters.py`. If artifacts and default comparison are
enough, skip `_results.py`.

| File | Responsibility | Keep out |
| --- | --- | --- |
| `__init__.py` | Public plugin class, settings model, class variables, lifecycle wiring, and small orchestration methods. | Long MDAnalysis loops, large plot functions, bulky formatters. |
| `_mda.py` | MDAnalysis job functions, `AnalysisBase`-compatible workers, collectors, sidecar writing, and helpers that translate runtime results into `ReplicateArtifact` objects. | Public plugin registration, CLI formatting, plot rendering. |
| `_plotters.py` | Artifact-only plot helpers that load `ConditionArtifact`, `ComparisonArtifact`, and registered sidecars, then write figures. | Trajectory loading, MDAnalysis job execution, recomputation. |
| `_results.py` | Pydantic models or typed result helpers when the plugin returns structured outputs beyond plain artifact payloads. | Raw MDAnalysis `Results` objects or large frame-by-frame arrays that belong in sidecars. |
| `_formatters.py` | Substantial text/table formatting used by `format()`. | Scientific computation, artifact writes, plotting side effects. |

## Keep `__init__.py` as the public facade

After the split, `__init__.py` should still be the place where discovery finds
the plugin class. Keep the public class readable and import private helpers from
the same package with relative imports.

```python
from typing import ClassVar

from pydantic import BaseModel, Field

from polyzymd.analyses.base import Analysis, PlotContext
from polyzymd.analyses.mda import MDACollectorContext, MDAReplicateJobContext

from ._formatters import format_solvent_shell
from ._mda import SolventShellCollector, build_solvent_shell_jobs
from ._plotters import plot_solvent_shell_summary


class SolventShellSettings(BaseModel):
    """Settings for solvent-shell analysis."""

    selection: str = Field(default="protein and chainid A")


class SolventShellAnalysis(Analysis):
    """Analyze solvent-shell behavior around the protein."""

    name: ClassVar[str] = "solvent_shell"
    Settings: ClassVar[type[BaseModel]] = SolventShellSettings

    def build_mda_jobs(self, ctx: MDAReplicateJobContext):
        return build_solvent_shell_jobs(ctx)

    def build_mda_collector(self, ctx: MDACollectorContext):
        del ctx
        return SolventShellCollector()

    def plot(self, ctx: PlotContext):
        return plot_solvent_shell_summary(ctx)

    def format(self, result, output_format: str = "text") -> str:
        return format_solvent_shell(result, output_format=output_format)
```

This file wires together the lifecycle but does not hide heavy computation or
plotting details inside the class body.

## Move trajectory work to `_mda.py`

Put trajectory-native details in `_mda.py`: job construction, function-adapter
jobs, `AnalysisBase`-compatible workers, collectors, and sidecar writes. Import
heavy dependencies lazily inside functions or methods that need them.

```python
from polyzymd.analyses.mda import (
    MDAAnalysisJob,
    MDACollectorContext,
    MDAJobResult,
    MDAReplicateJobContext,
    ReplicateArtifact,
    frame_selection_payload,
)


def calculate_shell_counts(universe, *, selection: str, **frame_kwargs):
    """Calculate compact per-replicate values for one trajectory."""

    import numpy as np

    atoms = universe.select_atoms(selection)
    values: list[float] = []
    frames = frame_kwargs.get("frames")
    iterator = universe.trajectory[frames] if frames is not None else universe.trajectory

    for _ts in iterator:
        values.append(float(np.asarray(atoms.positions).shape[0]))

    return {"mean_count": float(np.mean(values)), "n_frames": len(values)}


def build_solvent_shell_jobs(ctx: MDAReplicateJobContext) -> list[MDAAnalysisJob]:
    settings = ctx.settings
    return [
        MDAAnalysisJob.from_function(
            name="solvent_shell_counts",
            function=calculate_shell_counts,
            universe=ctx.universe,
            frame_selection=ctx.frame_selection,
            universe_policy=ctx.universe_policy,
            function_kwargs={"selection": settings.selection},
        )
    ]


class SolventShellCollector:
    """Collect one completed job into a replicate artifact."""

    def __call__(
        self,
        ctx: MDACollectorContext,
        completed_jobs: list[MDAJobResult],
    ) -> ReplicateArtifact:
        job = completed_jobs[0]
        mean_count = float(job.results["mean_count"])
        return ReplicateArtifact(
            analysis_name=ctx.analysis_name,
            condition_label=ctx.condition_label,
            replicate=ctx.replicate,
            payload={"metrics": {"mean_shell_count": mean_count}},
            provenance={"frame_selection": frame_selection_payload(ctx.frame_selection)},
            metadata={"result_kind": "solvent_shell_replicate"},
            warnings=list(ctx.warnings),
        )
```

The collector returns a durable artifact. It should not serialize raw MDAnalysis
`Results` objects. If the job produces arrays or event tables, write registered
sidecars as shown in {doc}`sidecars`.

## Move artifact-only plotting to `_plotters.py`

Plot helpers should read cached artifacts and registered sidecars only. They
should not load trajectories, create `MDAAnalysisJob` objects, or rerun compute
work.

```python
from pathlib import Path

from polyzymd.analyses.base import PlotContext
from polyzymd.analyses.mda import ArtifactStore


def plot_solvent_shell_summary(ctx: PlotContext) -> list[Path]:
    """Plot from cached condition artifacts and sidecars."""

    import matplotlib.pyplot as plt

    output_path = ctx.output_dir / "solvent_shell_summary.png"
    labels: list[str] = []
    values: list[float] = []

    for condition in ctx.conditions:
        analysis_dir = ctx.analysis_dirs.get(condition.label)
        if analysis_dir is None:
            continue
        artifact = ArtifactStore(analysis_dir).read_condition_result()
        metric = artifact.payload.get("metrics", {}).get("mean_shell_count")
        if metric is None:
            continue
        labels.append(condition.label)
        values.append(float(metric["mean"]))

    if not values:
        return []

    fig, ax = plt.subplots()
    ax.bar(labels, values)
    ax.set_ylabel("Mean shell count")
    fig.tight_layout()
    fig.savefig(output_path)
    plt.close(fig)
    return [output_path]
```

Keep all heavy plotting dependencies inside the plotting function. This preserves
fast imports for users who only need configuration, discovery, or API docs.

## Add `_results.py` only for structured outputs

Use `_results.py` when the plugin has a custom Pydantic output model or typed
helpers that make aggregation and comparison clearer. Do not add a result model
only to wrap one scalar already stored under `payload["metrics"]`.

Good uses for `_results.py` include:

- custom comparison outputs with several tables or ranked entries;
- aggregate summaries that need schema validation beyond artifact payloads;
- small typed models that point to sidecars by reference; and
- compatibility loaders for an established plugin-specific result format.

Large arrays, per-frame matrices, and event streams should remain in registered
sidecars, with `_results.py` storing only validated summaries or sidecar
references.

## Add `_formatters.py` only for substantial formatting

Use `_formatters.py` when `format()` has enough branches or table-building logic
to obscure the `Analysis` subclass. Keep `format()` in `__init__.py` as a small
delegation method and put the formatting implementation in `_formatters.py`.

```python
def format_solvent_shell(result, *, output_format: str = "text") -> str:
    """Format solvent-shell comparison output for the CLI."""

    if output_format == "json":
        return result.model_dump_json(indent=2)
    return "Solvent-shell comparison complete."
```

Do not compute new scientific results in a formatter. It should only render data
that aggregation or comparison already produced.

## Migrate in small steps

1. **Start from a passing single-file plugin.** Run the focused plugin tests
   before the split so you know whether a later failure came from the refactor.
2. **Create the package directory.** Replace `solvent_shell.py` with
   `solvent_shell/__init__.py` and keep the public `Analysis` subclass name,
   `name`, and `Settings` unchanged.
3. **Move trajectory helpers first.** Put job functions, workers, collectors, and
   sidecar writes in `_mda.py`. Keep public imports from
   `polyzymd.analyses.base`, `polyzymd.analyses.mda`, and documented
   `polyzymd.analyses.shared` utilities.
4. **Move plot helpers next.** Put plot data loading and rendering in
   `_plotters.py`. Verify plots read artifacts and sidecars only.
5. **Move structured models only if needed.** Add `_results.py` when the plugin
   has a custom result contract. Otherwise, keep relying on canonical artifacts.
6. **Move formatting only if needed.** Add `_formatters.py` when CLI rendering is
   substantial enough to deserve its own file.
7. **Run focused tests and docs checks.** Use the plugin test file and the docs
   build before opening a pull request.

Example validation commands:

```bash
pixi run -e build pytest tests/analyses/plugins/test_solvent_shell.py -v
pixi run -e build make -C docs clean html
```

## Migration checklist

Before you consider the package split complete, check that:

- `__init__.py` still exposes exactly one public plugin class for discovery.
- The plugin's `name` and `Settings` are unchanged unless you intentionally made
  a user-facing configuration change.
- All imports from PolyzyMD use public facades or documented shared utilities.
- No contributor code imports from another plugin's private helper modules.
- Heavy dependencies such as MDAnalysis, NumPy-heavy routines, or matplotlib are
  imported lazily inside functions that need them when practical.
- Collectors return `ReplicateArtifact` objects and write large arrays or tables
  as registered sidecars.
- Aggregation and comparison consume artifacts or validated structured outputs,
  not raw runtime containers.
- `plot()` delegates to helpers that load cached artifacts and sidecars only.
- `_results.py` and `_formatters.py` exist only if they have clear current
  responsibilities.
- Focused plugin tests pass after the file move.

## Success state

You have converted a plugin into an advanced package when the public plugin class
is easier to read, private helpers have focused responsibilities, imports stay on
public PolyzyMD facades, heavy dependencies remain lazy, and plotting is still
artifact-only.
