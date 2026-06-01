# Extend PolyzyMD with MDAnalysis-native analyses

This guide shows the supported contributor workflow for adding a new analysis
plugin. PolyzyMD treats each trajectory-native analysis as an MDAnalysis-style
analysis, then lifts it from one trajectory to condition/replicate ensembles.
MDAnalysis owns the per-trajectory `Universe`, `AtomGroup`, frame iteration,
`AnalysisBase`, and `Results` layer. PolyzyMD owns discovery, cache provenance,
replicate artifacts, condition aggregation, cross-condition statistics, CLI
wiring, and plotting from cached artifacts.

Start with the scaffold unless you are updating a built-in analysis. The
generated code matches the current MDAnalysis-native extension layer and is the
smallest working example of the contract.

## What you will build

A contributor plugin normally provides:

1. a Pydantic settings model for user-facing options;
2. one or more `MDAAnalysisJob` objects in `build_mda_jobs()`;
3. a collector from completed MDAnalysis jobs to a `ReplicateArtifact`;
4. optional custom aggregation when `payload["metrics"]` is not enough;
5. default artifact comparison by implementing `extract_metrics()` over canonical
   `ConditionArtifact` payloads, or a custom `compare()`;
6. optional `plot()` code that reads artifacts and sidecars only.

The default path is intentionally small. A single-file plugin can wrap a
function with `MDAAnalysisJob.from_function()`. More advanced plugins can build
direct `AnalysisBase`-compatible objects in a package-level `_mda.py` helper.

## Start with the scaffold

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
`src/polyzymd/analyses/<name>/` both participate in discovery. No registry
edits, decorators, or bootstrap imports are needed.

The scaffold implementation still accepts the historical `style="measurement"`
token for programmatic compatibility. That token only means the default
single-file MDAnalysis-native scaffold; it does not refer to the removed
scalar-measurement API and is not exposed as a `--style` CLI choice.

Use an advanced package when your plugin needs a direct `AnalysisBase` helper,
sidecars, multiple metrics, custom aggregation, or custom plotting:

```bash
pixi run -e build polyzymd new-analysis solvent_shell --advanced
pixi run -e build polyzymd new-analysis solvent_shell --style dict
```

`--advanced` and `--style dict` generate MDAnalysis-native packages that store
dict metrics in framework-owned `ReplicateArtifact` and `ConditionArtifact`
files through `ArtifactStore`.

## Public imports

Contributor plugins should import from public facades only:

```python
from polyzymd.analyses.base import Analysis, MetricValue
from polyzymd.analyses.mda import (
    MDAAnalysisJob,
    MDACollectorContext,
    MDAReplicateJobContext,
    ReplicateArtifact,
)
```

`polyzymd.analyses.base` is the stable contributor surface for the `Analysis`
base class, lifecycle contexts, `MetricValue`, and comparison result models.
`polyzymd.analyses.mda` is the stable MDAnalysis extension-layer surface for
jobs, frame selection, artifacts, stores, aggregation, and comparison helpers.
Do not import private `_framework/` modules from contributor plugins.

## Minimal function-adapter plugin

Use `MDAAnalysisJob.from_function()` when one function can compute a
replicate-level result from an already-loaded MDAnalysis `Universe`. The
function receives MDAnalysis-style frame-selection kwargs and should return a
strict JSON-compatible object.

```python
from __future__ import annotations

from collections.abc import Sequence
from typing import Any, ClassVar

from pydantic import BaseModel, Field

from polyzymd.analyses.base import Analysis
from polyzymd.analyses.exceptions import PluginContractError
from polyzymd.analyses.mda import (
    MDAAnalysisJob,
    MDACollectorContext,
    MDAJobResult,
    MDAReplicateJobContext,
    ReplicateArtifact,
    frame_selection_payload,
    strict_json_payload,
)

METRIC_NAME = "mean_shell_count"


class SolventShellSettings(BaseModel):
    """Settings for solvent shell analysis."""

    selection: str = Field(default="protein and name CA", min_length=1)
    scale: float = Field(default=1.0, gt=0.0)


def measure_solvent_shell(
    universe: Any,
    *,
    settings: SolventShellSettings,
    start: int | None = None,
    stop: int | None = None,
    step: int | None = None,
    frames: Sequence[Any] | None = None,
) -> dict[str, Any]:
    """Return a JSON-compatible replicate result."""

    atoms = universe.select_atoms(settings.selection)
    n_frames = _selected_frame_count(universe, start=start, stop=stop, step=step, frames=frames)
    mean_shell_count = float(len(atoms)) * settings.scale
    return {"metrics": {METRIC_NAME: mean_shell_count}, "n_frames": n_frames}


class SolventShellArtifactCollector:
    """Collect completed jobs into one replicate artifact."""

    def __call__(
        self,
        ctx: MDACollectorContext,
        completed_jobs: Sequence[MDAJobResult],
    ) -> ReplicateArtifact:
        if len(completed_jobs) != 1:
            raise PluginContractError("solvent_shell expects exactly one MDA job.")
        job = completed_jobs[0]
        result_payload = strict_json_payload(job.results, analysis_name=ctx.analysis_name)
        metrics = result_payload.get("metrics")
        if not isinstance(metrics, dict):
            raise PluginContractError("Job results must include a metrics mapping.")
        metadata = {"result_kind": "solvent_shell_replicate"}
        if ctx.settings_fingerprint is not None:
            metadata["settings_fingerprint"] = ctx.settings_fingerprint
        return ReplicateArtifact(
            analysis_name=ctx.analysis_name,
            condition_label=ctx.condition_label,
            replicate=ctx.replicate,
            payload={"metrics": {name: float(value) for name, value in metrics.items()}},
            provenance={
                "source": "solvent_shell_function_adapter",
                "frame_selection": frame_selection_payload(ctx.frame_selection),
                "universe_policy": strict_json_payload(
                    ctx.universe_policy.as_dict(), analysis_name=ctx.analysis_name
                ),
            },
            metadata=metadata,
            warnings=list(ctx.warnings),
        )


class SolventShellAnalysis(Analysis):
    """Solvent shell analysis backed by an MDAnalysis-compatible job."""

    name: ClassVar[str] = "solvent_shell"
    Settings: ClassVar[type[BaseModel]] = SolventShellSettings

    def build_mda_jobs(self, ctx: MDAReplicateJobContext) -> Sequence[MDAAnalysisJob]:
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

    def build_mda_collector(self, ctx: MDACollectorContext) -> SolventShellArtifactCollector:
        del ctx
        return SolventShellArtifactCollector()
```

The helper `_selected_frame_count()` is generated by the scaffold and handles
both `start`/`stop`/`step` and explicit `frames` selectors. Replace the
placeholder calculation with your scientific logic.

Each replicate declares explicit scalar values in
`ReplicateArtifact.payload["metrics"]`; the default aggregator combines those
values into a canonical `ConditionArtifact` with per-metric `mean`, `std`, `sem`,
`n`, and replicate `values`. During comparison, `extract_metrics()` receives the
canonical condition artifact when a plugin needs to customize metric descriptors;
simple scalar artifact plugins can use the framework default directly.

## Use direct AnalysisBase-compatible jobs for frame algorithms

Use a direct MDAnalysis `AnalysisBase` subclass when the calculation needs to do
work on every frame, store arrays/events, or use MDAnalysis internal backends.
Put the helper in `_mda.py` for package plugins and keep heavy imports lazy:

```python
def build_solvent_shell_analysis(universe, *, settings):
    from MDAnalysis.analysis.base import AnalysisBase

    class SolventShellAnalysisBase(AnalysisBase):
        def __init__(self, universe, *, settings):
            self.universe = universe
            self.settings = settings
            self._atoms = universe.select_atoms(settings.selection)
            self._counts = []
            super().__init__(universe.trajectory)

        def _prepare(self):
            self._counts = []

        def _single_frame(self):
            self._counts.append(len(self._atoms))

        def _conclude(self):
            n_frames = len(self._counts)
            mean_shell_count = float(sum(self._counts) / n_frames) if n_frames else 0.0
            self.results.metrics = {"mean_shell_count": mean_shell_count}
            self.results.n_frames = n_frames

    return SolventShellAnalysisBase(universe, settings=settings)
```

Then build the job with an analysis factory:

```python
MDAAnalysisJob(
    name="solvent_shell",
    analysis_factory=lambda: build_solvent_shell_analysis(ctx.universe, settings=ctx.settings),
    frame_selection=ctx.frame_selection,
    universe_policy=ctx.universe_policy,
)
```

Non-MDAnalysis kernels are acceptable only when exposed through an
`AnalysisBase`-compatible object so PolyzyMD can keep one job/artifact lifecycle.

## FrameSelection and backend policy

PolyzyMD resolves each replicate's analysis window into a `FrameSelection`.
`FrameSelection.run_kwargs()` maps directly to MDAnalysis `run()` keyword
arguments:

| FrameSelection field | MDAnalysis meaning |
|----------------------|-------------------|
| `start` | first frame index |
| `stop` | exclusive final frame index |
| `step` | stride |
| `frames` | explicit integer frame list or boolean mask |

Do not mix `frames` with `start`/`stop`/`step`. If your analysis has reference
construction, autocorrelation, or variance-based subsampling requirements,
record the policy in artifact provenance so stale caches can be rejected.

MDAnalysis internal parallel backends are opt-in per job through
`MDABackendPolicy`. The default policy forwards no backend kwargs, because
PolyzyMD normally parallelizes over analyses, conditions, and replicates. Avoid
nested parallelism on HPC unless the scheduler configuration explicitly reserves
cores for each replicate job.

## Artifact contract

Collectors must map raw MDAnalysis outputs into PolyzyMD artifacts. The standard
objects are:

| Object | Scope | Contributor responsibility |
|--------|-------|----------------------------|
| `ReplicateArtifact` | one analysis on one replicate | include JSON-compatible payload, provenance, warnings, and sidecar refs |
| `ConditionArtifact` | aggregated replicates for one condition | use the default aggregator or a custom plugin reducer |
| `ComparisonArtifact` | cross-condition comparison | produced by the default MDA comparison path or custom `compare()` |
| `ArtifactStore` | filesystem persistence | write/read `result.json`, manifests, and sidecars safely |

For simple scalar plugins, put one replicate-level value per metric in
`payload["metrics"]`:

```python
payload={"metrics": {"mean_shell_count": 12.5}}
```

The default aggregator reads `payload["metrics"]` or
`payload["replicate_metrics"]`, validates finite scalar values, and produces a
condition artifact whose `payload["metrics"]` contains `mean`, `std`, `sem`,
`n`, and replicate `values` for each metric. It does not reduce arrays, event
tables, or frame-level values; those are scientific choices that belong in your
plugin's custom aggregation.

Large arrays, event tables, and profile data must be sidecars, not large JSON
fields. Use `ArtifactStore` to register sidecars so size and SHA-256 hashes are
validated on load.

## Raw MDAnalysis Results mapping rule

Never store raw MDAnalysis `Results` or `ResultsGroup` objects in artifact
payloads, provenance, metadata, sidecar refs, or comparison outputs. They are
runtime containers, not cache schemas. A collector must convert them first:

- scalar metrics to JSON numbers in `payload["metrics"]`;
- small summaries to JSON-compatible dictionaries/lists;
- arrays, profiles, and event tables to sidecars plus JSON metadata;
- labels, selections, frame policy, transformations, and software versions to
  provenance or metadata.

Artifact validation rejects raw MDAnalysis results recursively. This keeps cache
files stable across MDAnalysis versions and makes aggregation/comparison
independent of trajectory loading.

## Aggregation and comparison

Use the default aggregation path when each replicate artifact declares explicit
scalar metrics. The default MDA comparison path consumes `ConditionArtifact`
objects, validates that replicate IDs and aggregate statistics match the active
comparison, then delegates statistics to PolyzyMD's scalar comparison engine.

For the all-`ConditionArtifact` workflow, choose stable metric keys in each
`ReplicateArtifact.payload["metrics"]`. If the framework default comparison is
not expressive enough, implement `extract_metrics(summary)` to read the
canonical `ConditionArtifact` payload and return `MetricValue` objects with
thoughtful `higher_is_better`, unit, and direction metadata because those fields
drive ranking and CLI interpretation.

Override `aggregate()` when the replicate payload contains arrays or events that
need scientific reduction before comparison. Override `compare()` only when the
default scalar comparison cannot represent the output, such as per-residue
hypothesis families, multi-run tables, or custom statistical models.

## Plot from artifacts only

Plotting must not load trajectories, rebuild universes, rerun MDAnalysis jobs,
or scan non-artifact cache filenames. `plot(ctx)` should load canonical aggregate
or comparison artifacts and validated sidecars, then write figures to
`ctx.output_dir` and return their paths.

Small plugins can keep plotting inline. Extract `_plotters.py` when a plugin has
several figure types, sidecar loaders, or enough plotting code that lifecycle
wiring becomes hard to review.

## Trajectory selections and lazy imports

Heavy scientific dependencies such as MDAnalysis, OpenMM, OpenFF, MDTraj,
ParmEd, and PDBFixer must be imported lazily inside functions or methods. Do
not import them at module level in contributor plugins.

Selection strings are passed to MDAnalysis `Universe.select_atoms()` unless your
plugin documents an explicit wrapper. Follow the PolyzyMD chain convention:

| Chain | Role |
|-------|------|
| A | protein/enzyme |
| B | substrate/small molecule |
| C | polymer/conjugate |
| D+ | solvent, ions, and other molecules |

Common selection examples are:

- `protein`
- `protein and name CA`
- `chainid A`
- `chainid C`
- `resname SBM`
- `protein and (resid 77 or resid 156)`

Validate selections on representative topologies because available attributes
depend on topology format. GROMACS/GRO files may not preserve chain IDs; failures
should report the selection, condition, replicate, and topology source.

## Testing checklist

Generated MDAnalysis-native tests cover the scaffold contract:

- discovery and class variables;
- settings defaults and validation;
- `build_mda_jobs()` with fake universes or fake `AnalysisBase` objects;
- collector output as `ReplicateArtifact`;
- default artifact aggregation from explicit `payload["metrics"]`.

Production plugins should add tests for the behavior they customize:

- `FrameSelection` behavior, explicit `frames`, and backend policy when relevant;
- collector mapping from raw `results` to JSON payloads and sidecars;
- aggregation over replicate artifacts, including stale/missing sidecars;
- default `ConditionArtifact` comparison, `extract_metrics()` over canonical
  condition payloads, or custom comparison;
- plotting from artifacts without trajectory loading.

Run plugin tests through the pixi environment:

```bash
pixi run -e build pytest tests/analyses/plugins/test_<name>.py -v
```

## Further reading

PolyzyMD intentionally follows MDAnalysis idioms for trajectory-native work.
Before writing a complex plugin, review:

- MDAnalysis custom trajectory analysis tutorial:
  <https://userguide.mdanalysis.org/stable/examples/analysis/custom_trajectory_analysis.html>
- Michaud-Agrawal, N., Denning, E. J., Woolf, T. B., & Beckstein, O. (2011).
  MDAnalysis: A toolkit for the analysis of molecular dynamics simulations.
  *Journal of Computational Chemistry*, 32(10), 2319-2327.
- Gowers, R. J., Linke, M., Barnoud, J., Reddy, T. J. E., Melo, M. N., Seyler,
  S. L., Domański, J., Dotson, D. L., Buchoux, S., Kenney, I. M., & Beckstein,
  O. (2016). MDAnalysis: A Python package for the rapid analysis of molecular
  dynamics simulations. *Proceedings of the 15th Python in Science Conference*,
  98-105.

## Style checklist

- Use NumPy-style docstrings for new classes and methods.
- Keep imports ordered stdlib, third-party, local.
- Keep heavy scientific dependencies lazy.
- Use `X | None` annotations rather than `Optional[X]`.
- Run Ruff and Black checks on modified Python files.

```bash
pixi run -e build ruff check src/polyzymd/analyses/<name>.py tests/analyses/plugins/test_<name>.py
pixi run -e build black src/polyzymd/analyses/<name>.py tests/analyses/plugins/test_<name>.py --check
```
