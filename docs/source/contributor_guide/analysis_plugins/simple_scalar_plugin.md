# Build a simple scalar analysis plugin

This tutorial turns the scaffold pattern into one small trajectory-native
analysis. You will sketch a plugin that computes one scalar value per replicate,
stores it in a `ReplicateArtifact`, lets PolyzyMD aggregate replicate metrics
into a `ConditionArtifact`, and lets the MDAnalysis artifact comparison path
build the default scalar comparison from that artifact.

The metric in this tutorial is **teaching-only**: the mean x-coordinate of the
selected protein atoms. It is not validated science and should not be used to
interpret enzyme-polymer simulations. It exists only because it is small enough
to show the plugin lifecycle in one sitting.

## Before you start

You should already have completed {doc}`first_scaffold`. Keep using a scratch
plugin while you learn. The snippets below are meant to replace the scaffolded
placeholder pieces in a tutorial plugin such as `solvent_shell`; they are not a
request to add this metric to PolyzyMD.

You also need the architecture vocabulary from {doc}`architecture`:

- `build_mda_jobs()` creates `MDAAnalysisJob` objects.
- A collector maps completed jobs to a `ReplicateArtifact`.
- Default aggregation reads scalar replicate metrics from
  `payload["metrics"]` and writes a `ConditionArtifact`.
- Default comparison reads the canonical `ConditionArtifact` payload; use
  `extract_metrics()` only to customize metric descriptors for that payload.

## 1. Define a tiny settings model

Start with one atom selection and one constant metric key. Use lowercase
`chainid` in MDAnalysis selections. In the PolyzyMD chain convention,
`chainid A` is the protein.

```python
from typing import ClassVar

from pydantic import BaseModel, Field

from polyzymd.analyses.base import Analysis
from polyzymd.analyses.mda import (
    MDAAnalysisJob,
    MDACollectorContext,
    MDAJobResult,
    MDAReplicateJobContext,
    ReplicateArtifact,
    frame_selection_payload,
    strict_json_payload,
)


METRIC_NAME: str = "mean_protein_x_nm"


class MeanProteinXSettings(BaseModel):
    """Settings for the teaching-only scalar metric."""

    selection: str = Field(
        default="protein and chainid A",
        description="MDAnalysis selection for atoms to summarize.",
    )
```

This imports PolyzyMD symbols only from the public contributor facades:
`polyzymd.analyses.base` and `polyzymd.analyses.mda`. The metric key is a
module constant instead of a setting so every replicate and condition writes the
same scalar name.

## 2. Write a function job that returns JSON-compatible data

`MDAAnalysisJob.from_function()` adapts a plain function into the MDAnalysis job
shape. The function receives the loaded universe plus frame-selection keyword
arguments from PolyzyMD.

Return a JSON-compatible dictionary. Do not return or store raw MDAnalysis
`Results` objects in artifacts.

```python
def calculate_mean_x(
    universe,
    *,
    selection: str,
    start=None,
    stop=None,
    step=None,
    frames=None,
    **_frame_kwargs,
) -> dict[str, float | int | str]:
    """Calculate a teaching-only scalar for one replicate."""

    atoms = universe.select_atoms(selection)
    values: list[float] = []

    iterator = (
        universe.trajectory[frames]
        if frames is not None
        else universe.trajectory[start:stop:step]
    )
    for _ts in iterator:
        if len(atoms) == 0:
            raise ValueError(f"Selection matched no atoms: {selection!r}")
        # MDAnalysis positions are commonly in Angstrom; divide by 10 for nm.
        values.append(float(atoms.positions[:, 0].mean() / 10.0))

    mean_x = sum(values) / len(values) if values else 0.0
    return {
        "mean_x_nm": mean_x,
        "n_frames": len(values),
        "selection": selection,
    }
```

The important boundary is the return type. The function can use MDAnalysis while
it is running, but its output is already reduced to primitive values that can be
validated and serialized.

## 3. Build the function job

In your `Analysis` subclass, use the framework-provided context. Do not load the
configuration again, and do not pass a backend policy to `from_function()` for a
function-adapter job.

```python
class MeanProteinXAnalysis(Analysis):
    """Teaching-only scalar analysis plugin."""

    name: ClassVar[str] = "mean_protein_x"
    Settings: ClassVar[type[BaseModel]] = MeanProteinXSettings

    def build_mda_jobs(self, ctx: MDAReplicateJobContext) -> list[MDAAnalysisJob]:
        settings = ctx.settings
        return [
            MDAAnalysisJob.from_function(
                name="mean_protein_x",
                function=calculate_mean_x,
                universe=ctx.universe,
                frame_selection=ctx.frame_selection,
                universe_policy=ctx.universe_policy,
                function_kwargs={"selection": settings.selection},
            )
        ]
```

This job runs once for a replicate. The function handles the selected frame
iteration internally and returns a compact result dictionary.

## 4. Collect the job result into a ReplicateArtifact

The collector is where you translate completed job output into the durable
replicate contract. Put the scalar that should be aggregated under
`payload["metrics"]`.

```python
class MeanProteinXCollector:
    """Collect one completed function job into a replicate artifact."""

    def __call__(
        self,
        ctx: MDACollectorContext,
        completed_jobs: list[MDAJobResult],
    ) -> ReplicateArtifact:
        if len(completed_jobs) != 1:
            raise ValueError(f"Expected one job, got {len(completed_jobs)}")

        job = completed_jobs[0]
        result = dict(job.results)
        mean_x = float(result["mean_x_nm"])

        metadata = {"result_kind": "teaching_scalar"}
        if ctx.settings_fingerprint is not None:
            metadata["settings_fingerprint"] = ctx.settings_fingerprint

        return ReplicateArtifact(
            analysis_name=ctx.analysis_name,
            condition_label=ctx.condition_label,
            replicate=ctx.replicate,
            payload={
                "selection": result["selection"],
                "n_frames": int(result["n_frames"]),
                "metrics": {METRIC_NAME: mean_x},
            },
            provenance={
                "source": "mean_protein_x_teaching_function",
                "frame_selection": frame_selection_payload(ctx.frame_selection),
                "universe_policy": strict_json_payload(
                    ctx.universe_policy.as_dict(),
                    analysis_name=ctx.analysis_name,
                ),
            },
            metadata=metadata,
            warnings=list(ctx.warnings),
        )

```

Then add the collector hook to the `MeanProteinXAnalysis` class above:

```python
    def build_mda_collector(self, ctx: MDACollectorContext):
        del ctx
        return MeanProteinXCollector()
```

Because `payload["metrics"]` contains one finite scalar, the default aggregation
path can combine replicate artifacts without a custom `aggregate()` method.

## 5. Let default aggregation build the ConditionArtifact

For this simple scalar plugin, do not override `aggregate()`. PolyzyMD's default
MDAnalysis aggregation reads each replicate artifact's `payload["metrics"]` and
creates a `ConditionArtifact` with a payload shaped like this:

```python
{
    "metrics": {
        "mean_protein_x_nm": {
            "name": "mean_protein_x_nm",
            "values": [1.2, 1.3, 1.1],
            "mean": 1.2,
            "sem": 0.0577,
            "std": 0.1,
            "n": 3,
        }
    },
    "replicate_metrics": {
        "1": {"mean_protein_x_nm": 1.2},
        "2": {"mean_protein_x_nm": 1.3},
        "3": {"mean_protein_x_nm": 1.1},
    },
    "n_replicates": 3,
}
```

The exact numbers will differ. The stable idea is that aggregation computes
replicate values, mean, standard deviation, and SEM for each named scalar metric.

## 6. Let artifact comparison use the condition metrics

For this simple MDAnalysis artifact tutorial, the default comparison reads the
canonical `ConditionArtifact.payload["metrics"]`, builds the internal
`MetricValue` inputs from the stored means, SEMs, and replicate values, and
returns a comparison artifact. Implement `extract_metrics()` only when you need
to customize metric direction, labels, or units from the same canonical payload.

Because the teaching metric has no scientific direction, treat any comparison
output as a lifecycle check only. For a real plugin, define metric direction and
labels through the supported comparison contract only after the metric has a
validated interpretation.

## Success state

You have the pieces of a simple scalar plugin when:

- `build_mda_jobs()` returns one `MDAAnalysisJob.from_function()` job.
- The function returns a JSON-compatible dictionary, not raw MDAnalysis results.
- The collector returns a `ReplicateArtifact` with a finite scalar under
  `payload["metrics"]`.
- You do not override `aggregate()` because default aggregation can create the
  `ConditionArtifact` from replicate metrics.
- Default comparison reads the `ConditionArtifact` metrics directly, or your
  `extract_metrics()` reads that canonical payload to add metric metadata.

## Common mistakes

- **Treating the teaching metric as science.** The x-coordinate example is only a
  lifecycle exercise. Replace it with a validated quantity for real work.
- **Using the wrong chain-selection spelling.** Use lowercase `chainid`, for
  example `protein and chainid A`.
- **Passing backend settings to `MDAAnalysisJob.from_function()`.** Function jobs
  use the default adapter path; do not pass `backend_policy` in this tutorial.
- **Serializing raw MDAnalysis results.** Store primitive values in the artifact
  payload. Use {doc}`sidecars` later for arrays or tables.
- **Hiding aggregation inside the collector.** The collector should describe one
  replicate. Let the default aggregation stage combine replicates.
- **Importing private framework modules.** Contributor examples should use
  `polyzymd.analyses.base` and `polyzymd.analyses.mda`.

## What to read next

- {doc}`../extending_analyses` for the full implementation guide.
- {doc}`sidecars` for the next practical guide when your plugin needs array or
  table outputs.
- {doc}`../../api/analyses_base` for `Analysis` and `MetricValue` API details.
- {doc}`../../api/analyses_mda` for `MDAAnalysisJob`, collectors, and artifact
  models.
