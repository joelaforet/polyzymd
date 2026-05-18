# Store arrays and tables with sidecars

This how-to is for contributors whose analysis produces arrays, per-frame
profiles, event tables, or other data that should not be stored directly in an
artifact JSON payload.

Use sidecars when you need durable data for aggregation, comparison, or plots,
but the data is too large or too structured for the small JSON summary in a
`ReplicateArtifact` or `ConditionArtifact`.

## Decide what goes in the payload and what goes in a sidecar

Keep the artifact payload small and JSON-compatible. Put enough information in
the payload to identify and summarize the result; put bulky data in sidecars and
register those sidecars on the artifact.

| Data | Put it in | Why |
| --- | --- | --- |
| Scalar metrics used by default aggregation | `payload["metrics"]` | Aggregation can combine finite replicate scalars directly. |
| Counts, labels, selections, dimensions, and compact summaries | `payload` or `metadata` | These make the artifact readable without opening large files. |
| Per-frame arrays, residue-by-frame matrices, distance matrices | NPZ sidecar | Arrays stay typed and compact, with a stored hash and size. |
| Event rows, contact tables, or per-frame tabular records | CSV or other table sidecar | Tables stay streamable and do not bloat JSON. |
| Raw MDAnalysis `Results` objects | Neither | Map them to JSON primitives and sidecars before constructing the artifact. |

Sidecar paths are artifact-store-relative paths such as
`sidecars/pair_distances.npz`. Do not store absolute paths in artifact payloads.
Use the `ArtifactSidecarRef` returned by `ArtifactStore`; it records the
relative path, SHA-256 digest, size, media type, and metadata.

## Write an NPZ array sidecar from a collector

In a collector, use the `artifact_store` provided by `MDACollectorContext`. The
store writes the file under the replicate artifact directory and returns a
registered sidecar reference.

```python
import numpy as np

from polyzymd.analyses.mda import (
    MDACollectorContext,
    MDAJobResult,
    ReplicateArtifact,
    frame_selection_payload,
    strict_json_payload,
)


class DistanceMatrixCollector:
    """Collect one completed distance job into a sidecar-backed artifact."""

    def __call__(
        self,
        ctx: MDACollectorContext,
        completed_jobs: list[MDAJobResult],
    ) -> ReplicateArtifact:
        if len(completed_jobs) != 1:
            raise ValueError(f"Expected one job, got {len(completed_jobs)}")

        job = completed_jobs[0]
        distances = np.asarray(job.results["distances_nm"], dtype=np.float64)
        frames = np.asarray(job.results["frames"], dtype=np.int64)

        sidecar = ctx.artifact_store.write_npz_sidecar(
            "sidecars/distance_matrix.npz",
            distances_nm=distances,
            frames=frames,
            metadata={
                "kind": "distance_matrix",
                "layout": "pair_x_frame",
                "shape": list(distances.shape),
            },
        )

        return ReplicateArtifact(
            analysis_name=ctx.analysis_name,
            condition_label=ctx.condition_label,
            replicate=ctx.replicate,
            payload={
                "n_pairs": int(distances.shape[0]),
                "n_frames": int(distances.shape[1]),
                "mean_distance_nm": float(np.nanmean(distances)),
                "metrics": {"mean_distance_nm": float(np.nanmean(distances))},
            },
            sidecars=[sidecar],
            provenance={
                "source": "distance_matrix_job",
                "frame_selection": frame_selection_payload(ctx.frame_selection),
                "universe_policy": strict_json_payload(
                    ctx.universe_policy.as_dict(),
                    analysis_name=ctx.analysis_name,
                ),
            },
            metadata={"result_kind": "distance_matrix_replicate"},
            warnings=list(ctx.warnings),
        )
```

The payload contains dimensions, a scalar summary, and optional default
aggregation metrics. The complete per-pair, per-frame array stays in the NPZ
sidecar.

## Write and register a CSV table sidecar

For tabular outputs, write the file under a store-relative `sidecars/` path and
then call `register_sidecar()`. The example below uses Python's standard `csv`
module, but the same pattern works for any writer that creates a file under the
artifact store root.

```python
import csv

from polyzymd.analyses.mda import MDACollectorContext, ReplicateArtifact


def write_contact_events(
    ctx: MDACollectorContext,
    rows: list[dict[str, int | str | float]],
) -> ReplicateArtifact:
    sidecar_path = "sidecars/contact_events.csv"
    csv_path = ctx.artifact_store.resolve_sidecar(sidecar_path)
    csv_path.parent.mkdir(parents=True, exist_ok=True)

    with csv_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["frame", "protein_resid", "polymer_atom", "distance_nm"],
        )
        writer.writeheader()
        writer.writerows(rows)

    sidecar = ctx.artifact_store.register_sidecar(
        sidecar_path,
        media_type="text/csv",
        metadata={
            "kind": "contact_events",
            "n_rows": len(rows),
            "columns": ["frame", "protein_resid", "polymer_atom", "distance_nm"],
        },
    )

    return ReplicateArtifact(
        analysis_name=ctx.analysis_name,
        condition_label=ctx.condition_label,
        replicate=ctx.replicate,
        payload={
            "n_contact_events": len(rows),
            "event_table_sidecar": sidecar.path,
        },
        sidecars=[sidecar],
        provenance={"source": "contact_event_collector"},
        metadata={"result_kind": "contact_events_replicate"},
        warnings=list(ctx.warnings),
    )
```

Registering the sidecar is the important step. It makes the CSV part of the
artifact contract instead of an arbitrary file that later code has to discover.

## Load sidecars through the artifact reference

Aggregation, comparison, and plotting code should load sidecars from the
artifact record. Do not scan directories for files that "look right".

```python
import csv

from polyzymd.analyses.mda import ArtifactSidecarRef, ArtifactStore, ReplicateArtifact


def sidecar_by_kind(artifact: ReplicateArtifact, kind: str) -> ArtifactSidecarRef:
    matches = [ref for ref in artifact.sidecars if ref.metadata.get("kind") == kind]
    if len(matches) != 1:
        raise ValueError(f"Expected one {kind!r} sidecar, found {len(matches)}")
    return matches[0]


def load_distance_matrix(run_dir):
    store = ArtifactStore(run_dir)
    artifact = store.read_replicate_result()
    sidecar = sidecar_by_kind(artifact, "distance_matrix")

    with store.load_npz_sidecar(sidecar) as data:
        distances_nm = data["distances_nm"].copy()
        frames = data["frames"].copy()

    return artifact, distances_nm, frames


def load_contact_events(run_dir):
    store = ArtifactStore(run_dir)
    artifact = store.read_replicate_result()
    sidecar = sidecar_by_kind(artifact, "contact_events")
    csv_path = store.validate_sidecar(sidecar)

    with csv_path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))

    return artifact, rows
```

`load_npz_sidecar()` and `validate_sidecar()` check the stored path, file size,
and SHA-256 digest before returning data. If the sidecar is missing or stale,
the store raises an artifact-store error instead of silently loading the wrong
file.

## Use sidecars for artifact-only plotting

`plot()` must read cached artifacts and registered sidecars only. It should not
reload trajectories, re-run MDAnalysis jobs, or infer outputs by walking the
filesystem.

A typical plot method should:

1. Locate each condition's artifact directory from the `PlotContext` data the
   orchestrator provides.
2. Read `result.json` with `ArtifactStore.read_replicate_result()`,
   `read_condition_result()`, or a plugin helper built on those methods.
3. Select the needed `ArtifactSidecarRef` from `artifact.sidecars` by metadata or
   a documented payload key.
4. Validate and load the sidecar with `ArtifactStore`.
5. Render figures from those cached values.

This keeps plotting reproducible: the same artifact and sidecar bytes that were
aggregated or compared are the bytes used for the figure.

## Avoid these anti-patterns

```python
# Wrong: huge JSON payloads make artifacts slow to read, diff, and validate.
ReplicateArtifact(
    analysis_name="my_analysis",
    condition_label="Polymer",
    replicate=1,
    payload={"all_frame_distances": distances.tolist()},
)
```

```python
# Wrong: raw MDAnalysis Results objects are runtime containers, not artifacts.
ReplicateArtifact(
    analysis_name="my_analysis",
    condition_label="Polymer",
    replicate=1,
    payload={"results": job.results},
)
```

```python
# Wrong: arbitrary file discovery bypasses artifact validation and provenance.
for path in run_dir.glob("*.csv"):
    rows = path.read_text()
```

Use a compact payload plus registered sidecars instead: the artifact says which
files belong to it, and `ArtifactStore` validates those files before loading.

## Success state

You have implemented sidecar handling when:

- bulky arrays or tables are written under a store-relative `sidecars/` path;
- every sidecar is registered and included in `artifact.sidecars`;
- the JSON payload contains summaries, dimensions, labels, and sidecar keys, not
  large frame-by-frame data;
- aggregation and plotting load sidecars through `ArtifactStore` and the
  artifact's sidecar references; and
- no contributor code imports private framework modules.
