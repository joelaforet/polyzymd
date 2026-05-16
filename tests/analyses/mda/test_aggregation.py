"""Tests for MDAnalysis replicate artifact aggregation."""

from __future__ import annotations

from pathlib import Path

import pytest

from polyzymd.analyses.mda.aggregation import (
    MDAAggregationContext,
    MDAAggregationError,
    aggregate_replicate_artifacts,
    aggregate_replicate_artifacts_from_disk,
)
from polyzymd.analyses.mda.artifacts import ArtifactSidecarRef, ReplicateArtifact
from polyzymd.analyses.mda.store import ArtifactStore


def _frame_selection() -> dict[str, int | None]:
    """Return deterministic frame-selection provenance for test artifacts."""

    return {"start": 0, "stop": 10, "step": 1, "n_frames_selected": 10}


def _artifact(
    replicate: int,
    metrics: dict[str, object] | None = None,
    *,
    analysis_name: str = "rmsd",
    condition_label: str = "PEG",
    settings_fingerprint: str = "abc123",
    frame_selection: dict[str, object] | None = None,
    sidecars: list[ArtifactSidecarRef] | None = None,
) -> ReplicateArtifact:
    """Build a replicate artifact with scalar metrics and provenance."""

    return ReplicateArtifact(
        analysis_name=analysis_name,
        condition_label=condition_label,
        replicate=replicate,
        payload={"metrics": metrics or {"mean_rmsd": 1.0 + replicate / 10.0}},
        provenance={"frame_selection": frame_selection or _frame_selection()},
        metadata={"settings_fingerprint": settings_fingerprint},
        sidecars=sidecars or [],
    )


def _context(
    replicates: tuple[int, ...] = (1, 2, 3),
    *,
    allow_partial: bool = False,
    min_replicates: int = 2,
) -> MDAAggregationContext:
    """Build a deterministic aggregation context."""

    return MDAAggregationContext(
        analysis_name="rmsd",
        condition_label="PEG",
        expected_replicates=replicates,
        settings_fingerprint="abc123",
        min_replicates=min_replicates,
        allow_partial=allow_partial,
    )


def test_aggregates_replicate_metrics_without_pooling_frames() -> None:
    """Aggregation should summarize one scalar per replicate per metric."""

    artifacts = [
        _artifact(1, {"mean_rmsd": 1.0}),
        _artifact(2, {"mean_rmsd": 1.2}),
        _artifact(3, {"mean_rmsd": 1.4}),
    ]

    condition = aggregate_replicate_artifacts(artifacts, _context())

    metric = condition.payload["metrics"]["mean_rmsd"]
    assert condition.artifact_type == "condition"
    assert condition.replicates == [1, 2, 3]
    assert condition.payload["replicate_metrics"] == {
        "1": {"mean_rmsd": 1.0},
        "2": {"mean_rmsd": 1.2},
        "3": {"mean_rmsd": 1.4},
    }
    assert metric["values"] == [1.0, 1.2, 1.4]
    assert metric["mean"] == pytest.approx(1.2)
    assert metric["std"] == pytest.approx(0.2)
    assert metric["sem"] == pytest.approx(0.2 / 3**0.5)
    assert condition.metadata["settings_fingerprint"] == "abc123"
    assert condition.provenance["frame_selection"] == _frame_selection()


def test_uses_replicate_metrics_payload_alias() -> None:
    """The default policy should accept the explicit replicate_metrics alias."""

    artifacts = [
        ReplicateArtifact(
            analysis_name="rmsd",
            condition_label="PEG",
            replicate=1,
            payload={"replicate_metrics": {"mean_rmsd": 1.0}},
            provenance={"frame_selection": _frame_selection()},
            metadata={"settings_fingerprint": "abc123"},
        ),
        ReplicateArtifact(
            analysis_name="rmsd",
            condition_label="PEG",
            replicate=2,
            payload={"replicate_metrics": {"mean_rmsd": 1.2}},
            provenance={"frame_selection": _frame_selection()},
            metadata={"settings_fingerprint": "abc123"},
        ),
    ]

    condition = aggregate_replicate_artifacts(artifacts, _context((1, 2)))

    assert condition.payload["metrics"]["mean_rmsd"]["n"] == 2


@pytest.mark.parametrize("bad_value", [[1.0, 2.0], {"event": 1}, "1.0", float("nan")])
def test_default_metric_policy_rejects_non_scalar_or_non_finite_metrics(bad_value: object) -> None:
    """Default aggregation should reject arrays, events, strings, and non-finite scalars."""

    artifacts = [_artifact(1, {"mean_rmsd": bad_value}), _artifact(2, {"mean_rmsd": 1.2})]

    with pytest.raises(MDAAggregationError, match="metric 'mean_rmsd'"):
        aggregate_replicate_artifacts(artifacts, _context((1, 2)))


def test_rejects_identity_and_settings_mismatches() -> None:
    """Aggregation should fail when artifact provenance does not match the context."""

    with pytest.raises(MDAAggregationError, match="analysis mismatch"):
        aggregate_replicate_artifacts(
            [_artifact(1, analysis_name="rg"), _artifact(2)],
            _context((1, 2)),
        )

    with pytest.raises(MDAAggregationError, match="settings fingerprint mismatch"):
        aggregate_replicate_artifacts(
            [_artifact(1, settings_fingerprint="old"), _artifact(2)],
            _context((1, 2)),
        )


def test_rejects_duplicate_unexpected_missing_and_incompatible_frame_selection() -> None:
    """Aggregation should enforce replicate identity and compatible frame provenance."""

    with pytest.raises(MDAAggregationError, match="Duplicate replicate"):
        aggregate_replicate_artifacts([_artifact(1), _artifact(1)], _context((1, 2)))
    with pytest.raises(MDAAggregationError, match="Unexpected replicate"):
        aggregate_replicate_artifacts([_artifact(1), _artifact(3)], _context((1, 2)))
    with pytest.raises(MDAAggregationError, match="missing replicate"):
        aggregate_replicate_artifacts([_artifact(1)], _context((1, 2)))
    with pytest.raises(MDAAggregationError, match="incompatible frame-selection"):
        aggregate_replicate_artifacts(
            [_artifact(1), _artifact(2, frame_selection={"start": 1, "stop": 10})],
            _context((1, 2)),
        )


def test_allow_partial_records_skipped_replicates_and_enforces_minimum() -> None:
    """Partial aggregation should be explicit and still require enough replicates."""

    partial = aggregate_replicate_artifacts(
        [_artifact(1), _artifact(3)],
        _context((1, 2, 3), allow_partial=True, min_replicates=2),
    )
    assert partial.replicates == [1, 3]
    assert partial.skipped_replicates == [{"replicate": 2, "reason": "missing artifact"}]

    with pytest.raises(MDAAggregationError, match="need at least 2"):
        aggregate_replicate_artifacts(
            [_artifact(1)],
            _context((1, 2, 3), allow_partial=True, min_replicates=2),
        )


def test_disk_helper_loads_artifacts_and_records_source_hashes(tmp_path: Path) -> None:
    """Disk aggregation should read artifacts through per-run stores only."""

    analysis_dir = tmp_path / "analysis" / "PEG" / "rmsd"
    for replicate in (1, 2):
        store = ArtifactStore(analysis_dir / f"run_{replicate}")
        store.write_replicate_result(_artifact(replicate))

    condition = aggregate_replicate_artifacts_from_disk(analysis_dir, _context((1, 2)))

    assert condition.replicates == [1, 2]
    assert [entry["replicate"] for entry in condition.source_replicates] == [1, 2]
    assert condition.source_replicates[0]["artifact"]["path"] == "result.json"
    assert len(condition.source_replicates[0]["artifact"]["sha256"]) == 64


def test_disk_helper_rejects_unexpected_discovered_artifacts(tmp_path: Path) -> None:
    """Disk aggregation should reject extra run_N artifacts outside the request."""

    analysis_dir = tmp_path / "analysis" / "PEG" / "rmsd"
    for replicate in (1, 2, 3):
        store = ArtifactStore(analysis_dir / f"run_{replicate}")
        store.write_replicate_result(_artifact(replicate))

    with pytest.raises(MDAAggregationError, match="unexpected replicate artifact"):
        aggregate_replicate_artifacts_from_disk(analysis_dir, _context((1, 2)))


def test_disk_helper_rejects_embedded_replicate_id_mismatch(tmp_path: Path) -> None:
    """Disk aggregation should reject artifacts whose payload identity lies."""

    analysis_dir = tmp_path / "analysis" / "PEG" / "rmsd"
    ArtifactStore(analysis_dir / "run_1").write_replicate_result(_artifact(1))
    ArtifactStore(analysis_dir / "run_2").write_replicate_result(_artifact(1))

    with pytest.raises(MDAAggregationError, match="embedded replicate ID mismatch"):
        aggregate_replicate_artifacts_from_disk(analysis_dir, _context((1, 2)))


def test_disk_helper_validates_sidecars_and_partial_missing(tmp_path: Path) -> None:
    """Disk aggregation should validate sidecars and record missing partial inputs."""

    analysis_dir = tmp_path / "analysis" / "PEG" / "rmsd"
    store = ArtifactStore(analysis_dir / "run_1")
    sidecar_path = store.root / "arrays" / "values.txt"
    sidecar_path.parent.mkdir(parents=True)
    sidecar_path.write_text("fresh")
    sidecar = store.register_sidecar("arrays/values.txt")
    store.write_replicate_result(_artifact(1, sidecars=[sidecar]))
    sidecar_path.write_text("stale")

    with pytest.raises(MDAAggregationError, match="stale sidecar"):
        aggregate_replicate_artifacts_from_disk(analysis_dir, _context((1,), min_replicates=1))

    sidecar_path.write_text("fresh")
    partial = aggregate_replicate_artifacts_from_disk(
        analysis_dir,
        _context((1, 2), allow_partial=True, min_replicates=1),
    )
    assert partial.skipped_replicates == [
        {
            "replicate": 2,
            "reason": "missing artifact",
            "path": str(analysis_dir / "run_2" / "result.json"),
        }
    ]


def test_disk_helper_allow_partial_records_stale_sidecars(tmp_path: Path) -> None:
    """Partial disk aggregation should skip stale sidecars and enforce minimum replicates."""

    analysis_dir = tmp_path / "analysis" / "PEG" / "rmsd"
    stale_store = ArtifactStore(analysis_dir / "run_1")
    sidecar_path = stale_store.root / "arrays" / "values.txt"
    sidecar_path.parent.mkdir(parents=True)
    sidecar_path.write_text("fresh")
    sidecar = stale_store.register_sidecar("arrays/values.txt")
    stale_store.write_replicate_result(_artifact(1, sidecars=[sidecar]))
    sidecar_path.write_text("stale")
    ArtifactStore(analysis_dir / "run_2").write_replicate_result(_artifact(2))

    partial = aggregate_replicate_artifacts_from_disk(
        analysis_dir,
        _context((1, 2), allow_partial=True, min_replicates=1),
    )

    assert partial.replicates == [2]
    assert partial.skipped_replicates[0]["replicate"] == 1
    assert partial.skipped_replicates[0]["path"] == str(analysis_dir / "run_1" / "result.json")
    assert "stale sidecar" in partial.skipped_replicates[0]["reason"]

    with pytest.raises(MDAAggregationError, match="need at least 2"):
        aggregate_replicate_artifacts_from_disk(
            analysis_dir,
            _context((1, 2), allow_partial=True, min_replicates=2),
        )
