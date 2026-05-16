"""Tests for MDAnalysis extension-layer artifact models and store."""

from __future__ import annotations

from pathlib import Path

import pytest
from pydantic import ValidationError

from polyzymd.analyses.mda.artifacts import (
    MDA_ARTIFACT_SCHEMA_VERSION,
    ArtifactManifest,
    ArtifactSidecarRef,
    ComparisonArtifact,
    ConditionArtifact,
    ReplicateArtifact,
    raw_mdanalysis_results_path,
)
from polyzymd.analyses.mda.store import ArtifactStore, ArtifactStoreError

VALID_SHA256 = "0" * 64


class FakeMDAnalysisResults(dict):
    """Import-light stand-in for ``MDAnalysis.analysis.results.Results``."""

    __module__ = "MDAnalysis.analysis.results"


def test_replicate_artifact_json_roundtrip() -> None:
    """Replicate envelopes should roundtrip through JSON with sidecar refs."""

    sidecar = ArtifactSidecarRef(
        path="arrays/rmsd.npz",
        sha256=VALID_SHA256,
        size_bytes=12,
        media_type="application/x-npz",
        metadata={"array": "rmsd"},
    )
    artifact = ReplicateArtifact(
        analysis_name="rmsd",
        condition_label="PEG",
        replicate=2,
        payload={"mean_rmsd": 1.25},
        sidecars=[sidecar],
        custom_future_field={"kept": True},
    )

    loaded = ReplicateArtifact.model_validate_json(artifact.model_dump_json())

    assert loaded.schema_version == MDA_ARTIFACT_SCHEMA_VERSION
    assert loaded.artifact_type == "replicate"
    assert loaded.analysis_name == "rmsd"
    assert loaded.condition_label == "PEG"
    assert loaded.replicate == 2
    assert loaded.payload == {"mean_rmsd": 1.25}
    assert loaded.sidecars == [sidecar]
    assert loaded.model_extra == {"custom_future_field": {"kept": True}}
    assert loaded.model_dump()["artifact_type"] == "replicate"


def test_condition_and_comparison_artifact_roundtrip() -> None:
    """Condition and comparison envelopes should keep generic payload dictionaries."""

    condition = ConditionArtifact(
        analysis_name="rmsd",
        condition_label="PEG",
        replicates=[1, 2],
        payload={"replicate_values": [1.0, 1.2]},
    )
    comparison = ComparisonArtifact(
        analysis_name="rmsd",
        conditions=["Control", "PEG"],
        control_label="Control",
        effective_control="Control",
        payload={"pairwise": [{"condition": "PEG", "p_value": 0.2}]},
    )

    assert ConditionArtifact.model_validate_json(condition.model_dump_json()) == condition
    assert ComparisonArtifact.model_validate_json(comparison.model_dump_json()) == comparison


@pytest.mark.parametrize(
    ("field", "kwargs"),
    [
        ("payload", {"payload": {"raw": FakeMDAnalysisResults(value=1)}}),
        ("provenance", {"provenance": {"raw": FakeMDAnalysisResults(value=1)}}),
        ("metadata", {"metadata": {"raw": FakeMDAnalysisResults(value=1)}}),
    ],
)
def test_replicate_artifacts_reject_raw_mdanalysis_results(
    field: str,
    kwargs: dict[str, object],
) -> None:
    """Artifact dictionaries should not accept raw MDAnalysis Results."""

    with pytest.raises(ValidationError) as exc_info:
        ReplicateArtifact(analysis_name="rmsd", condition_label="PEG", replicate=1, **kwargs)
    message = str(exc_info.value)
    assert field in message
    assert "raw MDAnalysis Results" in message


@pytest.mark.parametrize(
    "kwargs",
    [
        {"payload": {"raw": FakeMDAnalysisResults(value=1)}},
        {"provenance": {"raw": FakeMDAnalysisResults(value=1)}},
        {"metadata": {"raw": FakeMDAnalysisResults(value=1)}},
    ],
)
def test_replicate_artifact_raw_results_error_mentions_results(kwargs: dict[str, object]) -> None:
    """Artifact raw-results errors should explain the mapping requirement."""

    with pytest.raises(ValidationError, match="raw MDAnalysis Results"):
        ReplicateArtifact(analysis_name="rmsd", condition_label="PEG", replicate=1, **kwargs)


def test_sidecar_and_manifest_metadata_reject_raw_mdanalysis_results() -> None:
    """Raw MDAnalysis Results should be rejected outside payload fields too."""

    with pytest.raises(ValidationError, match="raw MDAnalysis Results"):
        ArtifactSidecarRef(
            path="arrays/raw.npz",
            sha256=VALID_SHA256,
            size_bytes=0,
            metadata={"raw": FakeMDAnalysisResults(value=1)},
        )
    with pytest.raises(ValidationError, match="raw MDAnalysis Results"):
        ArtifactManifest(analysis_name="rmsd", provenance={"raw": FakeMDAnalysisResults(value=1)})


def test_artifact_extra_fields_reject_raw_mdanalysis_results() -> None:
    """Allowed future extension fields should still reject raw Results."""

    with pytest.raises(ValidationError, match="raw MDAnalysis Results"):
        ReplicateArtifact(
            analysis_name="rmsd",
            condition_label="PEG",
            replicate=1,
            raw=FakeMDAnalysisResults(value=1),
        )
    with pytest.raises(ValidationError, match="raw MDAnalysis Results"):
        ConditionArtifact(
            analysis_name="rmsd",
            condition_label="PEG",
            raw=FakeMDAnalysisResults(value=1),
        )
    with pytest.raises(ValidationError, match="raw MDAnalysis Results"):
        ComparisonArtifact(analysis_name="rmsd", raw=FakeMDAnalysisResults(value=1))


def test_sidecar_and_manifest_extra_fields_reject_raw_mdanalysis_results() -> None:
    """Sidecar and manifest extension fields should reject raw Results."""

    with pytest.raises(ValidationError, match="raw MDAnalysis Results"):
        ArtifactSidecarRef(
            path="arrays/raw.npz",
            sha256=VALID_SHA256,
            size_bytes=0,
            raw=FakeMDAnalysisResults(value=1),
        )
    with pytest.raises(ValidationError, match="raw MDAnalysis Results"):
        ArtifactManifest(analysis_name="rmsd", raw=FakeMDAnalysisResults(value=1))


def test_raw_mdanalysis_results_path_reports_nested_location() -> None:
    """The raw-results helper should report a useful nested path."""

    path = raw_mdanalysis_results_path({"outer": [{"inner": FakeMDAnalysisResults()}]})

    assert path == "$.outer[0].inner"


def test_raw_mdanalysis_results_path_detects_model_extra_bypass() -> None:
    """The raw-results helper should inspect Pydantic model extras."""

    artifact = ReplicateArtifact.model_construct(
        analysis_name="rmsd",
        condition_label="PEG",
        replicate=1,
        raw=FakeMDAnalysisResults(value=1),
    )

    assert raw_mdanalysis_results_path(artifact) == "$.raw"


@pytest.mark.parametrize("bad_path", ["/abs/file.npz", "../file.npz", "safe/../file.npz"])
def test_sidecar_refs_reject_absolute_and_parent_paths(bad_path: str) -> None:
    """Sidecar references should never store absolute or parent-traversal paths."""

    with pytest.raises(ValidationError):
        ArtifactSidecarRef(path=bad_path, sha256=VALID_SHA256, size_bytes=0)


@pytest.mark.parametrize("bad_sha", ["", "abc", "g" * 64, "1" * 63, "1" * 65])
def test_sidecar_refs_require_sha256(bad_sha: str) -> None:
    """Sidecar references should require SHA-256 hex digests."""

    with pytest.raises(ValidationError):
        ArtifactSidecarRef(path="arrays/data.npz", sha256=bad_sha, size_bytes=0)


def test_store_writes_and_reads_replicate_result_and_manifest(tmp_path: Path) -> None:
    """ArtifactStore should roundtrip replicate result and manifest JSON."""

    store = ArtifactStore(tmp_path)
    artifact = ReplicateArtifact(
        analysis_name="rg",
        condition_label="Polymer",
        replicate=1,
        payload={"mean_rg": 2.5},
    )
    manifest = ArtifactManifest(
        analysis_name="rg",
        artifact_id="rg-Polymer-1",
        inputs={"trajectory": "traj.xtc"},
    )

    result_path = store.write_replicate_result(artifact)
    manifest_path = store.write_manifest(manifest)

    assert result_path == tmp_path / "result.json"
    assert manifest_path == tmp_path / "manifest.json"
    assert store.read_replicate_result() == artifact
    assert store.read_manifest() == manifest


@pytest.mark.parametrize(
    "payload",
    ["{not json", '{"artifact_type": "replicate", "analysis_name": "rmsd"}'],
)
def test_store_read_replicate_result_wraps_invalid_json_and_schema(
    tmp_path: Path,
    payload: str,
) -> None:
    """ArtifactStore should include artifact paths for invalid JSON or schema."""

    store = ArtifactStore(tmp_path)
    artifact_path = tmp_path / "result.json"
    artifact_path.write_text(payload)

    with pytest.raises(
        ArtifactStoreError, match=r"Failed to validate replicate artifact .*result\.json"
    ):
        store.read_replicate_result()


def test_store_writes_npz_sidecar_lazily_and_validates(tmp_path: Path) -> None:
    """NPZ writing should lazy-import NumPy and produce a valid relative sidecar ref."""

    import numpy as np

    store = ArtifactStore(tmp_path)

    ref = store.write_npz_sidecar("arrays/values.npz", values=np.array([1.0, 2.0, 3.0]))

    assert ref.path == "arrays/values.npz"
    assert ref.size_bytes > 0
    assert len(ref.sha256) == 64
    assert store.resolve_sidecar(ref) == tmp_path / "arrays" / "values.npz"
    assert store.validate_sidecar(ref) == tmp_path / "arrays" / "values.npz"
    with np.load(store.resolve_sidecar(ref)) as data:
        assert data["values"].tolist() == [1.0, 2.0, 3.0]


def test_register_sidecar_records_hash_and_size(tmp_path: Path) -> None:
    """Registering an existing sidecar should stream hash and size metadata."""

    store = ArtifactStore(tmp_path)
    sidecar_path = tmp_path / "events" / "raw.bin"
    sidecar_path.parent.mkdir(parents=True)
    sidecar_path.write_bytes(b"abcdef")

    ref = store.register_sidecar("events/raw.bin", media_type="application/octet-stream")

    assert ref.path == "events/raw.bin"
    assert ref.size_bytes == 6
    assert ref.media_type == "application/octet-stream"
    assert store.validate_sidecar(ref) == sidecar_path


@pytest.mark.parametrize(
    "bad_path", ["/absolute/result.json", "../result.json", "nested/../x.json"]
)
def test_store_rejects_absolute_and_parent_paths(tmp_path: Path, bad_path: str) -> None:
    """Store operations should reject unsafe JSON and sidecar paths."""

    store = ArtifactStore(tmp_path)
    artifact = ReplicateArtifact(analysis_name="rmsd", condition_label="A", replicate=1)

    with pytest.raises(ArtifactStoreError):
        store.write_replicate_result(artifact, bad_path)
    with pytest.raises(ArtifactStoreError):
        store.resolve_sidecar(bad_path)


def test_register_sidecar_rejects_missing_file(tmp_path: Path) -> None:
    """Registering a missing sidecar should fail explicitly."""

    store = ArtifactStore(tmp_path)

    with pytest.raises(ArtifactStoreError, match="Sidecar does not exist"):
        store.register_sidecar("missing/data.npz")


def test_validate_sidecar_rejects_missing_file(tmp_path: Path) -> None:
    """Validating a missing sidecar reference should fail explicitly."""

    store = ArtifactStore(tmp_path)
    ref = ArtifactSidecarRef(path="missing/data.npz", sha256=VALID_SHA256, size_bytes=0)

    with pytest.raises(ArtifactStoreError, match="Missing sidecar"):
        store.validate_sidecar(ref)


def test_validate_sidecar_rejects_hash_mismatch(tmp_path: Path) -> None:
    """Sidecar validation should detect stale or corrupted bytes."""

    store = ArtifactStore(tmp_path)
    path = tmp_path / "arrays" / "values.bin"
    path.parent.mkdir(parents=True)
    path.write_bytes(b"old")
    ref = store.register_sidecar("arrays/values.bin")

    path.write_bytes(b"new")

    with pytest.raises(ArtifactStoreError, match="SHA-256 mismatch"):
        store.validate_sidecar(ref)


def test_validate_sidecar_rejects_size_mismatch(tmp_path: Path) -> None:
    """Sidecar validation should check byte size before hashing."""

    store = ArtifactStore(tmp_path)
    path = tmp_path / "arrays" / "values.bin"
    path.parent.mkdir(parents=True)
    path.write_bytes(b"old")
    ref = store.register_sidecar("arrays/values.bin")

    path.write_bytes(b"longer")

    with pytest.raises(ArtifactStoreError, match="size mismatch"):
        store.validate_sidecar(ref)
