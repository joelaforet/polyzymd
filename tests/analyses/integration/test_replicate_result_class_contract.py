"""Regression tests for canonical MDA replicate artifact contracts."""

from __future__ import annotations

from pathlib import Path

import pytest

from polyzymd.analyses.mda import ArtifactStore, ReplicateArtifact


class TestReplicateArtifactContract:
    """Replicate artifacts round-trip through the artifact store."""

    @pytest.mark.parametrize("plugin_name", ["artifact_contract"])
    def test_replicate_artifact_roundtrip(self, plugin_name: str, tmp_path: Path) -> None:
        """Canonical replicate artifacts should roundtrip through the artifact store."""

        artifact = ReplicateArtifact(
            analysis_name=plugin_name,
            condition_label="Artifact Contract",
            replicate=1,
            payload={"metrics": {"smoke": 1.0}, "replicate_metrics": {"smoke": 1.0}},
            provenance={"source": "contract_test"},
            metadata={"settings_fingerprint": "contract"},
        )

        store = ArtifactStore(tmp_path)
        store.write_replicate_result(artifact)
        loaded = store.read_replicate_result()

        assert loaded == artifact
        assert loaded.analysis_name == plugin_name
        assert loaded.replicate == 1
