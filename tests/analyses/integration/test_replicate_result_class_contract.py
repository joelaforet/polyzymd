"""Regression tests for canonical MDA replicate artifact contracts."""

from __future__ import annotations

from pathlib import Path

import pytest

from polyzymd.analyses import get_analysis, list_analyses
from polyzymd.analyses.mda import ArtifactStore, ReplicateArtifact

MDA_ARTIFACT_PLUGINS = {
    "rmsf",
    "rmsd",
    "rg",
    "catalytic_triad",
    "distances",
    "secondary_structure",
    "contacts",
    "sasa",
    "hydrogen_bonds",
}


class TestReplicateArtifactContract:
    """Compute-stage plugins use canonical artifacts instead of typed result classes."""

    @pytest.mark.parametrize("plugin_name", sorted(MDA_ARTIFACT_PLUGINS))
    def test_migrated_plugin_replicate_result_class_is_none(self, plugin_name: str) -> None:
        """Migrated built-ins should not declare legacy replicate result classes."""

        cls = get_analysis(plugin_name)
        assert cls.ReplicateResultClass is None

    def test_all_compute_stage_plugins_are_mda_artifact_plugins(self) -> None:
        """All discovered compute plugins should be represented by the artifact contract."""

        discovered_compute_plugins = {
            name for name, cls in list_analyses().items() if cls.has_compute_stage
        }
        assert discovered_compute_plugins == MDA_ARTIFACT_PLUGINS

    @pytest.mark.parametrize("plugin_name", sorted(MDA_ARTIFACT_PLUGINS))
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
