"""Tests for framework cache identity utilities."""

from __future__ import annotations

from pathlib import Path

import pytest
from pydantic import BaseModel

from polyzymd.analyses._framework.cache_identity import (
    compute_cache_identity,
    extract_settings_fingerprint_from_path,
    settings_fingerprint,
    validate_settings_fingerprint,
)


class SimpleSettings(BaseModel):
    """Minimal settings model for fingerprint tests."""

    cutoff: float = 4.5
    selection: str = "protein"


class NestedSettings(BaseModel):
    """Settings model with nested dict/list structures."""

    grouping: dict[str, list[int]]
    metadata: dict[str, dict[str, float]]


class TestSettingsFingerprint:
    """Tests for settings_fingerprint canonicalization."""

    def test_identical_settings_produce_same_fingerprint(self):
        s1 = SimpleSettings(cutoff=4.5, selection="protein")
        s2 = SimpleSettings(cutoff=4.5, selection="protein")
        assert settings_fingerprint(s1) == settings_fingerprint(s2)

    def test_different_values_produce_different_fingerprint(self):
        s1 = SimpleSettings(cutoff=4.5)
        s2 = SimpleSettings(cutoff=5.0)
        assert settings_fingerprint(s1) != settings_fingerprint(s2)

    def test_fingerprint_is_8_hex_chars(self):
        s = SimpleSettings()
        fp = settings_fingerprint(s)
        assert len(fp) == 8
        assert all(c in "0123456789abcdef" for c in fp)

    def test_fingerprint_stable_across_construction_order(self):
        """Dict ordering should not matter with sorted key serialization."""
        s1 = SimpleSettings(cutoff=4.5, selection="protein")
        s2 = SimpleSettings(selection="protein", cutoff=4.5)
        assert settings_fingerprint(s1) == settings_fingerprint(s2)

    def test_fingerprint_stable_for_nested_structures_across_dict_order(self):
        """Nested dict/list content order should be canonicalized for fingerprints."""
        s1 = NestedSettings(
            grouping={"group_a": [1, 2], "group_b": [3, 4]},
            metadata={
                "weights": {"alpha": 0.2, "beta": 0.8},
                "limits": {"min": 0.0, "max": 1.0},
            },
        )
        s2 = NestedSettings(
            metadata={
                "limits": {"max": 1.0, "min": 0.0},
                "weights": {"beta": 0.8, "alpha": 0.2},
            },
            grouping={"group_b": [3, 4], "group_a": [1, 2]},
        )

        assert settings_fingerprint(s1) == settings_fingerprint(s2)

    def test_fingerprint_changes_when_nested_value_changes(self):
        """Changing nested values should change the settings fingerprint."""
        original = NestedSettings(
            grouping={"group_a": [1, 2], "group_b": [3, 4]},
            metadata={"weights": {"alpha": 0.2, "beta": 0.8}},
        )
        changed = NestedSettings(
            grouping={"group_a": [1, 2], "group_b": [3, 4]},
            metadata={"weights": {"alpha": 0.25, "beta": 0.75}},
        )

        assert settings_fingerprint(original) != settings_fingerprint(changed)


class TestCacheIdentity:
    """Tests for unified cache identity helper."""

    def test_compute_cache_identity_stable_for_same_inputs(self):
        """Cache identity should be deterministic for identical inputs."""
        settings = SimpleSettings(cutoff=4.5, selection="protein")
        a = compute_cache_identity(
            config_hash="abc123",
            settings=settings,
            cache_params={"equilibration": "10ns", "replicate": 1},
        )
        b = compute_cache_identity(
            config_hash="abc123",
            settings=settings,
            cache_params={"replicate": 1, "equilibration": "10ns"},
        )
        assert a == b
        assert len(a) == 12

    def test_compute_cache_identity_changes_when_settings_change(self):
        """Changing settings should change cache identity."""
        low = SimpleSettings(cutoff=4.0)
        high = SimpleSettings(cutoff=4.5)
        low_id = compute_cache_identity(config_hash="abc123", settings=low)
        high_id = compute_cache_identity(config_hash="abc123", settings=high)
        assert low_id != high_id

    def test_compute_cache_identity_requires_settings_or_fingerprint(self):
        """Helper should reject calls without settings identity input."""
        with pytest.raises(ValueError, match="Provide either settings or settings_fp"):
            compute_cache_identity(config_hash="abc123")


class TestSettingsFingerprintValidation:
    """Tests for settings fingerprint extraction and validation."""

    def test_extract_settings_fingerprint_from_path(self):
        """Extract helper should parse embedded settings fingerprints."""
        path = Path("/tmp/contacts_eq10ns_cut4.5_s1a2b3c4d_rep1.json")
        assert extract_settings_fingerprint_from_path(path) == "1a2b3c4d"

    def test_extract_settings_fingerprint_from_path_returns_none_when_absent(self):
        """Non-canonical paths without fingerprint should return None."""
        path = Path("/tmp/rmsf_eq10ns.json")
        assert extract_settings_fingerprint_from_path(path) is None

    def test_validate_settings_fingerprint_accepts_match(self):
        """Matching fingerprints should be accepted without warnings."""
        settings = SimpleSettings()
        current = settings_fingerprint(settings)
        assert validate_settings_fingerprint(current, settings, warn=False)

    def test_validate_settings_fingerprint_rejects_mismatch_with_warning(self):
        """Mismatched fingerprints should force recompute path."""
        settings = SimpleSettings(cutoff=4.5)
        with pytest.warns(UserWarning, match="Cached settings fingerprint mismatch"):
            valid = validate_settings_fingerprint("deadbeef", settings)
        assert valid is False

    def test_validate_settings_fingerprint_rejects_missing_by_default(self):
        """Missing fingerprint should be rejected by default."""
        settings = SimpleSettings(cutoff=4.5)
        with pytest.warns(UserWarning, match="rejecting cache"):
            valid = validate_settings_fingerprint(None, settings)
        assert valid is False

    def test_validate_settings_fingerprint_always_rejects_missing(self):
        """Missing fingerprints should be rejected without an opt-in path."""
        settings = SimpleSettings(cutoff=4.5)
        with pytest.warns(UserWarning, match="rejecting cache"):
            valid = validate_settings_fingerprint(None, settings)
        assert valid is False


class TestPluginFingerprintAgreement:
    """A plugin's own cache tag must match what the framework stamps on artifacts."""

    @staticmethod
    def _analyses_with_a_private_cache_tag() -> list[tuple[str, object]]:
        """Return registered analyses that compute their own settings cache tag."""

        from polyzymd.analyses.discovery import get_analysis, list_analyses

        found = []
        for name in list_analyses():
            analysis = get_analysis(name)()
            if hasattr(analysis, "_make_settings_cache_tag"):
                found.append((name, analysis))
        return found

    def test_private_cache_tag_matches_aggregate_fingerprint(self):
        """Every plugin that folds extra identity in must also override the hook.

        ``_stamp_replicate_identity`` writes ``aggregate_settings_fingerprint``
        onto a replicate artifact, while aggregation compares against the
        plugin's own tag. A plugin that folds a version into one and not the
        other rejects every artifact it just wrote.
        """

        checked = 0
        for name, analysis in self._analyses_with_a_private_cache_tag():
            settings = analysis.Settings()
            assert analysis.aggregate_settings_fingerprint(
                settings
            ) == analysis._make_settings_cache_tag(
                settings
            ), f"{name}: aggregate_settings_fingerprint and _make_settings_cache_tag disagree"
            checked += 1
        assert checked > 0

    def test_stamped_rmsd_artifact_passes_aggregation(self, tmp_path: Path):
        """An RMSD artifact stamped by the lifecycle is accepted by aggregation."""

        from polyzymd.analyses._framework.lifecycle import _stamp_replicate_identity
        from polyzymd.analyses.mda import ReplicateArtifact
        from polyzymd.analyses.rmsd import RMSDAnalysis, RMSDSettings
        from polyzymd.analyses.rmsd._mda import _validate_and_order_artifacts
        from polyzymd.analyses.shared.autocorrelation import AUTOCORRELATION_ESTIMATOR_VERSION

        del tmp_path
        analysis = RMSDAnalysis()
        settings = RMSDSettings()
        artifact = ReplicateArtifact(
            analysis_name="rmsd",
            condition_label="Control",
            replicate=1,
            payload={"runs": [{"label": run.label} for run in settings.runs]},
            metadata={
                "autocorrelation_estimator_version": AUTOCORRELATION_ESTIMATOR_VERSION,
            },
        )
        _stamp_replicate_identity(artifact, analysis, settings, "10ns")

        ordered = _validate_and_order_artifacts(
            condition_label="Control",
            expected_replicates=[1],
            run_labels=[run.label for run in settings.runs],
            settings_fingerprint=analysis._make_settings_cache_tag(settings),
            artifacts=[artifact],
        )
        assert [a.replicate for a in ordered] == [1]
