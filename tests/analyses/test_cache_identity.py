"""Tests for framework cache identity utilities."""

from __future__ import annotations

from pathlib import Path

import pytest
from pydantic import BaseModel

from polyzymd.analyses._framework.cache_identity import (
    _WARNED_CONFIG_HASH_MISMATCHES,
    compute_cache_identity,
    extract_settings_fingerprint_from_path,
    settings_fingerprint,
    validate_config_hash,
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


class TestConfigHashValidation:
    """Tests for config hash mismatch diagnostics."""

    def test_config_hash_mismatch_warns_once_for_same_pair(self, monkeypatch):
        """Repeated mismatch checks should not spam identical warnings."""
        _WARNED_CONFIG_HASH_MISMATCHES.clear()
        config = object()
        monkeypatch.setattr(
            "polyzymd.analyses._framework.cache_identity.compute_config_hash",
            lambda current_config: "currenthash",
        )

        with pytest.warns(UserWarning, match="CONFIG HASH MISMATCH") as warnings_record:
            first = validate_config_hash("storedhash", config)
            second = validate_config_hash("storedhash", config)

        assert first is False
        assert second is False
        assert len(warnings_record) == 1
        assert "--recompute" in str(warnings_record[0].message)
        assert "clear the analysis cache" in str(warnings_record[0].message)
        _WARNED_CONFIG_HASH_MISMATCHES.clear()


class TestSettingsFingerprintValidation:
    """Tests for settings fingerprint extraction and validation."""

    def test_extract_settings_fingerprint_from_path(self):
        """Extract helper should parse embedded settings fingerprints."""
        path = Path("/tmp/contacts_eq10ns_cut4.5_s1a2b3c4d_rep1.json")
        assert extract_settings_fingerprint_from_path(path) == "1a2b3c4d"

    def test_extract_settings_fingerprint_from_path_returns_none_when_absent(self):
        """Legacy paths without fingerprint should return None."""
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

    def test_validate_settings_fingerprint_allows_legacy_missing_with_warning(self):
        """Missing fingerprint should remain backward-compatible."""
        settings = SimpleSettings(cutoff=4.5)
        with pytest.warns(UserWarning, match="missing settings fingerprint"):
            valid = validate_settings_fingerprint(None, settings)
        assert valid is True

    def test_validate_settings_fingerprint_can_reject_missing_strictly(self):
        """Strict callers should be able to reject missing fingerprints."""
        settings = SimpleSettings(cutoff=4.5)
        with pytest.warns(UserWarning, match="rejecting cache without strict validation"):
            valid = validate_settings_fingerprint(None, settings, allow_missing=False)
        assert valid is False
