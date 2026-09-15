"""Tests for framework cache identity utilities."""

from __future__ import annotations

from pathlib import Path

import pytest
from pydantic import BaseModel

from polyzymd.analyses._framework.cache_identity import (
    settings_fingerprint,
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
        other rejects every artifact it just wrote. Every built-in plugin now
        runs on the observable contract, where the runner owns the identity
        block, so the list this walks is empty; the check stays as a guard for
        a plugin that brings a private tag back.
        """

        for name, analysis in self._analyses_with_a_private_cache_tag():
            settings = analysis.Settings()
            assert analysis.aggregate_settings_fingerprint(
                settings
            ) == analysis._make_settings_cache_tag(
                settings
            ), f"{name}: aggregate_settings_fingerprint and _make_settings_cache_tag disagree"
