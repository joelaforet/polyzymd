"""Tests for config hash utilities."""

from __future__ import annotations

from pydantic import BaseModel

from polyzymd.analyses.shared.config_hash import settings_fingerprint


class SimpleSettings(BaseModel):
    """Minimal settings model for fingerprint tests."""

    cutoff: float = 4.5
    selection: str = "protein"


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
