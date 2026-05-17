"""Compatibility smoke tests for shared cache identity imports."""

from __future__ import annotations

from polyzymd.analyses._framework import cache_identity
from polyzymd.analyses.shared import config_hash


def test_shared_config_hash_facade_reexports_framework_helpers() -> None:
    """Legacy shared imports should point at framework cache helpers."""
    assert config_hash.compute_config_hash is cache_identity.compute_config_hash
    assert config_hash.validate_config_hash is cache_identity.validate_config_hash
    assert config_hash.settings_fingerprint is cache_identity.settings_fingerprint
