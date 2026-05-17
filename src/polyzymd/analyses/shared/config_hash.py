"""Compatibility facade for framework-owned cache identity helpers.

The implementation moved to :mod:`polyzymd.analyses._framework.cache_identity`.
This module remains temporarily available for existing imports.
"""

from polyzymd.analyses._framework.cache_identity import (
    _WARNED_CONFIG_HASH_MISMATCHES,
    SETTINGS_FINGERPRINT_PATTERN,
    compute_cache_identity,
    compute_config_hash,
    extract_settings_fingerprint_from_path,
    settings_fingerprint,
    validate_config_hash,
    validate_settings_fingerprint,
)

__all__ = [
    "SETTINGS_FINGERPRINT_PATTERN",
    "_WARNED_CONFIG_HASH_MISMATCHES",
    "compute_cache_identity",
    "compute_config_hash",
    "extract_settings_fingerprint_from_path",
    "settings_fingerprint",
    "validate_config_hash",
    "validate_settings_fingerprint",
]
