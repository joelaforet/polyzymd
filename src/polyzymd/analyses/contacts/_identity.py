"""Contacts-domain settings identity helpers."""

from __future__ import annotations

import hashlib
import json
from typing import Any

from pydantic import BaseModel

CONTACTS_SETTINGS_FINGERPRINT_DOMAIN = "contacts-v2"
CONTACTS_DETECTION_FINGERPRINT_DOMAIN = "contacts-detection-v1"
CONTACTS_MDA_IMPLEMENTATION_VERSION = "contacts-mda-sparse-v1"
CONTACTS_GROUPING_CLASSIFIER_VERSION = "contacts-grouping-v1"
CONTACT_SEMANTICS_VERSION = "any_atom_residue_pair-v1"
CONTACTS_PBC_POLICY = "box=timestep.dimensions"


def contacts_settings_fingerprint(settings: BaseModel) -> str:
    """Compute the contacts-domain settings fingerprint.

    The contacts-domain identity includes only settings that alter the contact
    trajectory result itself. Downstream comparison-only and plotting settings
    are intentionally excluded so those workflows can change without
    invalidating contact caches.

    Parameters
    ----------
    settings : BaseModel
        Contacts-compatible settings model.

    Returns
    -------
    str
        First eight hexadecimal characters of the SHA-256 digest.
    """

    serialized = json.dumps(contacts_identity_payload(settings), sort_keys=True)
    return hashlib.sha256(serialized.encode("utf-8")).hexdigest()[:8]


def contacts_detection_fingerprint(settings: BaseModel) -> str:
    """Compute the sparse contact-detection settings fingerprint.

    Parameters
    ----------
    settings : BaseModel
        Contacts-compatible settings model.

    Returns
    -------
    str
        First eight hexadecimal characters of the SHA-256 digest.
    """

    serialized = json.dumps(contacts_detection_identity_payload(settings), sort_keys=True)
    return hashlib.sha256(serialized.encode("utf-8")).hexdigest()[:8]


def contacts_detection_identity_payload(settings: BaseModel) -> dict[str, Any]:
    """Return identity fields that alter sparse contact detection."""

    polymer_selection = _normalized_text(getattr(settings, "polymer_selection", "chainid C"))
    protein_selection = _normalized_text(getattr(settings, "protein_selection", "chainid A"))
    polymer_types = normalize_polymer_types(getattr(settings, "polymer_types", None))
    grouping = _normalized_text(getattr(settings, "grouping", "aa_class"))
    return {
        "domain": CONTACTS_DETECTION_FINGERPRINT_DOMAIN,
        "protein_selection": protein_selection,
        "polymer_selection": polymer_selection,
        "polymer_types": polymer_types,
        "effective_polymer_selection": effective_polymer_selection(
            polymer_selection,
            polymer_types,
        ),
        "cutoff": {"value": float(getattr(settings, "cutoff", 4.5)), "units": "angstrom"},
        "grouping": grouping,
        "grouping_classifier_version": CONTACTS_GROUPING_CLASSIFIER_VERSION,
        "pbc_policy": CONTACTS_PBC_POLICY,
        "contact_semantics": CONTACT_SEMANTICS_VERSION,
        "implementation_version": CONTACTS_MDA_IMPLEMENTATION_VERSION,
    }


def contacts_identity_payload(settings: BaseModel) -> dict[str, Any]:
    """Return the normalized settings payload used for contact cache identity.

    Parameters
    ----------
    settings : BaseModel
        Contacts-compatible settings model.

    Returns
    -------
    dict[str, Any]
        Canonical JSON-compatible identity payload.
    """

    polymer_selection = _normalized_text(getattr(settings, "polymer_selection", "chainid C"))
    protein_selection = _normalized_text(getattr(settings, "protein_selection", "chainid A"))
    polymer_types = normalize_polymer_types(getattr(settings, "polymer_types", None))

    return {
        "protein_selection": protein_selection,
        "polymer_selection": polymer_selection,
        "polymer_types": polymer_types,
        "effective_polymer_selection": effective_polymer_selection(
            polymer_selection,
            polymer_types,
        ),
        "cutoff": float(getattr(settings, "cutoff", 4.5)),
        "grouping": _normalized_text(getattr(settings, "grouping", "aa_class")),
        "compute_residence_times": bool(getattr(settings, "compute_residence_times", True)),
    }


def contacts_settings_fingerprint_candidates(settings: BaseModel) -> tuple[str, ...]:
    """Return current and legacy fingerprints accepted for contact caches.

    Parameters
    ----------
    settings : BaseModel
        Contacts-compatible settings model.

    Returns
    -------
    tuple[str, ...]
        Current contacts-domain fingerprint. When residence times are enabled,
        legacy fingerprints that predate the residence-time identity flag are
        appended for compatibility.
    """

    current = contacts_settings_fingerprint(settings)
    candidates = [current]
    if not bool(getattr(settings, "compute_residence_times", True)):
        return tuple(candidates)

    legacy_payload = dict(contacts_identity_payload(settings))
    legacy_payload.pop("compute_residence_times", None)
    legacy = hashlib.sha256(json.dumps(legacy_payload, sort_keys=True).encode("utf-8")).hexdigest()[
        :8
    ]
    if legacy not in candidates:
        candidates.append(legacy)

    from polyzymd.analyses.shared.config_hash import settings_fingerprint

    full_settings = settings_fingerprint(settings)
    if full_settings not in candidates:
        candidates.append(full_settings)

    return tuple(candidates)


def normalize_polymer_types(polymer_types: Any) -> list[str]:
    """Normalize the polymer type filter for stable cache identity.

    Parameters
    ----------
    polymer_types : Any
        Configured polymer type filter.

    Returns
    -------
    list[str]
        Sorted unique non-empty polymer residue names.
    """

    if polymer_types is None:
        return []
    return sorted(
        {str(polymer_type).strip() for polymer_type in polymer_types if str(polymer_type).strip()}
    )


def effective_polymer_selection(polymer_selection: str, polymer_types: list[str]) -> str:
    """Return the effective polymer selection encoded by contact identity.

    Parameters
    ----------
    polymer_selection : str
        Base polymer atom selection.
    polymer_types : list[str]
        Normalized polymer residue-name filter.

    Returns
    -------
    str
        Effective MDAnalysis selection string.
    """

    if not polymer_types:
        return polymer_selection
    return f"({polymer_selection}) and (resname {' '.join(polymer_types)})"


def _normalized_text(value: Any) -> str:
    """Normalize scalar text values for stable contact identity.

    Parameters
    ----------
    value : Any
        Value to coerce to normalized text.

    Returns
    -------
    str
        Stripped string representation.
    """

    return str(value).strip()
