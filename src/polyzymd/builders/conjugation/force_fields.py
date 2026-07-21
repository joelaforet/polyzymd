"""Attachment-scoped force-field resolution for conjugate construction."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any, Literal

ForceFieldSourceKind = Literal["glycam", "openff"]
ConjugateForceFieldRoute = Literal[
    "standard_interchange",
    "native_exact",
    "mixed_overlay",
]

_GLYCAM_ALIASES = {"glycam06", "glycam_06", "glycam-06", "glycam", "glycam06j"}


@dataclass(frozen=True)
class ResolvedAttachmentForceField:
    """Resolved force-field source for one enabled attachment.

    The resolver records the configured source without carrying heavyweight force-field
    objects. The result is safe to serialize through :meth:`to_dict`.
    """

    attachment_name: str
    source: str
    source_kind: ForceFieldSourceKind

    @property
    def is_glycam(self) -> bool:
        """Return whether the attachment is assigned to canonical GLYCAM."""

        return self.source_kind == "glycam"

    def to_dict(self) -> dict[str, Any]:
        """Return a JSON-serializable representation."""

        return {
            "attachment_name": self.attachment_name,
            "source": self.source,
            "source_kind": self.source_kind,
        }


@dataclass(frozen=True)
class ResolvedConjugateForceFields:
    """Resolved force-field route and attachment assignments."""

    route: ConjugateForceFieldRoute
    attachments: tuple[ResolvedAttachmentForceField, ...]
    protein_force_field: str

    @property
    def has_glycam(self) -> bool:
        """Return whether any attachment uses GLYCAM."""

        return any(attachment.is_glycam for attachment in self.attachments)

    @property
    def has_generic(self) -> bool:
        """Return whether any attachment uses an OpenFF-compatible source."""

        return any(not attachment.is_glycam for attachment in self.attachments)

    def to_dict(self) -> dict[str, Any]:
        """Return a JSON-serializable route manifest."""

        return {
            "route": self.route,
            "protein_force_field": self.protein_force_field,
            "attachments": [attachment.to_dict() for attachment in self.attachments],
        }


def resolve_conjugate_force_fields(config: Any) -> ResolvedConjugateForceFields:
    """Resolve enabled conjugation attachments to force-field sources and a route.

    Every enabled attachment must explicitly declare its owner. Canonical GLYCAM
    is requested with ``glycam06``. OpenFF OFFXML names and filesystem paths use
    the generic Interchange route. Missing, blank, and unknown values fail closed.
    """

    force_field = getattr(config, "force_field", None)
    protein_force_field = str(getattr(force_field, "protein", ""))

    conjugation = getattr(config, "conjugation", None)
    attachments = tuple(
        attachment
        for attachment in tuple(getattr(conjugation, "attachments", ()) or ())
        if getattr(conjugation, "enabled", True)
        if getattr(attachment, "enabled", True)
    )
    resolved = tuple(_resolve_attachment_force_field(attachment) for attachment in attachments)
    route = _route_from_resolved_attachments(resolved)
    return ResolvedConjugateForceFields(
        route=route,
        attachments=resolved,
        protein_force_field=protein_force_field,
    )


def native_glycam_only_route(config: Any) -> bool:
    """Return whether the config uses the exact native GLYCAM-only route."""

    return resolve_conjugate_force_fields(config).route == "native_exact"


def mixed_overlay_route(config: Any) -> bool:
    """Return whether the config requires the mixed GLYCAM/OpenFF overlay route."""

    return resolve_conjugate_force_fields(config).route == "mixed_overlay"


def _resolve_attachment_force_field(attachment: Any) -> ResolvedAttachmentForceField:
    """Resolve one attachment force-field value."""

    moiety = getattr(attachment, "moiety", None)
    raw_value = getattr(moiety, "force_field", None)
    attachment_name = str(getattr(attachment, "name", "attachment"))
    if raw_value is None:
        raise ValueError(
            f"Enabled conjugation attachment {attachment_name!r} must explicitly declare "
            "moiety.force_field"
        )
    value = str(raw_value)
    source, source_kind = _canonical_source(value)
    return ResolvedAttachmentForceField(
        attachment_name=attachment_name,
        source=source,
        source_kind=source_kind,
    )


def _canonical_source(value: str) -> tuple[str, ForceFieldSourceKind]:
    """Canonicalize one configured force-field source or raise."""

    normalized = value.strip()
    if not normalized:
        raise ValueError("Conjugation moiety force_field must not be blank")
    label = normalized.lower().replace("_", "-")
    alias_key = normalized.lower().replace("-", "_")
    if label in _GLYCAM_ALIASES or alias_key in _GLYCAM_ALIASES:
        return "glycam06", "glycam"
    if normalized.lower().endswith(".offxml") or Path(normalized).suffix.lower() == ".offxml":
        return normalized, "openff"
    raise ValueError(
        "Unknown conjugation moiety force_field "
        f"{value!r}; use canonical 'glycam06' or an OpenFF .offxml source"
    )


def _route_from_resolved_attachments(
    attachments: tuple[ResolvedAttachmentForceField, ...],
) -> ConjugateForceFieldRoute:
    """Return the construction route implied by resolved attachments."""

    if not attachments:
        return "standard_interchange"
    has_glycam = any(attachment.is_glycam for attachment in attachments)
    has_generic = any(not attachment.is_glycam for attachment in attachments)
    if has_glycam and has_generic:
        return "mixed_overlay"
    if has_glycam:
        return "native_exact"
    return "standard_interchange"
