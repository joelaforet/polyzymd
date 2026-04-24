"""Attachment site normalization helpers for conjugation planning."""

from __future__ import annotations

from pydantic import BaseModel, Field, field_validator

from polyzymd.builders.conjugation.exceptions import ConjugationNotImplementedError
from polyzymd.builders.conjugation.mechanisms import ReactionMechanism, SiteAtomRule
from polyzymd.config.schema import ConjugationSiteConfig


class AttachmentSite(BaseModel):
    """Normalized explicit atom attachment site."""

    chain_id: str = Field(..., min_length=1, max_length=1, description="Normalized chain ID")
    residue_name: str = Field(..., min_length=1, description="Normalized residue name")
    residue_number: int = Field(..., description="Residue number in the input structure")
    atom_name: str = Field(..., min_length=1, description="Normalized reactive atom name")
    selector: str | None = Field(None, description="Original selector expression, if any")
    selector_mode: str = Field("explicit_atom", description="Normalized selector mode")

    @field_validator("chain_id")
    @classmethod
    def normalize_chain_id(cls, value: str) -> str:
        """Normalize and validate a one-character chain ID."""
        normalized = value.strip().upper()
        if len(normalized) != 1 or not normalized.isalpha():
            raise ValueError("Attachment site chain ID must be a single alphabetic character")
        return normalized

    @field_validator("residue_name", "atom_name")
    @classmethod
    def normalize_labels(cls, value: str) -> str:
        """Normalize residue and atom labels for mechanism matching."""
        normalized = value.strip().upper()
        if not normalized:
            raise ValueError("Attachment site residue and atom labels must be non-empty strings")
        return normalized


def normalize_attachment_site(config: ConjugationSiteConfig) -> AttachmentSite:
    """Normalize an explicit atom site from conjugation configuration.

    Parameters
    ----------
    config : ConjugationSiteConfig
        User configuration for one attachment site.

    Returns
    -------
    AttachmentSite
        Normalized explicit atom site.

    Raises
    ------
    ValueError
        If required explicit site fields are missing.
    ConjugationNotImplementedError
        If a future selector-only mode is requested.
    """
    has_selector = bool(config.selector and config.selector.strip())
    explicit_values = {
        "residue_name": config.residue_name,
        "residue_number": config.residue_number,
        "atom_name": config.atom_name,
    }
    missing = [name for name, value in explicit_values.items() if value is None]

    if has_selector and missing:
        raise ConjugationNotImplementedError(
            "Selector-based conjugation site discovery is not implemented yet. Provide explicit "
            "residue_name, residue_number, and atom_name fields for this phase."
        )
    if missing:
        missing_text = ", ".join(missing)
        raise ValueError(
            "Explicit conjugation sites require residue_name, residue_number, and atom_name; "
            f"missing: {missing_text}"
        )

    return AttachmentSite(
        chain_id=config.chain_id,
        residue_name=config.residue_name,
        residue_number=config.residue_number,
        atom_name=config.atom_name,
        selector=config.selector,
    )


def match_site_rule(site: AttachmentSite, mechanism: ReactionMechanism) -> SiteAtomRule:
    """Return the mechanism rule matching an attachment site.

    Parameters
    ----------
    site : AttachmentSite
        Normalized attachment site.
    mechanism : ReactionMechanism
        Mechanism declaration to validate against.

    Returns
    -------
    SiteAtomRule
        Matching allowed site rule.

    Raises
    ------
    ValueError
        If the site does not satisfy any mechanism rule.
    """
    for rule in mechanism.allowed_sites:
        residue_matches = rule.residue_name == "ANY" or rule.residue_name == site.residue_name
        if residue_matches and rule.atom_name == site.atom_name:
            return rule

    allowed = ", ".join(
        f"{rule.residue_name}:{rule.atom_name}" for rule in mechanism.allowed_sites
    )
    raise ValueError(
        f"Site {site.chain_id}:{site.residue_name}{site.residue_number}:{site.atom_name} is not "
        f"allowed for mechanism '{mechanism.identifier}'. Allowed rules: {allowed}"
    )
