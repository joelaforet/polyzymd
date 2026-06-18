"""Moiety descriptor normalization for conjugation planning."""

from __future__ import annotations

from pathlib import Path
from typing import Any

from pydantic import BaseModel, Field, field_validator

from polyzymd.config.schema import ConjugationMoietyConfig


class MoietyDescriptor(BaseModel):
    """Normalized declarative description of an attachment moiety."""

    name: str = Field(..., min_length=1, description="Moiety identifier")
    role: str = Field(..., min_length=1, description="Normalized component role")
    kind: str = Field(..., min_length=1, description="Normalized descriptor kind")
    residue_name: str | None = Field(None, description="Optional moiety residue name")
    input_path: Path | None = Field(None, description="Optional structure or template input path")
    smiles: str | None = Field(None, description="Optional declarative SMILES string")
    metadata: dict[str, Any] = Field(
        default_factory=dict, description="Serializable descriptor notes"
    )

    @field_validator("name")
    @classmethod
    def normalize_name(cls, value: str) -> str:
        """Normalize a moiety name while preserving user-facing capitalization."""
        normalized = value.strip()
        if not normalized:
            raise ValueError("Moiety names must be non-empty strings")
        return normalized

    @field_validator("role", "kind")
    @classmethod
    def normalize_classifier(cls, value: str) -> str:
        """Normalize classifier labels for deterministic metadata."""
        normalized = value.strip().lower()
        if not normalized:
            raise ValueError("Moiety classifier labels must be non-empty strings")
        return normalized

    @field_validator("residue_name")
    @classmethod
    def normalize_residue_name(cls, value: str | None) -> str | None:
        """Normalize optional moiety residue names."""
        if value is None:
            return None
        normalized = value.strip().upper()
        if not normalized:
            raise ValueError("Moiety residue names must be non-empty strings when provided")
        return normalized


def normalize_moiety_descriptor(config: ConjugationMoietyConfig) -> MoietyDescriptor:
    """Normalize a conjugation moiety config into a descriptor.

    Parameters
    ----------
    config : ConjugationMoietyConfig
        User configuration for one attachment moiety.

    Returns
    -------
    MoietyDescriptor
        Normalized declarative moiety descriptor.
    """
    role = config.role.strip().lower()
    kind = _classify_moiety_kind(config, role)
    metadata: dict[str, Any] = {
        "source": kind,
        "declared_role": config.role,
        "has_smiles": config.smiles is not None,
        "has_input_path": config.input_path is not None,
    }
    if config.input_path is not None:
        metadata["input_path"] = str(config.input_path)
        metadata["input_suffix"] = config.input_path.suffix.lower()
    if config.smiles is not None:
        metadata["smiles"] = config.smiles
    if config.residue_name is not None:
        metadata["residue_name"] = config.residue_name.upper()

    return MoietyDescriptor(
        name=config.name,
        role=role,
        kind=kind,
        residue_name=config.residue_name,
        input_path=config.input_path,
        smiles=config.smiles,
        metadata=metadata,
    )


def _classify_moiety_kind(config: ConjugationMoietyConfig, role: str) -> str:
    """Classify a moiety descriptor using available config fields."""
    if config.input_path is not None:
        return "file"
    if config.smiles is not None:
        return "smiles"
    if role in {"glycan", "polymer", "ptm", "residue", "moiety"}:
        return f"declarative_{role}"
    if config.residue_name is not None:
        return "declarative_residue"
    return "declarative_moiety"
