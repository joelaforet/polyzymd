"""Artifact envelope models for MDAnalysis extension-layer results."""

from __future__ import annotations

from pathlib import PurePosixPath
from typing import Any, Literal

from pydantic import BaseModel, ConfigDict, Field, field_validator

MDA_ARTIFACT_SCHEMA_VERSION: str = "1"


def validate_artifact_relative_path(value: str) -> str:
    """Validate an artifact-relative POSIX path string.

    Parameters
    ----------
    value : str
        Candidate path stored in a JSON artifact.

    Returns
    -------
    str
        Normalized relative path using POSIX separators.

    Raises
    ------
    ValueError
        Raised when the path is empty, absolute, or contains parent traversal.
    """

    if not isinstance(value, str):
        raise TypeError("sidecar path must be a string")
    path = PurePosixPath(value)
    if not value or path.is_absolute():
        raise ValueError("sidecar path must be relative")
    if any(part in {"", ".."} for part in path.parts):
        raise ValueError("sidecar path must not contain empty or parent-traversal parts")
    if str(path) == ".":
        raise ValueError("sidecar path must reference a file")
    return str(path)


class ArtifactSidecarRef(BaseModel):
    """Relative reference to a sidecar file owned by an artifact store."""

    path: str = Field(description="Artifact-store-relative POSIX path")
    sha256: str = Field(description="SHA-256 hex digest of the sidecar bytes")
    size_bytes: int = Field(ge=0, description="Sidecar size in bytes")
    media_type: str | None = Field(default=None, description="Optional media or content type")
    metadata: dict[str, Any] = Field(default_factory=dict, description="Sidecar-specific metadata")

    model_config = ConfigDict(extra="allow")

    @field_validator("path")
    @classmethod
    def _validate_path(cls, value: str) -> str:
        """Validate the stored sidecar path.

        Parameters
        ----------
        value : str
            Candidate sidecar path.

        Returns
        -------
        str
            Normalized relative POSIX sidecar path.
        """

        return validate_artifact_relative_path(value)

    @field_validator("sha256")
    @classmethod
    def _validate_sha256(cls, value: str) -> str:
        """Validate that a sidecar ref uses a SHA-256 digest.

        Parameters
        ----------
        value : str
            Candidate digest string.

        Returns
        -------
        str
            Lowercase SHA-256 hex digest.

        Raises
        ------
        ValueError
            Raised when the digest is not a 64-character hexadecimal SHA-256 value.
        """

        digest = value.lower()
        if len(digest) != 64 or any(char not in "0123456789abcdef" for char in digest):
            raise ValueError("sidecar hashes must be SHA-256 hex digests")
        return digest


class ArtifactManifest(BaseModel):
    """Manifest for one artifact directory and its sidecar files."""

    schema_version: str = Field(default=MDA_ARTIFACT_SCHEMA_VERSION)
    analysis_name: str
    artifact_id: str | None = None
    artifact_type: str = "manifest"
    polyzymd_version: str | None = None
    mdanalysis_version: str | None = None
    inputs: dict[str, Any] = Field(default_factory=dict)
    provenance: dict[str, Any] = Field(default_factory=dict)
    sidecars: list[ArtifactSidecarRef] = Field(default_factory=list)
    metadata: dict[str, Any] = Field(default_factory=dict)

    model_config = ConfigDict(extra="allow")

    @field_validator("schema_version")
    @classmethod
    def _validate_schema_version(cls, value: str) -> str:
        """Validate the artifact schema version.

        Parameters
        ----------
        value : str
            Candidate schema version.

        Returns
        -------
        str
            Supported schema version.
        """

        if value != MDA_ARTIFACT_SCHEMA_VERSION:
            raise ValueError(
                f"Unsupported MDA artifact schema version {value!r}; "
                f"expected {MDA_ARTIFACT_SCHEMA_VERSION!r}"
            )
        return value


class ArtifactEnvelope(BaseModel):
    """Extensible JSON envelope for MDAnalysis extension-layer artifacts."""

    schema_version: str = Field(default=MDA_ARTIFACT_SCHEMA_VERSION)
    artifact_type: str = "artifact"
    analysis_name: str
    payload: dict[str, Any] = Field(default_factory=dict)
    sidecars: list[ArtifactSidecarRef] = Field(default_factory=list)
    provenance: dict[str, Any] = Field(default_factory=dict)
    metadata: dict[str, Any] = Field(default_factory=dict)
    warnings: list[str] = Field(default_factory=list)

    model_config = ConfigDict(extra="allow")

    @field_validator("schema_version")
    @classmethod
    def _validate_schema_version(cls, value: str) -> str:
        """Validate the artifact schema version.

        Parameters
        ----------
        value : str
            Candidate schema version.

        Returns
        -------
        str
            Supported schema version.
        """

        if value != MDA_ARTIFACT_SCHEMA_VERSION:
            raise ValueError(
                f"Unsupported MDA artifact schema version {value!r}; "
                f"expected {MDA_ARTIFACT_SCHEMA_VERSION!r}"
            )
        return value


class ReplicateArtifact(ArtifactEnvelope):
    """Result artifact produced for one replicate trajectory."""

    artifact_type: Literal["replicate"] = "replicate"

    condition_label: str
    replicate: int = Field(ge=1)


class ConditionArtifact(ArtifactEnvelope):
    """Aggregated artifact produced for one simulation condition."""

    artifact_type: Literal["condition"] = "condition"

    condition_label: str
    replicates: list[int] = Field(default_factory=list)


class ComparisonArtifact(ArtifactEnvelope):
    """Cross-condition comparison artifact."""

    artifact_type: Literal["comparison"] = "comparison"

    conditions: list[str] = Field(default_factory=list)
    control_label: str | None = None
    effective_control: str | None = None
