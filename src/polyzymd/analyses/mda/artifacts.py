"""Artifact envelope models for MDAnalysis extension-layer results."""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from pathlib import PurePosixPath
from typing import Any, Literal

from pydantic import BaseModel, ConfigDict, Field, field_validator, model_validator

MDA_ARTIFACT_SCHEMA_VERSION: str = "1"
_RAW_MDA_RESULTS_MODULE = "MDAnalysis.analysis.results"
_RAW_MDA_RESULTS_GROUP_CLASS_NAME = "ResultsGroup"


def is_raw_mdanalysis_results(value: Any) -> bool:
    """Return whether a value is an MDAnalysis raw results container.

    Detection is intentionally import-light and relies on module/class metadata
    so artifact validation never imports MDAnalysis.

    Parameters
    ----------
    value : Any
        Candidate value to inspect.

    Returns
    -------
    bool
        ``True`` when ``value`` looks like an MDAnalysis ``Results`` object.
    """

    value_type = type(value)
    class_name = getattr(value_type, "__name__", "")
    return getattr(value_type, "__module__", "") == _RAW_MDA_RESULTS_MODULE and (
        class_name.endswith("Results") or class_name == _RAW_MDA_RESULTS_GROUP_CLASS_NAME
    )


def raw_mdanalysis_results_path(value: Any) -> str | None:
    """Return the nested path to raw MDAnalysis results, if present.

    Parameters
    ----------
    value : Any
        Candidate artifact field value.

    Returns
    -------
    str or None
        Human-readable nested path to the first raw results container, or
        ``None`` when no raw results are present.
    """

    return _raw_mdanalysis_results_path(value, path="$", seen=set())


def reject_raw_mdanalysis_results(value: Any, *, field_name: str) -> Any:
    """Reject raw MDAnalysis ``Results`` objects in artifact fields.

    Parameters
    ----------
    value : Any
        Candidate artifact field value.
    field_name : str
        Name used in validation diagnostics.

    Returns
    -------
    Any
        The original value when no raw results are found.

    Raises
    ------
    ValueError
        Raised when raw MDAnalysis results are found recursively.
    """

    raw_path = raw_mdanalysis_results_path(value)
    if raw_path is not None:
        raise ValueError(
            f"{field_name} must not contain raw MDAnalysis Results at {raw_path}; "
            "map Results to JSON primitives or sidecar artifacts first"
        )
    return value


def _raw_mdanalysis_results_path(value: Any, *, path: str, seen: set[int]) -> str | None:
    """Recursively find raw MDAnalysis results without importing MDAnalysis.

    Parameters
    ----------
    value : Any
        Candidate value.
    path : str
        Current human-readable traversal path.
    seen : set of int
        Object IDs already visited to avoid cycles.

    Returns
    -------
    str or None
        Path to the first raw results object, if found.
    """

    if is_raw_mdanalysis_results(value):
        return path
    value_id = id(value)
    if value_id in seen:
        return None
    if isinstance(value, (Mapping, Sequence, BaseModel)) and not isinstance(
        value, (str, bytes, bytearray)
    ):
        seen.add(value_id)
    if isinstance(value, BaseModel):
        model_data: dict[str, Any] = dict(vars(value))
        model_extra = getattr(value, "model_extra", None)
        if isinstance(model_extra, Mapping):
            model_data.update(model_extra)
        return _raw_mdanalysis_results_path(model_data, path=path, seen=seen)
    if isinstance(value, Mapping):
        for key, item in value.items():
            key_path = _raw_mdanalysis_results_path(key, path=f"{path}.<key>", seen=seen)
            if key_path is not None:
                return key_path
            nested_path = f"{path}.{key}" if isinstance(key, str) else f"{path}[{key!r}]"
            item_path = _raw_mdanalysis_results_path(item, path=nested_path, seen=seen)
            if item_path is not None:
                return item_path
        return None
    if isinstance(value, Sequence) and not isinstance(value, (str, bytes, bytearray)):
        for index, item in enumerate(value):
            item_path = _raw_mdanalysis_results_path(item, path=f"{path}[{index}]", seen=seen)
            if item_path is not None:
                return item_path
    return None


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

    @model_validator(mode="before")
    @classmethod
    def _reject_raw_results_anywhere(cls, value: Any) -> Any:
        """Reject raw MDAnalysis results across sidecar fields and extras.

        Parameters
        ----------
        value : Any
            Candidate model input.

        Returns
        -------
        Any
            Original input when valid.
        """

        return reject_raw_mdanalysis_results(value, field_name="sidecar reference")

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

    @field_validator("metadata", mode="before")
    @classmethod
    def _reject_raw_results_metadata(cls, value: Any) -> Any:
        """Reject raw MDAnalysis results in sidecar metadata.

        Parameters
        ----------
        value : Any
            Candidate sidecar metadata.

        Returns
        -------
        Any
            Original metadata when valid.
        """

        return reject_raw_mdanalysis_results(value, field_name="sidecar metadata")

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

    @model_validator(mode="before")
    @classmethod
    def _reject_raw_results_anywhere(cls, value: Any) -> Any:
        """Reject raw MDAnalysis results across manifest fields and extras.

        Parameters
        ----------
        value : Any
            Candidate model input.

        Returns
        -------
        Any
            Original input when valid.
        """

        return reject_raw_mdanalysis_results(value, field_name="manifest")

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

    @field_validator("inputs", "provenance", "metadata", mode="before")
    @classmethod
    def _reject_raw_results_dicts(cls, value: Any) -> Any:
        """Reject raw MDAnalysis results in manifest dictionaries.

        Parameters
        ----------
        value : Any
            Candidate manifest dictionary field.

        Returns
        -------
        Any
            Original value when valid.
        """

        return reject_raw_mdanalysis_results(value, field_name="manifest field")


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

    @model_validator(mode="before")
    @classmethod
    def _reject_raw_results_anywhere(cls, value: Any) -> Any:
        """Reject raw MDAnalysis results across artifact fields and extras.

        Parameters
        ----------
        value : Any
            Candidate model input.

        Returns
        -------
        Any
            Original input when valid.
        """

        return reject_raw_mdanalysis_results(value, field_name="artifact")

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

    @field_validator("payload", "provenance", "metadata", mode="before")
    @classmethod
    def _reject_raw_results_dicts(cls, value: Any) -> Any:
        """Reject raw MDAnalysis results in artifact dictionaries.

        Parameters
        ----------
        value : Any
            Candidate artifact dictionary field.

        Returns
        -------
        Any
            Original value when valid.
        """

        return reject_raw_mdanalysis_results(value, field_name="artifact field")


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
    source_replicates: list[dict[str, Any]] = Field(default_factory=list)
    skipped_replicates: list[dict[str, Any]] = Field(default_factory=list)


class ComparisonArtifact(ArtifactEnvelope):
    """Cross-condition comparison artifact."""

    artifact_type: Literal["comparison"] = "comparison"

    conditions: list[str] = Field(default_factory=list)
    control_label: str | None = None
    effective_control: str | None = None
