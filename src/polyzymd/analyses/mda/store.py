"""Filesystem artifact store for MDAnalysis extension-layer outputs."""

from __future__ import annotations

import hashlib
from pathlib import Path, PurePosixPath
from typing import Any, Sequence

from pydantic import ValidationError

from polyzymd.analyses.exceptions import AggregateValidationError
from polyzymd.analyses.mda.artifacts import (
    ArtifactManifest,
    ArtifactSidecarRef,
    ComparisonArtifact,
    ConditionArtifact,
    ReplicateArtifact,
    stamp_software_versions,
    validate_artifact_relative_path,
)
from polyzymd.analyses.mda.base import MDAnalysisExtensionError

_HASH_CHUNK_SIZE = 1024 * 1024


class ArtifactStoreError(MDAnalysisExtensionError):
    """Error raised when artifact-store path or validation checks fail."""


class ArtifactStore:
    """Store JSON artifacts and relative sidecars under one root directory."""

    def __init__(self, root: str | Path) -> None:
        """Create an artifact store rooted at ``root``.

        Parameters
        ----------
        root : str or Path
            Directory that owns artifact JSON files and sidecars.
        """

        self.root = Path(root).expanduser().resolve()

    def write_replicate_result(
        self,
        artifact: ReplicateArtifact,
        path: str | Path = "result.json",
    ) -> Path:
        """Write a replicate artifact JSON file.

        Parameters
        ----------
        artifact : ReplicateArtifact
            Replicate artifact envelope to serialize.
        path : str or Path, optional
            Store-relative JSON path, by default ``"result.json"``.

        Returns
        -------
        Path
            Absolute path to the written JSON file.
        """

        return self._write_json_model(artifact, path)

    def read_replicate_result(self, path: str | Path = "result.json") -> ReplicateArtifact:
        """Read a replicate artifact JSON file.

        Parameters
        ----------
        path : str or Path, optional
            Store-relative JSON path, by default ``"result.json"``.

        Returns
        -------
        ReplicateArtifact
            Deserialized replicate artifact.
        """

        resolved_path = self._resolve_relative_path(path)
        try:
            return ReplicateArtifact.model_validate_json(resolved_path.read_text())
        except OSError as exc:
            raise ArtifactStoreError(
                f"Failed to read replicate artifact {resolved_path}: {exc}"
            ) from exc
        except ValidationError as exc:
            raise ArtifactStoreError(
                f"Failed to validate replicate artifact {resolved_path}: {exc}"
            ) from exc

    def write_condition_result(
        self,
        artifact: ConditionArtifact,
        path: str | Path = "result.json",
    ) -> Path:
        """Write a condition artifact JSON file.

        Parameters
        ----------
        artifact : ConditionArtifact
            Condition artifact envelope to serialize.
        path : str or Path, optional
            Store-relative JSON path, by default ``"result.json"``.

        Returns
        -------
        Path
            Absolute path to the written JSON file.
        """

        return self._write_json_model(artifact, path)

    def read_condition_result(self, path: str | Path = "result.json") -> ConditionArtifact:
        """Read a condition artifact JSON file.

        Parameters
        ----------
        path : str or Path, optional
            Store-relative JSON path, by default ``"result.json"``.

        Returns
        -------
        ConditionArtifact
            Deserialized condition artifact.
        """

        resolved_path = self._resolve_relative_path(path)
        try:
            artifact = ConditionArtifact.model_validate_json(resolved_path.read_text())
        except (OSError, ValidationError) as exc:
            raise ArtifactStoreError(
                f"Failed to validate condition artifact {resolved_path}: {exc}"
            ) from exc
        # Every canonical aggregate read funnels through here, including the
        # plugin ``_load_aggregated_result`` overrides that bypass the framework
        # loader, so this is where the staleness check cannot be routed around.
        validate_aggregate_not_outdated(
            artifact,
            analysis_name=artifact.analysis_name,
            source=resolved_path,
        )
        return artifact

    def write_comparison_result(
        self,
        artifact: ComparisonArtifact,
        path: str | Path = "result.json",
    ) -> Path:
        """Write a comparison artifact JSON file.

        Parameters
        ----------
        artifact : ComparisonArtifact
            Comparison artifact envelope to serialize.
        path : str or Path, optional
            Store-relative JSON path, by default ``"result.json"``.

        Returns
        -------
        Path
            Absolute path to the written JSON file.
        """

        return self._write_json_model(artifact, path)

    def read_comparison_result(self, path: str | Path = "result.json") -> ComparisonArtifact:
        """Read a comparison artifact JSON file.

        Parameters
        ----------
        path : str or Path, optional
            Store-relative JSON path, by default ``"result.json"``.

        Returns
        -------
        ComparisonArtifact
            Deserialized comparison artifact.
        """

        resolved_path = self._resolve_relative_path(path)
        try:
            return ComparisonArtifact.model_validate_json(resolved_path.read_text())
        except (OSError, ValidationError) as exc:
            raise ArtifactStoreError(
                f"Failed to validate comparison artifact {resolved_path}: {exc}"
            ) from exc

    def write_manifest(
        self,
        manifest: ArtifactManifest,
        path: str | Path = "manifest.json",
    ) -> Path:
        """Write an artifact manifest JSON file.

        Parameters
        ----------
        manifest : ArtifactManifest
            Manifest to serialize.
        path : str or Path, optional
            Store-relative JSON path, by default ``"manifest.json"``.

        Returns
        -------
        Path
            Absolute path to the written manifest.
        """

        return self._write_json_model(manifest, path)

    def read_manifest(self, path: str | Path = "manifest.json") -> ArtifactManifest:
        """Read an artifact manifest JSON file.

        Parameters
        ----------
        path : str or Path, optional
            Store-relative JSON path, by default ``"manifest.json"``.

        Returns
        -------
        ArtifactManifest
            Deserialized manifest.
        """

        resolved_path = self._resolve_relative_path(path)
        try:
            return ArtifactManifest.model_validate_json(resolved_path.read_text())
        except OSError as exc:
            raise ArtifactStoreError(
                f"Failed to read artifact manifest {resolved_path}: {exc}"
            ) from exc

    def register_sidecar(
        self,
        path: str | Path,
        *,
        media_type: str | None = None,
        metadata: dict[str, Any] | None = None,
    ) -> ArtifactSidecarRef:
        """Register an existing store-relative sidecar file.

        Parameters
        ----------
        path : str or Path
            Store-relative sidecar path.
        media_type : str or None, optional
            Optional sidecar content type.
        metadata : dict[str, Any] or None, optional
            Sidecar metadata to copy into the reference.

        Returns
        -------
        ArtifactSidecarRef
            Relative sidecar reference with streamed SHA-256 and size metadata.
        """

        resolved_path = self._resolve_relative_path(path)
        if not resolved_path.is_file():
            raise ArtifactStoreError(f"Sidecar does not exist: {resolved_path}")
        stat = resolved_path.stat()
        return ArtifactSidecarRef(
            path=self._relative_ref_for(resolved_path),
            sha256=self._sha256_file(resolved_path),
            size_bytes=stat.st_size,
            media_type=media_type,
            metadata=dict(metadata or {}),
        )

    def write_npz_sidecar(
        self,
        path: str | Path,
        *,
        compressed: bool = True,
        media_type: str | None = "application/x-npz",
        metadata: dict[str, Any] | None = None,
        **arrays: Any,
    ) -> ArtifactSidecarRef:
        """Write NumPy arrays to an NPZ sidecar and register it.

        Parameters
        ----------
        path : str or Path
            Store-relative sidecar path.
        compressed : bool, optional
            Whether to use ``numpy.savez_compressed``, by default True.
        media_type : str or None, optional
            Sidecar media type stored in the reference.
        metadata : dict[str, Any] or None, optional
            Additional sidecar metadata.
        **arrays : Any
            Named arrays forwarded to NumPy.

        Returns
        -------
        ArtifactSidecarRef
            Registered NPZ sidecar reference.
        """

        import numpy as np

        resolved_path = self._resolve_relative_path(path)
        resolved_path.parent.mkdir(parents=True, exist_ok=True)
        writer = np.savez_compressed if compressed else np.savez
        try:
            writer(resolved_path, **arrays)
        except OSError as exc:
            raise ArtifactStoreError(f"Failed to write NPZ sidecar {resolved_path}: {exc}") from exc
        return self.register_sidecar(
            resolved_path.relative_to(self.root), media_type=media_type, metadata=metadata
        )

    def resolve_sidecar(self, ref: ArtifactSidecarRef | str | Path) -> Path:
        """Resolve a sidecar reference to an absolute path under the store root.

        Parameters
        ----------
        ref : ArtifactSidecarRef or str or Path
            Sidecar reference or store-relative path.

        Returns
        -------
        Path
            Absolute sidecar path under the store root.
        """

        path = ref.path if isinstance(ref, ArtifactSidecarRef) else ref
        return self._resolve_relative_path(path)

    def validate_sidecar(self, ref: ArtifactSidecarRef) -> Path:
        """Validate sidecar existence, size, path containment, and SHA-256.

        Parameters
        ----------
        ref : ArtifactSidecarRef
            Sidecar reference to validate.

        Returns
        -------
        Path
            Absolute sidecar path when validation passes.
        """

        resolved_path = self.resolve_sidecar(ref)
        if not resolved_path.is_file():
            raise ArtifactStoreError(f"Missing sidecar: {resolved_path}")
        size_bytes = resolved_path.stat().st_size
        if size_bytes != ref.size_bytes:
            raise ArtifactStoreError(
                f"Sidecar size mismatch for {ref.path}: expected {ref.size_bytes}, got {size_bytes}"
            )
        actual_sha256 = self._sha256_file(resolved_path)
        if actual_sha256 != ref.sha256:
            raise ArtifactStoreError(
                f"Sidecar SHA-256 mismatch for {ref.path}: expected {ref.sha256}, got {actual_sha256}"
            )
        return resolved_path

    def load_npz_sidecar(self, ref: ArtifactSidecarRef) -> Any:
        """Validate and open an NPZ sidecar.

        Parameters
        ----------
        ref : ArtifactSidecarRef
            Sidecar reference to validate before loading.

        Returns
        -------
        Any
            Open ``numpy.load`` handle. Callers should use it as a context
            manager to close the underlying file promptly.
        """

        import numpy as np

        resolved_path = self.validate_sidecar(ref)
        try:
            return np.load(resolved_path)
        except (OSError, ValueError) as exc:
            raise ArtifactStoreError(f"Failed to load NPZ sidecar {resolved_path}: {exc}") from exc

    def source_artifact_ref(self, path: str | Path = "result.json") -> dict[str, Any]:
        """Return a hashed reference for an artifact JSON file.

        Parameters
        ----------
        path : str or Path, optional
            Store-relative artifact path, by default ``"result.json"``.

        Returns
        -------
        dict[str, Any]
            Relative path, SHA-256 digest, and size metadata for the artifact.
        """

        resolved_path = self._resolve_relative_path(path)
        if not resolved_path.is_file():
            raise ArtifactStoreError(f"Artifact source does not exist: {resolved_path}")
        stat = resolved_path.stat()
        return {
            "path": self._relative_ref_for(resolved_path),
            "sha256": self._sha256_file(resolved_path),
            "size_bytes": stat.st_size,
        }

    def _write_json_model(self, model: Any, path: str | Path) -> Path:
        """Write a Pydantic-like model to a store-relative JSON path.

        Parameters
        ----------
        model : Any
            Object exposing ``model_dump_json``.
        path : str or Path
            Store-relative output path.

        Returns
        -------
        Path
            Absolute output path.

        Notes
        -----
        Every artifact the framework persists goes through here, so this is
        where the running software versions are recorded. Stamping at the
        collector instead left comparison artifacts and direct compute returns
        with null versions.
        """

        stamp_software_versions(model)
        resolved_path = self._resolve_relative_path(path)
        resolved_path.parent.mkdir(parents=True, exist_ok=True)
        try:
            resolved_path.write_text(model.model_dump_json(indent=2), encoding="utf-8")
        except OSError as exc:
            raise ArtifactStoreError(
                f"Failed to write artifact JSON {resolved_path}: {exc}"
            ) from exc
        return resolved_path

    def _resolve_relative_path(self, path: str | Path) -> Path:
        """Resolve a store-relative path and reject root escape attempts.

        Parameters
        ----------
        path : str or Path
            Candidate store-relative path.

        Returns
        -------
        Path
            Absolute path under the store root.
        """

        relative_path = self._normalize_relative_path(path)
        resolved_path = (self.root / relative_path).resolve(strict=False)
        if not resolved_path.is_relative_to(self.root):
            raise ArtifactStoreError(f"Artifact path escapes store root: {path}")
        return resolved_path

    @staticmethod
    def _normalize_relative_path(path: str | Path) -> Path:
        """Normalize a candidate relative path for filesystem use.

        Parameters
        ----------
        path : str or Path
            Candidate relative path.

        Returns
        -------
        Path
            Relative filesystem path.
        """

        if isinstance(path, Path) and path.is_absolute():
            raise ArtifactStoreError(f"Artifact paths must be relative: {path}")
        candidate = PurePosixPath(str(path).replace("\\", "/"))
        try:
            relative = validate_artifact_relative_path(str(candidate))
        except (TypeError, ValueError) as exc:
            raise ArtifactStoreError(str(exc)) from exc
        return Path(relative)

    def _relative_ref_for(self, path: Path) -> str:
        """Return a validated POSIX reference for an absolute store path.

        Parameters
        ----------
        path : Path
            Absolute path that must live under the store root.

        Returns
        -------
        str
            Store-relative POSIX reference.
        """

        resolved_path = path.resolve(strict=False)
        if not resolved_path.is_relative_to(self.root):
            raise ArtifactStoreError(f"Artifact path escapes store root: {path}")
        relative_path = resolved_path.relative_to(self.root).as_posix()
        return validate_artifact_relative_path(relative_path)

    @staticmethod
    def _sha256_file(path: Path) -> str:
        """Stream a file into a SHA-256 digest.

        Parameters
        ----------
        path : Path
            File to hash.

        Returns
        -------
        str
            Lowercase SHA-256 hexadecimal digest.
        """

        digest = hashlib.sha256()
        with path.open("rb") as handle:
            for chunk in iter(lambda: handle.read(_HASH_CHUNK_SIZE), b""):
                digest.update(chunk)
        return digest.hexdigest()


def validate_aggregate_not_outdated(
    result: Any,
    *,
    analysis_name: str,
    source: str | Path | None,
    expected_replicates: Sequence[int] | None = None,
) -> None:
    """Reject an aggregate that is older than the replicate results it covers.

    An aggregate summarises the ``run_N/result.json`` files beside it. If one
    of those files was written after the aggregate, the aggregate describes a
    replicate result that no longer exists.

    Call this only for an aggregate that was just read from ``source``. An
    aggregate that is about to be written there is newer than every replicate
    result by construction, even though the file still on disk is not.

    Parameters
    ----------
    result : Any
        Coerced aggregate result.
    expected_replicates : sequence of int or None
        Replicate IDs requested by the caller, used when the aggregate records
        none itself.
    analysis_name : str
        Analysis name used in diagnostics.
    source : str or Path or None
        Path the aggregate was loaded from. The check is skipped for aggregates
        that were not loaded from a file.
    """

    if not isinstance(source, Path) or not source.is_file():
        return
    condition_dir = source.parent.parent
    recorded = [
        entry
        for entry in getattr(result, "source_replicates", None) or []
        if isinstance(entry, dict) and entry.get("fingerprint") and "replicate" in entry
    ]
    if recorded:
        # Content, not modification times, decides: a copied, unzipped or cloned
        # study gets new modification times in no particular order.
        from polyzymd.analyses.identity import file_fingerprint

        changed = []
        for entry in recorded:
            replicate_path = condition_dir / f"run_{int(entry['replicate'])}" / source.name
            if (
                replicate_path.is_file()
                and file_fingerprint(replicate_path) != entry["fingerprint"]
            ):
                changed.append(str(replicate_path))
        if changed:
            raise _validation_error(
                analysis_name,
                source,
                "replicate result(s) changed after they were aggregated ("
                + ", ".join(changed)
                + "); rerun with --recompute to rebuild it",
            )
        return
    replicates = _replicates_from_result(result)
    if replicates is None:
        replicates = tuple(expected_replicates) if expected_replicates is not None else ()
    aggregate_mtime_ns = source.stat().st_mtime_ns
    newer: list[str] = []
    for replicate in replicates:
        replicate_path = condition_dir / f"run_{int(replicate)}" / source.name
        if not replicate_path.is_file():
            continue
        if replicate_path.stat().st_mtime_ns > aggregate_mtime_ns:
            newer.append(str(replicate_path))
    if newer:
        raise _validation_error(
            analysis_name,
            source,
            "aggregate is older than the replicate result(s) it summarizes ("
            + ", ".join(newer)
            + "); rerun with --recompute to rebuild it",
        )


def _replicates_from_result(result: Any) -> tuple[int, ...] | None:
    """Return aggregate replicate IDs when the result declares them."""

    for field_name in ("replicates", "replicate_ids"):
        value = _field_value(result, field_name)
        if value is not None:
            try:
                return tuple(int(rep) for rep in value)
            except (TypeError, ValueError) as exc:
                raise AggregateValidationError(
                    "Aggregate result has invalid replicate identity. Recompute the analysis "
                    "or clear stale analysis cache files."
                ) from exc
    for field_name in ("replicates", "replicate_ids"):
        value = _metadata_value(result, field_name)
        if value is not None:
            try:
                return tuple(int(rep) for rep in value)
            except (TypeError, ValueError) as exc:
                raise AggregateValidationError(
                    "Aggregate result metadata has invalid replicate identity. Recompute the "
                    "analysis or clear stale analysis cache files."
                ) from exc
    return None


def _validation_error(
    analysis_name: str,
    source: str | Path | None,
    detail: str,
) -> AggregateValidationError:
    """Build a user-actionable aggregate validation error."""

    source_text = f" at {source}" if source is not None else ""
    return AggregateValidationError(
        f"{analysis_name}: stale or incompatible aggregated result{source_text}: {detail}. "
        "Recompute the analysis with --recompute or clear stale analysis cache files."
    )


def _field_value(result: Any, field_name: str) -> Any:
    """Return a top-level result field from models or mappings."""

    if isinstance(result, dict):
        return result.get(field_name)
    return getattr(result, field_name, None)


def _metadata_value(result: Any, field_name: str) -> Any:
    """Return a metadata field from models or mappings."""

    metadata = _field_value(result, "metadata")
    if isinstance(metadata, dict):
        return metadata.get(field_name)
    return None
