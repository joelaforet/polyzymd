"""Filesystem artifact store for MDAnalysis extension-layer outputs."""

from __future__ import annotations

import hashlib
from pathlib import Path, PurePosixPath
from typing import Any

from polyzymd.analyses.mda.artifacts import (
    ArtifactManifest,
    ArtifactSidecarRef,
    ReplicateArtifact,
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
        """

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
