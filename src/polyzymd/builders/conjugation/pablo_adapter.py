"""Lazy OpenFF Pablo adapter foundation for covalent modification ingestion."""

from __future__ import annotations

import importlib
import importlib.metadata
from pathlib import Path
from typing import Any

from pydantic import BaseModel, Field

from polyzymd.builders.conjugation.exceptions import (
    ConjugationNotImplementedError,
    PabloIngestionError,
)

SUPPORTED_STRUCTURE_SUFFIXES = frozenset({".pdb", ".cif", ".mmcif", ".pdbx"})


class PabloAvailability(BaseModel):
    """OpenFF Pablo import availability details."""

    available: bool
    version: str | None = None
    module_path: str | None = None
    error: str | None = None
    warnings: list[str] = Field(default_factory=list)


class PabloStructurePreflight(BaseModel):
    """Preflight diagnostics for a structure intended for Pablo ingestion."""

    intended_mode: str
    path: Path
    suffix: str
    pablo: PabloAvailability
    inspection_attempted: bool = False
    inspection_implemented: bool = False
    ingestion_implemented: bool = False
    warnings: list[str] = Field(default_factory=list)


class PabloIngestor:
    """Adapter boundary for future CCD/Pablo ingestion.

    The class intentionally keeps OpenFF Pablo imports inside methods so that
    importing PolyzyMD does not require the full chemistry stack.
    """

    def __init__(self, policy: Any) -> None:
        """Initialize a Pablo ingestion adapter.

        Parameters
        ----------
        policy : Any
            CCD/Pablo policy config from :class:`ConjugationConfig`.
        """
        self._policy = policy

    def ingest_existing(self, pdb_path: Path | str | None) -> Any:
        """Ingest a pre-conjugated PDB with OpenFF Pablo.

        Parameters
        ----------
        pdb_path : Path, str, or None
            Path to a pre-conjugated PDB structure.

        Returns
        -------
        Any
            Future OpenFF topology or ingestion result.

        Raises
        ------
        ConjugationNotImplementedError
            Always raised in the Phase 0-1 skeleton.
        """
        if pdb_path is not None and Path(pdb_path).exists():
            self.preflight_structure(pdb_path)
        else:
            self.probe_available()

        raise ConjugationNotImplementedError(
            "OpenFF Pablo ingestion for conjugation mode 'ingest_existing' is not implemented "
            "in the current conjugation foundation. Use preflight_structure() for file and "
            "Pablo availability diagnostics while CCD-aware residue and crosslink wiring are "
            "implemented."
        )

    def check_available(self) -> PabloAvailability:
        """Import OpenFF Pablo lazily and return version/path details.

        Returns
        -------
        PabloAvailability
            Successful Pablo import details.

        Raises
        ------
        PabloIngestionError
            If OpenFF Pablo is not importable in the active environment.
        """
        availability = self.probe_available()
        if not availability.available:
            raise PabloIngestionError(
                "OpenFF Pablo is not importable in this environment. Use the "
                "'conjugation-py312' pixi environment for conjugation ingestion preflight. "
                f"Original import error: {availability.error}"
            )
        return availability

    def probe_available(self) -> PabloAvailability:
        """Probe OpenFF Pablo availability without raising on import failure.

        Returns
        -------
        PabloAvailability
            Pablo import details or a structured import error.
        """
        try:
            pablo_module = importlib.import_module("openff.pablo")
        except ImportError as exc:
            return PabloAvailability(
                available=False,
                error=str(exc),
                warnings=[
                    "OpenFF Pablo is unavailable; install/use the conjugation-py312 pixi "
                    "environment before production ingestion"
                ],
            )

        version = _detect_pablo_version(pablo_module)
        module_path = getattr(pablo_module, "__file__", None)
        warnings = [
            "Pablo import succeeded, but PolyzyMD has not implemented production ingestion yet"
        ]
        return PabloAvailability(
            available=True,
            version=version,
            module_path=str(module_path) if module_path is not None else None,
            warnings=warnings,
        )

    def preflight_structure(self, path: Path | str) -> PabloStructurePreflight:
        """Validate a structure path and report Pablo readiness for future ingestion.

        The preflight intentionally does not construct an OpenFF topology. Pablo's
        production ingestion API will be wired in a later chemistry phase after
        residue definitions, crosslinks, and charge policy are finalized.

        Parameters
        ----------
        path : Path or str
            Structure path intended for future Pablo ingestion.

        Returns
        -------
        PabloStructurePreflight
            File validation and Pablo availability diagnostics.

        Raises
        ------
        PabloIngestionError
            If the path is missing or has an unsupported suffix.
        """
        structure_path = Path(path)
        _validate_structure_path(structure_path)

        availability = self.probe_available()
        warnings = [
            "Pablo structure inspection is not wired yet; this report is preflight-only",
            *availability.warnings,
        ]
        return PabloStructurePreflight(
            intended_mode="ingest_existing",
            path=structure_path,
            suffix=structure_path.suffix.lower(),
            pablo=availability,
            warnings=warnings,
        )


def _validate_structure_path(path: Path) -> None:
    """Validate a structure path for Pablo preflight."""
    if not path.exists():
        raise PabloIngestionError(f"Structure file does not exist: {path}")
    if not path.is_file():
        raise PabloIngestionError(f"Structure path is not a file: {path}")
    if path.suffix.lower() not in SUPPORTED_STRUCTURE_SUFFIXES:
        supported = ", ".join(sorted(SUPPORTED_STRUCTURE_SUFFIXES))
        raise PabloIngestionError(
            f"Unsupported structure suffix '{path.suffix}'. Expected one of: {supported}"
        )


def _detect_pablo_version(pablo_module: Any) -> str | None:
    """Return the OpenFF Pablo version when available."""
    try:
        return importlib.metadata.version("openff-pablo")
    except importlib.metadata.PackageNotFoundError:
        return getattr(pablo_module, "__version__", None)
