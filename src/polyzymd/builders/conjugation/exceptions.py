"""Exceptions for covalent modification builders."""

from __future__ import annotations


class ConjugationError(RuntimeError):
    """Base class for covalent modification workflow errors."""


class ConjugationNotImplementedError(ConjugationError):
    """Raised when a requested covalent modification path is not implemented."""


class PabloIngestionError(ConjugationError):
    """Raised when CCD/Pablo ingestion cannot proceed."""
