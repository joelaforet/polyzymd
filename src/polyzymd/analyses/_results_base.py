"""Base classes for analysis results.

All analysis result types inherit from BaseAnalysisResult, which provides:
- Config hash validation
- JSON serialization/deserialization
- Standard metadata fields

Design Principles
-----------------
1. Results are immutable after creation
2. All results store the config hash for validation
3. Results can be saved/loaded from JSON
4. Timestamps and versions are tracked for reproducibility
"""

from __future__ import annotations

import json
from abc import ABC, abstractmethod
from datetime import datetime
from pathlib import Path
from typing import Any, ClassVar, Self

from pydantic import BaseModel, Field


class BaseAnalysisResult(BaseModel, ABC):
    """Base class for all analysis result types.

    All analysis results should inherit from this class to ensure
    consistent serialization, metadata tracking, and config validation.

    Serialization Contract
    ----------------------
    Subclasses get ``save()`` and ``load()`` for free.  ``save()`` writes
    the Pydantic model to JSON via ``model_dump(mode="json")``, and
    ``load()`` reconstructs it via ``model_validate()``.  This means:

    - All fields must be JSON-serializable (Pydantic handles ``datetime``,
      ``Path``, enums, and nested ``BaseModel`` objects automatically).
    - For NumPy arrays, convert to ``list`` in a ``field_serializer`` or
      store as a plain list field and convert at access time.
    - For large binary data (e.g., per-frame SASA arrays), use a separate
      NPZ + JSON sidecar pattern instead of inheriting from this class.

    Adding a New Result Type
    ------------------------
    1. **Inherit from** ``BaseAnalysisResult`` and set ``analysis_type``
       as a ``ClassVar[str]``:

       .. code-block:: python

           from typing import ClassVar
           from polyzymd.analysis.results.base import BaseAnalysisResult

           class MyResult(BaseAnalysisResult):
               analysis_type: ClassVar[str] = "my_analysis"
               # your fields here
               score: float
               labels: list[str]

               def summary(self) -> str:
                   return f"MyResult: score={self.score:.3f}"

    2. **Implement** ``summary()`` — a human-readable one-liner for logs
       and CLI output.

    3. **Do NOT reimplement** ``save()`` or ``load()`` — the inherited
       versions handle JSON serialization automatically.

    4. **Nested data objects** should inherit from ``pydantic.BaseModel``
       (not ``BaseAnalysisResult``).  Only top-level result containers
       need the full ``BaseAnalysisResult`` machinery.

    5. **Metadata fields are optional** — ``config_hash``, ``equilibration_time``,
       and ``selection_string`` have sensible defaults for analyzers that
       operate without full config context (e.g., low-level ``ContactAnalyzer``).

    6. **Backward compatibility** — if migrating from a legacy format, use a
       ``model_validator(mode="before")`` to remap old field names.  See
       ``ContactResult`` for an example mapping ``analysis_timestamp`` →
       ``created_at``.

    Attributes
    ----------
    analysis_type : str
        Type of analysis (e.g., "rmsf", "distances").
    config_hash : str
        Hash of config at time of analysis for validation.
        Defaults to ``"unknown"`` for contexts without config access.
    created_at : datetime
        Timestamp when result was created.
    polyzymd_version : str
        Version of PolyzyMD used for analysis.
    replicate : int | None
        Replicate number (None for aggregated results).
    equilibration_time : float
        Time skipped for equilibration. Defaults to ``0.0`` for
        low-level analyzers that operate without config context.
    equilibration_unit : str
        Unit of equilibration time (e.g., "ns", "ps").
    selection_string : str
        MDAnalysis selection string used. Defaults to ``""`` for
        contexts where selections are implicit.

    See Also
    --------
    BaseAnalyzer : Base class for analyzer implementations.
    AggregatedResultMixin : Mixin for multi-replicate aggregated results.
    """

    # Class variable - subclasses should override
    analysis_type: ClassVar[str] = "base"

    # Metadata
    config_hash: str = Field(
        default="unknown",
        description="SHA-256 hash of config for cache validation",
    )
    created_at: datetime = Field(
        default_factory=datetime.now, description="Timestamp of result creation"
    )
    polyzymd_version: str = Field(
        default="unknown", description="PolyzyMD version used for analysis"
    )

    # Analysis parameters
    replicate: int | None = Field(
        default=None, description="Replicate number (1-indexed), None for aggregated"
    )
    equilibration_time: float = Field(
        default=0.0,
        description="Time skipped for equilibration",
    )
    equilibration_unit: str = Field(default="ns", description="Unit of equilibration time")
    selection_string: str = Field(
        default="",
        description="MDAnalysis selection string used",
    )

    # Correlation time info (if autocorrelation was computed)
    correlation_time: float | None = Field(default=None, description="Estimated correlation time τ")
    correlation_time_unit: str | None = Field(default=None, description="Unit of correlation time")
    n_independent_frames: int | None = Field(
        default=None, description="Number of independent frames used"
    )

    model_config = {"extra": "forbid"}

    @classmethod
    def get_analysis_type(cls) -> str:
        """Get the analysis type for this result class."""
        return cls.analysis_type

    def save(self, filepath: str | Path) -> Path:
        """Save result to JSON file.

        Parameters
        ----------
        filepath : str or Path
            Output file path. Parent directories will be created.

        Returns
        -------
        Path
            Path to saved file
        """
        filepath = Path(filepath)
        filepath.parent.mkdir(parents=True, exist_ok=True)

        with open(filepath, "w") as f:
            json.dump(self.model_dump(mode="json"), f, indent=2, default=str)

        return filepath

    @classmethod
    def load(cls, filepath: str | Path) -> Self:
        """Load result from JSON file.

        Parameters
        ----------
        filepath : str or Path
            Path to JSON file

        Returns
        -------
        Self
            Loaded result instance
        """
        filepath = Path(filepath)
        with open(filepath) as f:
            data = json.load(f)

        return cls.model_validate(data)

    @abstractmethod
    def summary(self) -> str:
        """Return a human-readable summary of the result."""
        pass

    def _format_equilibration(self) -> str:
        """Format equilibration time for display."""
        return f"{self.equilibration_time}{self.equilibration_unit}"

    def _format_correlation_time(self) -> str:
        """Format correlation time for display."""
        if self.correlation_time is None:
            return "not computed"
        return f"{self.correlation_time:.2f}{self.correlation_time_unit}"


class AggregatedResultMixin:
    """Mixin for aggregated (multi-replicate) results.

    Provides common methods for results that combine data from multiple replicates.

    Note: Child classes must define `replicates: list[int]` and `n_replicates: int`
    fields themselves. This mixin only provides utility methods.
    """

    # These are expected to be defined by child classes as Pydantic fields.
    # We declare them here only for type checking purposes (not as Fields).
    replicates: list[int]
    n_replicates: int

    @property
    def replicate_range(self) -> str:
        """Format replicate list as range string."""
        reps = sorted(self.replicates)  # type: ignore
        if len(reps) == 0:
            return "none"
        if len(reps) == 1:
            return str(reps[0])

        # Check if consecutive
        if reps == list(range(reps[0], reps[-1] + 1)):
            return f"{reps[0]}-{reps[-1]}"

        # Non-consecutive, list them
        return ",".join(map(str, reps))


def get_polyzymd_version() -> str:
    """Get current PolyzyMD version."""
    try:
        from importlib.metadata import version

        return version("polyzymd")
    except Exception:
        return "unknown"
