"""Pydantic configuration model for exposure dynamics analysis."""

from __future__ import annotations

from pydantic import BaseModel, Field, model_validator
from typing_extensions import Self


class ExposureConfig(BaseModel):
    """Configuration for residue exposure dynamics analysis.

    Residues are classified into three stability categories based on the
    fraction of frames in which they are exposed (rSASA > exposure_threshold):

    - ``stably_exposed``: exposed in >= ``transient_upper`` of frames
    - ``stably_buried``: exposed in <= ``transient_lower`` of frames
    - ``transient``: exposure fraction falls between the two thresholds

    Parameters
    ----------
    transient_lower : float
        Lower bound of transient exposure fraction.  Residues with
        exposure_fraction <= transient_lower are classified as
        "stably_buried".  Default 0.2 (20 %).
    transient_upper : float
        Upper bound of transient exposure fraction.  Residues with
        exposure_fraction >= transient_upper are classified as
        "stably_exposed".  Default 0.8 (80 %).
    min_event_length : int
        Minimum number of consecutive frames required for a chaperone or
        refolding event to be counted.  Filters out single-frame noise.
        Default 1 (all events counted).
    polymer_resnames : list[str]
        Residue names of polymer monomers used to look up contacts.
        E.g. ["SBM", "EGM"].  If empty, all polymer types are included.
    """

    transient_lower: float = Field(
        default=0.2,
        ge=0.0,
        lt=1.0,
        description="Exposure fraction threshold below which residue is stably buried",
    )
    transient_upper: float = Field(
        default=0.8,
        gt=0.0,
        le=1.0,
        description="Exposure fraction threshold above which residue is stably exposed",
    )
    min_event_length: int = Field(
        default=1,
        ge=1,
        description="Minimum frames for an event to be counted",
    )
    polymer_resnames: list[str] = Field(
        default_factory=list,
        description="Polymer monomer resnames to include (empty = all)",
    )

    @model_validator(mode="after")
    def validate_thresholds(self) -> Self:
        if self.transient_lower >= self.transient_upper:
            raise ValueError(
                f"transient_lower ({self.transient_lower}) must be less than "
                f"transient_upper ({self.transient_upper})"
            )
        return self
