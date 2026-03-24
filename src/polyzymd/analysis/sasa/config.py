"""Pydantic configuration model for SASA trajectory analysis."""

from __future__ import annotations

from pydantic import BaseModel, Field


class SASAConfig(BaseModel):
    """Configuration for trajectory SASA computation.

    Parameters
    ----------
    exposure_threshold : float
        Relative SASA threshold for classifying a residue as "exposed".
        A residue is exposed if (SASA / max_SASA) > threshold.
        Default 0.2 (20% of theoretical maximum).
    probe_radius_nm : float
        Probe radius in nanometers for shrake_rupley. Default 0.14 nm (water).
    n_sphere_points : int
        Number of sphere points for shrake_rupley. Higher = more accurate but slower.
        Default 960 (MDTraj default).
    chain_id : str
        Chain ID to extract as protein. Default "A" (chain convention: A=protein).
    cache_sasa : bool
        Whether to cache SASA arrays to disk to avoid recomputation.
        Default True.
    """

    exposure_threshold: float = Field(
        default=0.2,
        gt=0.0,
        lt=1.0,
        description="Relative SASA threshold to classify residue as exposed (0-1)",
    )
    probe_radius_nm: float = Field(
        default=0.14,
        gt=0.0,
        description="Probe radius in nanometers",
    )
    n_sphere_points: int = Field(
        default=960,
        ge=100,
        description="Number of sphere points for SASA calculation",
    )
    chain_id: str = Field(
        default="A",
        description="Chain ID of the protein (chain convention: A=protein)",
    )
    cache_sasa: bool = Field(
        default=True,
        description="Cache SASA arrays to disk",
    )
