"""Surface exposure filtering using SASA calculations.

This module provides tools for determining which protein residues are
surface-exposed based on solvent-accessible surface area (SASA) calculations.
Surface exposure filtering is critical for binding preference analysis because
buried residues cannot participate in polymer-protein contacts.

Uses `rust_sasa_python` for fast SASA computation on initial PDB structures.

The surface exposure threshold determines which residues are considered
"exposed" based on their relative SASA (actual SASA / max theoretical SASA).
A threshold of 0.2 means residues with at least 20% of their theoretical
maximum surface area exposed are considered surface-accessible.

Examples
--------
>>> from polyzymd.analysis.contacts.surface_exposure import SurfaceExposureFilter
>>> filter = SurfaceExposureFilter(threshold=0.2)
>>> result = filter.calculate("enzyme.pdb")
>>> print(f"Found {result.exposed_count} surface-exposed residues")
>>> print(result.exposed_by_aa_class())
{'aromatic': 12, 'polar': 25, 'nonpolar': 18, ...}
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from pathlib import Path
from typing import TYPE_CHECKING

from polyzymd.analysis.common.aa_classification import MAX_ASA_TABLE, get_aa_class

if TYPE_CHECKING:
    pass

logger = logging.getLogger(__name__)


@dataclass
class ResidueExposure:
    """SASA data for a single protein residue.

    Attributes
    ----------
    resid : int
        Residue ID (1-indexed, as in PDB)
    resname : str
        3-letter residue name (e.g., "ALA", "PHE")
    chain_id : str
        Chain identifier from PDB
    sasa : float
        Calculated SASA in Angstrom^2
    max_sasa : float
        Maximum theoretical SASA for this residue type (Tien et al. 2013)
    relative_sasa : float
        Ratio of actual to maximum SASA (0.0 to ~1.0+)
    is_exposed : bool
        True if relative_sasa > threshold
    aa_class : str
        Amino acid classification (aromatic, polar, etc.)
    """

    resid: int
    resname: str
    chain_id: str
    sasa: float
    max_sasa: float
    relative_sasa: float
    is_exposed: bool
    aa_class: str


@dataclass
class SurfaceExposureResult:
    """Surface exposure analysis result for a protein structure.

    Contains SASA data for all residues and provides convenience methods
    for filtering and grouping by exposure status and amino acid class.

    Attributes
    ----------
    residue_exposures : list[ResidueExposure]
        SASA data for each residue
    threshold : float
        Relative SASA threshold used (e.g., 0.2 for 20%)
    pdb_path : str
        Path to the PDB file analyzed
    probe_radius : float
        Probe radius used for SASA calculation
    n_points : int
        Number of points used for SASA calculation
    """

    residue_exposures: list[ResidueExposure] = field(default_factory=list)
    threshold: float = 0.2
    pdb_path: str = ""
    probe_radius: float = 1.4
    n_points: int = 1000

    @property
    def exposed_resids(self) -> set[int]:
        """Set of residue IDs that are surface-exposed."""
        return {r.resid for r in self.residue_exposures if r.is_exposed}

    @property
    def buried_resids(self) -> set[int]:
        """Set of residue IDs that are buried (not surface-exposed)."""
        return {r.resid for r in self.residue_exposures if not r.is_exposed}

    @property
    def all_resids(self) -> set[int]:
        """Set of all residue IDs."""
        return {r.resid for r in self.residue_exposures}

    @property
    def exposed_count(self) -> int:
        """Number of surface-exposed residues."""
        return len(self.exposed_resids)

    @property
    def total_count(self) -> int:
        """Total number of residues analyzed."""
        return len(self.residue_exposures)

    def exposed_by_aa_class(self) -> dict[str, int]:
        """Count of exposed residues per amino acid class.

        Returns
        -------
        dict[str, int]
            Mapping from AA class to count of exposed residues.
            E.g., {'aromatic': 12, 'polar': 25, ...}
        """
        counts: dict[str, int] = {}
        for r in self.residue_exposures:
            if r.is_exposed:
                counts[r.aa_class] = counts.get(r.aa_class, 0) + 1
        return counts

    def total_by_aa_class(self) -> dict[str, int]:
        """Count of all residues per amino acid class (exposed + buried).

        Returns
        -------
        dict[str, int]
            Mapping from AA class to total count.
        """
        counts: dict[str, int] = {}
        for r in self.residue_exposures:
            counts[r.aa_class] = counts.get(r.aa_class, 0) + 1
        return counts

    def exposed_in_selection(self, resids: set[int]) -> set[int]:
        """Get exposed residues within a specific selection.

        Parameters
        ----------
        resids : set[int]
            Set of residue IDs to filter

        Returns
        -------
        set[int]
            Intersection of resids with exposed_resids
        """
        return self.exposed_resids & resids

    def get_exposure(self, resid: int) -> ResidueExposure | None:
        """Get exposure data for a specific residue.

        Parameters
        ----------
        resid : int
            Residue ID to look up

        Returns
        -------
        ResidueExposure or None
            Exposure data, or None if resid not found
        """
        for r in self.residue_exposures:
            if r.resid == resid:
                return r
        return None

    def to_dict(self) -> dict:
        """Serialize to dictionary for JSON storage."""
        return {
            "threshold": self.threshold,
            "pdb_path": self.pdb_path,
            "probe_radius": self.probe_radius,
            "n_points": self.n_points,
            "residue_exposures": [
                {
                    "resid": r.resid,
                    "resname": r.resname,
                    "chain_id": r.chain_id,
                    "sasa": r.sasa,
                    "max_sasa": r.max_sasa,
                    "relative_sasa": r.relative_sasa,
                    "is_exposed": r.is_exposed,
                    "aa_class": r.aa_class,
                }
                for r in self.residue_exposures
            ],
        }

    @classmethod
    def from_dict(cls, data: dict) -> "SurfaceExposureResult":
        """Deserialize from dictionary."""
        exposures = [
            ResidueExposure(
                resid=r["resid"],
                resname=r["resname"],
                chain_id=r["chain_id"],
                sasa=r["sasa"],
                max_sasa=r["max_sasa"],
                relative_sasa=r["relative_sasa"],
                is_exposed=r["is_exposed"],
                aa_class=r["aa_class"],
            )
            for r in data.get("residue_exposures", [])
        ]
        return cls(
            residue_exposures=exposures,
            threshold=data.get("threshold", 0.2),
            pdb_path=data.get("pdb_path", ""),
            probe_radius=data.get("probe_radius", 1.4),
            n_points=data.get("n_points", 1000),
        )


class SurfaceExposureFilter:
    """Compute surface exposure from initial PDB structure.

    Uses `rust_sasa_python` for fast SASA calculation. Residues are classified
    as "exposed" if their relative SASA (actual/max theoretical) exceeds
    the threshold.

    Parameters
    ----------
    threshold : float
        Relative SASA threshold (0-1). Residues with
        SASA/maxSASA > threshold are considered exposed.
        Default 0.2 (20% of max theoretical).
    probe_radius : float
        Probe radius for SASA calculation in Angstroms.
        Default 1.4 (water molecule radius).
    n_points : int
        Number of points for SASA calculation.
        Higher = more accurate but slower.
        Default 1000.

    Examples
    --------
    >>> filter = SurfaceExposureFilter(threshold=0.2)
    >>> result = filter.calculate("enzyme.pdb")
    >>> exposed = result.exposed_resids
    >>> print(f"Found {len(exposed)} exposed residues")

    >>> # Get exposed residues by AA class
    >>> by_class = result.exposed_by_aa_class()
    >>> print(f"Exposed aromatics: {by_class.get('aromatic', 0)}")
    """

    def __init__(
        self,
        threshold: float = 0.2,
        probe_radius: float = 1.4,
        n_points: int = 1000,
    ):
        if not 0.0 < threshold < 1.0:
            logger.warning(
                f"Surface exposure threshold {threshold} is outside typical range (0-1). "
                "This may classify too few or too many residues as exposed."
            )
        self.threshold = threshold
        self.probe_radius = probe_radius
        self.n_points = n_points

    def calculate(self, pdb_path: Path | str) -> SurfaceExposureResult:
        """Calculate surface exposure for all protein residues.

        Parameters
        ----------
        pdb_path : Path or str
            Path to PDB file (typically the enzyme PDB from config).
            Should contain only the protein (not polymer/solvent).

        Returns
        -------
        SurfaceExposureResult
            Surface exposure data for all residues

        Raises
        ------
        FileNotFoundError
            If PDB file does not exist
        ImportError
            If rust_sasa_python is not installed
        """
        # Lazy import to avoid dependency issues
        try:
            import rust_sasa_python as sasa
        except ImportError as e:
            raise ImportError(
                "rust_sasa_python is required for surface exposure analysis. "
                "Install with: pip install rust-sasa-python"
            ) from e

        pdb_path = Path(pdb_path)
        if not pdb_path.exists():
            raise FileNotFoundError(f"PDB file not found: {pdb_path}")

        logger.info(f"Calculating SASA for {pdb_path.name} (threshold={self.threshold})")

        # Create calculator and configure
        calculator = sasa.SASACalculator(str(pdb_path))
        calculator = calculator.with_probe_radius(self.probe_radius)
        calculator = calculator.with_n_points(self.n_points)

        # Calculate per-residue SASA
        residue_sasa_values = calculator.calculate_residue()

        exposures = []
        for res_data in residue_sasa_values:
            resname = res_data.residue_name
            resid = res_data.residue_number
            chain_id = getattr(res_data, "chain_id", "A")  # Default to A if not present
            sasa_value = res_data.sasa

            # Get max SASA from table, fallback to 200 A^2 for unknown residues
            max_sasa = MAX_ASA_TABLE.get(resname.upper(), 200.0)

            # Calculate relative SASA
            relative_sasa = sasa_value / max_sasa if max_sasa > 0 else 0.0

            # Determine if exposed
            is_exposed = relative_sasa > self.threshold

            exposures.append(
                ResidueExposure(
                    resid=resid,
                    resname=resname,
                    chain_id=chain_id,
                    sasa=sasa_value,
                    max_sasa=max_sasa,
                    relative_sasa=relative_sasa,
                    is_exposed=is_exposed,
                    aa_class=get_aa_class(resname),
                )
            )

        result = SurfaceExposureResult(
            residue_exposures=exposures,
            threshold=self.threshold,
            pdb_path=str(pdb_path),
            probe_radius=self.probe_radius,
            n_points=self.n_points,
        )

        logger.info(
            f"SASA analysis complete: {result.exposed_count}/{result.total_count} "
            f"residues are surface-exposed (>{self.threshold * 100:.0f}% relative SASA)"
        )

        return result
