"""Configuration schema for analysis.yaml files.

This module defines the YAML schema for analysis.yaml files that
configure which analyses to run for a simulation.

The analysis.yaml must live alongside config.yaml (same directory)
to maintain the config.yaml as the single source of truth.
"""

from __future__ import annotations

from pathlib import Path
from typing import Optional

import yaml
from pydantic import BaseModel, Field, field_validator


# =============================================================================
# Individual Analysis Configurations
# =============================================================================


class RMSFConfig(BaseModel):
    """Configuration for RMSF analysis.

    Attributes
    ----------
    enabled : bool
        Whether to run RMSF analysis
    selection : str
        MDAnalysis selection string for atoms to analyze
    reference_mode : str
        Reference structure mode: 'centroid', 'average', or 'frame'
    reference_frame : int, optional
        Frame number if reference_mode is 'frame'
    """

    enabled: bool = True
    selection: str = "protein and name CA"
    reference_mode: str = "centroid"
    reference_frame: Optional[int] = None


class DistancePairConfig(BaseModel):
    """Configuration for a single distance pair.

    Attributes
    ----------
    label : str
        Human-readable label for this pair
    selection_a : str
        First atom/point selection
    selection_b : str
        Second atom/point selection
    """

    label: str
    selection_a: str
    selection_b: str


class DistancesConfig(BaseModel):
    """Configuration for distance analysis.

    Attributes
    ----------
    enabled : bool
        Whether to run distance analysis
    pairs : list[DistancePairConfig]
        List of atom pairs to measure distances between
    """

    enabled: bool = False
    pairs: list[DistancePairConfig] = Field(default_factory=list)


class TriadPairConfig(BaseModel):
    """Configuration for a catalytic triad distance pair.

    Attributes
    ----------
    label : str
        Human-readable label for this pair (e.g., "Asp133-His156")
    selection_a : str
        First atom/point selection
    selection_b : str
        Second atom/point selection
    """

    label: str
    selection_a: str
    selection_b: str


class CatalyticTriadConfig(BaseModel):
    """Configuration for catalytic triad analysis.

    Attributes
    ----------
    enabled : bool
        Whether to run triad analysis
    name : str
        Name of the catalytic triad/active site
    threshold : float
        Distance threshold for contact/H-bond analysis (Angstroms)
    pairs : list[TriadPairConfig]
        List of atom pairs that define the triad
    """

    enabled: bool = False
    name: str = "catalytic_triad"
    threshold: float = 3.5
    pairs: list[TriadPairConfig] = Field(default_factory=list)

    @field_validator("pairs", mode="after")
    @classmethod
    def validate_pairs_if_enabled(cls, v: list[TriadPairConfig], info) -> list[TriadPairConfig]:
        """Warn if enabled but no pairs defined."""
        # Note: We can't access 'enabled' here easily, so validation
        # happens at runtime in the CLI
        return v


# =============================================================================
# Main Analysis Configuration
# =============================================================================


class AnalysisDefaults(BaseModel):
    """Default parameters applied to all analyses.

    Attributes
    ----------
    equilibration_time : str
        Time to skip for equilibration (e.g., "10ns", "5000ps")
    """

    equilibration_time: str = "10ns"


class AnalysisConfig(BaseModel):
    """Schema for analysis.yaml configuration files.

    This configuration defines which analyses to run for a simulation.
    It must be located in the same directory as the simulation's config.yaml.

    Attributes
    ----------
    replicates : list[int]
        List of replicate numbers to analyze
    defaults : AnalysisDefaults
        Default parameters for all analyses
    rmsf : RMSFConfig
        RMSF analysis configuration
    distances : DistancesConfig
        Distance analysis configuration
    catalytic_triad : CatalyticTriadConfig
        Catalytic triad analysis configuration

    Examples
    --------
    >>> config = AnalysisConfig.from_yaml("analysis.yaml")
    >>> if config.rmsf.enabled:
    ...     print("RMSF analysis enabled")
    """

    replicates: list[int] = Field(default_factory=lambda: [1, 2, 3])
    defaults: AnalysisDefaults = Field(default_factory=AnalysisDefaults)
    rmsf: RMSFConfig = Field(default_factory=RMSFConfig)
    distances: DistancesConfig = Field(default_factory=DistancesConfig)
    catalytic_triad: CatalyticTriadConfig = Field(default_factory=CatalyticTriadConfig)

    @field_validator("replicates", mode="before")
    @classmethod
    def ensure_list(cls, v: list[int] | int) -> list[int]:
        """Ensure replicates is a list."""
        if isinstance(v, int):
            return [v]
        return list(v)

    @classmethod
    def from_yaml(cls, path: Path | str) -> "AnalysisConfig":
        """Load analysis config from YAML file.

        Parameters
        ----------
        path : Path or str
            Path to analysis.yaml file

        Returns
        -------
        AnalysisConfig
            Loaded and validated configuration

        Raises
        ------
        FileNotFoundError
            If the config file doesn't exist
        ValidationError
            If the config is invalid
        """
        path = Path(path)
        if not path.exists():
            raise FileNotFoundError(f"Analysis config not found: {path}")

        with open(path) as f:
            data = yaml.safe_load(f)

        # Handle empty file
        if data is None:
            data = {}

        return cls(**data)

    def to_yaml(self, path: Path | str) -> None:
        """Save analysis config to YAML file.

        Parameters
        ----------
        path : Path or str
            Output path for analysis.yaml
        """
        path = Path(path)
        data = self.model_dump(exclude_none=True)

        with open(path, "w") as f:
            yaml.dump(data, f, default_flow_style=False, sort_keys=False)

    def get_enabled_analyses(self) -> list[str]:
        """Get list of enabled analysis types.

        Returns
        -------
        list[str]
            Names of enabled analyses (e.g., ['rmsf', 'catalytic_triad'])
        """
        enabled = []
        if self.rmsf.enabled:
            enabled.append("rmsf")
        if self.distances.enabled:
            enabled.append("distances")
        if self.catalytic_triad.enabled:
            enabled.append("catalytic_triad")
        return enabled

    def validate_config(self) -> list[str]:
        """Validate the analysis configuration.

        Returns
        -------
        list[str]
            List of error/warning messages (empty if valid)
        """
        issues = []

        # Check replicates
        if not self.replicates:
            issues.append("No replicates specified")

        # Check distances pairs if enabled
        if self.distances.enabled and not self.distances.pairs:
            issues.append("Distance analysis enabled but no pairs defined")

        # Check triad pairs if enabled
        if self.catalytic_triad.enabled and not self.catalytic_triad.pairs:
            issues.append("Catalytic triad analysis enabled but no pairs defined")

        return issues


def generate_analysis_template(eq_time: str = "10ns") -> str:
    """Generate a template analysis.yaml file.

    Parameters
    ----------
    eq_time : str
        Default equilibration time

    Returns
    -------
    str
        YAML template content
    """
    return f'''\
# ============================================================================
# PolyzyMD Analysis Configuration
# ============================================================================
# Configures which analyses to run for this simulation.
# Must be in same directory as config.yaml.
#
# Run all enabled analyses: polyzymd analyze run
# Docs: https://polyzymd.readthedocs.io/en/latest/
# ============================================================================

# Which replicates to analyze
replicates: [1, 2, 3]

# Default parameters (override per-analysis if needed)
# Note: Set equilibration_time appropriately for your simulation!
defaults:
  equilibration_time: "{eq_time}"

# ============================================================================
# RMSF Analysis
# ============================================================================
rmsf:
  enabled: true
  selection: "protein and name CA"
  reference_mode: "centroid"  # centroid, average, or frame

# ============================================================================
# Distance Analysis
# ============================================================================
# distances:
#   enabled: true
#   pairs:
#     - label: "Ser77-Substrate"
#       selection_a: "resid 77 and name OG"
#       selection_b: "resname RBY and name C1"

# ============================================================================
# Catalytic Triad Analysis
# ============================================================================
# For cross-condition comparison, define triad in comparison.yaml instead.
#
# catalytic_triad:
#   enabled: true
#   name: "LipA_catalytic_triad"
#   threshold: 3.5
#   pairs:
#     - label: "Asp133-His156"
#       selection_a: "midpoint(resid 133 and name OD1 OD2)"
#       selection_b: "resid 156 and name ND1"
#     - label: "His156-Ser77"
#       selection_a: "resid 156 and name NE2"
#       selection_b: "resid 77 and name OG"
'''
