"""Configuration schema for comparison projects.

This module defines the YAML schema for comparison.yaml files that
specify which simulation conditions to compare.

Also includes catalytic triad/active site configuration for distance
and H-bond analysis.
"""

from __future__ import annotations

from datetime import datetime
from pathlib import Path
from typing import Optional

import yaml
from pydantic import BaseModel, Field, field_validator


# =============================================================================
# Catalytic Triad / Active Site Configuration
# =============================================================================


class TriadPairConfig(BaseModel):
    """Configuration for one distance pair in a catalytic triad/active site.

    Each pair defines two atom selections between which distances will be
    measured. Selections support special syntax:
    - Standard MDAnalysis: "resid 77 and name OG"
    - Midpoint: "midpoint(resid 133 and name OD1 OD2)"
    - Center of mass: "com(resid 50-75)"

    Attributes
    ----------
    label : str
        Human-readable label for this pair (e.g., "Asp133-His156")
    selection_a : str
        First atom/point selection
    selection_b : str
        Second atom/point selection

    Examples
    --------
    >>> pair = TriadPairConfig(
    ...     label="Asp133-His156",
    ...     selection_a="midpoint(resid 133 and name OD1 OD2)",
    ...     selection_b="resid 156 and name ND1",
    ... )
    """

    label: str = Field(..., description="Human-readable label for this pair")
    selection_a: str = Field(..., description="First atom/point selection")
    selection_b: str = Field(..., description="Second atom/point selection")


class CatalyticTriadConfig(BaseModel):
    """Configuration for catalytic triad/active site distance analysis.

    Defines multiple distance pairs that together comprise an active site
    or catalytic machinery. Computes:
    - Per-pair distance distributions and statistics
    - Simultaneous contact fraction (all pairs below threshold at same time)

    Attributes
    ----------
    name : str
        Name of the triad/active site (e.g., "LipA_catalytic_triad")
    pairs : list[TriadPairConfig]
        Distance pairs to monitor
    threshold : float
        Distance threshold for contact/H-bond analysis (Angstroms)
    description : str, optional
        Description of the active site

    Examples
    --------
    In comparison.yaml:

    ```yaml
    catalytic_triad:
      name: "LipA_catalytic_triad"
      description: "Ser-His-Asp catalytic triad of Lipase A"
      threshold: 3.5
      pairs:
        - label: "Asp133-His156"
          selection_a: "midpoint(resid 133 and name OD1 OD2)"
          selection_b: "resid 156 and name ND1"
        - label: "His156-Ser77"
          selection_a: "resid 156 and name NE2"
          selection_b: "resid 77 and name OG"
    ```
    """

    name: str = Field(..., description="Name of the catalytic triad/active site")
    pairs: list[TriadPairConfig] = Field(..., description="Distance pairs to monitor")
    threshold: float = Field(
        default=3.5, description="Distance threshold for contact analysis (Angstroms)"
    )
    description: Optional[str] = Field(default=None, description="Description of the active site")

    @field_validator("pairs", mode="after")
    @classmethod
    def validate_pairs(cls, v: list[TriadPairConfig]) -> list[TriadPairConfig]:
        """Ensure at least one pair is defined."""
        if len(v) == 0:
            raise ValueError("At least one distance pair must be defined")
        return v

    @property
    def n_pairs(self) -> int:
        """Number of distance pairs."""
        return len(self.pairs)

    def get_pair_selections(self) -> list[tuple[str, str]]:
        """Get list of (selection_a, selection_b) tuples for DistanceCalculator."""
        return [(p.selection_a, p.selection_b) for p in self.pairs]

    def get_pair_labels(self) -> list[str]:
        """Get list of pair labels."""
        return [p.label for p in self.pairs]


# =============================================================================
# Comparison Configuration
# =============================================================================


class ConditionConfig(BaseModel):
    """Configuration for one condition in a comparison.

    Attributes
    ----------
    label : str
        Display name for this condition (e.g., "No Polymer", "100% SBMA")
    config : Path
        Path to the simulation's config.yaml file
    replicates : list[int]
        List of replicate numbers to include in the analysis
    """

    label: str
    config: Path
    replicates: list[int]

    @field_validator("config", mode="before")
    @classmethod
    def resolve_path(cls, v: str | Path) -> Path:
        """Convert string paths to Path objects."""
        return Path(v)

    @field_validator("replicates", mode="before")
    @classmethod
    def ensure_list(cls, v: list[int] | int) -> list[int]:
        """Ensure replicates is a list."""
        if isinstance(v, int):
            return [v]
        return list(v)


class AnalysisDefaults(BaseModel):
    """Default parameters for comparison analyses.

    These can be overridden on the command line.

    Attributes
    ----------
    equilibration_time : str
        Time to skip for equilibration (e.g., "10ns", "5000ps")
    selection : str
        MDAnalysis selection string for RMSF calculation
    reference_mode : str
        Reference structure mode for alignment
    """

    equilibration_time: str = "10ns"
    selection: str = "protein and name CA"
    reference_mode: str = "centroid"


class ComparisonConfig(BaseModel):
    """Schema for comparison.yaml configuration files.

    A comparison config defines multiple simulation conditions to compare,
    along with default analysis parameters and optional catalytic triad
    configuration.

    Attributes
    ----------
    name : str
        Name of the comparison project
    description : str, optional
        Description of what is being compared
    control : str, optional
        Label of the control condition for relative comparisons
    conditions : list[ConditionConfig]
        List of conditions to compare
    defaults : AnalysisDefaults
        Default analysis parameters
    catalytic_triad : CatalyticTriadConfig, optional
        Configuration for catalytic triad/active site analysis

    Examples
    --------
    >>> config = ComparisonConfig.from_yaml("comparison.yaml")
    >>> print(config.name)
    "Polymer Stabilization Study"
    >>> for cond in config.conditions:
    ...     print(f"{cond.label}: {cond.config}")
    >>> if config.catalytic_triad:
    ...     print(f"Triad: {config.catalytic_triad.name}")
    """

    name: str
    description: Optional[str] = None
    control: Optional[str] = None
    conditions: list[ConditionConfig]
    defaults: AnalysisDefaults = AnalysisDefaults()
    catalytic_triad: Optional[CatalyticTriadConfig] = None

    @classmethod
    def from_yaml(cls, path: Path | str) -> "ComparisonConfig":
        """Load comparison config from YAML file.

        Parameters
        ----------
        path : Path or str
            Path to comparison.yaml file

        Returns
        -------
        ComparisonConfig
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
            raise FileNotFoundError(f"Comparison config not found: {path}")

        with open(path) as f:
            data = yaml.safe_load(f)

        # Resolve relative paths relative to the config file location
        config_dir = path.parent.resolve()
        if "conditions" in data:
            for cond in data["conditions"]:
                if "config" in cond:
                    cond_path = Path(cond["config"])
                    if not cond_path.is_absolute():
                        cond["config"] = str(config_dir / cond_path)

        return cls(**data)

    def to_yaml(self, path: Path | str) -> None:
        """Save comparison config to YAML file.

        Parameters
        ----------
        path : Path or str
            Output path for comparison.yaml
        """
        path = Path(path)

        # Convert to dict, handling Path objects
        data = self.model_dump()
        for cond in data["conditions"]:
            cond["config"] = str(cond["config"])

        with open(path, "w") as f:
            yaml.dump(data, f, default_flow_style=False, sort_keys=False)

    def get_condition(self, label: str) -> ConditionConfig:
        """Get a condition by its label.

        Parameters
        ----------
        label : str
            The condition label to find

        Returns
        -------
        ConditionConfig
            The matching condition

        Raises
        ------
        KeyError
            If no condition with that label exists
        """
        for cond in self.conditions:
            if cond.label == label:
                return cond
        raise KeyError(f"Condition '{label}' not found in: {[c.label for c in self.conditions]}")

    def validate_config(self) -> list[str]:
        """Validate the comparison configuration.

        Returns
        -------
        list[str]
            List of error messages (empty if valid)
        """
        errors = []

        # Check minimum conditions
        if len(self.conditions) < 2:
            errors.append("Need at least 2 conditions to compare")

        # Check for duplicate labels
        labels = [c.label for c in self.conditions]
        if len(labels) != len(set(labels)):
            errors.append("Duplicate condition labels found")

        # Check control label exists
        if self.control and self.control not in labels:
            errors.append(f"Control '{self.control}' not in conditions: {labels}")

        # Check config files exist
        for cond in self.conditions:
            if not cond.config.exists():
                errors.append(f"Config not found for '{cond.label}': {cond.config}")

        return errors


def generate_comparison_template(name: str, eq_time: str = "10ns") -> str:
    """Generate a template comparison.yaml file.

    Parameters
    ----------
    name : str
        Project name
    eq_time : str
        Default equilibration time

    Returns
    -------
    str
        YAML template content
    """
    return f'''\
# ============================================================================
# PolyzyMD Comparison Configuration
# ============================================================================
# Compare analyses across simulation conditions (e.g., polymer, temperature).
# Docs: https://polyzymd.readthedocs.io/en/latest/
# ============================================================================

name: "{name}"
description: "Comparison of simulation conditions"

# Control condition for relative comparisons (null if none)
control: null

# ============================================================================
# Conditions
# ============================================================================
# Each condition points to a simulation's config.yaml file.
conditions:
  - label: "Condition A"
    config: "../path/to/condition_a/config.yaml"
    replicates: [1, 2, 3]

  - label: "Condition B"
    config: "../path/to/condition_b/config.yaml"
    replicates: [1, 2, 3]

# ============================================================================
# Defaults
# ============================================================================
# Note: If equilibration_time is 0ns, you'll get a warning - set appropriately!
defaults:
  equilibration_time: "{eq_time}"
  selection: "protein and name CA"
  reference_mode: "centroid"

# ============================================================================
# Catalytic Triad (for polyzymd compare triad)
# ============================================================================
# Define your enzyme's catalytic machinery ONCE here - it's shared across
# all conditions since the enzyme is the same.
#
# Run: polyzymd compare triad -c comparison.yaml
#
# catalytic_triad:
#   name: "enzyme_catalytic_triad"
#   threshold: 3.5  # Angstroms (H-bond cutoff)
#   pairs:
#     - label: "Asp-His"
#       selection_a: "midpoint(resid 133 and name OD1 OD2)"
#       selection_b: "resid 156 and name ND1"
#     - label: "His-Ser"
#       selection_a: "resid 156 and name NE2"
#       selection_b: "resid 77 and name OG"
'''
