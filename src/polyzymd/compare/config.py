"""Configuration schema for comparison projects.

This module defines the YAML schema for comparison.yaml files that
specify which simulation conditions to compare.

The schema has two main sections:
- analysis_settings: Defines WHAT analyses to run (shared across conditions)
- comparison_settings: Defines HOW to compare (statistical parameters)

Both sections use a registry-based approach for extensibility. New analysis
types can be added by registering with AnalysisSettingsRegistry and
ComparisonSettingsRegistry (see polyzymd.compare.settings).
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Any

import yaml
from pydantic import BaseModel, Field, field_validator, model_validator

from polyzymd.analysis.core.registry import (
    AnalysisSettingsRegistry,
    BaseAnalysisSettings,
    BaseComparisonSettings,
    ComparisonSettingsRegistry,
)

# Import settings to trigger registration
from polyzymd.compare.settings import (  # noqa: F401
    BindingFreeEnergyAnalysisSettings,
    BindingFreeEnergyComparisonSettings,
    CatalyticTriadAnalysisSettings,
    CatalyticTriadComparisonSettings,
    ContactsAnalysisSettings,
    ContactsComparisonSettings,
    DistancePairSettings,
    DistancesAnalysisSettings,
    DistancesComparisonSettings,
    RMSFAnalysisSettings,
    RMSFComparisonSettings,
    TriadPairSettings,
)

# Backward-compatible aliases for analysis module
# The analysis/triad module still imports these old names
CatalyticTriadConfig = CatalyticTriadAnalysisSettings
TriadPairConfig = TriadPairSettings

logger = logging.getLogger(__name__)


# ============================================================================
# Condition Configuration
# ============================================================================


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
        Time to skip for equilibration (e.g., "10ns", "5000ps").
        Shared across all analyses (RMSF, contacts, triad).
    """

    equilibration_time: str = "10ns"


# ============================================================================
# Dynamic Settings Containers
# ============================================================================


class AnalysisSettingsContainer(BaseModel):
    """Container for analysis settings (WHAT to analyze).

    Uses dynamic attribute access to support any registered analysis type
    without hardcoding field names.
    """

    model_config = {"extra": "allow"}

    def __init__(self, **data: Any):
        """Initialize with dynamic analysis settings.

        Parameters
        ----------
        **data : Any
            Analysis settings keyed by analysis type name.
        """
        # Parse each setting using the registry
        parsed_settings: dict[str, BaseAnalysisSettings] = {}
        for key, value in data.items():
            if value is None:
                continue
            key_lower = key.lower()
            if AnalysisSettingsRegistry.is_registered(key_lower):
                settings_class = AnalysisSettingsRegistry.get(key_lower)
                if isinstance(value, dict):
                    parsed_settings[key_lower] = settings_class(**value)
                elif isinstance(value, BaseAnalysisSettings):
                    parsed_settings[key_lower] = value
                else:
                    raise ValueError(
                        f"Invalid value for {key}: expected dict or {settings_class.__name__}"
                    )
            else:
                logger.warning(f"Unknown analysis type '{key}' - skipping")

        super().__init__(**parsed_settings)

    def get(self, analysis_type: str) -> BaseAnalysisSettings | None:
        """Get settings for a specific analysis type.

        Parameters
        ----------
        analysis_type : str
            Analysis type identifier (e.g., "rmsf", "contacts").

        Returns
        -------
        BaseAnalysisSettings or None
            Settings for the analysis type, or None if not configured.
        """
        return getattr(self, analysis_type.lower(), None)

    def get_enabled_analyses(self) -> list[str]:
        """Get list of enabled analysis types.

        Returns
        -------
        list[str]
            Names of configured analyses (presence implies enabled).

        Notes
        -----
        Uses actual model data from comparison.yaml rather than relying on
        a registry. This makes comparison.yaml the source of truth for which
        analyses are enabled.
        """
        # Use actual model data - if an analysis section exists and has a value, it's enabled
        return [key for key, value in self.model_dump().items() if value is not None]

    def to_analysis_yaml_dict(self, replicates: list[int], eq_time: str) -> dict[str, Any]:
        """Convert to analysis.yaml-compatible dictionary.

        Parameters
        ----------
        replicates : list[int]
            Replicate numbers for the analysis.yaml.
        eq_time : str
            Equilibration time for the analysis.yaml.

        Returns
        -------
        dict[str, Any]
            Dictionary suitable for writing to analysis.yaml.
        """
        result: dict[str, Any] = {
            "replicates": replicates,
            "defaults": {"equilibration_time": eq_time},
        }
        for analysis_type in self.get_enabled_analyses():
            settings = self.get(analysis_type)
            if settings is not None:
                result[analysis_type] = settings.to_analysis_yaml_dict()
        return result


class ComparisonSettingsContainer(BaseModel):
    """Container for comparison settings (HOW to compare).

    Uses dynamic attribute access to support any registered comparison type.
    Each analysis type in analysis_settings must have a corresponding entry
    here (can be empty dict) to enable comparison.
    """

    model_config = {"extra": "allow"}

    def __init__(self, **data: Any):
        """Initialize with dynamic comparison settings.

        Parameters
        ----------
        **data : Any
            Comparison settings keyed by analysis type name.
        """
        # Parse each setting using the registry
        parsed_settings: dict[str, BaseComparisonSettings] = {}
        for key, value in data.items():
            if value is None:
                continue
            key_lower = key.lower()
            if ComparisonSettingsRegistry.is_registered(key_lower):
                settings_class = ComparisonSettingsRegistry.get(key_lower)
                if isinstance(value, dict):
                    parsed_settings[key_lower] = settings_class(**value)
                elif isinstance(value, BaseComparisonSettings):
                    parsed_settings[key_lower] = value
                else:
                    raise ValueError(
                        f"Invalid value for {key}: expected dict or {settings_class.__name__}"
                    )
            else:
                logger.warning(f"Unknown comparison type '{key}' - skipping")

        super().__init__(**parsed_settings)

    def get(self, analysis_type: str) -> BaseComparisonSettings | None:
        """Get settings for a specific comparison type.

        Parameters
        ----------
        analysis_type : str
            Analysis type identifier (e.g., "rmsf", "contacts").

        Returns
        -------
        BaseComparisonSettings or None
            Comparison settings, or None if not configured.
        """
        return getattr(self, analysis_type.lower(), None)

    def get_enabled_comparisons(self) -> list[str]:
        """Get list of enabled comparison types.

        Returns
        -------
        list[str]
            Names of configured comparisons.
        """
        enabled = []
        for analysis_type in ComparisonSettingsRegistry.list_available():
            if self.get(analysis_type) is not None:
                enabled.append(analysis_type)
        return enabled


# ============================================================================
# Plot Settings Configuration
# ============================================================================


class RMSFPlotSettings(BaseModel):
    """RMSF-specific plot customization.

    Attributes
    ----------
    show_error : bool
        Show error bands/bars on plots (default True)
    highlight_residues : list[int]
        Residue numbers to highlight with vertical lines (e.g., active site)
    figsize_profile : tuple[float, float]
        Figure size for per-residue profile plots
    figsize_comparison : tuple[float, float]
        Figure size for bar comparison plots
    """

    show_error: bool = True
    highlight_residues: list[int] = Field(default_factory=list)
    figsize_profile: tuple[float, float] = (14, 4)
    figsize_comparison: tuple[float, float] = (8, 6)


class TriadPlotSettings(BaseModel):
    """Triad-specific plot customization.

    Attributes
    ----------
    generate_kde_panel : bool
        Generate multi-row KDE panel plot (default True)
    generate_bars : bool
        Generate grouped threshold bar chart (default True)
    generate_2d_kde : bool
        Generate 2D joint KDE plot (default False, more specialized)
    threshold_line_color : str
        Color for threshold vertical line
    kde_fill_alpha : float
        Transparency for KDE fill (0-1)
    figsize_kde_panel : tuple[float, float] | None
        Figure size for KDE panel (auto-calculated if None)
    figsize_bars : tuple[float, float]
        Figure size for bar chart
    """

    generate_kde_panel: bool = True
    generate_bars: bool = True
    generate_2d_kde: bool = False
    threshold_line_color: str = "red"
    kde_fill_alpha: float = 0.7
    figsize_kde_panel: tuple[float, float] | None = None
    figsize_bars: tuple[float, float] = (10, 6)


class DistancesPlotSettings(BaseModel):
    """Distance analysis plot customization.

    Attributes
    ----------
    show_threshold : bool
        Show threshold line on distribution plots
    use_kde : bool
        Use KDE instead of histogram for distributions
    figsize : tuple[float, float]
        Default figure size for distance plots
    """

    show_threshold: bool = True
    use_kde: bool = True
    figsize: tuple[float, float] = (10, 6)


class ContactsPlotSettings(BaseModel):
    """Contacts analysis plot customization.

    Attributes
    ----------
    figsize : tuple[float, float]
        Default figure size for contact plots
    generate_enrichment_heatmap : bool
        Generate binding preference enrichment heatmap (default True)
    generate_enrichment_bars : bool
        Generate binding preference bar charts (default True)
    figsize_enrichment_heatmap : tuple[float, float] | None
        Figure size for enrichment heatmap (auto-calculated if None)
    figsize_enrichment_bars : tuple[float, float]
        Figure size for enrichment bar charts
    enrichment_colormap : str
        Colormap for enrichment heatmap (diverging recommended)
    show_enrichment_error : bool
        Show error bars on enrichment bar charts (default True)
    generate_system_coverage_heatmap : bool
        Generate system coverage enrichment heatmap (default True)
    generate_system_coverage_bars : bool
        Generate system coverage bar charts (default True)
    figsize_system_coverage_heatmap : tuple[float, float] | None
        Figure size for system coverage heatmap (auto-calculated if None)
    figsize_system_coverage_bars : tuple[float, float]
        Figure size for system coverage bar charts
    show_system_coverage_error : bool
        Show error bars on system coverage bar charts (default True)
    generate_user_partition_bars : bool
        Generate user-defined partition bar charts (default True)
    figsize_user_partition_bars : tuple[float, float]
        Figure size for user-defined partition bar charts
    show_user_partition_error : bool
        Show error bars on user-defined partition bar charts (default True)
    """

    figsize: tuple[float, float] = (10, 8)
    generate_enrichment_heatmap: bool = True
    generate_enrichment_bars: bool = True
    figsize_enrichment_heatmap: tuple[float, float] | None = None
    figsize_enrichment_bars: tuple[float, float] = (10, 6)
    enrichment_colormap: str = "RdBu_r"  # Diverging: red=high, blue=low
    show_enrichment_error: bool = True

    # System coverage plot settings
    generate_system_coverage_heatmap: bool = True
    generate_system_coverage_bars: bool = True
    figsize_system_coverage_heatmap: tuple[float, float] | None = None
    figsize_system_coverage_bars: tuple[float, float] = (10, 6)
    show_system_coverage_error: bool = True

    # User-defined partition plot settings
    generate_user_partition_bars: bool = True
    figsize_user_partition_bars: tuple[float, float] = (10, 6)
    show_user_partition_error: bool = True


class BFEPlotSettings(BaseModel):
    """Binding free energy plot customization.

    Attributes
    ----------
    generate_heatmap : bool
        Generate ΔΔG heatmap (rows = AA groups, columns = conditions). Default True.
    generate_bars : bool
        Generate ΔΔG grouped bar chart (one bar per condition per AA group). Default True.
    figsize_heatmap : tuple[float, float] | None
        Figure size for ΔΔG heatmap (auto-calculated if None).
    figsize_bars : tuple[float, float]
        Figure size for ΔΔG bar charts.
    colormap : str
        Diverging colormap for heatmap (default "RdBu_r": red = avoidance, blue = preference).
    show_error_bars : bool
        Show SEM error bars on bar charts. Default True.
    annotate_heatmap : bool
        Annotate each heatmap cell with its ΔΔG value. Default True.
    """

    generate_heatmap: bool = True
    generate_bars: bool = True
    figsize_heatmap: tuple[float, float] | None = None
    figsize_bars: tuple[float, float] = (10, 6)
    colormap: str = "RdBu_r"
    show_error_bars: bool = True
    annotate_heatmap: bool = True


class PlotSettings(BaseModel):
    """Global plot settings for comparison.yaml.

    Controls plot generation for all analyses. Can be customized
    per-analysis type using nested settings.

    Attributes
    ----------
    output_dir : Path
        Directory for generated plots (relative to comparison.yaml)
    format : str
        Image format: "png", "pdf", or "svg"
    dpi : int
        Resolution for raster formats (PNG)
    style : str
        Plot style preset: "publication", "presentation", or "minimal"
    color_palette : str
        Seaborn/matplotlib color palette name

    Examples
    --------
    In comparison.yaml:

    .. code-block:: yaml

        plot_settings:
          output_dir: "figures/"
          format: "png"
          dpi: 300
          style: "publication"

          rmsf:
            highlight_residues: [77, 133, 156]

          triad:
            generate_2d_kde: true
    """

    output_dir: Path = Field(default=Path("figures/"))
    format: str = Field(default="png", pattern="^(png|pdf|svg)$")
    dpi: int = Field(default=300, ge=50, le=600)
    style: str = Field(default="publication", pattern="^(publication|presentation|minimal)$")
    color_palette: str = "tab10"

    # Per-analysis plot settings
    rmsf: RMSFPlotSettings = Field(default_factory=RMSFPlotSettings)
    triad: TriadPlotSettings = Field(default_factory=TriadPlotSettings)
    distances: DistancesPlotSettings = Field(default_factory=DistancesPlotSettings)
    contacts: ContactsPlotSettings = Field(default_factory=ContactsPlotSettings)
    binding_free_energy: BFEPlotSettings = Field(default_factory=BFEPlotSettings)

    @field_validator("output_dir", mode="before")
    @classmethod
    def resolve_output_dir(cls, v: str | Path) -> Path:
        """Convert string paths to Path objects."""
        return Path(v)


# ============================================================================
# Main Comparison Configuration
# ============================================================================


class ComparisonConfig(BaseModel):
    """Schema for comparison.yaml configuration files.

    A comparison config defines multiple simulation conditions to compare,
    along with analysis settings and comparison-specific parameters.

    The schema follows a three-section pattern:
    - analysis_settings: WHAT to analyze (shared across conditions)
    - comparison_settings: HOW to compare (statistical parameters)
    - plot_settings: HOW to visualize (plot customization)

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
        Default analysis parameters (equilibration_time)
    analysis_settings : AnalysisSettingsContainer
        Analysis parameters (WHAT to analyze)
    comparison_settings : ComparisonSettingsContainer
        Comparison parameters (HOW to compare)
    plot_settings : PlotSettings
        Plot customization (HOW to visualize)

    Examples
    --------
    >>> config = ComparisonConfig.from_yaml("comparison.yaml")
    >>> print(config.name)
    "Polymer Stabilization Study"
    >>> for cond in config.conditions:
    ...     print(f"{cond.label}: {cond.config}")
    >>> print("Enabled analyses:", config.analysis_settings.get_enabled_analyses())
    >>> rmsf_settings = config.analysis_settings.get("rmsf")
    >>> if rmsf_settings:
    ...     print(f"RMSF selection: {rmsf_settings.selection}")
    """

    name: str
    description: str | None = None
    control: str | None = None
    conditions: list[ConditionConfig]
    defaults: AnalysisDefaults = Field(default_factory=AnalysisDefaults)
    analysis_settings: AnalysisSettingsContainer = Field(default_factory=AnalysisSettingsContainer)
    comparison_settings: ComparisonSettingsContainer = Field(
        default_factory=ComparisonSettingsContainer
    )
    plot_settings: PlotSettings = Field(default_factory=PlotSettings)
    source_path: Path | None = Field(default=None, exclude=True)

    @field_validator("analysis_settings", mode="before")
    @classmethod
    def parse_analysis_settings(cls, v: Any) -> AnalysisSettingsContainer:
        """Parse analysis_settings from dict or container."""
        if v is None:
            return AnalysisSettingsContainer()
        if isinstance(v, dict):
            return AnalysisSettingsContainer(**v)
        return v

    @field_validator("comparison_settings", mode="before")
    @classmethod
    def parse_comparison_settings(cls, v: Any) -> ComparisonSettingsContainer:
        """Parse comparison_settings from dict or container."""
        if v is None:
            return ComparisonSettingsContainer()
        if isinstance(v, dict):
            return ComparisonSettingsContainer(**v)
        return v

    @model_validator(mode="after")
    def validate_comparison_coverage(self) -> "ComparisonConfig":
        """Validate that comparison_settings covers all analysis_settings.

        Each analysis type in analysis_settings must have a corresponding
        entry in comparison_settings (can be empty {}).
        """
        enabled_analyses = self.analysis_settings.get_enabled_analyses()
        enabled_comparisons = self.comparison_settings.get_enabled_comparisons()

        missing = set(enabled_analyses) - set(enabled_comparisons)
        if missing:
            raise ValueError(
                f"Missing comparison_settings for: {sorted(missing)}. "
                f"Add 'comparison_settings.{list(missing)[0]}: {{}}' to enable comparison."
            )
        return self

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

        config = cls(**data)
        config.source_path = path.resolve()
        return config

    def to_yaml(self, path: Path | str) -> None:
        """Save comparison config to YAML file.

        Parameters
        ----------
        path : Path or str
            Output path for comparison.yaml
        """
        path = Path(path)

        # Convert to dict, handling Path objects and nested containers
        data = self.model_dump(mode="json")
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

        # Check analysis/comparison coverage
        enabled_analyses = self.analysis_settings.get_enabled_analyses()
        enabled_comparisons = self.comparison_settings.get_enabled_comparisons()
        missing = set(enabled_analyses) - set(enabled_comparisons)
        if missing:
            errors.append(
                f"Missing comparison_settings for: {sorted(missing)}. "
                f"Add comparison_settings entries for these analyses."
            )

        return errors

    def generate_analysis_yaml(self, condition: ConditionConfig) -> str:
        """Generate analysis.yaml content for a specific condition.

        Parameters
        ----------
        condition : ConditionConfig
            The condition to generate analysis.yaml for.

        Returns
        -------
        str
            YAML content for the analysis.yaml file.
        """
        data = self.analysis_settings.to_analysis_yaml_dict(
            replicates=condition.replicates,
            eq_time=self.defaults.equilibration_time,
        )
        return yaml.dump(data, default_flow_style=False, sort_keys=False)

    def generate_analysis_yaml_for_all(self) -> dict[str, str]:
        """Generate analysis.yaml content for all conditions.

        Returns
        -------
        dict[str, str]
            Dictionary mapping condition labels to analysis.yaml content.
        """
        return {cond.label: self.generate_analysis_yaml(cond) for cond in self.conditions}


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
    return f"""\
# ============================================================================
# PolyzyMD Comparison Configuration
# ============================================================================
# Compare analyses across simulation conditions (e.g., polymer, temperature).
# Docs: https://polyzymd.readthedocs.io/en/latest/
# ============================================================================

name: "{name}"
description: "Comparison of simulation conditions"

# Control condition for relative comparisons.
# Must match one of the 'label' values in 'conditions' below, or null if none.
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
defaults:
  equilibration_time: "{eq_time}"

# ============================================================================
# Analysis Settings (WHAT to analyze - applied to all conditions)
# ============================================================================
# Define which analyses to run. Presence of a section enables that analysis.
# Running `polyzymd compare analyze` will run these for each condition.

analysis_settings:
  # RMSF Analysis
  rmsf:
    selection: "protein and name CA"
    reference_mode: "centroid"  # centroid, average, or frame
    # reference_frame: 500      # Required if reference_mode is "frame"

  # Catalytic Triad / Active Site Distances
  #
  # IMPORTANT: Always use "protein and resid X" for protein residues!
  # Residue numbers restart per chain. Without "protein and", your selection
  # may match atoms from polymer or water chains, causing incorrect distances.
  #
  # catalytic_triad:
  #   name: "enzyme_catalytic_triad"
  #   threshold: 3.5  # Angstroms (H-bond cutoff)
  #   pairs:
  #     - label: "Asp-His"
  #       selection_a: "midpoint(protein and resid 133 and name OD1 OD2)"
  #       selection_b: "protein and resid 156 and name ND1"
  #     - label: "His-Ser"
  #       selection_a: "protein and resid 156 and name NE2"
  #       selection_b: "protein and resid 77 and name OG"

  # Distance Analysis (general inter-atomic distances)
  #
  # IMPORTANT: Always use "protein and resid X" for protein residues!
  # See warning above in catalytic_triad section.
  #
  # distances:
  #   threshold: 3.5
  #   pairs:
  #     - label: "Ser77-Substrate"
  #       selection_a: "protein and resid 77 and name OG"
  #       selection_b: "resname RBY and name C1"

  # Polymer-Protein Contact Analysis
  # contacts:
  #   polymer_selection: "chainID C"
  #   protein_selection: "protein"
  #   cutoff: 4.5
  #   grouping: "aa_class"  # aa_class, secondary_structure, or none
  #   compute_residence_times: true
  #
  #   # Binding Preference Analysis (enrichment by residue group)
  #   # --------------------------------------------------------
  #   # Computes which residue types (aromatic, polar, etc.) are preferentially
  #   # contacted by the polymer, normalized by surface exposure.
  #   #
  #   # IMPORTANT: Place your enzyme PDB in the structures/ directory!
  #   # The path is relative to this comparison.yaml file.
  #   #
  #   # compute_binding_preference: true
  #   # surface_exposure_threshold: 0.2  # 20% relative SASA = surface exposed
  #   # enzyme_pdb_for_sasa: "structures/enzyme.pdb"
  #   # include_default_aa_groups: true  # aromatic, polar, nonpolar, charged
  #   #
  #   # Custom protein groups (residue IDs from your enzyme):
  #   # protein_groups:
  #   #   catalytic_triad: [77, 133, 156]
  #   #   lid_helix_5: [141, 142, 143, 144, 145]
  #   #   lid_helix_10: [281, 282, 283, 284, 285]
  #   #
  #   # User-defined partitions for system coverage plots
  #   # -------------------------------------------------
  #   # Groups within a partition MUST be mutually exclusive (no overlapping
  #   # residues). A "rest_of_protein" element is automatically added if the
  #   # groups don't cover all surface-exposed residues. One plot is generated
  #   # per partition.
  #   #
  #   # protein_partitions:
  #   #   lid_helices:          # partition name (becomes plot title)
  #   #     - lid_helix_5       # must be defined in protein_groups above
  #   #     - lid_helix_10

  # Exposure Dynamics Analysis (chaperone-like polymer activity)
  # Requires contacts analysis to be run first for each condition.
  # Run: polyzymd compare exposure
  #
  # exposure:
  #   exposure_threshold: 0.20     # fraction SASA defining 'exposed' residue
  #   transient_lower: 0.20        # lower bound: residue must reach this to be transient
  #   transient_upper: 0.80        # upper bound: residue must also reach this threshold
  #   min_event_length: 1          # minimum consecutive frames to count as an event
  #   protein_chain: "A"           # chain ID for the protein
  #   protein_selection: "protein" # MDAnalysis selection for protein
  #   polymer_selection: "chainID C"  # MDAnalysis selection for polymer
  #   # polymer_resnames: [SBM, EGM]  # optional: residue names for enrichment analysis
  #   probe_radius_nm: 0.14        # SASA probe radius (nm)
  #   n_sphere_points: 960         # number of sphere points for SASA computation

  # Binding Free Energy Analysis (ΔΔG via Boltzmann inversion)
  # Requires contacts analysis with compute_binding_preference: true to be run first.
  # Converts binding preference probabilities into ΔΔG = -k_B·T·ln(contact_share / expected_share).
  # Run: polyzymd compare binding-free-energy
  #
  # binding_free_energy:
  #   units: kcal/mol              # energy units: kcal/mol or kJ/mol
  #   surface_exposure_threshold: 0.2  # minimum relative SASA to be considered surface-exposed
  #   # protein_partitions: null   # optional: restrict to user-defined AA partitions

# ============================================================================
# Comparison Settings (HOW to compare - statistical parameters)
# ============================================================================
# Each analysis in analysis_settings MUST have a corresponding entry here.
# Use empty {{}} for analyses with no comparison-specific parameters.

comparison_settings:
  rmsf: {{}}  # No comparison-specific parameters

  # catalytic_triad: {{}}

  # distances: {{}}

  # contacts:
  #   fdr_alpha: 0.05           # FDR for Benjamini-Hochberg correction
  #   min_effect_size: 0.5      # Cohen's d threshold (0.2=small, 0.5=medium, 0.8=large)
  #   top_residues: 10          # Number of top residues to show in console

  # exposure: {{}}              # No comparison-specific parameters for exposure

  # binding_free_energy:
  #   fdr_alpha: 0.05           # FDR for Benjamini-Hochberg correction
"""
