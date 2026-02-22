"""Result models for binding free energy comparison analysis.

Physics background
------------------
In the NPT ensemble (constant pressure, as used in all polyzymd simulations)
the correct thermodynamic potential is the **Gibbs free energy** G.

The quantity computed here is a **selectivity free energy difference** (ΔΔG)
that measures how much more (or less) favorable it is for a polymer to contact
a given group of protein residues compared to what would be expected if the
polymer contacted each exposed surface residue in proportion to that residue
group's share of the total solvent-exposed protein surface.

Concretely: if aromatic residues make up 10% of the solvent-exposed surface
but receive 20% of the polymer's contacts, the polymer preferentially contacts
aromatic residues. The reference (expected) distribution is simply proportional
to surface availability — not any property of the polymer itself.

    ΔΔG_j = -k_B·T · ln(contact_share_j / expected_share_j)

where:
    contact_share_j  = (contact frames involving residues in group j) /
                       (total contact frames across all protein residues)
                       — the observed fraction of polymer contacts directed
                       at group j
    expected_share_j = (number of solvent-exposed residues in group j) /
                       (total number of solvent-exposed protein residues)
                       — the fraction of the protein surface belonging to
                       group j; this is the reference assuming contacts are
                       distributed purely by surface area
    k_B              = Boltzmann constant (0.0019872041 kcal mol⁻¹ K⁻¹)
    T                = simulation temperature in Kelvin

Note: contact_share / expected_share = enrichment_ratio = enrichment + 1
(where enrichment is the existing dimensionless enrichment score from binding
preference analysis). So ΔΔG = -kT·ln(enrichment + 1), and the two
representations are mathematically equivalent; ΔΔG simply puts the enrichment
score on a physically meaningful energy scale.

Sign convention:
    ΔΔG < 0  →  preferential contact (observed > surface-availability reference)
    ΔΔG > 0  →  contact avoidance (observed < surface-availability reference)
    ΔΔG = 0  →  contacts match the surface-availability reference exactly

Uncertainty propagation
-----------------------
When multiple independent replicates are available, two uncertainty estimates
are reported:

1. **Between-replicate SEM on ΔΔG** (primary, used for pairwise statistics):
   ΔΔG is computed independently for each replicate, and the SEM is taken
   directly across those values. This is the most statistically sound approach
   for independent replicates and is the quantity used in t-tests.

2. **Delta-method propagation** (analytical approximation, stored for reference):
   For the mean contact_share and its SEM, uncertainty is propagated through
   the logarithm using first-order error propagation (Taylor 1997, ch. 3;
   Bevington & Robinson 2003, ch. 3):

       σ(ΔΔG) ≈ k_B·T · √[(σ_cs / cs)² + (σ_es / es)²]

   where σ_cs = SEM of contact_share across replicates, and σ_es ≈ 0 because
   expected_share is computed from a single static PDB structure (no replicate
   variance). This simplifies to σ(ΔΔG) ≈ k_B·T · (σ_cs / cs).

   References:
   - Taylor, J. R. (1997). *An Introduction to Error Analysis*, 2nd ed.
     University Science Books. (Ch. 3: Error propagation for functions of
     one or more variables)
   - Bevington, P. R. & Robinson, D. K. (2003). *Data Reduction and Error
     Analysis for the Physical Sciences*, 3rd ed. McGraw-Hill. (Ch. 3)
   - Wikipedia: Delta method,
     https://en.wikipedia.org/wiki/Delta_method

Temperature handling:
    ΔΔG computed at temperature T is NOT directly comparable to ΔΔG at
    temperature T'. The probability ratios and kT scale factors both differ.
    Pairwise statistical comparisons are only computed between conditions
    sharing the same simulation temperature. Cross-temperature ΔΔG values
    are displayed with their respective temperatures labeled clearly.
"""

from __future__ import annotations

from datetime import datetime
from pathlib import Path
from typing import Optional

from pydantic import BaseModel, Field

from polyzymd import __version__


class FreeEnergyEntry(BaseModel):
    """Free energy analysis for one (polymer_type, protein_group) pair in one condition.

    Stores both the ΔΔG value and the raw probability quantities used to compute
    it, enabling reproducibility and downstream verification.

    Attributes
    ----------
    polymer_type : str
        Polymer residue type (e.g., "SBM", "EGM").
    protein_group : str
        Protein amino acid group label (e.g., "aromatic", "charged_positive").
    partition_name : str
        Name of the partition this group belongs to (e.g., "aa_class").
    contact_share : float
        Observed fraction of polymer contacts directed at this group.
        This is P_obs in ΔΔG = -kT·ln(P_obs / P_ref).
    expected_share : float
        Surface-availability-weighted reference fraction.
        This is P_ref in ΔΔG = -kT·ln(P_obs / P_ref).
    enrichment_ratio : float
        contact_share / expected_share (= enrichment + 1).
        Stored for traceability; ΔΔG = -kT·ln(enrichment_ratio).
    delta_G : float | None
        ΔΔG in the configured units. None when contact_share = 0 or
        expected_share = 0 (log undefined or reference missing).
    delta_G_uncertainty : float | None
        σ(ΔΔG) from delta-method error propagation. None if delta_G is None
        or if SEM data is unavailable (single replicate).
    delta_G_per_replicate : list[float]
        Per-replicate ΔΔG values used for cross-condition statistics.
    units : str
        Energy units ("kcal/mol" or "kJ/mol").
    temperature_K : float
        Simulation temperature in Kelvin (used as kT denominator).
    n_replicates : int
        Number of replicates with valid data for this entry.
    n_exposed_in_group : int
        Number of surface-exposed residues in this group (used for expected_share).
    """

    polymer_type: str
    protein_group: str
    partition_name: str = "aa_class"
    contact_share: float
    expected_share: float
    enrichment_ratio: float = Field(
        description="contact_share / expected_share; ΔΔG = -kT·ln(this value)"
    )
    delta_G: Optional[float] = Field(
        default=None,
        description="ΔΔG = -k_B·T·ln(contact_share / expected_share) in units",
    )
    delta_G_uncertainty: Optional[float] = Field(
        default=None,
        description="σ(ΔΔG) from delta-method: k_B·T·√[(σ_cs/cs)²+(σ_es/es)²]",
    )
    delta_G_per_replicate: list[float] = Field(
        default_factory=list,
        description="Per-replicate ΔΔG values for statistical testing",
    )
    units: str = "kcal/mol"
    temperature_K: float
    n_replicates: int = 0
    n_exposed_in_group: int = 0


class FreeEnergyConditionSummary(BaseModel):
    """Free energy summary for one simulation condition.

    Aggregates FreeEnergyEntry objects across all (polymer_type, protein_group)
    pairs for a single condition, together with condition metadata.

    Attributes
    ----------
    label : str
        Display name for this condition.
    config_path : str
        Path to the SimulationConfig YAML used.
    temperature_K : float
        Simulation temperature in Kelvin.
    n_replicates : int
        Number of replicates in this condition.
    units : str
        Energy units ("kcal/mol" or "kJ/mol").
    entries : list[FreeEnergyEntry]
        All (polymer_type, protein_group) ΔΔG entries.
    polymer_types : list[str]
        Polymer residue types present.
    protein_groups : list[str]
        Protein group labels analyzed.
    """

    label: str
    config_path: str
    temperature_K: float
    n_replicates: int
    units: str = "kcal/mol"
    entries: list[FreeEnergyEntry] = Field(default_factory=list)
    polymer_types: list[str] = Field(default_factory=list)
    protein_groups: list[str] = Field(default_factory=list)

    @property
    def primary_metric_value(self) -> float:
        """Mean ΔΔG across all valid entries (for BaseConditionSummary compatibility)."""
        vals = [e.delta_G for e in self.entries if e.delta_G is not None]
        return float(sum(vals) / len(vals)) if vals else 0.0

    @property
    def primary_metric_sem(self) -> float:
        """Mean σ(ΔΔG) across all valid entries."""
        vals = [e.delta_G_uncertainty for e in self.entries if e.delta_G_uncertainty is not None]
        return float(sum(vals) / len(vals)) if vals else 0.0

    def get_entry(self, polymer_type: str, protein_group: str) -> Optional[FreeEnergyEntry]:
        """Get the FreeEnergyEntry for a (polymer_type, protein_group) pair.

        Parameters
        ----------
        polymer_type : str
            Polymer type.
        protein_group : str
            AA group label.

        Returns
        -------
        FreeEnergyEntry or None
        """
        for e in self.entries:
            if e.polymer_type == polymer_type and e.protein_group == protein_group:
                return e
        return None


class FreeEnergyPairwiseEntry(BaseModel):
    """Pairwise ΔΔG comparison between two conditions for one (polymer, group) pair.

    Statistics are only computed when both conditions share the same simulation
    temperature. If temperatures differ, all stat fields are None and the
    ``cross_temperature`` flag is set to True.

    Attributes
    ----------
    polymer_type : str
        Polymer residue type.
    protein_group : str
        Protein group label.
    condition_a : str
        Label of the first condition.
    condition_b : str
        Label of the second condition.
    temperature_a_K : float
        Temperature of condition A in Kelvin.
    temperature_b_K : float
        Temperature of condition B in Kelvin.
    cross_temperature : bool
        True when temperatures differ — statistics are suppressed.
    delta_G_a : float | None
        ΔΔG for condition A.
    delta_G_b : float | None
        ΔΔG for condition B.
    delta_delta_G : float | None
        Difference: ΔΔG_B − ΔΔG_A. Positive → B has less favorable selectivity.
    t_statistic : float | None
        T-test statistic (None for cross-temperature pairs).
    p_value : float | None
        Two-tailed p-value (None for cross-temperature pairs).
    """

    polymer_type: str
    protein_group: str
    condition_a: str
    condition_b: str
    temperature_a_K: float
    temperature_b_K: float
    cross_temperature: bool = False
    delta_G_a: Optional[float] = None
    delta_G_b: Optional[float] = None
    delta_delta_G: Optional[float] = None
    t_statistic: Optional[float] = None
    p_value: Optional[float] = None


class BindingFreeEnergyResult(BaseModel):
    """Complete binding free energy comparison result.

    This is the main output from BindingFreeEnergyComparator.compare().

    Physics summary
    ---------------
    Formula: ΔΔG = -k_B·T · ln(contact_share / expected_share)

    Uncertainty: σ(ΔΔG) = k_B·T · √[(σ_cs/cs)² + (σ_es/es)²]

    Temperature note: pairwise statistics are suppressed between conditions
    at different temperatures. The ``mixed_temperatures`` flag indicates this
    occurred. Each condition's temperature is stored in its summary.

    Attributes
    ----------
    name : str
        Name of the comparison project.
    units : str
        Energy units ("kcal/mol" or "kJ/mol").
    formula : str
        Human-readable formula string (for documentation/output).
    mixed_temperatures : bool
        True if conditions span more than one simulation temperature.
    temperature_groups : dict[float, list[str]]
        Mapping of temperature (K) to condition labels at that temperature.
    conditions : list[FreeEnergyConditionSummary]
        Summary for each condition.
    pairwise_comparisons : list[FreeEnergyPairwiseEntry]
        All pairwise comparisons (cross-T pairs have stats suppressed).
    polymer_types : list[str]
        All polymer types found.
    protein_groups : list[str]
        All protein groups analyzed.
    surface_exposure_threshold : float | None
        SASA threshold used (from binding preference settings).
    equilibration_time : str
        Equilibration time used.
    created_at : datetime
        When the analysis was run.
    polyzymd_version : str
        Version of polyzymd used.
    """

    name: str
    units: str = "kcal/mol"
    formula: str = "ΔΔG = -k_B·T · ln(contact_share / expected_share)"
    mixed_temperatures: bool = False
    temperature_groups: dict[str, list[str]] = Field(
        default_factory=dict,
        description="temperature_K (as str key) → list of condition labels",
    )
    conditions: list[FreeEnergyConditionSummary] = Field(default_factory=list)
    pairwise_comparisons: list[FreeEnergyPairwiseEntry] = Field(default_factory=list)
    polymer_types: list[str] = Field(default_factory=list)
    protein_groups: list[str] = Field(default_factory=list)
    surface_exposure_threshold: Optional[float] = None
    equilibration_time: str = ""
    created_at: datetime = Field(default_factory=datetime.now)
    polyzymd_version: str = Field(default_factory=lambda: __version__)

    def save(self, path: Path | str) -> Path:
        """Save result to JSON file.

        Parameters
        ----------
        path : Path or str
            Output path.

        Returns
        -------
        Path
            Path to the saved file.
        """
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(self.model_dump_json(indent=2))
        return path

    @classmethod
    def load(cls, path: Path | str) -> "BindingFreeEnergyResult":
        """Load result from JSON file.

        Parameters
        ----------
        path : Path or str
            Path to JSON file.

        Returns
        -------
        BindingFreeEnergyResult
            Loaded result.
        """
        path = Path(path)
        return cls.model_validate_json(path.read_text())

    def get_condition(self, label: str) -> FreeEnergyConditionSummary:
        """Get a condition summary by label.

        Parameters
        ----------
        label : str
            Condition label.

        Returns
        -------
        FreeEnergyConditionSummary

        Raises
        ------
        KeyError
            If not found.
        """
        for cond in self.conditions:
            if cond.label == label:
                return cond
        raise KeyError(f"Condition '{label}' not found")
