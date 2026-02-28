"""Result models for polymer affinity score comparison analysis.

The polymer affinity score is a comparative metric that quantifies the total
strength of polymer-protein interactions by summing per-contact free energy
contributions weighted by the number of simultaneous contacts.

Physics
-------
For each (polymer_type, protein_group) pair, the affinity score is:

    S_{p,g} = N_{p,g} × ΔΔG_{p,g}

where:
    N_{p,g}   = mean number of simultaneous contacts per frame
              = mean_contact_fraction × n_exposed_in_group
    ΔΔG_{p,g} = -ln(contact_share / expected_share)  [in units of k_bT]

The total affinity score for a polymer type is:

    S_p = Σ_g S_{p,g}

The total affinity score for a condition is:

    S = Σ_p S_p

Independence assumption
-----------------------
This formulation assumes contacts are thermodynamically independent — each
contact contributes the same free energy regardless of what other contacts
exist simultaneously. This is the standard polyvalent binding approximation
(Mammen et al., Angew. Chem. Int. Ed. 1998, 37, 2754).

The absolute values are NOT rigorous thermodynamic binding free energies.
However, the *relative differences* between polymer compositions are meaningful
as a comparative scoring function, analogous to scoring functions in molecular
docking or MM/PBSA decomposition.

Sign convention
---------------
    S < 0  →  net favorable polymer-protein interaction
    S > 0  →  net unfavorable (avoidance dominates)
    S = 0  →  contacts match the surface-availability reference

Interpretation
--------------
More negative total score → stronger net polymer-protein affinity. When
combined with structural stability metrics (RMSF, triad contacts), the
affinity score helps rank polymer compositions by total interaction strength.

Uncertainty propagation
-----------------------
Per-replicate scores are computed independently:

    S_rep = N_rep × ΔΔG_rep

where N_rep = contact_fraction_rep × n_exposed_in_group, and
ΔΔG_rep = -ln(enrichment_rep + 1). The mean and SEM are taken across
replicates. This approach naturally captures the covariance between N and ΔΔG.

When per-replicate data is unavailable, analytical error propagation is used:

    σ(S) = √[(N·σ_ΔΔG)² + (ΔΔG·σ_N)²]
"""

from __future__ import annotations

from datetime import datetime
from pathlib import Path
from typing import Optional

from pydantic import BaseModel, Field

from polyzymd import __version__


class AffinityScoreEntry(BaseModel):
    """Affinity score for one (polymer_type, protein_group) pair in one condition.

    Stores both the composite score and its constituent quantities for
    reproducibility and downstream verification.

    Attributes
    ----------
    polymer_type : str
        Polymer residue type (e.g., "SBM", "EGM").
    protein_group : str
        Protein amino acid group label (e.g., "aromatic", "charged_positive").
    partition_name : str
        Name of the partition this group belongs to (e.g., "aa_class").
    n_contacts : float
        Mean number of simultaneous contacts per frame.
        Computed as mean_contact_fraction * n_exposed_in_group.
    delta_G_per_contact : float | None
        Per-contact selectivity free energy in kT.
        Computed as -ln(contact_share / expected_share).
    affinity_score : float | None
        Composite score: n_contacts * delta_G_per_contact (kT).
        More negative = stronger favorable interaction.
    affinity_score_uncertainty : float | None
        Uncertainty on affinity_score. From replicate SEM when available,
        otherwise from analytical error propagation.
    affinity_score_per_replicate : list[float]
        Per-replicate affinity scores for statistical testing.
    mean_contact_fraction : float
        Mean per-residue contact fraction in this group (from binding preference).
    n_exposed_in_group : int
        Number of surface-exposed residues in this group.
    contact_share : float
        Observed fraction of polymer contacts directed at this group.
    expected_share : float
        Surface-availability reference fraction.
    temperature_K : float
        Simulation temperature in Kelvin.
    n_replicates : int
        Number of replicates with valid data.
    """

    polymer_type: str
    protein_group: str
    partition_name: str = "aa_class"
    n_contacts: float = Field(description="Mean simultaneous contacts per frame = mcf * n_exposed")
    delta_G_per_contact: Optional[float] = Field(
        default=None,
        description="Per-contact ΔΔG = -ln(contact_share / expected_share) [kT]",
    )
    affinity_score: Optional[float] = Field(
        default=None,
        description="N_contacts × ΔΔG_per_contact [kT]; negative = favorable",
    )
    affinity_score_uncertainty: Optional[float] = Field(
        default=None,
        description="σ(affinity_score) from replicate SEM or error propagation",
    )
    affinity_score_per_replicate: list[float] = Field(
        default_factory=list,
        description="Per-replicate affinity scores for statistical testing",
    )
    mean_contact_fraction: float = Field(
        default=0.0,
        description="Mean per-residue contact fraction in group",
    )
    n_exposed_in_group: int = Field(
        default=0,
        description="Surface-exposed residues in group",
    )
    contact_share: float = Field(
        default=0.0,
        description="Observed fraction of polymer contacts to this group",
    )
    expected_share: float = Field(
        default=0.0,
        description="Expected share based on surface availability",
    )
    temperature_K: float = 0.0
    n_replicates: int = 0


class PolymerTypeScore(BaseModel):
    """Aggregated affinity score for one polymer type across all protein groups.

    The score is the sum of per-group affinity scores:
        S_p = Σ_g (N_g × ΔΔG_g)

    Attributes
    ----------
    polymer_type : str
        Polymer residue type (e.g., "SBM", "EGM").
    total_score : float
        Sum of affinity scores across all protein groups (kT).
    total_score_uncertainty : float | None
        Uncertainty on total_score.
    total_score_per_replicate : list[float]
        Per-replicate total scores for statistical testing.
    total_n_contacts : float
        Total mean simultaneous contacts per frame across all groups.
    group_contributions : list[AffinityScoreEntry]
        Breakdown by protein group (for detail reporting).
    """

    polymer_type: str
    total_score: float = Field(
        description="Σ_g (N_g × ΔΔG_g) summed over protein groups [kT]",
    )
    total_score_uncertainty: Optional[float] = Field(
        default=None,
        description="σ(total_score) from replicate SEM or quadrature sum",
    )
    total_score_per_replicate: list[float] = Field(
        default_factory=list,
        description="Per-replicate total scores for statistical testing",
    )
    total_n_contacts: float = Field(
        default=0.0,
        description="Total mean simultaneous contacts per frame",
    )
    group_contributions: list[AffinityScoreEntry] = Field(
        default_factory=list,
        description="Per-group breakdown of affinity scores",
    )


class AffinityScoreConditionSummary(BaseModel):
    """Affinity score summary for one simulation condition.

    Aggregates scores at three levels: per (polymer_type, protein_group),
    per polymer_type, and total condition score.

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
    total_score : float
        Grand total affinity score across all polymer types and groups (kT).
    total_score_uncertainty : float | None
        Uncertainty on total_score.
    total_score_per_replicate : list[float]
        Per-replicate grand total scores for pairwise statistics.
    total_n_contacts : float
        Total mean simultaneous contacts per frame (all types, all groups).
    polymer_type_scores : list[PolymerTypeScore]
        Per-polymer-type score breakdown.
    entries : list[AffinityScoreEntry]
        All (polymer_type, protein_group) entries.
    polymer_types : list[str]
        Polymer types present.
    protein_groups : list[str]
        Protein groups analyzed.
    """

    label: str
    config_path: str
    temperature_K: float
    n_replicates: int = 0
    total_score: float = Field(
        default=0.0,
        description="Grand total affinity score [kT]",
    )
    total_score_uncertainty: Optional[float] = Field(
        default=None,
        description="σ(total_score)",
    )
    total_score_per_replicate: list[float] = Field(
        default_factory=list,
        description="Per-replicate grand total scores",
    )
    total_n_contacts: float = Field(
        default=0.0,
        description="Total mean simultaneous contacts per frame",
    )
    polymer_type_scores: list[PolymerTypeScore] = Field(
        default_factory=list,
    )
    entries: list[AffinityScoreEntry] = Field(default_factory=list)
    polymer_types: list[str] = Field(default_factory=list)
    protein_groups: list[str] = Field(default_factory=list)

    @property
    def primary_metric_value(self) -> float:
        """Total affinity score (for ranking compatibility)."""
        return self.total_score

    @property
    def primary_metric_sem(self) -> float:
        """Uncertainty on total affinity score."""
        return self.total_score_uncertainty if self.total_score_uncertainty is not None else 0.0


class AffinityScorePairwiseEntry(BaseModel):
    """Pairwise affinity score comparison between two conditions.

    Compares total affinity scores. Statistics are suppressed for
    cross-temperature pairs.

    Attributes
    ----------
    condition_a : str
        Label of the first condition (typically control or reference).
    condition_b : str
        Label of the second condition.
    temperature_a_K : float
        Temperature of condition A in Kelvin.
    temperature_b_K : float
        Temperature of condition B in Kelvin.
    cross_temperature : bool
        True when temperatures differ (statistics suppressed).
    score_a : float
        Total affinity score for condition A (kT).
    score_b : float
        Total affinity score for condition B (kT).
    delta_score : float | None
        Difference: score_B - score_A (kT).
        Negative = B has stronger affinity than A.
    t_statistic : float | None
        T-test statistic (None for cross-temperature pairs).
    p_value : float | None
        Two-tailed p-value (None for cross-temperature pairs).
    """

    condition_a: str
    condition_b: str
    temperature_a_K: float
    temperature_b_K: float
    cross_temperature: bool = False
    score_a: float = 0.0
    score_b: float = 0.0
    delta_score: Optional[float] = None
    t_statistic: Optional[float] = None
    p_value: Optional[float] = None


class PolymerAffinityScoreResult(BaseModel):
    """Complete polymer affinity score comparison result.

    This is the main output from PolymerAffinityScoreComparator.compare().

    The polymer affinity score quantifies total polymer-protein interaction
    strength as a comparative metric. It is computed by summing per-contact
    selectivity free energies weighted by the number of simultaneous contacts:

        S = Σ_{p,g} N_{p,g} × ΔΔG_{p,g}

    where the sum runs over all (polymer_type, protein_group) pairs.

    .. important::

       This quantity assumes contact independence and should be interpreted
       as a relative affinity score, not a rigorous thermodynamic binding
       free energy. See the module docstring for details.

    Attributes
    ----------
    name : str
        Name of the comparison project.
    methodology : str
        Human-readable description of the scoring methodology.
    mixed_temperatures : bool
        True if conditions span more than one simulation temperature.
    temperature_groups : dict[str, list[str]]
        Mapping of temperature (K, as str) to condition labels.
    conditions : list[AffinityScoreConditionSummary]
        Summary for each condition.
    pairwise_comparisons : list[AffinityScorePairwiseEntry]
        All pairwise comparisons.
    polymer_types : list[str]
        All polymer types found.
    protein_groups : list[str]
        All protein groups analyzed.
    surface_exposure_threshold : float | None
        SASA threshold used (from settings).
    equilibration_time : str
        Equilibration time used.
    created_at : datetime
        When the analysis was run.
    polyzymd_version : str
        Version of polyzymd used.
    """

    name: str
    methodology: str = (
        "Polymer Affinity Score: S = Σ (N_contacts × ΔΔG_per_contact) [kT]. "
        "N_contacts = mean_contact_fraction × n_exposed_in_group. "
        "ΔΔG_per_contact = -ln(contact_share / expected_share). "
        "More negative = stronger net polymer-protein affinity. "
        "Assumes contact independence; interpret as comparative scoring metric."
    )
    mixed_temperatures: bool = False
    temperature_groups: dict[str, list[str]] = Field(
        default_factory=dict,
        description="temperature_K (as str key) → list of condition labels",
    )
    conditions: list[AffinityScoreConditionSummary] = Field(default_factory=list)
    pairwise_comparisons: list[AffinityScorePairwiseEntry] = Field(default_factory=list)
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
    def load(cls, path: Path | str) -> "PolymerAffinityScoreResult":
        """Load result from JSON file.

        Parameters
        ----------
        path : Path or str
            Path to JSON file.

        Returns
        -------
        PolymerAffinityScoreResult
            Loaded result.
        """
        path = Path(path)
        return cls.model_validate_json(path.read_text())

    def get_condition(self, label: str) -> AffinityScoreConditionSummary | None:
        """Look up a condition summary by label.

        Parameters
        ----------
        label : str
            Condition display name.

        Returns
        -------
        AffinityScoreConditionSummary or None
        """
        for c in self.conditions:
            if c.label == label:
                return c
        return None

    def get_ranking(self) -> list[AffinityScoreConditionSummary]:
        """Return conditions ranked by total affinity score (most negative first).

        Returns
        -------
        list[AffinityScoreConditionSummary]
            Conditions sorted by total_score ascending.
        """
        return sorted(
            [c for c in self.conditions if c.entries],
            key=lambda c: c.total_score,
        )
