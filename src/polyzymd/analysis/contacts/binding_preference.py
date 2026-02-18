"""Binding preference metrics for polymer-protein contacts.

This module computes enrichment ratios that answer the scientific question:
"Does polymer type X preferentially bind amino acid class Y?"

Enrichment Calculation
----------------------
For each (polymer_type, protein_group) pair:

    contact_share = Σ(contact_frames for exposed residues in group) /
                    Σ(contact_frames for all exposed residues)

    expected_share = n_exposed_residues_in_group / n_total_exposed_residues

    enrichment = contact_share / expected_share

Interpretation:
- Enrichment > 1.0 = preferential binding (polymer contacts this group
  more than expected by random chance)
- Enrichment = 1.0 = neutral (contact frequency matches surface availability)
- Enrichment < 1.0 = avoidance (polymer contacts this group less than expected)

The contact-frame weighting ensures that residues contacted for longer durations
contribute proportionally more to the enrichment calculation. A residue contacted
for 60% of the simulation contributes 60x more than one contacted for 1 frame.

Examples
--------
>>> from polyzymd.analysis.contacts.binding_preference import compute_binding_preference
>>> from polyzymd.analysis.contacts.surface_exposure import SurfaceExposureFilter
>>>
>>> # Load contact results and compute surface exposure
>>> contact_result = ContactResult.load("contacts.json")
>>> exposure_filter = SurfaceExposureFilter(threshold=0.2)
>>> surface_exposure = exposure_filter.calculate("enzyme.pdb")
>>>
>>> # Define protein groups (resolved to residue IDs)
>>> protein_groups = {
...     "aromatic": {12, 45, 67, 89},  # resids of aromatic residues
...     "charged": {23, 34, 56, 78},
... }
>>>
>>> # Compute binding preference
>>> result = compute_binding_preference(
...     contact_result,
...     surface_exposure,
...     protein_groups,
... )
>>>
>>> # Check enrichment
>>> print(result.get_enrichment("SBM", "aromatic"))
1.45  # SBMA prefers aromatic residues
>>> print(result.enrichment_matrix())
{'SBM': {'aromatic': 1.45, 'charged': 0.82}, 'EGM': {...}}
"""

from __future__ import annotations

import json
import logging
from pathlib import Path
from typing import TYPE_CHECKING, Any

import numpy as np
from pydantic import BaseModel, Field

from polyzymd.analysis.common.aa_classification import DEFAULT_AA_CLASS_SELECTIONS

if TYPE_CHECKING:
    from MDAnalysis.core.universe import Universe

    from polyzymd.analysis.config import ContactsConfig
    from polyzymd.analysis.contacts.results import ContactResult
    from polyzymd.analysis.contacts.surface_exposure import SurfaceExposureResult

logger = logging.getLogger(__name__)


class BindingPreferenceEntry(BaseModel):
    """Single entry in the binding preference matrix.

    Represents the binding preference metrics for one
    (polymer_type, protein_group) combination.

    Attributes
    ----------
    polymer_type : str
        Polymer residue type (e.g., "SBM", "EGM")
    protein_group : str
        Protein group label (e.g., "aromatic", "charged_positive")
    total_contact_frames : int
        Sum of contact frames across all exposed residues in this group.
        This is the numerator for contact_share calculation.
    mean_contact_fraction : float
        Average per-residue contact fraction within this group.
        (total_contact_frames / n_frames) / n_exposed_in_group
    n_residues_in_group : int
        Total residues in this protein group (exposed + buried)
    n_exposed_in_group : int
        Surface-exposed residues in this group (used for enrichment)
    n_residues_contacted : int
        Number of exposed residues that had at least one contact
    contact_share : float
        Fraction of this polymer's total contacts that went to this group.
        = total_contact_frames / total_polymer_contact_frames
    expected_share : float
        Expected share based on surface availability.
        = n_exposed_in_group / n_total_exposed
    enrichment_ratio : float | None
        contact_share / expected_share.
        None if denominator is zero.
    """

    polymer_type: str
    protein_group: str
    total_contact_frames: int = Field(
        description="Sum of contact frames for all exposed residues in group"
    )
    mean_contact_fraction: float = Field(
        description="Average per-residue contact fraction in group"
    )
    n_residues_in_group: int = Field(description="Total residues in this protein group")
    n_exposed_in_group: int = Field(description="Surface-exposed residues in group")
    n_residues_contacted: int = Field(
        default=0, description="Exposed residues with at least one contact"
    )
    contact_share: float = Field(
        default=0.0, description="Fraction of polymer's contacts to this group"
    )
    expected_share: float = Field(default=0.0, description="Expected share by surface availability")
    enrichment_ratio: float | None = Field(
        default=None, description="Enrichment vs expected (contact_share/expected_share)"
    )


class BindingPreferenceResult(BaseModel):
    """Complete binding preference analysis result.

    Provides enrichment-normalized metrics for polymer-protein
    binding preferences, answering questions like:

    - "Does SBMA preferentially bind aromatic residues?"
    - "How does EGMA's preference for charged residues compare to SBMA?"
    - "Which amino acid class does this polymer type prefer?"

    Attributes
    ----------
    entries : list[BindingPreferenceEntry]
        All (polymer_type × protein_group) combinations
    n_frames : int
        Total frames analyzed
    total_exposed_residues : int
        Number of surface-exposed residues considered
    surface_exposure_threshold : float
        SASA threshold used for surface filtering
    protein_groups_used : dict[str, str]
        Mapping of group name to MDAnalysis selection string
    polymer_types_used : dict[str, str]
        Mapping of polymer type name to MDAnalysis selection string
    schema_version : int
        Version for forward compatibility
    """

    entries: list[BindingPreferenceEntry] = Field(default_factory=list)
    n_frames: int = 0
    total_exposed_residues: int = 0
    surface_exposure_threshold: float | None = None
    protein_groups_used: dict[str, str] = Field(default_factory=dict)
    polymer_types_used: dict[str, str] = Field(default_factory=dict)
    metadata: dict[str, Any] = Field(default_factory=dict)
    schema_version: int = 1

    def to_dataframe(self) -> "pd.DataFrame":
        """Convert to pandas DataFrame for analysis/plotting.

        Returns
        -------
        pd.DataFrame
            Columns: polymer_type, protein_group, total_contact_frames,
            mean_contact_fraction, n_residues_in_group, n_exposed_in_group,
            n_residues_contacted, contact_share, expected_share, enrichment_ratio
        """
        import pandas as pd

        return pd.DataFrame([e.model_dump() for e in self.entries])

    def enrichment_matrix(self) -> dict[str, dict[str, float]]:
        """Get enrichment as nested dict: {polymer_type: {protein_group: value}}.

        Useful for quick lookups and programmatic comparison.

        Returns
        -------
        dict[str, dict[str, float]]
            Nested mapping of enrichment values.
            Missing/invalid values are returned as 0.0.

        Examples
        --------
        >>> matrix = result.enrichment_matrix()
        >>> print(matrix["SBM"]["aromatic"])
        1.45
        """
        result: dict[str, dict[str, float]] = {}
        for entry in self.entries:
            if entry.polymer_type not in result:
                result[entry.polymer_type] = {}
            result[entry.polymer_type][entry.protein_group] = (
                entry.enrichment_ratio if entry.enrichment_ratio is not None else 0.0
            )
        return result

    def contact_fraction_matrix(self) -> dict[str, dict[str, float]]:
        """Get mean contact fractions as nested dict.

        Returns
        -------
        dict[str, dict[str, float]]
            Nested mapping: {polymer_type: {protein_group: mean_frac}}
        """
        result: dict[str, dict[str, float]] = {}
        for entry in self.entries:
            if entry.polymer_type not in result:
                result[entry.polymer_type] = {}
            result[entry.polymer_type][entry.protein_group] = entry.mean_contact_fraction
        return result

    def contact_share_matrix(self) -> dict[str, dict[str, float]]:
        """Get contact shares as nested dict.

        Returns
        -------
        dict[str, dict[str, float]]
            Nested mapping: {polymer_type: {protein_group: contact_share}}
        """
        result: dict[str, dict[str, float]] = {}
        for entry in self.entries:
            if entry.polymer_type not in result:
                result[entry.polymer_type] = {}
            result[entry.polymer_type][entry.protein_group] = entry.contact_share
        return result

    def get_enrichment(self, polymer_type: str, protein_group: str) -> float | None:
        """Get enrichment for a specific (polymer_type, protein_group) pair.

        Parameters
        ----------
        polymer_type : str
            Polymer type name
        protein_group : str
            Protein group name

        Returns
        -------
        float or None
            Enrichment ratio, or None if pair not found
        """
        for entry in self.entries:
            if entry.polymer_type == polymer_type and entry.protein_group == protein_group:
                return entry.enrichment_ratio
        return None

    def get_entry(self, polymer_type: str, protein_group: str) -> BindingPreferenceEntry | None:
        """Get the full entry for a (polymer_type, protein_group) pair.

        Parameters
        ----------
        polymer_type : str
            Polymer type name
        protein_group : str
            Protein group name

        Returns
        -------
        BindingPreferenceEntry or None
            Full entry, or None if not found
        """
        for entry in self.entries:
            if entry.polymer_type == polymer_type and entry.protein_group == protein_group:
                return entry
        return None

    def polymer_types(self) -> list[str]:
        """Get list of polymer types in this result."""
        return sorted(set(e.polymer_type for e in self.entries))

    def protein_groups(self) -> list[str]:
        """Get list of protein groups in this result."""
        return sorted(set(e.protein_group for e in self.entries))

    def save(self, path: str | Path) -> None:
        """Save to JSON file.

        Parameters
        ----------
        path : str or Path
            Output path for JSON file
        """
        Path(path).write_text(json.dumps(self.model_dump(), indent=2))
        logger.info(f"Saved binding preference result to {path}")

    @classmethod
    def load(cls, path: str | Path) -> "BindingPreferenceResult":
        """Load from JSON file.

        Parameters
        ----------
        path : str or Path
            Path to JSON file

        Returns
        -------
        BindingPreferenceResult
            Loaded result
        """
        data = json.loads(Path(path).read_text())
        return cls.model_validate(data)


def compute_binding_preference(
    contact_result: "ContactResult",
    surface_exposure: "SurfaceExposureResult",
    protein_groups: dict[str, set[int]],
    polymer_types: list[str] | None = None,
    protein_group_selections: dict[str, str] | None = None,
    polymer_type_selections: dict[str, str] | None = None,
) -> BindingPreferenceResult:
    """Compute binding preference from contact results.

    This function computes enrichment ratios for each (polymer_type, protein_group)
    combination, answering: "Does this polymer type preferentially bind this
    protein group compared to random chance?"

    The enrichment calculation accounts for:
    1. Surface exposure (only exposed residues are considered)
    2. Contact duration (contact frames are summed, not binary counts)

    Parameters
    ----------
    contact_result : ContactResult
        Raw contact analysis results from trajectory
    surface_exposure : SurfaceExposureResult
        Surface exposure data for filtering buried residues
    protein_groups : dict[str, set[int]]
        Mapping of group name to residue IDs in that group.
        Only surface-exposed residues in each group are counted.
        Example: {"aromatic": {12, 45, 67}, "charged": {23, 34}}
    polymer_types : list[str], optional
        If provided, only compute for these polymer residue types.
        If None, all polymer types found in contacts are used.
    protein_group_selections : dict[str, str], optional
        Original MDAnalysis selections (for metadata/reproducibility)
    polymer_type_selections : dict[str, str], optional
        Original MDAnalysis selections (for metadata/reproducibility)

    Returns
    -------
    BindingPreferenceResult
        Enrichment-normalized binding preference metrics

    Notes
    -----
    The enrichment ratio is computed as:

        enrichment = contact_share / expected_share

    where:
        contact_share = polymer_contacts_to_group / polymer_total_contacts
        expected_share = n_exposed_in_group / n_total_exposed

    This normalizes by surface availability, so enrichment > 1 indicates
    that the polymer contacts this group more than expected by chance.
    """
    exposed_resids = surface_exposure.exposed_resids
    n_frames = contact_result.n_frames
    total_exposed = len(exposed_resids)

    logger.info(
        f"Computing binding preference: {total_exposed} exposed residues, "
        f"{len(protein_groups)} protein groups"
    )

    # Filter protein_groups to only exposed residues
    exposed_groups: dict[str, set[int]] = {}
    for group_name, resids in protein_groups.items():
        exposed_groups[group_name] = resids & exposed_resids
        n_orig = len(resids)
        n_exposed = len(exposed_groups[group_name])
        logger.debug(f"  {group_name}: {n_exposed}/{n_orig} residues are exposed")

    # Collect contact frame counts per (polymer_type, protein_group)
    # Structure: {polymer_type: {protein_group: {"frames": int, "residues": set}}}
    contact_data: dict[str, dict[str, dict[str, Any]]] = {}

    # Also track total contacts per polymer type (for contact_share denominator)
    total_contacts_by_polymer: dict[str, int] = {}

    for rc in contact_result.residue_contacts:
        resid = rc.protein_resid
        if resid not in exposed_resids:
            continue  # Skip buried residues

        # Determine which protein groups this residue belongs to
        residue_groups = [gname for gname, gresids in exposed_groups.items() if resid in gresids]

        if not residue_groups:
            # Residue is exposed but not in any defined group
            continue

        # Get contacts by polymer type for this residue
        contacts_by_type = rc.contacts_by_polymer_type(n_frames)

        for poly_type, frac in contacts_by_type.items():
            if polymer_types and poly_type not in polymer_types:
                continue

            # Convert fraction to frame count
            contact_frames = int(round(frac * n_frames))

            if contact_frames == 0:
                continue

            # Initialize polymer type if needed
            if poly_type not in contact_data:
                contact_data[poly_type] = {}
                total_contacts_by_polymer[poly_type] = 0

            # Add to total contacts for this polymer
            total_contacts_by_polymer[poly_type] += contact_frames

            # Add to each group this residue belongs to
            for gname in residue_groups:
                if gname not in contact_data[poly_type]:
                    contact_data[poly_type][gname] = {
                        "total_frames": 0,
                        "residues_contacted": set(),
                    }
                contact_data[poly_type][gname]["total_frames"] += contact_frames
                contact_data[poly_type][gname]["residues_contacted"].add(resid)

    # Build result entries with enrichment calculations
    entries = []

    for poly_type in sorted(contact_data.keys()):
        total_poly_contacts = total_contacts_by_polymer.get(poly_type, 0)

        for group_name in sorted(protein_groups.keys()):
            n_total_in_group = len(protein_groups.get(group_name, set()))
            n_exposed_in_group = len(exposed_groups.get(group_name, set()))

            # Get contact data for this (polymer, group) pair
            gdata = contact_data[poly_type].get(group_name, {})
            contact_frames = gdata.get("total_frames", 0)
            residues_contacted = gdata.get("residues_contacted", set())
            n_residues_contacted = len(residues_contacted)

            # Calculate mean contact fraction
            if n_exposed_in_group > 0 and n_frames > 0:
                # Per-residue average: total_frames / (n_frames * n_exposed)
                mean_frac = contact_frames / (n_frames * n_exposed_in_group)
            else:
                mean_frac = 0.0

            # Calculate contact share and expected share
            if total_poly_contacts > 0:
                contact_share = contact_frames / total_poly_contacts
            else:
                contact_share = 0.0

            if total_exposed > 0:
                expected_share = n_exposed_in_group / total_exposed
            else:
                expected_share = 0.0

            # Calculate enrichment ratio
            if expected_share > 0 and contact_share > 0:
                enrichment = contact_share / expected_share
            elif expected_share > 0 and contact_share == 0:
                enrichment = 0.0  # No contacts = zero enrichment
            else:
                enrichment = None  # Can't compute (no exposed residues in group)

            entries.append(
                BindingPreferenceEntry(
                    polymer_type=poly_type,
                    protein_group=group_name,
                    total_contact_frames=contact_frames,
                    mean_contact_fraction=mean_frac,
                    n_residues_in_group=n_total_in_group,
                    n_exposed_in_group=n_exposed_in_group,
                    n_residues_contacted=n_residues_contacted,
                    contact_share=contact_share,
                    expected_share=expected_share,
                    enrichment_ratio=enrichment,
                )
            )

    result = BindingPreferenceResult(
        entries=entries,
        n_frames=n_frames,
        total_exposed_residues=total_exposed,
        surface_exposure_threshold=surface_exposure.threshold,
        protein_groups_used=protein_group_selections or {},
        polymer_types_used=polymer_type_selections or {},
    )

    # Log summary
    polymer_types_found = result.polymer_types()
    logger.info(
        f"Binding preference computed: {len(entries)} entries for "
        f"{len(polymer_types_found)} polymer types × {len(protein_groups)} groups"
    )

    return result


def aggregate_binding_preference(
    results: list[BindingPreferenceResult],
) -> "AggregatedBindingPreferenceResult":
    """Aggregate binding preference across replicates.

    Computes mean ± SEM for enrichment ratios across multiple replicates.

    Parameters
    ----------
    results : list[BindingPreferenceResult]
        Binding preference results from multiple replicates

    Returns
    -------
    AggregatedBindingPreferenceResult
        Aggregated results with mean and SEM
    """
    if not results:
        raise ValueError("No results to aggregate")

    # Collect all (polymer_type, protein_group) pairs
    all_pairs: set[tuple[str, str]] = set()
    for r in results:
        for e in r.entries:
            all_pairs.add((e.polymer_type, e.protein_group))

    entries = []
    for poly_type, prot_group in sorted(all_pairs):
        # Collect values from each replicate
        enrichments = []
        contact_fractions = []
        contact_shares = []
        contact_frames = []

        for r in results:
            entry = r.get_entry(poly_type, prot_group)
            if entry is not None:
                if entry.enrichment_ratio is not None:
                    enrichments.append(entry.enrichment_ratio)
                contact_fractions.append(entry.mean_contact_fraction)
                contact_shares.append(entry.contact_share)
                contact_frames.append(entry.total_contact_frames)

        # Compute statistics
        n_reps = len(enrichments)
        if n_reps > 0:
            mean_enrichment = float(np.mean(enrichments))
            sem_enrichment = (
                float(np.std(enrichments, ddof=1) / np.sqrt(n_reps)) if n_reps > 1 else 0.0
            )
        else:
            mean_enrichment = None
            sem_enrichment = None

        mean_contact_fraction = float(np.mean(contact_fractions)) if contact_fractions else 0.0
        sem_contact_fraction = (
            float(np.std(contact_fractions, ddof=1) / np.sqrt(len(contact_fractions)))
            if len(contact_fractions) > 1
            else 0.0
        )

        mean_contact_share = float(np.mean(contact_shares)) if contact_shares else 0.0

        # Get group metadata from first result
        first_entry = results[0].get_entry(poly_type, prot_group)
        n_exposed = first_entry.n_exposed_in_group if first_entry else 0
        n_total = first_entry.n_residues_in_group if first_entry else 0
        expected_share = first_entry.expected_share if first_entry else 0.0

        entries.append(
            AggregatedBindingPreferenceEntry(
                polymer_type=poly_type,
                protein_group=prot_group,
                mean_enrichment=mean_enrichment,
                sem_enrichment=sem_enrichment,
                mean_contact_fraction=mean_contact_fraction,
                sem_contact_fraction=sem_contact_fraction,
                mean_contact_share=mean_contact_share,
                expected_share=expected_share,
                n_exposed_in_group=n_exposed,
                n_residues_in_group=n_total,
                n_replicates=n_reps,
                per_replicate_enrichments=enrichments,
            )
        )

    return AggregatedBindingPreferenceResult(
        entries=entries,
        n_replicates=len(results),
        total_exposed_residues=results[0].total_exposed_residues if results else 0,
        surface_exposure_threshold=results[0].surface_exposure_threshold if results else None,
        protein_groups_used=results[0].protein_groups_used if results else {},
        polymer_types_used=results[0].polymer_types_used if results else {},
    )


class AggregatedBindingPreferenceEntry(BaseModel):
    """Aggregated binding preference for one (polymer_type, protein_group) pair.

    Contains mean ± SEM across replicates.
    """

    polymer_type: str
    protein_group: str
    mean_enrichment: float | None = Field(description="Mean enrichment ratio across replicates")
    sem_enrichment: float | None = Field(description="Standard error of enrichment")
    mean_contact_fraction: float = Field(description="Mean per-residue contact fraction")
    sem_contact_fraction: float = Field(
        default=0.0, description="Standard error of contact fraction"
    )
    mean_contact_share: float = Field(default=0.0, description="Mean contact share")
    expected_share: float = Field(default=0.0, description="Expected share by surface availability")
    n_exposed_in_group: int = Field(default=0, description="Surface-exposed residues in group")
    n_residues_in_group: int = Field(default=0, description="Total residues in group")
    n_replicates: int = Field(default=0, description="Number of replicates with valid data")
    per_replicate_enrichments: list[float] = Field(
        default_factory=list, description="Enrichment values from each replicate"
    )


class AggregatedBindingPreferenceResult(BaseModel):
    """Binding preference aggregated across replicates.

    Contains mean ± SEM for all metrics across multiple replicates.
    """

    entries: list[AggregatedBindingPreferenceEntry] = Field(default_factory=list)
    n_replicates: int = 0
    total_exposed_residues: int = 0
    surface_exposure_threshold: float | None = None
    protein_groups_used: dict[str, str] = Field(default_factory=dict)
    polymer_types_used: dict[str, str] = Field(default_factory=dict)
    schema_version: int = 1

    def to_dataframe(self) -> "pd.DataFrame":
        """Convert to pandas DataFrame."""
        import pandas as pd

        return pd.DataFrame([e.model_dump() for e in self.entries])

    def enrichment_matrix(self) -> dict[str, dict[str, float]]:
        """Get mean enrichment as nested dict."""
        result: dict[str, dict[str, float]] = {}
        for entry in self.entries:
            if entry.polymer_type not in result:
                result[entry.polymer_type] = {}
            result[entry.polymer_type][entry.protein_group] = (
                entry.mean_enrichment if entry.mean_enrichment is not None else 0.0
            )
        return result

    def get_entry(
        self, polymer_type: str, protein_group: str
    ) -> AggregatedBindingPreferenceEntry | None:
        """Get entry for a (polymer_type, protein_group) pair."""
        for entry in self.entries:
            if entry.polymer_type == polymer_type and entry.protein_group == protein_group:
                return entry
        return None

    def polymer_types(self) -> list[str]:
        """Get list of polymer types."""
        return sorted(set(e.polymer_type for e in self.entries))

    def protein_groups(self) -> list[str]:
        """Get list of protein groups."""
        return sorted(set(e.protein_group for e in self.entries))

    def save(self, path: str | Path) -> None:
        """Save to JSON file."""
        Path(path).write_text(json.dumps(self.model_dump(), indent=2))
        logger.info(f"Saved aggregated binding preference to {path}")

    @classmethod
    def load(cls, path: str | Path) -> "AggregatedBindingPreferenceResult":
        """Load from JSON file."""
        data = json.loads(Path(path).read_text())
        return cls.model_validate(data)


# =============================================================================
# Selection Resolution Helpers
# =============================================================================


def resolve_protein_group_selections(
    universe: "Universe",
    protein_group_selections: dict[str, str] | None = None,
) -> dict[str, set[int]]:
    """Resolve protein group MDAnalysis selections to residue IDs.

    Converts MDAnalysis selection strings into sets of residue IDs at
    analysis time. This allows the user to define groups with flexible
    selections while we work with concrete residue IDs internally.

    Parameters
    ----------
    universe : Universe
        MDAnalysis Universe with loaded topology
    protein_group_selections : dict[str, str], optional
        Mapping of group name to MDAnalysis selection string.
        If None, uses DEFAULT_AA_CLASS_SELECTIONS.

    Returns
    -------
    dict[str, set[int]]
        Mapping of group name to set of residue IDs

    Examples
    --------
    >>> selections = {"aromatic": "protein and resname PHE TRP TYR HIS"}
    >>> groups = resolve_protein_group_selections(universe, selections)
    >>> print(groups["aromatic"])
    {12, 45, 67, 89, ...}  # Set of aromatic residue IDs
    """
    if protein_group_selections is None:
        protein_group_selections = DEFAULT_AA_CLASS_SELECTIONS.copy()

    resolved: dict[str, set[int]] = {}

    for group_name, selection in protein_group_selections.items():
        try:
            atoms = universe.select_atoms(selection)
            resids = set(atoms.residues.resids)
            resolved[group_name] = resids
            logger.debug(
                f"Resolved protein group '{group_name}': {len(resids)} residues ({selection})"
            )
        except Exception as e:
            logger.warning(
                f"Failed to resolve protein group '{group_name}' with selection '{selection}': {e}"
            )
            resolved[group_name] = set()

    return resolved


def resolve_protein_groups_from_surface_exposure(
    surface_exposure: "SurfaceExposureResult",
    include_default_aa_groups: bool = True,
    custom_protein_groups: dict[str, list[int]] | None = None,
) -> dict[str, set[int]]:
    """Resolve protein groups from surface exposure data without Universe.

    This function derives protein group → residue ID mappings using only the
    surface exposure result (which already contains resid, resname, and aa_class
    for each residue). This allows binding preference computation without
    requiring an MDAnalysis Universe at comparison time.

    The function supports:
    - Default AA class groups (aromatic, polar, nonpolar, charged_+, charged_-)
      derived from the aa_class field in surface exposure data
    - Custom user-defined groups specified as resid lists
    - Override behavior: if a custom group has the same name as a default,
      the custom definition takes precedence

    Parameters
    ----------
    surface_exposure : SurfaceExposureResult
        Surface exposure analysis result containing residue data
    include_default_aa_groups : bool, default True
        If True, include default AA class groupings (aromatic, polar, etc.)
        derived from surface exposure data
    custom_protein_groups : dict[str, list[int]], optional
        User-defined protein groups as {group_name: [resid1, resid2, ...]}.
        If a group name matches a default AA class, it overrides that default.

    Returns
    -------
    dict[str, set[int]]
        Mapping of group name to set of residue IDs

    Examples
    --------
    >>> from polyzymd.analysis.contacts.surface_exposure import SurfaceExposureFilter
    >>> filter = SurfaceExposureFilter(threshold=0.2)
    >>> surface_result = filter.calculate("enzyme.pdb")
    >>> # Get default AA groups + custom active_site group
    >>> groups = resolve_protein_groups_from_surface_exposure(
    ...     surface_result,
    ...     include_default_aa_groups=True,
    ...     custom_protein_groups={"active_site": [77, 133, 156]}
    ... )
    >>> print(groups.keys())
    dict_keys(['aromatic', 'polar', 'nonpolar', 'charged_positive', 'charged_negative', 'active_site'])
    """
    resolved: dict[str, set[int]] = {}
    all_valid_resids = surface_exposure.all_resids

    # Step 1: Build default AA class groups from surface exposure data
    if include_default_aa_groups:
        # Group residues by their aa_class from surface exposure
        aa_class_groups: dict[str, set[int]] = {}
        for res_exp in surface_exposure.residue_exposures:
            aa_class = res_exp.aa_class
            if aa_class and aa_class != "unknown":
                if aa_class not in aa_class_groups:
                    aa_class_groups[aa_class] = set()
                aa_class_groups[aa_class].add(res_exp.resid)

        for group_name, resids in aa_class_groups.items():
            resolved[group_name] = resids
            logger.debug(
                f"Default AA group '{group_name}': {len(resids)} residues from surface exposure"
            )

    # Step 2: Add custom protein groups (with override behavior)
    if custom_protein_groups:
        for group_name, resid_list in custom_protein_groups.items():
            # Convert to set and validate
            requested_resids = set(resid_list)
            valid_resids = requested_resids & all_valid_resids
            invalid_resids = requested_resids - all_valid_resids

            # Warn about invalid resids
            if invalid_resids:
                logger.warning(
                    f"Custom protein group '{group_name}': {len(invalid_resids)} resids not found "
                    f"in enzyme structure and will be ignored: {sorted(invalid_resids)}"
                )

            # Override if same name as default (user intent takes precedence)
            if group_name in resolved:
                logger.info(f"Custom protein group '{group_name}' overrides default AA class group")

            resolved[group_name] = valid_resids
            logger.debug(f"Custom protein group '{group_name}': {len(valid_resids)} valid residues")

    return resolved


def resolve_polymer_type_selections(
    universe: "Universe",
    polymer_type_selections: dict[str, str] | None = None,
) -> list[str]:
    """Resolve polymer type selections and return list of polymer types.

    If explicit selections are provided, validates them and returns the keys.
    If None, auto-detects polymer types from the standard polymer chain (C).

    Parameters
    ----------
    universe : Universe
        MDAnalysis Universe with loaded topology
    polymer_type_selections : dict[str, str], optional
        Mapping of type name to MDAnalysis selection string.
        If None, auto-detects from "chainID C".

    Returns
    -------
    list[str]
        List of polymer type names (resnames or selection keys)

    Examples
    --------
    >>> # Auto-detect from chain C
    >>> polymer_types = resolve_polymer_type_selections(universe, None)
    >>> print(polymer_types)
    ['SBM', 'EGM']

    >>> # Explicit selections
    >>> selections = {"SBMA": "chainID C and resname SBM"}
    >>> polymer_types = resolve_polymer_type_selections(universe, selections)
    >>> print(polymer_types)
    ['SBMA']
    """
    if polymer_type_selections is not None:
        # Validate selections and return keys
        valid_types = []
        for type_name, selection in polymer_type_selections.items():
            try:
                atoms = universe.select_atoms(selection)
                if len(atoms) > 0:
                    valid_types.append(type_name)
                    logger.debug(f"Polymer type '{type_name}': {len(atoms)} atoms ({selection})")
                else:
                    logger.warning(
                        f"Polymer type '{type_name}' selection matched no atoms: {selection}"
                    )
            except Exception as e:
                logger.warning(f"Invalid selection for polymer type '{type_name}': {e}")
        return valid_types

    # Auto-detect from chain C
    try:
        polymer_atoms = universe.select_atoms("chainID C")
        if len(polymer_atoms) == 0:
            logger.warning("No atoms found in chainID C for polymer auto-detection")
            return []

        # Get unique resnames
        resnames = set(polymer_atoms.residues.resnames)
        logger.debug(f"Auto-detected polymer types from chain C: {resnames}")
        return sorted(resnames)

    except Exception as e:
        logger.warning(f"Failed to auto-detect polymer types: {e}")
        return []


def compute_binding_preference_from_config(
    contact_result: "ContactResult",
    universe: "Universe",
    enzyme_pdb_path: Path | str,
    config: "ContactsConfig",
) -> BindingPreferenceResult:
    """Compute binding preference using ContactsConfig settings.

    This is a convenience function that orchestrates the full binding
    preference computation using configuration from analysis.yaml.
    It:
    1. Calculates surface exposure from enzyme PDB
    2. Resolves protein group selections to residue IDs
    3. Computes binding preference with enrichment ratios

    Parameters
    ----------
    contact_result : ContactResult
        Contact analysis results from trajectory
    universe : Universe
        MDAnalysis Universe (used for resolving selections)
    enzyme_pdb_path : Path or str
        Path to enzyme PDB file for SASA calculation
    config : ContactsConfig
        Configuration from analysis.yaml with binding preference settings

    Returns
    -------
    BindingPreferenceResult
        Binding preference with enrichment ratios

    Notes
    -----
    This function requires rust_sasa_python to be installed for
    surface exposure calculation. Install with:
        pip install rust-sasa-python

    Examples
    --------
    >>> from polyzymd.analysis.config import ContactsConfig
    >>> config = ContactsConfig(
    ...     compute_binding_preference=True,
    ...     surface_exposure_threshold=0.2,
    ... )
    >>> result = compute_binding_preference_from_config(
    ...     contact_result, universe, "enzyme.pdb", config
    ... )
    >>> print(result.enrichment_matrix())
    {'SBM': {'aromatic': 1.45, 'polar': 0.82, ...}, ...}
    """
    from polyzymd.analysis.contacts.surface_exposure import SurfaceExposureFilter

    logger.info("Computing binding preference from config...")

    # Step 1: Calculate surface exposure
    exposure_filter = SurfaceExposureFilter(threshold=config.surface_exposure_threshold)
    surface_exposure = exposure_filter.calculate(enzyme_pdb_path)

    logger.info(
        f"Surface exposure: {surface_exposure.exposed_count}/{surface_exposure.total_count} "
        f"residues exposed (threshold={config.surface_exposure_threshold})"
    )

    # Step 2: Resolve protein group selections to residue IDs
    protein_groups = resolve_protein_group_selections(universe, config.protein_group_selections)

    for group_name, resids in protein_groups.items():
        exposed_in_group = surface_exposure.exposed_in_selection(resids)
        logger.debug(f"  {group_name}: {len(resids)} residues, {len(exposed_in_group)} exposed")

    # Step 3: Get polymer types
    polymer_types = resolve_polymer_type_selections(universe, config.polymer_type_selections)
    logger.info(f"Polymer types for binding preference: {polymer_types}")

    # Step 4: Compute binding preference
    result = compute_binding_preference(
        contact_result=contact_result,
        surface_exposure=surface_exposure,
        protein_groups=protein_groups,
        polymer_types=polymer_types if polymer_types else None,
        protein_group_selections=config.protein_group_selections,
        polymer_type_selections=config.polymer_type_selections,
    )

    return result
