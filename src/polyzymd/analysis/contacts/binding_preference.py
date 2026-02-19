"""Binding preference metrics for polymer-protein contacts.

This module computes enrichment ratios that answer the scientific question:
"Does polymer type X preferentially bind amino acid class Y?"

Enrichment Calculation (Zero-Centered)
--------------------------------------
For each (polymer_type, protein_group) pair:

    contact_share = Σ(contact_frames for exposed residues in group) /
                    Σ(contact_frames for all exposed residues)

Enrichment is normalized by polymer composition and **centered at zero**:

    enrichment = (contact_share / expected_share) - 1

Interpretation (applies to both normalization methods):
- enrichment > 0: Preferential binding (more contacts than expected)
    - +0.5 means "50% more contacts than expected"
- enrichment = 0: Neutral (contact frequency matches random chance)
- enrichment < 0: Avoidance (fewer contacts than expected)
    - -0.3 means "30% fewer contacts than expected"
- enrichment = -1: Complete avoidance (no contacts at all)

Dual Normalization
------------------
Two normalization methods are provided to distinguish chemical affinity
from geometric/steric effects:

**Residue-based** (default, matches experimental concentration ratios):

    expected_by_residue = polymer_residue_count / total_polymer_residues
    enrichment_by_residue = (contact_share / expected_by_residue) - 1

**Atom-based** (accounts for monomer size differences):

    expected_by_atoms = polymer_heavy_atoms / total_polymer_heavy_atoms
    enrichment_by_atoms = (contact_share / expected_by_atoms) - 1

Interpretation of dual metrics:
- Both positive → Strong evidence of chemical preference
- Both negative → Strong evidence of avoidance
- Residue positive, Atom ~0 → Enrichment explained by larger monomer size
- Residue ~0, Atom positive → Smaller monomer "punches above its weight"

The contact-frame weighting ensures that residues contacted for longer durations
contribute proportionally more to the enrichment calculation. A residue contacted
for 60% of the simulation contributes 60x more than one contacted for 1 frame.

Examples
--------
>>> from polyzymd.analysis.contacts.binding_preference import (
...     compute_binding_preference,
...     extract_polymer_composition,
... )
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
>>> # Extract polymer composition from trajectory
>>> polymer_composition = extract_polymer_composition(universe)
>>>
>>> # Compute binding preference
>>> result = compute_binding_preference(
...     contact_result,
...     surface_exposure,
...     protein_groups,
...     polymer_composition,
... )
>>>
>>> # Check enrichment (default: by residue)
>>> print(result.get_enrichment("SBM", "aromatic"))
0.45  # SBMA has 45% more contacts than expected with aromatic residues
>>>
>>> # Compare normalizations
>>> print(result.enrichment_matrix(by="residue"))  # Default
{'SBM': {'aromatic': 0.45, 'charged': -0.18}, 'EGM': {...}}
>>> print(result.enrichment_matrix(by="atoms"))    # Size-adjusted
{'SBM': {'aromatic': 0.12, 'charged': -0.25}, 'EGM': {...}}
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


class PolymerComposition(BaseModel):
    """Polymer composition data extracted from trajectory.

    Contains counts of residues and heavy atoms for each polymer type,
    enabling normalization of enrichment ratios by polymer availability.

    This data is used for dual normalization of binding preference:

    - **Residue-based**: Normalizes by number of polymer residues per type.
      Matches the experimental viewpoint where concentrations are specified
      in terms of monomer units.

    - **Atom-based**: Normalizes by number of heavy atoms (non-hydrogen) per type.
      Accounts for differences in monomer size, since larger monomers have
      more surface area and thus more opportunity for contacts.

    Attributes
    ----------
    residue_counts : dict[str, int]
        Number of residues per polymer type (e.g., {"SBM": 50, "EGM": 50})
    heavy_atom_counts : dict[str, int]
        Number of heavy atoms (non-hydrogen) per polymer type.
        Heavy atoms are defined as all atoms with element != 'H'.

    Examples
    --------
    >>> composition = PolymerComposition(
    ...     residue_counts={"SBM": 50, "EGM": 50},
    ...     heavy_atom_counts={"SBM": 750, "EGM": 400},
    ... )
    >>> composition.total_residues
    100
    >>> composition.residue_fraction("SBM")
    0.5
    >>> composition.heavy_atom_fraction("SBM")
    0.652  # SBM has larger monomers
    """

    residue_counts: dict[str, int] = Field(
        default_factory=dict,
        description="Number of residues per polymer type",
    )
    heavy_atom_counts: dict[str, int] = Field(
        default_factory=dict,
        description="Number of heavy atoms (non-H) per polymer type",
    )

    @property
    def total_residues(self) -> int:
        """Total polymer residues across all types."""
        return sum(self.residue_counts.values())

    @property
    def total_heavy_atoms(self) -> int:
        """Total heavy atoms across all polymer types."""
        return sum(self.heavy_atom_counts.values())

    def residue_fraction(self, polymer_type: str) -> float:
        """Fraction of residues that are this polymer type.

        Parameters
        ----------
        polymer_type : str
            Polymer type name (e.g., "SBM")

        Returns
        -------
        float
            Fraction in range [0, 1], or 0.0 if type not found
        """
        total = self.total_residues
        if total == 0:
            return 0.0
        return self.residue_counts.get(polymer_type, 0) / total

    def heavy_atom_fraction(self, polymer_type: str) -> float:
        """Fraction of heavy atoms that are this polymer type.

        Parameters
        ----------
        polymer_type : str
            Polymer type name (e.g., "SBM")

        Returns
        -------
        float
            Fraction in range [0, 1], or 0.0 if type not found
        """
        total = self.total_heavy_atoms
        if total == 0:
            return 0.0
        return self.heavy_atom_counts.get(polymer_type, 0) / total

    def polymer_types(self) -> list[str]:
        """Get sorted list of polymer types in this composition."""
        return sorted(set(self.residue_counts.keys()) | set(self.heavy_atom_counts.keys()))


class BindingPreferenceEntry(BaseModel):
    """Single entry in the binding preference matrix.

    Represents the binding preference metrics for one
    (polymer_type, protein_group) combination.

    Enrichment Interpretation (centered at zero)
    --------------------------------------------
    - enrichment > 0: Preferential binding (more contacts than expected)
        - +0.5 means "50% more contacts than expected"
    - enrichment = 0: Neutral (contact frequency matches random chance)
    - enrichment < 0: Avoidance (fewer contacts than expected)
        - -0.3 means "30% fewer contacts than expected"
    - enrichment = -1: Complete avoidance (no contacts at all)

    Dual Normalization
    ------------------
    Two normalization methods are provided to distinguish chemical
    affinity from geometric/steric effects:

    - **enrichment_by_residue**: Normalized by residue count. Matches
      experimental concentration ratios. Use this for direct comparison
      with experimental expectations.

    - **enrichment_by_atoms**: Normalized by heavy atom count. Accounts
      for monomer size differences. Use this to identify true chemical
      affinity vs. effects explained by larger monomers having more
      contact surface area.

    Interpretation of agreement/disagreement:
    - Both positive: Strong evidence of chemical preference
    - Both negative: Strong evidence of avoidance
    - Residue positive, Atom ~0: Enrichment explained by larger monomer size
    - Residue ~0, Atom positive: Smaller monomer "punches above its weight"

    Attributes
    ----------
    polymer_type : str
        Polymer residue type (e.g., "SBM", "EGM")
    protein_group : str
        Protein group label (e.g., "aromatic", "charged_positive")
    total_contact_frames : int
        Sum of contact frames across all exposed residues in this group.
    mean_contact_fraction : float
        Average per-residue contact fraction within this group.
    n_residues_in_group : int
        Total residues in this protein group (exposed + buried)
    n_exposed_in_group : int
        Surface-exposed residues in this group (used for enrichment)
    n_residues_contacted : int
        Number of exposed residues that had at least one contact
    contact_share : float
        Fraction of this polymer's total contacts that went to this group.
    polymer_residue_count : int
        Number of residues of this polymer type in the system.
    total_polymer_residues : int
        Total polymer residues across all types.
    expected_share_by_residue : float
        Expected contact share based on polymer residue count.
    enrichment_by_residue : float | None
        Enrichment normalized by residue count: (contact_share / expected) - 1
    polymer_heavy_atom_count : int
        Number of heavy atoms (non-hydrogen) for this polymer type.
    total_polymer_heavy_atoms : int
        Total heavy atoms across all polymer types.
    expected_share_by_atoms : float
        Expected contact share based on heavy atom count.
    enrichment_by_atoms : float | None
        Enrichment normalized by heavy atom count: (contact_share / expected) - 1
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

    # Residue-based normalization (matches experimental ratios)
    polymer_residue_count: int = Field(
        default=0,
        description="Number of residues of this polymer type in the system",
    )
    total_polymer_residues: int = Field(
        default=0,
        description="Total polymer residues across all types",
    )
    expected_share_by_residue: float = Field(
        default=0.0,
        description="Expected contact share based on polymer residue count",
    )
    enrichment_by_residue: float | None = Field(
        default=None,
        description="Enrichment by residue count: (contact_share / expected) - 1",
    )

    # Atom-based normalization (accounts for monomer size)
    polymer_heavy_atom_count: int = Field(
        default=0,
        description="Number of heavy atoms (non-H) for this polymer type",
    )
    total_polymer_heavy_atoms: int = Field(
        default=0,
        description="Total heavy atoms across all polymer types",
    )
    expected_share_by_atoms: float = Field(
        default=0.0,
        description="Expected contact share based on heavy atom count",
    )
    enrichment_by_atoms: float | None = Field(
        default=None,
        description="Enrichment by heavy atoms: (contact_share / expected) - 1",
    )


class BindingPreferenceResult(BaseModel):
    """Complete binding preference analysis result.

    Provides enrichment-normalized metrics for polymer-protein
    binding preferences, answering questions like:

    - "Does SBMA preferentially bind aromatic residues?"
    - "How does EGMA's preference for charged residues compare to SBMA?"
    - "Which amino acid class does this polymer type prefer?"

    Enrichment values are centered at zero:
    - enrichment > 0: Preferential binding
    - enrichment = 0: Neutral (random chance)
    - enrichment < 0: Avoidance

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
    polymer_composition : PolymerComposition
        Polymer composition data (residue/atom counts per type)
    schema_version : int
        Version for forward compatibility. Version 2 introduces
        dual normalization (by residue and by heavy atoms).
    """

    entries: list[BindingPreferenceEntry] = Field(default_factory=list)
    n_frames: int = 0
    total_exposed_residues: int = 0
    surface_exposure_threshold: float | None = None
    protein_groups_used: dict[str, str] = Field(default_factory=dict)
    polymer_types_used: dict[str, str] = Field(default_factory=dict)
    polymer_composition: PolymerComposition | None = Field(
        default=None,
        description="Polymer composition data (residue/atom counts per type)",
    )
    metadata: dict[str, Any] = Field(default_factory=dict)
    schema_version: int = 2  # Version 2: dual normalization

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

    def enrichment_matrix(self, by: str = "residue") -> dict[str, dict[str, float]]:
        """Get enrichment as nested dict: {polymer_type: {protein_group: value}}.

        Enrichment values are centered at zero:
        - > 0: Preferential binding (more contacts than expected)
        - = 0: Neutral (matches random chance)
        - < 0: Avoidance (fewer contacts than expected)

        Parameters
        ----------
        by : str, default "residue"
            Normalization method to use:

            - "residue": Normalize by polymer residue count. Matches
              experimental concentration ratios. This is the default
              because it aligns with how experimentalists typically
              describe polymer compositions.

            - "atoms": Normalize by heavy atom count. Accounts for
              monomer size differences. Use this to distinguish true
              chemical affinity from geometric/steric effects.

        Returns
        -------
        dict[str, dict[str, float]]
            Nested mapping of enrichment values.
            Missing/invalid values are returned as 0.0.

        Examples
        --------
        >>> # Default: residue-based enrichment
        >>> matrix = result.enrichment_matrix()
        >>> print(matrix["SBM"]["aromatic"])
        0.45  # 45% more contacts than expected

        >>> # Atom-based enrichment (accounts for monomer size)
        >>> matrix_atoms = result.enrichment_matrix(by="atoms")
        >>> print(matrix_atoms["SBM"]["aromatic"])
        0.12  # Only 12% enrichment when accounting for size

        Notes
        -----
        To change the default normalization in your analysis, explicitly
        pass the ``by`` parameter. For comparing both normalizations:

        >>> residue_enr = result.enrichment_matrix(by="residue")
        >>> atom_enr = result.enrichment_matrix(by="atoms")
        """
        if by not in ("residue", "atoms"):
            raise ValueError(f"by must be 'residue' or 'atoms', got {by!r}")

        result: dict[str, dict[str, float]] = {}
        for entry in self.entries:
            if entry.polymer_type not in result:
                result[entry.polymer_type] = {}

            if by == "atoms":
                value = entry.enrichment_by_atoms
            else:
                value = entry.enrichment_by_residue

            result[entry.polymer_type][entry.protein_group] = value if value is not None else 0.0
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

    def get_enrichment(
        self, polymer_type: str, protein_group: str, by: str = "residue"
    ) -> float | None:
        """Get enrichment for a specific (polymer_type, protein_group) pair.

        Parameters
        ----------
        polymer_type : str
            Polymer type name
        protein_group : str
            Protein group name
        by : str, default "residue"
            Normalization method: "residue" or "atoms"

        Returns
        -------
        float or None
            Enrichment value (centered at zero), or None if pair not found
        """
        for entry in self.entries:
            if entry.polymer_type == polymer_type and entry.protein_group == protein_group:
                if by == "atoms":
                    return entry.enrichment_by_atoms
                return entry.enrichment_by_residue
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
    polymer_composition: PolymerComposition,
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
    3. Polymer composition (normalization by residue count or heavy atom count)

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
    polymer_composition : PolymerComposition
        Polymer composition data (residue and heavy atom counts per type).
        Used for dual normalization of enrichment ratios.
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
        Binding preference metrics with dual enrichment normalization

    Notes
    -----
    Enrichment Calculation (centered at zero)
    -----------------------------------------
    For each (polymer_type, protein_group) pair:

        contact_share = polymer_contacts_to_group / polymer_total_contacts

    Two normalization methods are computed:

    1. **Residue-based** (matches experimental concentration ratios):

        expected_by_residue = polymer_residue_count / total_polymer_residues
        enrichment_by_residue = (contact_share / expected_by_residue) - 1

    2. **Atom-based** (accounts for monomer size via heavy atoms):

        expected_by_atoms = polymer_heavy_atoms / total_polymer_heavy_atoms
        enrichment_by_atoms = (contact_share / expected_by_atoms) - 1

    Interpretation (both methods):
    - enrichment > 0: Preferential binding (more contacts than expected)
    - enrichment = 0: Neutral (matches random chance)
    - enrichment < 0: Avoidance (fewer contacts than expected)
    - enrichment = -1: Complete avoidance (no contacts at all)

    When to use which metric:
    - enrichment_by_residue: Direct comparison to experimental concentration ratios
    - enrichment_by_atoms: Reveals true chemical affinity vs. geometric/steric effects
    """
    exposed_resids = surface_exposure.exposed_resids
    n_frames = contact_result.n_frames
    total_exposed = len(exposed_resids)

    logger.info(
        f"Computing binding preference: {total_exposed} exposed residues, "
        f"{len(protein_groups)} protein groups"
    )
    logger.info(
        f"Polymer composition: {polymer_composition.residue_counts} residues, "
        f"{polymer_composition.heavy_atom_counts} heavy atoms"
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

    # Helper function for enrichment calculation (centered at zero)
    def calc_enrichment(contact_share: float, expected: float) -> float | None:
        """Calculate enrichment: (observed / expected) - 1, centered at zero."""
        if expected > 0 and contact_share > 0:
            return (contact_share / expected) - 1
        elif expected > 0 and contact_share == 0:
            return -1.0  # Complete avoidance
        else:
            return None  # Cannot compute (no polymer of this type)

    # Get polymer composition totals
    total_poly_residues = polymer_composition.total_residues
    total_poly_atoms = polymer_composition.total_heavy_atoms

    # Build result entries with enrichment calculations
    entries = []

    for poly_type in sorted(contact_data.keys()):
        total_poly_contacts = total_contacts_by_polymer.get(poly_type, 0)

        # Get polymer composition for this type
        poly_res_count = polymer_composition.residue_counts.get(poly_type, 0)
        poly_atom_count = polymer_composition.heavy_atom_counts.get(poly_type, 0)

        # Calculate expected shares by each normalization method
        if total_poly_residues > 0:
            expected_share_by_residue = poly_res_count / total_poly_residues
        else:
            expected_share_by_residue = 0.0

        if total_poly_atoms > 0:
            expected_share_by_atoms = poly_atom_count / total_poly_atoms
        else:
            expected_share_by_atoms = 0.0

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

            # Calculate contact share
            if total_poly_contacts > 0:
                contact_share = contact_frames / total_poly_contacts
            else:
                contact_share = 0.0

            # Calculate enrichment (centered at zero) for both normalizations
            enrichment_by_residue = calc_enrichment(contact_share, expected_share_by_residue)
            enrichment_by_atoms = calc_enrichment(contact_share, expected_share_by_atoms)

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
                    # Residue-based normalization
                    polymer_residue_count=poly_res_count,
                    total_polymer_residues=total_poly_residues,
                    expected_share_by_residue=expected_share_by_residue,
                    enrichment_by_residue=enrichment_by_residue,
                    # Atom-based normalization
                    polymer_heavy_atom_count=poly_atom_count,
                    total_polymer_heavy_atoms=total_poly_atoms,
                    expected_share_by_atoms=expected_share_by_atoms,
                    enrichment_by_atoms=enrichment_by_atoms,
                )
            )

    result = BindingPreferenceResult(
        entries=entries,
        n_frames=n_frames,
        total_exposed_residues=total_exposed,
        surface_exposure_threshold=surface_exposure.threshold,
        protein_groups_used=protein_group_selections or {},
        polymer_types_used=polymer_type_selections or {},
        polymer_composition=polymer_composition,
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

    Computes mean ± SEM for both residue-based and atom-based enrichment
    ratios across multiple replicates.

    Parameters
    ----------
    results : list[BindingPreferenceResult]
        Binding preference results from multiple replicates

    Returns
    -------
    AggregatedBindingPreferenceResult
        Aggregated results with mean and SEM for both normalization methods
    """
    if not results:
        raise ValueError("No results to aggregate")

    # Helper function for computing mean and SEM
    def _compute_stats(values: list[float]) -> tuple[float | None, float | None]:
        """Compute mean and SEM from a list of values."""
        n = len(values)
        if n == 0:
            return None, None
        mean_val = float(np.mean(values))
        sem_val = float(np.std(values, ddof=1) / np.sqrt(n)) if n > 1 else 0.0
        return mean_val, sem_val

    # Collect all (polymer_type, protein_group) pairs
    all_pairs: set[tuple[str, str]] = set()
    for r in results:
        for e in r.entries:
            all_pairs.add((e.polymer_type, e.protein_group))

    entries = []
    for poly_type, prot_group in sorted(all_pairs):
        # Collect values from each replicate
        enrichments_by_residue = []
        enrichments_by_atoms = []
        contact_fractions = []
        contact_shares = []

        for r in results:
            entry = r.get_entry(poly_type, prot_group)
            if entry is not None:
                if entry.enrichment_by_residue is not None:
                    enrichments_by_residue.append(entry.enrichment_by_residue)
                if entry.enrichment_by_atoms is not None:
                    enrichments_by_atoms.append(entry.enrichment_by_atoms)
                contact_fractions.append(entry.mean_contact_fraction)
                contact_shares.append(entry.contact_share)

        # Compute statistics for both normalizations
        mean_enr_res, sem_enr_res = _compute_stats(enrichments_by_residue)
        mean_enr_atoms, sem_enr_atoms = _compute_stats(enrichments_by_atoms)
        mean_contact_fraction, sem_contact_fraction = _compute_stats(contact_fractions)
        mean_contact_share = float(np.mean(contact_shares)) if contact_shares else 0.0

        # Get group metadata from first result
        first_entry = results[0].get_entry(poly_type, prot_group)
        n_exposed = first_entry.n_exposed_in_group if first_entry else 0
        n_total = first_entry.n_residues_in_group if first_entry else 0
        expected_by_residue = first_entry.expected_share_by_residue if first_entry else 0.0
        expected_by_atoms = first_entry.expected_share_by_atoms if first_entry else 0.0

        entries.append(
            AggregatedBindingPreferenceEntry(
                polymer_type=poly_type,
                protein_group=prot_group,
                # Residue-based enrichment
                mean_enrichment_by_residue=mean_enr_res,
                sem_enrichment_by_residue=sem_enr_res,
                per_replicate_enrichments_by_residue=enrichments_by_residue,
                # Atom-based enrichment
                mean_enrichment_by_atoms=mean_enr_atoms,
                sem_enrichment_by_atoms=sem_enr_atoms,
                per_replicate_enrichments_by_atoms=enrichments_by_atoms,
                # Contact metrics
                mean_contact_fraction=mean_contact_fraction if mean_contact_fraction else 0.0,
                sem_contact_fraction=sem_contact_fraction if sem_contact_fraction else 0.0,
                mean_contact_share=mean_contact_share,
                # Expected shares
                expected_share_by_residue=expected_by_residue,
                expected_share_by_atoms=expected_by_atoms,
                # Group metadata
                n_exposed_in_group=n_exposed,
                n_residues_in_group=n_total,
                n_replicates=len(enrichments_by_residue),
            )
        )

    return AggregatedBindingPreferenceResult(
        entries=entries,
        n_replicates=len(results),
        total_exposed_residues=results[0].total_exposed_residues if results else 0,
        surface_exposure_threshold=results[0].surface_exposure_threshold if results else None,
        protein_groups_used=results[0].protein_groups_used if results else {},
        polymer_types_used=results[0].polymer_types_used if results else {},
        polymer_composition=results[0].polymer_composition if results else None,
    )


class AggregatedBindingPreferenceEntry(BaseModel):
    """Aggregated binding preference for one (polymer_type, protein_group) pair.

    Contains mean ± SEM across replicates for both residue-based and
    atom-based enrichment normalizations.

    Enrichment values are centered at zero:
    - > 0: Preferential binding
    - = 0: Neutral (random chance)
    - < 0: Avoidance
    """

    polymer_type: str
    protein_group: str

    # Residue-based enrichment (matches experimental ratios)
    mean_enrichment_by_residue: float | None = Field(
        default=None,
        description="Mean residue-normalized enrichment across replicates",
    )
    sem_enrichment_by_residue: float | None = Field(
        default=None,
        description="Standard error of residue-normalized enrichment",
    )
    per_replicate_enrichments_by_residue: list[float] = Field(
        default_factory=list,
        description="Residue-normalized enrichment values from each replicate",
    )

    # Atom-based enrichment (accounts for monomer size)
    mean_enrichment_by_atoms: float | None = Field(
        default=None,
        description="Mean atom-normalized enrichment across replicates",
    )
    sem_enrichment_by_atoms: float | None = Field(
        default=None,
        description="Standard error of atom-normalized enrichment",
    )
    per_replicate_enrichments_by_atoms: list[float] = Field(
        default_factory=list,
        description="Atom-normalized enrichment values from each replicate",
    )

    # Contact metrics
    mean_contact_fraction: float = Field(
        default=0.0,
        description="Mean per-residue contact fraction",
    )
    sem_contact_fraction: float = Field(
        default=0.0,
        description="Standard error of contact fraction",
    )
    mean_contact_share: float = Field(
        default=0.0,
        description="Mean contact share",
    )

    # Expected shares (from polymer composition)
    expected_share_by_residue: float = Field(
        default=0.0,
        description="Expected contact share based on polymer residue count",
    )
    expected_share_by_atoms: float = Field(
        default=0.0,
        description="Expected contact share based on heavy atom count",
    )

    # Group metadata
    n_exposed_in_group: int = Field(
        default=0,
        description="Surface-exposed residues in group",
    )
    n_residues_in_group: int = Field(
        default=0,
        description="Total residues in group",
    )
    n_replicates: int = Field(
        default=0,
        description="Number of replicates with valid data",
    )


class AggregatedBindingPreferenceResult(BaseModel):
    """Binding preference aggregated across replicates.

    Contains mean ± SEM for all metrics across multiple replicates,
    with dual enrichment normalization (by residue and by heavy atoms).

    Enrichment values are centered at zero:
    - > 0: Preferential binding
    - = 0: Neutral (random chance)
    - < 0: Avoidance
    """

    entries: list[AggregatedBindingPreferenceEntry] = Field(default_factory=list)
    n_replicates: int = 0
    total_exposed_residues: int = 0
    surface_exposure_threshold: float | None = None
    protein_groups_used: dict[str, str] = Field(default_factory=dict)
    polymer_types_used: dict[str, str] = Field(default_factory=dict)
    polymer_composition: PolymerComposition | None = Field(
        default=None,
        description="Polymer composition data (residue/atom counts per type)",
    )
    schema_version: int = 2  # Version 2: dual normalization

    def to_dataframe(self) -> "pd.DataFrame":
        """Convert to pandas DataFrame."""
        import pandas as pd

        return pd.DataFrame([e.model_dump() for e in self.entries])

    def enrichment_matrix(self, by: str = "residue") -> dict[str, dict[str, float]]:
        """Get mean enrichment as nested dict.

        Parameters
        ----------
        by : str, default "residue"
            Normalization method: "residue" or "atoms"

        Returns
        -------
        dict[str, dict[str, float]]
            Nested mapping: {polymer_type: {protein_group: mean_enrichment}}
        """
        if by not in ("residue", "atoms"):
            raise ValueError(f"by must be 'residue' or 'atoms', got {by!r}")

        result: dict[str, dict[str, float]] = {}
        for entry in self.entries:
            if entry.polymer_type not in result:
                result[entry.polymer_type] = {}

            if by == "atoms":
                value = entry.mean_enrichment_by_atoms
            else:
                value = entry.mean_enrichment_by_residue

            result[entry.polymer_type][entry.protein_group] = value if value is not None else 0.0
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


def extract_polymer_composition(
    universe: "Universe",
    polymer_type_selections: dict[str, str] | None = None,
) -> PolymerComposition:
    """Extract polymer composition (residue and heavy atom counts) from trajectory.

    This function counts the number of residues and heavy atoms for each
    polymer type in the system. The composition data is used for dual
    normalization of binding preference enrichment ratios.

    Parameters
    ----------
    universe : Universe
        MDAnalysis Universe with loaded topology
    polymer_type_selections : dict[str, str], optional
        Mapping of type name to MDAnalysis selection string.
        If None, auto-detects from "chainID C" using unique resnames.

    Returns
    -------
    PolymerComposition
        Composition data with residue_counts and heavy_atom_counts per type.
        Heavy atoms are defined as all atoms with element != 'H'.

    Notes
    -----
    For auto-detection (polymer_type_selections=None), the function:
    1. Selects all atoms in chainID C
    2. Groups residues by resname
    3. Counts residues and heavy atoms per resname

    For explicit selections, each selection string defines which atoms
    belong to that polymer type. The function counts residues (unique
    resids) and heavy atoms within each selection.

    Examples
    --------
    >>> # Auto-detect from chain C
    >>> composition = extract_polymer_composition(universe)
    >>> print(composition.residue_counts)
    {'SBM': 50, 'EGM': 50}
    >>> print(composition.heavy_atom_counts)
    {'SBM': 750, 'EGM': 400}

    >>> # Explicit selections
    >>> selections = {"SBMA": "chainID C and resname SBM"}
    >>> composition = extract_polymer_composition(universe, selections)
    """
    residue_counts: dict[str, int] = {}
    heavy_atom_counts: dict[str, int] = {}

    if polymer_type_selections is not None:
        # Use explicit selections
        for type_name, selection in polymer_type_selections.items():
            try:
                atoms = universe.select_atoms(selection)
                if len(atoms) == 0:
                    logger.warning(
                        f"Polymer type '{type_name}' selection matched no atoms: {selection}"
                    )
                    continue

                # Count unique residues
                n_residues = len(atoms.residues)
                residue_counts[type_name] = n_residues

                # Count heavy atoms (non-hydrogen)
                # MDAnalysis stores element info - filter out hydrogens
                try:
                    heavy_atoms = atoms.select_atoms("not element H")
                    n_heavy = len(heavy_atoms)
                except Exception:
                    # Fallback: count atoms with mass > 1.1 (heavier than H)
                    n_heavy = sum(1 for a in atoms if a.mass > 1.1)

                heavy_atom_counts[type_name] = n_heavy

                logger.debug(
                    f"Polymer type '{type_name}': {n_residues} residues, {n_heavy} heavy atoms"
                )

            except Exception as e:
                logger.warning(f"Failed to extract composition for '{type_name}': {e}")

    else:
        # Auto-detect from chain C
        try:
            polymer_atoms = universe.select_atoms("chainID C")
            if len(polymer_atoms) == 0:
                logger.warning("No atoms found in chainID C for polymer composition")
                return PolymerComposition()

            # Group by resname
            resnames = set(polymer_atoms.residues.resnames)

            for resname in resnames:
                # Select atoms of this resname within chain C
                type_atoms = universe.select_atoms(f"chainID C and resname {resname}")

                # Count residues
                n_residues = len(type_atoms.residues)
                residue_counts[resname] = n_residues

                # Count heavy atoms
                try:
                    heavy_atoms = type_atoms.select_atoms("not element H")
                    n_heavy = len(heavy_atoms)
                except Exception:
                    # Fallback: count atoms with mass > 1.1
                    n_heavy = sum(1 for a in type_atoms if a.mass > 1.1)

                heavy_atom_counts[resname] = n_heavy

                logger.debug(
                    f"Auto-detected polymer '{resname}': {n_residues} residues, "
                    f"{n_heavy} heavy atoms"
                )

        except Exception as e:
            logger.warning(f"Failed to extract polymer composition: {e}")
            return PolymerComposition()

    composition = PolymerComposition(
        residue_counts=residue_counts,
        heavy_atom_counts=heavy_atom_counts,
    )

    logger.info(
        f"Polymer composition extracted: {composition.total_residues} total residues, "
        f"{composition.total_heavy_atoms} total heavy atoms across {len(residue_counts)} types"
    )

    return composition


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

    # Step 4: Extract polymer composition for dual normalization
    polymer_composition = extract_polymer_composition(universe, config.polymer_type_selections)

    # Step 5: Compute binding preference
    result = compute_binding_preference(
        contact_result=contact_result,
        surface_exposure=surface_exposure,
        protein_groups=protein_groups,
        polymer_composition=polymer_composition,
        polymer_types=polymer_types if polymer_types else None,
        protein_group_selections=config.protein_group_selections,
        polymer_type_selections=config.polymer_type_selections,
    )

    return result
