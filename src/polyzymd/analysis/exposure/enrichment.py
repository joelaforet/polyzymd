"""Dynamic chaperone enrichment with dual residue/atom normalization.

Enrichment measures whether a polymer type contacts transiently-exposed
residues more (or less) than would be expected by chance.

Formula (time-averaged, Option B from the design spec)
-------------------------------------------------------
For a polymer type P and residue group G (e.g. "aromatic"):

    observed[t]  = fraction of residues in G ∩ {exposed at t} that P contacts
    expected[t]  = |G ∩ {exposed at t}| / |{all protein residues exposed at t}|

    enrichment   = mean_t(observed[t]) / mean_t(expected[t]) - 1

A value of 0 means the polymer contacts transiently exposed residues in
proportion to their abundance among all exposed residues (random).
Positive enrichment → preferential contact with that group when exposed.

Dual normalization (Issue #33 design decision)
----------------------------------------------
- **Residue-based** (primary): each residue counts as 1 regardless of size.
  Matches the experimentalist framing ("90:10 SBMA:EGMA by monomer fraction").
- **Atom-based** (secondary): weights each residue by its number of heavy atoms.
  Accounts for the fact that larger monomers geometrically expose more atoms.

``size_contribution`` measures the fraction of enrichment attributable to
atom count differences vs. chemical preference.

Scope
-----
This module operates on a **single trajectory** (one replicate).  Aggregation
across replicates is done upstream by the comparator.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from typing import TYPE_CHECKING

import numpy as np
from numpy.typing import NDArray

if TYPE_CHECKING:
    from polyzymd.analysis.contacts.results import ContactResult
    from polyzymd.analysis.sasa.trajectory import SASATrajectoryResult

logger = logging.getLogger(__name__)

# Small constant to avoid division by zero
_EPS = 1e-12


# ---------------------------------------------------------------------------
# Per-polymer-type, per-group enrichment entry
# ---------------------------------------------------------------------------


@dataclass
class GroupEnrichmentEntry:
    """Enrichment of one polymer type toward one amino-acid group.

    Attributes
    ----------
    polymer_type : str
        Polymer monomer resname (e.g. "SBM", "EGM").
    aa_group : str
        Amino-acid group label (e.g. "aromatic", "charged_positive").
    n_residues_in_group : int
        Number of protein residues in the group.
    enrichment_residue : float
        Enrichment relative to residue-level expectation.
    enrichment_atom : float
        Enrichment relative to atom-level expectation.
    size_contribution : float
        (enrichment_atom - enrichment_residue) / (|enrichment_residue| + eps).
        Positive → atoms-per-residue drives enrichment above the residue signal.
    observed_contact_fraction : float
        Mean fraction of exposed-group residues contacted by this polymer type.
    expected_contact_fraction_residue : float
        Mean expected fraction (residue-based).
    expected_contact_fraction_atom : float
        Mean expected fraction (atom-based).
    n_frames_with_exposed : int
        Number of frames where >= 1 group residue was exposed.
    """

    polymer_type: str
    aa_group: str
    n_residues_in_group: int
    enrichment_residue: float
    enrichment_atom: float
    size_contribution: float
    observed_contact_fraction: float
    expected_contact_fraction_residue: float
    expected_contact_fraction_atom: float
    n_frames_with_exposed: int


@dataclass
class ChaperoneEnrichmentResult:
    """Dynamic chaperone enrichment for all polymer types and AA groups.

    Attributes
    ----------
    entries : list[GroupEnrichmentEntry]
        One entry per (polymer_type, aa_group) pair.
    polymer_types : list[str]
        All polymer types present.
    aa_groups : list[str]
        All AA groups present in the protein.
    n_frames : int
        Number of frames analyzed.
    trajectory_path : str
        Source trajectory path (provenance).
    """

    entries: list[GroupEnrichmentEntry] = field(default_factory=list)
    polymer_types: list[str] = field(default_factory=list)
    aa_groups: list[str] = field(default_factory=list)
    n_frames: int = 0
    trajectory_path: str = ""

    def get(self, polymer_type: str, aa_group: str) -> GroupEnrichmentEntry | None:
        """Return the entry for a specific (polymer_type, aa_group) pair."""
        for e in self.entries:
            if e.polymer_type == polymer_type and e.aa_group == aa_group:
                return e
        return None

    def enrichment_matrix_residue(self) -> NDArray[np.float64]:
        """Return enrichment matrix (polymer_types × aa_groups), residue-based.

        Rows correspond to ``self.polymer_types``, columns to ``self.aa_groups``.
        Missing pairs are filled with NaN.
        """
        mat = np.full((len(self.polymer_types), len(self.aa_groups)), np.nan)
        pt_idx = {pt: i for i, pt in enumerate(self.polymer_types)}
        ag_idx = {ag: j for j, ag in enumerate(self.aa_groups)}
        for e in self.entries:
            i, j = pt_idx.get(e.polymer_type, -1), ag_idx.get(e.aa_group, -1)
            if i >= 0 and j >= 0:
                mat[i, j] = e.enrichment_residue
        return mat

    def enrichment_matrix_atom(self) -> NDArray[np.float64]:
        """Return enrichment matrix (polymer_types × aa_groups), atom-based."""
        mat = np.full((len(self.polymer_types), len(self.aa_groups)), np.nan)
        pt_idx = {pt: i for i, pt in enumerate(self.polymer_types)}
        ag_idx = {ag: j for j, ag in enumerate(self.aa_groups)}
        for e in self.entries:
            i, j = pt_idx.get(e.polymer_type, -1), ag_idx.get(e.aa_group, -1)
            if i >= 0 and j >= 0:
                mat[i, j] = e.enrichment_atom
        return mat


# ---------------------------------------------------------------------------
# Atom count lookup (approximate heavy-atom counts per residue type)
# ---------------------------------------------------------------------------

# Approximate heavy-atom counts (from standard amino-acid structures).
# Used only for relative weighting; exact values are not critical.
_HEAVY_ATOM_COUNTS: dict[str, int] = {
    "ALA": 5,
    "ARG": 11,
    "ASN": 8,
    "ASP": 8,
    "CYS": 6,
    "GLN": 9,
    "GLU": 9,
    "GLY": 4,
    "HIS": 10,
    "ILE": 8,
    "LEU": 8,
    "LYS": 9,
    "MET": 8,
    "PHE": 11,
    "PRO": 7,
    "SER": 6,
    "THR": 7,
    "TRP": 14,
    "TYR": 12,
    "VAL": 7,
    # Common non-standard / protonation states
    "HID": 10,
    "HIE": 10,
    "HIP": 10,
    "CYX": 6,
    "LYN": 9,
    "ASH": 8,
    "GLH": 9,
}
_DEFAULT_HEAVY_ATOMS = 8  # fallback


def _heavy_atoms(resname: str) -> int:
    return _HEAVY_ATOM_COUNTS.get(resname.upper(), _DEFAULT_HEAVY_ATOMS)


# ---------------------------------------------------------------------------
# Core computation
# ---------------------------------------------------------------------------


def compute_chaperone_enrichment(
    sasa_result: "SASATrajectoryResult",
    contact_result: "ContactResult",
    polymer_resnames: list[str] | None = None,
) -> ChaperoneEnrichmentResult:
    """Compute time-averaged chaperone enrichment for a single trajectory.

    Parameters
    ----------
    sasa_result : SASATrajectoryResult
        Per-frame SASA result (protein chain only).
    contact_result : ContactResult
        Contact analysis result for the same trajectory.
    polymer_resnames : list[str], optional
        Subset of polymer monomer resnames to include.  If None, all polymer
        types present in ``contact_result`` are used.

    Returns
    -------
    ChaperoneEnrichmentResult
    """
    n_frames = sasa_result.n_frames
    n_residues = sasa_result.n_residues

    # ------------------------------------------------------------------ #
    # Step 1: Build per-residue arrays                                     #
    # ------------------------------------------------------------------ #

    resids = list(sasa_result.resids)
    resnames_list = sasa_result.resnames
    aa_classes = sasa_result.aa_classes
    exposed_mask = sasa_result.exposed_mask_per_frame()  # (n_frames, n_residues)

    # Heavy atom counts per residue — for atom-weighted expected
    atom_weights = np.array([_heavy_atoms(rn) for rn in resnames_list], dtype=np.float64)

    # ------------------------------------------------------------------ #
    # Step 2: Build per-polymer-type contact arrays (n_frames, n_residues) #
    # ------------------------------------------------------------------ #

    # Discover polymer types present in the contact result
    all_polymer_types: set[str] = set()
    for rc in contact_result.residue_contacts:
        for sc in rc.segment_contacts:
            all_polymer_types.add(sc.polymer_resname)

    if polymer_resnames is not None:
        polymer_types = sorted(set(polymer_resnames) & all_polymer_types)
    else:
        polymer_types = sorted(all_polymer_types)

    if not polymer_types:
        logger.warning("No polymer types found in contact result; returning empty enrichment")
        return ChaperoneEnrichmentResult(
            n_frames=n_frames,
            trajectory_path=sasa_result.trajectory_path,
        )

    # (n_polymer_types, n_frames, n_residues) contact boolean arrays
    # Build residue-index lookup for fast access
    resid_to_idx: dict[int, int] = {int(rid): i for i, rid in enumerate(resids)}

    # contact_by_type[pt_idx, frame, res_idx] = True if contact
    contact_by_type = np.zeros((len(polymer_types), n_frames, n_residues), dtype=bool)

    for rc in contact_result.residue_contacts:
        res_idx = resid_to_idx.get(rc.protein_resid)
        if res_idx is None:
            continue
        for sc in rc.segment_contacts:
            if sc.polymer_resname not in polymer_types:
                continue
            pt_idx = polymer_types.index(sc.polymer_resname)
            arr = sc.to_binary_array(n_frames)
            contact_by_type[pt_idx, :, res_idx] |= arr

    # ------------------------------------------------------------------ #
    # Step 3: Group residues by AA class                                   #
    # ------------------------------------------------------------------ #

    aa_groups_set: set[str] = set(aa_classes)
    aa_groups = sorted(aa_groups_set)

    # group_mask[group_idx, res_idx] = True
    group_masks: dict[str, NDArray[np.bool_]] = {}
    for ag in aa_groups:
        mask = np.array([c == ag for c in aa_classes], dtype=bool)
        group_masks[ag] = mask

    # ------------------------------------------------------------------ #
    # Step 4: Compute enrichment per (polymer_type, aa_group)              #
    # ------------------------------------------------------------------ #

    entries: list[GroupEnrichmentEntry] = []

    # Total exposed residues per frame: shape (n_frames,)
    total_exposed_per_frame = exposed_mask.sum(axis=1).astype(np.float64)  # (n_frames,)
    # Atom-weighted total exposed per frame
    total_exposed_atom_per_frame = exposed_mask.astype(np.float64) @ atom_weights  # (n_frames,)

    for pt_idx, ptype in enumerate(polymer_types):
        contact_arr = contact_by_type[pt_idx]  # (n_frames, n_residues)

        for ag in aa_groups:
            gmask = group_masks[ag]  # (n_residues,)
            n_in_group = int(gmask.sum())
            if n_in_group == 0:
                continue

            # Group-exposed mask: (n_frames, n_in_group)
            group_exposed = exposed_mask[:, gmask]  # (n_frames, n_group)
            # Group contacts: (n_frames, n_group)
            group_contact = contact_arr[:, gmask]

            # Frames where at least one group residue is exposed
            frames_with_exposed = group_exposed.any(axis=1)  # (n_frames,)
            n_frames_with_exposed = int(frames_with_exposed.sum())

            if n_frames_with_exposed == 0:
                continue  # this group never exposed → skip

            # --- Observed (same for both normalizations) ---
            # observed[t] = fraction of (exposed group residues) that are contacted
            group_exposed_float = group_exposed.astype(np.float64)
            group_contact_and_exposed = (group_contact & group_exposed).astype(np.float64)

            n_exposed_in_group_per_frame = group_exposed_float.sum(axis=1)  # (n_frames,)
            n_contacted_in_group_per_frame = group_contact_and_exposed.sum(axis=1)

            # Only count frames where group has exposed residues
            valid = n_exposed_in_group_per_frame > 0
            obs_per_frame = np.where(
                valid,
                n_contacted_in_group_per_frame / np.maximum(n_exposed_in_group_per_frame, _EPS),
                0.0,
            )
            mean_observed = float(obs_per_frame[valid].mean()) if valid.any() else 0.0

            # --- Expected: residue-based ---
            # expected_residue[t] = n_exposed_in_group[t] / n_total_exposed[t]
            exp_residue_per_frame = np.where(
                total_exposed_per_frame > 0,
                n_exposed_in_group_per_frame / np.maximum(total_exposed_per_frame, _EPS),
                0.0,
            )
            valid_for_exp = (total_exposed_per_frame > 0) & valid
            mean_expected_residue = (
                float(exp_residue_per_frame[valid_for_exp].mean()) if valid_for_exp.any() else 0.0
            )

            # --- Expected: atom-based ---
            group_atom_weights = atom_weights[gmask]  # (n_group,)
            n_exposed_atom_in_group_per_frame = (
                group_exposed_float @ group_atom_weights
            )  # (n_frames,)
            exp_atom_per_frame = np.where(
                total_exposed_atom_per_frame > 0,
                n_exposed_atom_in_group_per_frame / np.maximum(total_exposed_atom_per_frame, _EPS),
                0.0,
            )
            mean_expected_atom = (
                float(exp_atom_per_frame[valid_for_exp].mean()) if valid_for_exp.any() else 0.0
            )

            # --- Enrichment ---
            enrichment_residue = (
                mean_observed / max(mean_expected_residue, _EPS) - 1.0
                if mean_expected_residue > _EPS
                else 0.0
            )
            enrichment_atom = (
                mean_observed / max(mean_expected_atom, _EPS) - 1.0
                if mean_expected_atom > _EPS
                else 0.0
            )

            # size_contribution: how much atom-based enrichment exceeds residue-based
            size_contribution = (enrichment_atom - enrichment_residue) / (
                abs(enrichment_residue) + _EPS
            )

            entries.append(
                GroupEnrichmentEntry(
                    polymer_type=ptype,
                    aa_group=ag,
                    n_residues_in_group=n_in_group,
                    enrichment_residue=enrichment_residue,
                    enrichment_atom=enrichment_atom,
                    size_contribution=size_contribution,
                    observed_contact_fraction=mean_observed,
                    expected_contact_fraction_residue=mean_expected_residue,
                    expected_contact_fraction_atom=mean_expected_atom,
                    n_frames_with_exposed=n_frames_with_exposed,
                )
            )

    logger.debug(
        f"Enrichment computed: {len(entries)} (polymer_type, aa_group) pairs, "
        f"{len(polymer_types)} polymer types, {len(aa_groups)} AA groups"
    )

    return ChaperoneEnrichmentResult(
        entries=entries,
        polymer_types=polymer_types,
        aa_groups=aa_groups,
        n_frames=n_frames,
        trajectory_path=sasa_result.trajectory_path,
    )
