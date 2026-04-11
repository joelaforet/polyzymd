"""Helpers for resolving protein and polymer selections."""

from __future__ import annotations

import logging
from typing import TYPE_CHECKING

from polyzymd.analyses.shared.aa_classification import DEFAULT_AA_CLASS_SELECTIONS

if TYPE_CHECKING:
    from MDAnalysis.core.universe import Universe

    from polyzymd.analyses.shared.surface_exposure import SurfaceExposureResult

logger = logging.getLogger(__name__)


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
        except (ValueError, AttributeError, TypeError) as e:
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
    >>> from polyzymd.analyses.shared.surface_exposure import SurfaceExposureFilter
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
    polymer_chain: str = "C",
) -> list[str]:
    """Resolve polymer type selections and return list of polymer types.

    If explicit selections are provided, validates them and returns the keys.
    If None, auto-detects polymer types from the standard polymer chain.

    Parameters
    ----------
    universe : Universe
        MDAnalysis Universe with loaded topology
    polymer_type_selections : dict[str, str], optional
        Mapping of type name to MDAnalysis selection string.
        If None, auto-detects from ``chainID <polymer_chain>``.
    polymer_chain : str
        Chain ID for auto-detection when *polymer_type_selections* is None.
        Defaults to ``"C"`` (PolyzyMD chain convention).

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
            except (ValueError, AttributeError, TypeError) as e:
                logger.warning(f"Invalid selection for polymer type '{type_name}': {e}")
        return valid_types

    # Auto-detect from polymer chain
    try:
        polymer_atoms = universe.select_atoms(f"chainID {polymer_chain}")
        if len(polymer_atoms) == 0:
            logger.warning(f"No atoms found in chainID {polymer_chain} for polymer auto-detection")
            return []

        # Get unique resnames
        resnames = set(polymer_atoms.residues.resnames)
        logger.debug(f"Auto-detected polymer types from chain {polymer_chain}: {resnames}")
        return sorted(resnames)

    except (ValueError, AttributeError, TypeError) as e:
        logger.warning(f"Failed to auto-detect polymer types: {e}")
        return []


