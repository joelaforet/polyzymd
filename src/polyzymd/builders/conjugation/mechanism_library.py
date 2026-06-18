"""Built-in declarative reaction mechanisms for conjugation planning."""

from __future__ import annotations

from typing import Final

from polyzymd.builders.conjugation.mechanisms import ReactionMechanism

BUILTIN_MECHANISM_DEFINITIONS: Final[tuple[dict, ...]] = (
    {
        "identifier": "nhs_lys_amide",
        "display_name": "NHS ester to lysine amide",
        "allowed_sites": [
            {
                "residue_name": "LYS",
                "atom_name": "NZ",
                "role": "nucleophile",
                "rationale": "Lysine NZ attacks the activated NHS ester carbonyl",
            }
        ],
        "moiety_reactive_group": {
            "group": "nhs_ester",
            "anchor_atom": "C1",
            "leaving_atoms": ["O_NHS"],
            "rationale": "The acyl carbon becomes the amide carbonyl in the final linkage",
        },
        "graph_edits": {
            "delete_site_atoms": [],
            "delete_moiety_atoms": ["O_NHS"],
            "remove_site_hydrogens": ["HZ1"],
            "remove_moiety_protons": [],
            "add_bonds": [{"site_atom": "NZ", "moiety_atom": "C1", "order": 1}],
        },
        "charge_patch_hint": {
            "strategy": "local_junction_patch",
            "patch_radius_bonds": 3,
            "rationale": "Amide formation changes the local lysine and linker charge model",
        },
        "rationale": "Represents NHS ester conjugation to lysine as an amide linkage without "
        "executing graph surgery in this phase",
        "notes": [
            "Atom labels are placeholders until moiety templates define canonical linker names",
            "Full proton and charge reconciliation is deferred to a later phase",
        ],
    },
    {
        "identifier": "n_glycosidic_asn",
        "display_name": "Asparagine N-glycosidic placeholder",
        "allowed_sites": [
            {
                "residue_name": "ASN",
                "atom_name": "ND2",
                "role": "glycosylation_site",
                "rationale": "N-linked glycans attach through the ASN side-chain amide nitrogen",
            }
        ],
        "moiety_reactive_group": {
            "group": "glycan_root",
            "anchor_atom": "C1",
            "leaving_atoms": ["O1"],
            "rationale": "The root sugar anomeric carbon anchors the glycosidic bond",
        },
        "graph_edits": {
            "delete_site_atoms": [],
            "delete_moiety_atoms": ["O1"],
            "remove_site_hydrogens": ["HD21"],
            "remove_moiety_protons": [],
            "add_bonds": [{"site_atom": "ND2", "moiety_atom": "C1", "order": 1}],
        },
        "charge_patch_hint": {
            "strategy": "defer_to_ccd_residue_templates",
            "rationale": "CCD-backed glycan residue templates should define most junction charges",
        },
        "rationale": "Captures an ASN-linked glycan root as a declarative future graph edit",
        "notes": [
            "This placeholder does not validate sequon context yet",
            "CCD/Pablo residue naming will be reconciled in a later implementation phase",
        ],
    },
    {
        "identifier": "residue_replacement",
        "display_name": "Residue replacement placeholder",
        "allowed_sites": [
            {
                "residue_name": "ANY",
                "atom_name": "CA",
                "role": "replacement_anchor",
                "rationale": "Replacement-style PTMs are anchored by residue identity, not one bond edit",
            }
        ],
        "moiety_reactive_group": {
            "group": "replacement_residue",
            "anchor_atom": "CA",
            "leaving_atoms": [],
            "rationale": "The replacement residue supplies the complete post-translational state",
        },
        "graph_edits": {
            "delete_site_atoms": [],
            "delete_moiety_atoms": [],
            "remove_site_hydrogens": [],
            "remove_moiety_protons": [],
            "add_bonds": [],
        },
        "charge_patch_hint": {
            "strategy": "template_replacement",
            "rationale": "Residue replacement should use a complete residue template charge model",
        },
        "rationale": "Represents PTM-style residue replacement as a declarative placeholder",
        "notes": [
            "No executable atom deletion or bond creation is defined in this phase",
            "Future code should validate backbone compatibility before replacement",
        ],
    },
)


def load_builtin_mechanisms() -> dict[str, ReactionMechanism]:
    """Load all built-in mechanism declarations.

    Returns
    -------
    dict[str, ReactionMechanism]
        Built-in mechanisms keyed by normalized identifier.
    """
    mechanisms = [
        ReactionMechanism.model_validate(definition) for definition in BUILTIN_MECHANISM_DEFINITIONS
    ]
    return {mechanism.identifier: mechanism for mechanism in mechanisms}


def list_builtin_mechanisms() -> tuple[str, ...]:
    """List available built-in mechanism identifiers.

    Returns
    -------
    tuple[str, ...]
        Sorted mechanism identifiers.
    """
    return tuple(sorted(load_builtin_mechanisms()))


def get_builtin_mechanism(name: str) -> ReactionMechanism:
    """Return one built-in mechanism by identifier.

    Parameters
    ----------
    name : str
        Mechanism identifier from a conjugation attachment config.

    Returns
    -------
    ReactionMechanism
        Matching built-in mechanism declaration.

    Raises
    ------
    KeyError
        If the mechanism identifier is unknown.
    """
    normalized = name.strip().lower()
    mechanisms = load_builtin_mechanisms()
    if normalized not in mechanisms:
        available = ", ".join(sorted(mechanisms))
        raise KeyError(f"Unknown conjugation mechanism '{name}'. Available mechanisms: {available}")
    return mechanisms[normalized]
