"""Tests for generic covalent reaction role exploration helpers."""

from __future__ import annotations

import inspect
from types import SimpleNamespace

import pytest

from polyzymd.builders.conjugation.reactions import _roles as reaction_roles
from polyzymd.builders.conjugation.reactions._roles import (
    AtomMappedReaction,
    AtomRoleSpec,
    PdbAtomIdentity,
    ReactionParticipant,
    atom_mapped_reaction_from_mechanism_config,
    derive_bond_changes,
    resolve_reaction_roles_from_identity_map,
    validate_mapped_smarts,
)


def _nhs_like_reaction() -> AtomMappedReaction:
    """Build a small NHS-Lys-like mapped reaction without PDB atom-name assumptions."""
    return AtomMappedReaction(
        name="generic_nhs_amine_poc",
        reactant_smarts=("[N:1]([H:2])", "[C:3](=[O:4])[O:5][N:6]"),
        product_smarts=("[N:1][C:3](=[O:4])",),
        participants=(
            ReactionParticipant(name="target amine", role="site", reactant_index=0),
            ReactionParticipant(name="activated ester", role="moiety", reactant_index=1),
        ),
        atom_roles=(
            AtomRoleSpec(map_number=1, participant="site", role="linking", label="site_link"),
            AtomRoleSpec(map_number=2, participant="site", role="leaving", label="site_proton"),
            AtomRoleSpec(map_number=3, participant="moiety", role="linking", label="moiety_link"),
            AtomRoleSpec(
                map_number=4,
                participant="moiety",
                role="geometry_anchor",
                label="carbonyl_anchor",
            ),
            AtomRoleSpec(
                map_number=5, participant="moiety", role="leaving", label="leaving_oxygen"
            ),
        ),
    )


def test_nhs_like_mapped_smarts_derives_added_removed_bonds_and_leaving_maps():
    """Mapped SMARTS should expose bond changes independent of atom names."""
    reaction = _nhs_like_reaction()

    validation = validate_mapped_smarts(reaction.reactant_smarts, reaction.product_smarts)
    changes = derive_bond_changes(reaction.reactant_smarts, reaction.product_smarts)

    assert validation.missing_product_maps == (2, 5, 6)
    assert [change.atom_maps for change in changes.added_bonds] == [(1, 3)]
    assert (3, 5) in {change.atom_maps for change in changes.removed_bonds}
    assert changes.order_changes == ()


def test_rdkit_reaction_parse_detects_missing_product_atom_maps():
    """Native reaction parsing should retain missing-map diagnostics."""
    pytest.importorskip("rdkit.Chem.rdChemReactions")

    validation = validate_mapped_smarts(("[N:1]([H:2])",), ("[N:1]",))

    assert validation.reactant_maps == (1, 2)
    assert validation.product_maps == (1,)
    assert validation.missing_product_maps == (2,)


def test_bond_change_detection_matches_fallback_for_added_removed_and_order_changed(
    monkeypatch,
):
    """RDKit-backed bond diagnostics should match the existing fallback behavior."""
    reactants = ("[C:1]=[O:2].[N:3][H:4]",)
    products = ("[O:2]-[C:1][N:3]",)
    normal = derive_bond_changes(reactants, products)

    def raise_value_error(smarts_entries):
        raise ValueError("force fallback")

    monkeypatch.setattr(reaction_roles, "_rdkit_mapped_bond_orders", raise_value_error)
    fallback = derive_bond_changes(reactants, products)

    assert normal.model_dump() == fallback.model_dump()
    assert [change.atom_maps for change in normal.added_bonds] == [(1, 3)]
    assert [change.atom_maps for change in normal.removed_bonds] == [(3, 4)]
    assert [
        (change.atom_maps, change.reactant_order, change.product_order)
        for change in normal.order_changes
    ] == [((1, 2), 2.0, 1.0)]


def test_rdkit_reaction_parse_preserves_participant_indices_and_map_evidence_order():
    """RDKit parsing should not mutate participant order or sorted map evidence."""
    pytest.importorskip("rdkit.Chem.rdChemReactions")
    reaction = AtomMappedReaction(
        name="participant_order",
        reactant_smarts=("[C:10]", "[N:1].[O:5]"),
        product_smarts=("[N:1][C:10]=[O:5]",),
        participants=(
            ReactionParticipant(name="moiety", role="moiety", reactant_index=1),
            ReactionParticipant(name="site", role="site", reactant_index=0),
        ),
    )

    validation = validate_mapped_smarts(reaction.reactant_smarts, reaction.product_smarts)

    assert [participant.reactant_index for participant in reaction.participants] == [1, 0]
    assert validation.reactant_maps == (1, 5, 10)
    assert validation.product_maps == (1, 5, 10)


def test_mapped_smarts_validation_fallback_still_handles_expected_rdkit_parse_errors(
    monkeypatch,
):
    """Expected RDKit reaction parse errors should fall back to text diagnostics."""

    def raise_value_error(reactant_smarts, product_smarts):
        raise ValueError("invalid reaction query")

    monkeypatch.setattr(reaction_roles, "_rdkit_reaction_mapped_atom_numbers", raise_value_error)

    validation = validate_mapped_smarts(("[N:1][H:2]",), ("[N:1]",))

    assert validation.reactant_maps == (1, 2)
    assert validation.product_maps == (1,)
    assert validation.missing_product_maps == (2,)


def test_mapped_smarts_validation_does_not_swallow_unexpected_runtime_errors(monkeypatch):
    """Unexpected reaction parsing runtime errors should propagate."""

    def raise_runtime_error(reactant_smarts, product_smarts):
        raise RuntimeError("unexpected reaction parser failure")

    monkeypatch.setattr(reaction_roles, "_rdkit_reaction_mapped_atom_numbers", raise_runtime_error)

    with pytest.raises(RuntimeError, match="unexpected reaction parser failure"):
        validate_mapped_smarts(("[N:1][H:2]",), ("[N:1]",))


def test_resolve_reaction_roles_preserves_geometry_anchor_and_concrete_identities():
    """Role resolution should keep role labels and supplied PDB identities generic."""
    reaction = _nhs_like_reaction()
    identities = {
        1: PdbAtomIdentity(chain_id="A", residue_name="AAA", residue_number=10, atom_name="N1"),
        2: PdbAtomIdentity(chain_id="A", residue_name="AAA", residue_number=10, atom_name="H1"),
        3: PdbAtomIdentity(chain_id="B", residue_name="BBB", residue_number=1, atom_name="C1"),
        4: PdbAtomIdentity(chain_id="B", residue_name="BBB", residue_number=1, atom_name="O1"),
        5: PdbAtomIdentity(chain_id="B", residue_name="BBB", residue_number=1, atom_name="O2"),
    }

    plan = resolve_reaction_roles_from_identity_map(reaction, identities)

    assert plan.inferred_leaving_maps == (2, 5, 6)
    assert plan.resolved_added_bonds[0].change.atom_maps == (1, 3)
    anchor_roles = [role for role in plan.resolved_roles if role.spec.role == "geometry_anchor"]
    assert anchor_roles[0].spec.map_number == 4
    assert anchor_roles[0].pdb_identity == identities[4]
    dumped = plan.model_dump_json()
    assert "LYX" not in dumped
    assert "NHX" not in dumped
    assert "C047" not in dumped
    assert "O020" not in dumped


def test_reaction_role_framework_does_not_hardcode_poc_residue_or_atom_names():
    """Framework logic should stay generic; POC names belong in notebooks/diagnostics."""
    source = inspect.getsource(reaction_roles)

    for forbidden in ("LYX", "NHX", "NHSX", "C047", "O020", "H11", "H13", "HZ2", "HZ3"):
        assert forbidden not in source


def test_mechanism_config_adapter_derives_role_and_bond_artifacts_without_atom_names():
    """Config-style reaction SMARTS should preflight without PDB atom-name assumptions."""
    mechanism = SimpleNamespace(
        name="generic_amide",
        reaction_smarts="[N:1]([H:2]).[C:3](=[O:4])[O:5]>>[N:1][C:3](=[O:4])",
        atom_roles=[
            {"map_number": 1, "participant": "site", "role": "linking", "label": "site"},
            {"map_number": 2, "participant": "site", "role": "leaving"},
            {"map_number": 3, "participant": "moiety", "role": "linking"},
            {"map_number": 4, "participant": "moiety", "role": "geometry_anchor"},
            {"map_number": 5, "participant": "moiety", "role": "leaving"},
        ],
    )

    reaction = atom_mapped_reaction_from_mechanism_config(mechanism)
    plan = resolve_reaction_roles_from_identity_map(
        reaction,
        {},
        require_required_identities=False,
    )

    assert [role.map_number for role in reaction.atom_roles] == [1, 2, 3, 4, 5]
    assert [change.atom_maps for change in plan.bond_changes.added_bonds] == [(1, 3)]
    assert plan.inferred_leaving_maps == (2, 5)
    assert all(role.pdb_identity is None for role in plan.resolved_roles)


def test_identity_map_resolver_boundary_requires_supplied_required_identities():
    """The generic resolver should not pretend to match structures by itself."""
    reaction = _nhs_like_reaction()

    with pytest.raises(ValueError, match="No PDB identity was provided"):
        resolve_reaction_roles_from_identity_map(reaction, {})


def test_pablo_strict_validation_can_require_balanced_product_maps():
    """The Pablo path can opt into all-reactant-maps-present validation."""
    reaction = _nhs_like_reaction()

    with pytest.raises(ValueError, match="missing reactant atom-map numbers"):
        validate_mapped_smarts(
            reaction.reactant_smarts,
            reaction.product_smarts,
            require_all_reactant_maps_in_products=True,
        )


def test_mapped_bond_order_fallback_does_not_swallow_programming_errors(monkeypatch):
    """Unexpected RDKit helper errors should not silently use the fallback parser."""

    def raise_programming_error(smarts_entries):
        raise TypeError("bad helper call")

    monkeypatch.setattr(reaction_roles, "_rdkit_mapped_bond_orders", raise_programming_error)

    with pytest.raises(TypeError, match="bad helper call"):
        reaction_roles._mapped_bond_orders(("[N:1][C:2]",))


def test_mapped_bond_order_fallback_does_not_swallow_runtime_errors(monkeypatch):
    """Unexpected runtime errors should propagate instead of using fallback parsing."""

    def raise_runtime_error(smarts_entries):
        raise RuntimeError("unexpected helper failure")

    monkeypatch.setattr(reaction_roles, "_rdkit_mapped_bond_orders", raise_runtime_error)

    with pytest.raises(RuntimeError, match="unexpected helper failure"):
        reaction_roles._mapped_bond_orders(("[N:1][C:2]",))


def test_mapped_bond_order_fallback_still_handles_value_errors(monkeypatch):
    """Expected RDKit parsing errors should still fall back to the lightweight parser."""

    def raise_value_error(smarts_entries):
        raise ValueError("invalid query")

    monkeypatch.setattr(reaction_roles, "_rdkit_mapped_bond_orders", raise_value_error)

    assert reaction_roles._mapped_bond_orders(("[N:1][C:2]",)) == {(1, 2): 1.0}
