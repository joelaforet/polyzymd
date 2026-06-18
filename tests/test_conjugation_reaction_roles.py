"""Tests for generic covalent reaction role exploration helpers."""

from __future__ import annotations

import inspect
from types import SimpleNamespace

import pytest

from polyzymd.builders.conjugation import reaction_roles
from polyzymd.builders.conjugation.pablo_reaction import (
    PabloReactionRequest,
    explore_pablo_residue_reaction,
)
from polyzymd.builders.conjugation.reaction_roles import (
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


def test_pablo_reaction_helper_passes_two_product_residue_names_and_counts_products():
    """Pablo exploration should pass through two product names and summarize products."""
    calls = []

    class FakeResidueDefinition:
        @staticmethod
        def react(**kwargs):
            calls.append(kwargs)
            assert kwargs["product_residue_names"] == ("LYX", "NHX")
            lys_product = SimpleNamespace(
                residue_name="LYX",
                atoms=(SimpleNamespace(name="N1", leaving=False),),
                linking_bond=None,
                crosslink="fake-crosslink",
            )
            nhs_product = SimpleNamespace(
                residue_name="NHX",
                atoms=(
                    SimpleNamespace(name="C1", leaving=False),
                    SimpleNamespace(name="O2", leaving=True),
                ),
                linking_bond="fake-linking-bond",
                crosslink=None,
            )
            return [(lys_product, nhs_product)]

    fake_pablo = SimpleNamespace(ResidueDefinition=FakeResidueDefinition)
    request = PabloReactionRequest(
        reactant_smarts=("[N:1][H:2]", "[C:3](=[O:4])[O:5]"),
        product_smarts=("[N:1][C:3](=[O:4])", "[O:5][H:2]"),
        product_residue_names=("LYX", "NHX"),
    )

    diagnostic = explore_pablo_residue_reaction(
        reactants=(object(), object()),
        request=request,
        pablo_module=fake_pablo,
    )

    assert calls
    assert diagnostic.success is True
    assert diagnostic.product_set_count == 1
    assert [product.residue_name for product in diagnostic.products[0]] == ["LYX", "NHX"]
    assert diagnostic.products[0][1].leaving_atom_names == ("O2",)


def test_pablo_reaction_helper_skips_gracefully_when_pablo_is_missing(monkeypatch):
    """Missing Pablo should return a diagnostic instead of raising by default."""
    import polyzymd.builders.conjugation.pablo_reaction as pablo_reaction

    def fail_import(name: str):
        assert name == "openff.pablo"
        raise ImportError("no pablo")

    monkeypatch.setattr(pablo_reaction.importlib, "import_module", fail_import)
    request = PabloReactionRequest(
        reactant_smarts=("[N:1]",),
        product_smarts=("[N:1]",),
        product_residue_names=("LYX",),
    )

    diagnostic = explore_pablo_residue_reaction(reactants=(), request=request)

    assert diagnostic.success is False
    assert diagnostic.pablo_available is False
    assert diagnostic.error_type == "ImportError"
    assert diagnostic.warnings
