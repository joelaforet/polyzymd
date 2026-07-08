"""Tests for config-driven conjugation system workflow helpers."""

from __future__ import annotations

import json
import logging
from pathlib import Path
from types import SimpleNamespace

import pytest

from polyzymd.builders.conjugation._linkage import (
    ExplicitLinkageContract,
    LinkageBond,
    PabloCrosslinkRequirement,
    PdbAtomSelector,
    ReactiveEndpoint,
    ResolvedAttachmentPlan,
)
from polyzymd.builders.conjugation.construction import ModifierConstructionSettings
from polyzymd.builders.conjugation.pablo.ingestion import PabloAvailability, PabloIngestionResult
from polyzymd.builders.conjugation.pablo.parameterization import InterchangeParameterizationResult
from polyzymd.builders.conjugation.placement import PackmolModifierPlacementResult
from polyzymd.builders.conjugation.polymer import (
    GeneratedMoietyFragment,
    GeneratedPolymerFragment,
    MultiResidueGenerationResult,
    PolymerFragmentAtom,
    PolymerFragmentResidue,
    PolymerMonomerRecipe,
    PolymerRecipe,
)
from polyzymd.builders.conjugation.structure.pdb import (
    CrosslinkedPdbAssemblyResult,
    PdbAtomRecord,
    PlacedPolymerFragment,
)
from polyzymd.builders.conjugation.structure.preparation import (
    ProteinCanonicalizationResult,
    ProteinCanonicalizationSettings,
)
from polyzymd.builders.conjugation.system_workflow import (
    _apply_pdb_atom_names_to_topology,
    _build_direct_solvated_system,
    _build_solvated_system,
    _local_minimization_settings_for_product,
    _policy_with_resolved_crosslink,
    _prepared_protein_pdb_path,
    _require_supported_coordinate_backend,
    _restore_pdb_atom_name_fields,
)
from polyzymd.builders.conjugation.validation import ValidationStatus
from polyzymd.config.schema import (
    ConjugationCcdPabloPolicyConfig,
    ConjugationChainPolicyConfig,
)


class _AtomDouble:
    def __init__(self, name: str) -> None:
        self.name = name
        self.metadata: dict[str, str] = {}


class _TopologyDouble:
    def __init__(self, atoms: list[_AtomDouble]) -> None:
        self._atoms = atoms

    def atoms(self) -> tuple[_AtomDouble, ...]:
        return tuple(self._atoms)


def test_restore_pdb_atom_name_fields_updates_only_template_atom_count(tmp_path: Path):
    """Conjugate template names should be restored without touching solvent names."""
    template = tmp_path / "linked.pdb"
    target = tmp_path / "solvated.pdb"
    template.write_text(
        "".join(
            [
                _pdb_line(1, " N  ", "ALA", "A", 1),
                _pdb_line(2, " CA ", "ALA", "A", 1),
            ]
        ),
        encoding="utf-8",
    )
    target.write_text(
        "".join(
            [
                _pdb_line(1, " N1x", "ALA", "A", 1),
                _pdb_line(2, " C1x", "ALA", "A", 1),
                _pdb_line(3, " O  ", "HOH", "D", 1),
            ]
        ),
        encoding="utf-8",
    )

    restored = _restore_pdb_atom_name_fields(target, template)
    lines = target.read_text(encoding="utf-8").splitlines()

    assert restored == 2
    assert lines[0][12:16] == " N  "
    assert lines[1][12:16] == " CA "
    assert lines[2][12:16] == " O  "


def test_apply_pdb_atom_names_to_topology_uses_same_order_template(tmp_path: Path):
    """OpenFF topology atom names should be reset from the linked PDB template."""
    template = tmp_path / "linked.pdb"
    template.write_text(
        "".join(
            [
                _pdb_line(1, " N  ", "ALA", "A", 1),
                _pdb_line(2, " CA ", "ALA", "A", 1),
            ]
        ),
        encoding="utf-8",
    )
    atoms = [_AtomDouble("N1x"), _AtomDouble("C1x")]

    _apply_pdb_atom_names_to_topology(_TopologyDouble(atoms), template)

    assert [atom.name for atom in atoms] == ["N", "CA"]
    assert [atom.metadata["atom_name"] for atom in atoms] == ["N", "CA"]


def test_policy_with_resolved_crosslink_uses_product_state_leaving_atoms():
    """Generated Pablo policies should not remove atoms from product-state PDBs."""
    requirement = PabloCrosslinkRequirement(
        residues=("LYX", "NHX"),
        linking_atoms=("NZ", "C047"),
        leaving_atoms=(("H11", "H13"), ("O020",)),
        bond_order=1,
    )
    policy = ConjugationCcdPabloPolicyConfig(crosslinks=[])

    generated = _policy_with_resolved_crosslink(
        policy,
        SimpleNamespace(pablo_crosslink_requirement=requirement),
    )

    assert generated.crosslinks[0].leaving_atoms == ((), ())


def test_system_workflow_settings_enable_public_product_state_defaults():
    """Public workflow defaults should use the generic conjugate protocol."""
    import polyzymd.builders.conjugation.system_workflow as workflow_module

    settings = workflow_module.ConjugatedPolymerSystemSettings()

    assert settings.canonicalize_source_protein_hydrogens is True
    assert settings.use_product_state_pablo_library is True
    assert settings.run_relaxation is True
    assert settings.run_product_state_local_minimization is False
    assert settings.protein_canonicalization.ph == pytest.approx(7.0)


def test_relaxed_conjugate_pdb_prefers_final_conjugate_relaxation_artifact(tmp_path):
    """Downstream solvation should consume the Stage B relaxed PDB."""
    import polyzymd.builders.conjugation.system_workflow as workflow_module

    minimized = tmp_path / "conjugate_minimized.pdb"
    relaxed = tmp_path / "conjugate_relaxed.pdb"
    minimized.write_text("END\n", encoding="utf-8")
    relaxed.write_text("END\n", encoding="utf-8")
    construction = SimpleNamespace(
        local_minimization=None,
        relaxation=SimpleNamespace(
            minimized_pdb_path=minimized,
            relaxed_pdb_path=relaxed,
        ),
    )

    assert workflow_module._relaxed_conjugate_pdb(construction) == relaxed


def test_build_solvated_system_uses_conjugation_final_interchange(monkeypatch, tmp_path: Path):
    """Config-driven final solvation should use the conjugation-specific helper."""
    import polyzymd.builders.conjugation.system_workflow as workflow_module

    config = SimpleNamespace(
        substrate=None,
        polymers=SimpleNamespace(enabled=True),
        solvent=SimpleNamespace(),
    )
    fake_builder = _NoKwargInterchangeBuilder(tmp_path / "solvated.pdb")
    product_library = _fake_product_state_library()
    calls = {}

    monkeypatch.setattr(
        workflow_module.SystemBuilder,
        "from_config",
        classmethod(lambda cls, config: fake_builder),
    )
    monkeypatch.setattr(workflow_module, "_build_and_pack_free_polymers", lambda *a, **k: None)

    def fake_final_interchange(builder, *, product_state_pablo_library, settings=None):
        calls["builder"] = builder
        calls["product_state_pablo_library"] = product_state_pablo_library
        calls["settings"] = settings
        builder._interchange = object()
        return builder._interchange

    monkeypatch.setattr(
        workflow_module,
        "create_final_conjugated_interchange",
        fake_final_interchange,
    )

    result = _build_solvated_system(
        config,
        relaxed_conjugate_topology=object(),
        working_dir=tmp_path,
        polymer_seed=123,
        create_interchange=True,
        product_state_pablo_library=product_library,
    )

    assert result is fake_builder
    assert fake_builder.create_interchange_calls == 0
    assert calls["builder"] is fake_builder
    assert calls["product_state_pablo_library"] is product_library


def test_build_direct_solvated_system_uses_conjugation_final_interchange(
    monkeypatch, tmp_path: Path
):
    """Direct final solvation should use the conjugation-specific helper."""
    import polyzymd.builders.conjugation.system_workflow as workflow_module

    fake_builder = _NoKwargInterchangeBuilder(tmp_path / "solvated.pdb")
    product_library = _fake_product_state_library()
    calls = {}
    monkeypatch.setattr(workflow_module, "SystemBuilder", lambda: fake_builder)
    monkeypatch.setattr(
        workflow_module,
        "create_final_conjugated_interchange",
        lambda builder, *, product_state_pablo_library, settings=None: calls.setdefault(
            "args", (builder, product_state_pablo_library, settings)
        ),
    )

    result = _build_direct_solvated_system(
        relaxed_conjugate_topology=object(),
        working_dir=tmp_path,
        create_interchange=True,
        product_state_pablo_library=product_library,
        padding=1.7,
        box_shape="truncated_octahedron",
    )

    assert result is fake_builder
    assert fake_builder.create_interchange_calls == 0
    assert fake_builder.solvation_settings == {
        "padding": 1.7,
        "box_shape": "truncated_octahedron",
    }
    assert calls["args"][0] is fake_builder
    assert calls["args"][1] is product_library


def test_prepared_protein_path_canonicalizes_to_construction_dir(monkeypatch, tmp_path: Path):
    """Config workflow should use a canonicalized source-protein PDB by default."""
    import polyzymd.builders.conjugation.system_workflow as workflow_module

    source = tmp_path / "protein.pdb"
    source.write_text("END\n", encoding="utf-8")
    calls = {}

    def fake_canonicalize(pdb_path, output_path=None, *, settings=None):
        calls["pdb_path"] = Path(pdb_path)
        calls["output_path"] = Path(output_path)
        calls["settings"] = settings
        Path(output_path).write_text("END\n", encoding="utf-8")
        return ProteinCanonicalizationResult(
            input_path=Path(pdb_path),
            output_path=Path(output_path),
            ph=settings.ph,
            force_field_name=settings.force_field_name,
            atom_count=1,
            hydrogen_count=0,
        )

    monkeypatch.setattr(workflow_module, "canonicalize_protein_hydrogens", fake_canonicalize)
    settings = workflow_module.ConjugatedPolymerSystemSettings(
        protein_canonicalization=ProteinCanonicalizationSettings(ph=6.5)
    )

    prepared_path, result = _prepared_protein_pdb_path(
        source,
        output_dir=tmp_path / "construction",
        settings=settings,
    )

    assert prepared_path == tmp_path / "construction" / "source_protein_canonical_pH6.5.pdb"
    assert result is not None
    assert calls["pdb_path"] == source
    assert calls["settings"].ph == pytest.approx(6.5)
    assert prepared_path.with_suffix(".canonicalization.json").exists()


def test_local_minimization_settings_use_product_atom_identities(tmp_path: Path):
    """Local minimization selectors should use generic product link atoms."""
    from polyzymd.builders.conjugation.local_minimization import LocalMinimizationSettings

    product = tmp_path / "product.pdb"
    product.write_text(
        "".join(
            [
                _pdb_line(10, " NZ ", "LYX", "A", 23, element="N"),
                _pdb_line(20, "C047", "NHX", "C", 5, element="C"),
                _pdb_line(21, "OX1 ", "NHX", "C", 5, element="O"),
                _pdb_line(22, "OY2 ", "NHX", "C", 5, element="O"),
                "CONECT   10   20\n",
                "CONECT   20   21\n",
            ]
        ),
        encoding="utf-8",
    )
    requirement = PabloCrosslinkRequirement(
        residues=("LYX", "NHX"),
        linking_atoms=("NZ", "C047"),
        leaving_atoms=(("H11",), ("O021",)),
        bond_order=1,
    )

    settings = _local_minimization_settings_for_product(
        product,
        base_settings=LocalMinimizationSettings(),
        requirements=(requirement,),
    )

    linkage = settings.linkages[0]
    assert linkage.protein_link_selector.serial is None
    assert linkage.protein_link_selector.chain_id == "A"
    assert linkage.protein_link_selector.residue_name == "LYX"
    assert linkage.protein_link_selector.residue_number == 23
    assert linkage.modifier_link_selector.chain_id == "C"
    assert linkage.modifier_link_selector.residue_number == 5
    assert linkage.modifier_link_selector.atom_name == "C047"


def test_local_minimization_settings_preserve_explicit_selectors(tmp_path: Path):
    """Explicit local minimization selectors should override product-state inference."""
    from polyzymd.builders.conjugation.local_minimization import (
        LocalLinkageAtomSelector,
        LocalLinkageSelectors,
        LocalMinimizationSettings,
    )

    settings = LocalMinimizationSettings(
        linkages=(
            LocalLinkageSelectors(
                protein_link_selector=LocalLinkageAtomSelector(
                    chain_id="A",
                    residue_name="LYX",
                    residue_number=23,
                    atom_name="NZ",
                ),
                modifier_link_selector=LocalLinkageAtomSelector(
                    chain_id="C",
                    residue_name="NHX",
                    residue_number=5,
                    atom_name="USER",
                ),
                target_bond_length_angstrom=1.33,
            ),
        ),
    )
    requirement = PabloCrosslinkRequirement(
        residues=("LYX", "NHX"),
        linking_atoms=("NZ", "C047"),
        leaving_atoms=(("H11",), ("O021",)),
        bond_order=1,
    )

    derived = _local_minimization_settings_for_product(
        tmp_path / "missing-product.pdb",
        base_settings=settings,
        requirements=(requirement,),
    )

    assert derived is settings
    assert derived.linkages[0].modifier_link_selector.atom_name == "USER"


def test_local_minimization_settings_do_not_require_modifier_oxygen_anchor(tmp_path: Path):
    """Generic local minimization should not require an auxiliary oxygen anchor."""
    from polyzymd.builders.conjugation.local_minimization import LocalMinimizationSettings

    product = tmp_path / "product.pdb"
    product.write_text(
        "".join(
            [
                _pdb_line(10, " NZ ", "LYX", "A", 23, element="N"),
                _pdb_line(20, "C047", "NHX", "C", 5, element="C"),
                _pdb_line(21, "ODEF", "NHX", "C", 5, element="O"),
            ]
        ),
        encoding="utf-8",
    )
    requirement = PabloCrosslinkRequirement(
        residues=("LYX", "NHX"),
        linking_atoms=("NZ", "C047"),
        leaving_atoms=(("H11",), ("O021",)),
        bond_order=1,
    )

    settings = _local_minimization_settings_for_product(
        product,
        base_settings=LocalMinimizationSettings(),
        requirements=(requirement,),
        product_state_pablo_library=SimpleNamespace(definitions=()),
    )

    assert settings.linkages[0].modifier_link_selector.atom_name == "C047"


def test_nhs_lys_shim_uses_product_state_local_minimization(
    monkeypatch,
    tmp_path: Path,
):
    """The legacy NHS-Lys shim should preserve single-site local minimization."""
    import polyzymd.builders.conjugation.pablo.product as product_pablo_module
    import polyzymd.builders.conjugation.system_workflow as workflow_module

    requirement = PabloCrosslinkRequirement(
        residues=("LYX", "NHX"),
        linking_atoms=("NZ", "C047"),
        leaving_atoms=(("H11",), ("O020",)),
        bond_order=1,
    )
    resolved_plan = _resolved_plan(requirement)
    product_library = SimpleNamespace(residue_library=object())
    production_template = SimpleNamespace(partial_charges=(0.0,))
    calls = {}

    placed = PlacedPolymerFragment(
        atoms=(
            PdbAtomRecord(
                serial=2,
                atom_name="C047",
                residue_name="NHX",
                chain_id="C",
                residue_number=5,
                x=1,
                y=0,
                z=0,
            ),
        ),
        reactive_atom_name="C047",
    )
    placement = PackmolModifierPlacementResult(
        placed_modifier=placed,
        packmol_input_path=tmp_path / "packmol.inp",
        packmol_output_path=tmp_path / "packmol.pdb",
        protein_sterics_pdb_path=tmp_path / "protein.pdb",
        modifier_pdb_path=tmp_path / "modifier.pdb",
        packmol_input_text="",
        target_bond_length_angstrom=1.33,
        placed_bond_length_angstrom=1.33,
        min_modifier_protein_distance_angstrom=2.0,
    )

    def fake_place(protein_pdb_path, modifiers, plans, output_dir, *, settings=None):
        calls["placed_protein_path"] = Path(protein_pdb_path)
        calls["modifier_count"] = len(modifiers)
        calls["plan_count"] = len(plans)
        return (placement,)

    def fake_write(protein_pdb_path, polymer_fragment, attachment, output_path, options):
        calls["assembly_protein_path"] = Path(protein_pdb_path)
        Path(output_path).write_text(
            "".join(
                [
                    _pdb_line(10, " NZ ", "LYX", "A", 23, element="N"),
                    _pdb_line(20, "C047", "NHX", "C", 5, element="C"),
                    _pdb_line(21, "OX1 ", "NHX", "C", 5, element="O"),
                    "CONECT   10   20\n",
                    "CONECT   20   21\n",
                ]
            ),
            encoding="utf-8",
        )
        return CrosslinkedPdbAssemblyResult(
            output_path=Path(output_path),
            protein_atom_count=1,
            polymer_atom_count=2,
            removed_protein_atom_count=1,
            removed_polymer_atom_count=1,
            removed_atom_serials=(11, 21),
            removed_atom_names=("H11", "O020"),
            added_conect_pair=(10, 20),
        )

    def fake_product_library(**kwargs):
        calls["product_pdb"] = Path(kwargs["product_pdb"])
        calls["source_protein_pdb"] = Path(kwargs["source_protein_pdb"])
        calls["polymer_sdf"] = Path(kwargs["polymer_sdf"])
        return product_library

    class FakeIngestor:
        def __init__(self, policy):
            self.policy = policy

        def ingest_structure(
            self, path, *, chain_policy=None, output_dir=None, residue_library=None
        ):
            calls["residue_library"] = residue_library
            return PabloIngestionResult(
                success=True,
                path=Path(path),
                suffix=".pdb",
                pablo=PabloAvailability(available=True),
                topology=SimpleNamespace(molecules=(object(),)),
            )

    def fake_parameterize(
        topology,
        *,
        settings=None,
        charge_from_molecules=None,
        require_charge_templates=False,
    ):
        calls["charge_template_count"] = len(charge_from_molecules or ())
        calls["require_charge_templates"] = require_charge_templates
        return InterchangeParameterizationResult(
            success=True,
            interchange=object(),
            force_field_names=("fake.offxml",),
            topology_type="FakeTopology",
        )

    def fake_validation(
        interchange,
        output_dir,
        *,
        settings=None,
        crosslinked_pdb_path=None,
        attachment_specs=(),
    ):
        calls["validation"] = calls.get("validation", 0) + 1
        calls["validation_crosslinked_pdb_path"] = crosslinked_pdb_path
        calls["validation_attachment_specs"] = attachment_specs
        return SimpleNamespace(success=True)

    def fake_local_minimization(pdb_path, output_dir, **kwargs):
        calls["local_minimization"] = calls.get("local_minimization", 0) + 1
        calls["local_settings"] = kwargs["settings"]
        calls["local_policy"] = kwargs["pablo_policy"]
        calls["local_product_library"] = kwargs["product_state_pablo_library"]
        relaxed = Path(output_dir) / "local-relaxed.pdb"
        relaxed.write_text("END\n", encoding="utf-8")
        return SimpleNamespace(
            success=False,
            relaxed_pdb_path=relaxed,
            blocker=None,
        )

    monkeypatch.setattr(workflow_module, "place_modifiers_with_resolved_plans", fake_place)
    monkeypatch.setattr(
        workflow_module, "placed_fragment_from_resolved_plan", lambda fragment, plan: fragment
    )
    monkeypatch.setattr(workflow_module, "write_crosslinked_pdb", fake_write)
    monkeypatch.setattr(
        product_pablo_module, "build_product_state_pablo_library", fake_product_library
    )
    monkeypatch.setattr(
        workflow_module,
        "build_conjugate_charge_templates",
        lambda topology, library: (production_template,),
    )
    monkeypatch.setattr(workflow_module, "PabloIngestor", FakeIngestor)
    monkeypatch.setattr(
        workflow_module, "create_interchange_from_pablo_topology", fake_parameterize
    )
    monkeypatch.setattr(
        workflow_module,
        "run_post_crosslink_local_minimization",
        fake_local_minimization,
    )

    policy = ConjugationCcdPabloPolicyConfig(
        crosslinks=[
            {
                "residues": ("LYX", "NHX"),
                "linking_atoms": ("NZ", "C047"),
                "leaving_atoms": ((), ()),
                "bond_order": 1,
            }
        ]
    )
    prepared_protein = tmp_path / "canonical.pdb"
    prepared_protein.write_text("END\n", encoding="utf-8")
    polymer_sdf = tmp_path / "polymer.sdf"
    polymer_sdf.write_text("", encoding="utf-8")

    spec = SimpleNamespace(
        attachment_id="nhs_lys_attachment_01",
        attachment_index=1,
        reaction_name="nhs_lys",
        generated_fragment=object(),
        resolved_plan=resolved_plan,
        source_sidecars={"sdf": polymer_sdf},
        fragment=SimpleNamespace(source_kind="polymer"),
    )

    construction, topology = workflow_module._construct_conjugate_from_specs(
        protein_pdb_path=prepared_protein,
        specs=(spec,),
        ccd_pablo_policy=policy,
        chain_policy=None,
        output_dir=tmp_path / "construction",
        settings=ModifierConstructionSettings(run_relaxation=False),
        use_product_state_pablo_library=True,
        use_conjugate_relaxation=False,
        run_product_state_local_minimization=True,
        local_minimization_settings=workflow_module._default_local_minimization_settings(),
    )

    assert topology is not None
    assert construction.local_minimization is not None
    assert construction.local_minimization.relaxed_pdb_path.name == "local-relaxed.pdb"
    assert construction.local_minimization.success is False
    assert construction.relaxation is None
    assert calls["placed_protein_path"] == prepared_protein
    assert calls["modifier_count"] == 1
    assert calls["plan_count"] == 1
    assert calls["assembly_protein_path"] == prepared_protein
    assert calls["source_protein_pdb"] == prepared_protein
    assert calls["polymer_sdf"] == polymer_sdf
    assert calls["residue_library"] is product_library.residue_library
    assert calls["charge_template_count"] == 1
    assert calls["require_charge_templates"] is True
    assert calls.get("validation", 0) == 0
    assert calls["local_minimization"] == 1
    assert calls["local_settings"].linkages[0].protein_link_selector.residue_name == "LYX"
    assert calls["local_settings"].linkages[0].modifier_link_selector.residue_name == "NHX"
    assert calls["local_policy"].crosslinks[0].leaving_atoms == ((), ())
    assert calls["local_product_library"] is product_library
    assert any(
        "local minimization completed and wrote a relaxed PDB" in item
        for item in construction.diagnostics
    )


def test_multi_modifier_construction_places_parameterizes_and_relaxes_once(
    monkeypatch,
    tmp_path: Path,
):
    """Multi-site construction should jointly place and parameterize the product once."""
    import polyzymd.builders.conjugation.system_workflow as workflow_module

    plans = (
        _generic_resolved_plan(residue_number=42, modifier_residue_number=1, modifier_atom="C001"),
        _generic_resolved_plan(residue_number=87, modifier_residue_number=2, modifier_atom="C001"),
    )
    modifiers = tuple(
        _generated_fragment(residue_name="NAG", residue_number=index + 1) for index in range(2)
    )
    moieties = tuple(
        _moiety_fragment(residue_name="NAG", residue_number=index + 1) for index in range(2)
    )
    calls = {"placements": 0, "parameterize": 0, "validation": 0, "local_minimization": 0}

    def fake_place(protein_pdb_path, modifiers_arg, plans_arg, output_dir, *, settings=None):
        calls["placements"] += 1
        calls["placed_modifier_count"] = len(modifiers_arg)
        calls["placed_plan_count"] = len(plans_arg)
        return tuple(
            PackmolModifierPlacementResult(
                placed_modifier=modifier.to_placed_fragment(),
                packmol_input_path=Path(output_dir) / "packmol.inp",
                packmol_output_path=Path(output_dir) / "packmol.pdb",
                protein_sterics_pdb_path=Path(output_dir) / "protein.pdb",
                modifier_pdb_path=Path(output_dir) / f"modifier_{index}.pdb",
                packmol_input_text="",
                target_bond_length_angstrom=plan.target_bond_length_angstrom,
                placed_bond_length_angstrom=plan.target_bond_length_angstrom,
                min_modifier_protein_distance_angstrom=2.0,
            )
            for index, (modifier, plan) in enumerate(
                zip(modifiers_arg, plans_arg, strict=True),
                start=1,
            )
        )

    def fake_write(protein_pdb_path, polymer_fragment, attachment, output_path, options):
        calls["fragment_count"] = len(polymer_fragment)
        calls["attachment_count"] = len(attachment)
        calls["protein_products"] = tuple(item.protein_target_resname for item in attachment)
        Path(output_path).write_text(
            "".join(
                [
                    _pdb_line(1, "ND2 ", "ASX", "A", 42, element="N"),
                    _pdb_line(2, "C001", "NAG", "C", 1, element="C"),
                    _pdb_line(3, "ND2 ", "ASX", "A", 87, element="N"),
                    _pdb_line(4, "C001", "NAG", "C", 2, element="C"),
                    "CONECT    1    2\n",
                    "CONECT    2    1\n",
                    "CONECT    3    4\n",
                    "CONECT    4    3\n",
                ]
            ),
            encoding="utf-8",
        )
        return CrosslinkedPdbAssemblyResult(
            output_path=Path(output_path),
            protein_atom_count=10,
            polymer_atom_count=2,
            removed_protein_atom_count=2,
            removed_polymer_atom_count=0,
            removed_atom_serials=(),
            removed_atom_names=(),
            added_conect_pair=(1, 2),
            added_conect_pairs=((1, 2), (3, 4)),
        )

    class FakeIngestor:
        def __init__(self, policy):
            self.policy = policy

        def ingest_structure(
            self, path, *, chain_policy=None, output_dir=None, residue_library=None
        ):
            return PabloIngestionResult(
                success=True,
                path=Path(path),
                suffix=".pdb",
                pablo=PabloAvailability(available=True),
                topology=SimpleNamespace(molecules=(object(),)),
            )

    def fake_parameterize(
        topology,
        *,
        settings=None,
        charge_from_molecules=None,
        require_charge_templates=False,
    ):
        calls["parameterize"] += 1
        calls["require_charge_templates"] = require_charge_templates
        return InterchangeParameterizationResult(
            success=True,
            interchange=object(),
            force_field_names=("fake.offxml",),
            topology_type="FakeTopology",
        )

    def fake_validation(
        interchange,
        output_dir,
        *,
        settings=None,
        crosslinked_pdb_path=None,
        attachment_specs=(),
    ):
        calls["validation"] += 1
        calls["validation_crosslinked_pdb_path"] = crosslinked_pdb_path
        calls["validation_attachment_specs"] = attachment_specs
        return SimpleNamespace(success=True)

    def fake_local_minimization(pdb_path, output_dir, **kwargs):
        calls["local_minimization"] += 1
        calls["local_linkage_count"] = len(kwargs["settings"].linkages)
        calls["local_settings"] = kwargs["settings"]
        relaxed = Path(output_dir) / "local-relaxed.pdb"
        relaxed.write_text("END\n", encoding="utf-8")
        return SimpleNamespace(success=True, relaxed_pdb_path=relaxed, blocker=None)

    monkeypatch.setattr(workflow_module, "place_modifiers_with_resolved_plans", fake_place)
    monkeypatch.setattr(
        workflow_module, "placed_fragment_from_resolved_plan", lambda fragment, plan: fragment
    )
    monkeypatch.setattr(workflow_module, "write_crosslinked_pdb", fake_write)
    monkeypatch.setattr(workflow_module, "PabloIngestor", FakeIngestor)
    monkeypatch.setattr(
        workflow_module, "create_interchange_from_pablo_topology", fake_parameterize
    )
    monkeypatch.setattr(
        workflow_module,
        "run_post_crosslink_local_minimization",
        fake_local_minimization,
    )

    policy = ConjugationCcdPabloPolicyConfig(
        crosslinks=[
            {
                "residues": ("ASX", "NAG"),
                "linking_atoms": ("ND2", "C001"),
                "leaving_atoms": ((), ()),
                "bond_order": 1,
            }
        ]
    )

    specs = tuple(
        SimpleNamespace(
            attachment_id=f"attachment_{index:02d}",
            attachment_index=index,
            reaction_name="n_glycosylation",
            generated_fragment=modifier,
            source_fragment=moiety,
            resolved_plan=plan,
            source_sidecars={},
            fragment=SimpleNamespace(source_kind="moiety"),
        )
        for index, (modifier, moiety, plan) in enumerate(
            zip(modifiers, moieties, plans, strict=True),
            start=1,
        )
    )

    construction, topology = workflow_module._construct_conjugate_from_specs(
        protein_pdb_path=tmp_path / "protein.pdb",
        specs=specs,
        ccd_pablo_policy=policy,
        output_dir=tmp_path / "construction",
        chain_policy=None,
        settings=ModifierConstructionSettings(run_relaxation=False),
        use_product_state_pablo_library=False,
        run_product_state_local_minimization=True,
        local_minimization_settings=object(),
    )

    assert topology is not None
    assert len(construction.resolved_plans) == 2
    assert calls["placements"] == 1
    assert calls["placed_modifier_count"] == 2
    assert calls["placed_plan_count"] == 2
    assert calls["fragment_count"] == 2
    assert calls["attachment_count"] == 2
    assert calls["protein_products"] == ("ASX", "ASX")
    assert calls["parameterize"] == 1
    assert calls["local_minimization"] == 1
    assert calls["local_linkage_count"] == 2
    assert calls["local_settings"].linkages[0].protein_link_selector.residue_name == "ASX"
    assert calls["local_settings"].linkages[1].modifier_link_selector.residue_name == "NAG"
    assert calls["validation"] == 0
    assert construction.local_minimization is not None


def test_single_generic_local_minimization_request_runs_minimizer(
    monkeypatch,
    tmp_path: Path,
):
    """Single non-NHS attachments should use the generic local minimization path."""
    import polyzymd.builders.conjugation.system_workflow as workflow_module

    plan = _generic_resolved_plan(
        residue_number=42,
        modifier_residue_number=1,
        modifier_atom="C001",
    )
    modifier = _generated_fragment(residue_name="NAG", residue_number=1)
    spec = SimpleNamespace(
        attachment_id="n_gly_attachment_01",
        attachment_index=1,
        reaction_name="n_glycosylation",
        generated_fragment=modifier,
        resolved_plan=plan,
        source_sidecars={},
        fragment=SimpleNamespace(source_kind="moiety"),
    )
    calls = {"validation": 0, "local_minimization": 0}

    def fake_place(protein_pdb_path, modifiers_arg, plans_arg, output_dir, *, settings=None):
        return (
            PackmolModifierPlacementResult(
                placed_modifier=modifiers_arg[0].to_placed_fragment(),
                packmol_input_path=Path(output_dir) / "packmol.inp",
                packmol_output_path=Path(output_dir) / "packmol.pdb",
                protein_sterics_pdb_path=Path(output_dir) / "protein.pdb",
                modifier_pdb_path=Path(output_dir) / "modifier.pdb",
                packmol_input_text="",
                target_bond_length_angstrom=plan.target_bond_length_angstrom,
                placed_bond_length_angstrom=plan.target_bond_length_angstrom,
                min_modifier_protein_distance_angstrom=2.0,
            ),
        )

    def fake_write(protein_pdb_path, polymer_fragment, attachment, output_path, options):
        Path(output_path).write_text("END\n", encoding="utf-8")
        return CrosslinkedPdbAssemblyResult(
            output_path=Path(output_path),
            protein_atom_count=10,
            polymer_atom_count=1,
            removed_protein_atom_count=1,
            removed_polymer_atom_count=0,
            removed_atom_serials=(),
            removed_atom_names=(),
            added_conect_pair=(1, 2),
            added_conect_pairs=((1, 2),),
        )

    class FakeIngestor:
        def __init__(self, policy):
            self.policy = policy

        def ingest_structure(
            self, path, *, chain_policy=None, output_dir=None, residue_library=None
        ):
            return PabloIngestionResult(
                success=True,
                path=Path(path),
                suffix=".pdb",
                pablo=PabloAvailability(available=True),
                topology=SimpleNamespace(molecules=(object(),)),
            )

    def fake_parameterize(
        topology,
        *,
        settings=None,
        charge_from_molecules=None,
        require_charge_templates=False,
    ):
        calls["require_charge_templates"] = require_charge_templates
        return InterchangeParameterizationResult(
            success=True,
            interchange=object(),
            force_field_names=("fake.offxml",),
            topology_type="FakeTopology",
        )

    def fake_validation(
        interchange,
        output_dir,
        *,
        settings=None,
        crosslinked_pdb_path=None,
        attachment_specs=(),
    ):
        calls["validation"] += 1
        calls["validation_crosslinked_pdb_path"] = crosslinked_pdb_path
        calls["validation_attachment_specs"] = attachment_specs
        return SimpleNamespace(success=True)

    def fake_local_minimization(*args, **kwargs):
        calls["local_minimization"] += 1
        relaxed = tmp_path / "local-relaxed.pdb"
        relaxed.write_text("END\n", encoding="utf-8")
        return SimpleNamespace(success=True, relaxed_pdb_path=relaxed)

    monkeypatch.setattr(workflow_module, "place_modifiers_with_resolved_plans", fake_place)
    monkeypatch.setattr(
        workflow_module, "placed_fragment_from_resolved_plan", lambda fragment, plan: fragment
    )
    monkeypatch.setattr(workflow_module, "write_crosslinked_pdb", fake_write)
    monkeypatch.setattr(workflow_module, "PabloIngestor", FakeIngestor)
    monkeypatch.setattr(
        workflow_module, "create_interchange_from_pablo_topology", fake_parameterize
    )
    monkeypatch.setattr(
        workflow_module,
        "run_post_crosslink_local_minimization",
        fake_local_minimization,
    )
    monkeypatch.setattr(
        workflow_module,
        "_local_minimization_settings_for_product",
        lambda *args, **kwargs: object(),
    )

    policy = ConjugationCcdPabloPolicyConfig(
        crosslinks=[
            {
                "residues": ("ASX", "NAG"),
                "linking_atoms": ("ND2", "C001"),
                "leaving_atoms": ((), ()),
                "bond_order": 1,
            }
        ]
    )

    construction, topology = workflow_module._construct_conjugate_from_specs(
        protein_pdb_path=tmp_path / "protein.pdb",
        specs=(spec,),
        ccd_pablo_policy=policy,
        output_dir=tmp_path / "construction",
        chain_policy=None,
        settings=ModifierConstructionSettings(run_relaxation=False),
        use_product_state_pablo_library=False,
        use_conjugate_relaxation=False,
        run_product_state_local_minimization=True,
        local_minimization_settings=object(),
    )

    assert topology is not None
    assert construction.local_minimization is not None
    assert construction.relaxation is None
    assert calls["local_minimization"] == 1
    assert calls["validation"] == 0


@pytest.mark.parametrize(
    (
        "run_relaxation",
        "use_conjugate_relaxation",
        "expect_relaxation",
        "relaxation_evidence_status",
    ),
    (
        (False, False, False, ValidationStatus.PASS),
        (True, True, True, ValidationStatus.PASS),
        (True, True, True, ValidationStatus.FAIL),
    ),
)
def test_relaxation_receives_product_path_and_attachment_specs(
    monkeypatch,
    tmp_path: Path,
    run_relaxation: bool,
    use_conjugate_relaxation: bool,
    expect_relaxation: bool,
    relaxation_evidence_status: ValidationStatus,
):
    """Conjugate relaxation should use the current product-aware signature."""
    import polyzymd.builders.conjugation.system_workflow as workflow_module

    plan = _generic_resolved_plan(
        residue_number=42,
        modifier_residue_number=1,
        modifier_atom="C001",
    )
    modifier = _generated_fragment(residue_name="NAG", residue_number=1)
    spec = SimpleNamespace(
        attachment_id="n_gly_attachment_01",
        attachment_index=1,
        reaction_name="n_glycosylation",
        generated_fragment=modifier,
        resolved_plan=plan,
        source_sidecars={},
        fragment=SimpleNamespace(source_kind="moiety"),
    )
    calls = {}

    def fake_place(protein_pdb_path, modifiers_arg, plans_arg, output_dir, *, settings=None):
        return (
            PackmolModifierPlacementResult(
                placed_modifier=modifiers_arg[0].to_placed_fragment(),
                packmol_input_path=Path(output_dir) / "packmol.inp",
                packmol_output_path=Path(output_dir) / "packmol.pdb",
                protein_sterics_pdb_path=Path(output_dir) / "protein.pdb",
                modifier_pdb_path=Path(output_dir) / "modifier.pdb",
                packmol_input_text="",
                target_bond_length_angstrom=plan.target_bond_length_angstrom,
                placed_bond_length_angstrom=plan.target_bond_length_angstrom,
                min_modifier_protein_distance_angstrom=2.0,
            ),
        )

    def fake_write(protein_pdb_path, polymer_fragment, attachment, output_path, options):
        Path(output_path).write_text("END\n", encoding="utf-8")
        return CrosslinkedPdbAssemblyResult(
            output_path=Path(output_path),
            protein_atom_count=10,
            polymer_atom_count=1,
            removed_protein_atom_count=1,
            removed_polymer_atom_count=0,
            removed_atom_serials=(),
            removed_atom_names=(),
            added_conect_pair=(1, 2),
            added_conect_pairs=((1, 2),),
        )

    class FakeIngestor:
        def __init__(self, policy):
            self.policy = policy

        def ingest_structure(
            self, path, *, chain_policy=None, output_dir=None, residue_library=None
        ):
            return PabloIngestionResult(
                success=True,
                path=Path(path),
                suffix=".pdb",
                pablo=PabloAvailability(available=True),
                topology=SimpleNamespace(molecules=(object(),)),
            )

    def fake_relax_conjugate(
        interchange,
        output_dir,
        *,
        product_pdb_path=None,
        attachment_specs=(),
        assembly=None,
        settings=None,
    ):
        calls["relaxation_product_pdb_path"] = product_pdb_path
        calls["relaxation_attachment_specs"] = attachment_specs
        return SimpleNamespace(success=True, relaxed_pdb_path=Path(output_dir) / "relaxed.pdb")

    class RaisingInterchange:
        """Interchange fake that must not be converted by validation reporting."""

        def to_openmm_system(self):
            """Raise if validation attempts to create its own OpenMM system."""
            raise AssertionError("validation must not call to_openmm_system")

    def fake_validation_report(**kwargs):
        calls["validation_kwargs"] = kwargs
        return SimpleNamespace(
            report_path=Path(kwargs["output_dir"]) / "validation.json",
            relaxation_evidence=SimpleNamespace(status=relaxation_evidence_status),
        )

    monkeypatch.setattr(workflow_module, "place_modifiers_with_resolved_plans", fake_place)
    monkeypatch.setattr(
        workflow_module, "placed_fragment_from_resolved_plan", lambda fragment, plan: fragment
    )
    monkeypatch.setattr(workflow_module, "write_crosslinked_pdb", fake_write)
    monkeypatch.setattr(workflow_module, "PabloIngestor", FakeIngestor)
    monkeypatch.setattr(
        workflow_module,
        "create_interchange_from_pablo_topology",
        lambda topology, **kwargs: InterchangeParameterizationResult(
            success=True,
            interchange=RaisingInterchange(),
            force_field_names=("fake.offxml",),
            topology_type="FakeTopology",
        ),
    )
    monkeypatch.setattr(workflow_module, "relax_conjugate", fake_relax_conjugate)
    monkeypatch.setattr(
        workflow_module,
        "build_conjugate_validation_report",
        fake_validation_report,
    )

    construction_kwargs = {
        "protein_pdb_path": tmp_path / "protein.pdb",
        "specs": (spec,),
        "ccd_pablo_policy": ConjugationCcdPabloPolicyConfig(
            crosslinks=[
                {
                    "residues": ("ASX", "NAG"),
                    "linking_atoms": ("ND2", "C001"),
                    "leaving_atoms": ((), ()),
                    "bond_order": 1,
                }
            ]
        ),
        "output_dir": tmp_path / "construction",
        "chain_policy": None,
        "settings": ModifierConstructionSettings(run_relaxation=run_relaxation),
        "use_product_state_pablo_library": False,
        "use_conjugate_relaxation": use_conjugate_relaxation,
        "run_product_state_local_minimization": False,
    }

    if relaxation_evidence_status == ValidationStatus.FAIL:
        with pytest.raises(RuntimeError, match="relaxation evidence failed"):
            workflow_module._construct_conjugate_from_specs(**construction_kwargs)
        assert ("relaxation_product_pdb_path" in calls) is expect_relaxation
        assert "interchange" not in calls["validation_kwargs"]
        return

    construction, _topology = workflow_module._construct_conjugate_from_specs(**construction_kwargs)

    assert ("relaxation_product_pdb_path" in calls) is expect_relaxation
    if expect_relaxation:
        assert calls["relaxation_product_pdb_path"] == construction.crosslinked_pdb_path
        assert calls["relaxation_attachment_specs"] == (spec,)
    assert "interchange" not in calls["validation_kwargs"]


def test_construction_final_interchange_uses_strict_charge_bridge(
    monkeypatch,
    tmp_path: Path,
):
    """Final Interchange should pass product and standard charge templates strictly."""
    import polyzymd.builders.conjugation.pablo.charge_bridge as charge_bridge_module
    import polyzymd.builders.conjugation.system_workflow as workflow_module
    from polyzymd.builders.conjugation.final_interchange import (
        create_final_conjugated_interchange,
    )

    plan = _generic_resolved_plan(
        residue_number=42,
        modifier_residue_number=1,
        modifier_atom="C001",
    )
    modifier = _generated_fragment(residue_name="NAG", residue_number=1)
    spec = SimpleNamespace(
        attachment_id="n_gly_attachment_01",
        attachment_index=1,
        reaction_name="n_glycosylation",
        generated_fragment=modifier,
        resolved_plan=plan,
        source_sidecars={},
        fragment=SimpleNamespace(source_kind="moiety"),
    )
    product_template = _charged_molecule_like("product", residue_name="NAG")
    product_template.properties = {"polyzymd_charge_provenance": "production:test"}
    standard_template = _charged_molecule_like("water", residue_name="HOH")
    product_library = SimpleNamespace(
        definitions=(SimpleNamespace(residue_name="NAG"),),
        residue_names=("NAG",),
        residue_library=object(),
    )

    def fake_place(protein_pdb_path, modifiers_arg, plans_arg, output_dir, *, settings=None):
        return (
            PackmolModifierPlacementResult(
                placed_modifier=modifiers_arg[0].to_placed_fragment(),
                packmol_input_path=Path(output_dir) / "packmol.inp",
                packmol_output_path=Path(output_dir) / "packmol.pdb",
                protein_sterics_pdb_path=Path(output_dir) / "protein.pdb",
                modifier_pdb_path=Path(output_dir) / "modifier.pdb",
                packmol_input_text="",
                target_bond_length_angstrom=plan.target_bond_length_angstrom,
                placed_bond_length_angstrom=plan.target_bond_length_angstrom,
                min_modifier_protein_distance_angstrom=2.0,
            ),
        )

    def fake_write(protein_pdb_path, polymer_fragment, attachment, output_path, options):
        Path(output_path).write_text("END\n", encoding="utf-8")
        return CrosslinkedPdbAssemblyResult(
            output_path=Path(output_path),
            protein_atom_count=10,
            polymer_atom_count=1,
            removed_protein_atom_count=1,
            removed_polymer_atom_count=0,
            removed_atom_serials=(),
            removed_atom_names=(),
            added_conect_pair=(1, 2),
            added_conect_pairs=((1, 2),),
        )

    class FakeIngestor:
        def __init__(self, policy):
            self.policy = policy

        def ingest_structure(
            self, path, *, chain_policy=None, output_dir=None, residue_library=None
        ):
            return PabloIngestionResult(
                success=True,
                path=Path(path),
                suffix=".pdb",
                pablo=PabloAvailability(available=True),
                topology=SimpleNamespace(molecules=(product_template,)),
            )

    intermediate = {}

    def fake_parameterize(
        topology,
        *,
        settings=None,
        charge_from_molecules=None,
        require_charge_templates=False,
    ):
        intermediate["charge_from_molecules"] = tuple(charge_from_molecules or ())
        intermediate["require_charge_templates"] = require_charge_templates
        return InterchangeParameterizationResult(
            success=True,
            interchange=object(),
            force_field_names=("fake.offxml",),
            topology_type="FakeTopology",
        )

    monkeypatch.setattr(workflow_module, "place_modifiers_with_resolved_plans", fake_place)
    monkeypatch.setattr(
        workflow_module, "placed_fragment_from_resolved_plan", lambda fragment, plan: fragment
    )
    monkeypatch.setattr(workflow_module, "write_crosslinked_pdb", fake_write)
    monkeypatch.setattr(workflow_module, "PabloIngestor", FakeIngestor)
    monkeypatch.setattr(
        workflow_module,
        "_product_state_pablo_library_for_specs",
        lambda **kwargs: product_library,
    )
    monkeypatch.setattr(
        charge_bridge_module,
        "build_product_state_charge_bridge",
        lambda **kwargs: kwargs["product_state_pablo_library"],
    )
    monkeypatch.setattr(
        workflow_module, "create_interchange_from_pablo_topology", fake_parameterize
    )

    construction, _topology = workflow_module._construct_conjugate_from_specs(
        protein_pdb_path=tmp_path / "protein.pdb",
        specs=(spec,),
        ccd_pablo_policy=ConjugationCcdPabloPolicyConfig(
            crosslinks=[
                {
                    "residues": ("ASX", "NAG"),
                    "linking_atoms": ("ND2", "C001"),
                    "leaving_atoms": ((), ()),
                    "bond_order": 1,
                }
            ]
        ),
        output_dir=tmp_path / "construction",
        chain_policy=None,
        settings=ModifierConstructionSettings(run_relaxation=False),
        use_product_state_pablo_library=True,
        run_product_state_local_minimization=False,
    )
    captured = {}
    builder = SimpleNamespace(
        _solvated_topology=SimpleNamespace(molecules=(product_template, standard_template)),
        _interchange=None,
        collect_standard_charge_templates=lambda: (standard_template,),
    )

    def fake_final_parameterizer(topology, **kwargs):
        captured["kwargs"] = kwargs
        return SimpleNamespace(success=True, interchange=object())

    create_final_conjugated_interchange(
        builder,
        product_state_pablo_library=construction.product_state_pablo_library,
        parameterizer=fake_final_parameterizer,
    )

    templates = tuple(captured["kwargs"]["charge_from_molecules"])
    assert construction.product_state_pablo_library.charge_templates == (product_template,)
    assert len(intermediate["charge_from_molecules"]) == 1
    assert intermediate["charge_from_molecules"][0] is not product_template
    assert intermediate["require_charge_templates"] is True
    assert len(templates) == 2
    assert templates[0] is not product_template
    assert templates[1] is standard_template
    assert captured["kwargs"]["require_charge_templates"] is True


def test_direct_n_gly_path_builds_specs_before_construction(monkeypatch, tmp_path: Path):
    """Direct N-gly requests should resolve one build spec per enabled attachment."""
    import polyzymd.builders.conjugation.system_workflow as workflow_module

    source = tmp_path / "protein.pdb"
    source.write_text("END\n", encoding="utf-8")
    relaxed = tmp_path / "relaxed.pdb"
    relaxed.write_text("END\n", encoding="utf-8")
    attachments = tuple(_direct_n_gly_attachment(index) for index in (1, 2))
    built_specs = []
    calls = {}
    real_spec_builder = workflow_module.attachment_spec_from_generated_polymer_plan

    class FakeReaction:
        coordinate_backend_mechanism = "n_glycosylation"

        @staticmethod
        def resolve_plan(protein_pdb_path, site, fragment, *, settings=None):
            return _generic_resolved_plan(
                residue_number=site.residue_number,
                modifier_residue_number=fragment.residues[0].residue_number,
                modifier_atom="C001",
            )

    def counting_spec_builder(*args, **kwargs):
        spec = real_spec_builder(*args, **kwargs)
        built_specs.append(spec)
        return spec

    def fake_construct(**kwargs):
        calls["construct_spec_count"] = len(built_specs)
        calls["spec_count"] = len(kwargs["specs"])
        calls["attachment_specs"] = kwargs["specs"]
        calls["resolved_plan_count"] = len([spec.resolved_plan for spec in kwargs["specs"]])
        calls["run_local_minimization"] = kwargs["run_product_state_local_minimization"]
        calls["local_minimization_settings"] = kwargs["local_minimization_settings"]
        return (
            SimpleNamespace(
                relaxation=SimpleNamespace(minimized_pdb_path=relaxed, relaxed_pdb_path=relaxed),
                crosslinked_pdb_path=tmp_path / "crosslinked.pdb",
                diagnostics=("fake construction",),
            ),
            object(),
        )

    monkeypatch.setattr(workflow_module, "get_reaction", lambda name: FakeReaction)

    def fake_resolve_source(attachment, *, attachment_index, **kwargs):
        fragment = _moiety_fragment(
            residue_name=attachment.moiety.residue_name,
            residue_number=attachment_index,
        )
        return SimpleNamespace(
            fragment=fragment,
            source_fragment=fragment,
            source_kind="smiles",
            sidecars={},
            generation=None,
            reactive_sequence_index=None,
            reactive_selector=None,
        )

    monkeypatch.setattr(workflow_module, "resolve_moiety_source", fake_resolve_source)
    monkeypatch.setattr(
        workflow_module,
        "attachment_spec_from_generated_polymer_plan",
        counting_spec_builder,
    )
    monkeypatch.setattr(
        workflow_module,
        "_prepared_protein_pdb_path",
        lambda protein_pdb_path, **kwargs: (Path(protein_pdb_path), None),
    )
    monkeypatch.setattr(workflow_module, "_construct_conjugate_from_specs", fake_construct)
    monkeypatch.setattr(
        workflow_module, "topology_with_pdb_positions", lambda topology, path: topology
    )
    monkeypatch.setattr(workflow_module, "_restore_pdb_atom_name_fields", lambda *args: 0)
    monkeypatch.setattr(
        workflow_module,
        "_build_direct_solvated_system",
        lambda **kwargs: _FakeSystemBuilder(tmp_path / "solvated.pdb"),
    )

    result = workflow_module.build_direct_moiety_conjugate(
        protein_pdb_path=source,
        attachments=attachments,
        output_dir=tmp_path / "out",
        settings=workflow_module.ConjugatedPolymerSystemSettings(
            run_product_state_local_minimization=True
        ),
    )

    assert len(built_specs) == 2
    assert calls["construct_spec_count"] == 2
    assert calls["spec_count"] == 2
    assert calls["resolved_plan_count"] == 2
    assert calls["attachment_specs"] == tuple(built_specs)
    assert calls["run_local_minimization"] is True
    assert calls["local_minimization_settings"] is not None
    assert result.workflow_json_path.exists()
    payload = json.loads(result.workflow_json_path.read_text(encoding="utf-8"))
    assert payload["workflow_json_path"] == str(result.workflow_json_path)
    assert payload["artifact_paths"]["workflow_json"] == str(result.workflow_json_path)


def test_config_nhs_lys_path_builds_specs_before_shared_construction(
    monkeypatch,
    tmp_path: Path,
    caplog: pytest.LogCaptureFixture,
):
    """Config NHS-Lys construction should resolve all enabled specs then construct once."""
    import polyzymd.builders.conjugation.system_workflow as workflow_module

    source = tmp_path / "protein.pdb"
    source.write_text("END\n", encoding="utf-8")
    polymer_pdb = tmp_path / "polymer.pdb"
    polymer_pdb.write_text(
        _pdb_line(1, "C001", "NHS", "C", 1, element="C") + "END\n",
        encoding="utf-8",
    )
    polymer_sdf = tmp_path / "polymer.sdf"
    polymer_sdf.write_text("", encoding="utf-8")
    relaxed = tmp_path / "relaxed.pdb"
    relaxed.write_text("END\n", encoding="utf-8")
    attachments = (
        _config_nhs_attachment("nhs_polymer_1", residue_number=23),
        _config_nhs_attachment("nhs_polymer_2", residue_number=45),
    )
    config = SimpleNamespace(
        enzyme=SimpleNamespace(pdb_path=source),
        conjugation=SimpleNamespace(
            enabled=True,
            attachments=attachments,
            ccd_pablo=ConjugationCcdPabloPolicyConfig(),
            chain_policy=ConjugationChainPolicyConfig(),
        ),
    )
    generation = MultiResidueGenerationResult(
        recipe_name="test",
        sequence="AA",
        cache_directory=tmp_path / "cache",
        monomer_group_path=tmp_path / "monomers.json",
        pdb_path=polymer_pdb,
        sdf_path=polymer_sdf,
        object_type="FakePolymer",
        atom_count=1,
    )
    built_specs = []
    calls = {}
    real_spec_builder = workflow_module.attachment_spec_from_generated_polymer_plan

    class FakeReaction:
        name = "nhs_lys"
        coordinate_backend_mechanism = "nhs_lys_amide"

        @staticmethod
        def resolve_plan(protein_pdb_path, site, fragment, *, settings=None):
            calls.setdefault("reaction_fragments", []).append(fragment)
            return _resolved_plan(
                PabloCrosslinkRequirement(
                    residues=("LYX", "NHX"),
                    linking_atoms=("NZ", "C047"),
                    leaving_atoms=((), ()),
                    bond_order=1,
                )
            )

    def counting_spec_builder(*args, **kwargs):
        spec = real_spec_builder(*args, **kwargs)
        built_specs.append(spec)
        return spec

    def fake_construct(**kwargs):
        calls["construct_spec_count"] = len(built_specs)
        calls["specs"] = kwargs["specs"]
        calls["run_local_minimization"] = kwargs["run_product_state_local_minimization"]
        calls["local_minimization_settings"] = kwargs["local_minimization_settings"]
        return (
            SimpleNamespace(
                relaxation=SimpleNamespace(minimized_pdb_path=relaxed, relaxed_pdb_path=relaxed)
            ),
            object(),
        )

    monkeypatch.setattr(workflow_module, "get_reaction", lambda name: FakeReaction())
    monkeypatch.setattr(
        workflow_module,
        "resolve_moiety_source",
        lambda *args, **kwargs: SimpleNamespace(
            fragment=_generated_fragment(residue_name="NHS", residue_number=1),
            source_fragment=_generated_fragment(residue_name="NHS", residue_number=1),
            source_kind="polymer",
            sidecars={"sdf": polymer_sdf},
            generation=generation,
            reactive_sequence_index=0,
            reactive_selector={"residue_name": "NHS", "residue_number": 1},
        ),
    )
    monkeypatch.setattr(
        workflow_module,
        "attachment_spec_from_generated_polymer_plan",
        counting_spec_builder,
    )
    monkeypatch.setattr(workflow_module, "_construct_conjugate_from_specs", fake_construct)
    monkeypatch.setattr(workflow_module, "_relaxed_conjugate_pdb", lambda construction: relaxed)
    monkeypatch.setattr(
        workflow_module,
        "topology_with_pdb_positions",
        lambda topology, path, **kwargs: topology,
    )
    monkeypatch.setattr(
        workflow_module,
        "_build_solvated_system",
        lambda *args, **kwargs: _FakeSystemBuilder(tmp_path / "solvated.pdb"),
    )
    monkeypatch.setattr(workflow_module, "_restore_pdb_atom_name_fields", lambda *args: 0)

    settings = workflow_module.ConjugatedPolymerSystemSettings(
        canonicalize_source_protein_hydrogens=False,
        preserve_reference_atom_names=False,
        run_product_state_local_minimization=True,
    )
    caplog.set_level(logging.INFO, logger=workflow_module.__name__)
    result = workflow_module.build_conjugated_polymer_system_from_config(
        config,
        output_dir=tmp_path / "out",
        settings=settings,
    )

    assert len(built_specs) == 2
    assert calls["construct_spec_count"] == 2
    assert calls["specs"] == tuple(built_specs)
    assert calls["run_local_minimization"] is True
    assert calls["local_minimization_settings"] is settings.local_minimization
    assert result.generated_sequences == ("AA", "AA")
    assert result.reactive_sequence_indices == (0, 0)
    assert result.modifiers == tuple(spec.generated_fragment for spec in built_specs)
    assert result.workflow_json_path.exists()
    payload = json.loads(result.workflow_json_path.read_text(encoding="utf-8"))
    assert payload["workflow_json_path"] == str(result.workflow_json_path)
    assert payload["artifact_paths"]["workflow_json"] == str(result.workflow_json_path)
    messages = [record.getMessage() for record in caplog.records]
    addition_messages = [message for message in messages if message.startswith("Addition ")]
    assert "Enabled conjugation attachment count: 2" in messages
    assert addition_messages == [
        "Addition 1/2: adding nhs_polymer_1 to A:LYS23:NZ",
        "Addition 2/2: adding nhs_polymer_2 to A:LYS45:NZ",
    ]
    assert any("Generating conjugate polymer/moiety" in message for message in messages)
    assert any("Preparing relaxed conjugate topology" in message for message in messages)
    assert any("Saved conjugation workflow JSON" in message for message in messages)


def test_config_nhs_lys_path_still_accepts_one_attachment(monkeypatch, tmp_path: Path):
    """The existing single-attachment config workflow remains valid."""
    import polyzymd.builders.conjugation.system_workflow as workflow_module

    source = tmp_path / "protein.pdb"
    source.write_text("END\n", encoding="utf-8")
    attachment = _config_nhs_attachment("nhs_polymer", residue_number=23)
    config = SimpleNamespace(
        enzyme=SimpleNamespace(pdb_path=source),
        conjugation=SimpleNamespace(
            enabled=True,
            attachments=(attachment,),
            ccd_pablo=ConjugationCcdPabloPolicyConfig(),
            chain_policy=ConjugationChainPolicyConfig(),
        ),
    )
    generation = MultiResidueGenerationResult(
        recipe_name="test",
        sequence="AA",
        cache_directory=tmp_path / "cache",
        monomer_group_path=tmp_path / "monomers.json",
        pdb_path=tmp_path / "polymer.pdb",
        sdf_path=tmp_path / "polymer.sdf",
        object_type="FakePolymer",
        atom_count=1,
    )
    generation.pdb_path.write_text(_pdb_line(1, "C001", "NHS", "C", 1), encoding="utf-8")
    generation.sdf_path.write_text("", encoding="utf-8")
    relaxed = tmp_path / "relaxed.pdb"
    relaxed.write_text("END\n", encoding="utf-8")

    class FakeReaction:
        name = "nhs_lys"
        coordinate_backend_mechanism = "nhs_lys_amide"

        @staticmethod
        def resolve_plan(protein_pdb_path, site, fragment, *, settings=None):
            return _resolved_plan(
                PabloCrosslinkRequirement(
                    residues=("LYX", "NHX"),
                    linking_atoms=("NZ", "C047"),
                    leaving_atoms=((), ()),
                    bond_order=1,
                )
            )

    monkeypatch.setattr(workflow_module, "get_reaction", lambda name: FakeReaction())
    monkeypatch.setattr(
        workflow_module,
        "resolve_moiety_source",
        lambda *args, **kwargs: SimpleNamespace(
            fragment=_generated_fragment(residue_name="NHS", residue_number=1),
            source_fragment=_generated_fragment(residue_name="NHS", residue_number=1),
            source_kind="polymer",
            sidecars={"sdf": generation.sdf_path},
            generation=generation,
            reactive_sequence_index=0,
            reactive_selector={"residue_name": "NHS", "residue_number": 1},
        ),
    )
    monkeypatch.setattr(
        workflow_module,
        "_construct_conjugate_from_specs",
        lambda **kwargs: (
            SimpleNamespace(
                relaxation=SimpleNamespace(minimized_pdb_path=relaxed, relaxed_pdb_path=relaxed)
            ),
            object(),
        ),
    )
    monkeypatch.setattr(workflow_module, "_relaxed_conjugate_pdb", lambda construction: relaxed)
    monkeypatch.setattr(
        workflow_module, "topology_with_pdb_positions", lambda topology, path, **kwargs: topology
    )
    monkeypatch.setattr(
        workflow_module,
        "_build_solvated_system",
        lambda *args, **kwargs: _FakeSystemBuilder(tmp_path / "solvated.pdb"),
    )
    monkeypatch.setattr(workflow_module, "_restore_pdb_atom_name_fields", lambda *args: 0)

    result = workflow_module.build_conjugated_polymer_system_from_config(
        config,
        output_dir=tmp_path / "out",
        settings=workflow_module.ConjugatedPolymerSystemSettings(
            canonicalize_source_protein_hydrogens=False,
            preserve_reference_atom_names=False,
        ),
    )

    assert result.generated_sequence == "AA"
    assert result.workflow_json_path.exists()


def test_attachment_spec_uses_generic_reaction_and_smiles_provider(monkeypatch, tmp_path: Path):
    """Generic spec resolution should use the provider and reaction template for any source."""
    import polyzymd.builders.conjugation.system_workflow as workflow_module

    source_fragment = _moiety_fragment(residue_name="NAG", residue_number=1)
    generated_fragment = _generated_fragment(residue_name="NAG", residue_number=1)
    plan = _generic_resolved_plan(
        residue_number=42,
        modifier_residue_number=1,
        modifier_atom="C001",
    )
    calls = {}

    class FakeReaction:
        name = "n_glycosylation"
        coordinate_backend_mechanism = "n_glycosylation"

        @staticmethod
        def resolve_plan(protein_pdb_path, site, fragment, *, settings=None):
            calls["protein_pdb_path"] = protein_pdb_path
            calls["site"] = site
            calls["fragment"] = fragment
            return plan

    attachment = SimpleNamespace(
        name="glycan",
        enabled=True,
        site=SimpleNamespace(chain_id="A", residue_name="ASN", residue_number=42, atom_name="ND2"),
        moiety=SimpleNamespace(name="nag", smiles="CO", residue_name="NAG"),
        mechanism=SimpleNamespace(name="n_glycosylation", reaction_smarts=None),
    )

    monkeypatch.setattr(workflow_module, "get_reaction", lambda name: FakeReaction())
    monkeypatch.setattr(
        workflow_module,
        "resolve_moiety_source",
        lambda *args, **kwargs: SimpleNamespace(
            fragment=generated_fragment,
            source_fragment=source_fragment,
            source_kind="smiles",
            sidecars={"sdf": tmp_path / "nag.sdf"},
            generation=None,
            reactive_sequence_index=None,
            reactive_selector=None,
        ),
    )
    monkeypatch.setattr(
        workflow_module,
        "generated_fragment_for_resolved_source",
        lambda source, resolved_plan: generated_fragment,
    )

    spec, generation, reactive_index, reactive_selector, reaction = (
        workflow_module._build_attachment_spec(
            attachment,
            attachment_index=1,
            protein_pdb_path=tmp_path / "protein.pdb",
            artifact_dir=tmp_path,
            workflow_settings=workflow_module.ConjugatedPolymerSystemSettings(),
        )
    )

    assert spec.reaction_name == "n_glycosylation"
    assert spec.fragment.source_kind == "moiety"
    assert spec.source_fragment is source_fragment
    assert spec.generated_fragment is generated_fragment
    assert generation is None
    assert reactive_index == 0
    assert reactive_selector == {}
    assert reaction.name == "n_glycosylation"
    assert calls["fragment"] is source_fragment


def test_config_nhs_lys_mixed_unsupported_mechanisms_fail_clearly(tmp_path: Path):
    """Mixed config workflows should not silently enter the NHS-Lys construction path."""
    import polyzymd.builders.conjugation.system_workflow as workflow_module

    source = tmp_path / "protein.pdb"
    source.write_text("END\n", encoding="utf-8")
    config = SimpleNamespace(
        enzyme=SimpleNamespace(pdb_path=source),
        conjugation=SimpleNamespace(
            enabled=True,
            attachments=(
                _config_nhs_attachment("nhs_polymer", residue_number=23),
                SimpleNamespace(
                    name="custom",
                    enabled=True,
                    moiety=SimpleNamespace(name="custom", polymer_recipe=object()),
                    mechanism=SimpleNamespace(name="custom", reaction_smarts=None),
                ),
            ),
            ccd_pablo=ConjugationCcdPabloPolicyConfig(),
            chain_policy=ConjugationChainPolicyConfig(),
        ),
    )

    with pytest.raises(NotImplementedError, match="coordinate surgery only"):
        workflow_module.build_conjugated_polymer_system_from_config(
            config,
            output_dir=tmp_path / "out",
            settings=workflow_module.ConjugatedPolymerSystemSettings(
                canonicalize_source_protein_hydrogens=False
            ),
        )


def test_coordinate_backend_gate_allows_explicit_nhs_lys_mechanism():
    """The system workflow should only route the named implemented backend to NHS-Lys."""
    attachment = SimpleNamespace(
        mechanism=SimpleNamespace(name="nhs_lys_amide", reaction_smarts=None)
    )

    _require_supported_coordinate_backend(attachment)


def test_coordinate_backend_gate_preflights_generic_smarts_then_blocks():
    """Generic reaction SMARTS should not silently enter the NHS-Lys coordinate path."""
    attachment = SimpleNamespace(
        mechanism=SimpleNamespace(
            name="generic_amide",
            reaction_smarts="[N:1]([H:2]).[C:3](=[O:4])[O:5]>>[N:1][C:3](=[O:4])",
            atom_roles=[
                {"map_number": 1, "participant": "site", "role": "linking"},
                {"map_number": 2, "participant": "site", "role": "leaving"},
                {"map_number": 3, "participant": "moiety", "role": "linking"},
                {"map_number": 5, "participant": "moiety", "role": "leaving"},
            ],
        )
    )

    with pytest.raises(NotImplementedError, match="generic SMARTS preflight") as excinfo:
        _require_supported_coordinate_backend(attachment)

    message = str(excinfo.value)
    assert "1 added" in message
    assert "coordinate surgery only for mechanism 'nhs_lys_amide'" in message


def test_coordinate_backend_gate_blocks_unspecified_mechanisms_without_smarts():
    """Unsupported mechanisms without SMARTS should fail before polymer generation."""
    attachment = SimpleNamespace(mechanism=SimpleNamespace(name="custom", reaction_smarts=None))

    with pytest.raises(NotImplementedError, match="currently implements coordinate surgery only"):
        _require_supported_coordinate_backend(attachment)


class _FakeSystemBuilder:
    def __init__(self, solvated_topology):
        self.solvated_topology = solvated_topology
        self.interchange = None

    def save_topology(self, path):
        Path(path).parent.mkdir(parents=True, exist_ok=True)
        Path(path).write_text("END\n", encoding="utf-8")


class _NoKwargInterchangeBuilder(_FakeSystemBuilder):
    """Fake builder whose final Interchange method rejects stale kwargs."""

    def __init__(self, solvated_topology):
        super().__init__(solvated_topology)
        self._combined_topology = object()
        self._solvent_builder = SimpleNamespace(
            solvated_topology=solvated_topology,
            solvate_from_config=lambda topology, solvent: None,
        )
        self.create_interchange_calls = 0

    def combine_solutes(self):
        """Record solute combination without building real OpenFF objects."""
        self.combined_solutes = True

    def solvate(self, *, padding, box_shape):
        """Record direct solvation without building real OpenFF objects."""
        self.solvation_settings = {"padding": padding, "box_shape": box_shape}

    def create_interchange(self):
        """Record final Interchange creation while accepting no kwargs."""
        self.create_interchange_calls += 1
        raise AssertionError("generic SystemBuilder.create_interchange must not be called")


def _fake_product_state_library():
    """Build product-state provenance for workflow helper tests."""
    return SimpleNamespace(
        definitions=(SimpleNamespace(residue_name="LYX"),),
        residue_names=("LYX",),
    )


def _charged_molecule_like(smiles: str, *, residue_name: str) -> SimpleNamespace:
    """Build a charged molecule-like object with residue metadata."""
    return SimpleNamespace(
        partial_charges=(0.0,),
        atoms=(SimpleNamespace(metadata={"residue_name": residue_name, "atom_name": "C001"}),),
        to_smiles=lambda mapped=False: smiles,
    )


def _direct_n_gly_attachment(index: int) -> SimpleNamespace:
    return SimpleNamespace(
        name=f"glycan_{index}",
        enabled=True,
        site=SimpleNamespace(residue_number=40 + index),
        moiety=SimpleNamespace(
            smiles="CO",
            residue_name="NAG",
            name=f"nag_{index}",
        ),
        mechanism=SimpleNamespace(name="n_glycosylation"),
    )


def _config_nhs_attachment(
    name: str = "nhs_polymer", *, residue_number: int = 23
) -> SimpleNamespace:
    recipe = PolymerRecipe(
        name="test_nhs_polymer",
        monomers=(
            PolymerMonomerRecipe(
                label="A",
                name="NHS",
                residue_name="NHS",
                smiles="C",
                probability=1.0,
            ),
        ),
        length=2,
        reactive_monomer_index=0,
        fixed_sequence="AA",
    )
    return SimpleNamespace(
        name=name,
        enabled=True,
        moiety=SimpleNamespace(name=name, polymer_recipe=recipe),
        mechanism=SimpleNamespace(name="nhs_lys_amide", reaction_smarts=None),
        site=SimpleNamespace(
            chain_id="A",
            residue_name="LYS",
            residue_number=residue_number,
            insertion_code="",
            atom_name="NZ",
        ),
    )


def _pdb_line(
    serial: int,
    atom_name_field: str,
    residue_name: str,
    chain_id: str,
    residue_number: int,
    *,
    element: str = "C",
) -> str:
    return (
        f"ATOM  {serial:5d} {atom_name_field[:4]} {residue_name:>3} {chain_id}"
        f"{residue_number:4d}       0.000   0.000   0.000  1.00  0.00"
        f"          {element:>2}  \n"
    )


def _resolved_plan(requirement: PabloCrosslinkRequirement) -> ResolvedAttachmentPlan:
    protein_selector = PdbAtomSelector(
        chain_id="A",
        residue_name="LYS",
        residue_number=23,
        atom_name="NZ",
    )
    modifier_selector = PdbAtomSelector(
        chain_id="C",
        residue_name="NHS",
        residue_number=5,
        atom_name="C047",
    )
    return ResolvedAttachmentPlan(
        contract=ExplicitLinkageContract(
            protein_endpoint=ReactiveEndpoint(
                participant="protein",
                selector=protein_selector,
                product_residue_name="LYX",
                leaving_atom_names=("H11",),
            ),
            modifier_endpoint=ReactiveEndpoint(
                participant="modifier",
                selector=modifier_selector,
                product_residue_name="NHX",
                leaving_atom_names=("O020",),
            ),
            bond=LinkageBond(
                protein_atom_selector=protein_selector,
                modifier_atom_selector=modifier_selector,
            ),
        ),
        protein_link_atom=PdbAtomRecord(
            serial=10,
            atom_name="NZ",
            residue_name="LYX",
            chain_id="A",
            residue_number=23,
            x=0,
            y=0,
            z=0,
        ),
        modifier_link_atom=PdbAtomRecord(
            serial=20,
            atom_name="C047",
            residue_name="NHX",
            chain_id="C",
            residue_number=5,
            x=1,
            y=0,
            z=0,
        ),
        protein_leaving_atoms=(
            PdbAtomRecord(
                serial=11,
                atom_name="H11",
                residue_name="LYS",
                chain_id="A",
                residue_number=23,
                x=0,
                y=1,
                z=0,
            ),
        ),
        modifier_leaving_atoms=(
            PdbAtomRecord(
                serial=21,
                atom_name="O020",
                residue_name="NHS",
                chain_id="C",
                residue_number=5,
                x=1,
                y=1,
                z=0,
            ),
        ),
        protein_product_residue_name="LYX",
        modifier_product_residue_name="NHX",
        pablo_crosslink_requirement=requirement,
    )


def _generic_resolved_plan(
    *,
    residue_number: int,
    modifier_residue_number: int,
    modifier_atom: str,
) -> ResolvedAttachmentPlan:
    requirement = PabloCrosslinkRequirement(
        residues=("ASX", "NAG"),
        linking_atoms=("ND2", modifier_atom),
        leaving_atoms=((), ()),
        bond_order=1,
    )
    protein_selector = PdbAtomSelector(
        chain_id="A",
        residue_name="ASN",
        residue_number=residue_number,
        atom_name="ND2",
    )
    modifier_selector = PdbAtomSelector(
        chain_id="C",
        residue_name="NAG",
        residue_number=modifier_residue_number,
        atom_name=modifier_atom,
        atom_serial=1,
        atom_index=0,
    )
    contract = ExplicitLinkageContract(
        protein_endpoint=ReactiveEndpoint(
            participant="protein",
            selector=protein_selector,
            product_residue_name="ASX",
        ),
        modifier_endpoint=ReactiveEndpoint(
            participant="modifier",
            selector=modifier_selector,
            product_residue_name="NAG",
        ),
        bond=LinkageBond(
            protein_atom_selector=protein_selector,
            modifier_atom_selector=modifier_selector,
            protein_atom_name="ND2",
            modifier_atom_name=modifier_atom,
            target_bond_length_angstrom=1.45,
        ),
        mechanism_name="n_glycosylation",
    )
    return ResolvedAttachmentPlan(
        contract=contract,
        protein_link_atom=PdbAtomRecord(
            serial=residue_number,
            atom_index=residue_number,
            atom_name="ND2",
            residue_name="ASN",
            chain_id="A",
            residue_number=residue_number,
            x=0,
            y=0,
            z=0,
            element="N",
        ),
        modifier_link_atom=PdbAtomRecord(
            serial=1,
            atom_index=0,
            atom_name=modifier_atom,
            residue_name="NAG",
            chain_id="C",
            residue_number=modifier_residue_number,
            x=1,
            y=0,
            z=0,
            element="C",
        ),
        protein_product_residue_name="ASX",
        modifier_product_residue_name="NAG",
        pablo_crosslink_requirement=requirement,
        target_bond_length_angstrom=1.45,
    )


def _generated_fragment(*, residue_name: str, residue_number: int) -> GeneratedPolymerFragment:
    atom = PolymerFragmentAtom(
        atom_index=0,
        serial=1,
        atom_name="C001",
        residue_name=residue_name,
        residue_number=residue_number,
        x=1.0,
        y=0.0,
        z=0.0,
        element="C",
    )
    residue = PolymerFragmentResidue(
        sequence_index=0,
        name=residue_name.lower(),
        residue_name=residue_name,
        residue_number=residue_number,
    )
    return GeneratedPolymerFragment(
        atoms=(atom,),
        residues=(residue,),
        reactive_atom_serial=1,
        reactive_atom_index=0,
        name=f"{residue_name.lower()}_{residue_number}",
    )


def _moiety_fragment(*, residue_name: str, residue_number: int) -> GeneratedMoietyFragment:
    generated = _generated_fragment(residue_name=residue_name, residue_number=residue_number)
    return GeneratedMoietyFragment(
        atoms=generated.atoms,
        residues=generated.residues,
        residue_name=residue_name,
        name=generated.name,
    )
