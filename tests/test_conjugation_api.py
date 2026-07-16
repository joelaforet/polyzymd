"""Tests for the Phase 1 conjugation facade skeleton."""

from __future__ import annotations

import inspect
import json
from pathlib import Path
from types import SimpleNamespace

import pytest

from polyzymd.builders.conjugation import (
    ConjugateBuildRequest,
    ConjugationEngine,
    ConjugationResult,
    build_conjugate,
    build_conjugate_from_config,
)
from polyzymd.builders.conjugation._linkage import (
    ExplicitLinkageContract,
    LinkageBond,
    PabloCrosslinkRequirement,
    PdbAtomSelector,
    ReactiveEndpoint,
    ResolvedAttachmentPlan,
)
from polyzymd.builders.conjugation.polymer import (
    GeneratedMoietyFragment,
    PolymerFragmentAtom,
    PolymerFragmentResidue,
)
from polyzymd.builders.conjugation.reactions import (
    NhsLysReaction,
    ReactionTemplate,
    get_reaction,
    list_reactions,
)
from polyzymd.config.schema import ConjugationAttachmentConfig


def test_new_conjugation_modules_import():
    """New Phase 1 modules should import without optional simulation dependencies."""
    import polyzymd.builders.conjugation.api as api
    import polyzymd.builders.conjugation.engine as engine
    import polyzymd.builders.conjugation.reactions as reactions

    assert api.build_conjugate_from_config is build_conjugate_from_config
    assert engine.ConjugationEngine is ConjugationEngine
    assert reactions.NhsLysReaction is NhsLysReaction


def test_linkage_and_pablo_product_import_without_polymer_package_cycle():
    """Linkage and Pablo modules should import without polymer package cycles."""
    import polyzymd.builders.conjugation._linkage as linkage
    from polyzymd.builders.conjugation.pablo import product

    assert linkage.PabloCrosslinkRequirement is PabloCrosslinkRequirement
    assert callable(product.build_product_state_pablo_library)


def test_reaction_library_exposes_nhs_lys_template():
    """The initial built-in reaction registry should expose NHS-Lys."""
    import polyzymd.builders.conjugation.reactions.library as library

    reactions = list_reactions()

    assert "nhs_lys" in reactions
    assert get_reaction("nhs_lys") is NhsLysReaction
    assert get_reaction("nhs_lys_amide") is NhsLysReaction
    assert issubclass(get_reaction("nhs_lys"), ReactionTemplate)
    assert NhsLysReaction.identifiers() == ("nhs_lys", "nhs_lys_amide")

    registry_source = inspect.getsource(library)
    for forbidden in ("LYX", "NHX", "NZ", "C047", "O020", "H11", "H13"):
        assert forbidden not in registry_source


def test_nhs_lys_reaction_template_settings_and_metadata_defaults():
    """NHS-Lys-specific defaults should live on the reaction template."""
    settings = NhsLysReaction.Settings()

    assert settings.product_residue_names == ("LYX", "NHX")
    assert settings.target_atom_name == "NZ"
    assert settings.bond_order == 1
    assert settings.target_bond_length_angstrom == pytest.approx(1.33)
    assert settings.max_nz_hydrogens_to_remove == 2
    assert NhsLysReaction.product_residue_names(settings) == ("LYX", "NHX")
    assert "[N:1]" in NhsLysReaction.reaction_smarts()

    role_model = NhsLysReaction.build_role_model()
    assert role_model.name == "nhs_lys_amide"
    assert [role.map_number for role in role_model.atom_roles] == [1, 2, 3, 4, 5]
    assert [participant.role for participant in role_model.participants] == ["site", "moiety"]


def test_direct_build_conjugate_is_explicitly_pending():
    """Direct non-config construction should fail with a clear engine message."""
    with pytest.raises(NotImplementedError, match="molecule/topology inputs"):
        build_conjugate(protein_mol=object(), moiety_mol=object())


def test_build_conjugate_facade_delegates_to_engine(monkeypatch, tmp_path):
    """The direct public facade should route supported inputs through the engine."""
    import polyzymd.builders.conjugation.api as api_module

    expected = ConjugationResult(output_dir=tmp_path / "out")
    calls = {}

    class FakeEngine:
        def __init__(self, *, settings=None):
            calls["settings"] = settings

        def build(self, request=None, **kwargs):
            calls["request"] = request
            calls["kwargs"] = kwargs
            return expected

    monkeypatch.setattr(api_module, "ConjugationEngine", FakeEngine)

    config_path = tmp_path / "config.yaml"
    result = build_conjugate(config_path, output_dir=tmp_path / "out", free_polymer_seed=5)

    assert result is expected
    assert calls == {
        "settings": None,
        "request": config_path,
        "kwargs": {"output_dir": tmp_path / "out", "free_polymer_seed": 5},
    }


def test_build_conjugate_from_config_facade_delegates_to_engine(monkeypatch, tmp_path):
    """The config facade should keep orchestration centralized in the engine."""
    import polyzymd.builders.conjugation.api as api_module

    expected = ConjugationResult(output_dir=tmp_path / "out")
    calls = {}

    class FakeEngine:
        def __init__(self, *, settings=None):
            calls["settings"] = settings

        def build_from_config(self, config, *, output_dir=None, free_polymer_seed=None):
            calls["config"] = config
            calls["output_dir"] = output_dir
            calls["free_polymer_seed"] = free_polymer_seed
            return expected

    monkeypatch.setattr(api_module, "ConjugationEngine", FakeEngine)

    config_path = tmp_path / "config.yaml"
    result = build_conjugate_from_config(config_path, output_dir=tmp_path / "out")

    assert result is expected
    assert calls == {
        "settings": None,
        "config": config_path,
        "output_dir": tmp_path / "out",
        "free_polymer_seed": None,
    }


def test_build_conjugate_from_config_path_delegates(monkeypatch, tmp_path):
    """Config paths should delegate to the path-based workflow."""
    import polyzymd.builders.conjugation.system_workflow as workflow_module

    expected = _workflow_result(tmp_path)
    calls = {}

    def fake_build_from_path(
        config_path, *, output_dir=None, settings=None, free_polymer_seed=None
    ):
        calls["config_path"] = config_path
        calls["output_dir"] = output_dir
        calls["settings"] = settings
        calls["free_polymer_seed"] = free_polymer_seed
        return expected

    monkeypatch.setattr(
        workflow_module,
        "build_conjugated_polymer_system_from_config_path",
        fake_build_from_path,
    )

    config_path = tmp_path / "config.yaml"
    output_dir = tmp_path / "out"
    result = build_conjugate_from_config(config_path, output_dir=output_dir, free_polymer_seed=7)

    assert isinstance(result, ConjugationResult)
    assert result.status == "completed"
    assert result.output_dir == expected.output_dir
    assert result.config_path == config_path
    assert result.workflow_json_path == expected.workflow_json_path
    assert result.artifact_paths["solvated_pdb"] == expected.solvated_pdb_path
    assert calls == {
        "config_path": config_path,
        "output_dir": output_dir,
        "settings": None,
        "free_polymer_seed": 7,
    }


def test_build_conjugate_from_config_object_delegates(monkeypatch, tmp_path):
    """In-memory configs should delegate to the object-based workflow."""
    import polyzymd.builders.conjugation.system_workflow as workflow_module

    expected = _workflow_result(tmp_path)
    config = object()
    calls = {}

    def fake_build_from_config(config_obj, *, output_dir, settings=None, free_polymer_seed=None):
        calls["config"] = config_obj
        calls["output_dir"] = output_dir
        calls["settings"] = settings
        calls["free_polymer_seed"] = free_polymer_seed
        return expected

    monkeypatch.setattr(
        workflow_module,
        "build_conjugated_polymer_system_from_config",
        fake_build_from_config,
    )

    output_dir = tmp_path / "out"
    result = build_conjugate_from_config(config, output_dir=output_dir, free_polymer_seed=11)

    assert isinstance(result, ConjugationResult)
    assert result.status == "completed"
    assert result.output_dir == expected.output_dir
    assert result.config_path is None
    assert calls == {
        "config": config,
        "output_dir": output_dir,
        "settings": None,
        "free_polymer_seed": 11,
    }


def test_engine_requires_output_dir_for_in_memory_config():
    """The delegated object workflow still requires an output directory."""
    with pytest.raises(ValueError, match="output_dir is required"):
        ConjugationEngine().build_from_config(object())


def test_engine_build_from_request_delegates_to_config_workflow(monkeypatch, tmp_path):
    """Engine request inputs should use the current config-driven workflow."""
    import polyzymd.builders.conjugation.system_workflow as workflow_module

    expected = _workflow_result(tmp_path)
    calls = {}

    def fake_build_from_path(
        config_path, *, output_dir=None, settings=None, free_polymer_seed=None
    ):
        calls["config_path"] = config_path
        calls["output_dir"] = output_dir
        calls["settings"] = settings
        calls["free_polymer_seed"] = free_polymer_seed
        return expected

    monkeypatch.setattr(
        workflow_module,
        "build_conjugated_polymer_system_from_config_path",
        fake_build_from_path,
    )

    request = ConjugateBuildRequest(
        config_path=tmp_path / "config.yaml",
        output_dir=tmp_path / "out",
        free_polymer_seed=17,
    )

    result = ConjugationEngine().build(request)

    assert result.output_dir == expected.output_dir
    assert calls == {
        "config_path": tmp_path / "config.yaml",
        "output_dir": tmp_path / "out",
        "settings": None,
        "free_polymer_seed": 17,
    }


def test_direct_request_builds_two_smiles_n_glycosylation_plans_once(
    monkeypatch,
    tmp_path,
):
    """Public direct requests should plan all SMILES moieties before one construction run."""
    import polyzymd.builders.conjugation.system_workflow as workflow_module

    protein_path = tmp_path / "protein.pdb"
    protein_path.write_text("END\n", encoding="utf-8")
    attachments = (
        _n_gly_attachment("glycan_1", residue_number=42),
        _n_gly_attachment("glycan_2", residue_number=87),
    )
    calls = {"moieties": [], "plans": [], "construct": 0, "solvate": 0}

    def fake_resolve_source(
        attachment,
        *,
        attachment_index,
        output_dir,
        protein_pdb_path,
        random_seed=None,
        **kwargs,
    ):
        calls["moieties"].append(
            (
                attachment.moiety.smiles,
                attachment.moiety.residue_name,
                attachment.moiety.name,
                Path(output_dir),
                random_seed,
            )
        )
        fragment = _moiety_fragment(
            attachment.moiety.residue_name,
            name=attachment.moiety.name,
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

    class FakeReaction:
        coordinate_backend_mechanism = "n_glycosylation"

        @staticmethod
        def settings_from_attachment(attachment):
            return SimpleNamespace(target_atom_name="ND2")

        @staticmethod
        def resolve_plan(protein, site, fragment, *, settings=None):
            calls["plans"].append((Path(protein), site.residue_number, site.atom_name, settings))
            return _resolved_plan(
                residue_number=site.residue_number,
                modifier_residue_name=fragment.residue_name,
                modifier_atom_name="C001",
            )

    def fake_get_reaction(name):
        assert name == "n_glycosylation"
        return FakeReaction

    def fake_construct(**kwargs):
        calls["construct"] += 1
        assert len(kwargs["specs"]) == 2
        assert len([spec.resolved_plan for spec in kwargs["specs"]]) == 2
        crosslinked = Path(kwargs["output_dir"]) / "assembled_crosslinked.pdb"
        relaxed = Path(kwargs["output_dir"]) / "conjugate_relaxed.pdb"
        crosslinked.write_text("END\n", encoding="utf-8")
        relaxed.write_text("END\n", encoding="utf-8")
        construction = SimpleNamespace(
            crosslinked_pdb_path=crosslinked,
            relaxation=SimpleNamespace(relaxed_pdb_path=relaxed),
        )
        return construction, object()

    class FakeSolvatedBuilder:
        interchange = None
        solvated_topology = object()

        def save_topology(self, path):
            Path(path).write_text("END\n", encoding="utf-8")

    def fake_solvate(**kwargs):
        calls["solvate"] += 1
        assert kwargs["working_dir"] == tmp_path / "out"
        assert kwargs["create_interchange"] is False
        return FakeSolvatedBuilder()

    monkeypatch.setattr(workflow_module, "resolve_moiety_source", fake_resolve_source)
    monkeypatch.setattr(workflow_module, "get_reaction", fake_get_reaction)
    monkeypatch.setattr(
        workflow_module,
        "_prepared_protein_pdb_path",
        lambda path, *, output_dir, settings: (Path(path), None),
    )
    monkeypatch.setattr(workflow_module, "_construct_conjugate_from_specs", fake_construct)
    monkeypatch.setattr(
        workflow_module, "topology_with_pdb_positions", lambda topology, path: topology
    )
    monkeypatch.setattr(workflow_module, "_build_direct_solvated_system", fake_solvate)
    monkeypatch.setattr(
        workflow_module, "_restore_pdb_atom_name_fields", lambda target, template: 0
    )

    result = ConjugationEngine().build(
        ConjugateBuildRequest(
            protein_pdb_path=protein_path,
            attachments=attachments,
            output_dir=tmp_path / "out",
            free_polymer_seed=19,
        )
    )

    assert result.crosslinked_conjugate_pdb_path == (
        tmp_path / "out" / "conjugate-construction" / "assembled_crosslinked.pdb"
    )
    assert result.relaxed_conjugate_pdb_path == (
        tmp_path / "out" / "conjugate-construction" / "conjugate_relaxed.pdb"
    )
    assert result.solvated_pdb_path == tmp_path / "out" / "solvated_conjugate_free_polymers.pdb"
    assert [entry[1] for entry in calls["moieties"]] == ["NAG", "NAG"]
    assert [entry[2] for entry in calls["plans"]] == [None, None]
    assert calls["construct"] == 1
    assert calls["solvate"] == 1


def test_smiles_moiety_residue_name_validation_rejects_non_one_residue_code(tmp_path):
    """SMILES moieties should require a three-character PDB-safe residue code."""
    with pytest.raises(ValueError, match="SMILES moiety residue names"):
        ConjugateBuildRequest(
            protein_pdb_path=tmp_path / "protein.pdb",
            output_dir=tmp_path / "out",
            attachments=(
                {
                    "name": "bad_glycan",
                    "site": {"chain_id": "A", "residue_name": "ASN", "residue_number": 42},
                    "moiety": {"name": "glycan", "smiles": "CO", "residue_name": "GLCN"},
                    "mechanism": {"name": "n_glycosylation"},
                },
            ),
        )


def test_conjugation_result_collects_workflow_output_paths(tmp_path):
    """The public result should carry status and useful workflow artifact paths."""
    workflow_result = _workflow_result(tmp_path)

    result = ConjugationResult.from_workflow_result(
        workflow_result,
        config_path=tmp_path / "config.yaml",
    )

    assert result.status == "completed"
    assert result.output_dir == tmp_path / "workflow-out"
    assert result.config_path == tmp_path / "config.yaml"
    assert result.crosslinked_conjugate_pdb_path == tmp_path / "crosslinked.pdb"
    assert result.relaxed_conjugate_pdb_path == tmp_path / "relaxed.pdb"
    assert result.solvated_pdb_path == tmp_path / "solvated.pdb"
    assert result.workflow_json_path == tmp_path / "workflow.json"
    assert result.final_interchange_created is False
    assert result.artifact_paths == {
        "crosslinked_conjugate_pdb": tmp_path / "crosslinked.pdb",
        "relaxed_conjugate_pdb": tmp_path / "relaxed.pdb",
        "solvated_pdb": tmp_path / "solvated.pdb",
        "workflow_json": tmp_path / "workflow.json",
    }


def test_conjugation_result_requires_final_interchange():
    """Final Interchange access should fail clearly when export state is unavailable."""
    result = ConjugationResult()

    with pytest.raises(RuntimeError, match="final Interchange"):
        result.require_final_interchange()

    interchange = object()
    assert (
        ConjugationResult(final_interchange=interchange).require_final_interchange() is interchange
    )


def test_conjugation_result_get_component_info_delegates():
    """Component metadata should come from the retained system builder."""
    info = {"protein": 1}
    builder = SimpleNamespace(get_component_info=lambda: info)

    assert ConjugationResult(system_builder=builder).get_component_info() == info

    with pytest.raises(RuntimeError, match="system builder"):
        ConjugationResult().get_component_info()


def test_conjugation_result_serialization_excludes_heavy_fields(tmp_path):
    """Serialized result payloads should omit in-memory topology and builder objects."""
    result = ConjugationResult(
        output_dir=tmp_path / "out",
        final_interchange=object(),
        system_builder=object(),
        relaxed_conjugate_topology=object(),
        solvated_topology=object(),
        construction=object(),
    )

    payload = result.model_dump()
    output_path = result.save(tmp_path / "result.json")
    saved_payload = json.loads(output_path.read_text(encoding="utf-8"))

    for field in (
        "final_interchange",
        "system_builder",
        "relaxed_conjugate_topology",
        "solvated_topology",
        "construction",
    ):
        assert field not in payload
        assert field not in saved_payload


def test_package_facade_exports_public_api_only():
    """The package facade should expose only the stable public API."""
    import polyzymd.builders.conjugation as conjugation

    expected = {
        "ConjugateBuildRequest",
        "ConjugationResult",
        "ConjugationEngine",
        "build_conjugate",
        "build_conjugate_from_config",
    }

    assert set(conjugation.__all__) == expected
    for name in expected:
        assert hasattr(conjugation, name)

    assert not hasattr(conjugation, "build_conjugated_polymer_system_from_config")
    assert not hasattr(conjugation, "list_reactions")


def _workflow_result(tmp_path: Path) -> ConjugationResult:
    return ConjugationResult(
        output_dir=tmp_path / "workflow-out",
        construction=SimpleNamespace(
            crosslinked_pdb_path=tmp_path / "crosslinked.pdb",
            relaxation=SimpleNamespace(relaxed_pdb_path=tmp_path / "relaxed.pdb"),
        ),
        relaxed_conjugate_pdb_path=tmp_path / "relaxed.pdb",
        solvated_pdb_path=tmp_path / "solvated.pdb",
        workflow_json_path=tmp_path / "workflow.json",
        final_interchange_created=False,
    )


def _n_gly_attachment(name: str, *, residue_number: int) -> ConjugationAttachmentConfig:
    return ConjugationAttachmentConfig.model_validate(
        {
            "name": name,
            "site": {"chain_id": "A", "residue_name": "ASN", "residue_number": residue_number},
            "moiety": {
                "name": name,
                "smiles": "CO",
                "residue_name": "NAG",
            },
            "mechanism": {"name": "n_glycosylation"},
        }
    )


def _moiety_fragment(residue_name: str, *, name: str) -> GeneratedMoietyFragment:
    atom = PolymerFragmentAtom(
        atom_index=0,
        serial=1,
        atom_name="C001",
        residue_name=residue_name,
        residue_number=1,
        x=0.0,
        y=0.0,
        z=0.0,
        element="C",
    )
    return GeneratedMoietyFragment(
        atoms=(atom,),
        bonds=(),
        bond_orders=(),
        residues=(
            PolymerFragmentResidue(
                sequence_index=0,
                name=name,
                residue_name=residue_name,
                residue_number=1,
            ),
        ),
        residue_name=residue_name,
        name=name,
    )


def _resolved_plan(
    *,
    residue_number: int,
    modifier_residue_name: str,
    modifier_atom_name: str,
) -> ResolvedAttachmentPlan:
    protein_selector = PdbAtomSelector(
        chain_id="A",
        residue_name="ASN",
        residue_number=residue_number,
        atom_name="ND2",
    )
    modifier_selector = PdbAtomSelector(
        chain_id="C",
        residue_name=modifier_residue_name,
        residue_number=1,
        atom_name=modifier_atom_name,
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
            product_residue_name=modifier_residue_name,
        ),
        bond=LinkageBond(
            protein_atom_name="ND2",
            modifier_atom_name=modifier_atom_name,
            target_bond_length_angstrom=1.45,
        ),
        mechanism_name="n_glycosylation",
    )
    protein_atom = _pdb_atom_record(
        serial=10 + residue_number,
        atom_index=residue_number,
        atom_name="ND2",
        residue_name="ASN",
        chain_id="A",
        residue_number=residue_number,
        element="N",
    )
    modifier_atom = _pdb_atom_record(
        serial=1,
        atom_index=0,
        atom_name=modifier_atom_name,
        residue_name=modifier_residue_name,
        chain_id="C",
        residue_number=1,
        element="C",
    )
    requirement = PabloCrosslinkRequirement(
        residues=("ASX", modifier_residue_name),
        linking_atoms=("ND2", modifier_atom_name),
        leaving_atoms=((), ()),
        bond_order=1,
    )
    return ResolvedAttachmentPlan(
        contract=contract,
        protein_link_atom=protein_atom,
        modifier_link_atom=modifier_atom,
        protein_product_residue_name="ASX",
        modifier_product_residue_name=modifier_residue_name,
        pablo_crosslink_requirement=requirement,
        target_bond_length_angstrom=1.45,
    )


def _pdb_atom_record(**kwargs) -> object:
    from polyzymd.builders.conjugation.structure.pdb import PdbAtomRecord

    defaults = {"x": 0.0, "y": 0.0, "z": 0.0}
    defaults.update(kwargs)
    return PdbAtomRecord(**defaults)
