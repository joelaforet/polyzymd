"""Tests for config-driven conjugation system workflow helpers."""

from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace

import pytest

from polyzymd.builders.conjugation._assembly import (
    ModifierConstructionSettings,
    PackmolModifierPlacementResult,
)
from polyzymd.builders.conjugation._linkage import (
    ExplicitLinkageContract,
    LinkageBond,
    PabloCrosslinkRequirement,
    PdbAtomSelector,
    ReactiveEndpoint,
    ResolvedAttachmentPlan,
)
from polyzymd.builders.conjugation.pablo.ingestion import PabloAvailability, PabloIngestionResult
from polyzymd.builders.conjugation.pablo.parameterization import InterchangeParameterizationResult
from polyzymd.builders.conjugation.reactions.nhs_lys import NhsLysReaction
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
    _construct_nhs_lys_modifier_linked_protein,
    _local_minimization_settings_for_product,
    _nhs_lys_linker_from_attachment,
    _policy_with_resolved_crosslink,
    _prepared_protein_pdb_path,
    _require_supported_coordinate_backend,
    _restore_pdb_atom_name_fields,
)
from polyzymd.config.schema import ConjugationCcdPabloPolicyConfig


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
    """Public workflow defaults should mirror the known-working product-state POC path."""
    import polyzymd.builders.conjugation.system_workflow as workflow_module

    settings = workflow_module.ConjugatedPolymerSystemSettings()

    assert settings.canonicalize_source_protein_hydrogens is True
    assert settings.use_product_state_pablo_library is True
    assert settings.run_product_state_local_minimization is True
    assert settings.protein_canonicalization.ph == pytest.approx(7.0)


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
    """Local minimization selectors should use product atoms, not leaving atoms."""
    from polyzymd.builders.conjugation._relaxation import LocalMinimizationSettings

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
        requirement=requirement,
    )

    assert settings.nz_selector.serial is None
    assert settings.nz_selector.chain_id == "A"
    assert settings.nz_selector.residue_name == "LYX"
    assert settings.nz_selector.residue_number == 23
    assert settings.c047_selector.chain_id == "C"
    assert settings.c047_selector.residue_number == 5
    assert settings.o020_selector.chain_id == "C"
    assert settings.o020_selector.residue_number == 5
    assert settings.o020_selector.atom_name == "OX1"
    assert settings.o020_selector.atom_name != requirement.leaving_atoms[1][0]


def test_local_minimization_settings_preserve_explicit_selectors(tmp_path: Path):
    """Explicit local minimization selectors should override product-state inference."""
    from polyzymd.builders.conjugation._relaxation import (
        CrosslinkAtomSelector,
        LocalMinimizationSettings,
    )

    settings = LocalMinimizationSettings(
        nz_selector=CrosslinkAtomSelector(
            chain_id="A",
            residue_name="LYX",
            residue_number=23,
            atom_name="NZ",
        ),
        c047_selector=CrosslinkAtomSelector(
            chain_id="C",
            residue_name="NHX",
            residue_number=5,
            atom_name="C047",
        ),
        o020_selector=CrosslinkAtomSelector(
            chain_id="C",
            residue_name="NHX",
            residue_number=5,
            atom_name="USER",
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
        requirement=requirement,
    )

    assert derived is settings
    assert derived.o020_selector.atom_name == "USER"


def test_local_minimization_settings_can_use_product_residue_definition(tmp_path: Path):
    """Product residue definitions can identify the carbonyl oxygen when CONECT is absent."""
    from polyzymd.builders.conjugation._relaxation import LocalMinimizationSettings

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
    product_library = SimpleNamespace(
        definitions=(
            SimpleNamespace(
                residue_name="NHX",
                atoms=(
                    SimpleNamespace(name="C047", symbol="C"),
                    SimpleNamespace(name="ODEF", symbol="O"),
                ),
                bonds=(SimpleNamespace(atom1="C047", atom2="ODEF"),),
            ),
        )
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
        requirement=requirement,
        product_state_pablo_library=product_library,
    )

    assert settings.o020_selector.atom_name == "ODEF"


def test_construct_uses_product_state_library_and_local_minimization(
    monkeypatch,
    tmp_path: Path,
):
    """NHS-Lys public construction should pass product-state definitions to Pablo."""
    import polyzymd.builders.conjugation._relaxation as local_min_module
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

    def fake_place(*args, **kwargs):
        calls["placed_protein_path"] = Path(args[0])
        return placement

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

    def fake_parameterize(topology, *, settings=None, charge_from_molecules=None):
        calls["charge_template_count"] = len(charge_from_molecules or ())
        return InterchangeParameterizationResult(
            success=True,
            interchange=object(),
            force_field_names=("fake.offxml",),
            topology_type="FakeTopology",
        )

    def fake_local_minimize(pdb_path, output_dir=None, **kwargs):
        calls["local_product_library"] = kwargs["product_state_pablo_library"]
        calls["local_requirement"] = kwargs["pablo_crosslink_requirement"]
        relaxed = Path(output_dir) / "relaxed.pdb"
        relaxed.write_text(Path(pdb_path).read_text(encoding="utf-8"), encoding="utf-8")
        return SimpleNamespace(success=True, relaxed_pdb_path=relaxed, blocker=None)

    monkeypatch.setattr(workflow_module, "place_modifier_with_packmol", fake_place)
    monkeypatch.setattr(
        workflow_module, "placed_fragment_from_resolved_plan", lambda fragment, plan: fragment
    )
    monkeypatch.setattr(workflow_module, "write_crosslinked_pdb", fake_write)
    monkeypatch.setattr(
        product_pablo_module, "build_product_state_pablo_library", fake_product_library
    )
    monkeypatch.setattr(workflow_module, "PabloIngestor", FakeIngestor)
    monkeypatch.setattr(
        workflow_module, "build_formal_charge_smoke_template", lambda molecule: object()
    )
    monkeypatch.setattr(
        workflow_module, "create_interchange_from_pablo_topology", fake_parameterize
    )
    monkeypatch.setattr(
        local_min_module, "run_post_crosslink_local_minimization", fake_local_minimize
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

    construction, topology = _construct_nhs_lys_modifier_linked_protein(
        protein_pdb_path=tmp_path / "raw.pdb",
        prepared_protein_pdb_path=prepared_protein,
        modifier=object(),
        polymer_sdf_path=polymer_sdf,
        linker=object(),
        resolved_plan=resolved_plan,
        ccd_pablo_policy=policy,
        chain_policy=None,
        output_dir=tmp_path / "construction",
        settings=ModifierConstructionSettings(run_smoke=True),
        use_product_state_pablo_library=True,
        run_product_state_local_minimization=True,
        local_minimization_settings=workflow_module._default_local_minimization_settings(),
    )

    assert topology is not None
    assert construction.local_minimization is not None
    assert construction.smoke is None
    assert calls["placed_protein_path"] == prepared_protein
    assert calls["assembly_protein_path"] == prepared_protein
    assert calls["source_protein_pdb"] == prepared_protein
    assert calls["polymer_sdf"] == polymer_sdf
    assert calls["residue_library"] is product_library.residue_library
    assert calls["local_product_library"] is product_library
    assert calls["local_requirement"].leaving_atoms == ((), ())
    assert calls["charge_template_count"] == 1


def test_coordinate_backend_gate_allows_explicit_nhs_lys_mechanism():
    """The system workflow should only route the named implemented backend to NHS-Lys."""
    attachment = SimpleNamespace(
        mechanism=SimpleNamespace(name="nhs_lys_amide", reaction_smarts=None)
    )

    _require_supported_coordinate_backend(attachment)


def test_nhs_lys_workflow_linker_uses_reaction_template_defaults():
    """Workflow linker defaults should come from the NHS-Lys reaction template."""
    attachment = SimpleNamespace(
        site=SimpleNamespace(
            chain_id="A",
            residue_name=None,
            residue_number=23,
            insertion_code="",
            atom_name=None,
        ),
        mechanism=SimpleNamespace(
            product_residues=SimpleNamespace(site=None, moiety=None),
            bond=SimpleNamespace(site_atom=None, order=1, target_bond_length_angstrom=1.33),
        ),
    )

    linker = _nhs_lys_linker_from_attachment(attachment)
    defaults = NhsLysReaction.Settings()

    assert linker.target_residue_name == defaults.source_site_residue_name
    assert linker.target_atom_name == defaults.target_atom_name
    assert linker.lysine_target_resname == defaults.product_site_residue_name
    assert linker.modifier_target_resname == defaults.product_moiety_residue_name


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
