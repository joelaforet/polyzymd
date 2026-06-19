"""Tests for integrated modifier-linking construction orchestration."""

from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace

import pytest

from polyzymd.builders.conjugation._assembly import (
    construct_explicit_pdb_linkage,
    construct_modifier_linked_protein,
)
from polyzymd.builders.conjugation._linkage import (
    ExplicitLinkageContract,
    LinkageBond,
    MissingPabloCrosslinkError,
    NhsLysModifierLinker,
    PdbAtomSelector,
    ReactiveEndpoint,
)
from polyzymd.builders.conjugation._relaxation import VacuumSmokeResult
from polyzymd.builders.conjugation.builder import CovalentModificationBuilder
from polyzymd.builders.conjugation.exceptions import ConjugationNotImplementedError
from polyzymd.builders.conjugation.pablo.ingestion import PabloAvailability, PabloIngestionResult
from polyzymd.builders.conjugation.pablo.parameterization import InterchangeParameterizationResult
from polyzymd.builders.conjugation.polymer import GeneratedPolymerFragment
from polyzymd.builders.conjugation.structure.pdb import PdbAtomRecord
from polyzymd.config.schema import (
    ConjugationAttachmentConfig,
    ConjugationConfig,
    ConjugationMechanismConfig,
    ConjugationMode,
    ConjugationMoietyConfig,
    ConjugationSiteConfig,
)


def test_missing_crosslink_config_fails_before_packmol(tmp_path: Path):
    """Construction should require explicit LYX/NHX Pablo crosslinks before placement."""
    protein_path = _protein_pdb(tmp_path)
    modifier = _generated_modifier()
    linker = NhsLysModifierLinker(target_residue_number=23)

    def forbidden_packmol(*args, **kwargs):
        """Fail if placement is attempted before crosslink validation."""
        raise AssertionError("Packmol should not run when ccd_pablo.crosslinks is missing")

    with pytest.raises(MissingPabloCrosslinkError, match="ccd_pablo.crosslinks"):
        construct_modifier_linked_protein(
            protein_pdb_path=protein_path,
            modifier=modifier,
            linker=linker,
            ccd_pablo_policy=SimpleNamespace(crosslinks=[]),
            output_dir=tmp_path / "out",
            run_packmol_func=forbidden_packmol,
        )


def test_integrated_helper_orchestrates_placement_assembly_pablo_parameterization_smoke(
    tmp_path: Path,
):
    """Integrated helper should produce expected artifacts using mocked heavy boundaries."""
    protein_path = _protein_pdb(tmp_path)
    modifier = _generated_modifier()
    linker = NhsLysModifierLinker(target_residue_number=23)
    policy = SimpleNamespace(
        crosslinks=[
            SimpleNamespace(
                residues=("LYX", "NHX"),
                linking_atoms=("NZ", "RC"),
                leaving_atoms=(("HZ2", "HZ3"), ("LG",)),
                bond_order=1,
            )
        ]
    )
    calls: list[str] = []

    def fake_run_packmol(input_text: str, work_dir: Path) -> Path:
        """Write a Packmol-like concatenated output."""
        calls.append("packmol")
        output_path = work_dir / "packmol_output.pdb"
        protein_lines = [
            line
            for line in (work_dir / "protein_fixed_sterics.pdb").read_text().splitlines(True)
            if line.startswith(("ATOM", "HETATM"))
        ]
        modifier_lines = [
            line
            for line in (work_dir / "modifier_retained.pdb").read_text().splitlines(True)
            if line.startswith(("ATOM", "HETATM"))
        ]
        output_path.write_text("".join([*protein_lines, *modifier_lines, "END\n"]))
        assert "movebadrandom" in input_text
        return output_path

    topology = object()

    class FakePabloIngestor:
        """Small Pablo ingestion fake."""

        def ingest_structure(self, path, *, chain_policy=None, output_dir=None):
            """Return a successful Pablo result."""
            calls.append("pablo")
            assert Path(path).name == "assembled_crosslinked.pdb"
            return PabloIngestionResult(
                success=True,
                path=Path(path),
                suffix=".pdb",
                pablo=PabloAvailability(available=True),
                topology=topology,
            )

    interchange = object()

    def fake_parameterizer(observed_topology, *, settings):
        """Return an injected Interchange result."""
        calls.append("parameterizer")
        assert observed_topology is topology
        return InterchangeParameterizationResult(
            success=True,
            interchange=interchange,
            force_field_names=settings.force_field_names,
            topology_type="FakeTopology",
        )

    def fake_smoke_runner(observed_interchange, output_dir, *, settings):
        """Return an injected smoke result and write its JSON artifact."""
        calls.append("smoke")
        assert observed_interchange is interchange
        smoke_json = Path(output_dir) / settings.smoke_json_name
        smoke_json.write_text("{}\n")
        return VacuumSmokeResult(
            success=True,
            output_dir=Path(output_dir),
            smoke_json_path=smoke_json,
            platform_name="Mock",
            restrained_atom_count=3,
            energy_before_min_kj_mol=1.0,
            energy_after_min_kj_mol=0.5,
            energy_before_nvt_kj_mol=0.5,
            energy_after_nvt_kj_mol=0.4,
            nvt_steps=10,
        )

    result = construct_modifier_linked_protein(
        protein_pdb_path=protein_path,
        modifier=modifier,
        linker=linker,
        ccd_pablo_policy=policy,
        output_dir=tmp_path / "out",
        run_packmol_func=fake_run_packmol,
        pablo_ingestor=FakePabloIngestor(),
        parameterizer=fake_parameterizer,
        smoke_runner=fake_smoke_runner,
    )

    assert calls == ["packmol", "pablo", "parameterizer", "smoke"]
    assert result.crosslinked_pdb_path.exists()
    assert result.crosslink_validation.residues == ("LYX", "NHX")
    assert result.placement.packmol_input_path.exists()
    assert result.assembly.added_conect_pair[0] != result.assembly.added_conect_pair[1]
    assert result.pablo.topology is topology
    assert result.parameterization.interchange is interchange
    assert result.smoke is not None
    assert result.smoke.smoke_json_path.exists()


def test_explicit_pdb_linkage_orchestrates_generic_modifier_pdb(tmp_path: Path):
    """Generic construction should consume a user-supplied multi-residue modifier PDB."""
    protein_path = _protein_pdb(tmp_path)
    modifier_path = _modifier_pdb(tmp_path)
    contract = _explicit_contract()
    policy = SimpleNamespace(
        crosslinks=[
            SimpleNamespace(
                residues=("LYX", "NHX"),
                linking_atoms=("NZ", "RC"),
                leaving_atoms=(("HZ2", "HZ3"), ("LG",)),
                bond_order=1,
            )
        ]
    )
    calls: list[str] = []

    def fake_run_packmol(input_text: str, work_dir: Path) -> Path:
        """Write a Packmol-like concatenated output for the generic path."""
        calls.append("packmol")
        output_path = work_dir / "packmol_output.pdb"
        protein_lines = [
            line
            for line in (work_dir / "protein_fixed_sterics.pdb").read_text().splitlines(True)
            if line.startswith(("ATOM", "HETATM"))
        ]
        modifier_lines = [
            line
            for line in (work_dir / "modifier_retained.pdb").read_text().splitlines(True)
            if line.startswith(("ATOM", "HETATM"))
        ]
        output_path.write_text("".join([*protein_lines, *modifier_lines, "END\n"]))
        assert "atoms 2" in input_text
        return output_path

    topology = object()

    class FakePabloIngestor:
        """Small Pablo ingestion fake for generic explicit construction."""

        def ingest_structure(self, path, *, chain_policy=None, output_dir=None):
            """Return a successful Pablo result."""
            calls.append("pablo")
            return PabloIngestionResult(
                success=True,
                path=Path(path),
                suffix=".pdb",
                pablo=PabloAvailability(available=True),
                topology=topology,
            )

    interchange = object()

    def fake_parameterizer(observed_topology, *, settings):
        """Return an injected Interchange result."""
        calls.append("parameterizer")
        assert observed_topology is topology
        return InterchangeParameterizationResult(
            success=True,
            interchange=interchange,
            force_field_names=settings.force_field_names,
            topology_type="FakeTopology",
        )

    def fake_smoke_runner(observed_interchange, output_dir, *, settings):
        """Return an injected smoke result."""
        calls.append("smoke")
        smoke_json = Path(output_dir) / settings.smoke_json_name
        smoke_json.write_text("{}\n")
        return VacuumSmokeResult(
            success=True,
            output_dir=Path(output_dir),
            smoke_json_path=smoke_json,
            platform_name="Mock",
            restrained_atom_count=3,
            energy_before_min_kj_mol=1.0,
            energy_after_min_kj_mol=0.5,
            energy_before_nvt_kj_mol=0.5,
            energy_after_nvt_kj_mol=0.4,
            nvt_steps=10,
        )

    result = construct_explicit_pdb_linkage(
        protein_pdb_path=protein_path,
        modifier_pdb_path=modifier_path,
        contract=contract,
        ccd_pablo_policy=policy,
        output_dir=tmp_path / "generic-out",
        run_packmol_func=fake_run_packmol,
        pablo_ingestor=FakePabloIngestor(),
        parameterizer=fake_parameterizer,
        smoke_runner=fake_smoke_runner,
    )

    assert calls == ["packmol", "pablo", "parameterizer", "smoke"]
    assert result.resolved_plan.modifier_link_atom.residue_number == 2
    assert result.crosslink_validation.linking_atoms == ("NZ", "RC")
    assert result.assembly.removed_atom_names == ("HZ2", "HZ3", "LG")


def test_explicit_pdb_linkage_missing_crosslink_fails_before_packmol(tmp_path: Path):
    """Generic construction should validate Pablo crosslinks before Packmol."""
    protein_path = _protein_pdb(tmp_path)
    modifier_path = _modifier_pdb(tmp_path)

    def forbidden_packmol(*args, **kwargs):
        """Fail if placement runs before crosslink validation."""
        raise AssertionError("Packmol should not run when crosslink config is missing")

    with pytest.raises(MissingPabloCrosslinkError, match="ccd_pablo.crosslinks"):
        construct_explicit_pdb_linkage(
            protein_pdb_path=protein_path,
            modifier_pdb_path=modifier_path,
            contract=_explicit_contract(),
            ccd_pablo_policy=SimpleNamespace(crosslinks=[]),
            output_dir=tmp_path / "generic-out",
            run_packmol_func=forbidden_packmol,
        )


def test_builder_does_not_use_v1_direct_bridge_as_core_acceptance(tmp_path: Path):
    """The builder should keep direct OpenFF handoff outside the core path."""
    protein_path = _protein_pdb(tmp_path)
    modifier = _generated_modifier()
    plan = NhsLysModifierLinker(target_residue_number=23).resolve_plan(protein_path, modifier)
    topology = object()

    def fake_direct_bridge(**kwargs):
        """Fail if the quarantined direct bridge is invoked."""
        raise AssertionError("Direct OpenFF bridge should remain quarantined")

    config = ConjugationConfig(
        enabled=True,
        mode=ConjugationMode.CONSTRUCT,
        attachments=[
            ConjugationAttachmentConfig(
                name="lys23-acb",
                site=ConjugationSiteConfig(
                    chain_id="A",
                    residue_name="LYS",
                    residue_number=23,
                    atom_name="NZ",
                ),
                moiety=ConjugationMoietyConfig(
                    name="NHS-linker",
                    role="moiety",
                    smiles="CC(=O)ON1C(=O)CCC1=O",
                ),
                mechanism=ConjugationMechanismConfig(name="nhs_lys_amide"),
            )
        ],
    )
    builder = CovalentModificationBuilder(config, output_dir=tmp_path / "direct")

    with pytest.raises(ConjugationNotImplementedError, match="covalent graph surgery"):
        builder.build(
            topology,
            context={
                "protein_pdb_path": protein_path,
                "conjugation_placed_modifier": modifier.to_placed_fragment(),
                "conjugation_resolved_plan": plan,
                "conjugation_direct_bridge": fake_direct_bridge,
            },
        )


def _protein_pdb(tmp_path: Path) -> Path:
    """Create a small lysine-containing protein PDB."""
    path = tmp_path / "protein.pdb"
    path.write_text(
        _pdb_atom(1, "N", "LYS", "A", 23, 0.0, 0.0, 0.0, element="N")
        + _pdb_atom(2, "CA", "LYS", "A", 23, 1.0, 0.0, 0.0)
        + _pdb_atom(3, "CE", "LYS", "A", 23, 1.5, 0.0, 0.0)
        + _pdb_atom(4, "NZ", "LYS", "A", 23, 2.0, 0.0, 0.0, element="N")
        + _pdb_atom(5, "HZ1", "LYS", "A", 23, 2.0, 0.7, 0.0, element="H")
        + _pdb_atom(6, "HZ2", "LYS", "A", 23, 2.0, -0.7, 0.0, element="H")
        + _pdb_atom(7, "HZ3", "LYS", "A", 23, 2.0, 0.0, 0.7, element="H")
        + _pdb_atom(8, "N", "ALA", "A", 24, 4.0, 0.0, 0.0, element="N")
        + "END\n"
    )
    return path


def _modifier_pdb(tmp_path: Path) -> Path:
    """Create a user-supplied multi-residue modifier PDB."""
    path = tmp_path / "modifier.pdb"
    path.write_text(
        _pdb_atom(101, "C1", "SBM", "C", 1, 5.0, 0.0, 0.0, record="HETATM")
        + _pdb_atom(102, "RC", "NHS", "C", 2, 3.3, 0.0, 0.0, record="HETATM")
        + _pdb_atom(103, "O1", "NHS", "C", 2, 3.8, 0.5, 0.0, element="O", record="HETATM")
        + _pdb_atom(104, "LG", "NHS", "C", 2, 4.2, 1.0, 0.0, element="O", record="HETATM")
        + _pdb_atom(105, "LG", "EGP", "C", 3, 5.2, 1.0, 0.0, element="O", record="HETATM")
        + "CONECT  101  102\n"
        + "CONECT  102  101  103  104\n"
        + "CONECT  103  102\n"
        + "CONECT  104  102\n"
        + "END\n"
    )
    return path


def _explicit_contract() -> ExplicitLinkageContract:
    """Build a generic explicit linkage contract for integrated tests."""
    protein_selector = PdbAtomSelector(
        chain_id="A",
        residue_name="LYS",
        residue_number=23,
        atom_name="NZ",
    )
    modifier_selector = PdbAtomSelector(
        chain_id="C",
        residue_name="NHS",
        residue_number=2,
        atom_name="RC",
    )
    return ExplicitLinkageContract(
        protein_endpoint=ReactiveEndpoint(
            participant="protein",
            selector=protein_selector,
            product_residue_name="LYX",
            leaving_atom_names=("HZ2", "HZ3"),
        ),
        modifier_endpoint=ReactiveEndpoint(
            participant="modifier",
            selector=modifier_selector,
            product_residue_name="NHX",
            leaving_atom_names=("LG",),
        ),
        bond=LinkageBond(protein_atom_name="NZ", modifier_atom_name="RC", bond_order=1),
        mechanism_name="explicit_linkage",
    )


def _generated_modifier() -> GeneratedPolymerFragment:
    """Create a small generated modifier with a leaving group."""
    atoms = (
        PdbAtomRecord(
            serial=101,
            atom_index=0,
            atom_name="C1",
            residue_name="SB1",
            chain_id="Z",
            residue_number=1,
            x=5.0,
            y=0.0,
            z=0.0,
            element="C",
            record_name="HETATM",
        ),
        PdbAtomRecord(
            serial=102,
            atom_index=1,
            atom_name="RC",
            residue_name="NHS",
            chain_id="Z",
            residue_number=2,
            x=3.3,
            y=0.0,
            z=0.0,
            element="C",
            record_name="HETATM",
        ),
        PdbAtomRecord(
            serial=103,
            atom_index=2,
            atom_name="O1",
            residue_name="NHS",
            chain_id="Z",
            residue_number=2,
            x=3.8,
            y=0.5,
            z=0.0,
            element="O",
            record_name="HETATM",
        ),
        PdbAtomRecord(
            serial=104,
            atom_index=3,
            atom_name="LG",
            residue_name="NHS",
            chain_id="Z",
            residue_number=2,
            x=4.2,
            y=1.0,
            z=0.0,
            element="O",
            record_name="HETATM",
        ),
    )
    return GeneratedPolymerFragment.from_atom_records(
        atoms,
        bonds=((101, 102), (102, 103), (102, 104)),
        reactive_atom_serial=102,
        reactive_atom_index=1,
        reactive_atom_name="RC",
        leaving_atom_serials=(104,),
        leaving_atom_indices=(3,),
        leaving_atom_names=("LG",),
        name="modifier",
    )


def _pdb_atom(
    serial: int,
    atom_name: str,
    residue_name: str,
    chain_id: str,
    residue_number: int,
    x_coord: float,
    y_coord: float,
    z_coord: float,
    *,
    element: str = "C",
    record: str = "ATOM",
) -> str:
    """Format one PDB atom line for tests."""
    return (
        f"{record:<6}{serial:5d} {atom_name:<4} {residue_name:>3} {chain_id:1}"
        f"{residue_number:4d}    {x_coord:8.3f}{y_coord:8.3f}{z_coord:8.3f}"
        f"  1.00  0.00          {element:>2}\n"
    )
