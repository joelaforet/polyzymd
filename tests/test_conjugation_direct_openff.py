"""Tests for the direct OpenFF linkage bridge."""

from __future__ import annotations

import builtins
from pathlib import Path

import pytest

import polyzymd.builders.conjugation.direct_openff as direct_openff
from polyzymd.builders.conjugation.contracts import (
    ExplicitLinkageContract,
    LinkageBond,
    PdbAtomSelector,
    ReactiveEndpoint,
    resolve_explicit_linkage_contract,
)
from polyzymd.builders.conjugation.direct_openff import (
    _normalize_product_formal_charges,
    build_direct_openff_linkage,
)
from polyzymd.builders.conjugation.polymer_fragment import GeneratedPolymerFragment

POC_PROTEIN_PATH = (
    Path(__file__).resolve().parents[1]
    / "src"
    / "polyzymd"
    / "builders"
    / "conjugation"
    / "poc"
    / "data"
    / "NH3_terminal_His_proton_updated.pdb"
)


def test_direct_linkage_proves_bond_removes_leaving_atoms_and_writes_summary(tmp_path: Path):
    """The direct bridge should emit a linked PDB with a proven NZ-RC bond."""
    protein_path = _protein_pdb(tmp_path)
    modifier = _generated_modifier()
    plan = resolve_explicit_linkage_contract(protein_path, modifier, _explicit_contract())

    result = build_direct_openff_linkage(
        protein_pdb_path=protein_path,
        modifier=modifier,
        resolved_plan=plan,
        output_dir=tmp_path / "direct",
        build_openff_topology=False,
    )

    assert result.linked_pdb_path.exists()
    assert result.summary_json_path.exists()
    assert result.linked_bond.protein_atom_name == "NZ"
    assert result.linked_bond.modifier_atom_name == "RC"
    assert result.linked_bond.conect_present is True
    assert result.finite_coordinates is True
    assert result.removed_atom_names == ("HZ2", "HZ3", "LG")
    text = result.linked_pdb_path.read_text(encoding="utf-8")
    assert " HZ2 " not in text
    assert " HZ3 " not in text
    assert " LG  NHX" not in text
    assert "CONECT" in text


def test_direct_linkage_can_skip_heavy_openff_imports(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
):
    """Disabling OpenFF topology construction should avoid heavy imports."""
    protein_path = _protein_pdb(tmp_path)
    modifier = _generated_modifier()
    plan = resolve_explicit_linkage_contract(protein_path, modifier, _explicit_contract())
    original_import = builtins.__import__

    def guarded_import(name, *args, **kwargs):
        """Fail if the disabled path tries to import OpenFF or RDKit."""
        if name.startswith(("openff", "rdkit")):
            raise AssertionError(f"Unexpected heavy import: {name}")
        return original_import(name, *args, **kwargs)

    monkeypatch.setattr(builtins, "__import__", guarded_import)

    result = build_direct_openff_linkage(
        protein_pdb_path=protein_path,
        modifier=modifier,
        resolved_plan=plan,
        output_dir=tmp_path / "direct",
        build_openff_topology=False,
    )

    assert result.topology is None
    assert any("disabled" in limitation for limitation in result.limitations)


def test_product_formal_charge_normalization_marks_tetravalent_nitrogen():
    """The direct topology handoff should repair quaternary nitrogen valence."""
    chem = pytest.importorskip("rdkit.Chem")
    rw_mol = chem.RWMol()
    nitrogen_index = rw_mol.AddAtom(chem.Atom("N"))
    for _ in range(4):
        carbon_index = rw_mol.AddAtom(chem.Atom("C"))
        rw_mol.AddBond(nitrogen_index, carbon_index, chem.BondType.SINGLE)

    _normalize_product_formal_charges(rw_mol)

    assert rw_mol.GetAtomWithIdx(nitrogen_index).GetFormalCharge() == 1


def test_direct_linkage_materializes_topology_and_charge_template_with_sulfonate(
    tmp_path: Path,
):
    """Direct artifacts should expose topology and charges before GPU solvation."""
    pytest.importorskip("openff.toolkit")
    pytest.importorskip("rdkit")
    if not POC_PROTEIN_PATH.exists():
        pytest.fail(f"Required lightweight POC protein fixture is missing: {POC_PROTEIN_PATH}")

    modifier = _sulfonate_modifier()
    plan = resolve_explicit_linkage_contract(
        POC_PROTEIN_PATH,
        modifier,
        _poc_explicit_contract(),
    )

    result = build_direct_openff_linkage(
        protein_pdb_path=POC_PROTEIN_PATH,
        modifier=modifier,
        resolved_plan=plan,
        output_dir=tmp_path / "direct-sulfonate",
    )

    assert result.topology is not None
    assert result.molecule is not None
    assert result.charge_template is not None
    assert result.linked_bond.conect_present is True
    assert not any("component-graph conversion failed" in warning for warning in result.warnings)


def test_direct_openff_conversion_failure_returns_best_effort_warning(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
):
    """Known RDKit/OpenFF conversion failures should preserve linked artifacts."""
    chem = pytest.importorskip("rdkit.Chem")
    pytest.importorskip("openff.toolkit")
    protein_path = _protein_pdb(tmp_path)
    modifier = _generated_modifier()
    plan = resolve_explicit_linkage_contract(protein_path, modifier, _explicit_contract())

    def fail_conversion(*args, **kwargs):
        """Raise the RDKit conversion exception class handled as best-effort."""
        raise chem.rdchem.MolSanitizeException("bad product valence")

    monkeypatch.setattr(
        direct_openff,
        "_build_component_graph_openff_molecule",
        fail_conversion,
    )

    result = build_direct_openff_linkage(
        protein_pdb_path=protein_path,
        modifier=modifier,
        resolved_plan=plan,
        output_dir=tmp_path / "direct-conversion-failure",
    )

    assert result.linked_pdb_path.exists()
    assert result.topology is None
    assert result.molecule is None
    assert any("component-graph conversion failed" in warning for warning in result.warnings)
    assert any("No production-ready OpenFF topology" in item for item in result.limitations)


def test_direct_openff_unexpected_graph_error_propagates_with_context(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
):
    """Unexpected graph assembly errors should not be downgraded to warnings."""
    pytest.importorskip("rdkit.Chem")
    pytest.importorskip("openff.toolkit")
    protein_path = _protein_pdb(tmp_path)
    modifier = _generated_modifier()
    plan = resolve_explicit_linkage_contract(protein_path, modifier, _explicit_contract())

    def fail_unexpectedly(*args, **kwargs):
        """Raise a data-integrity error from the direct assembly path."""
        raise ValueError("linked metadata atom count mismatch")

    monkeypatch.setattr(
        direct_openff,
        "_build_component_graph_openff_molecule",
        fail_unexpectedly,
    )

    with pytest.raises(
        RuntimeError,
        match="Unexpected direct OpenFF topology construction failure",
    ):
        build_direct_openff_linkage(
            protein_pdb_path=protein_path,
            modifier=modifier,
            resolved_plan=plan,
            output_dir=tmp_path / "direct-unexpected-error",
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
        + "END\n",
        encoding="utf-8",
    )
    return path


def _generated_modifier() -> GeneratedPolymerFragment:
    """Create a placed ACB-like modifier fragment with a leaving group."""
    lines = [
        _pdb_atom(101, "C1", "SBM", "C", 1, 3.3, 0.0, 0.0, record="HETATM"),
        _pdb_atom(102, "RC", "NHS", "C", 2, 3.33, 0.0, 0.0, record="HETATM"),
        _pdb_atom(103, "O1", "NHS", "C", 2, 3.8, 0.5, 0.0, element="O", record="HETATM"),
        _pdb_atom(104, "LG", "NHS", "C", 2, 4.2, 1.0, 0.0, element="O", record="HETATM"),
        _pdb_atom(105, "C2", "EGP", "C", 3, 5.2, 1.0, 0.0, record="HETATM"),
    ]
    return GeneratedPolymerFragment.from_pdb_lines(
        lines,
        bonds=((101, 102), (102, 103), (102, 104), (103, 105)),
        reactive_atom_serial=102,
        leaving_atom_serials=(104,),
        sequence="ACB",
        name="ACB",
    )


def _explicit_contract() -> ExplicitLinkageContract:
    """Build an explicit linkage contract for direct bridge tests."""
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


def _poc_explicit_contract() -> ExplicitLinkageContract:
    """Build an explicit Lys23 linkage contract for the lightweight POC protein."""
    protein_selector = PdbAtomSelector(
        chain_id="A",
        residue_name="LYS",
        residue_number=23,
        atom_name="NZ",
    )
    modifier_selector = PdbAtomSelector(
        chain_id="C",
        residue_name="NHX",
        residue_number=1,
        atom_name="RC",
    )
    return ExplicitLinkageContract(
        protein_endpoint=ReactiveEndpoint(
            participant="protein",
            selector=protein_selector,
            product_residue_name="LYX",
            leaving_atom_names=("H11", "H13"),
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


def _sulfonate_modifier() -> GeneratedPolymerFragment:
    """Create a minimal modifier that reproduces the SBMA sulfonate graph."""
    lines = [
        _pdb_atom(101, "RC", "NHX", "C", 1, 6.0, 9.1, 6.6, element="C", record="HETATM"),
        _pdb_atom(102, "O000", "NHX", "C", 1, 6.2, 9.7, 7.6, element="O", record="HETATM"),
        _pdb_atom(103, "LG", "NHX", "C", 1, 5.4, 8.2, 6.6, element="O", record="HETATM"),
        _pdb_atom(104, "C001", "SB1", "C", 2, 7.3, 9.0, 5.9, element="C", record="HETATM"),
        _pdb_atom(105, "H001", "SB1", "C", 2, 7.5, 10.0, 5.5, element="H", record="HETATM"),
        _pdb_atom(106, "H002", "SB1", "C", 2, 7.3, 8.3, 5.0, element="H", record="HETATM"),
        _pdb_atom(107, "C002", "SB1", "C", 2, 8.4, 8.6, 6.9, element="C", record="HETATM"),
        _pdb_atom(108, "H003", "SB1", "C", 2, 8.1, 7.8, 7.6, element="H", record="HETATM"),
        _pdb_atom(109, "H004", "SB1", "C", 2, 8.8, 9.5, 7.4, element="H", record="HETATM"),
        _pdb_atom(110, "S000", "SB1", "C", 2, 9.7, 8.0, 6.0, element="S", record="HETATM"),
        _pdb_atom(111, "O006", "SB1", "C", 2, 10.8, 7.8, 6.9, element="O", record="HETATM"),
        _pdb_atom(112, "O007", "SB1", "C", 2, 9.3, 6.8, 5.3, element="O", record="HETATM"),
        _pdb_atom(113, "O008", "SB1", "C", 2, 10.1, 9.1, 5.0, element="O", record="HETATM"),
    ]
    return GeneratedPolymerFragment.from_pdb_lines(
        lines,
        bond_orders=(
            (101, 102, 2.0),
            (101, 103, 1.0),
            (101, 104, 1.0),
            (104, 105, 1.0),
            (104, 106, 1.0),
            (104, 107, 1.0),
            (107, 108, 1.0),
            (107, 109, 1.0),
            (107, 110, 1.0),
            (110, 111, 2.0),
            (110, 112, 2.0),
            (110, 113, 1.0),
        ),
        reactive_atom_serial=101,
        leaving_atom_serials=(103,),
        sequence="SB",
        name="sulfonate-direct-smoke",
    )


def _pdb_atom(
    serial: int,
    atom_name: str,
    residue_name: str,
    chain: str,
    residue_number: int,
    x: float,
    y: float,
    z: float,
    *,
    element: str = "C",
    record: str = "ATOM",
) -> str:
    """Format a minimal PDB atom record."""
    return (
        f"{record:<6}{serial:5d} {atom_name:<4} {residue_name:>3} {chain}"
        f"{residue_number:4d}    {x:8.3f}{y:8.3f}{z:8.3f}"
        f"  1.00  0.00          {element:>2}\n"
    )
