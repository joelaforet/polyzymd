"""Tests for the specs-first conjugation construction boundary."""

from __future__ import annotations

import importlib.util
from pathlib import Path
from types import SimpleNamespace

import pytest

from polyzymd.builders.conjugation import system_workflow as workflow_module
from polyzymd.builders.conjugation._linkage import (
    ExplicitLinkageContract,
    LinkageBond,
    PdbAtomSelector,
    ReactiveEndpoint,
    resolve_explicit_linkage_contract,
)
from polyzymd.builders.conjugation._specs import attachment_spec_from_generated_polymer_plan
from polyzymd.builders.conjugation.polymer import GeneratedPolymerFragment
from polyzymd.builders.conjugation.structure.pdb import PdbAtomRecord


def test_legacy_assembly_module_and_single_linker_apis_are_removed():
    """Construction should no longer expose obsolete single-linker entry points."""
    assert importlib.util.find_spec("polyzymd.builders.conjugation._assembly") is None
    assert not hasattr(workflow_module, "construct_modifier_linked_protein")
    assert not hasattr(workflow_module, "construct_explicit_pdb_linkage")


def test_attachment_build_spec_is_the_construction_boundary(monkeypatch, tmp_path: Path):
    """Specs-first construction should consume resolved ``AttachmentBuildSpec`` objects."""
    protein_path = _protein_pdb(tmp_path)
    modifier = _generated_modifier()
    plan = resolve_explicit_linkage_contract(
        protein_path,
        modifier,
        _explicit_contract(),
    )
    spec = attachment_spec_from_generated_polymer_plan(
        modifier,
        None,
        plan,
        attachment_config=SimpleNamespace(name="spec-only"),
        attachment_index=1,
        reaction_name="explicit_linkage",
    )
    observed: dict[str, object] = {}

    def fake_construct(**kwargs):
        """Capture the construction input without running heavy boundaries."""
        observed.update(kwargs)
        return object(), object()

    monkeypatch.setattr(workflow_module, "_construct_conjugate_from_specs", fake_construct)
    result = workflow_module._construct_conjugate_from_specs(
        protein_pdb_path=protein_path,
        specs=(spec,),
        ccd_pablo_policy=SimpleNamespace(crosslinks=[]),
        output_dir=tmp_path / "out",
        chain_policy=None,
        settings=workflow_module.ModifierConstructionSettings(),
        use_product_state_pablo_library=False,
    )

    assert result == (observed and result[0], result[1])
    assert observed["specs"] == (spec,)


def test_construct_conjugate_from_specs_rejects_empty_specs(tmp_path: Path):
    """Specs-first construction should require at least one attachment spec."""
    with pytest.raises(ValueError, match="at least one attachment spec"):
        workflow_module._construct_conjugate_from_specs(
            protein_pdb_path=tmp_path / "protein.pdb",
            specs=(),
            ccd_pablo_policy=SimpleNamespace(crosslinks=[]),
            output_dir=tmp_path / "out",
            chain_policy=None,
            settings=workflow_module.ModifierConstructionSettings(),
            use_product_state_pablo_library=False,
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


def _explicit_contract() -> ExplicitLinkageContract:
    """Build a generic explicit linkage contract for specs-first tests."""
    return ExplicitLinkageContract(
        protein_endpoint=ReactiveEndpoint(
            participant="protein",
            selector=PdbAtomSelector(
                chain_id="A",
                residue_name="LYS",
                residue_number=23,
                atom_name="NZ",
            ),
            product_residue_name="LYX",
            leaving_atom_names=("HZ2", "HZ3"),
        ),
        modifier_endpoint=ReactiveEndpoint(
            participant="modifier",
            selector=PdbAtomSelector(
                chain_id="C",
                residue_name="NHS",
                residue_number=2,
                atom_name="RC",
            ),
            product_residue_name="NHX",
            leaving_atom_names=("LG",),
        ),
        bond=LinkageBond(
            protein_atom_name="NZ",
            modifier_atom_name="RC",
            bond_order=1,
            target_bond_length_angstrom=1.33,
        ),
        mechanism_name="explicit_linkage",
    )


def _generated_modifier() -> GeneratedPolymerFragment:
    """Create a small generated modifier fragment."""
    return GeneratedPolymerFragment.from_atom_records(
        (
            PdbAtomRecord(
                record_name="HETATM",
                serial=101,
                atom_name="C1",
                residue_name="SBM",
                chain_id="C",
                residue_number=1,
                x=5.0,
                y=0.0,
                z=0.0,
                element="C",
                atom_index=0,
            ),
            PdbAtomRecord(
                record_name="HETATM",
                serial=102,
                atom_name="RC",
                residue_name="NHS",
                chain_id="C",
                residue_number=2,
                x=3.3,
                y=0.0,
                z=0.0,
                element="C",
                atom_index=1,
            ),
            PdbAtomRecord(
                record_name="HETATM",
                serial=103,
                atom_name="LG",
                residue_name="NHS",
                chain_id="C",
                residue_number=2,
                x=4.2,
                y=1.0,
                z=0.0,
                element="O",
                atom_index=2,
            ),
        ),
        bonds=((101, 102), (102, 103)),
        reactive_atom_serial=102,
        reactive_atom_index=1,
        leaving_atom_serials=(103,),
        leaving_atom_indices=(2,),
        name="mock_modifier",
    )


def _pdb_atom(
    serial: int,
    name: str,
    residue: str,
    chain: str,
    residue_number: int,
    x_coord: float,
    y_coord: float,
    z_coord: float,
    *,
    element: str = "C",
) -> str:
    """Format a PDB atom line."""
    return (
        f"ATOM  {serial:5d} {name:<4s} {residue:>3s} {chain:1s}{residue_number:4d}    "
        f"{x_coord:8.3f}{y_coord:8.3f}{z_coord:8.3f}  1.00  0.00          {element:>2s}\n"
    )
