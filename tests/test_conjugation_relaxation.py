"""Tests for conjugate relaxation helpers."""

from __future__ import annotations

import json
import math
import os
import shutil
import sys
from types import SimpleNamespace

import numpy as np
import pytest
from pydantic import ValidationError

import polyzymd.builders.conjugation.relaxation._openmm_workflows as relaxation_workflows
import polyzymd.builders.conjugation.relaxation.openmm as relaxation
from polyzymd.builders.conjugation.relaxation._diagnostics import (
    validate_finite_energy,
    validate_finite_positions,
)
from polyzymd.builders.conjugation.relaxation._linkages import resolve_product_linkage_pairs
from polyzymd.builders.conjugation.relaxation._openmm_system import (
    _add_linkage_anchor_restraints,
    _freeze_protein_chain_masses,
    _remove_barostats,
    _restore_particle_masses,
    _run_fixed_product_md,
    _write_relaxed_product_pdb,
)
from polyzymd.builders.conjugation.relaxation.geometry import (
    positions_to_numpy,
    replace_pdb_coordinates,
)
from polyzymd.builders.conjugation.relaxation.models import (
    ConjugateRelaxationDiagnostics,
    ConjugateRelaxationSettings,
    ProductLinkage,
)
from polyzymd.builders.conjugation.relaxation.openmm import relax_conjugate


class _TopologyDouble:
    """Minimal OpenMM topology-like object for restraint selection tests."""

    def __init__(self, atoms: tuple[object, ...]) -> None:
        """Store atom-like objects returned by ``atoms()``."""
        self._atoms = atoms

    def atoms(self) -> tuple[object, ...]:
        """Return atom-like objects in topology order."""
        return self._atoms

    def bonds(self) -> tuple[tuple[object, object], ...]:
        """Return atom-like topology bonds."""
        return getattr(self, "_bonds", ())


class _AtomDouble:
    """Minimal OpenMM atom-like object for restraint selection tests."""

    def __init__(self, index: int, name: str, element: str, chain_id: str) -> None:
        """Store the atom attributes used by restraint inference."""
        self.index = index
        self.name = name
        self.element = type("ElementDouble", (), {"symbol": element})()
        chain = type("ChainDouble", (), {"id": chain_id})()
        self.residue = type(
            "ResidueDouble",
            (),
            {"chain": chain, "name": "RES", "id": "1"},
        )()


def _pdb_line(
    serial: int | None,
    atom_name: str,
    residue_name: str,
    chain_id: str,
    residue_number: int,
    x: float,
    *,
    element: str | None = None,
) -> str:
    """Build a compact PDB atom line for product-linkage tests."""
    serial_field = "" if serial is None else str(serial)
    atom_element = element or atom_name[0]
    return (
        f"ATOM  {serial_field:>5s} {atom_name:<4s} {residue_name:>3s} {chain_id:1s}"
        f"{residue_number:4d}    {x:8.3f}{0.0:8.3f}{0.0:8.3f}  1.00  0.00"
        f"          {atom_element:>2s}\n"
    )


def _topology_atom(
    index: int,
    atom_name: str,
    residue_name: str,
    chain_id: str,
    residue_id: str,
) -> object:
    """Create a topology atom double with product identity fields."""
    chain = type("ChainDouble", (), {"id": chain_id})()
    residue = type("ResidueDouble", (), {"chain": chain, "name": residue_name, "id": residue_id})()
    return type("TopologyAtomDouble", (), {"index": index, "name": atom_name, "residue": residue})()


def test_validate_finite_energy_accepts_finite_values():
    """Finite energies should pass OpenMM relaxation checks."""
    validate_finite_energy(-123.4, label="test_energy")


def test_write_relaxed_product_pdb_preserves_identity_and_connectivity(tmp_path):
    """Relaxed coordinate export should preserve product PDB identity records."""
    product = tmp_path / "product.pdb"
    relaxed = tmp_path / "relaxed.pdb"
    product.write_text(
        "".join(
            (
                _pdb_line(1, "CA", "ALA", "A", 1, 0.0, element="C"),
                "TER       2      ALA A   1\n",
                _pdb_line(3, "C1", "PEG", "C", 2, 1.0, element="C"),
                "LINK         CA  ALA A   1                 C1  PEG C   2     1555   1555"
                "  1.50\n",
                "CONECT    1    3\n",
                "END\n",
            )
        ),
        encoding="utf-8",
    )
    final_positions_nm = np.array([[0.25, -0.5, 1.0], [1.5, 2.0, -0.25]])
    unit = SimpleNamespace(nanometer=1.0)

    _write_relaxed_product_pdb(
        product,
        final_positions_nm,
        relaxed,
        unit,
        expected_atom_count=2,
    )

    original_lines = product.read_text(encoding="utf-8").splitlines()
    lines = relaxed.read_text(encoding="utf-8").splitlines()
    assert lines[0][:30] == original_lines[0][:30]
    assert lines[2][:30] == original_lines[2][:30]
    assert lines[0][21] == "A"
    assert lines[2][21] == "C"
    assert lines[1].startswith("TER")
    assert lines[3].startswith("LINK")
    assert lines[4] == "CONECT    1    3"
    assert lines[0][30:54] == "   2.500  -5.000  10.000"
    assert lines[2][30:54] == "  15.000  20.000  -2.500"


def test_replace_pdb_coordinates_accepts_valid_fixed_width_edge_fields(tmp_path):
    """Fixed-width coordinate replacement should accept exactly fitting edge fields."""
    product = tmp_path / "product.pdb"
    relaxed = tmp_path / "relaxed.pdb"
    product.write_text(_pdb_line(1, "CA", "ALA", "A", 1, 0.0, element="C"), encoding="utf-8")

    replace_pdb_coordinates(product, np.array([[-999.999, 9999.999, 0.0]]), relaxed)

    line = relaxed.read_text(encoding="utf-8").splitlines()[0]
    assert line[30:54] == "-999.9999999.999   0.000"
    assert all(len(line[start : start + 8]) == 8 for start in (30, 38, 46))


@pytest.mark.parametrize("coordinate", [-1000.0, 9999.9996])
def test_replace_pdb_coordinates_rejects_fixed_width_overflow(tmp_path, coordinate):
    """Coordinate replacement should reject raw and rounded PDB field overflow."""
    product = tmp_path / "product.pdb"
    relaxed = tmp_path / "relaxed.pdb"
    product.write_text(_pdb_line(1, "CA", "ALA", "A", 1, 0.0, element="C"), encoding="utf-8")

    with pytest.raises(ValueError, match="8-character PDB field"):
        replace_pdb_coordinates(product, np.array([[coordinate, 0.0, 0.0]]), relaxed)


def test_write_relaxed_product_pdb_rejects_atom_count_mismatch(tmp_path):
    """Relaxed coordinate export should require one final position per product atom."""
    product = tmp_path / "product.pdb"
    product.write_text(_pdb_line(1, "CA", "ALA", "A", 1, 0.0, element="C"), encoding="utf-8")

    with pytest.raises(ValueError, match="atom count"):
        _write_relaxed_product_pdb(
            product,
            np.array([[0.0, 0.0, 0.0], [0.1, 0.0, 0.0]]),
            tmp_path / "relaxed.pdb",
            SimpleNamespace(nanometer=1.0),
            expected_atom_count=1,
        )


def test_write_relaxed_product_pdb_rejects_nonfinite_positions(tmp_path):
    """Relaxed coordinate export should reject non-finite final coordinates."""
    product = tmp_path / "product.pdb"
    product.write_text(_pdb_line(1, "CA", "ALA", "A", 1, 0.0, element="C"), encoding="utf-8")

    with pytest.raises(ValueError, match="non-finite"):
        _write_relaxed_product_pdb(
            product,
            np.array([[math.nan, 0.0, 0.0]]),
            tmp_path / "relaxed.pdb",
            SimpleNamespace(nanometer=1.0),
            expected_atom_count=1,
        )


def test_resolve_product_linkage_pairs_uses_generic_assembly_metadata(tmp_path):
    """Arbitrary residue and atom names should resolve from CONECT metadata."""
    product = tmp_path / "product.pdb"
    product.write_text(
        "".join(
            (
                _pdb_line(1, "Q1", "ABC", "A", 7, 0.0),
                _pdb_line(2, "Z9", "MNO", "C", 42, 1.4),
                "CONECT    1    2\nEND\n",
            )
        ),
        encoding="utf-8",
    )
    topology = _TopologyDouble(
        (
            _topology_atom(0, "Q1", "ABC", "A", "7"),
            _topology_atom(1, "Z9", "MNO", "C", "42"),
        )
    )
    spec = type(
        "SpecDouble",
        (),
        {
                "target_bond_length_angstrom": 1.4,
                "attachment_id": "x",
                "attachment_index": 1,
                "endpoint_provenance": {},
            "protein_product_residue_name": "ABC",
            "modifier_product_residue_name": "MNO",
            "protein_link_atom": SimpleNamespace(
                serial=1,
                atom_index=0,
                atom_name="Q1",
                residue_name="ABC",
                chain_id="A",
                residue_number=7,
                insertion_code="",
            ),
            "modifier_link_atom": SimpleNamespace(
                serial=2,
                atom_index=1,
                atom_name="Z9",
                residue_name="MNO",
                chain_id="C",
                residue_number=42,
                insertion_code="",
            ),
        },
    )()
    assembly = type("AssemblyDouble", (), {"added_conect_pairs": ((1, 2),)})()

    pairs = resolve_product_linkage_pairs(
        topology,
        product_pdb_path=product,
        attachment_specs=(spec,),
        assembly=assembly,
    )

    assert pairs[0].protein_atom_index == 0
    assert pairs[0].modifier_atom_index == 1
    assert pairs[0].target_bond_length_angstrom == pytest.approx(1.4)


def test_resolve_product_linkage_pairs_disambiguates_duplicate_moieties(tmp_path):
    """Identical duplicated moieties should resolve distinct product linkage pairs."""
    product = tmp_path / "product.pdb"
    product.write_text(
        "".join(
            (
                _pdb_line(1, "Q1", "ABC", "A", 7, 0.0),
                _pdb_line(2, "Z9", "MNO", "C", 42, 1.4),
                _pdb_line(3, "Q1", "ABC", "A", 8, 3.0),
                _pdb_line(4, "Z9", "MNO", "C", 43, 4.4),
                "END\n",
            )
        ),
        encoding="utf-8",
    )
    topology = _TopologyDouble(
        (
            _topology_atom(0, "Q1", "ABC", "A", "7"),
            _topology_atom(1, "Z9", "MNO", "C", "42"),
            _topology_atom(2, "Q1", "ABC", "A", "8"),
            _topology_atom(3, "Z9", "MNO", "C", "43"),
        )
    )
    specs = (
        _generic_product_mapped_relaxation_spec(
            attachment_id="x1",
            attachment_index=1,
            protein_serial=1,
            protein_residue_number=7,
            modifier_serial=2,
            modifier_residue_number=42,
        ),
        _generic_product_mapped_relaxation_spec(
            attachment_id="x2",
            attachment_index=2,
            protein_serial=3,
            protein_residue_number=8,
            modifier_serial=4,
            modifier_residue_number=43,
        ),
    )
    assembly = type("AssemblyDouble", (), {"added_conect_pairs": ((1, 2), (3, 4))})()

    pairs = resolve_product_linkage_pairs(
        topology,
        product_pdb_path=product,
        attachment_specs=specs,
        assembly=assembly,
    )

    assert [(pair.protein_atom_index, pair.modifier_atom_index) for pair in pairs] == [
        (0, 1),
        (2, 3),
    ]


def test_resolve_product_linkage_pairs_validates_distinct_product_mapped_specs(tmp_path):
    """Product-mapped specs should keep distinct ordered serial pairs and targets."""
    product = tmp_path / "product.pdb"
    product.write_text(
        "".join(
            (
                _pdb_line(1, "ND2", "ASX", "A", 42, 0.0),
                _pdb_line(2, "C1", "NAG", "C", 10, 1.4),
                _pdb_line(3, "ND2", "ASX", "A", 87, 3.0),
                _pdb_line(4, "C1", "NAG", "C", 20, 4.6),
                "END\n",
            )
        ),
        encoding="utf-8",
    )
    topology = _TopologyDouble(
        (
            _topology_atom(0, "ND2", "ASX", "A", "42"),
            _topology_atom(1, "C1", "NAG", "C", "10"),
            _topology_atom(2, "ND2", "ASX", "A", "87"),
            _topology_atom(3, "C1", "NAG", "C", "20"),
        )
    )
    assembly = SimpleNamespace(added_conect_pairs=((1, 2), (3, 4)))

    pairs = resolve_product_linkage_pairs(
        topology,
        product_pdb_path=product,
        attachment_specs=_product_mapped_relaxation_specs(),
        assembly=assembly,
    )

    assert [(pair.protein_serial, pair.modifier_serial) for pair in pairs] == [(1, 2), (3, 4)]
    assert [pair.target_bond_length_angstrom for pair in pairs] == [
        pytest.approx(1.4),
        pytest.approx(1.6),
    ]


def test_resolve_product_linkage_pairs_rejects_swapped_product_mapped_pair(tmp_path):
    """A positional assembly pair cannot be swapped away from exact spec endpoints."""
    product = tmp_path / "product.pdb"
    product.write_text(
        "".join(
            (
                _pdb_line(1, "ND2", "ASX", "A", 42, 0.0),
                _pdb_line(2, "C1", "NAG", "C", 10, 1.4),
                _pdb_line(3, "ND2", "ASX", "A", 87, 3.0),
                _pdb_line(4, "C1", "NAG", "C", 20, 4.6),
                "END\n",
            )
        ),
        encoding="utf-8",
    )
    topology = _TopologyDouble(
        (
            _topology_atom(0, "ND2", "ASX", "A", "42"),
            _topology_atom(1, "C1", "NAG", "C", "10"),
            _topology_atom(2, "ND2", "ASX", "A", "87"),
            _topology_atom(3, "C1", "NAG", "C", "20"),
        )
    )
    assembly = SimpleNamespace(added_conect_pairs=((2, 1), (3, 4)))

    with pytest.raises(RuntimeError, match=r"attachment_id='first'.*observed_pair=\(2, 1\)"):
        resolve_product_linkage_pairs(
            topology,
            product_pdb_path=product,
            attachment_specs=_product_mapped_relaxation_specs(),
            assembly=assembly,
        )


def test_resolve_product_linkage_pairs_accepts_valid_mapped_serials(tmp_path):
    """Mapped source serials should resolve when the product atom identity matches."""
    product = tmp_path / "product.pdb"
    product.write_text(
        "".join(
            (
                _pdb_line(10, "ND2", "ASX", "A", 42, 0.0, element="N"),
                _pdb_line(20, "C1", "NAG", "C", 10, 1.4, element="C"),
                "END\n",
            )
        ),
        encoding="utf-8",
    )
    topology = _TopologyDouble(
        (
            _topology_atom(0, "ND2", "ASX", "A", "42"),
            _topology_atom(1, "C1", "NAG", "C", "10"),
        )
    )

    pairs = resolve_product_linkage_pairs(
        topology,
        product_pdb_path=product,
        attachment_specs=(
            _single_product_mapped_relaxation_spec(protein_serial=10, modifier_serial=20),
        ),
    )

    assert (pairs[0].protein_serial, pairs[0].modifier_serial) == (10, 20)
    assert (pairs[0].protein_atom_index, pairs[0].modifier_atom_index) == (0, 1)


def test_resolve_product_linkage_pairs_falls_back_from_stale_serial_collision(tmp_path):
    """A stale serial pointing at another product atom should use unique identity fallback."""
    product = tmp_path / "product.pdb"
    product.write_text(
        "".join(
            (
                _pdb_line(9, "CA", "GLY", "A", 1, -1.0, element="C"),
                _pdb_line(10, "ND2", "ASX", "A", 42, 0.0, element="N"),
                _pdb_line(20, "C1", "NAG", "C", 10, 1.4, element="C"),
                "END\n",
            )
        ),
        encoding="utf-8",
    )
    topology = _TopologyDouble(
        (
            _topology_atom(0, "CA", "GLY", "A", "1"),
            _topology_atom(1, "ND2", "ASX", "A", "42"),
            _topology_atom(2, "C1", "NAG", "C", "10"),
        )
    )

    pairs = resolve_product_linkage_pairs(
        topology,
        product_pdb_path=product,
        attachment_specs=(
            _single_product_mapped_relaxation_spec(protein_serial=9, modifier_serial=20),
        ),
    )

    assert (pairs[0].protein_serial, pairs[0].modifier_serial) == (10, 20)
    assert (pairs[0].protein_atom_index, pairs[0].modifier_atom_index) == (1, 2)


def test_resolve_product_linkage_pairs_falls_back_when_serial_absent(tmp_path):
    """Missing source serial metadata should resolve through one exact identity match."""
    product = tmp_path / "product.pdb"
    product.write_text(
        "".join(
            (
                _pdb_line(10, "ND2", "ASX", "A", 42, 0.0, element="N"),
                _pdb_line(20, "C1", "NAG", "C", 10, 1.4, element="C"),
                "END\n",
            )
        ),
        encoding="utf-8",
    )
    topology = _TopologyDouble(
        (
            _topology_atom(0, "ND2", "ASX", "A", "42"),
            _topology_atom(1, "C1", "NAG", "C", "10"),
        )
    )

    pairs = resolve_product_linkage_pairs(
        topology,
        product_pdb_path=product,
        attachment_specs=(
            _single_product_mapped_relaxation_spec(protein_serial=None, modifier_serial=None),
        ),
    )

    assert (pairs[0].protein_serial, pairs[0].modifier_serial) == (10, 20)


def test_resolve_product_linkage_pairs_rejects_duplicate_product_serial(tmp_path):
    """Duplicate product serials should be rejected instead of accepting a wrong atom."""
    product = tmp_path / "product.pdb"
    product.write_text(
        "".join(
            (
                _pdb_line(10, "ND2", "ASX", "A", 42, 0.0, element="N"),
                _pdb_line(10, "CA", "GLY", "A", 1, -1.0, element="C"),
                _pdb_line(20, "C1", "NAG", "C", 10, 1.4, element="C"),
                "END\n",
            )
        ),
        encoding="utf-8",
    )
    topology = _TopologyDouble(
        (
            _topology_atom(0, "ND2", "ASX", "A", "42"),
            _topology_atom(1, "CA", "GLY", "A", "1"),
            _topology_atom(2, "C1", "NAG", "C", "10"),
        )
    )

    with pytest.raises(RuntimeError, match="duplicate serials=10"):
        resolve_product_linkage_pairs(
            topology,
            product_pdb_path=product,
            attachment_specs=(
                _single_product_mapped_relaxation_spec(protein_serial=10, modifier_serial=20),
            ),
        )


def test_resolve_product_linkage_pairs_rejects_unrelated_duplicate_product_serial(tmp_path):
    """Duplicate product serials should be rejected before serial-to-topology mapping."""
    product = tmp_path / "product.pdb"
    product.write_text(
        "".join(
            (
                _pdb_line(10, "ND2", "ASX", "A", 42, 0.0, element="N"),
                _pdb_line(20, "C1", "NAG", "C", 10, 1.4, element="C"),
                _pdb_line(99, "CA", "GLY", "A", 1, -1.0, element="C"),
                _pdb_line(99, "CB", "GLY", "A", 1, -1.2, element="C"),
                "END\n",
            )
        ),
        encoding="utf-8",
    )
    topology = _TopologyDouble(
        (
            _topology_atom(0, "ND2", "ASX", "A", "42"),
            _topology_atom(1, "C1", "NAG", "C", "10"),
            _topology_atom(2, "CA", "GLY", "A", "1"),
            _topology_atom(3, "CB", "GLY", "A", "1"),
        )
    )

    with pytest.raises(RuntimeError, match="duplicate serials=99"):
        resolve_product_linkage_pairs(
            topology,
            product_pdb_path=product,
            attachment_specs=(
                _single_product_mapped_relaxation_spec(protein_serial=10, modifier_serial=20),
            ),
        )


def test_resolve_product_linkage_pairs_rejects_ambiguous_product_identity(tmp_path):
    """Ambiguous identity fallback should reject instead of guessing one atom."""
    product = tmp_path / "product.pdb"
    product.write_text(
        "".join(
            (
                _pdb_line(10, "ND2", "ASX", "A", 42, 0.0, element="N"),
                _pdb_line(11, "ND2", "ASX", "A", 42, 0.2, element="N"),
                _pdb_line(20, "C1", "NAG", "C", 10, 1.4, element="C"),
                "END\n",
            )
        ),
        encoding="utf-8",
    )
    topology = _TopologyDouble(
        (
            _topology_atom(0, "ND2", "ASX", "A", "42"),
            _topology_atom(1, "ND2", "ASX", "A", "42"),
            _topology_atom(2, "C1", "NAG", "C", "10"),
        )
    )

    with pytest.raises(RuntimeError, match="identity is ambiguous"):
        resolve_product_linkage_pairs(
            topology,
            product_pdb_path=product,
            attachment_specs=(
                _single_product_mapped_relaxation_spec(protein_serial=None, modifier_serial=20),
            ),
        )


def test_resolve_product_linkage_pairs_rejects_multi_attachment_unresolved_pair(tmp_path):
    """Multi-attachment products should not accept positional pairs without endpoint proof."""
    product, topology = _write_mapped_multi_attachment_product(tmp_path)
    specs = _product_mapped_relaxation_specs()
    specs[0].protein_link_atom.residue_number = 999
    assembly = SimpleNamespace(added_conect_pairs=((1, 2), (3, 4)))

    with pytest.raises(RuntimeError, match="expected pair could not be resolved"):
        resolve_product_linkage_pairs(
            topology,
            product_pdb_path=product,
            attachment_specs=specs,
            assembly=assembly,
        )


def test_resolve_product_linkage_pairs_rejects_multi_attachment_missing_pair(tmp_path):
    """Multi-attachment products should require assembly pairs for every attachment."""
    product, topology = _write_mapped_multi_attachment_product(tmp_path)
    assembly = SimpleNamespace(added_conect_pairs=((1, 2),))

    with pytest.raises(RuntimeError, match="assembly metadata"):
        resolve_product_linkage_pairs(
            topology,
            product_pdb_path=product,
            attachment_specs=_product_mapped_relaxation_specs(),
            assembly=assembly,
        )


def test_resolve_product_linkage_pairs_rejects_multi_attachment_invalid_pair(tmp_path):
    """Multi-attachment products should reject malformed assembly pair metadata."""
    product, topology = _write_mapped_multi_attachment_product(tmp_path)
    assembly = SimpleNamespace(added_conect_pairs=((1, 2), "bad"))

    with pytest.raises(RuntimeError, match="observed_pair='bad'"):
        resolve_product_linkage_pairs(
            topology,
            product_pdb_path=product,
            attachment_specs=_product_mapped_relaxation_specs(),
            assembly=assembly,
        )


def test__freeze_protein_chain_masses_supports_chain_a_and_multiple_chains():
    """Temporary zero masses should support chain A and multiple configured chains."""

    class SystemDouble:
        """Minimal system double with mutable particle masses."""

        def __init__(self) -> None:
            self.masses = [12.0, 14.0, 16.0, 32.0]

        def getParticleMass(self, index: int) -> float:
            """Return one particle mass."""
            return self.masses[index]

        def setParticleMass(self, index: int, mass: float) -> None:
            """Set one particle mass."""
            self.masses[index] = mass

    topology = _TopologyDouble(
        (
            _AtomDouble(0, "CA", "C", "A"),
            _AtomDouble(1, "CB", "C", "A"),
            _AtomDouble(2, "C1", "C", "C"),
            _AtomDouble(3, "S1", "S", "D"),
        )
    )
    system = SystemDouble()

    frozen, original = _freeze_protein_chain_masses(
        system,
        topology,
        type("Unit", (), {"dalton": 1.0})(),
        chain_ids=("A", "D"),
    )
    assert frozen == (0, 1, 3)
    assert system.masses == [0.0, 0.0, 16.0, 0.0]

    _restore_particle_masses(system, original)
    assert system.masses == [12.0, 14.0, 16.0, 32.0]


def test__freeze_protein_chain_masses_supports_configurable_chains():
    """Temporary zero masses should support generic configured protein chains."""

    class SystemDouble:
        """Minimal system double with mutable particle masses."""

        def __init__(self) -> None:
            """Store particle masses."""
            self.masses = [12.0, 14.0, 16.0]

        def getParticleMass(self, index: int) -> float:
            """Return one particle mass."""
            return self.masses[index]

        def setParticleMass(self, index: int, mass: float) -> None:
            """Set one particle mass."""
            self.masses[index] = mass

    topology = _TopologyDouble(
        (
            _AtomDouble(0, "CA", "C", "A"),
            _AtomDouble(1, "CB", "C", "B"),
            _AtomDouble(2, "C1", "C", "C"),
        )
    )
    system = SystemDouble()

    frozen, original = _freeze_protein_chain_masses(
        system,
        topology,
        type("Unit", (), {"dalton": 1.0})(),
        chain_ids=("B",),
    )

    assert frozen == (1,)
    assert system.masses == [12.0, 0.0, 16.0]
    _restore_particle_masses(system, original)
    assert system.masses == [12.0, 14.0, 16.0]


def test_temporary_anchor_restraint_count_matches_linkages():
    """Temporary anchor restraints should add one bond per attachment."""

    class ForceDouble:
        """CustomBondForce test double."""

        def __init__(self, expression: str) -> None:
            self.expression = expression
            self.bonds = []

        def addPerBondParameter(self, name: str) -> None:
            """Accept per-bond parameter declarations."""

        def addGlobalParameter(self, name: str, value: float) -> None:
            """Accept global parameter declarations."""

        def addBond(self, atom_i: int, atom_j: int, parameters: list[float]) -> None:
            """Store one anchor bond."""
            self.bonds.append((atom_i, atom_j, parameters))

    class SystemDouble:
        """System test double collecting added forces."""

        def __init__(self) -> None:
            self.forces = []

        def addForce(self, force: ForceDouble) -> None:
            """Store one added force."""
            self.forces.append(force)

    openmm = type("OpenMMDouble", (), {"CustomBondForce": ForceDouble})()
    unit = type("UnitDouble", (), {"kilojoule_per_mole": 1.0, "nanometer": 1.0})()
    pairs = (
        ProductLinkage(
            protein_atom_index=0, modifier_atom_index=1, target_bond_length_angstrom=1.3
        ),
        ProductLinkage(
            protein_atom_index=2, modifier_atom_index=3, target_bond_length_angstrom=1.4
        ),
    )
    system = SystemDouble()

    count = _add_linkage_anchor_restraints(system, pairs, 100.0, openmm, unit)

    assert count == 2
    assert len(system.forces[0].bonds) == 2


def test_conjugate_relaxation_stage_b_uses_fresh_fixed_system(monkeypatch, tmp_path):
    """Stage A should use normal masses and Stage B should use a fresh frozen system."""
    calls: list[tuple[str, str, tuple[float, ...]]] = []
    openmm, openmm_app = _install_fake_openmm(monkeypatch, calls)
    monkeypatch.setattr(
        relaxation_workflows,
        "_select_platform",
        lambda *_args: SimpleNamespace(getName=lambda: "CPU"),
    )
    monkeypatch.setattr(relaxation_workflows, "_assign_force_groups", lambda _system: {})
    monkeypatch.setattr(
        relaxation_workflows, "_force_group_labels", lambda _system, *, existing_labels: {}
    )
    monkeypatch.setattr(relaxation_workflows, "_force_group_energies", lambda *_args: {})
    monkeypatch.setattr(relaxation_workflows, "_add_linkage_anchor_restraints", lambda *_args: 0)
    monkeypatch.setattr(
        relaxation_workflows, "_write_relaxed_product_pdb", lambda *_args, **_kwargs: None
    )

    topology = _relaxation_topology()
    interchange = _RelaxationInterchange(topology)

    result = relax_conjugate(
        interchange,
        tmp_path,
        product_pdb_path=tmp_path / "product.pdb",
        attachment_specs=(),
        settings=ConjugateRelaxationSettings(md_steps=5),
    )

    assert result.success is True
    assert len(interchange.systems) == 2
    assert interchange.systems[0] is not interchange.systems[1]
    assert calls == [
        ("stage_a", "init", (12.0, 16.0)),
        ("stage_b", "init", (0.0, 16.0)),
        ("stage_b", "step", (0.0, 16.0)),
    ]
    assert interchange.systems[0].masses == [12.0, 16.0]
    assert interchange.systems[1].masses == [12.0, 16.0]
    assert sys.modules["openmm"] is openmm
    assert sys.modules["openmm.app"] is openmm_app
    assert openmm.minimize_calls[0]["maxIterations"] == 0


def test_conjugate_relaxation_uses_uncapped_minimization(monkeypatch, tmp_path):
    """Stage A minimization should pass the OpenMM uncapped iteration sentinel."""
    calls: list[tuple[str, str, tuple[float, ...]]] = []
    openmm, _openmm_app = _install_fake_openmm(monkeypatch, calls)
    monkeypatch.setattr(
        relaxation_workflows,
        "_select_platform",
        lambda *_args: SimpleNamespace(getName=lambda: "CPU"),
    )
    monkeypatch.setattr(relaxation_workflows, "_assign_force_groups", lambda _system: {})
    monkeypatch.setattr(
        relaxation_workflows, "_force_group_labels", lambda _system, *, existing_labels: {}
    )
    monkeypatch.setattr(relaxation_workflows, "_force_group_energies", lambda *_args: {})
    monkeypatch.setattr(relaxation_workflows, "_add_linkage_anchor_restraints", lambda *_args: 0)
    monkeypatch.setattr(
        relaxation_workflows, "_write_relaxed_product_pdb", lambda *_args, **_kwargs: None
    )

    relax_conjugate(
        _RelaxationInterchange(_relaxation_topology()),
        tmp_path,
        product_pdb_path=tmp_path / "product.pdb",
        attachment_specs=(),
        settings=ConjugateRelaxationSettings(md_steps=1),
    )

    assert openmm.minimize_calls
    assert openmm.minimize_calls[0]["maxIterations"] == 0


def test_conjugate_relaxation_restores_masses_after_stage_b_error(monkeypatch, tmp_path):
    """Temporary Stage B zero masses should be restored when MD fails."""
    calls: list[tuple[str, str, tuple[float, ...]]] = []
    _install_fake_openmm(monkeypatch, calls)
    monkeypatch.setattr(
        relaxation_workflows,
        "_select_platform",
        lambda *_args: SimpleNamespace(getName=lambda: "CPU"),
    )
    monkeypatch.setattr(relaxation_workflows, "_assign_force_groups", lambda _system: {})
    monkeypatch.setattr(
        relaxation_workflows, "_force_group_labels", lambda _system, *, existing_labels: {}
    )
    monkeypatch.setattr(relaxation_workflows, "_force_group_energies", lambda *_args: {})
    monkeypatch.setattr(relaxation_workflows, "_add_linkage_anchor_restraints", lambda *_args: 0)
    monkeypatch.setattr(
        relaxation_workflows, "_write_relaxed_product_pdb", lambda *_args, **_kwargs: None
    )

    def fail_md(*_args, **_kwargs):
        raise RuntimeError("backend instability")

    monkeypatch.setattr(relaxation_workflows, "_run_fixed_product_md", fail_md)
    interchange = _RelaxationInterchange(_relaxation_topology())

    with pytest.raises(RuntimeError, match="backend instability"):
        relax_conjugate(
            interchange,
            tmp_path,
            product_pdb_path=tmp_path / "product.pdb",
            attachment_specs=(),
            settings=ConjugateRelaxationSettings(md_steps=5),
        )

    assert interchange.systems[1].masses == [12.0, 16.0]


def test_conjugate_relaxation_stage_b_tolerance_violation_fails(monkeypatch, tmp_path):
    """Stage B tolerance violations should fail relaxation evidence."""
    calls: list[tuple[str, str, tuple[float, ...]]] = []
    _install_fake_openmm(monkeypatch, calls)
    monkeypatch.setattr(
        relaxation_workflows,
        "_select_platform",
        lambda *_args: SimpleNamespace(getName=lambda: "CPU"),
    )
    monkeypatch.setattr(relaxation_workflows, "_assign_force_groups", lambda _system: {})
    monkeypatch.setattr(
        relaxation_workflows, "_force_group_labels", lambda _system, *, existing_labels: {}
    )
    monkeypatch.setattr(relaxation_workflows, "_force_group_energies", lambda *_args: {})
    monkeypatch.setattr(relaxation_workflows, "_add_linkage_anchor_restraints", lambda *_args: 0)
    monkeypatch.setattr(
        relaxation_workflows, "_write_relaxed_product_pdb", lambda *_args, **_kwargs: None
    )

    def moving_md(*_args, **_kwargs):
        """Return final coordinates with an immobilized-protein displacement violation."""
        return -1.0, -2.0, np.array([[0.1, 0.0, 0.0], [0.15, 0.0, 0.0]])

    monkeypatch.setattr(relaxation_workflows, "_run_fixed_product_md", moving_md)

    with pytest.raises(RuntimeError, match="Stage B protein RMSD"):
        relax_conjugate(
            _RelaxationInterchange(_relaxation_topology()),
            tmp_path,
            product_pdb_path=tmp_path / "product.pdb",
            attachment_specs=(),
            settings=ConjugateRelaxationSettings(md_steps=5, max_protein_rmsd_angstrom=0.05),
        )

    payload = json.loads((tmp_path / "conjugate_relaxation.json").read_text(encoding="utf-8"))
    assert payload["success"] is False
    assert payload["stage_b_success"] is False
    assert any("Stage B protein RMSD" in warning for warning in payload["warnings"])


def test_conjugate_relaxation_settings_reject_zero_md_steps():
    """Conjugate relaxation must run a real Stage B MD segment."""
    with pytest.raises(ValidationError, match="greater than 0"):
        ConjugateRelaxationSettings(md_steps=0)


def test_frozen_product_md_propagates_programming_errors(monkeypatch):
    """Retry handling should not mask programming errors as instability."""

    class BuggySimulation:
        def __init__(self, *_args):
            raise TypeError("bad fake API")

    warnings: list[str] = []
    openmm = SimpleNamespace(
        LangevinMiddleIntegrator=lambda *_args: object(),
        OpenMMException=RuntimeError,
    )
    app = SimpleNamespace(Simulation=BuggySimulation)
    unit = SimpleNamespace(kelvin=1.0, picosecond=1.0, femtosecond=1.0, nanometer=1.0)

    with pytest.raises(TypeError, match="bad fake API"):
        _run_fixed_product_md(
            _relaxation_topology(),
            _RelaxationSystem("stage_b"),
            np.zeros((2, 3)),
            ConjugateRelaxationSettings(md_steps=5),
            openmm,
            app,
            unit,
            object(),
            warnings,
        )

    assert warnings == []


def test_remove_barostats_removes_only_barostat_forces():
    """Vacuum relaxation should discard barostat forces from transient systems."""

    class MonteCarloBarostat:
        """Barostat marker double."""

    class HarmonicBondForce:
        """Non-barostat marker double."""

    class SystemDouble:
        """System double with removable forces."""

        def __init__(self) -> None:
            self.forces = [HarmonicBondForce(), MonteCarloBarostat(), HarmonicBondForce()]

        def getNumForces(self) -> int:
            """Return force count."""
            return len(self.forces)

        def getForce(self, index: int) -> object:
            """Return force by index."""
            return self.forces[index]

        def removeForce(self, index: int) -> None:
            """Remove force by index."""
            self.forces.pop(index)

    system = SystemDouble()

    assert _remove_barostats(system) == 1
    assert [type(force).__name__ for force in system.forces] == [
        "HarmonicBondForce",
        "HarmonicBondForce",
    ]


def test_conjugate_relaxation_diagnostics_serializes(tmp_path):
    """Conjugate diagnostics should write JSON evidence."""
    path = tmp_path / "diagnostics.json"
    diagnostics = ConjugateRelaxationDiagnostics(
        success=True,
        stage_a_success=True,
        stage_b_success=True,
        frozen_atom_count=2,
        temporary_anchor_count=1,
        stage_a_energy_after_min_kj_mol=-10.0,
        stage_b_energy_after_md_kj_mol=-11.0,
        stage_a_protein_rmsd_from_initial_angstrom=0.2,
        stage_b_protein_rmsd_from_stage_a_angstrom=0.0,
        linkage_distances_angstrom=(1.4,),
        stage_b_linkage_distances_angstrom=(1.4,),
        final_relaxed_pdb_path=tmp_path / "relaxed.pdb",
    )

    diagnostics.write_json(path)

    payload = json.loads(path.read_text(encoding="utf-8"))
    assert payload["success"] is True
    assert payload["stage_a_success"] is True
    assert payload["stage_b_success"] is True
    assert payload["temporary_anchor_count"] == 1
    assert payload["barostat_used"] is False
    assert payload["final_relaxed_pdb_path"] == str(tmp_path / "relaxed.pdb")
    assert payload["json_path"] == str(path)


@pytest.mark.parametrize("value", [math.nan, math.inf, -math.inf])
def test_validate_finite_energy_rejects_nonfinite_values(value: float):
    """NaN and infinite energies should fail hard."""
    with pytest.raises(RuntimeError, match="non-finite test_energy"):
        validate_finite_energy(value, label="test_energy")


def test_validate_finite_positions_accepts_numpy_arrays():
    """Finite coordinate arrays should pass OpenMM relaxation checks."""
    span_nm = validate_finite_positions(np.zeros((3, 3)), label="positions")

    assert span_nm == 0.0


def test_positions_to_numpy_logs_expected_conversion_fallback(caplog):
    """Expected unit conversion failures should fall back with a warning."""

    class RawArrayPositions:
        """Position container with a broken unit API and raw array fallback."""

        def value_in_unit(self, unit):
            """Raise an expected conversion error for the unit-aware path."""
            raise TypeError(f"cannot convert to {unit}")

        def __array__(self, dtype=None):
            """Return raw nanometer coordinates for NumPy coercion."""
            return np.asarray([[0.0, 0.0, 0.0]], dtype=dtype)

    unit_module = type("UnitModule", (), {"nanometer": "nanometer"})()

    with caplog.at_level("WARNING"):
        array = positions_to_numpy(RawArrayPositions(), unit_module)

    np.testing.assert_allclose(array, [[0.0, 0.0, 0.0]])
    assert "Falling back to raw np.asarray()" in caplog.text
    assert "RawArrayPositions" in caplog.text


def test_validate_finite_positions_propagates_unexpected_conversion_errors():
    """Unexpected position conversion failures should not be hidden."""

    class BrokenPositions:
        """Position container with an unexpected conversion failure."""

        def value_in_unit(self, unit):
            """Raise an unexpected error from the unit-aware path."""
            raise RuntimeError(f"unexpected conversion failure for {unit}")

    unit_module = type("UnitModule", (), {"nanometer": "nanometer"})()

    with pytest.raises(RuntimeError, match="unexpected conversion failure"):
        validate_finite_positions(BrokenPositions(), unit_module, label="positions")


def test_validate_finite_positions_rejects_nonfinite_arrays():
    """NaN coordinates should fail hard before downstream use."""
    positions = np.array([[0.0, 0.0, 0.0], [np.nan, 1.0, 2.0]])

    with pytest.raises(RuntimeError, match="non-finite positions"):
        validate_finite_positions(positions, label="positions")


def test_validate_finite_positions_rejects_unrealistic_span():
    """Expanded post-relaxation coordinates should fail before solvation."""
    positions = np.array([[0.0, 0.0, 0.0], [51.0, 1.0, 1.0]])

    with pytest.raises(RuntimeError, match="unrealistic coordinate span"):
        validate_finite_positions(positions, label="equilibrated_positions", max_span_nm=50.0)


class _RelaxationSystem:
    """Mutable system double for staged conjugate relaxation tests."""

    def __init__(self, label: str) -> None:
        """Store a label and particle masses."""
        self.label = label
        self.masses = [12.0, 16.0]

    def getNumForces(self) -> int:  # noqa: N802 - OpenMM API compatibility
        """Return no forces for the minimal relaxation path."""
        return 0

    def getParticleMass(self, index: int) -> float:  # noqa: N802 - OpenMM API compatibility
        """Return one particle mass."""
        return self.masses[index]

    def setParticleMass(self, index: int, mass: float) -> None:  # noqa: N802
        """Set one particle mass."""
        self.masses[index] = float(mass)


class _RelaxationInterchange:
    """Interchange double returning fresh OpenMM systems."""

    def __init__(self, topology: _TopologyDouble) -> None:
        """Store topology and emitted systems."""
        self.topology = topology
        self.positions = np.array([[0.0, 0.0, 0.0], [0.15, 0.0, 0.0]])
        self.systems: list[_RelaxationSystem] = []

    def to_openmm_topology(self) -> _TopologyDouble:
        """Return the fake topology."""
        return self.topology

    def to_openmm_system(self) -> _RelaxationSystem:
        """Return a fresh system for each staged OpenMM conversion."""
        label = "stage_a" if not self.systems else "stage_b"
        system = _RelaxationSystem(label)
        self.systems.append(system)
        return system


class _RelaxationState:
    """OpenMM state double with energy and positions."""

    def __init__(self, positions: np.ndarray) -> None:
        """Store state positions."""
        self._positions = np.asarray(positions, dtype=float)

    def getPotentialEnergy(self) -> float:  # noqa: N802 - OpenMM API compatibility
        """Return a finite energy."""
        return -1.0

    def getPositions(self, asNumpy: bool = False) -> np.ndarray:  # noqa: N802, FBT001, FBT002
        """Return the stored positions."""
        return self._positions.copy()


class _RelaxationContext:
    """OpenMM context double storing positions."""

    def __init__(self) -> None:
        """Initialize empty positions."""
        self.positions = np.zeros((2, 3))

    def setPositions(self, positions: np.ndarray) -> None:  # noqa: N802
        """Store context positions."""
        self.positions = np.asarray(positions, dtype=float)

    def setVelocities(self, _velocities: np.ndarray) -> None:  # noqa: N802
        """Accept velocity initialization."""

    def getState(self, **_kwargs) -> _RelaxationState:  # noqa: N802
        """Return a finite state for requested data."""
        return _RelaxationState(self.positions)


def _relaxation_topology() -> _TopologyDouble:
    """Build a minimal two-atom protein-conjugate topology."""
    return _TopologyDouble(
        (
            _AtomDouble(0, "CA", "C", "A"),
            _AtomDouble(1, "C1", "C", "C"),
        )
    )


def _install_fake_openmm(monkeypatch, calls: list[tuple[str, str, tuple[float, ...]]]):
    """Install fake OpenMM modules for staged relaxation behavior tests."""

    class FakeSimulation:
        """OpenMM Simulation double recording stage masses."""

        def __init__(self, _topology, system: _RelaxationSystem, _integrator, _platform) -> None:
            self.system = system
            self.context = _RelaxationContext()
            calls.append((system.label, "init", tuple(system.masses)))

        def step(self, _steps: int) -> None:
            """Record masses used during MD integration."""
            calls.append((self.system.label, "step", tuple(self.system.masses)))

    minimize_calls = []

    def minimize(*_args, **kwargs) -> None:
        """Record minimization kwargs."""

        minimize_calls.append(kwargs)

    openmm = SimpleNamespace(
        VerletIntegrator=lambda *_args: object(),
        LangevinMiddleIntegrator=lambda *_args: object(),
        LocalEnergyMinimizer=SimpleNamespace(minimize=minimize),
        OpenMMException=RuntimeError,
        minimize_calls=minimize_calls,
    )
    openmm_app = SimpleNamespace(
        Simulation=FakeSimulation,
        PDBFile=SimpleNamespace(writeFile=lambda *_args, **_kwargs: None),
    )
    openmm_unit = SimpleNamespace(
        picoseconds=1.0,
        picosecond=1.0,
        femtosecond=1.0,
        kelvin=1.0,
        kilojoule_per_mole=1.0,
        nanometer=1.0,
        dalton=1.0,
    )
    openmm.app = openmm_app
    openmm.unit = openmm_unit
    monkeypatch.setitem(sys.modules, "openmm", openmm)
    monkeypatch.setitem(sys.modules, "openmm.app", openmm_app)
    monkeypatch.setitem(sys.modules, "openmm.unit", openmm_unit)
    return openmm, openmm_app


def _write_mapped_multi_attachment_product(tmp_path: Path) -> tuple[Path, _TopologyDouble]:
    """Write a two-attachment product PDB and matching topology double."""
    product = tmp_path / "product.pdb"
    product.write_text(
        "".join(
            (
                _pdb_line(1, "ND2", "ASX", "A", 42, 0.0),
                _pdb_line(2, "C1", "NAG", "C", 10, 1.4),
                _pdb_line(3, "ND2", "ASX", "A", 87, 3.0),
                _pdb_line(4, "C1", "NAG", "C", 20, 4.6),
                "END\n",
            )
        ),
        encoding="utf-8",
    )
    topology = _TopologyDouble(
        (
            _topology_atom(0, "ND2", "ASX", "A", "42"),
            _topology_atom(1, "C1", "NAG", "C", "10"),
            _topology_atom(2, "ND2", "ASX", "A", "87"),
            _topology_atom(3, "C1", "NAG", "C", "20"),
        )
    )
    return product, topology


def _generic_product_mapped_relaxation_spec(
    *,
    attachment_id: str,
    attachment_index: int,
    protein_serial: int,
    protein_residue_number: int,
    modifier_serial: int,
    modifier_residue_number: int,
) -> SimpleNamespace:
    """Return one generic mapped spec for duplicate-moiety linkage tests."""
    return SimpleNamespace(
        attachment_id=attachment_id,
        attachment_index=attachment_index,
        endpoint_provenance={},
        protein_link_atom=SimpleNamespace(
            serial=protein_serial,
            chain_id="A",
            residue_name="ABC",
            residue_number=protein_residue_number,
            insertion_code="",
            atom_name="Q1",
        ),
        modifier_link_atom=SimpleNamespace(
            serial=modifier_serial,
            chain_id="C",
            residue_name="MNO",
            residue_number=modifier_residue_number,
            insertion_code="",
            atom_name="Z9",
        ),
        protein_product_residue_name="ABC",
        modifier_product_residue_name="MNO",
        target_bond_length_angstrom=1.4,
    )


def _product_mapped_relaxation_specs() -> tuple[SimpleNamespace, SimpleNamespace]:
    """Return two mapped attachment specs with distinct targets and IDs."""
    return (
        SimpleNamespace(
            attachment_id="first",
            attachment_index=1,
            endpoint_provenance={},
            protein_link_atom=SimpleNamespace(
                serial=1,
                chain_id="A",
                residue_name="ASX",
                residue_number=42,
                insertion_code="",
                atom_name="ND2",
            ),
            modifier_link_atom=SimpleNamespace(
                serial=2,
                chain_id="C",
                residue_name="NAG",
                residue_number=10,
                insertion_code="",
                atom_name="C1",
            ),
            protein_product_residue_name="ASX",
            modifier_product_residue_name="NAG",
            target_bond_length_angstrom=1.4,
        ),
        SimpleNamespace(
            attachment_id="second",
            attachment_index=2,
            endpoint_provenance={},
            protein_link_atom=SimpleNamespace(
                serial=3,
                chain_id="A",
                residue_name="ASX",
                residue_number=87,
                insertion_code="",
                atom_name="ND2",
            ),
            modifier_link_atom=SimpleNamespace(
                serial=4,
                chain_id="C",
                residue_name="NAG",
                residue_number=20,
                insertion_code="",
                atom_name="C1",
            ),
            protein_product_residue_name="ASX",
            modifier_product_residue_name="NAG",
            target_bond_length_angstrom=1.6,
        ),
    )


def _single_product_mapped_relaxation_spec(
    *,
    protein_serial: int | None,
    modifier_serial: int | None,
) -> SimpleNamespace:
    """Return one mapped attachment spec for product-linkage resolution tests."""
    return SimpleNamespace(
        attachment_id="single",
        attachment_index=1,
        endpoint_provenance={},
        protein_link_atom=SimpleNamespace(
            serial=protein_serial,
            chain_id="A",
            residue_name="ASX",
            residue_number=42,
            insertion_code="",
            atom_name="ND2",
            element="N",
        ),
        modifier_link_atom=SimpleNamespace(
            serial=modifier_serial,
            chain_id="C",
            residue_name="NAG",
            residue_number=10,
            insertion_code="",
            atom_name="C1",
            element="C",
        ),
        protein_product_residue_name="ASX",
        modifier_product_residue_name="NAG",
        target_bond_length_angstrom=1.4,
    )


@pytest.mark.slow
@pytest.mark.conjugation_stack
def test_opt_in_conjugation_stack_validation_requirements_are_available():
    """Opt-in guard for stack availability before the integrated validation test.

    The physics acceptance test lives in
    ``tests/test_conjugation_integrated_validation.py`` and should be run under::

        module load slurm/blanca
        salloc ...
        PYTHONNOUSERSITE=1 POLYZYMD_RUN_CONJUGATION_PABLO_VALIDATION=1 \
          pixi run -e sim-cuda-12-4 pytest \
          tests/test_conjugation_integrated_validation.py -v
    """
    if os.environ.get("POLYZYMD_RUN_CONJUGATION_PABLO_VALIDATION") != "1":
        pytest.skip(
            "Set POLYZYMD_RUN_CONJUGATION_PABLO_VALIDATION=1 to check the Pablo validation stack"
        )
    if shutil.which("packmol") is None:
        pytest.skip("Packmol binary is not available on PATH")
    pytest.importorskip("polymerist")
    pytest.importorskip("openff.toolkit")
    pytest.importorskip("openff.interchange")
    pytest.importorskip("openff.pablo")
    pytest.importorskip("openmm")
