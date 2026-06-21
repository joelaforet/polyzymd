"""Tests for conjugation Interchange parameterization helpers."""

from __future__ import annotations

import logging
from types import SimpleNamespace

import numpy as np
import pytest

from polyzymd.builders.conjugation.pablo.parameterization import (
    DEFAULT_CONJUGATION_FORCE_FIELD_NAMES,
    _formal_charge_value,
    _topology_positions_as_angstrom,
    create_interchange_from_openff_topology,
    create_interchange_from_pablo_topology,
    deduplicate_charge_templates,
    load_combined_smirnoff_force_field,
    set_topology_positions_from_pdb,
    suppress_openff_library_charge_info,
)


def test_load_combined_smirnoff_force_field_uses_default_ff14sb_and_sage():
    """The conjugation force field loader should use ff14SB plus Sage."""

    class FakeForceField:
        """Small ForceField stand-in that records constructor arguments."""

        def __init__(self, *names: str):
            self.names = names

    force_field = load_combined_smirnoff_force_field(force_field_cls=FakeForceField)

    assert force_field.names == DEFAULT_CONJUGATION_FORCE_FIELD_NAMES
    assert "ff14sb_off_impropers_0.0.4.offxml" in force_field.names
    assert "openff-2.0.0.offxml" in force_field.names


def test_create_interchange_from_openff_topology_passes_charge_templates():
    """Interchange creation should forward deduplicated charge templates."""
    topology = SimpleNamespace(box_vectors="box")
    interchange = SimpleNamespace(box=None)
    template = _template("[H:1][O:2][H:3]")
    duplicate = _template("[H:1][O:2][H:3]")
    captured: dict[str, object] = {}

    class FakeForceField:
        """Small ForceField stand-in that captures create_interchange kwargs."""

        def create_interchange(self, observed_topology, **kwargs):
            captured["topology"] = observed_topology
            captured["kwargs"] = kwargs
            return interchange

    result = create_interchange_from_openff_topology(
        topology,
        force_field=FakeForceField(),
        charge_from_molecules=(template, duplicate),
    )

    assert result.interchange is interchange
    assert interchange.box == "box"
    assert captured["topology"] is topology
    assert captured["kwargs"] == {"charge_from_molecules": [template]}


def test_create_interchange_from_openff_topology_requires_charge_templates():
    """The direct conjugation path should not allow implicit full-protein charging."""
    with pytest.raises(ValueError, match="requires explicit charged templates"):
        create_interchange_from_openff_topology(
            SimpleNamespace(),
            force_field=SimpleNamespace(create_interchange=lambda topology, **kwargs: object()),
            charge_from_molecules=(),
        )


def test_suppress_openff_library_charge_info_restores_specific_logger_level():
    """Scoped OpenFF log suppression should avoid broad logging side effects."""
    nonbonded_logger = logging.getLogger("openff.interchange.smirnoff._nonbonded")
    polyzymd_logger = logging.getLogger("polyzymd.builders.conjugation.pablo.parameterization")
    previous_nonbonded_level = nonbonded_logger.level
    previous_polyzymd_level = polyzymd_logger.level
    nonbonded_logger.setLevel(logging.INFO)
    polyzymd_logger.setLevel(logging.INFO)
    try:
        with suppress_openff_library_charge_info():
            assert nonbonded_logger.level == logging.WARNING
            assert polyzymd_logger.level == logging.INFO

        assert nonbonded_logger.level == logging.INFO
        assert polyzymd_logger.level == logging.INFO
    finally:
        nonbonded_logger.setLevel(previous_nonbonded_level)
        polyzymd_logger.setLevel(previous_polyzymd_level)


def test_create_interchange_from_pablo_topology_can_require_charge_templates():
    """Pablo topology helper should forward strict charge-template requirements."""
    with pytest.raises(ValueError, match="requires explicit charged templates"):
        create_interchange_from_pablo_topology(
            SimpleNamespace(),
            force_field=SimpleNamespace(create_interchange=lambda topology, **kwargs: object()),
            charge_from_molecules=(),
            require_charge_templates=True,
        )


def test_deduplicate_charge_templates_rejects_uncharged_template():
    """Charge templates must already contain partial charges."""
    with pytest.raises(ValueError, match="partial_charges"):
        deduplicate_charge_templates((_template("CC", partial_charges=None),))


def test_set_topology_positions_from_pdb_preserves_angstrom_magnitudes(tmp_path):
    """PDB Angstrom coordinates should become nanometer OpenFF positions."""
    topology = _water_topology_with_positions(
        np.array(
            [
                [0.0, 0.0, 0.0],
                [1.0, 0.0, 0.0],
                [0.0, 1.0, 0.0],
            ]
        )
    )
    pdb_coordinates = np.array(
        [
            [10.0, 20.0, 30.0],
            [11.5, 20.0, 30.0],
            [10.0, 22.5, 30.0],
        ]
    )
    pdb_path = _write_water_pdb(tmp_path, pdb_coordinates)

    result = set_topology_positions_from_pdb(topology, pdb_path)

    assert result is topology
    np.testing.assert_allclose(topology.get_positions().m_as("angstrom"), pdb_coordinates)
    np.testing.assert_allclose(topology.get_positions().m_as("nanometer"), pdb_coordinates / 10.0)


def test_set_topology_positions_from_pdb_rescales_double_converted_coordinates(tmp_path):
    """Absurd 100x PDB spans should be corrected against existing topology spans."""
    reference_coordinates = np.array(
        [
            [1.0, 2.0, 3.0],
            [2.5, 2.0, 3.0],
            [1.0, 4.0, 3.0],
        ]
    )
    topology = _water_topology_with_positions(reference_coordinates)
    pdb_path = _write_water_pdb(tmp_path, reference_coordinates * 100.0)

    set_topology_positions_from_pdb(topology, pdb_path)

    np.testing.assert_allclose(topology.get_positions().m_as("angstrom"), reference_coordinates)


def test_topology_positions_expected_failure_warns_and_returns_none(caplog):
    """Expected topology position API failures should skip scale inference with a warning."""

    class BrokenTopology:
        """Topology stand-in with unavailable positions."""

        def get_positions(self):
            """Raise an expected API failure."""
            raise ValueError("positions are unavailable")

    with caplog.at_level("WARNING"):
        result = _topology_positions_as_angstrom(BrokenTopology())

    assert result is None
    assert "Could not read topology positions" in caplog.text
    assert "positions are unavailable" in caplog.text


def test_formal_charge_scalar_fallback_warns(caplog):
    """Formal charge scalar fallback should warn after unit conversion fails."""

    class ScalarFallbackCharge:
        """Charge container with broken unit conversion and valid scalar fallback."""

        def m_as(self, unit):
            """Raise an expected unit conversion failure."""
            raise TypeError(f"cannot convert to {unit}")

        def __float__(self):
            """Return the scalar formal charge."""
            return -1.0

    with caplog.at_level("WARNING"):
        charge = _formal_charge_value(ScalarFallbackCharge())

    assert charge == -1.0
    assert "Falling back to scalar formal charge conversion" in caplog.text
    assert "cannot convert" in caplog.text


def test_formal_charge_total_failure_raises_value_error():
    """Formal charge conversion should fail clearly when all paths fail."""

    class BrokenCharge:
        """Charge container with no usable conversion path."""

        def m_as(self, unit):
            """Raise an expected unit conversion failure."""
            raise TypeError(f"cannot convert to {unit}")

        def __float__(self):
            """Raise an expected scalar conversion failure."""
            raise TypeError("not scalar")

    with pytest.raises(ValueError, match="unit-aware conversion or scalar fallback"):
        _formal_charge_value(BrokenCharge())


def _template(smiles: str, *, partial_charges: object | None = (0.0,)):
    """Build a molecule-like template for helper tests."""

    def to_smiles(*, mapped: bool = False) -> str:
        """Return a stable mapped or unmapped SMILES string."""
        return smiles if mapped else smiles.replace(":", "")

    return SimpleNamespace(partial_charges=partial_charges, to_smiles=to_smiles)


def _water_topology_with_positions(coordinates_angstrom: np.ndarray):
    """Build a water topology with reference Angstrom coordinates."""
    pytest.importorskip("openff.toolkit")
    from openff.toolkit import Molecule, Topology
    from openff.units import Quantity

    molecule = Molecule.from_smiles("O")
    topology = Topology.from_molecules([molecule])
    topology.set_positions(Quantity(coordinates_angstrom, "angstrom"))
    return topology


def _write_water_pdb(tmp_path, coordinates_angstrom: np.ndarray):
    """Write a same-order three-atom PDB with supplied Angstrom coordinates."""
    path = tmp_path / "water.pdb"
    atom_names = ("O", "H1", "H2")
    elements = ("O", "H", "H")
    lines = [
        _pdb_atom(index, name, *coords, element=element)
        for index, (name, coords, element) in enumerate(
            zip(atom_names, coordinates_angstrom, elements, strict=True),
            start=1,
        )
    ]
    path.write_text("".join(lines) + "END\n", encoding="utf-8")
    return path


def _pdb_atom(
    index: int, name: str, x_coord: float, y_coord: float, z_coord: float, *, element: str
):
    """Format a minimal fixed-width PDB atom line."""
    return (
        f"ATOM  {index:5d} {name:<4} HOH A   1    "
        f"{x_coord:8.3f}{y_coord:8.3f}{z_coord:8.3f}  1.00  0.00          {element:>2}\n"
    )
