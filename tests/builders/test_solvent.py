"""Tests for solvent composition counting helpers."""

import numpy as np
import pytest

from polyzymd.builders.solvent import AVOGADRO_CONSTANT, SolventBuilder
from polyzymd.config.schema import CoSolventSpec, SolventConfig


def test_orthorhombic_padding_preserves_independent_axis_lengths() -> None:
    """Rectangular boxes should add padding per face without cubic expansion."""
    from openff.units import Quantity

    from polyzymd.utils.boxvectors import pad_box_vectors_uniform

    bounding_box = Quantity(np.diag([10.0, 20.0, 30.0]), "angstrom")
    padded = pad_box_vectors_uniform(bounding_box, Quantity(1.0, "nanometer"))

    np.testing.assert_allclose(padded.m_as("angstrom"), np.diag([30.0, 40.0, 50.0]))


@pytest.mark.parametrize("shape", ["orthorhombic", "cube"])
def test_rectangular_shape_labels_use_identity_transform(shape) -> None:
    """Canonical and legacy rectangular labels should share box geometry."""
    import openff.packmol as packmol

    matrix = SolventBuilder()._get_box_shape_matrix(shape)

    np.testing.assert_allclose(matrix, packmol.UNIT_CUBE)


def test_concentration_count_uses_liters_and_avogadro() -> None:
    """Concentration counts should use volume in liters and Avogadro's constant."""
    count = SolventBuilder._count_concentration_molecules(
        concentration_molar=2.0,
        box_volume_liters=1.0e-22,
    )

    assert count == round(2.0 * 1.0e-22 * AVOGADRO_CONSTANT)


def test_mole_fraction_counts_replace_water_from_mass_budget() -> None:
    """Mole-fraction co-solvents should share the neutral solvent mass budget."""
    water_count, cosolvent_counts = SolventBuilder._calculate_mole_fraction_counts(
        neutral_solvent_mass=2400.0,
        water_mass=18.0,
        cosolvent_mole_fractions=[("dmso", 0.10, 78.0)],
    )

    assert cosolvent_counts == [("dmso", 10)]
    assert water_count == 90


def test_neutral_solvent_mass_subtracts_actual_final_ions() -> None:
    """Mass budgeting should reserve the actual final neutralized ion counts."""
    neutral_solvent_mass = SolventBuilder._calculate_neutral_solvent_mass(
        solvent_mass=5000.0,
        na_count=51,
        cl_count=36,
        na_mass=23.0,
        cl_mass=35.5,
    )

    assert neutral_solvent_mass == 5000.0 - (51 * 23.0 + 36 * 35.5)
    assert neutral_solvent_mass != 5000.0 - (43 * (23.0 + 35.5))


def test_config_translation_uses_mole_fraction(monkeypatch) -> None:
    """solvate_from_config should pass mole_fraction into builder composition."""
    captured = {}

    def fake_solvate(self, *, topology, composition, **kwargs):
        """Capture the composition without running the heavy build path."""
        captured["topology"] = topology
        captured["composition"] = composition
        captured["kwargs"] = kwargs
        return topology

    monkeypatch.setattr(SolventBuilder, "solvate", fake_solvate)

    config = SolventConfig(
        co_solvents=[
            CoSolventSpec(name="dmso", mole_fraction=0.10),
            CoSolventSpec(name="urea", concentration=2.0),
        ]
    )
    topology = object()

    result = SolventBuilder().solvate_from_config(topology, config)

    assert result is topology
    co_solvents = captured["composition"].co_solvents
    assert co_solvents[0].mole_fraction == 0.10
    assert co_solvents[0].concentration is None
    assert co_solvents[1].mole_fraction is None
    assert co_solvents[1].concentration == 2.0


def test_neutralizing_ion_counts_use_target_and_tie_toward_more_ions() -> None:
    """Neutralization should target a concentration and prefer more ions on ties."""
    na_count, cl_count = SolventBuilder._calculate_ion_counts(
        nacl_to_add=43,
        solute_charge=-15,
        neutralize=True,
    )

    assert (na_count, cl_count) == (51, 36)
    assert -15 + na_count - cl_count == 0


@pytest.mark.parametrize(
    ("solute_charge", "expected_counts"),
    [
        (-4, (12, 8)),
        (5, (8, 13)),
        (0, (10, 10)),
    ],
)
def test_neutralizing_ion_counts_produce_zero_net_charge(
    solute_charge: int,
    expected_counts: tuple[int, int],
) -> None:
    """Neutralizing final ion counts should exactly cancel integer solute charge."""
    na_count, cl_count = SolventBuilder._calculate_ion_counts(
        nacl_to_add=10,
        solute_charge=solute_charge,
        neutralize=True,
    )

    assert (na_count, cl_count) == expected_counts
    assert solute_charge + na_count - cl_count == 0


def test_non_neutralizing_ion_counts_preserve_equal_salt_pairs() -> None:
    """Non-neutralizing ion counts should preserve equal NaCl pairs."""
    na_count, cl_count = SolventBuilder._calculate_ion_counts(
        nacl_to_add=43,
        solute_charge=-15,
        neutralize=False,
    )

    assert (na_count, cl_count) == (43, 43)


def test_charge_to_integer_tolerates_tiny_floating_noise() -> None:
    """Integer charge conversion should tolerate tiny floating-point noise."""
    assert SolventBuilder._charge_to_integer(-15.00000024) == -15


def test_charge_to_integer_rejects_true_fractional_charge() -> None:
    """Integer charge conversion should reject fractional net charge."""
    with pytest.raises(ValueError, match="must be an integer"):
        SolventBuilder._charge_to_integer(-15.25)
