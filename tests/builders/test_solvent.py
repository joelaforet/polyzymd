"""Tests for solvent composition counting helpers."""

from polyzymd.builders.solvent import AVOGADRO_CONSTANT, SolventBuilder
from polyzymd.config.schema import CoSolventSpec, SolventConfig


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
