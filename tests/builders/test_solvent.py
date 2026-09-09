"""Tests for solvent composition counting helpers."""

import pytest

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


def test_config_translation_forwards_packmol_seed(monkeypatch) -> None:
    """solvate_from_config should forward the Packmol seed to solvate()."""
    captured = {}

    def fake_solvate(self, *, topology, composition, **kwargs):
        captured["kwargs"] = kwargs
        return topology

    monkeypatch.setattr(SolventBuilder, "solvate", fake_solvate)
    SolventBuilder().solvate_from_config(object(), SolventConfig(), seed=4)
    assert captured["kwargs"]["seed"] == 4


def test_solvate_default_seed_is_none(monkeypatch) -> None:
    """Without a seed, solvate_from_config must not invent one."""
    captured = {}

    def fake_solvate(self, *, topology, composition, **kwargs):
        captured["kwargs"] = kwargs
        return topology

    monkeypatch.setattr(SolventBuilder, "solvate", fake_solvate)
    SolventBuilder().solvate_from_config(object(), SolventConfig())
    assert captured["kwargs"]["seed"] is None


# ---------------------------------------------------------------------------
# Deterministic periodic cell
# ---------------------------------------------------------------------------


class _FakeConformer:
    """Minimal stand-in for an OpenFF conformer quantity."""

    def __init__(self, coords):
        import numpy as np

        self._coords = np.asarray(coords, dtype=float)

    def m_as(self, _unit):
        return self._coords


class _FakeMolecule:
    """Molecule exposing a single conformer, as get_topology_bbox expects."""

    def __init__(self, coords):
        self.conformers = [_FakeConformer(coords)]
        self.n_conformers = 1


class _FakeTopology:
    """Topology exposing an iterable of molecules."""

    def __init__(self, *molecules):
        self.molecules = list(molecules)


def _solute_molecule():
    """A solute spanning 10 x 20 x 30 Angstrom."""
    return _FakeMolecule([[0.0, 0.0, 0.0], [10.0, 20.0, 30.0]])


def _legacy_box_vectors(topology, padding_nm):
    """The box the pre-fix code computed: bbox + 2*padding, shaped."""
    import openff.packmol as packmol
    from openff.units import Quantity

    from polyzymd.utils import boxvectors

    bbox = boxvectors.get_topology_bbox(topology)
    padded = boxvectors.pad_box_vectors_uniform(bbox, Quantity(padding_nm, "nanometer"))
    return packmol.RHOMBIC_DODECAHEDRON @ padded


def test_compute_box_vectors_matches_legacy_for_a_solute_only_build() -> None:
    """Control bundles must keep the box they have today (no extra padding)."""
    import numpy as np

    topology = _FakeTopology(_solute_molecule())
    computed = SolventBuilder().compute_box_vectors(topology, padding=1.2)
    legacy = _legacy_box_vectors(topology, 1.2)

    np.testing.assert_allclose(
        computed.m_as("nanometer"), legacy.m_as("nanometer"), rtol=0, atol=1e-12
    )


def test_compute_box_vectors_reserves_room_for_polymers() -> None:
    """With polymers the cell grows by 2 * packing padding in every direction."""
    import numpy as np

    topology = _FakeTopology(_solute_molecule())
    without = SolventBuilder().compute_box_vectors(topology, padding=1.2)
    with_polymers = SolventBuilder().compute_box_vectors(topology, padding=1.2, extra_padding=2.0)
    grown = _legacy_box_vectors(topology, 1.2 + 2.0)

    np.testing.assert_allclose(
        with_polymers.m_as("nanometer"), grown.m_as("nanometer"), rtol=0, atol=1e-12
    )
    assert np.linalg.det(with_polymers.m_as("nanometer")) > np.linalg.det(without.m_as("nanometer"))


def test_compute_box_vectors_is_independent_of_polymer_positions() -> None:
    """Two replicates whose chains landed elsewhere must share one box."""
    import numpy as np

    solute = _solute_molecule()
    replicate_a = _FakeTopology(solute, _FakeMolecule([[40.0, 40.0, 40.0]]))
    replicate_b = _FakeTopology(solute, _FakeMolecule([[-60.0, 90.0, 15.0]]))

    builder = SolventBuilder()
    box_a = builder.compute_box_vectors(_FakeTopology(solute), padding=1.2, extra_padding=2.0)
    box_b = builder.compute_box_vectors(_FakeTopology(solute), padding=1.2, extra_padding=2.0)
    np.testing.assert_array_equal(box_a.m_as("nanometer"), box_b.m_as("nanometer"))

    # ...whereas deriving the box from the packed topology does not.
    packed_a = builder.compute_box_vectors(replicate_a, padding=1.2)
    packed_b = builder.compute_box_vectors(replicate_b, padding=1.2)
    assert not np.allclose(packed_a.m_as("nanometer"), packed_b.m_as("nanometer"))


def test_compute_box_vectors_from_config_uses_config_padding() -> None:
    """The config wrapper must forward padding, shape and the extra padding."""
    import numpy as np

    topology = _FakeTopology(_solute_molecule())
    config = SolventConfig()
    computed = SolventBuilder().compute_box_vectors_from_config(
        topology, config, extra_padding_nm=2.0
    )
    expected = SolventBuilder().compute_box_vectors(
        topology,
        padding=config.box.padding,
        box_shape=config.box.shape.value,
        extra_padding=2.0,
    )
    np.testing.assert_array_equal(computed.m_as("nanometer"), expected.m_as("nanometer"))


def _real_solute_topology():
    """A real one-molecule OpenFF topology with coordinates (methane)."""
    from openff.toolkit import Molecule, Topology

    molecule = Molecule.from_smiles("C")
    molecule.generate_conformers(n_conformers=1)
    return Topology.from_molecules([molecule])


def test_solvate_with_precomputed_box_skips_centring(monkeypatch) -> None:
    """A topology already framed in the brick must not be moved again."""
    import numpy as np

    import polyzymd.utils.packmol as packmol_utils

    captured = {}

    def fake_solvate_with_packmol(**kwargs):
        captured.update(kwargs)
        return kwargs["solute"]

    def fail_if_centred(self, topology, box_vecs):
        raise AssertionError("an already-framed topology must not be re-centred")

    monkeypatch.setattr(packmol_utils, "solvate_with_packmol", fake_solvate_with_packmol)
    monkeypatch.setattr(SolventBuilder, "_center_topology_in_box", fail_if_centred)

    topology = _real_solute_topology()
    box = SolventBuilder().compute_box_vectors(topology, padding=1.2, extra_padding=2.0)

    SolventBuilder().solvate(topology=topology, box_vectors=box)

    assert captured["center_solute"] is False
    np.testing.assert_array_equal(captured["box_vectors"].m_as("nanometer"), box.m_as("nanometer"))


def test_solvate_without_box_centres_and_derives_the_box(monkeypatch) -> None:
    """Solute-only builds keep the legacy compute-centre-solvate behaviour."""
    import numpy as np

    import polyzymd.utils.packmol as packmol_utils

    captured = {}
    centred = []

    def fake_solvate_with_packmol(**kwargs):
        captured.update(kwargs)
        return kwargs["solute"]

    monkeypatch.setattr(packmol_utils, "solvate_with_packmol", fake_solvate_with_packmol)
    monkeypatch.setattr(
        SolventBuilder,
        "_center_topology_in_box",
        lambda self, topology, box_vecs: centred.append(box_vecs),
    )

    topology = _real_solute_topology()
    expected = SolventBuilder().compute_box_vectors(topology, padding=1.2)

    SolventBuilder().solvate(topology=topology, padding=1.2)

    assert captured["center_solute"] is True
    assert len(centred) == 1
    np.testing.assert_allclose(
        captured["box_vectors"].m_as("nanometer"), expected.m_as("nanometer")
    )
