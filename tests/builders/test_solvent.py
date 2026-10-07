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


@pytest.mark.parametrize(
    ("solute_charge", "expected_counts"),
    [
        (-4, (14, 10)),
        (5, (10, 15)),
        (0, (10, 10)),
    ],
)
def test_neutralizing_ions_are_added_on_top_of_the_salt(
    solute_charge: int,
    expected_counts: tuple[int, int],
) -> None:
    """Neutralizing ions are added to the requested salt pairs, which stay complete."""
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
    """Molecule exposing a single conformer, as get_topology_positions expects."""

    def __init__(self, coords):
        self.conformers = [_FakeConformer(coords)]
        self.n_conformers = 1


class _FakeTopology:
    """Topology exposing an iterable of molecules."""

    def __init__(self, *molecules):
        self.molecules = list(molecules)


def _solute_molecule():
    """A solute spanning 10 x 20 x 30 Angstrom (diameter sqrt(1400) Angstrom)."""
    return _FakeMolecule([[0.0, 0.0, 0.0], [10.0, 20.0, 30.0]])


def _elongated_solute():
    """A filled ellipsoid spanning 37.5 x 41.8 x 50 Angstrom, long along z like T4 lysozyme."""
    import numpy as np

    semi_axes = np.array([18.75, 20.9, 25.0])
    grid = np.mgrid[-25:26:3, -25:26:3, -25:26:3].reshape(3, -1).T.astype(float)
    inside = grid[np.sum((grid / semi_axes) ** 2, axis=1) <= 1.0]
    tips = np.vstack([np.diag(semi_axes), -np.diag(semi_axes)])
    return _FakeMolecule(np.vstack([inside, tips]))


def _closest_image_distance_nm(coords_angstrom, box_nm) -> float:
    """Smallest distance between the solute and any of its periodic copies."""
    import itertools

    import numpy as np
    from scipy.spatial.distance import cdist

    coords = np.asarray(coords_angstrom, dtype=float) / 10.0
    closest = np.inf
    for shift in itertools.product((-1, 0, 1), repeat=3):
        if shift == (0, 0, 0):
            continue
        image = coords + np.asarray(shift, dtype=float) @ box_nm
        closest = min(closest, float(cdist(coords, image).min()))
    return closest


def test_compute_box_vectors_reserves_room_for_polymers() -> None:
    """With polymers the cell edge grows by 2 * packing padding."""
    import numpy as np

    topology = _FakeTopology(_solute_molecule())
    without = SolventBuilder().compute_box_vectors(topology, padding=1.2)
    with_polymers = SolventBuilder().compute_box_vectors(topology, padding=1.2, extra_padding=2.0)

    edge = with_polymers.m_as("nanometer")[0, 0]
    assert edge == pytest.approx(without.m_as("nanometer")[0, 0] + 4.0)
    assert np.linalg.det(with_polymers.m_as("nanometer")) > np.linalg.det(without.m_as("nanometer"))


@pytest.mark.parametrize("shape", ["cube", "rhombic_dodecahedron"])
def test_solute_is_two_paddings_from_its_periodic_copies(shape: str) -> None:
    """Padding is the solute-to-box-edge distance: copies are 2 x padding apart."""
    solute = _elongated_solute()
    box = SolventBuilder().compute_box_vectors(_FakeTopology(solute), padding=1.2, box_shape=shape)

    distance = _closest_image_distance_nm(solute.conformers[0].m_as("angstrom"), box.m_as("nm"))
    assert distance >= 2 * 1.2 - 1e-9


def test_cube_box_is_cubic() -> None:
    """shape: cube gives three equal, orthogonal edges."""
    import numpy as np

    box = SolventBuilder().compute_box_vectors(
        _FakeTopology(_elongated_solute()), padding=1.2, box_shape="cube"
    )
    edge = box.m_as("nanometer")[0, 0]
    np.testing.assert_allclose(box.m_as("nanometer"), edge * np.eye(3))


def test_default_box_fits_an_elongated_solute_in_the_brick() -> None:
    """The default dodecahedron brick leaves at least the Packmol tolerance on every face."""
    import numpy as np

    solute = _elongated_solute()
    config = SolventConfig()
    box = SolventBuilder().compute_box_vectors_from_config(_FakeTopology(solute), config)

    coords = solute.conformers[0].m_as("angstrom") / 10.0
    extent = coords.max(axis=0) - coords.min(axis=0)
    clearance = (np.diagonal(box.m_as("nanometer")) - extent) / 2.0
    assert np.all(clearance >= config.box.tolerance / 10.0 - 1e-9), clearance


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
    """The config wrapper must forward padding, shape, tolerance and the extra padding."""
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
        tolerance=config.box.tolerance,
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


def _packmol_counts(monkeypatch, tmp_path, composition, runs: int = 1) -> list[dict]:
    """Solvate methane `runs` times with Packmol replaced; return what it was asked each time."""
    from openff.toolkit import Molecule, Topology

    import polyzymd.utils.packmol as packmol_utils
    from polyzymd.data import solvent_molecules

    calls: list[dict] = []

    def fake_solvate_with_packmol(**kwargs):
        calls.append(kwargs)
        return kwargs["solute"]

    monkeypatch.setattr(packmol_utils, "solvate_with_packmol", fake_solvate_with_packmol)
    monkeypatch.setattr(solvent_molecules, "_USER_CACHE_DIR", tmp_path)
    monkeypatch.setattr(solvent_molecules, "_loaded_molecules", {})
    solute = Molecule.from_smiles("C")
    solute.generate_conformers(n_conformers=1)
    results = []
    for _ in range(runs):
        SolventBuilder().solvate(Topology.from_molecules([solute]), composition, padding=1.5)
        names = ["water", "na", "cl", "cosolvent"]
        counts = dict(zip(names, calls[-1]["number_of_copies"]))
        counts["molecule"] = calls[-1]["molecules"][3]
        results.append(counts)
    return results


def _solvate_counts(
    monkeypatch, tmp_path, co_solvent_smiles: str, neutralize: bool = True, nacl: float = 0.0
) -> dict:
    """Solvate methane with one co-solvent at 0.5 M, Packmol replaced, and return what it was asked."""
    from polyzymd.builders.solvent import CoSolvent, SolventComposition

    cosolvent = CoSolvent(name="surf", smiles=co_solvent_smiles, concentration=0.5)
    composition = SolventComposition(
        co_solvents=[cosolvent], neutralize=neutralize, nacl_concentration=nacl
    )
    return _packmol_counts(monkeypatch, tmp_path, composition)[0]


class TestCoSolventCharge:
    def test_a_charged_cosolvent_is_neutralized(self, monkeypatch, tmp_path) -> None:
        """An acetate SMILES carries -1 each, which the Na+ count must balance."""
        counts = _solvate_counts(monkeypatch, tmp_path, "CC(=O)[O-]")
        assert counts["cosolvent"] > 0
        assert counts["na"] - counts["cl"] - counts["cosolvent"] == 0

    @pytest.mark.parametrize("neutralize", [True, False])
    def test_a_counter_ion_in_the_smiles_becomes_its_own_ion(
        self, monkeypatch, tmp_path, neutralize
    ) -> None:
        """'...[O-].[Na+]' packs acetate without Na and one Na+ ion per acetate."""
        counts = _solvate_counts(monkeypatch, tmp_path, "CC(=O)[O-].[Na+]", neutralize)
        assert counts["cosolvent"] > 0
        assert counts["na"] == counts["cosolvent"] and counts["cl"] == 0
        assert 11 not in [atom.atomic_number for atom in counts["molecule"].atoms]

    @pytest.mark.parametrize("smiles", ["CC(=O)[O-]", "CC(=O)[O-].[Na+]"])
    def test_requested_salt_is_kept_beside_a_charged_cosolvent(
        self, monkeypatch, tmp_path, smiles
    ) -> None:
        """Both spellings keep every requested Cl-, and the Na+ count adds one per acetate."""
        salt = _solvate_counts(monkeypatch, tmp_path, "CCO", nacl=0.5)["cl"]
        counts = _solvate_counts(monkeypatch, tmp_path, smiles, nacl=0.5)
        assert salt > 0 and counts["cl"] == salt
        assert counts["na"] == salt + counts["cosolvent"]

    def test_a_net_charge_without_neutralize_is_a_warning(
        self, monkeypatch, tmp_path, caplog
    ) -> None:
        """With neutralize off the build goes on, but says the system is charged."""
        import logging

        with caplog.at_level(logging.WARNING):
            _solvate_counts(monkeypatch, tmp_path, "CC(=O)[O-]", neutralize=False)
        assert any("net charge" in r.getMessage() for r in caplog.records)

    @pytest.mark.parametrize("neutralize", [True, False])
    def test_a_reused_composition_gives_the_same_counts(
        self, monkeypatch, tmp_path, neutralize
    ) -> None:
        """A second solvate() with the same composition adds the same counter-ions."""
        from polyzymd.builders.solvent import CoSolvent, SolventComposition

        cosolvent = CoSolvent(name="surf", smiles="CC(=O)[O-].[Na+]", count=3)
        composition = SolventComposition(
            co_solvents=[cosolvent], neutralize=neutralize, nacl_concentration=0.15
        )
        first, second = _packmol_counts(monkeypatch, tmp_path, composition, runs=2)
        assert first["na"] == first["cl"] + 3
        assert [second[k] for k in ("water", "na", "cl", "cosolvent")] == [
            first[k] for k in ("water", "na", "cl", "cosolvent")
        ]

    def test_a_mole_fraction_with_counter_ions_is_met(self, monkeypatch, tmp_path) -> None:
        """Counter-ion mass comes out of the budget before the co-solvent count is taken."""
        from polyzymd.builders.solvent import CoSolvent, SolventComposition

        cosolvent = CoSolvent(name="surf", smiles="CC(=O)[O-].[Na+]", mole_fraction=0.1)
        composition = SolventComposition(co_solvents=[cosolvent], nacl_concentration=0.0)
        counts = _packmol_counts(monkeypatch, tmp_path, composition)[0]
        assert counts["na"] == counts["cosolvent"] and counts["cl"] == 0
        achieved = counts["cosolvent"] / (counts["cosolvent"] + counts["water"])
        assert achieved == pytest.approx(0.1, abs=0.002)
